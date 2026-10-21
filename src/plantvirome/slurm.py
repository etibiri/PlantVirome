"""Generate sbatch scripts, one per step, chained with afterok dependencies."""
import shlex
import subprocess
from pathlib import Path

from .config import COHORT_STEPS, STEPS, Config


def render(cfg: Config, step: str, sample: str | None, config_path: Path,
           outdir: Path, logdir: Path) -> str:
    res = cfg.resources(step)
    name = f"pv_{step}" + (f"_{sample}" if sample else "")
    lines = ["#!/bin/bash",
             f"#SBATCH --job-name={name}",
             f"#SBATCH --cpus-per-task={res.cpus}",
             f"#SBATCH --mem={res.mem_gb}G",
             f"#SBATCH --time={res.time}",
             f"#SBATCH --output={logdir}/{name}.%j.log"]
    if cfg.slurm.partition:
        lines.append(f"#SBATCH --partition={cfg.slurm.partition}")
    if cfg.slurm.account:
        lines.append(f"#SBATCH --account={cfg.slurm.account}")
    lines += ["set -euo pipefail", *cfg.slurm.setup]
    cmd = ["plantvirome", step, "--config", str(config_path), "--outdir", str(outdir),
           "--threads", str(res.cpus)]
    if sample:
        cmd += ["--sample", sample]
    lines.append(shlex.join(cmd))
    return "\n".join(lines) + "\n"


def write_scripts(cfg: Config, samples: list[str], steps: list[str], config_path: Path,
                  outdir: Path) -> list[tuple[str, str | None, Path]]:
    """Write scripts in execution order; returns (step, sample, path) tuples."""
    unknown = set(steps) - set(STEPS)
    if unknown:
        raise ValueError(f"unknown steps: {sorted(unknown)}")
    jobdir = outdir / "slurm"
    logdir = jobdir / "logs"
    logdir.mkdir(parents=True, exist_ok=True)
    scripts = []
    for step in sorted(steps, key=STEPS.index):
        for sample in ([None] if step in COHORT_STEPS else samples):
            path = jobdir / (f"{step}.sh" if sample is None else f"{step}.{sample}.sh")
            path.write_text(render(cfg, step, sample, config_path, outdir, logdir))
            scripts.append((step, sample, path))
    return scripts


def submit(scripts: list[tuple[str, str | None, Path]]) -> None:
    """sbatch every script; a job waits for all jobs of the previous step (afterok),
    so per-sample chains stay simple and a cohort step waits for every sample."""
    prev: list[str] = []
    current_step, current = None, []
    for step, _sample, path in scripts:
        if step != current_step:
            prev, current, current_step = current, [], step
        cmd = ["sbatch", "--parsable"]
        if prev:
            cmd.append("--dependency=afterok:" + ":".join(prev))
        out = subprocess.run(cmd + [str(path)], capture_output=True, text=True, check=True)
        current.append(out.stdout.strip().split(";")[0])
