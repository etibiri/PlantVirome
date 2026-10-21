"""Step 1: adapter trimming (porechop_abi), filtering (chopper), read statistics."""
import gzip
import json
import shutil
from pathlib import Path

from ..config import Config
from ..logging import event
from ..runner import ToolError, mark_done, run, step_done


def find_reads(paths: list[Path]) -> list[Path]:
    files: list[Path] = []
    for p in paths:
        files += sorted(f for f in p.iterdir() if ".fastq" in f.name or ".fq" in f.name) \
            if p.is_dir() else [p]
    if not files:
        raise ToolError(f"no FASTQ files found in {paths}")
    return files


def concat(files: list[Path], dest: Path) -> None:
    """Byte-concatenate FASTQ files into one gzip file (gzip members concatenate
    validly; plain files are compressed first)."""
    with open(dest, "wb") as out:
        for f in files:
            if f.suffix == ".gz":
                with open(f, "rb") as fh:
                    shutil.copyfileobj(fh, out)
            else:
                with open(f, "rb") as fh, gzip.GzipFile(fileobj=out, mode="wb") as gz:
                    shutil.copyfileobj(fh, gz)


def read_lengths(fastq_gz: Path) -> list[int]:
    with gzip.open(fastq_gz, "rt") as fh:
        return [len(line) - 1 for i, line in enumerate(fh) if i % 4 == 1]


def length_stats(lengths: list[int]) -> dict:
    total = sum(lengths)
    n50, acc = 0, 0
    for length in sorted(lengths, reverse=True):
        acc += length
        if acc * 2 >= total:
            n50 = length
            break
    return {"reads": len(lengths), "bases": total, "n50": n50,
            "mean": round(total / len(lengths), 1) if lengths else 0,
            "max": max(lengths, default=0)}


def run_qc(cfg: Config, sampledir: Path, reads: list[Path] | None, threads: int,
           force: bool = False) -> None:
    out = sampledir / "qc"
    if step_done(out) and not force:
        event("qc skipped (done)", sample=sampledir.name)
        return
    out.mkdir(parents=True, exist_ok=True)
    files = find_reads(reads or [sampledir / "input"])
    raw, trimmed, clean = out / "raw.fastq.gz", out / "trimmed.fastq.gz", out / "reads.qc.fastq.gz"
    concat(files, raw)
    src = raw
    if cfg.qc.trim_adapters:
        run(["porechop_abi", "-abi", "--no_split", "-i", raw, "-o", trimmed, "-t", threads])
        src = trimmed
    cmd = ["chopper", "-q", cfg.qc.min_quality, "-l", cfg.qc.min_length, "-t", threads, "-i", src]
    if cfg.qc.max_length:
        cmd += ["--maxlength", cfg.qc.max_length]
    # chopper writes plain FASTQ to stdout; compress while writing.
    plain = out / "reads.qc.fastq"
    run(cmd, stdout=plain)
    with open(plain, "rb") as fh, gzip.open(clean, "wb") as gz:
        shutil.copyfileobj(fh, gz)
    plain.unlink()
    trimmed.unlink(missing_ok=True)
    stats = {"raw": length_stats(read_lengths(raw)), "filtered": length_stats(read_lengths(clean))}
    (out / "stats.json").write_text(json.dumps(stats, indent=2))
    raw.unlink()
    event("qc stats", sample=sampledir.name, **stats["filtered"])
    mark_done(out)
