import pytest
from pydantic import ValidationError
from typer.testing import CliRunner

from plantvirome.cli import app
from plantvirome.config import Config, load_config
from plantvirome.runner import ToolError, run


def test_default_config_loads(tmp_path):
    cfg = load_config("config/default.yaml")
    assert cfg.qc.min_length == 200


def test_unknown_key_rejected():
    with pytest.raises(ValidationError):
        Config.model_validate({"qcc": {}})


def test_run_missing_tool():
    with pytest.raises(ToolError):
        run(["definitely_not_a_tool"])


def test_submit_dry_run_chain(tmp_path):
    res = CliRunner().invoke(app, ["submit", "--config", "config/default.yaml",
                                   "--outdir", str(tmp_path), "--sample", "s1",
                                   "--sample", "s2", "--steps", "qc,diversity", "--dry-run"])
    assert res.exit_code == 0, res.output
    lines = res.output.strip().splitlines()
    assert len(lines) == 3  # qc x2 samples + 1 cohort diversity job
    script = (tmp_path / "slurm" / "qc.s1.sh").read_text()
    assert "--sample s1" in script and "#SBATCH --cpus-per-task=4" in script


def test_slurm_commands_match_cli():
    """Every command line that `submit` generates must use options the step accepts."""
    import re

    from plantvirome.config import COHORT_STEPS, STEPS
    from plantvirome.slurm import render
    cfg = load_config("config/default.yaml")
    from pathlib import Path
    for step in STEPS:
        sample = None if step in COHORT_STEPS else "s1"
        script = render(cfg, step, sample, Path("c.yaml"), Path("out"), Path("logs"))
        used = set(re.findall(r"--[a-z-]+", script.splitlines()[-1]))
        helptext = CliRunner().invoke(app, [step, "--help"]).output
        assert all(opt in helptext for opt in used), (step, used)
