import json
import random
import shutil

import pytest
from simulate import rand_seq, simulate, write_fastq

from plantvirome.config import Config
from plantvirome.steps.qc import length_stats, run_qc


def test_length_stats():
    assert length_stats([10, 10, 20, 60])["n50"] == 60
    assert length_stats([])["reads"] == 0


@pytest.mark.needs_tools
@pytest.mark.skipif(not (shutil.which("chopper") and shutil.which("porechop_abi")),
                    reason="chopper/porechop_abi not on PATH")
def test_qc_end_to_end(tmp_path):
    rng = random.Random(0)
    reads, _ = simulate([rand_seq(rng, 2700)], n_reads=100)
    (tmp_path / "s1" / "input").mkdir(parents=True)
    write_fastq(tmp_path / "s1" / "input" / "a.fastq.gz", reads)
    cfg = Config.model_validate({"qc": {"min_length": 1000, "min_quality": 10}})
    run_qc(cfg, tmp_path / "s1", None, 2)
    stats = json.loads((tmp_path / "s1" / "qc" / "stats.json").read_text())
    assert 0 < stats["filtered"]["reads"] <= stats["raw"]["reads"]
    assert (tmp_path / "s1" / "qc" / ".done").exists()
