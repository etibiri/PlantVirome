import json
import random
import shutil

import pytest
from simulate import rand_seq, simulate, write_fastq

from plantvirome.config import Config
from plantvirome.steps.hostdeplete import host_read_ids, run_hostdeplete


def test_union_coverage(tmp_path):
    # two overlapping hits cover 0-600 of a 1000-bp read -> host at 0.5, not at 0.7
    paf = tmp_path / "x.paf"
    paf.write_text("r1\t1000\t0\t400\t+\tc\t9\t0\t400\t380\t400\t60\n"
                   "r1\t1000\t300\t600\t+\tc\t9\t0\t300\t290\t300\t60\n")
    assert host_read_ids(paf, 0.5, 0.8) == {"r1"}
    assert host_read_ids(paf, 0.7, 0.8) == set()


@pytest.mark.skipif(not shutil.which("minimap2"), reason="minimap2 not on PATH")
@pytest.mark.parametrize("with_host_ref", [True, False])
def test_hostdeplete_separates_host(tmp_path, with_host_ref):
    rng = random.Random(0)
    host = rand_seq(rng, 20000)
    reads, truth = simulate([rand_seq(rng, 2700)], n_reads=150, host=host)
    qc = tmp_path / "s1" / "qc"
    qc.mkdir(parents=True)
    write_fastq(qc / "reads.qc.fastq.gz", reads)
    refs = []
    if with_host_ref:
        (tmp_path / "host.fa").write_text(f">h\n{host}\n")
        refs = [str(tmp_path / "host.fa")]
    cfg = Config.model_validate({"host": {"references": refs}})
    run_hostdeplete(cfg, tmp_path / "s1", 2)
    stats = json.loads((tmp_path / "s1" / "hostdeplete" / "stats.json").read_text())
    n_host = sum(v == "host" for v in truth.values())
    if with_host_ref:
        assert abs(stats["host"] - n_host) <= 0.1 * n_host  # 5% errors, near-perfect recall
    else:
        assert stats["host"] == 0 and stats["nonhost"] == len(reads)
