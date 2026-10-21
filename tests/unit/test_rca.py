import json
import random
import shutil

import pytest
from simulate import rand_seq, simulate, write_fastq

from plantvirome.config import Config
from plantvirome.steps.rca import classify, run_rca

R = {"readlen": 1000, "start": 0, "end": 900}


def test_classify():
    assert classify([], 0.8) == "untandem"
    assert classify([R], 0.8) == "clean"
    assert classify([{**R, "end": 500}], 0.8) == "partial"
    assert classify([R, R], 0.8) == "multi"


@pytest.mark.skipif(not shutil.which("TideHunter"), reason="TideHunter not on PATH")
def test_rca_units_and_chimeras(tmp_path):
    rng = random.Random(0)
    glen = 2700
    reads, truth = simulate([rand_seq(rng, glen), rand_seq(rng, 2400)], n_reads=120,
                            host_frac=0, chimera_frac=0.25, copies=(3, 8))
    hd = tmp_path / "s1" / "hostdeplete"
    hd.mkdir(parents=True)
    write_fastq(hd / "reads.nonhost.fastq.gz", reads)
    run_rca(Config(), tmp_path / "s1", 2)
    cls = dict(l.split("\t")[:2] for l in
               (tmp_path / "s1" / "rca" / "reads.tsv").read_text().splitlines()[1:])
    chim = [n for n, t in truth.items() if t == "chimera"]
    virus = [n for n, t in truth.items() if t != "chimera"]
    # chimeras must come out flagged (multi/partial) or untandem, almost never clean
    assert sum(cls[n] == "clean" for n in virus) >= 0.9 * len(virus)
    assert sum(cls[n] == "clean" for n in chim) <= 0.1 * len(chim)
    sizes = {len(l) for l in (tmp_path / "s1" / "rca" / "units.fasta").read_text().split("\n")[1::2] if l}
    assert all(abs(s - glen) < 0.03 * glen or abs(s - 2400) < 0.03 * 2400 for s in sizes)
