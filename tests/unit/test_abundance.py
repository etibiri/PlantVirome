import random
import shutil

import pytest
from simulate import rand_seq, simulate, write_fastq

from plantvirome.config import Config
from plantvirome.steps.abundance import assign, chimeras, depth_profile, run_abundance


def test_depth_profile():
    assert depth_profile([(0, 0, 0, "c", 0, 3), (0, 0, 0, "c", 2, 5)], 5) == [1, 1, 2, 1, 1]


def test_assign_and_chimeras():
    alns = {"r1": [(100, 0, 90, "A", 0, 90)],
            "r2": [(100, 0, 40, "A", 0, 40), (100, 40, 90, "B", 0, 50)],
            "r3": [(100, 0, 10, "A", 0, 10)]}
    assert set(assign(alns, 0.5)) == {"r1", "r2"}
    assert chimeras(alns, 0.2) == {("A", "B"): 1}


@pytest.mark.skipif(not shutil.which("minimap2"), reason="minimap2 not on PATH")
def test_abundance_ratio_and_chimeras(tmp_path):
    rng = random.Random(0)
    g0, g1 = rand_seq(rng, 2700), rand_seq(rng, 2400)
    reads, truth = simulate([g0, g0, g0, g1], n_reads=400, host_frac=0, chimera_frac=0.15,
                            copies=(3, 8))
    s = tmp_path / "s1"
    (s / "hostdeplete").mkdir(parents=True)
    (s / "assemble").mkdir()
    write_fastq(s / "hostdeplete" / "reads.nonhost.fastq.gz", reads)
    (s / "assemble" / "contigs.fasta").write_text(f">c0\n{g0}\n>c1\n{g1}\n")
    run_abundance(Config(), s, 2)
    rows = {l.split("\t")[0]: l.split("\t") for l in
            (s / "abundance" / "abundance.tsv").read_text().splitlines()[1:]}
    ratio = float(rows["c0"][3]) / float(rows["c1"][3])
    assert 2.0 < ratio < 4.5                      # truth: 3:1 reads
    assert rows["c0"][6] == "True" and float(rows["c0"][5]) > 0.99
    assert abs(float(rows["c0"][7]) + float(rows["c1"][7]) - 1) < 1e-3
    assert int(rows["c0"][8]) > 0                 # c0/c1 chimeras were injected
