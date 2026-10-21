import math
import random
import shutil

import numpy as np
import pytest
from simulate import mutate, rand_seq

from plantvirome.config import Config
from plantvirome.steps.diversity import alpha, beta, pcoa, run_diversity


def test_alpha_even_two():
    a = alpha(np.array([5.0, 5.0, 0.0]), np.array([10, 10, 0]))
    assert a["observed"] == 2 and a["shannon"] == pytest.approx(math.log(2))
    assert a["simpson"] == pytest.approx(0.5) and a["pielou"] == pytest.approx(1.0)


def test_beta_and_pcoa():
    t = np.array([[1.0, 1, 0], [1, 1, 1], [0, 0, 4]])
    b = beta(t)
    assert b["jaccard"][0, 1] == pytest.approx(1 / 3) and b["jaccard"][0, 2] == 1.0
    assert b["bray_curtis"][0, 0] == 0
    coords, expl = pcoa(b["bray_curtis"])
    assert coords.shape == (3, 2) and 0 < expl.sum() <= 1.0 + 1e-9
    # PCoA distances reproduce the input matrix when all axes are kept (3 samples -> 2 axes)
    full, _ = pcoa(b["bray_curtis"], axes=3)
    assert np.linalg.norm(full[0] - full[1]) == pytest.approx(b["bray_curtis"][0, 1], abs=1e-6)


@pytest.mark.skipif(not shutil.which("minimap2"), reason="minimap2 not on PATH")
def test_cohort(tmp_path):
    rng = random.Random(0)
    a, b, c = (rand_seq(rng, 2700), rand_seq(rng, 2400), rand_seq(rng, 2600))
    a2 = mutate(rng, a, 0.03)                  # same vOTU as `a` at 90% identity
    content = {"s1": {"s1_c1": a, "s1_c2": b}, "s2": {"s2_c1": a2, "s2_c2": b, "s2_c3": c},
               "s3": {"s3_c1": c}}
    for s, cs in content.items():
        (tmp_path / s / "assemble").mkdir(parents=True)
        (tmp_path / s / "abundance").mkdir()
        (tmp_path / s / "assemble" / "contigs.fasta").write_text("".join(f">{k}\n{v}\n" for k, v in cs.items()))
        (tmp_path / s / "abundance" / ".done").touch()
        (tmp_path / s / "abundance" / "abundance.tsv").write_text(
            "contig\tlength\treads\tmean_depth\tdepth_cv\tbreadth\tdetected\trel_abundance\tchimeric_reads\n"
            + "".join(f"{k}\t{len(v)}\t10\t100\t0.1\t1\tTrue\t0.5\t0\n" for k, v in cs.items()))
    run_diversity(Config(), tmp_path, None, 2)
    d = tmp_path / "cohort" / "diversity"
    table = (d / "votu_table.tsv").read_text().splitlines()
    assert len(table) == 1 + 3                  # a/a2, b, c
    alpha_rows = {l.split("\t")[0]: l.split("\t") for l in (d / "alpha.tsv").read_text().splitlines()[1:]}
    assert alpha_rows["s1"][1] == "2" and alpha_rows["s2"][1] == "3" and alpha_rows["s3"][1] == "1"
    jac = [l.split("\t") for l in (d / "beta_jaccard.tsv").read_text().splitlines()]
    assert float(jac[1][2]) == pytest.approx(0.3333, abs=1e-3)   # s1 vs s2
    assert float(jac[1][3]) == 1.0                                # s1 vs s3
    assert (d / "pcoa_bray_curtis.tsv").exists()
