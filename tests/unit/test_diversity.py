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


@pytest.mark.skipif(not shutil.which("minimap2"), reason="minimap2 not on PATH")
def test_blank_flags_contamination(tmp_path):
    """The same genome in a sample and in the blank: flagged when the sample is not >= 10x deeper."""
    rng = random.Random(1)
    a, b = rand_seq(rng, 2700), rand_seq(rng, 2400)
    content = {"s1": ({"s1_c1": a, "s1_c2": b}, {"s1_c1": 500, "s1_c2": 50}),    # blank depth is 10: a 50x -> ok, b 5x -> flagged
               "s2": ({"s2_c1": a}, {"s2_c1": 60}),                                # 6x: flagged
               "blank": ({"blank_c1": a, "blank_c2": b}, {"blank_c1": 10, "blank_c2": 10})}
    for s, (cs, dp) in content.items():
        (tmp_path / s / "assemble").mkdir(parents=True)
        (tmp_path / s / "abundance").mkdir()
        (tmp_path / s / "assemble" / "contigs.fasta").write_text("".join(f">{k}\n{v}\n" for k, v in cs.items()))
        (tmp_path / s / "abundance" / ".done").touch()
        (tmp_path / s / "abundance" / "abundance.tsv").write_text(
            "contig\tlength\treads\tmean_depth\tdepth_cv\tbreadth\tdetected\trel_abundance\tchimeric_reads\n"
            + "".join(f"{k}\t{len(v)}\t10\t{dp[k]}\t0.1\t1\tTrue\t0.5\t0\n" for k, v in cs.items()))
    cfg = Config.model_validate({"diversity": {"blanks": ["blank"], "blank_ratio": 10}})
    run_diversity(cfg, tmp_path, None, 2)
    d = tmp_path / "cohort" / "diversity"
    flags = {l.split("\t")[0]: l.split("\t")[6].strip() for l in (d / "contamination.tsv").read_text().splitlines()[1:]}
    assert flags == {"s1_c1": "ok", "s1_c2": "possible_contamination", "s2_c1": "possible_contamination"}
    alpha_rows = [l.split("\t")[0] for l in (d / "alpha.tsv").read_text().splitlines()[1:]]
    assert alpha_rows == ["s1", "s2"]                 # the blank is not analysed as a sample
    assert "blank" in (d / "votu_table.tsv").read_text().splitlines()[0]


@pytest.mark.skipif(not shutil.which("minimap2"), reason="minimap2 not on PATH")
def test_fragments_do_not_bridge_genomes(tmp_path):
    """A small piece shared by two genomes (a common region) must not merge them into one vOTU."""
    rng = random.Random(2)
    shared = rand_seq(rng, 250)
    a, b = rand_seq(rng, 1200) + shared + rand_seq(rng, 1250), rand_seq(rng, 1100) + shared + rand_seq(rng, 1350)
    content = {"s1": {"s1_c1": a, "s1_c2": b, "s1_c3": shared}}
    d = tmp_path / "s1"
    (d / "assemble").mkdir(parents=True)
    (d / "abundance").mkdir()
    (d / "assemble" / "contigs.fasta").write_text("".join(f">{k}\n{v}\n" for k, v in content["s1"].items()))
    (d / "abundance" / ".done").touch()
    (d / "abundance" / "abundance.tsv").write_text(
        "contig\tlength\treads\tmean_depth\tdepth_cv\tbreadth\tdetected\trel_abundance\tchimeric_reads\n"
        + "".join(f"{k}\t{len(v)}\t10\t100\t0.1\t1\tTrue\t0.5\t0\n" for k, v in content["s1"].items()))
    run_diversity(Config(), tmp_path, None, 2)
    votus = (tmp_path / "cohort" / "diversity" / "votus.tsv").read_text()
    assert "s1_c3" not in votus                                   # 250 bp: below min_contig_length
    assert len(set(l.split("\t")[1] for l in votus.splitlines()[1:])) == 2   # A and B stay apart
