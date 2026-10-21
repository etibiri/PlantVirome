import random
import shutil

import pytest
from simulate import mutate, rand_seq, write_fastq

from plantvirome.config import Config
from plantvirome.steps.variety import CS, depth_array, run_variety


def test_cs_regex_and_depth():
    assert [m.group(0) for m in CS.finditer(":5*ag+tt-cc:3")] == [":5", "*ag", "+tt", "-cc", ":3"]
    assert depth_array({1: 1, 4: -1}, 5) == [0, 1, 1, 1, 0]


def test_no_reference_skips(tmp_path):
    (tmp_path / "s1").mkdir()
    run_variety(Config(), tmp_path / "s1", 1)
    assert (tmp_path / "s1" / "variety" / ".done").exists()


@pytest.mark.skipif(not shutil.which("minimap2"), reason="minimap2 not on PATH")
def test_variants_and_partial_coverage(tmp_path):
    rng = random.Random(0)
    ref = rand_seq(rng, 5000)
    snv = {500: "", 1200: "", 2000: "", 3000: "", 4500: ""}      # 4500 is outside covered region
    sample = list(ref)
    for p in snv:
        sample[p] = next(b for b in "ACGT" if b != ref[p])
    sample = "".join(sample)
    reads = []
    for i in range(160):                                           # ~60x over positions 0-4000 only
        s = rng.randrange(0, 2500)
        reads.append((f"r{i}", mutate(rng, sample[s:s + 1500], 0.03)))
    (tmp_path / "ref.fa").write_text(f">chloro\n{ref}\n")
    (tmp_path / "s1" / "qc").mkdir(parents=True)
    write_fastq(tmp_path / "s1" / "qc" / "reads.qc.fastq.gz", reads)
    cfg = Config.model_validate({"host": {"variety_references": [str(tmp_path / "ref.fa")]}})
    run_variety(cfg, tmp_path / "s1", 2)
    out = tmp_path / "s1" / "variety"
    called = {int(l.split("\t")[1]) - 1 for l in (out / "variants.tsv").read_text().splitlines()[1:]}
    assert called == {500, 1200, 2000, 3000}                       # no false positives, 4500 unseen
    cov = (out / "coverage.tsv").read_text().splitlines()[1].split("\t")
    assert cov[-1] == "PARTIAL" and float(cov[4]) < 0.9
    cons = (out / "consensus.fasta").read_text().splitlines()[1]
    assert cons[4500] == "N" and cons[500] == sample[500]
