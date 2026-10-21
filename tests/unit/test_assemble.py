import random
import shutil

import pytest
from simulate import mutate, rand_seq, revcomp_seq, simulate, write_fastq

from plantvirome.config import Assemble, Config
from plantvirome.steps.assemble import cluster_units, consensus, read_fasta, run_assemble


def test_consensus_fixes_sub_ins_del():
    seed = "ACGTACGTAC"
    # three members agree: G->T at index 2, C deleted at 5, "GG" inserted after index 6
    out = consensus(seed, [(0, 10, ":2*gt:2-c:1+gg:3")] * 3)
    assert out == "ACTTAGGGTAC"
    # a lone dissenting member does not override the majority
    assert consensus(seed, [(0, 10, ":10")] * 2 + [(0, 10, ":2*gt:7")]) == seed


@pytest.mark.skipif(shutil.which("minimap2") is None, reason="minimap2 missing")
def test_cluster_units_rotation_strand_and_two_genomes(tmp_path):
    rng = random.Random(0)
    g1, g2 = rand_seq(rng, 2700), rand_seq(rng, 2300)
    units = {}
    for i in range(60):
        g = (g1, g2)[i % 3 == 0]
        k = rng.randrange(len(g))
        u = mutate(rng, g[k:] + g[:k], 0.02)
        units[f"u{i}"] = revcomp_seq(u) if i % 2 else u
    contigs = cluster_units(units, Assemble(), 2, tmp_path)
    assert [n for _, n in contigs] == [40, 20]
    # each consensus is the true genome up to rotation/strand
    for (seq, _), g in zip(contigs, (g1, g2)):
        assert abs(len(seq) - len(g)) <= 3
        assert g[200:260] in seq + seq or revcomp_seq(g[200:260]) in seq + seq


@pytest.mark.skipif(not (shutil.which("minimap2") and shutil.which("TideHunter")),
                    reason="minimap2/TideHunter missing")
def test_assemble_recovers_genomes(tmp_path):
    from plantvirome.steps.rca import run_rca
    rng = random.Random(0)
    genomes = [rand_seq(rng, 2700), rand_seq(rng, 2400)]
    reads, _ = simulate(genomes, n_reads=300, host_frac=0, chimera_frac=0, copies=(4, 8))
    hd = tmp_path / "s1" / "hostdeplete"
    hd.mkdir(parents=True)
    write_fastq(hd / "reads.nonhost.fastq.gz", reads)
    cfg = Config()
    run_rca(cfg, tmp_path / "s1", 2)
    run_assemble(cfg, tmp_path / "s1", 2)
    contigs = read_fasta(tmp_path / "s1" / "assemble" / "contigs.fasta")
    assert sorted(len(s) for s in contigs.values()) == pytest.approx([2400, 2700], abs=5)
    rows = (tmp_path / "s1" / "assemble" / "contigs.tsv").read_text().splitlines()[1:]
    assert all(r.split("\t")[3] == "True" for r in rows)


def test_too_few_units(tmp_path):
    rca = tmp_path / "s1" / "rca"
    rca.mkdir(parents=True)
    (rca / "units.fasta").write_text(">a\nACGT\n")
    run_assemble(Config(), tmp_path / "s1", 1)
    assert (tmp_path / "s1" / "assemble" / "contigs.fasta").read_text() == ""


def test_subsample_fastq(tmp_path):
    from plantvirome.steps.assemble import subsample_fastq
    reads = [(f"r{i}", "ACGT" * 250) for i in range(100)]          # 100 kb
    write_fastq(tmp_path / "in.fastq.gz", reads)
    total, kept = subsample_fastq(tmp_path / "in.fastq.gz", tmp_path / "a.fastq.gz", 30_000, 1)
    assert total == 100_000 and 15_000 < kept < 45_000
    _, kept2 = subsample_fastq(tmp_path / "in.fastq.gz", tmp_path / "b.fastq.gz", 30_000, 1)
    assert kept2 == kept                                              # deterministic for a seed
    assert subsample_fastq(tmp_path / "in.fastq.gz", tmp_path / "c.fastq.gz", None, 1) == (100_000, 100_000)


@pytest.mark.skipif(not all(shutil.which(t) for t in ("flye", "TideHunter")), reason="flye/TideHunter missing")
def test_flye_method_collapses_multimers(tmp_path):
    rng = random.Random(5)
    g1, g2 = rand_seq(rng, 2700), rand_seq(rng, 2400)
    reads, _ = simulate([g1, g1, g2], n_reads=150, host_frac=0, chimera_frac=0, copies=(3, 8))
    hd = tmp_path / "s1" / "hostdeplete"
    hd.mkdir(parents=True)
    write_fastq(hd / "reads.nonhost.fastq.gz", reads)
    cfg = Config.model_validate({"assemble": {"method": "flye"}})
    run_assemble(cfg, tmp_path / "s1", 2)
    contigs = read_fasta(tmp_path / "s1" / "assemble" / "contigs.fasta")
    assert sorted(len(s) for s in contigs.values()) == pytest.approx([2400, 2700], abs=30)
    rows = (tmp_path / "s1" / "assemble" / "contigs.tsv").read_text().splitlines()[1:]
    assert all(r.split("\t")[3] == "True" for r in rows)
