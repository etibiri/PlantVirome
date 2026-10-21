"""End to end on simulated data: reads -> submission files, two samples + cohort diversity."""
import random
import shutil

import pytest
from simulate import (make_genome, make_tiny_db, mutate, rand_protein, rand_seq, simulate,
                      write_fastq)
from typer.testing import CliRunner

from plantvirome.cli import app

TOOLS = ["porechop_abi", "chopper", "minimap2", "TideHunter", "diamond", "taxonkit"]


@pytest.mark.needs_tools
@pytest.mark.skipif(any(shutil.which(t) is None for t in TOOLS), reason="external tools missing")
def test_full_pipeline(tmp_path):
    rng = random.Random(7)
    p1, p2 = rand_protein(rng), rand_protein(rng)
    g1, g2 = make_genome(rng, p1), make_genome(rng, p2)
    host = rand_seq(rng, 20000)
    db, taxdump = make_tiny_db(tmp_path, {11: p1, 12: p2})
    (tmp_path / "host.fa").write_text(f">host\n{host}\n")

    for seed, (sample, genomes) in enumerate({"s1": [g1], "s2": [g1, g2]}.items()):
        reads, _ = simulate(genomes, n_reads=350, host_frac=0.2, chimera_frac=0.1, host=host,
                            seed=seed, copies=(3, 8))
        (tmp_path / "reads" / sample).mkdir(parents=True)
        write_fastq(tmp_path / "reads" / sample / "r.fastq.gz", reads)

    cfg = tmp_path / "config.yaml"
    cfg.write_text(f"""
threads: 2
host: {{name: test, references: [{tmp_path / 'host.fa'}], variety_references: [{tmp_path / 'host.fa'}]}}
databases: {{viral_protein: {db}, taxdump: {taxdump}}}
qc: {{min_length: 500}}
rca: {{min_period: 600, max_period: 3500}}
genbank: {{origin_motif: TAATATTAC, origin_offset: 7, metadata: {{host: Manihot esculenta}}}}
""")
    out = tmp_path / "results"
    runner = CliRunner()
    for sample in ("s1", "s2"):
        res = runner.invoke(app, ["run", "--config", str(cfg), "--outdir", str(out), "--sample", sample,
                                  "--reads", str(tmp_path / "reads" / sample)])
        assert res.exit_code == 0, res.output + str(res.exception)
    res = runner.invoke(app, ["diversity", "--config", str(cfg), "--outdir", str(out)])
    assert res.exit_code == 0, res.output + str(res.exception)

    fsa = {s: (out / s / "genbank" / "submission.fsa").read_text() for s in ("s1", "s2")}
    assert fsa["s1"].count(">") == 1 and fsa["s2"].count(">") == 2
    assert "[organism=GenA sp.]" in fsa["s1"] and "[organism=GenB sp.]" in fsa["s2"]
    assert all(rec.splitlines()[1].startswith("AC") for s in fsa.values() for rec in s.split(">")[1:])
    div = out / "cohort" / "diversity"
    alpha = {l.split("\t")[0]: l.split("\t")[1] for l in (div / "alpha.tsv").read_text().splitlines()[1:]}
    assert alpha == {"s1": "1", "s2": "2"}
    assert (out / "s1" / "variety" / "coverage.tsv").exists()
