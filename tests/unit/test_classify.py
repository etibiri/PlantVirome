import itertools
import random
import shutil

import pytest

from plantvirome.config import Config
from plantvirome.steps.classify import summarize, tier

# Standard genetic code (TCAG order); first codon seen per amino acid is used.
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {}
for aa, cod in zip(_AAS, ("".join(c) for c in itertools.product("TCAG", repeat=3))):
    if aa != "*":
        CODON.setdefault(aa, cod)
AA = list(CODON)


def test_tier():
    t = [(85.0, "close"), (60.0, "related"), (35.0, "distant")]
    assert [tier(p, t) for p in (95, 70, 40, 20)] == ["close", "related", "distant", "unclassified"]


def test_summarize_lca_inputs():
    h = lambda b, ids, s: {"pident": 90, "bits": b, "evalue": 1e-9, "qlen": 100, "span": s, "taxids": ids}
    s = summarize([h(100, {"1"}, (1, 50)), h(95, {"2"}, (1, 50)), h(10, {"3"}, (51, 100))], 10)
    assert s["taxids"] == {"1", "2"} and s["cover"] == 1.0


@pytest.mark.skipif(not (shutil.which("diamond") and shutil.which("taxonkit")),
                    reason="diamond/taxonkit missing")
def test_classify_tiny_db(tmp_path):
    import subprocess

    from plantvirome.steps.classify import run_classify
    rng = random.Random(0)
    prot = {"P1": "".join(rng.choice(AA) for _ in range(400)),
            "P2": "".join(rng.choice(AA) for _ in range(400))}
    mutate = lambda p, r: "".join(rng.choice(AA) if rng.random() < r else a for a in p)
    dna = lambda p: "".join(CODON[a] for a in p)
    tax = tmp_path / "tax"
    tax.mkdir()
    rows = [(1, 1, "no rank"), (2, 1, "superkingdom"), (10, 2, "family"), (11, 10, "genus"), (12, 10, "genus")]
    (tax / "nodes.dmp").write_text("".join(f"{a}\t|\t{b}\t|\t{r}\t|\t\t|\n" for a, b, r in rows))
    names = {1: "root", 2: "Viruses", 10: "FamA", 11: "GenA", 12: "GenB"}
    (tax / "names.dmp").write_text("".join(f"{k}\t|\t{v}\t|\t\t|\tscientific name\t|\n" for k, v in names.items()))
    for f in ("merged.dmp", "delnodes.dmp"):
        (tax / f).touch()
    (tmp_path / "db.faa").write_text("".join(f">{k}\n{v}\n" for k, v in prot.items()))
    (tmp_path / "map.tsv").write_text("accession\taccession.version\ttaxid\tgi\nP1\tP1\t11\t0\nP2\tP2\t12\t0\n")
    subprocess.run(["diamond", "makedb", "--in", tmp_path / "db.faa", "-d", tmp_path / "db",
                    "--taxonmap", tmp_path / "map.tsv", "--taxonnodes", tax / "nodes.dmp",
                    "--taxonnames", tax / "names.dmp"], check=True)
    asm = tmp_path / "s1" / "assemble"
    asm.mkdir(parents=True)
    (asm / "contigs.fasta").write_text(
        f">exact\n{dna(prot['P1'])}\n>diverged\n{dna(mutate(prot['P2'], 0.2))}\n"
        f">random\n{''.join(rng.choice('ACGT') for _ in range(1200))}\n")
    cfg = Config.model_validate({"databases": {"viral_protein": tmp_path / "db.dmnd", "taxdump": tax}})
    run_classify(cfg, tmp_path / "s1", 2)
    rows = {l.split("\t")[0]: l.split("\t") for l in
            (tmp_path / "s1" / "classify" / "taxonomy.tsv").read_text().splitlines()[1:]}
    assert rows["exact"][2:4] == ["11", "GenA"] and rows["exact"][-1] == "close"
    assert rows["diverged"][3] == "GenB" and rows["diverged"][-1] == "related"
    assert rows["random"][-1] == "unclassified"
