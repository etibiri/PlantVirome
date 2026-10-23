import random
import shutil

import pytest
from simulate import rand_seq

from plantvirome.steps.annotate import CODON, all_orfs, annotate, classify_genes, gene_reasons, merge_unnamed
from plantvirome.steps.download_refs import genbank_cds
from plantvirome.steps.genbank import revcomp

CODONS = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT"
          if a + b + c not in ("TAA", "TAG", "TGA", "ATG")]
GB = """VERSION     AB000001.1
FEATURES             Location/Qualifiers
     source          1..100
                     /organism="x"
     CDS             1..30
                     /gene="AC1"
                     /product="replication-associated
                     protein"
                     /translation="MKLV
                     QQ"
     CDS             complement(40..60)
                     /product="coat protein"
                     /translation="MA"
ORIGIN
        1 acgt
//
"""


def orf(rng, n):
    return "ATG" + "".join(rng.choice(CODONS) for _ in range(n)) + "TAA"


def test_genbank_cds_parses_multiline_and_missing_gene():
    assert genbank_cds(GB) == [("AB000001.1", "AC1", "replication-associated protein", "MKLVQQ"),
                               ("AB000001.1", "", "coat protein", "MA")]


def test_all_orfs_both_strands_with_coordinates():
    rng = random.Random(3)
    o = orf(rng, 100)
    g = rand_seq(rng, 300)
    g = g[:50] + "TAA" + o + g[50 + 3 + len(o):]       # stop before the ORF so the scan starts at its ATG
    for seq, strand in ((g, 1), (revcomp(g), -1)):
        hit = [x for x in all_orfs(seq, 60) if x["strand"] == strand and len(x["protein"]) == 101]
        assert hit and hit[0]["end"] - hit[0]["begin"] + 1 == len(o)
        assert seq[hit[0]["begin"] - 1:hit[0]["end"]] in (o, revcomp(o))


def test_classify_genes_states_and_duplicates():
    prots = [{"id": "a", "gene": "AV1", "key": "coat protein", "aa": 250},
             {"id": "b", "gene": "AV1", "key": "coat protein", "aa": 250},      # duplicate product
             {"id": "c", "gene": "AC1", "key": "rep", "aa": 358},
             {"id": "d", "gene": "AC2", "key": "trap", "aa": 134},
             {"id": "e", "gene": "AC5", "key": "ac5", "aa": 60}]               # too short to be required
    best = {"b": (90.0, 1.0), "c": (80.0, 0.4), "d": (20.0, 1.0)}
    r = classify_genes(prots, best, 100, 35.0)
    assert r == {"found": ["AV1"], "truncated": ["AC1"], "missing": ["AC2"]}
    assert gene_reasons({**r}) == ["missing_gene:AC2", "truncated_gene:AC1"]
    core = classify_genes(prots, best, 100, 35.0, core={"rep"})              # only the core products count
    assert core == {"found": [], "truncated": ["AC1"], "missing": []}


def test_merge_unnamed_keeps_only_non_overlapping_gene_finder_orfs():
    named = [{"begin": 100, "end": 400, "product": "x"}]
    found = [{"begin": 150, "end": 300}, {"begin": 500, "end": 700}]
    assert [(o["begin"], o["index"]) for o in merge_unnamed(named, found)] == [(100, 1), (500, 2)]


@pytest.mark.needs_tools
@pytest.mark.skipif(not shutil.which("diamond"), reason="diamond not on PATH")
def test_annotate_names_orf_and_flags_missing_gene(tmp_path):
    rng = random.Random(4)
    cp = "".join(rng.choice("ACDEFGHIKLNPQRSTVWY") for _ in range(149))
    rep = "".join(rng.choice("ACDEFGHIKLNPQRSTVWY") for _ in range(149))
    (tmp_path / "refs.faa").write_text(f">R1|AV1|coat_protein\nM{cp}\n>R1|AC1|replication_protein\nM{rep}\n")
    rev = {}
    for c in CODONS:
        rev.setdefault(CODON[c], []).append(c)
    dna = "ATG" + "".join(rng.choice(rev[a]) for a in cp) + "TAA"     
    res = annotate({"c1": rand_seq(rng, 40) + "TAA" + dna + rand_seq(rng, 40)}, tmp_path / "refs.faa",
                   tmp_path / "w", 1, 60)["c1"]
    assert [o["product"] for o in res["orfs"]] == ["coat protein"]
    assert res["found"] == ["AV1"] and res["missing"] == ["AC1"] and res["ref"] == "R1"
