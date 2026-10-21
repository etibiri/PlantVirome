import random

from simulate import rand_seq

from plantvirome.config import Config
from plantvirome.steps.genbank import (feature_table, find_orfs, orient_and_rotate, organism,
                                       revcomp, run_genbank)

MOTIF = "TAATATTAC"
CODONS = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT"
          if a + b + c not in ("TAA", "TAG", "TGA")]


def genome(rng, n=2700, orf=True):
    """Random genome with the origin motif at position 1000 and one long ORF."""
    g = rand_seq(rng, n)
    if orf:
        o = "ATG" + "".join(rng.choice(CODONS) for _ in range(300)) + "TAA"
        g = g[:1500] + o + g[1500 + len(o):]
    return g[:1000] + MOTIF + g[1000 + len(MOTIF):]


def test_rotate_both_strands_and_circular_junction():
    rng = random.Random(1)
    g = genome(rng, orf=False)
    want, note = orient_and_rotate(g, MOTIF, 7)
    assert note == "origin_ok" and want.startswith("AC") and len(want) == len(g)
    assert orient_and_rotate(revcomp(g), MOTIF, 7)[0] == want       # opposite strand
    k = 2000                                                        # motif split over the origin
    assert orient_and_rotate(g[k:] + g[:k], MOTIF, 7)[0] == want
    assert orient_and_rotate(g.replace(MOTIF, "A" * 9), MOTIF, 7)[1] == "origin_not_found"
    assert orient_and_rotate(g + MOTIF, MOTIF, 7)[1] == "origin_ambiguous"


def test_orfs_and_table():
    rng = random.Random(2)
    g = genome(rng)
    for seq, strand in ((g, 1), (revcomp(g), -1)):
        orfs = [o for o in find_orfs(seq, 60) if o["strand"] == strand]
        assert orfs and max(len(o["protein"]) for o in orfs) >= 300
        for i, o in enumerate(orfs, 1):
            o["index"] = i
        lines = feature_table("s", orfs).splitlines()
        a, b = map(lambda x: int(x.strip("<>")), lines[1].split("\t")[:2])
        assert (a < b) == (strand == 1)


def test_organism():
    assert organism({"name": "Cassava mosaic virus", "rank": "species"}) == ("Cassava mosaic virus", None)
    assert organism({"name": "Begomovirus", "rank": "genus"})[0] == "Begomovirus sp."
    assert organism({})[0] == "uncultured virus"


def test_run_genbank_status(tmp_path):
    rng = random.Random(3)
    good, junk = genome(rng), rand_seq(rng, 150)
    s = tmp_path / "s1"
    for d in ("assemble", "classify", "abundance"):
        (s / d).mkdir(parents=True)
    (s / "assemble" / "contigs.fasta").write_text(
        f">s1_c1\n{good[700:] + good[:700]}\n>s1_c2\n{junk}\n>s1_c3\n{rand_seq(rng, 2000)}\n")
    (s / "assemble" / "contigs.tsv").write_text(
        "contig\tlength\tcoverage\tcircular\tcirc_support\n"
        "s1_c1\t2700\t90\tTrue\t5\ns1_c2\t150\t90\tTrue\t5\ns1_c3\t2000\t90\tFalse\t0\n")
    (s / "classify" / "taxonomy.tsv").write_text(
        "contig\tlength\ttaxid\tname\trank\tlineage\tbest_pident\tquery_cover\tevalue\ttier\n"
        "s1_c1\t2700\t11\tCassava mosaic virus\tspecies\tx\t95\t0.9\t1e-50\tclose\n")
    (s / "abundance" / "abundance.tsv").write_text(
        "contig\tlength\treads\tmean_depth\tdepth_cv\tbreadth\tdetected\trel_abundance\tchimeric_reads\n"
        + "".join(f"s1_c{i}\t1\t50\t100\t0.1\t1\tTrue\t0.3\t0\n" for i in (1, 2, 3)))
    cfg = Config.model_validate({"genbank": {"origin_motif": MOTIF, "origin_offset": 7,
                                             "metadata": {"host": "Manihot esculenta"}}})
    run_genbank(cfg, s)
    rows = {l.split("\t")[0]: l.split("\t") for l in
            (s / "genbank" / "contigs.qc.tsv").read_text().splitlines()[1:]}
    assert rows["s1_c1"][8] == "submit" and int(rows["s1_c1"][9]) >= 1
    assert "too_short" in rows["s1_c2"][10] and rows["s1_c2"][8] == "reject"
    assert rows["s1_c3"][8] == "review"
    fsa = (s / "genbank" / "submission.fsa").read_text()
    assert fsa.count(">") == 1 and "[organism=Cassava mosaic virus]" in fsa
    assert "[topology=circular]" in fsa and "[host=Manihot esculenta]" in fsa
    assert fsa.splitlines()[1].startswith("AC")


def test_organism_prefers_species_call():
    tax = {"name": "Begomovirus", "rank": "genus"}
    assert organism(tax, {"call": "known_species", "species": "Real virus"}, "Begomovirus") == ("Real virus", None)
    assert organism(tax, {"call": "novel_species_candidate"}, "Begomovirus") == \
        ("Begomovirus sp.", "novel_species_candidate")
    assert organism(tax, {"call": "no_reference_match"}, "Begomovirus")[0] == "Begomovirus sp."
