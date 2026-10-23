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


def test_origin_with_mismatch():
    # non-palindromic motif; a random 9-mer has a 1-mismatch decoy about half the time
    # (-> ambiguous), so look for a genome where the damaged motif is the only hit
    motif = "ACGGTCATC"
    for seed in range(30):
        rng = random.Random(seed)
        g = rand_seq(rng, 2700)
        g = g[:1000] + motif + g[1009:]
        bad = g[:1000] + "ACGGTCATG" + g[1009:]
        if orient_and_rotate(bad, motif, 4)[1] != "origin_not_found":
            continue
        seq, note = orient_and_rotate(bad, motif, 4, max_mismatch=1)
        if note == "origin_approx":
            assert seq.startswith("TCATG")
            assert orient_and_rotate(g, motif, 4, max_mismatch=1)[1] == "origin_ok"   # exact hit wins
            return
        assert note == "origin_ambiguous"
    raise AssertionError("no seed gave a unique 1-mismatch hit")


def test_palindromic_motif_is_ambiguous_after_one_error():
    """TAATATTA is its own reverse complement, so a damaged begomovirus nonanucleotide matches
    on both strands at the same site: use orient_by_reference instead of origin_max_mismatch."""
    g = genome(random.Random(1), orf=False)
    bad = g.replace(MOTIF, "TAATATTAG")
    assert orient_and_rotate(bad, MOTIF, 7, max_mismatch=1)[1] == "origin_ambiguous"


def test_run_genbank_orients_by_reference(tmp_path):
    """Reverse-complemented contig with a damaged (palindromic) motif: the reference strand
    plus a 1-mismatch search on that strand only gives a unique, correctly oriented origin."""
    rng = random.Random(3)
    good = genome(rng)
    true_start = good[1007:] + good[:1007]                      # starts at "AC..."
    damaged = good.replace(MOTIF, "TAATATTAG")
    contig = revcomp(damaged)
    s = tmp_path / "s1"
    for d in ("assemble", "classify", "abundance", "species"):
        (s / d).mkdir(parents=True)
    (s / "assemble" / "contigs.fasta").write_text(f">s1_c1\n{contig}\n")
    (s / "assemble" / "contigs.tsv").write_text(
        "contig\tlength\tcoverage\tcircular\tcirc_support\ns1_c1\t2700\t90\tTrue\t5\n")
    (s / "abundance" / "abundance.tsv").write_text(
        "contig\tlength\treads\tmean_depth\tdepth_cv\tbreadth\tdetected\trel_abundance\tchimeric_reads\n"
        "s1_c1\t1\t50\t100\t0.1\t1\tTrue\t1\t0\n")
    (s / "species" / "species.tsv").write_text(
        "contig\tlength\tbest_ref\tspecies\ttaxid\tcomponent\tidentity\tcover\tcall\tcaveat\tref_strand\tref_origin\n"
        "s1_c1\t2700\tR1\tX virus\t1\tA\t0.95\t1.00\tknown_species\t\t-\t0\n")
    cfg = Config.model_validate({"genbank": {"origin_motif": MOTIF, "origin_offset": 7, "origin_max_mismatch": 1,
                                             "orient_by_reference": True}})
    run_genbank(cfg, s)
    seq = (s / "genbank" / "submission.fsa").read_text().splitlines()[1]
    assert seq.startswith("AG") and seq[2:50] == true_start[2:50]   # damaged base kept, origin right
    assert "origin_approx" in (s / "genbank" / "contigs.qc.tsv").read_text()


def test_organism_for_other_component_uses_nearest_reference():
    tax = {"name": "Begomovirus manihotisburkinafasoense", "rank": "species"}   # protein-level guess
    sp = {"call": "other_component", "species": "Begomovirus manihotis"}
    assert organism(tax, sp, "Begomovirus") == ("Begomovirus manihotis", "organism_from_nearest_component_reference")


def test_checkv_reasons_never_block_on_completeness():
    from plantvirome.steps.genbank import checkv_reasons
    ok = {"contamination": "0.0", "host_genes": "0", "viral_genes": "2", "completeness": "50.1"}
    assert checkv_reasons(ok, 10) == []
    assert checkv_reasons({**ok, "contamination": "NA"}, 10) == []
    assert checkv_reasons({**ok, "contamination": "25"}, 10) == ["checkv_contamination"]
    assert checkv_reasons({**ok, "host_genes": "3"}, 10) == ["checkv_host_genes"]
