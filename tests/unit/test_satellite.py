import random
import shutil

import pytest
from simulate import rand_seq

from plantvirome.config import Config
from plantvirome.steps.diversity import identified_contigs
from plantvirome.steps.genbank import run_genbank
from plantvirome.steps.satellite import protein_hints, run_satellite, satellite_call


def diverge(rng, seq, rate):
    out = [rng.choice([x for x in "ACGT" if x != b]) if rng.random() < rate else b for b in seq]
    k = rng.randrange(len(seq))
    return "".join(out[k:] + out[:k])


def test_satellite_call_mapping():
    assert satellite_call("known_species", True) == "known_satellite"
    assert satellite_call("novel_species_candidate", True) == "novel_satellite_candidate"
    assert satellite_call("known_species", False) == "satellite_related"


def make_sample(tmp_path, contigs, taxonomy=""):
    d = tmp_path / "s1"
    (d / "assemble").mkdir(parents=True)
    (d / "classify").mkdir()
    (d / "assemble" / "contigs.fasta").write_text("".join(f">{k}\n{v}\n" for k, v in contigs.items()))
    (d / "classify" / "taxonomy.tsv").write_text(
        "contig\tlength\ttaxid\tname\trank\tlineage\tbest_pident\tquery_cover\tevalue\ttier\n" + taxonomy)
    return d


@pytest.mark.skipif(shutil.which("minimap2") is None, reason="minimap2 missing")
def test_run_satellite(tmp_path):
    rng = random.Random(1)
    b1, b2, a1 = rand_seq(rng, 1350), rand_seq(rng, 1340), rand_seq(rng, 1380)
    (tmp_path / "beta.fa").write_text(f">B1.1|Beta_one|1|satellite\n{b1}\n>B2.1|Beta_two|2|satellite\n{b2}\n")
    (tmp_path / "alpha.fa").write_text(f">A1.1|Alpha_one|3|satellite\n{a1}\n")
    contigs = {"s1_c1": diverge(rng, b1, 0.03),        # beta, same species
               "s1_c2": diverge(rng, b2, 0.12),        # ~88 %: below the 91 % threshold
               "s1_c3": diverge(rng, a1, 0.05),        # alpha group has no threshold
               "s1_c4": rand_seq(rng, 1300),           # nothing
               "s1_c5": diverge(rng, b1, 0.02) + rand_seq(rng, 1400),   # 2.7 kb: outside the window
               "s1_c6": rand_seq(rng, 1250)}           # protein evidence only
    d = make_sample(tmp_path, contigs,
                    "s1_c6\t1250\t9\tAlphasatellitidae sp.\tfamily\tViruses;Alphasatellitidae\t45\t0.5\t1e-9\tdistant\n")
    cfg = Config.model_validate({"satellite": {"groups": {
        "betasatellite": {"reference": str(tmp_path / "beta.fa"), "species_threshold": 0.91},
        "alphasatellite": {"reference": str(tmp_path / "alpha.fa")}}}})
    run_satellite(cfg, d, 2)
    rows = {l.split("\t")[0]: l.split("\t") for l in
            (d / "satellite" / "satellites.tsv").read_text().splitlines()[1:]}
    assert rows["s1_c1"][2:5] == ["betasatellite", "B1.1", "Beta one"] and rows["s1_c1"][8] == "known_satellite"
    assert rows["s1_c2"][4] == "Beta two" and rows["s1_c2"][8] == "novel_satellite_candidate"
    assert rows["s1_c3"][2] == "alphasatellite" and rows["s1_c3"][8] == "satellite_related"
    assert rows["s1_c6"][8] == "protein_only_candidate" and "Alphasatellitidae" in rows["s1_c6"][14]
    assert "s1_c4" not in rows and "s1_c5" not in rows
    # c6 is identified through its protein tier (classify), not through the satellite file
    assert identified_contigs(d) == {"s1_c1", "s1_c2", "s1_c3", "s1_c6"}


def test_skipped_without_groups(tmp_path):
    (tmp_path / "s1").mkdir()
    run_satellite(Config(), tmp_path / "s1", 1)
    assert (tmp_path / "s1" / "satellite" / ".done").exists()


def test_protein_hints(tmp_path):
    f = tmp_path / "t.tsv"
    f.write_text("contig\tlength\ttaxid\tname\trank\tlineage\tbest_pident\tquery_cover\tevalue\ttier\n"
                 "a\t1\t1\tBegomovirus\tgenus\tViruses;Geminiviridae;Begomovirus\t90\t1\t0\tclose\n"
                 "b\t1\t2\tAlphasatellitidae\tfamily\tViruses;Alphasatellitidae\t50\t1\t0\tdistant\n")
    assert protein_hints(f) == {"b": "Alphasatellitidae"}


def test_genbank_names_known_satellite(tmp_path):
    """A satellite is named by the satellite step; orientation comes from its reference."""
    rng = random.Random(3)
    seq = rand_seq(rng, 1350)
    d = tmp_path / "s1"
    for sub in ("assemble", "classify", "abundance", "species", "satellite"):
        (d / sub).mkdir(parents=True)
    (d / "assemble" / "contigs.fasta").write_text(f">s1_c1\n{seq}\n")
    (d / "assemble" / "contigs.tsv").write_text("contig\tlength\tcoverage\tcircular\tcirc_support\ns1_c1\t1350\t90\tTrue\t5\n")
    (d / "abundance" / "abundance.tsv").write_text(
        "contig\tlength\treads\tmean_depth\tdepth_cv\tbreadth\tdetected\trel_abundance\tchimeric_reads\n"
        "s1_c1\t1\t50\t100\t0.1\t1\tTrue\t1\t0\n")
    (d / "species" / "species.tsv").write_text(
        "contig\tlength\tbest_ref\tspecies\ttaxid\tcomponent\tidentity\tcover\tcall\tcaveat\tref_strand\tref_origin\tref_length\tref_cover\n"
        "s1_c1\t1350\t\t\t\t\t0.000\t0.00\tno_reference_match\t\t\t\t\t\n")
    (d / "satellite" / "satellites.tsv").write_text(
        "contig\tlength\tgroup\tbest_ref\tspecies\ttaxid\tidentity\tcover\tcall\tcaveat\tref_strand\tref_origin\tref_length\tref_cover\tprotein_hint\n"
        "s1_c1\t1350\tbetasatellite\tB1.1\tCotton leaf curl Multan betasatellite\t1\t0.97\t1.00\tknown_satellite\t\t+\t100\t1350\t1.00\t\n")
    cfg = Config.model_validate({"genbank": {"min_length": 1000, "min_depth": 20,
                                             "origin_motif": "TAATATTAC", "orient_by_reference": True}})
    run_genbank(cfg, d)
    row = (d / "genbank" / "contigs.qc.tsv").read_text().splitlines()[1].split("\t")
    assert row[8] == "submit" and row[7] == "Cotton leaf curl Multan betasatellite"
    fsa = (d / "genbank" / "submission.fsa").read_text().splitlines()
    assert fsa[1] == seq[100:] + seq[:100]                      # starts where the reference starts
