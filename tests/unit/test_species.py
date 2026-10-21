import random
import shutil

import pytest
from simulate import rand_seq

from plantvirome.config import Config, Species
from plantvirome.steps.download_refs import component, fasta_records
from plantvirome.steps.species import call, run_species

REFS = {"R1": {"acc": "R1", "species": "Alpha virus", "taxid": "1", "component": "A"},
        "R2": {"acc": "R2", "species": "Alpha virus", "taxid": "1", "component": "B"}}
SP = Species(species_threshold=0.91, strain_threshold=0.94)


def test_component():
    assert component("X virus isolate Y segment DNA-A, complete sequence") == "A"
    assert component("X virus segment DNA-B, complete sequence") == "B"
    assert component("X virus complete genome") == "A"
    assert component("Cotton leaf curl Multan betasatellite, complete sequence") == "satellite"


def test_fasta_records():
    assert fasta_records(">A.1 t\nAC\nGT\n>B.1 u\nTT\n") == {"A.1": "ACGT", "B.1": "TT"}


def test_call_logic():
    c = lambda hits: call(hits, REFS, SP)
    assert c({}) ["call"] == "no_reference_match"
    assert c({"R1": (0.96, 1.0)})["call"] == "known_species" and c({"R1": (0.96, 1.0)})["caveat"] == ""
    assert c({"R1": (0.925, 1.0)})["caveat"] == "new_strain"
    assert c({"R1": (0.912, 1.0)})["caveat"] == "borderline,new_strain"
    assert c({"R1": (0.88, 0.95)})["call"] == "novel_species_candidate"
    partial = c({"R1": (0.93, 0.5)})            # high identity over half the contig is not enough
    assert partial["call"] == "novel_species_candidate" and partial["caveat"] == "partial_alignment"
    assert c({"R1": (0.9, 0.1)})["call"] == "no_reference_match"
    assert c({"R2": (0.95, 1.0)})["call"] == "other_component"
    # a better full-length DNA-B hit must not beat nothing; best identity wins among full hits
    assert c({"R1": (0.80, 1.0), "R2": (0.90, 1.0)})["ref"] == "R2"


def diverge(rng, seq, rate):
    """Substitutions at exactly the given rate (identity = 1 - rate), then rotate."""
    out = [rng.choice([x for x in "ACGT" if x != b]) if rng.random() < rate else b for b in seq]
    k = rng.randrange(len(seq))
    return "".join(out[k:] + out[:k])


@pytest.mark.skipif(shutil.which("minimap2") is None, reason="minimap2 missing")
def test_run_species(tmp_path):
    rng = random.Random(0)
    a, b, dnab = rand_seq(rng, 2750), rand_seq(rng, 2700), rand_seq(rng, 2650)
    (tmp_path / "refs.fasta").write_text(
        f">A1.1|Alpha_virus|1|A\n{a}\n>B1.1|Beta_virus|2|A\n{b}\n>A1B.1|Alpha_virus|1|B\n{dnab}\n")
    contigs = {"s1_c1": diverge(rng, a, 0.03),          # same species as Alpha
               "s1_c2": diverge(rng, b, 0.12),          # ~88 %: novel candidate near Beta
               "s1_c3": rand_seq(rng, 2600),            # unrelated
               "s1_c4": diverge(rng, dnab, 0.02)}       # DNA-B of Alpha
    d = tmp_path / "s1" / "assemble"
    d.mkdir(parents=True)
    (d / "contigs.fasta").write_text("".join(f">{k}\n{v}\n" for k, v in contigs.items()))
    cfg = Config.model_validate({"species": {"reference": str(tmp_path / "refs.fasta"),
                                             "taxon": "Testvirus", "strain_threshold": 0.94}})
    run_species(cfg, tmp_path / "s1", 2)
    rows = {l.split("\t")[0]: l.split("\t") for l in
            (tmp_path / "s1" / "species" / "species.tsv").read_text().splitlines()[1:]}
    assert rows["s1_c1"][3] == "Alpha virus" and rows["s1_c1"][8] == "known_species"
    assert rows["s1_c2"][3] == "Beta virus" and rows["s1_c2"][8] == "novel_species_candidate"
    assert rows["s1_c3"][8] == "no_reference_match"
    assert rows["s1_c4"][8] == "other_component" and rows["s1_c4"][5] == "B"


def test_skipped_without_reference(tmp_path):
    (tmp_path / "s1").mkdir()
    run_species(Config(), tmp_path / "s1", 1)
    assert (tmp_path / "s1" / "species" / ".done").exists()
