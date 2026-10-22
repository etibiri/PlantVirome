"""Step 5c: detection of satellite DNAs (alpha-, beta-, deltasatellites ...).

Satellites are circular, 1-1.4 kb, and share little with their helper virus, so they are
searched separately: contigs inside a group's size window are aligned to that group's
reference genomes (`download-refs --taxon Alphasatellitidae`, `Tolecusatellitidae`, ...) and
called like species in the species step:

  known_satellite            identity >= the group's species_threshold over >= min_cover
  novel_satellite_candidate  below it (nearest species reported)
  satellite_related          the group has no threshold: identity reported only
  protein_only_candidate     no nucleotide match, but the protein search (classify) put the
                             contig in a satellite family
  (no row)                   nothing satellite-like

Limits: nucleotide search finds satellites down to ~85 % identity to a reference (alignments
fragment below that), and alphasatellite species differ by >12 %, so a divergent alphasatellite
is found only through its Rep protein (protein_only_candidate) or not at all; a betasatellite's
single short protein rarely gives a protein hit. Size alone proves nothing: defective DNAs and
phi29 artifacts also fall in the 1-2 kb window. Satellite species thresholds are ICTV's and
vary by subfamily: set them per group in the config.
"""
import csv
import re
from pathlib import Path

from ..config import Config, Species
from ..logging import event
from ..runner import mark_done, run, step_done
from .assemble import read_fasta
from .species import call, parse_paf, read_refs, reference_origin

HEADER = ("contig\tlength\tgroup\tbest_ref\tspecies\ttaxid\tidentity\tcover\tcall\tcaveat\t"
          "ref_strand\tref_origin\tref_length\tref_cover\tprotein_hint\n")
SAT_LINEAGE = re.compile(r"satelli|Alphasatell|Tolecusatell|Deltasatell", re.IGNORECASE)


def satellite_call(species_call: str, has_threshold: bool) -> str:
    if species_call == "known_species":
        return "known_satellite" if has_threshold else "satellite_related"
    return "novel_satellite_candidate" if has_threshold else "satellite_related"


def protein_hints(classify_tsv: Path) -> dict[str, str]:
    """contig -> protein taxon, for contigs the protein search put in a satellite family."""
    if not classify_tsv.exists():
        return {}
    with open(classify_tsv) as fh:
        return {r["contig"]: r["name"] for r in csv.DictReader(fh, delimiter="\t")
                if SAT_LINEAGE.search(r["lineage"] + " " + r["name"])}


def run_satellite(cfg: Config, sampledir: Path, threads: int, force: bool = False) -> None:
    out = sampledir / "satellite"
    if step_done(out) and not force:
        event("satellite skipped (done)", sample=sampledir.name)
        return
    out.mkdir(parents=True, exist_ok=True)
    s = cfg.satellite
    if not s.groups:
        event("no satellite.groups configured; satellite step skipped", sample=sampledir.name)
        mark_done(out)
        return
    contigs = read_fasta(sampledir / "assemble" / "contigs.fasta")
    hints = protein_hints(sampledir / "classify" / "taxonomy.tsv")

    best: dict[str, tuple] = {}   # contig -> (sort key, row fields)
    for name, g in s.groups.items():
        sub = {k: v for k, v in contigs.items() if g.min_length <= len(v) <= g.max_length}
        if not sub:
            continue
        fa, paf = out / f"{name}.fa", out / f"{name}.paf"
        fa.write_text("".join(f">{k}\n{v}\n" for k, v in sub.items()))
        run(["minimap2", "-x", "asm20", "-c", "-N", 100, "-p", 0.3, "-t", threads, g.reference, fa],
            stdout=paf)
        refs = read_refs(g.reference, component="A")
        ref_lengths = {k: len(v) for k, v in read_fasta(g.reference).items()}
        hits, pieces = parse_paf(paf, {k: len(v) for k, v in sub.items()}, ref_lengths)
        fa.unlink()
        paf.unlink()
        sp = Species(species_threshold=g.species_threshold or 2.0, min_cover=s.min_cover,
                     min_partial_cover=s.min_partial_cover, min_ref_cover=s.min_ref_cover,
                     borderline_margin=s.borderline_margin)
        for contig, seq in sub.items():
            c = call(hits.get(contig, {}), refs, sp)
            if c["call"] == "no_reference_match":
                continue
            strand, pos = reference_origin(pieces[(contig, c["ref"])], len(seq))
            r = refs[c["ref"]]
            row = [contig, len(seq), name, r["acc"], r["species"], r["taxid"], f"{c['identity']:.3f}",
                   f"{c['cover']:.2f}", satellite_call(c["call"], g.species_threshold is not None),
                   c["caveat"], strand, pos, ref_lengths[c["ref"]], f"{c['ref_cover']:.2f}",
                   hints.get(contig, "")]
            key = (c["identity"], c["cover"])
            if contig not in best or key > best[contig][0]:
                best[contig] = (key, row)

    with open(out / "satellites.tsv", "w") as fh:
        fh.write(HEADER)
        for contig, (_, row) in best.items():
            fh.write("\t".join(map(str, row)) + "\n")
        for contig, taxon in hints.items():       # protein evidence without a nucleotide match
            if contig not in best:
                fh.write("\t".join(map(str, [contig, len(contigs[contig]), "", "", "", "", "", "",
                                             "protein_only_candidate", "", "", "", "", "", taxon])) + "\n")
    event("satellite done", sample=sampledir.name, nucleotide=len(best),
          protein_only=sum(1 for c in hints if c not in best))
    mark_done(out)
