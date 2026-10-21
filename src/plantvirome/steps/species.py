"""Step 5b: nucleotide-level species assignment against reference genomes.

Protein homology (`classify`) often stops at genus because related species share
near-identical proteins. Many virus genera define species by whole-genome pairwise
nucleotide identity instead (begomovirus DNA-A: 91 %), so each contig is aligned
to the references and the best match decides:

  known_species            identity >= species_threshold over >= min_cover of the contig
  novel_species_candidate  identity below it (nearest species reported)
  no_reference_match       no usable alignment to any reference

Identity here is matches / alignment columns from minimap2 (asm20), summed over the
pieces of a circular genome that wraps the contig ends. It tracks true identity well
from ~88 % up but reads ~0.5 % low at 91 % (gaps count) and the alignment fragments
below ~86 %, where coverage falls and the call becomes a partial-alignment caveat.
Treat calls within `borderline_margin` of the threshold, and every novel candidate,
as needing confirmation with a proper multiple alignment (e.g. SDT/MAFFT).
Reference headers are `accession|species_name|taxid|component` (underscores for spaces); only component `A`
(DNA-A, or a monopartite genome) is used for the threshold. Other components
(DNA-B, satellites) get `other_component` and no species decision.
"""
from collections import defaultdict
from pathlib import Path

from ..config import Config, Species
from ..logging import event
from ..runner import ToolError, mark_done, run, step_done
from .assemble import read_fasta

HEADER = "contig\tlength\tbest_ref\tspecies\ttaxid\tcomponent\tidentity\tcover\tcall\tcaveat\n"


def parse_paf(paf: Path, lengths: dict[str, int]) -> dict[str, dict[str, tuple[float, float]]]:
    """contig -> ref -> (identity, cover), summing all alignment pieces of the pair."""
    acc: dict[tuple[str, str], list[int]] = defaultdict(lambda: [0, 0, 0])
    with open(paf) as fh:
        for line in fh:
            f = line.split("\t")
            a = acc[(f[0], f[5])]
            a[0] += int(f[3]) - int(f[2])
            a[1] += int(f[9])
            a[2] += int(f[10])
    out: dict[str, dict[str, tuple[float, float]]] = defaultdict(dict)
    for (q, t), (aligned, matches, block) in acc.items():
        out[q][t] = (matches / block if block else 0.0, min(1.0, aligned / lengths[q]))
    return out


def call(hits: dict[str, tuple[float, float]], refs: dict[str, dict], sp: Species) -> dict:
    """Pick the best reference and make the species call for one contig."""
    if not hits:
        return {"ref": "", "identity": 0.0, "cover": 0.0, "call": "no_reference_match", "caveat": ""}
    full = {r: h for r, h in hits.items() if h[1] >= sp.min_cover}
    pool = full or hits
    # best = highest identity among full-length comparisons, else the most covered
    ref = max(pool, key=lambda r: (pool[r][0], pool[r][1]) if full else (pool[r][1], pool[r][0]))
    identity, cover = hits[ref]
    out = {"ref": ref, "identity": identity, "cover": cover, "caveat": ""}
    if cover < sp.min_partial_cover:
        return {**out, "call": "no_reference_match"}
    if refs[ref]["component"] != "A":
        return {**out, "call": "other_component", "caveat": f"component_{refs[ref]['component']}"}
    out["call"] = "known_species" if full and identity >= sp.species_threshold else "novel_species_candidate"
    caveats = []
    if not full:
        caveats.append("partial_alignment")
    if abs(identity - sp.species_threshold) < sp.borderline_margin:
        caveats.append("borderline")
    if out["call"] == "known_species" and sp.strain_threshold and identity < sp.strain_threshold:
        caveats.append("new_strain")
    out["caveat"] = ",".join(caveats)
    return out


def read_refs(path: Path) -> dict[str, dict]:
    refs = {}
    for header in read_fasta(path):
        acc, species, taxid, component = header.split("|")
        refs[header] = {"acc": acc, "species": species.replace("_", " "), "taxid": taxid, "component": component}
    return refs


def run_species(cfg: Config, sampledir: Path, threads: int, force: bool = False) -> None:
    out = sampledir / "species"
    if step_done(out) and not force:
        event("species skipped (done)", sample=sampledir.name)
        return
    out.mkdir(parents=True, exist_ok=True)
    sp = cfg.species
    if sp.reference is None:
        event("no species.reference configured; species step skipped", sample=sampledir.name)
        mark_done(out)
        return
    refs = read_refs(sp.reference)
    if not refs:
        raise ToolError(f"no reference sequences in {sp.reference}")
    contigs = sampledir / "assemble" / "contigs.fasta"
    lengths = {k: len(v) for k, v in read_fasta(contigs).items()}
    hits: dict[str, dict] = {}
    if lengths:
        paf = out / "contigs_vs_refs.paf"
        # many similar references: keep up to 100 alignments per contig, not just the top 5
        run(["minimap2", "-x", "asm20", "-c", "-N", 100, "-p", 0.3, "-t", threads,
             sp.reference, contigs], stdout=paf)
        hits = parse_paf(paf, lengths)
        paf.unlink()
    with open(out / "species.tsv", "w") as fh:
        fh.write(HEADER)
        for contig, length in lengths.items():
            c = call(hits.get(contig, {}), refs, sp)
            r = refs.get(c["ref"], {"species": "", "taxid": "", "component": "", "acc": ""})
            fh.write(f"{contig}\t{length}\t{r['acc']}\t{r['species']}\t{r['taxid']}\t{r['component']}\t"
                     f"{c['identity']:.3f}\t{c['cover']:.2f}\t{c['call']}\t{c['caveat']}\n")
    event("species done", sample=sampledir.name, contigs=len(lengths))
    mark_done(out)
