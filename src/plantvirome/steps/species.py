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

HEADER = ("contig\tlength\tbest_ref\tspecies\ttaxid\tcomponent\tidentity\tcover\tcall\tcaveat\t"
          "ref_strand\tref_origin\tref_length\tref_cover\n")


def parse_paf(paf: Path, lengths: dict[str, int], ref_lengths: dict[str, int]):
    """contig -> ref -> (identity, cover, ref_cover) summing all alignment pieces of the pair,
    and (contig, ref) -> [(contig_start, contig_end, ref_start, strand)] for the origin.
    `ref_lengths` gives ref_cover (fraction of the reference covered by the contig)."""
    acc: dict[tuple[str, str], list[int]] = defaultdict(lambda: [0, 0, 0, 0])
    pieces: dict[tuple[str, str], list[tuple]] = defaultdict(list)
    with open(paf) as fh:
        for line in fh:
            f = line.split("\t")
            a = acc[(f[0], f[5])]
            a[0] += int(f[3]) - int(f[2])
            a[1] += int(f[9])
            a[2] += int(f[10])
            a[3] += int(f[8]) - int(f[7])
            pieces[(f[0], f[5])].append((int(f[2]), int(f[3]), int(f[7]), f[4]))
    out: dict[str, dict[str, tuple[float, float, float]]] = defaultdict(dict)
    for (q, t), (aligned, matches, block, ref_aligned) in acc.items():
        out[q][t] = (matches / block if block else 0.0, min(1.0, aligned / lengths[q]),
                     min(1.0, ref_aligned / ref_lengths[t]))
    return out, pieces


def reference_origin(pieces: list[tuple], n: int) -> tuple[str, int]:
    """Where the reference's position 1 falls in the contig: ('+', pos) when the contig
    is already on the reference strand, ('-', pos) when its reverse complement is, with pos
    the 0-based start in that orientation. Uses the alignment piece nearest the reference
    start; exact when the contig and reference have no indels in between."""
    qs, qe, rs, strand = min(pieces, key=lambda p: p[2])
    if strand == "+":
        return "+", (qs - rs) % n
    return "-", (n - qe - rs) % n


def call(hits: dict[str, tuple[float, float, float]], refs: dict[str, dict], sp: Species) -> dict:
    """Pick the best reference and make the species call for one contig."""
    if not hits:
        return {"ref": "", "identity": 0.0, "cover": 0.0, "ref_cover": 0.0, "call": "no_reference_match",
                "caveat": ""}
    full = {r: h for r, h in hits.items() if h[1] >= sp.min_cover}
    pool = full or hits
    # best = highest identity among full-length comparisons, else the most covered
    ref = max(pool, key=lambda r: (pool[r][0], pool[r][1]) if full else (pool[r][1], pool[r][0]))
    identity, cover, ref_cover = hits[ref]
    out = {"ref": ref, "identity": identity, "cover": cover, "ref_cover": ref_cover, "caveat": ""}
    if cover < sp.min_partial_cover:
        return {**out, "call": "no_reference_match"}
    if refs[ref]["component"] != "A":
        caveats = [f"component_{refs[ref]['component']}"]
        if ref_cover < sp.min_ref_cover:
            caveats.append("partial_genome")
        return {**out, "call": "other_component", "caveat": ",".join(caveats)}
    out["call"] = "known_species" if full and identity >= sp.species_threshold else "novel_species_candidate"
    caveats = []
    if ref_cover < sp.min_ref_cover:
        caveats.append("partial_genome")
    if not full:
        caveats.append("partial_alignment")
    if abs(identity - sp.species_threshold) < sp.borderline_margin:
        caveats.append("borderline")
    if out["call"] == "known_species" and sp.strain_threshold and identity < sp.strain_threshold:
        caveats.append("new_strain")
    out["caveat"] = ",".join(caveats)
    return out


def read_refs(path: Path, component: str | None = None) -> dict[str, dict]:
    """Reference headers -> metadata. `component` overrides the header's component label
    (the satellite step treats every reference as a complete genome)."""
    refs = {}
    for header in read_fasta(path):
        acc, species, taxid, comp = header.split("|")
        component_ = component or comp
        refs[header] = {"acc": acc, "species": species.replace("_", " "), "taxid": taxid, "component": component_}
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
    ref_lengths = {k: len(v) for k, v in read_fasta(sp.reference).items()}
    if not refs:
        raise ToolError(f"no reference sequences in {sp.reference}")
    contigs = sampledir / "assemble" / "contigs.fasta"
    lengths = {k: len(v) for k, v in read_fasta(contigs).items()}
    hits: dict[str, dict] = {}
    pieces: dict = {}
    if lengths:
        paf = out / "contigs_vs_refs.paf"
        # many similar references: keep up to 100 alignments per contig, not just the top 5
        run(["minimap2", "-x", "asm20", "-c", "-N", 100, "-p", 0.3, "-t", threads,
             sp.reference, contigs], stdout=paf)
        hits, pieces = parse_paf(paf, lengths, ref_lengths)
        paf.unlink()
    with open(out / "species.tsv", "w") as fh:
        fh.write(HEADER)
        for contig, length in lengths.items():
            c = call(hits.get(contig, {}), refs, sp)
            r = refs.get(c["ref"], {"species": "", "taxid": "", "component": "", "acc": ""})
            fh.write(f"{contig}\t{length}\t{r['acc']}\t{r['species']}\t{r['taxid']}\t{r['component']}\t"
                     f"{c['identity']:.3f}\t{c['cover']:.2f}\t{c['call']}\t{c['caveat']}\t")
            if c["call"] != "no_reference_match":
                strand, pos = reference_origin(pieces[(contig, c["ref"])], length)
                fh.write(f"{strand}\t{pos}\t{ref_lengths[c['ref']]}\t{c['ref_cover']:.2f}\n")
            else:
                fh.write("\t\t\t\n")
    event("species done", sample=sampledir.name, contigs=len(lengths))
    mark_done(out)
