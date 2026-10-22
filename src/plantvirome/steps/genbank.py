"""Step 8: contig QC, origin rotation, ORF calling and GenBank-ready files.

Each contig gets status submit / review / reject with reasons (and non-blocking notes). Only "submit"
contigs enter submission.fsa and features.tbl (table2asn input). The organism
name and the "hypothetical protein" products are proposals: check them against
the literature and ICTV before submitting.
"""
import csv
from pathlib import Path

import pyrodigal_gv

from ..config import Config, GenBank
from ..logging import event
from ..runner import mark_done, step_done
from .assemble import read_fasta

COMP = str.maketrans("ACGTacgt", "TGCAtgca")
SEQID_MAX = 25  # GenBank limit


def revcomp(seq: str) -> str:
    return seq.translate(COMP)[::-1]


def orient_and_rotate(seq: str, motif: str, offset: int, max_mismatch: int = 0,
                      both_strands: bool = True) -> tuple[str, str]:
    """Rotate a circular sequence to start at motif[offset]; returns (seq, note).
    Looks for exact matches first, then up to max_mismatch substitutions (consensus errors);
    needs exactly one occurrence at the best level across both strands, else leaves it as is."""
    n, m = len(seq), len(motif)
    for k in range(max_mismatch + 1):
        hits = []
        for s in ((seq, revcomp(seq)) if both_strands else (seq,)):
            ext = s + s[:m - 1]
            hits += [(s, i) for i in range(n) if sum(a != b for a, b in zip(ext[i:i + m], motif)) == k]
        if hits:
            break
    else:
        return seq, "origin_not_found"
    if len(hits) > 1:
        return seq, "origin_ambiguous"
    s, pos = hits[0]
    start = (pos + offset) % n
    return s[start:] + s[:start], "origin_ok" if k == 0 else "origin_approx"


def rotate_to_reference(seq: str, strand: str, origin: int) -> str:
    """Orient and rotate a circular contig so it starts where its reference starts."""
    s = revcomp(seq) if strand == "-" else seq
    return s[origin:] + s[:origin]


def find_orfs(seq: str, min_aa: int) -> list[dict]:
    finder = pyrodigal_gv.ViralGeneFinder(meta=True)
    orfs = []
    for g in finder.find_genes(seq.encode()):
        protein = g.translate().rstrip("*")
        if len(protein) >= min_aa:
            orfs.append({"begin": g.begin, "end": g.end, "strand": g.strand, "protein": protein,
                         "partial_begin": g.partial_begin, "partial_end": g.partial_end})
    return orfs


def feature_table(seqid: str, orfs: list[dict]) -> str:
    """5-column feature table. Minus-strand features are written high -> low;
    partial ends get < / > (the lower coordinate is the 3' end on the minus strand)."""
    lines = [f">Feature {seqid}"]
    for o in orfs:
        plus = o["strand"] == 1
        low = f"{'<' if o['partial_begin'] and plus else '>' if o['partial_begin'] else ''}{o['begin']}"
        high = f"{'>' if o['partial_end'] and plus else '<' if o['partial_end'] else ''}{o['end']}"
        a, b = (low, high) if plus else (high, low)
        lines += [f"{a}\t{b}\tgene", "\t\t\tgene\torf" + f"{o['index']}",
                  f"{a}\t{b}\tCDS", "\t\t\tproduct\thypothetical protein"]
    return "\n".join(lines) + "\n"


def organism(row: dict, sp: dict | None = None, taxon: str | None = None) -> tuple[str, str | None]:
    """Proposed NCBI organism name, plus a warning if weak. A nucleotide-level species
    call (species step) wins over the protein-level classification."""
    if sp and sp.get("call") == "known_species":
        return sp["species"], None
    if sp and sp.get("call") == "other_component" and sp.get("species"):
        # DNA-B and satellites: protein hits cannot tell closely related species apart, so use
        # the nearest reference of that component; its cognate DNA-A decides the true species
        return sp["species"], "organism_from_nearest_component_reference"
    if sp and sp.get("call") == "novel_species_candidate" and taxon:
        return f"{taxon} sp.", "novel_species_candidate"
    name, rank = row.get("name", ""), row.get("rank", "")
    if rank == "species":
        return name, None
    if rank == "genus":
        return f"{name} sp.", "organism_proposed_from_genus"
    return "uncultured virus", "organism_unresolved"


def load_metadata(g: GenBank, sample: str) -> dict[str, str]:
    meta = dict(g.metadata)
    if g.metadata_csv:
        with open(g.metadata_csv) as fh:
            for r in csv.DictReader(fh):
                if r.pop("sample") == sample:
                    meta.update({k: v for k, v in r.items() if v})
    return meta


def _rows(path: Path, key: str) -> dict[str, dict]:
    if not path.exists():
        return {}
    with open(path) as fh:
        return {r[key]: r for r in csv.DictReader(fh, delimiter="\t")}


def run_genbank(cfg: Config, sampledir: Path, force: bool = False) -> None:
    out = sampledir / "genbank"
    if step_done(out) and not force:
        event("genbank skipped (done)", sample=sampledir.name)
        return
    out.mkdir(parents=True, exist_ok=True)
    g = cfg.genbank
    sample = sampledir.name
    seqs = read_fasta(sampledir / "assemble" / "contigs.fasta")
    info = _rows(sampledir / "assemble" / "contigs.tsv", "contig")
    tax = _rows(sampledir / "classify" / "taxonomy.tsv", "contig")
    ab = _rows(sampledir / "abundance" / "abundance.tsv", "contig")
    spec = _rows(sampledir / "species" / "species.tsv", "contig")
    meta = load_metadata(g, sample)
    # contigs whose vOTU is also in a negative control (written by the diversity step)
    contam = _rows(sampledir.parent / "cohort" / "diversity" / "contamination.tsv", "contig")

    qc_rows, fsa, tbl, faa = [], [], [], []
    for contig, seq in seqs.items():
        seq = seq.upper()
        circ = info.get(contig, {}).get("circular") == "True"
        t = tax.get(contig, {})
        a = ab.get(contig, {})
        sp = spec.get(contig)
        n_frac = sum(c not in "ACGT" for c in seq) / len(seq)
        depth = float(a.get("mean_depth", 0) or 0)
        # chimeric reads are not assigned to any contig, so compare with all reads touching it
        n_chim, n_reads = int(a.get("chimeric_reads", 0) or 0), int(a.get("reads", 0) or 0)
        chim = n_chim / max(1, n_reads + n_chim)
        seqid = f"{sample}_{contig.split('_')[-1]}" if not contig.startswith(sample) else contig
        reject, review, notes = [], [], []
        if len(seq) < g.min_length:
            reject.append("too_short")
        if n_frac > g.max_n_fraction:
            reject.append("too_many_ambiguous_bases")
        if depth < g.min_depth:
            reject.append("low_depth")
        if contam.get(contig, {}).get("flag") == "possible_contamination":
            review.append("possible_blank_contamination")
        if sp and "partial_genome" in sp.get("caveat", ""):
            review.append("shorter_than_reference")
        if sp and sp.get("call") == "no_reference_match":     # species step ran: no genome-like match
            review.append("no_reference_genome_match")
        if len(seqid) > SEQID_MAX:
            reject.append("seqid_too_long")
        # a nucleotide match to the reference genus is as good as a protein-tier hit
        if t.get("tier", "unclassified") not in g.allowed_tiers and not (
                sp and sp["call"] in ("known_species", "novel_species_candidate")):
            review.append("taxonomy_weak_or_missing")
        if not circ and not g.submit_partial:
            review.append("not_circular_partial_genome")
        if chim > g.max_chimeric_fraction:
            review.append("chimeric_reads_high")
        org, org_warn = organism(t, sp, cfg.species.taxon)
        if org_warn == "organism_unresolved":
            review.append(org_warn)
        elif org_warn:
            notes.append(org_warn)
        if sp and sp["caveat"]:
            notes.append("species_" + sp["caveat"].replace(",", ",species_"))
        note = ""
        if circ and not reject:
            ref_strand = sp.get("ref_strand") if g.orient_by_reference and sp else None
            if ref_strand:   # the reference fixes the strand; the motif then places the origin
                oriented = revcomp(seq) if ref_strand == "-" else seq
                note = "origin_not_found"
                if g.origin_motif:
                    rotated, note = orient_and_rotate(oriented, g.origin_motif, g.origin_offset,
                                                      g.origin_max_mismatch, both_strands=False)
                if note in ("origin_ok", "origin_approx"):
                    seq = rotated
                else:    # no usable motif: start where the reference starts
                    seq = rotate_to_reference(seq, ref_strand, int(sp["ref_origin"]))
                    note = "origin_from_reference"
                if note != "origin_ok":
                    notes.append(note)
            elif g.origin_motif:
                seq, note = orient_and_rotate(seq, g.origin_motif, g.origin_offset, g.origin_max_mismatch)
                if note == "origin_approx":
                    notes.append(note)
                elif note != "origin_ok":
                    review.append(note)
        status = "reject" if reject else "review" if review else "submit"
        n_orfs = 0
        if status == "submit":
            orfs = find_orfs(seq, g.min_orf_aa)
            for i, o in enumerate(orfs, 1):
                o["index"] = i
            n_orfs = len(orfs)
            mods = {"organism": org, "moltype": "genomic", "molecule": "DNA",
                    "topology": "circular" if circ else "linear", **meta}
            fsa.append(f">{seqid} " + " ".join(f"[{k}={v}]" for k, v in mods.items())
                       + f"\n{seq}\n")
            tbl.append(feature_table(seqid, orfs))
            faa += [f">{seqid}_orf{o['index']}\n{o['protein']}\n" for o in orfs]
        qc_rows.append([contig, seqid, len(seq), circ, f"{n_frac:.4f}", f"{depth:.1f}",
                        t.get("tier", "unclassified"), org, status, n_orfs,
                        ",".join(reject + review + notes)])
    (out / "submission.fsa").write_text("".join(fsa))
    (out / "features.tbl").write_text("".join(tbl))
    (out / "orfs.faa").write_text("".join(faa))
    with open(out / "contigs.qc.tsv", "w") as fh:
        fh.write("contig\tseqid\tlength\tcircular\tn_fraction\tmean_depth\ttier\torganism\t"
                 "status\torfs\treasons\n")
        fh.writelines("\t".join(map(str, r)) + "\n" for r in qc_rows)
    event("genbank done", sample=sample, submit=len(fsa), total=len(qc_rows))
    mark_done(out)
