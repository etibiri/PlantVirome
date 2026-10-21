"""Step 5: contig classification by protein homology (DIAMOND blastx) + LCA (taxonkit).

Protein-level search is used because divergent plant viruses share little
nucleotide identity. Contigs without a hit stay "unclassified": they may be
novel viruses, satellites without proteins, or host/organelle sequence.
"""
import subprocess
from collections import defaultdict
from pathlib import Path

from ..config import Config
from ..logging import event
from ..runner import ToolError, mark_done, run, step_done
from .assemble import read_fasta

HEADER = ["contig", "length", "taxid", "name", "rank", "lineage",
          "best_pident", "query_cover", "evalue", "tier"]
COLS = "qseqid sseqid pident length qstart qend qlen bitscore evalue staxids".split()


def parse_hits(path: Path) -> dict[str, list[dict]]:
    hits: dict[str, list[dict]] = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            f = dict(zip(COLS, line.rstrip("\n").split("\t")))
            hits[f["qseqid"]].append({
                "pident": float(f["pident"]), "bits": float(f["bitscore"]),
                "evalue": float(f["evalue"]), "qlen": int(f["qlen"]),
                "span": tuple(sorted((int(f["qstart"]), int(f["qend"])))),
                "taxids": {t for t in f["staxids"].split(";") if t and t != "0"}})
    return hits


def summarize(hits: list[dict], top_pct: float) -> dict:
    """Best hit stats, taxids of the top-scoring hits (for LCA), and query coverage."""
    best = max(hits, key=lambda h: h["bits"])
    keep = [h for h in hits if h["bits"] >= best["bits"] * (1 - top_pct / 100)]
    covered, end = 0, 0
    for s, e in sorted(h["span"] for h in hits):
        covered += max(0, e - max(s - 1, end))
        end = max(end, e)
    return {"pident": best["pident"], "evalue": best["evalue"],
            "cover": min(1.0, covered / best["qlen"]),
            "taxids": set().union(*(h["taxids"] for h in keep))}


def tier(pident: float, tiers: list[tuple[float, str]]) -> str:
    for min_id, label in sorted(tiers, reverse=True):
        if pident >= min_id:
            return label
    return "unclassified"


def taxonkit(args: list[str], lines: list[str], taxdump: Path) -> list[list[str]]:
    out = subprocess.run(["taxonkit", *args, "--data-dir", str(taxdump)], input="\n".join(lines) + "\n",
                         capture_output=True, text=True)
    if out.returncode:
        raise ToolError(f"taxonkit failed: {out.stderr[-1000:]}")
    return [l.split("\t") for l in out.stdout.splitlines()]


def run_classify(cfg: Config, sampledir: Path, threads: int, force: bool = False) -> None:
    out = sampledir / "classify"
    if step_done(out) and not force:
        event("classify skipped (done)", sample=sampledir.name)
        return
    out.mkdir(parents=True, exist_ok=True)
    try:
        db, taxdump = cfg.databases["viral_protein"], cfg.databases["taxdump"]
    except KeyError as e:
        raise ToolError(f"config databases.{e.args[0]} is required for classify") from None
    contigs = sampledir / "assemble" / "contigs.fasta"
    lengths = {k: len(v) for k, v in read_fasta(contigs).items()}
    c = cfg.classify
    tsv = out / "taxonomy.tsv"
    hits: dict[str, list[dict]] = {}
    if lengths:
        raw = out / "diamond.tsv"
        run(["diamond", "blastx", f"--{c.sensitivity}", "-d", db, "-q", contigs, "-o", raw,
             "--outfmt", "6", *COLS, "-e", c.evalue, "-p", threads, "--frameshift", 15,
             "--max-target-seqs", 25])
        hits = parse_hits(raw)
        if hits and not any(h["taxids"] for hs in hits.values() for h in hs):
            raise ToolError("DIAMOND hits carry no taxids: build the database with "
                            "--taxonmap and --taxonnodes")

    summ = {k: summarize(v, c.top_pct) for k, v in hits.items()}
    lca = {}
    if summ:  # one LCA per contig, then lineage/name/rank for each distinct taxid
        keys = list(summ)
        res = taxonkit(["lca"], [" ".join(summ[k]["taxids"]) for k in keys], taxdump)
        lca = dict(zip(keys, (r[1] for r in res)))
    info = {r[0]: r[1:] for r in taxonkit(["lineage", "-n", "-r"], sorted(set(lca.values())), taxdump)} \
        if lca else {}

    with open(tsv, "w") as fh:
        fh.write("\t".join(HEADER) + "\n")
        for contig, length in lengths.items():
            s = summ.get(contig)
            t = tier(s["pident"], c.tiers) if s else "unclassified"
            if s and t != "unclassified":
                taxid = lca[contig]
                lineage, name, rank = (info.get(taxid, ["", "", ""]) + ["", "", ""])[:3]
                row = [taxid, name, rank, lineage, f"{s['pident']:.1f}", f"{s['cover']:.2f}", f"{s['evalue']:.1e}"]
            else:
                row = ["0", "unclassified", "", "", f"{s['pident']:.1f}" if s else "", f"{s['cover']:.2f}" if s else "", ""]
            fh.write("\t".join([contig, str(length), *row, t]) + "\n")
    event("classify done", sample=sampledir.name, contigs=len(lengths), with_hits=len(summ))
    mark_done(out)
