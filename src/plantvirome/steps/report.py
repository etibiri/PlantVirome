"""Final step: one self-contained HTML report (summary tables + inline SVG figures).

Reads the tables the other steps already wrote; needs no plotting library.
Anything missing (e.g. no cohort diversity yet) is skipped, not an error.
"""
import csv
import statistics
from collections import Counter, defaultdict
from html import escape
from pathlib import Path

from .. import __version__
from ..logging import event

PALETTE = ["#4e79a7", "#f28e2b", "#e15759", "#76b7b2", "#59a14f", "#edc948", "#b07aa1", "#9c755f", "#bab0ac"]
STATUS_COLOR = {"submit": "#59a14f", "review": "#edc948", "reject": "#e15759"}


def read_tsv(path: Path) -> list[dict]:
    if not path.is_file():
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def table(rows: list[list], header: list[str]) -> str:
    head = "".join(f"<th>{escape(str(h))}</th>" for h in header)
    body = "".join("<tr>" + "".join(f"<td>{escape(str(c))}</td>" for c in r) + "</tr>" for r in rows)
    return f"<div class=scroll><table><tr>{head}</tr>{body}</table></div>"


def svg(w: int, h: int, inner: str) -> str:
    return (f'<svg viewBox="0 0 {w} {h}" width="100%" style="max-width:{w}px" '
            f'font-family="sans-serif" font-size="11">{inner}</svg>')


def hbars(items: list[tuple[str, float]], color: str = PALETTE[0]) -> str:
    """Horizontal bar chart of (label, value)."""
    if not items:
        return ""
    m = max(v for _, v in items) or 1
    inner = ""
    for i, (k, v) in enumerate(items):
        y = i * 20 + 4
        inner += (f'<text x="268" y="{y + 12}" text-anchor="end">{escape(str(k)[:44])}</text>'
                  f'<rect x="274" y="{y}" width="{max(v / m * 380, 1):.1f}" height="15" fill="{color}"/>'
                  f'<text x="{278 + v / m * 380:.1f}" y="{y + 12}">{v:g}</text>')
    return svg(760, len(items) * 20 + 8, inner)


def vbars(items: list[tuple[str, dict[str, float]]], colors: dict[str, str], legend: bool = True) -> str:
    """Vertical stacked bars: one bar per label, segments by key."""
    if not items:
        return ""
    w, h, top, bot = max(660, len(items) * 14 + 50), 260, 10, 70
    m = max(sum(d.values()) for _, d in items) or 1
    bw = (w - 40) / len(items)
    inner = f'<line x1="36" y1="{h - bot}" x2="{w}" y2="{h - bot}" stroke="#888"/>'
    inner += f'<text x="32" y="{top + 8}" text-anchor="end">{m:g}</text><text x="32" y="{h - bot}" text-anchor="end">0</text>'
    for i, (k, d) in enumerate(items):
        x, y = 38 + i * bw, h - bot
        for seg, v in d.items():
            hh = v / m * (h - bot - top)
            y -= hh
            inner += f'<rect x="{x:.1f}" y="{y:.1f}" width="{max(bw - 2, 1):.1f}" height="{hh:.1f}" fill="{colors.get(seg, "#999")}"/>'
        inner += (f'<text transform="translate({x + bw / 2:.1f},{h - bot + 8}) rotate(60)" '
                  f'font-size="9">{escape(str(k))}</text>')
    key = "".join(f'<span><i style="background:{c}"></i>{escape(s)}</span>' for s, c in colors.items())
    return svg(w, h, inner) + (f'<div class=legend>{key}</div>' if legend else "")


def histogram(values: list[float], bins: int = 30) -> str:
    if not values:
        return ""
    lo, hi = min(values), max(values)
    step = (hi - lo) / bins or 1
    counts = Counter(min(int((v - lo) / step), bins - 1) for v in values)
    items = [(f"{lo + i * step:.0f}", {"n": counts.get(i, 0)}) for i in range(bins)]
    return vbars(items, {"n": PALETTE[0]}, legend=False)


def scatter(points: list[tuple[str, float, float]], xl: str, yl: str) -> str:
    if not points:
        return ""
    w, h, pad = 560, 420, 40
    xs, ys = [p[1] for p in points], [p[2] for p in points]
    x0, x1, y0, y1 = min(xs), max(xs), min(ys), max(ys)
    sx = lambda v: pad + (v - x0) / ((x1 - x0) or 1) * (w - 2 * pad)
    sy = lambda v: h - pad - (v - y0) / ((y1 - y0) or 1) * (h - 2 * pad)
    inner = (f'<rect x="{pad}" y="{pad}" width="{w - 2 * pad}" height="{h - 2 * pad}" fill="none" stroke="#888"/>'
             f'<text x="{w / 2}" y="{h - 8}" text-anchor="middle">{escape(xl)}</text>'
             f'<text transform="translate(12,{h / 2}) rotate(-90)" text-anchor="middle">{escape(yl)}</text>')
    for s, x, y in points:
        inner += (f'<circle cx="{sx(x):.1f}" cy="{sy(y):.1f}" r="4" fill="{PALETTE[0]}" fill-opacity=".7">'
                  f'<title>{escape(s)}</title></circle>')
    return svg(w, h, inner)


def top(counter: Counter, n: int = 12) -> list[tuple[str, float]]:
    return [(k or "none", v) for k, v in counter.most_common(n)]


def build(outdir: Path) -> str:
    samples = sorted(p.parent.parent.name for p in outdir.glob("*/genbank/contigs.qc.tsv"))
    qc, species, ab = [], [], []
    for s in samples:
        for dst, rel in ((qc, "genbank/contigs.qc.tsv"), (species, "species/species.tsv"), (ab, "abundance/abundance.tsv")):
            dst += [{"sample": s, **r} for r in read_tsv(outdir / s / rel)]
    div = outdir / "cohort" / "diversity"
    alpha, pcoa, votus = read_tsv(div / "alpha.tsv"), read_tsv(div / "pcoa_bray_curtis.tsv"), read_tsv(div / "votu_table.tsv")
    contam = read_tsv(div / "contamination.tsv")

    status = Counter(r["status"] for r in qc)
    per_sample: dict[str, Counter] = defaultdict(Counter)
    for r in qc:
        per_sample[r["sample"]][r["status"]] += 1
    detected = {r["contig"] for r in ab if r["detected"] == "True"}

    sec = []
    sec.append("<h2>Overview</h2>" + table([
        ["Samples", len(samples)], ["Contigs", len(qc)],
        ["submit / review / reject", " / ".join(str(status[k]) for k in ("submit", "review", "reject"))],
        ["Contigs detected (abundance)", len(detected)], ["vOTUs", len(votus)],
        ["Possible contamination (vs blanks)", sum(r["flag"] != "ok" for r in contam)],
    ], ["Metric", "Value"]))

    sec.append("<h2>Contig quality</h2><h3>Status per sample</h3>"
               + vbars([(s, dict(per_sample[s])) for s in samples], STATUS_COLOR)
               + "<h3>Contig length (bp)</h3>" + histogram([float(r["length"]) for r in qc])
               + "<h3>Summary per sample</h3>" + table(
                   [[s, sum(per_sample[s].values())] + [per_sample[s][k] for k in ("submit", "review", "reject")]
                    for s in samples], ["sample", "contigs", "submit", "review", "reject"]))
    if any(r.get("checkv_quality") for r in qc):
        sec.append("<h3>CheckV quality</h3>" + hbars(top(Counter(r["checkv_quality"] for r in qc if r["checkv_quality"])), PALETTE[4])
                   + "<p class=note>Informational: CheckV compares with reference genomes and underestimates the "
                     "completeness of small circular genomes (e.g. begomovirus DNA-A/B). Only contamination "
                     "and host genes affect the status.</p>")
    reasons = Counter(x for r in qc for x in r["reasons"].split(",") if x)
    sec.append("<h3>Most frequent QC reasons</h3>" + hbars(top(reasons), PALETTE[1]))

    if species:
        sec.append("<h2>Species assignment</h2><h3>Call</h3>" + hbars(top(Counter(r["call"] for r in species)), PALETTE[3])
                   + "<h3>Best reference species</h3>" + hbars(top(Counter(r["species"] for r in species)), PALETTE[4])
                   + "<h3>Identity to best reference, by call</h3>"
                   + table([[c, n, f"{statistics.median(v):.3f}", f"{min(v):.3f}", f"{max(v):.3f}"]
                            for c, v in sorted(_by(species, "call", "identity").items()) for n in [len(v)]],
                           ["call", "contigs", "median identity", "min", "max"]))

    if ab:
        sec.append("<h2>Abundance</h2><h3>Mean depth per contig (log10)</h3>"
                   + histogram([_log10(float(r["mean_depth"])) for r in ab])
                   + table([["contigs", len(ab)], ["detected", len(detected)],
                            ["median depth", f"{statistics.median(float(r['mean_depth']) for r in ab):.1f}"]],
                           ["metric", "value"]))

    if alpha:
        sh = [float(r["shannon"]) for r in alpha]
        sec.append("<h2>Diversity</h2><h3>Shannon index per sample</h3>"
                   + vbars([(r["sample"], {"shannon": float(r["shannon"])}) for r in alpha], {"shannon": PALETTE[0]}, legend=False)
                   + table([["observed vOTUs", f"{statistics.median(float(r['observed']) for r in alpha):g}"],
                            ["Shannon (median)", f"{statistics.median(sh):.3f}"], ["Shannon (max)", f"{max(sh):.3f}"]],
                           ["alpha diversity", "median across samples"]))
    if pcoa:
        k = list(pcoa[0])
        sec.append(f"<h3>PCoA (Bray-Curtis)</h3>"
                   + scatter([(r["sample"], float(r[k[1]]), float(r[k[2]])) for r in pcoa], k[1], k[2]))
    if votus:
        sec.append("<h3>vOTUs by taxon</h3>" + hbars(top(Counter(r["taxon"] for r in votus)), PALETTE[6])
                   + table([[r["votu"], r["taxon"], r["tier"], r["length"], r["contigs"]] for r in votus[:50]],
                           ["vOTU", "taxon", "tier", "length", "contigs"])
                   + ("<p class=note>First 50 vOTUs; full table in cohort/diversity/votu_table.tsv.</p>" if len(votus) > 50 else ""))
    if contam:
        sec.append("<h3>Blank check</h3>" + hbars(top(Counter(r["flag"] for r in contam)), PALETTE[2]))
    return PAGE.format(version=__version__, outdir=escape(str(outdir)), body="\n".join(sec))


def _by(rows: list[dict], key: str, val: str) -> dict[str, list[float]]:
    out: dict[str, list[float]] = defaultdict(list)
    for r in rows:
        if r.get(val):
            out[r[key]].append(float(r[val]))
    return {k: v for k, v in out.items() if v}


def _log10(v: float) -> float:
    import math
    return math.log10(v) if v > 0 else 0.0


PAGE = """<!doctype html><html><head><meta charset=utf-8><title>plantvirome report</title>
<style>body{{font:14px/1.5 sans-serif;max-width:900px;margin:2em auto;padding:0 1em;color:#222}}
table{{border-collapse:collapse;margin:.5em 0}}th,td{{border:1px solid #ccc;padding:2px 8px;text-align:left}}
th{{background:#f0f0f0;position:sticky;top:0}}.scroll{{max-height:420px;overflow:auto;display:inline-block}}
.legend span{{margin-right:1em}}.legend i{{display:inline-block;width:10px;height:10px;margin-right:4px}}
.note{{color:#666}}h2{{border-bottom:2px solid #ddd;margin-top:2em}}</style></head><body>
<h1>plantvirome report</h1><p class=note>plantvirome {version} &middot; {outdir}</p>
{body}</body></html>"""


def run_report(outdir: Path, out: Path | None = None) -> Path:
    out = out or outdir / "report.html"
    out.write_text(build(outdir))
    event("report done", path=str(out))
    return out
