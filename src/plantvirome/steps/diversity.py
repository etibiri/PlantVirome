"""Step 7 (cohort): vOTU clustering, alpha and beta diversity, PCoA.

Features are vOTUs: contigs of all samples clustered by nucleotide similarity, so
novel viruses count even without a database hit; taxonomy is only an annotation.
Abundance is mean depth from the abundance step.

Limits (also stated in the README): RCA selectively amplifies small circular DNA
and misses most RNA viruses, so this is the diversity of the RCA-amplifiable
virome, not of the whole virome. Depth is semi-quantitative, Chao1 treats reads as
individuals (concatemer reads of unequal length violate that), and
the vOTU identity threshold sets what "richness" means.
"""
import csv
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.spatial.distance import pdist, squareform

from ..config import Config
from ..logging import event
from ..runner import ToolError, mark_done, run, step_done
from .assemble import read_fasta


def cluster(paf: Path, lengths: dict[str, int], min_identity: float, min_cover: float) -> dict[str, str]:
    """Union-find over contig pairs passing identity/coverage; returns contig -> vOTU id."""
    pairs: dict[tuple[str, str], list[int]] = defaultdict(lambda: [0, 0, 0])
    with open(paf) as fh:
        for line in fh:
            f = line.split("\t")
            q, t = f[0], f[5]
            if q == t:
                continue
            p = pairs[(q, t)]
            p[0] += int(f[3]) - int(f[2])
            p[1] += int(f[9])
            p[2] += int(f[10])
    parent = {c: c for c in lengths}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for (q, t), (aligned, matches, block) in pairs.items():
        if block and matches / block >= min_identity \
                and min(1.0, aligned / max(lengths[q], lengths[t])) >= min_cover:
            parent[find(q)] = find(t)
    groups: dict[str, list[str]] = defaultdict(list)
    for c in lengths:
        groups[find(c)].append(c)
    # Number vOTUs by decreasing size of their longest contig, deterministic.
    ordered = sorted(groups.values(), key=lambda g: (-max(lengths[c] for c in g), min(g)))
    return {c: f"votu_{i:03d}" for i, g in enumerate(ordered, 1) for c in g}


def alpha(depth: np.ndarray, reads: np.ndarray) -> dict:
    present = depth > 0
    p = depth[present] / depth.sum() if present.any() else np.array([])
    shannon = float(-(p * np.log(p)).sum()) if p.size else 0.0
    f1, f2 = int((reads == 1).sum()), int((reads == 2).sum())
    return {"observed": int(present.sum()), "shannon": shannon,
            "simpson": float(1 - (p ** 2).sum()) if p.size else 0.0,
            "pielou": shannon / math.log(p.size) if p.size > 1 else float("nan"),
            "chao1_reads": int(present.sum()) + f1 * (f1 - 1) / (2 * (f2 + 1))}


def beta(table: np.ndarray) -> dict[str, np.ndarray]:
    """Distance matrices between samples (rows of `table`: samples x vOTUs)."""
    rel = table / table.sum(axis=1, keepdims=True)
    pseudo = table[table > 0].min() / 2
    clr = np.log(table + pseudo)
    clr -= clr.mean(axis=1, keepdims=True)
    return {"bray_curtis": squareform(pdist(rel, "braycurtis")),
            "jaccard": squareform(pdist(table > 0, "jaccard")),
            "aitchison": squareform(pdist(clr, "euclidean"))}


def pcoa(dist: np.ndarray, axes: int = 2) -> tuple[np.ndarray, np.ndarray]:
    n = len(dist)
    j = np.eye(n) - 1 / n
    vals, vecs = np.linalg.eigh(-0.5 * j @ (dist ** 2) @ j)
    order = np.argsort(vals)[::-1]
    vals, vecs = vals[order], vecs[:, order]
    pos = np.clip(vals, 0, None)
    k = min(axes, n)
    return vecs[:, :k] * np.sqrt(pos[:k]), pos[:k] / pos.sum() if pos.sum() else pos[:k]


def _write_matrix(path: Path, names: list[str], m: np.ndarray) -> None:
    with open(path, "w") as fh:
        fh.write("\t" + "\t".join(names) + "\n")
        for n, row in zip(names, m):
            fh.write(n + "\t" + "\t".join(f"{v:.4f}" for v in row) + "\n")


def identified_contigs(sampledir: Path) -> set[str]:
    """Contigs with a protein-level tier or a nucleotide-level species call."""
    ok = set()
    for step, file, col, bad in (("classify", "taxonomy.tsv", "tier", {"unclassified"}),
                                 ("species", "species.tsv", "call", {"no_reference_match"}),
                                 ("satellite", "satellites.tsv", "call", {"protein_only_candidate"})):
        path = sampledir / step / file
        if path.exists():
            with open(path) as fh:
                ok |= {r["contig"] for r in csv.DictReader(fh, delimiter="\t") if r[col] not in bad}
    return ok


def write_contamination(path: Path, samples, blanks, real, ids, votu, contigs, depth, depth_t,
                        ratio: float) -> None:
    """Contigs whose vOTU is also in a blank: flag when sample depth < ratio x blank depth."""
    blank_depth = depth_t[blanks].sum(axis=0)
    with open(path, "w") as fh:
        fh.write("contig\tsample\tvotu\tcontig_depth\tvotu_blank_depth\tdepth_ratio\tflag\n")
        for c, v in sorted(votu.items()):
            s = contigs[c][0]
            b = blank_depth[ids.index(v)]
            if s in {samples[i] for i in blanks} or b == 0:
                continue
            sample_ratio = depth_t[samples.index(s), ids.index(v)] / b
            fh.write(f"{c}\t{s}\t{v}\t{depth[c]:.2f}\t{b:.2f}\t{sample_ratio:.1f}\t"
                     f"{'possible_contamination' if sample_ratio < ratio else 'ok'}\n")


def run_diversity(cfg: Config, outdir: Path, samples: list[str] | None, threads: int,
                  force: bool = False) -> None:
    out = outdir / "cohort" / "diversity"
    if step_done(out) and not force:
        event("diversity skipped (done)")
        return
    out.mkdir(parents=True, exist_ok=True)
    samples = samples or sorted(p.parent.parent.name for p in outdir.glob("*/abundance/.done"))
    if not samples:
        raise ToolError(f"no finished abundance steps found under {outdir}")

    contigs, lengths, depth, reads = {}, {}, {}, {}
    for s in samples:
        seqs = read_fasta(outdir / s / "assemble" / "contigs.fasta")
        keep = identified_contigs(outdir / s) if cfg.diversity.exclude_unclassified else None
        with open(outdir / s / "abundance" / "abundance.tsv") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if (r["detected"] == "True" and int(r["length"]) >= cfg.diversity.min_contig_length
                        and (keep is None or r["contig"] in keep)):
                    contigs[r["contig"]] = (s, seqs[r["contig"]])
                    lengths[r["contig"]] = int(r["length"])
                    depth[r["contig"]], reads[r["contig"]] = float(r["mean_depth"]), int(r["reads"])
    if not contigs:
        raise ToolError("no detected contigs in any sample")
    allfa = out / "all_contigs.fasta"
    allfa.write_text("".join(f">{c}\n{seq}\n" for c, (_, seq) in contigs.items()))
    paf = out / "all_vs_all.paf"
    run(["minimap2", "-x", "asm20", "-c", "-X", "-t", threads, allfa, allfa], stdout=paf)
    d = cfg.diversity
    votu = cluster(paf, lengths, d.cluster_identity, d.cluster_min_cover)
    paf.unlink()

    ids = sorted(set(votu.values()))
    col = {s: i for i, s in enumerate(samples)}
    depth_t = np.zeros((len(samples), len(ids)))
    reads_t = np.zeros_like(depth_t)
    for c, v in votu.items():
        i, j = col[contigs[c][0]], ids.index(v)
        depth_t[i, j] += depth[c]
        reads_t[i, j] += reads[c]

    rep = {v: max((c for c in votu if votu[c] == v), key=lambda c: lengths[c]) for v in ids}
    tax = {}
    for s in samples:
        tf = outdir / s / "classify" / "taxonomy.tsv"
        if tf.exists():
            with open(tf) as fh:
                tax.update({r["contig"]: r for r in csv.DictReader(fh, delimiter="\t")})
    with open(out / "votu_table.tsv", "w") as fh:
        fh.write("votu\trepresentative\tlength\tcontigs\ttaxon\ttier\tlineage\t" + "\t".join(samples) + "\n")
        for j, v in enumerate(ids):
            t = tax.get(rep[v], {})
            fh.write(f"{v}\t{rep[v]}\t{lengths[rep[v]]}\t{sum(1 for x in votu.values() if x == v)}\t"
                     f"{t.get('name', 'unclassified')}\t{t.get('tier', 'unclassified')}\t{t.get('lineage', '')}\t"
                     + "\t".join(f"{depth_t[i, j]:.2f}" for i in range(len(samples))) + "\n")
    with open(out / "votus.tsv", "w") as fh:
        fh.write("contig\tvotu\tsample\n")
        fh.writelines(f"{c}\t{v}\t{contigs[c][0]}\n" for c, v in sorted(votu.items()))
    allfa.unlink()

    blanks = [i for i, s in enumerate(samples) if s in d.blanks]
    real = [i for i in range(len(samples)) if i not in blanks]     # negative controls are not analysed
    if blanks:
        write_contamination(out / "contamination.tsv", samples, blanks, real, ids, votu, contigs,
                            depth, depth_t, d.blank_ratio)
    rows = [alpha(depth_t[i], reads_t[i]) for i in real]
    with open(out / "alpha.tsv", "w") as fh:
        keys = list(rows[0])
        fh.write("sample\t" + "\t".join(keys) + "\n")
        for s, r in zip([samples[i] for i in real], rows):
            fh.write(s + "\t" + "\t".join("NA" if isinstance(r[k], float) and math.isnan(r[k])
                                          else f"{r[k]:.4g}" for k in keys) + "\n")

    keep = [i for i in real if depth_t[i].sum() > 0]
    names = [samples[i] for i in keep]
    if len(keep) >= 2:
        dists = beta(depth_t[keep])
        for name, m in dists.items():
            _write_matrix(out / f"beta_{name}.tsv", names, m)
        if len(keep) >= 3:
            coords, expl = pcoa(dists["bray_curtis"])
            with open(out / "pcoa_bray_curtis.tsv", "w") as fh:
                fh.write("sample\t" + "\t".join(f"PC{k + 1}({e:.1%})" for k, e in enumerate(expl)) + "\n")
                for n, row in zip(names, coords):
                    fh.write(n + "\t" + "\t".join(f"{v:.4f}" for v in row) + "\n")
    else:
        event("fewer than 2 samples with detections; beta diversity skipped")
    event("diversity done", samples=len(samples), votus=len(ids))
    mark_done(out)
