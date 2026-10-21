"""Step 6: read mapping to contigs -> depth, breadth, relative abundance, chimeras.

Abundance is mean depth (aligned bases / contig length), not read counts: a
concatemer read carries many genome copies, so counting reads would depend on
read length. RCA amplification is uneven and favours small circular DNA, so
relative abundances are semi-quantitative at best; depth CV is reported to show
how uneven coverage is.
"""
import statistics
from collections import defaultdict
from pathlib import Path

from ..config import Config
from ..logging import event
from ..runner import mark_done, run, step_done
from .assemble import read_fasta


def parse_paf(path: Path, min_identity: float) -> dict[str, list[tuple]]:
    """query -> [(qlen, qstart, qend, contig, rstart, rend)] for alignments passing identity."""
    out: dict[str, list[tuple]] = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            f = line.split("\t")
            if int(f[10]) and int(f[9]) / int(f[10]) >= min_identity:
                out[f[0]].append((int(f[1]), int(f[2]), int(f[3]), f[5], int(f[7]), int(f[8])))
    return out


def assign(alns: dict[str, list[tuple]], min_cover: float) -> dict[str, list[tuple]]:
    """Assign each read to the contig with most aligned bases (if >= min_cover of the read)."""
    assigned: dict[str, list[tuple]] = {}
    for read, rows in alns.items():
        per_contig: dict[str, int] = defaultdict(int)
        for qlen, qs, qe, contig, *_ in rows:
            per_contig[contig] += qe - qs
        best = max(per_contig, key=per_contig.get)
        if per_contig[best] / rows[0][0] >= min_cover:
            assigned[read] = [r for r in rows if r[3] == best]
    return assigned


def chimeras(alns: dict[str, list[tuple]], min_frac: float) -> dict[tuple[str, str], int]:
    """Reads with >= min_frac of their length on each of two different contigs."""
    pairs: dict[tuple[str, str], int] = defaultdict(int)
    for rows in alns.values():
        per_contig: dict[str, int] = defaultdict(int)
        for qlen, qs, qe, contig, *_ in rows:
            per_contig[contig] += qe - qs
        big = sorted(c for c, n in per_contig.items() if n / rows[0][0] >= min_frac)
        for i, a in enumerate(big):
            for b in big[i + 1:]:
                pairs[(a, b)] += 1
    return pairs


def depth_profile(rows: list[tuple], length: int) -> list[int]:
    diff = [0] * (length + 1)
    for *_, rs, re_ in rows:
        diff[rs] += 1
        diff[min(re_, length)] -= 1
    depth, cur = [], 0
    for d in diff[:length]:
        cur += d
        depth.append(cur)
    return depth


def run_abundance(cfg: Config, sampledir: Path, threads: int, force: bool = False) -> None:
    out = sampledir / "abundance"
    if step_done(out) and not force:
        event("abundance skipped (done)", sample=sampledir.name)
        return
    out.mkdir(parents=True, exist_ok=True)
    contigs = sampledir / "assemble" / "contigs.fasta"
    reads = sampledir / "hostdeplete" / "reads.nonhost.fastq.gz"
    a = cfg.abundance
    lengths = {k: len(v) for k, v in read_fasta(contigs).items()}
    header = "contig\tlength\treads\tmean_depth\tdepth_cv\tbreadth\tdetected\trel_abundance\tchimeric_reads\n"
    alns: dict[str, list[tuple]] = {}
    if lengths:
        paf = out / "reads_vs_contigs.paf"
        run(["minimap2", "-c", "-x", "map-ont", "-t", threads, "--secondary=no", contigs, reads],
            stdout=paf)
        alns = parse_paf(paf, a.min_identity)
        paf.unlink()
    assigned = assign(alns, a.min_read_cover)
    pairs = chimeras(alns, a.chimera_min_frac)
    chim_count: dict[str, int] = defaultdict(int)
    for (x, y), n in pairs.items():
        chim_count[x] += n
        chim_count[y] += n

    per_contig = defaultdict(list)
    for rows in assigned.values():
        per_contig[rows[0][3]] += rows
    stats = {}
    for contig, length in lengths.items():
        depth = depth_profile(per_contig.get(contig, []), length)
        mean = sum(depth) / length
        stats[contig] = {
            "reads": sum(1 for r in assigned.values() if r[0][3] == contig),
            "mean": mean, "breadth": sum(d > 0 for d in depth) / length,
            "cv": statistics.pstdev(depth) / mean if mean else 0.0}
    total = sum(s["mean"] for s in stats.values() if s["breadth"] >= a.min_breadth)
    with open(out / "abundance.tsv", "w") as fh:
        fh.write(header)
        for contig, s in stats.items():
            det = s["breadth"] >= a.min_breadth
            fh.write(f"{contig}\t{lengths[contig]}\t{s['reads']}\t{s['mean']:.2f}\t{s['cv']:.2f}\t"
                     f"{s['breadth']:.3f}\t{det}\t{s['mean'] / total if det and total else 0:.4f}\t"
                     f"{chim_count.get(contig, 0)}\n")
    with open(out / "chimeras.tsv", "w") as fh:
        fh.write("contig_a\tcontig_b\treads\n")
        for (x, y), n in sorted(pairs.items()):
            fh.write(f"{x}\t{y}\t{n}\n")
    event("abundance done", sample=sampledir.name, contigs=len(lengths), mapped_reads=len(assigned))
    mark_done(out)
