"""Step 9 (optional): host variety markers from organelle reads.

All QC reads (host and non-host) are mapped to the references in
host.variety_references (e.g. chloroplast, mitochondrion, conserved genes).
SNVs are called from minimap2 cs strings: no extra tools, but a deliberately
simple caller (see limits). Coverage is ALWAYS reported, and flagged PARTIAL
unless most of each reference reaches min_depth: a variety call needs the
informative positions covered, not just a few variants.

Limits: SNVs only (Nanopore indels are too error-prone to call here); organelles
are treated as homoplasmic, so heteroplasmy is not assessed; copies of one RCA
concatemer come from a single molecule and are not independent evidence; identifying
a variety further needs a panel of known cultivar haplotypes, not provided here.
"""
import re
from collections import defaultdict
from pathlib import Path

from ..config import Config
from ..logging import event
from ..runner import mark_done, run, step_done
from .assemble import read_fasta

CS = re.compile(r":(\d+)|\*([a-z])([a-z])|\+([a-z]+)|-([a-z]+)")


def pileup(paf: Path, min_identity: float, min_mapq: int):
    """Per reference: depth diff array and {pos: {base: count}} for mismatches."""
    depth_diff: dict[str, dict[int, int]] = defaultdict(lambda: defaultdict(int))
    alts: dict[str, dict[int, dict[str, int]]] = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    with open(paf) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if int(f[10]) == 0 or int(f[9]) / int(f[10]) < min_identity or int(f[11]) < min_mapq:
                continue
            cs = next((x[5:] for x in f[12:] if x.startswith("cs:Z:")), None)
            if cs is None:
                continue
            ref, pos = f[5], int(f[7])
            depth_diff[ref][pos] += 1
            depth_diff[ref][int(f[8])] -= 1
            for m in CS.finditer(cs):
                if m.group(1):
                    pos += int(m.group(1))
                elif m.group(2):
                    alts[ref][pos][m.group(3).upper()] += 1
                    pos += 1
                elif m.group(5):
                    pos += len(m.group(5))
    return depth_diff, alts


def depth_array(diff: dict[int, int], length: int) -> list[int]:
    out, cur = [], 0
    for i in range(length):
        cur += diff.get(i, 0)
        out.append(cur)
    return out


def run_variety(cfg: Config, sampledir: Path, threads: int, force: bool = False) -> None:
    out = sampledir / "variety"
    if step_done(out) and not force:
        event("variety skipped (done)", sample=sampledir.name)
        return
    out.mkdir(parents=True, exist_ok=True)
    refs = cfg.host.variety_references
    if not refs:
        event("no host.variety_references configured; variety step skipped", sample=sampledir.name)
        mark_done(out)
        return
    v = cfg.variety
    panel = out / "refs.fa"
    panel.write_text("".join(f">{n}\n{s}\n" for r in refs for n, s in read_fasta(r).items()))
    seqs = read_fasta(panel)
    paf = out / "reads_vs_refs.paf"
    run(["minimap2", "-c", "--cs", "-x", "map-ont", "-t", threads, "--secondary=no", panel,
         sampledir / "qc" / "reads.qc.fastq.gz"], stdout=paf)
    depth_diff, alts = pileup(paf, v.min_identity, v.min_mapq)
    paf.unlink()

    cov_rows, var_rows, consensus = [], [], []
    for name, ref in seqs.items():
        depth = depth_array(depth_diff.get(name, {}), len(ref))
        cons = [c if d >= v.min_depth else "N" for c, d in zip(ref, depth)]
        n_var = 0
        for pos, counts in sorted(alts.get(name, {}).items()):
            alt, n = max(counts.items(), key=lambda kv: kv[1])
            if depth[pos] >= v.min_depth and n / depth[pos] >= v.min_alt_frac:
                var_rows.append(f"{name}\t{pos + 1}\t{ref[pos]}\t{alt}\t{depth[pos]}\t{n}\t{n / depth[pos]:.2f}\n")
                cons[pos] = alt
                n_var += 1
        breadth = sum(d >= v.min_depth for d in depth) / len(ref)
        cov_rows.append(f"{name}\t{len(ref)}\t{sum(depth) / len(ref):.1f}\t"
                        f"{sum(d > 0 for d in depth) / len(ref):.3f}\t{breadth:.3f}\t{n_var}\t"
                        f"{'OK' if breadth >= v.complete_breadth else 'PARTIAL'}\n")
        consensus.append(f">{name}\n{''.join(cons)}\n")
    (out / "coverage.tsv").write_text("reference\tlength\tmean_depth\tbreadth_1x\tbreadth_min_depth\t"
                                      "variants\tcoverage\n" + "".join(cov_rows))
    (out / "variants.tsv").write_text("reference\tpos\tref\talt\tdepth\talt_reads\talt_frac\n"
                                      + "".join(var_rows))
    (out / "consensus.fasta").write_text("".join(consensus))
    panel.unlink()
    event("variety done", sample=sampledir.name, refs=len(seqs), variants=len(var_rows))
    mark_done(out)
