"""Step 4: assemble RCA unit consensus sequences into circular genome contigs.

Every unit is one genome copy of a circular template, in random rotation and
strand. Generic long-read assemblers (Flye was tried) behave erratically on such
fixed-length, fully overlapping reads, so units are assembled directly:
  1. take a seed unit and map all remaining units to seed+seed (absorbs rotation),
  2. units aligning end to end (>= min_identity, >= min_cover) join the cluster,
  3. a per-position majority vote (substitutions, deletions, insertions) over the
     members gives the contig, in the seed's coordinates,
  4. repeat on the units that did not join.
Limits: closely related strains (above min_identity) collapse into one majority
consensus; divergent ones (below it) become separate contigs; a template
that is not a full circle (no concatemer) cannot be assembled here.
"""
import re
from collections import defaultdict
from pathlib import Path

from ..config import Assemble, Config
from ..logging import event
from ..runner import mark_done, run, step_done

CS = re.compile(r":(\d+)|\*([a-z])([a-z])|\+([a-z]+)|-([a-z]+)")


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif line:
            seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def consensus(seed: str, alns: list[tuple[int, int, str]]) -> str:
    """Majority-vote polish of `seed` from (ref_start, ref_end, cs) alignments to seed+seed."""
    n = len(seed)
    depth = [0] * n
    subs = [defaultdict(int) for _ in range(n)]
    dels = [0] * n
    ins = [defaultdict(int) for _ in range(n)]
    for rs, re_, cs in alns:
        for p in range(rs, re_):
            depth[p % n] += 1
        pos = rs
        for m in CS.finditer(cs):
            if m.group(1):
                pos += int(m.group(1))
            elif m.group(2):
                subs[pos % n][m.group(3).upper()] += 1
                pos += 1
            elif m.group(4):
                ins[(pos - 1) % n][m.group(4).upper()] += 1   # inserted after previous base
            else:
                for j in range(len(m.group(5))):
                    dels[(pos + j) % n] += 1
                pos += len(m.group(5))
    out = []
    for p in range(n):
        d = depth[p]
        if d and dels[p] > d / 2:
            pass
        else:
            counts = dict(subs[p])
            counts[seed[p]] = counts.get(seed[p], 0) + d - sum(subs[p].values()) - dels[p]
            out.append(max(counts, key=counts.get))
        if ins[p]:
            best, k = max(ins[p].items(), key=lambda kv: kv[1])
            if k > d / 2:
                out.append(best)
    return "".join(out)


def best_alignments(paf: Path, a: Assemble) -> dict[str, tuple[int, int, str]]:
    """Per unit, its longest alignment if it passes identity and coverage."""
    best: dict[str, tuple[int, tuple[int, int, str]]] = {}
    with open(paf) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            span = int(f[3]) - int(f[2])
            cs = next((x[5:] for x in f[12:] if x.startswith("cs:Z:")), None)
            if cs and int(f[10]) and int(f[9]) / int(f[10]) >= a.min_identity \
                    and span / int(f[1]) >= a.min_cover and span > best.get(f[0], (0,))[0]:
                best[f[0]] = (span, (int(f[7]), int(f[8]), cs))
    return {k: v[1] for k, v in best.items()}


def cluster_units(units: dict[str, str], a: Assemble, threads: int, work: Path) -> list[tuple[str, int]]:
    """Returns [(consensus, n_units)] sorted by decreasing n_units."""
    remaining = dict(units)
    contigs = []
    ref, qry, paf = work / "seed.fa", work / "units.fa", work / "seed.paf"
    while len(remaining) >= a.min_cluster_units:
        name, seed = next(iter(remaining.items()))
        ref.write_text(f">seed\n{seed}{seed}\n")
        qry.write_text("".join(f">{k}\n{v}\n" for k, v in remaining.items()))
        run(["minimap2", "-c", "--cs", "-x", "map-ont", "-t", threads, "--secondary=no", ref, qry],
            stdout=paf)
        members = best_alignments(paf, a)
        if len(members) >= a.min_cluster_units:
            contigs.append((consensus(seed, list(members.values())), len(members)))
            for k in members:
                remaining.pop(k, None)
        remaining.pop(name, None)  # the seed leaves in any case, so the loop ends
    for f in (ref, qry, paf):
        f.unlink(missing_ok=True)
    return sorted(contigs, key=lambda c: -c[1])


def run_assemble(cfg: Config, sampledir: Path, threads: int, force: bool = False) -> None:
    out = sampledir / "assemble"
    if step_done(out) and not force:
        event("assemble skipped (done)", sample=sampledir.name)
        return
    out.mkdir(parents=True, exist_ok=True)
    a = cfg.assemble
    units = read_fasta(sampledir / "rca" / "units.fasta")
    header = "contig\tlength\tcoverage\tcircular\tcirc_support\n"
    contigs = []
    if len(units) < a.min_units:
        event("too few units; empty assembly", sample=sampledir.name, units=len(units))
    else:
        contigs = cluster_units(units, a, threads, out)
    with open(out / "contigs.fasta", "w") as fa, open(out / "contigs.tsv", "w") as tsv:
        tsv.write(header)
        for i, (seq, n) in enumerate(contigs, 1):
            name = f"{sampledir.name}_c{i}"
            fa.write(f">{name}\n{seq}\n")
            tsv.write(f"{name}\t{len(seq)}\t{n}\t{n >= a.circ_min_support}\t{n}\n")
    event("assemble stats", sample=sampledir.name, contigs=len(contigs), units=len(units))
    mark_done(out)
