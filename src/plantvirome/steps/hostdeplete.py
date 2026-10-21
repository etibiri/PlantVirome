"""Step 2: soft host depletion against a (possibly partial) reference panel.

Reads are split into host / non-host files, never discarded. With no host
reference configured, all reads pass through as non-host.
"""
import gzip
import json
import shutil
from pathlib import Path

from ..config import Config
from ..logging import event
from ..runner import mark_done, run, step_done


def host_read_ids(paf: Path, min_cov: float, min_identity: float) -> set[str]:
    """Reads whose union of aligned query intervals covers >= min_cov of the read."""
    spans: dict[str, tuple[int, list[tuple[int, int]]]] = {}
    with open(paf) as fh:
        for line in fh:
            f = line.split("\t")
            matches, block = int(f[9]), int(f[10])
            if block == 0 or matches / block < min_identity:
                continue
            spans.setdefault(f[0], (int(f[1]), []))[1].append((int(f[2]), int(f[3])))
    hosts = set()
    for name, (qlen, ivs) in spans.items():
        covered, end = 0, 0
        for s, e in sorted(ivs):
            covered += max(0, e - max(s, end))
            end = max(end, e)
        if covered / qlen >= min_cov:
            hosts.add(name)
    return hosts


def split_fastq(src: Path, hosts: set[str], nonhost: Path, host: Path) -> tuple[int, int]:
    counts = [0, 0]
    with gzip.open(src, "rt") as fh, gzip.open(nonhost, "wt") as nh, gzip.open(host, "wt") as h:
        while header := fh.readline():
            rec = header + fh.readline() + fh.readline() + fh.readline()
            is_host = header[1:].split()[0] in hosts
            (h if is_host else nh).write(rec)
            counts[is_host] += 1
    return counts[0], counts[1]


def run_hostdeplete(cfg: Config, sampledir: Path, threads: int, force: bool = False) -> None:
    out = sampledir / "hostdeplete"
    if step_done(out) and not force:
        event("hostdeplete skipped (done)", sample=sampledir.name)
        return
    out.mkdir(parents=True, exist_ok=True)
    reads = sampledir / "qc" / "reads.qc.fastq.gz"
    nonhost, host = out / "reads.nonhost.fastq.gz", out / "reads.host.fastq.gz"
    refs = cfg.host.references
    if refs:
        panel = out / "host_panel.fa"
        with open(panel, "wb") as w:
            for r in refs:
                with open(r, "rb") as fh:
                    shutil.copyfileobj(fh, w)
        paf = out / "host.paf"
        run(["minimap2", "-c", "-x", "map-ont", "-t", threads, "--secondary=no", panel, reads],
            stdout=paf)
        hosts = host_read_ids(paf, cfg.hostdeplete.min_query_cov, cfg.hostdeplete.min_identity)
        n_nonhost, n_host = split_fastq(reads, hosts, nonhost, host)
        paf.unlink()
        panel.unlink()
    else:
        event("no host reference configured; passing all reads through",
              sample=sampledir.name)
        shutil.copyfile(reads, nonhost)
        with gzip.open(host, "wt"):
            pass
        with gzip.open(nonhost, "rt") as fh:
            n_nonhost, n_host = sum(1 for _ in fh) // 4, 0
    stats = {"host": n_host, "nonhost": n_nonhost,
             "host_fraction": round(n_host / max(1, n_host + n_nonhost), 3)}
    (out / "stats.json").write_text(json.dumps(stats, indent=2))
    event("hostdeplete stats", sample=sampledir.name, **stats)
    mark_done(out)
