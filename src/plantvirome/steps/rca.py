"""Step 3: RCA concatemer handling.

TideHunter finds tandem repeats in each read and calls an error-corrected
consensus of the repeat unit (one genome copy). Reads are then classified:
  clean      one repeat region covers >= min_read_cover of the read
  multi      >= 2 repeat regions in one read: putative chimera (template switching)
  partial    one region covering less than min_read_cover (unexplained flank)
  untandem   no concatemer found (single-copy or short read); kept aside, not used
Only "clean" units feed the assembler. "multi" and "partial" reads are putative
chimeras (TideHunter often reports just one of the two regions of a chimera, which
shows up as "partial"); their units go to units.flagged.fasta so they can be
inspected or rescued, never silently merged.
"""
import gzip
import json
from collections import defaultdict
from pathlib import Path

from ..config import Config
from ..logging import event
from ..runner import mark_done, run, step_done


def parse_tidehunter(path: Path) -> dict[str, list[dict]]:
    """Tabular TideHunter output (-f 2): name, rep, copies, readlen, start, end,
    unitlen, identity, fullcopies, unit_starts, consensus."""
    regions: dict[str, list[dict]] = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            regions[f[0]].append({"copies": float(f[2]), "readlen": int(f[3]),
                                  "start": int(f[4]), "end": int(f[5]),
                                  "unitlen": int(f[6]), "identity": float(f[7]),
                                  "seq": f[10]})
    return regions


def classify(regs: list[dict], min_cover: float) -> str:
    if not regs:
        return "untandem"
    if len(regs) > 1:
        return "multi"
    r = regs[0]
    return "clean" if (r["end"] - r["start"]) / r["readlen"] >= min_cover else "partial"


def run_rca(cfg: Config, sampledir: Path, threads: int, force: bool = False) -> None:
    out = sampledir / "rca"
    if step_done(out) and not force:
        event("rca skipped (done)", sample=sampledir.name)
        return
    out.mkdir(parents=True, exist_ok=True)
    reads = sampledir / "hostdeplete" / "reads.nonhost.fastq.gz"
    tab = out / "tidehunter.tsv"
    c = cfg.rca
    run(["TideHunter", "-f", 2, "-t", threads, "-c", c.min_copies, "-p", c.min_period,
         "-P", c.max_period, reads], stdout=tab)
    regions = parse_tidehunter(tab)

    classes: dict[str, str] = {}
    with gzip.open(reads, "rt") as fh:  # every input read gets a class
        for i, line in enumerate(fh):
            if i % 4 == 0:
                name = line[1:].split()[0]
                classes[name] = classify(regions.get(name, []), c.min_read_cover)

    with open(out / "units.fasta", "w") as ok, open(out / "units.flagged.fasta", "w") as flagged, \
            open(out / "reads.tsv", "w") as tsv:
        tsv.write("read\tclass\tregions\tcopies\tunit_length\tidentity\n")
        for name, cls in classes.items():
            regs = regions.get(name, [])
            for k, r in enumerate(regs):
                if cls != "untandem":
                    fa = ok if cls == "clean" else flagged
                    fa.write(f">{name}|r{k}|copies={r['copies']:g}|id={r['identity']}\n{r['seq']}\n")
            tsv.write(f"{name}\t{cls}\t{len(regs)}\t"
                      f"{','.join(f'{r['copies']:g}' for r in regs)}\t"
                      f"{','.join(str(r['unitlen']) for r in regs)}\t"
                      f"{','.join(str(r['identity']) for r in regs)}\n")
    counts = {k: sum(v == k for v in classes.values())
              for k in ("clean", "multi", "partial", "untandem")}
    (out / "stats.json").write_text(json.dumps(counts, indent=2))
    tab.unlink()
    event("rca stats", sample=sampledir.name, **counts)
    mark_done(out)
