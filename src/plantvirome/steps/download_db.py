"""Build the `classify` databases from NCBI RefSeq viral proteins and the NCBI taxdump.

Downloads viral.*.protein.faa.gz (sequences) and the matching .gpff.gz (the only
RefSeq viral file that carries the taxid of each protein), plus taxdump.tar.gz,
then runs `diamond makedb` with taxonomy. Provenance goes to download.json.
"""
import gzip
import json
import re
import shutil
import tarfile
import time
import urllib.request
from pathlib import Path

from ..logging import event
from ..runner import run

BASE_URL = "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/"
TAXDUMP_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
TAXDUMP_FILES = ("nodes.dmp", "names.dmp", "merged.dmp", "delnodes.dmp")
PART = re.compile(r'href="(viral(?:\.wp_protein)?\.\d+\.protein)\.faa\.gz"')


def discover_parts(index_html: str) -> list[str]:
    """Part names (without extension) of protein FASTA files listed in a directory index."""
    return sorted(set(PART.findall(index_html)))


def fetch(url: str, dest: Path) -> None:
    if dest.exists():
        event("already downloaded", file=dest.name)
        return
    tmp = dest.with_suffix(dest.suffix + ".part")
    event("downloading", url=url)
    with urllib.request.urlopen(url) as r, open(tmp, "wb") as w:
        shutil.copyfileobj(r, w, 1 << 20)
    tmp.rename(dest)


def gpff_taxids(path: Path):
    """Yield (protein accession.version, taxid) from a RefSeq .gpff(.gz) file."""
    version = None
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("VERSION"):
                version = line.split()[1]
            elif version and "/db_xref=\"taxon:" in line:
                yield version, line.split("taxon:")[1].split('"')[0]
                version = None


def build_db(faa_files: list[Path], gpff_files: list[Path], taxdump_tgz: Path,
             outdir: Path, threads: int) -> Path:
    """Write outdir/viral_protein.dmnd and outdir/taxdump/. Returns the .dmnd path."""
    tax = outdir / "taxdump"
    tax.mkdir(parents=True, exist_ok=True)
    with tarfile.open(taxdump_tgz) as tf:
        for name in TAXDUMP_FILES:
            with tf.extractfile(name) as src, open(tax / name, "wb") as dst:
                shutil.copyfileobj(src, dst)
    mapping = outdir / "prot.accession2taxid"
    n = 0
    with open(mapping, "w") as w:
        w.write("accession\taccession.version\ttaxid\tgi\n")
        for gpff in gpff_files:
            for version, taxid in gpff_taxids(gpff):
                w.write(f"{version.split('.')[0]}\t{version}\t{taxid}\t0\n")
                n += 1
    faa = outdir / "viral.faa.gz"
    with open(faa, "wb") as w:  # gzip members concatenate validly
        for f in faa_files:
            with open(f, "rb") as fh:
                shutil.copyfileobj(fh, w)
    run(["diamond", "makedb", "--in", faa, "-d", outdir / "viral_protein", "--threads", threads,
         "--taxonmap", mapping, "--taxonnodes", tax / "nodes.dmp", "--taxonnames", tax / "names.dmp"])
    mapping.unlink()
    faa.unlink()
    event("database built", proteins_with_taxid=n)
    return outdir / "viral_protein.dmnd"


def run_download_db(outdir: Path, threads: int, base_url: str = BASE_URL,
                    taxdump_url: str = TAXDUMP_URL) -> Path:
    raw = outdir / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(base_url) as r:
        parts = discover_parts(r.read().decode())
    if not parts:
        raise RuntimeError(f"no viral protein files found at {base_url}")
    faa, gpff = [], []
    for p in parts:
        for ext, bucket in (("faa.gz", faa), ("gpff.gz", gpff)):
            fetch(f"{base_url}{p}.{ext}", raw / f"{p}.{ext}")
            bucket.append(raw / f"{p}.{ext}")
    fetch(taxdump_url, raw / "taxdump.tar.gz")
    dmnd = build_db(faa, gpff, raw / "taxdump.tar.gz", outdir, threads)
    (outdir / "download.json").write_text(json.dumps(
        {"date": time.strftime("%Y-%m-%d"), "base_url": base_url, "taxdump_url": taxdump_url,
         "parts": parts}, indent=2))
    return dmnd
