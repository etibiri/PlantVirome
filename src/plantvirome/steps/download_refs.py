"""Fetch RefSeq reference genomes of a taxon from NCBI Entrez for the `species` step.

RefSeq holds one exemplar per species, which is what species demarcation compares
against. A contig that is >= threshold to a non-exemplar member but below it to the
exemplar would be called a novel candidate: confirm borderline calls.
"""
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

from ..logging import event

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
BATCH = 100


B_TITLE = re.compile(r"\bDNA[ _-]?B\b|\bDNA[ _-]?2\b|\b(?:segment|component)[ _-]?B\b|\bBR1\b|\bBL1\b",
                     re.IGNORECASE)


def component(title: str) -> str:
    """A = DNA-A or a monopartite genome (used for species demarcation); B = DNA-B.
    Derived from the RefSeq title, which names the second component inconsistently
    (DNA-B, DNA B, DNA_B, DNA2, segment B, gene names): records with no hint count as A,
    so a few unlabelled DNA-B references may still be taken for A."""
    if B_TITLE.search(title):
        return "B"
    if "satellite" in title.lower():
        return "satellite"
    return "A"


def _get(endpoint: str, **params) -> str:
    url = EUTILS + endpoint + "?" + urllib.parse.urlencode({"tool": "plantvirome", **params})
    time.sleep(0.4)  # NCBI allows 3 requests/second without an API key
    with urllib.request.urlopen(url, timeout=120) as r:
        return r.read().decode()


def fasta_records(text: str) -> dict[str, str]:
    recs: dict[str, list[str]] = {}
    for line in text.splitlines():
        if line.startswith(">"):
            acc = line[1:].split()[0]
            recs[acc] = []
        elif line.strip():
            recs[acc].append(line.strip())
    return {k: "".join(v) for k, v in recs.items()}


def genbank_cds(text: str) -> list[tuple[str, str, str, str]]:
    """(accession.version, gene, product, protein) of every CDS in GenBank flat-file text.
    The gene falls back to the product when the record has no /gene."""
    out = []
    for rec in text.split("\n//\n"):
        m = re.search(r"^VERSION\s+(\S+)", rec, re.M)
        feats = rec.partition("\nFEATURES")[2].partition("\nORIGIN")[0]
        for blk in re.split(r"\n {5}(?=\S)", feats)[1:]:
            if not blk.startswith("CDS") or not m:
                continue
            q = {k: re.sub(r"\s+", "" if k == "translation" else " ", v).strip()
                 for k, v in re.findall(r'/(\w+)="([^"]*)"', blk)}
            if q.get("translation"):
                out.append((m.group(1), q.get("gene", ""), q.get("product", ""), q["translation"]))
    return out


def run_download_proteins(refs_tsv: Path) -> Path:
    """Write refs.faa (CDS proteins of every accession of refs.tsv) next to it, for `annotate`.
    Headers: >ACC|gene|product_with_underscores."""
    with open(refs_tsv) as fh:
        acc = [line.split("\t")[0] for line in fh.read().splitlines()[1:]]
    out = refs_tsv.with_name("refs.faa")
    n = 0
    with open(out, "w") as fa:
        for i in range(0, len(acc), BATCH):
            for a, gene, product, prot in genbank_cds(_get("efetch.fcgi", db="nuccore", rettype="gb",
                                                          retmode="text", id=",".join(acc[i:i + BATCH]))):
                name = re.sub(r"[|\s]+", "_", product)
                fa.write(f">{a}|{re.sub(r'[|\s]+', '_', gene or product)}|{name}\n{prot}\n")
                n += 1
    event("reference proteins written", n=n, path=str(out))
    return out


def run_download_refs(taxon: str, outdir: Path, extra: dict[str, int | None] | None = None) -> Path:
    """`extra`: non-RefSeq accessions to add (e.g. isolates from INSDC) -> species taxid to
    label them with (None = the record's own organism, often a strain-level taxon)."""
    outdir.mkdir(parents=True, exist_ok=True)
    term = f'"{taxon}"[Organism:exp] AND srcdb_refseq[PROP] AND biomol_genomic[PROP]'
    ids = json.loads(_get("esearch.fcgi", db="nuccore", term=term, retmax=10000,
                          retmode="json"))["esearchresult"]["idlist"]
    if not ids:
        raise RuntimeError(f"no RefSeq genomes found for {taxon!r}")
    event("reference genomes found", taxon=taxon, n=len(ids))
    rows = []
    with open(outdir / "refs.fasta", "w") as fa:
        for i in range(0, len(ids), BATCH):
            batch = ",".join(ids[i:i + BATCH])
            summ = json.loads(_get("esummary.fcgi", db="nuccore", id=batch, retmode="json"))["result"]
            seqs = fasta_records(_get("efetch.fcgi", db="nuccore", id=batch, rettype="fasta",
                                      retmode="text"))
            for uid in summ["uids"]:
                r = summ[uid]
                acc = r["accessionversion"]
                species = re.sub(r"[|\s]+", " ", r["organism"]).strip()
                comp = component(r["title"])
                # no spaces in the header: minimap2 and FASTA readers keep only the first word
                fa.write(f">{acc}|{species.replace(' ', '_')}|{r['taxid']}|{comp}\n{seqs[acc]}\n")
                rows.append((acc, species, r["taxid"], comp, r["slen"], r["title"]))
        for acc, taxid in (extra or {}).items():
            summ = json.loads(_get("esummary.fcgi", db="nuccore", id=acc, retmode="json"))["result"]
            r = summ[summ["uids"][0]]
            acc = r["accessionversion"]
            name, tid = r["organism"], r["taxid"]
            if taxid:
                tax = json.loads(_get("esummary.fcgi", db="taxonomy", id=taxid, retmode="json"))["result"]
                name, tid = tax[str(taxid)]["scientificname"], taxid
            seq = fasta_records(_get("efetch.fcgi", db="nuccore", id=acc, rettype="fasta",
                                     retmode="text"))[acc]
            species, comp = re.sub(r"[|\s]+", " ", name).strip(), component(r["title"])
            fa.write(f">{acc}|{species.replace(' ', '_')}|{tid}|{comp}\n{seq}\n")
            rows.append((acc, species, tid, comp, r["slen"], r["title"]))
    with open(outdir / "refs.tsv", "w") as fh:
        fh.write("accession\tspecies\ttaxid\tcomponent\tlength\ttitle\n")
        fh.writelines("\t".join(map(str, r)) + "\n" for r in rows)
    (outdir / "download.json").write_text(json.dumps(
        {"date": time.strftime("%Y-%m-%d"), "taxon": taxon, "query": term, "extra": list(extra or {}),
         "records": len(rows)}, indent=2))
    run_download_proteins(outdir / "refs.tsv")
    return outdir / "refs.fasta"
