"""Gene annotation and genome-completeness QC of contigs, by homology to reference proteins.

Every ATG..stop ORF of the six frames is matched (DIAMOND blastp) to the CDS of the reference
genomes (`refs.faa`, from `download-proteins`). Unlike a gene finder this keeps overlapping genes,
which compact ssDNA genomes are full of. The reference record whose proteins are best recovered
gives the expected gene set; each expected protein is `found` (subject cover >= FULL), `truncated`
(frameshift or premature stop: Nanopore indels) or `missing`. Virus-agnostic: no gene list is built in.
Proteins are identified by their normalised product name, because the gene tags of GenBank
records are inconsistent. ORFs shorter than min_aa are never called; reference proteins shorter than expect_min_aa are not
expected (small accessory ORFs such as AC4/AC5 are too variable to demand).
"""
import re
from collections import defaultdict
from pathlib import Path

from ..logging import event
from ..runner import run
from .assemble import read_fasta

FULL = 0.8      # subject cover for a gene to count as complete
COLS = ["qseqid", "sseqid", "pident", "length", "slen", "bitscore", "qstart", "sstart"]
BASES = "TCAG"
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: AA[16 * i + 4 * j + k] for i, a in enumerate(BASES) for j, b in enumerate(BASES)
         for k, c in enumerate(BASES)}
COMP = str.maketrans("ACGT", "TGCA")


def read_proteins(faa: Path) -> dict[str, dict]:
    """id (`ACC|gene|product`) -> {acc, gene, product, key, aa}."""
    out, cur = {}, None
    for line in open(faa):
        if line.startswith(">"):
            cur = line[1:].strip()
            acc, gene, product = (cur.split("|") + ["", ""])[:3]
            product = product.replace("_", " ")
            out[cur] = {"acc": acc, "gene": gene or product, "product": product, "aa": 0,
                        "key": re.sub(r"[^a-z0-9]+", " ", re.sub(r"\(.*?\)", "", product.lower())).strip()}
        elif cur:
            out[cur]["aa"] += len(line.strip())
    return out


def all_orfs(seq: str, min_aa: int) -> list[dict]:
    """ATG..stop ORFs (from the first ATG after the previous stop) of both strands, 1-based coordinates
    on `seq`; ORFs without a stop are ignored."""
    n, orfs = len(seq), []
    for strand, s in ((1, seq), (-1, seq.translate(COMP)[::-1])):
        for f in range(3):
            start = None
            prot: list[str] = []
            for i in range(f, n - 2, 3):
                aa = CODON.get(s[i:i + 3], "X")
                if aa == "*":
                    if start is not None and len(prot) >= min_aa:
                        b, e = start + 1, i + 3
                        orfs.append({"begin": b if strand == 1 else n - e + 1, "end": e if strand == 1 else n - b + 1,
                                     "strand": strand, "protein": "".join(prot),
                                     "partial_begin": False, "partial_end": False})
                    start, prot = None, []
                elif start is None and aa == "M":
                    start, prot = i, ["M"]
                elif start is not None:
                    prot.append(aa)
    return orfs


def classify_genes(ref_proteins: list[dict], best: dict[str, tuple[float, float]], min_aa: int,
                   min_identity: float) -> dict[str, list[str]]:
    """Expected proteins of one reference -> found / truncated / missing (by gene name).
    `best`: protein id -> (identity %, subject cover) of its best ORF hit. Proteins with the same
    normalised product count once, in their best state."""
    rank = {"found": 0, "truncated": 1, "missing": 2}
    state: dict[str, tuple[str, str]] = {}
    for p in ref_proteins:
        if p["aa"] < min_aa:
            continue
        ident, cover = best.get(p["id"], (0.0, 0.0))
        s = "missing" if ident < min_identity else "found" if cover >= FULL else "truncated"
        if p["key"] not in state or rank[s] < rank[state[p["key"]][0]]:
            state[p["key"]] = (s, p["gene"])
    res: dict[str, list[str]] = {"found": [], "truncated": [], "missing": []}
    for s, gene in state.values():
        res[s].append(gene)
    return res


def annotate(contigs: dict[str, str], proteins: Path, work: Path, threads: int, min_aa: int,
             min_identity: float = 35.0, min_cover: float = 0.5, expect_min_aa: int = 100) -> dict[str, dict]:
    """contig -> {orfs: [named ORFs, begin-sorted, 1-based `index`], ref, found, truncated, missing}."""
    work.mkdir(parents=True, exist_ok=True)
    prots = read_proteins(proteins)
    by_acc: dict[str, list[dict]] = defaultdict(list)
    for pid, p in prots.items():
        by_acc[p["acc"]].append({**p, "id": pid})
    orfs = {c: all_orfs(s.upper(), min_aa) for c, s in contigs.items()}
    with open(work / "orfs.faa", "w") as fh:
        for c, os in orfs.items():
            for i, o in enumerate(os):
                o["id"] = i
                fh.write(f">{c}::{i}\n{o['protein']}\n")
    hits: dict[str, list[tuple]] = defaultdict(list)     # contig -> (orf id, pid, ident, cover, bits, qstart, sstart)
    if (work / "orfs.faa").stat().st_size:
        run(["diamond", "makedb", "--in", proteins, "-d", work / "refs", "--quiet"])
        run(["diamond", "blastp", "--sensitive", "-d", work / "refs", "-q", work / "orfs.faa",
             "-o", work / "orfs_vs_refs.tsv", "--outfmt", "6", *COLS, "-e", "1e-5", "-p", threads,
             "--max-target-seqs", 200])
        for line in open(work / "orfs_vs_refs.tsv"):
            q, s, pid, ln, slen, bits, qs, ss = line.split("\t")
            c, i = q.rsplit("::", 1)
            hits[c].append((int(i), s, float(pid), int(ln) / int(slen), float(bits), int(qs), int(ss)))
    out = {}
    for c, os in orfs.items():
        hs = hits.get(c, [])
        per_acc: dict[str, dict[str, tuple]] = defaultdict(dict)    # acc -> pid -> (ident, cover, bits)
        for _, s, ident, cover, bits, *_ in hs:
            if bits > per_acc[prots[s]["acc"]].get(s, (0, 0, 0))[2]:
                per_acc[prots[s]["acc"]][s] = (ident, cover, bits)
        scored = {}
        for acc, best in per_acc.items():
            res = classify_genes(by_acc[acc], {k: v[:2] for k, v in best.items()}, expect_min_aa, min_identity)
            scored[acc] = (len(res["found"]), len(res["truncated"]), sum(v[2] for v in best.values()), res)
        ref = max(scored, key=lambda a: scored[a][:3], default="")
        named = []
        for o in os:
            mine = [h for h in hs if h[0] == o["id"] and h[2] >= min_identity and h[3] >= min_cover]
            if not mine:
                continue
            _, s, ident, cover, _, qs, ss = max(mine, key=lambda h: (h[4], prots[h[1]]["acc"] == ref))
            o.update(gene=prots[s]["gene"], product=prots[s]["product"], identity=ident, cover=cover)
            if ss <= 3 and qs > 1 and "M" in o["protein"][qs - 1:qs + 5]:     # start at the aligned Met
                k = qs - 1 + o["protein"][qs - 1:qs + 5].index("M")
                o["protein"] = o["protein"][k:]
                o["begin" if o["strand"] == 1 else "end"] += 3 * k * o["strand"]
            named.append(o)
        named.sort(key=lambda o: (o["begin"], o["end"]))
        for i, o in enumerate(named, 1):
            o["index"] = i
        res = scored[ref][3] if ref else {"found": [], "truncated": [], "missing": []}
        out[c] = {"orfs": named, "ref": ref, **res}
    return out


def merge_unnamed(named: list[dict], found: list[dict]) -> list[dict]:
    """Named ORFs plus gene-finder ORFs that overlap none of them (e.g. genes absent from the references)."""
    keep = [o for o in found if not any(o["begin"] <= n["end"] and n["begin"] <= o["end"] for n in named)]
    merged = sorted(named + keep, key=lambda o: (o["begin"], o["end"]))
    for i, o in enumerate(merged, 1):
        o["index"] = i
    return merged


def gene_reasons(g: dict) -> list[str]:
    return [f"missing_gene:{x}" for x in g["missing"]] + [f"truncated_gene:{x}" for x in g["truncated"]]


def run_annotate(fasta: Path, proteins: Path, outdir: Path, threads: int = 4, min_aa: int = 60,
                 expect_min_aa: int = 100) -> None:
    """Standalone: annotate any FASTA of genomes; writes annotation.tsv, genome_qc.tsv, features.tbl."""
    from .genbank import feature_table, find_orfs
    outdir.mkdir(parents=True, exist_ok=True)
    seqs = {c: s.upper() for c, s in read_fasta(fasta).items()}
    res = annotate(seqs, proteins, outdir / "work", threads, min_aa, expect_min_aa=expect_min_aa)
    with open(outdir / "annotation.tsv", "w") as an, open(outdir / "genome_qc.tsv", "w") as qc, \
            open(outdir / "features.tbl", "w") as tbl:
        an.write("contig\torf\tbegin\tend\tstrand\taa\tgene\tproduct\tidentity\tcover\n")
        qc.write("contig\tlength\tref\tfound\ttruncated\tmissing\tstatus\n")
        for c, seq in seqs.items():
            r = res[c]
            orfs = merge_unnamed(r["orfs"], find_orfs(seq, min_aa))
            for o in orfs:
                an.write("\t".join(map(str, [c, o["index"], o["begin"], o["end"], o["strand"], len(o["protein"]),
                                             o.get("gene", ""), o.get("product", ""),
                                             f"{o['identity']:.1f}" if "identity" in o else "",
                                             f"{o['cover']:.2f}" if "cover" in o else ""])) + "\n")
            qc.write("\t".join([c, str(len(seq)), r["ref"], ",".join(r["found"]), ",".join(r["truncated"]),
                                ",".join(r["missing"]), "ok" if r["ref"] and not gene_reasons(r) else "review"]) + "\n")
            tbl.write(feature_table(c, orfs))
    event("annotate done", contigs=len(seqs), out=str(outdir))
