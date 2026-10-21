"""Synthetic RCA-Nanopore reads: circular genomes -> concatemers, chimeras, host reads."""
import gzip
import itertools
import random
import subprocess

BASES = "ACGT"


def rand_seq(rng, n):
    return "".join(rng.choice(BASES) for _ in range(n))


def mutate(rng, seq, rate=0.05):
    """Substitutions and 1-bp indels at roughly `rate` per base (Nanopore-like)."""
    out = []
    for b in seq:
        r = rng.random()
        if r < rate * 0.5:
            out.append(rng.choice(BASES))
        elif r < rate * 0.75:
            continue
        elif r < rate:
            out += [b, rng.choice(BASES)]
        else:
            out.append(b)
    return "".join(out)


def rca_read(rng, genome, copies):
    start = rng.randrange(len(genome))
    rot = genome[start:] + genome[:start]
    return mutate(rng, rot * copies)


def simulate(genomes, n_reads=200, chimera_frac=0.1, host_frac=0.2, host=None,
             seed=1, copies=(2, 8)):
    """Return (reads, truth): reads = [(name, seq)], truth maps name -> label."""
    rng = random.Random(seed)
    host = host or rand_seq(rng, 20000)
    reads, truth = [], {}
    for i in range(n_reads):
        r = rng.random()
        if r < host_frac:
            s = rng.randrange(len(host) - 3000)
            seq, label = mutate(rng, host[s:s + rng.randrange(500, 3000)]), "host"
        elif r < host_frac + chimera_frac:
            a, b = rng.sample(genomes, 2) if len(genomes) > 1 else (genomes[0], genomes[0])
            seq = rca_read(rng, a, rng.randint(2, 4)) + rca_read(rng, b, rng.randint(2, 4))
            label = "chimera"
        else:
            g = rng.randrange(len(genomes))
            seq, label = rca_read(rng, genomes[g], rng.randint(*copies)), f"virus{g}"
        name = f"read{i}"
        reads.append((name, seq))
        truth[name] = label
    return reads, truth


def write_fastq(path, reads, qual="I"):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt") as fh:
        for name, seq in reads:
            fh.write(f"@{name}\n{seq}\n+\n{qual * len(seq)}\n")


# Standard genetic code, first codon per amino acid (used to back-translate proteins).
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {}
for _aa, _cod in zip(_AAS, ("".join(c) for c in itertools.product("TCAG", repeat=3))):
    if _aa != "*":
        CODON.setdefault(_aa, _cod)
AA = list(CODON)


def rand_protein(rng, n=400):
    return "".join(rng.choice(AA) for _ in range(n))


def back_translate(protein):
    return "".join(CODON[a] for a in protein)


def make_genome(rng, protein, n=2700, motif="TAATATTAC"):
    """Random circular genome with `protein` encoded at 500 (ATG..TAA) and an origin motif at 2000."""
    orf = "ATG" + back_translate(protein) + "TAA"
    g = rand_seq(rng, n)
    g = g[:500] + orf + g[500 + len(orf):]
    return g[:2000] + motif + g[2000 + len(motif):]


def make_tiny_db(tmp_path, proteins_by_taxid):
    """DIAMOND db + NCBI-style taxdump. Taxa: 1 root, 2 Viruses, 10 FamA, 11 GenA, 12 GenB.
    Returns (db.dmnd path, taxdump dir)."""
    tax = tmp_path / "tax"
    tax.mkdir()
    rows = [(1, 1, "no rank"), (2, 1, "superkingdom"), (10, 2, "family"), (11, 10, "genus"), (12, 10, "genus")]
    (tax / "nodes.dmp").write_text("".join(f"{a}\t|\t{b}\t|\t{r}\t|\t\t|\n" for a, b, r in rows))
    names = {1: "root", 2: "Viruses", 10: "FamA", 11: "GenA", 12: "GenB"}
    (tax / "names.dmp").write_text("".join(f"{k}\t|\t{v}\t|\t\t|\tscientific name\t|\n" for k, v in names.items()))
    for f in ("merged.dmp", "delnodes.dmp"):
        (tax / f).touch()
    (tmp_path / "db.faa").write_text("".join(f">P{t}\n{p}\n" for t, p in proteins_by_taxid.items()))
    (tmp_path / "map.tsv").write_text("accession\taccession.version\ttaxid\tgi\n"
                                      + "".join(f"P{t}\tP{t}\t{t}\t0\n" for t in proteins_by_taxid))
    subprocess.run(["diamond", "makedb", "--in", tmp_path / "db.faa", "-d", tmp_path / "db",
                    "--taxonmap", tmp_path / "map.tsv", "--taxonnodes", tax / "nodes.dmp",
                    "--taxonnames", tax / "names.dmp"], check=True, capture_output=True)
    return tmp_path / "db.dmnd", tax


def revcomp_seq(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]
