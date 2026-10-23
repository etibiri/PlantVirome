# plantvirome

Virus-agnostic pipeline for **RCA (rolling-circle amplification) + Nanopore** plant virome data.
From basecalled reads to circular virus genomes, taxonomy, diversity indices and
GenBank-ready files. Nothing in the code is tied to a host or a virus: hosts, databases and
virus-specific parameters (genome size range, origin motif, species threshold) live in YAML.

Pipeline (one subcommand per step, each resumable through a `.done` marker):

```
qc → hostdeplete → rca → assemble → classify → species → satellite → abundance → genbank   (per sample)
                                               └→ diversity            (once, all samples)
                                    variety (optional, host organelle markers)
```

## Install

Requires Linux or macOS. External tools come from bioconda, the Python package from pip.

```bash
mamba env create -f environment.yml      # tools + `pip install -e .`
mamba activate plantvirome
pytest                                   # unit + end-to-end tests on simulated data (~2 min)
```

Container (Docker; Apptainer for the HPC):

```bash
docker build -f containers/Dockerfile -t plantvirome .
apptainer build plantvirome.sif containers/plantvirome.def
```

Tools: `minimap2`, `chopper`, `porechop_abi`, `TideHunter`, `DIAMOND`, `taxonkit`, `Flye` (only for `assemble.method: flye`), `Racon` (polishing)
(Python: `typer`, `pydantic`, `pyyaml`, `numpy`, `scipy`, `pyrodigal-gv`).
Tests that need a missing tool are skipped, not failed.

## Run

```bash
# one sample, locally (reads: FASTQ(.gz) files or a directory)
plantvirome run --config my.yaml --outdir results --sample s1 --reads fastq/s1/

# many samples on Slurm: writes results/slurm/<step>.<sample>.sh, chained with afterok
plantvirome submit --config my.yaml --outdir results --sample s1 --sample s2 --dry-run
plantvirome submit --config my.yaml --outdir results --sample s1 --sample s2

# cohort step after all samples are done (also part of `submit`)
plantvirome diversity --config my.yaml --outdir results
```

Individual steps: `plantvirome qc|hostdeplete|rca|assemble|classify|species|satellite|abundance|genbank|variety
--config … --outdir … --sample …` (`--force` reruns a finished step).
Logs are JSON lines (one object per event) on stderr; under Slurm they land in
`results/slurm/logs/`.

Without `--reads`, `qc` reads `results/<sample>/input/`.

### Databases

`classify` needs, in `databases`:

* `viral_protein` — a DIAMOND database **built with taxonomy**,
* `taxdump` — directory with NCBI `nodes.dmp`, `names.dmp`, `merged.dmp`, `delnodes.dmp`.

Build both from NCBI RefSeq viral proteins (about 320 MB download, 3 minutes, ~350 MB on disk):

```bash
plantvirome download-db --outdir db --threads 8      # add --keep-raw to keep the downloads
```

It prints the two lines to paste under `databases:` and records the source URLs and date in
`db/download.json` (cite it in your methods). Taxids come from the RefSeq `.gpff` files, so no
multi-GB accession2taxid file is needed. To use another protein set, build it yourself with
`diamond makedb --in x.faa -d x --taxonmap … --taxonnodes … --taxonnames …`.
Because related species share near-identical proteins, the LCA often stops at genus
(a conserved coat protein of a cassava mosaic virus classified as *Begomovirus*): species-level
calls need nucleotide comparison of whole genomes, which `plantvirome` does not do.

### Reference genomes for species calls (optional)

```bash
plantvirome download-refs --taxon Begomovirus --outdir db/begomovirus   # ~1 min, 725 RefSeq genomes
```

Fetches the RefSeq genomes of a taxon from NCBI (Entrez) into `refs.fasta`
(`accession|species|taxid|component`, component `A` = DNA-A or monopartite, `B` = DNA-B,
`satellite`), plus `refs.tsv` and `download.json`. Point `species.reference` at `refs.fasta`.
Without it the `species` step is skipped and nothing else changes.
`--extra ACC=TAXID` (repeatable) adds non-RefSeq records, labelled with the species taxid
`TAXID` (without it, the record's own, often strain-level, taxon), e.g.
`--extra AJ717531.1=62079 --extra AJ704962.1=62079` for EACMV isolates.

### Configuration

Start from `default.yaml`; presets: `begomovirus.yaml` (genome size range, 91 %
species threshold, `TAATATTAC` origin), `slurm.yaml` (partition, resources per step). They ship
with the package: `cp "$(plantvirome presets)"/default.yaml my.yaml` (source: `src/plantvirome/presets/`).
Unknown keys are rejected. Validate with `plantvirome validate-config my.yaml`.

| Section | Key parameters (default) |
|---|---|
| `host` | `references` (host FASTA, may be partial or empty), `variety_references` (organelles/markers) |
| `databases` | `viral_protein`, `taxdump` |
| `qc` | `min_length` 200, `max_length` none, `min_quality` 10, `trim_adapters` true |
| `hostdeplete` | `min_query_cov` 0.5, `min_identity` 0.8 |
| `rca` | `min_copies` 2, `min_period` 300, `max_period` 20000, `min_read_cover` 0.8 |
| `assemble` | `method` units/flye, `flye_mode` nano-hq, `flye_reads` nonhost/qc, `flye_max_bases` 20000000, `flye_collapse_multimers` true, `flye_min_coverage` 15; units method: `min_identity` 0.85, `min_cover` 0.9, `min_cluster_units` 3, `min_units` 20, `circ_min_support` 3; polishing: `polish_rounds` 2 (0 = off), `polish_medaka_model` null, `polish_max_bases` 50000000 |
| `classify` | `evalue` 1e-5, `sensitivity` sensitive, `top_pct` 10, `tiers` close ≥85 / related ≥60 / distant ≥35 (protein identity %) |
| `species` | `reference` (none = step skipped), `taxon`, `species_threshold` 0.91, `strain_threshold` none (begomovirus: 0.94), `min_cover` 0.8, `borderline_margin` 0.01 |
| `satellite` | `groups` {name: `reference`, `species_threshold`, `min_length`, `max_length`}, `min_cover` 0.8, `min_ref_cover` 0.8 |
| `abundance` | `min_identity` 0.8, `min_read_cover` 0.5, `min_breadth` 0.9, `chimera_min_frac` 0.2 |
| `diversity` | `cluster_identity` 0.90, `cluster_min_cover` 0.80, `exclude_unclassified` false, `blanks` [], `blank_ratio` 10 |
| `genbank` | `min_length` 200, `max_n_fraction` 0.01, `min_depth` 20, `allowed_tiers`, `submit_partial` false, `origin_motif`/`origin_offset`/`origin_max_mismatch` 0, `orient_by_reference` false, `min_orf_aa` 60, `metadata`, `metadata_csv` |
| `variety` | `min_depth` 10, `min_alt_frac` 0.8, `complete_breadth` 0.9 |
| `slurm` | `partition`, `account`, `setup` (shell lines), `default` and per-`steps` `cpus`/`mem_gb`/`time` |

## What each step does, and what to watch

Layout: `results/<sample>/<step>/`.

**qc** — `porechop_abi` adapter trimming (no internal splitting; set `trim_adapters: false`
for dorado output that already trimmed adapters and barcodes, since the ab-initio adapter search
takes many minutes per 100k reads) and `chopper` length/quality filter. Outputs `reads.qc.fastq.gz`, `stats.json`. *Do not set `max_length` low*: concatemers
of 10–30 kb are the signal. Do not raise `min_quality` much; later steps use the extra copies
to correct errors.

**hostdeplete** — maps reads to the host panel (`minimap2 -c`); a read is host if alignments
cover ≥ `min_query_cov` of it. Reads are split into `reads.host/nonhost.fastq.gz`, never
deleted. With an incomplete host genome, reads from missing regions stay in the non-host
set (classification removes them later); reads from endogenous viral elements in the host
panel are labelled host.

**rca** — `TideHunter` finds the tandem repeat in each read and calls an error-corrected
consensus of one genome copy ("unit"). Classes: `clean` (one repeat region covers ≥ 80 % of
the read), `multi`/`partial` (putative chimeras from template switching; units go to
`units.flagged.fasta`, never into assembly), `untandem` (no concatemer, not used). Set
`min_period`/`max_period` around the genome sizes you expect (satellites ~1.3 kb to
bipartite/ssDNA genomes ~3 kb). Chimeras between near-identical genomes can look clean.

**assemble** — two methods, selected with `assemble.method`:

* `units` (default) — units are clustered and polished directly (seed + majority vote; see the
  module docstring). Fast and deterministic. Contigs are circular by construction
  (`circ_support` = number of units). Strains closer than `min_identity` collapse into one
  majority consensus; samples with fewer than `min_units` clean units give an empty assembly.
* `flye` — de novo assembly of the **raw** host-depleted reads with `Flye --meta`
  (`flye_reads: qc` uses all QC'd reads). Raw concatemers make Flye assemble tandem multimers
  (a 2.7 kb genome comes out as a 13.5 kb "circular" contig), so such contigs are collapsed to
  one genome copy with TideHunter (`flye_collapse_multimers`; `circ_support` = copies found).
  It is slow on deep, redundant data (about 1.5 min for 4 Mb, 4 threads), so reads are randomly
  subsampled to `flye_max_bases` (20 Mb by default; `null` = all): **rare viruses can be lost
  when subsampling**. Chimeric reads are not filtered here, and Flye was unstable on the
  pre-corrected units, which is why `units` does not use it. Running both methods on the same
  sample and comparing the contigs is a good consistency check.

  *Polishing* (both methods): `polish_rounds` rounds (default 2) of `Racon` with the host-depleted
  reads (subsampled to `polish_max_bases`; minimap2 `map-ont`, one alignment per read, so a
  concatemer counts once per round), then, if `polish_medaka_model` is set (e.g.
  `r1041_e82_400bps_sup_v5.2.0`; it must match the basecalling model), one `Medaka` round.
  Circular contigs are padded with their own other end and cut back at an anchor k-mer, so
  the junction is polished too (Racon alone trims the ends, or with `--no-trimming` inserts
  bases at window edges); the origin moves by about 10 bp per round, which `genbank`
  re-fixes. Contigs under 300 bp, without read support, or whose length changes by more
  than 5 % stay unpolished. Pre-polish sequences are in `contigs.prepolish.fasta`.
  Medaka is **not** in `environment.yml` (heavy, and its conda solve was not tested with the
  other tools here): install it in its own environment and put `medaka_consensus`, `bgzip`
  and `tabix` on the `PATH`. On two real SUP contigs (barcode01, about 2.8 kb), Racon changed
  2 and 17 positions relative to the units consensus, and Medaka then reverted almost all of
  them (0 and 2 differences left versus the units consensus): Racon's edits looked like
  errors, not corrections. That is two contigs and no truth set, so treat it as a hint:
  for SUP data prefer `polish_rounds: 0` with Medaka, or no polishing, and confirm on your
  own data. Polishing cannot repair a consensus built from the wrong strain mix, or errors
  shared by all reads.

**classify** — `DIAMOND blastx` (frameshift-tolerant) + LCA over hits within `top_pct`
of the best (`taxonkit`). `taxonomy.tsv` has taxon, rank, lineage, best identity, query
coverage and a tier. Tiers are protein-identity bands, **not** ICTV species thresholds —
tune them per virus group. No hit ⇒ `unclassified` (novel virus, protein-free satellite,
host or organelle sequence).

**species** (optional) — nucleotide-level species assignment, because protein hits often stop
at genus. Each contig is aligned (`minimap2 -x asm20`) to the reference genomes and gets, in
`species.tsv`: `known_species` (identity ≥ `species_threshold` over ≥ 80 % of the contig; the
species of the best reference), `novel_species_candidate` (below it; the nearest species is
reported), `no_reference_match`, or `other_component` (best hit is DNA-B or a satellite: no species
decision, since the threshold applies to DNA-A). Caveats: `borderline` (within 1 point of the
threshold), `partial_alignment`, `new_strain` (below `strain_threshold`). Limits: identity is
alignment-based and reads ~0.5 points low at 91 %; alignments fragment below ~86 % identity, so
more divergent contigs are reported by coverage only; references are RefSeq exemplars, not
every member of a species. **Confirm every borderline and novel call with a multiple alignment
(e.g. SDT or MAFFT) and the current ICTV criteria** before naming anything. The species name is
the NCBI organism of the closest reference (can be a strain-level taxon, or a new ICTV binomial
such as *Begomovirus manihotis*). `genbank` uses a species call for the organism name.

**satellite** (optional) — satellite DNAs (alpha-, beta-, deltasatellites) are circular, 1 to 1.4 kb
and share little with the helper virus, so they get their own step. Build one reference set per
group (`plantvirome download-refs --taxon Alphasatellitidae --outdir db/alphasatellite`, likewise
`Tolecusatellitidae` for betasatellites and `Deltasatellite`) and list them under
`satellite.groups` with a length window and a species threshold. Contigs inside a group's window
are aligned to its references and called `known_satellite`, `novel_satellite_candidate`
(below the threshold), or `satellite_related` (group without a threshold); contigs the protein
search put in a satellite family but with no nucleotide match are `protein_only_candidate`.
`species.tsv`-style columns, including the reference start used by `genbank` for orientation.
ICTV species thresholds (check the current criteria): alphasatellites 88 % (Geminialphasatellitinae)
or 80 % (Nanoalphasatellitinae), betasatellites and deltasatellites 91 %. Limits: nucleotide search
finds satellites down to roughly 85 % identity to a reference; divergent alphasatellites show up
only through their Rep protein, and betasatellites (one short protein) mostly not at all.
A 1 to 2 kb circular contig is not evidence of a satellite: defective DNAs and phi29 artifacts
fall in the same window.

**abundance** — reads re-mapped to contigs; abundance is **mean depth** (aligned bases /
length), not read counts, because a concatemer carries many copies. Reports breadth,
depth CV, a `detected` flag, relative abundance among detected contigs, and contig pairs
joined by chimeric reads (`chimeras.tsv`). RCA amplifies unevenly and favours small circular
DNA, so abundances are semi-quantitative: use them to rank, not as viral loads.

**diversity** (cohort) — contigs of all samples are clustered into vOTUs (≥ `cluster_identity`
over ≥ `cluster_min_cover`), so novel viruses count without a database hit. Alpha: observed
vOTUs, Shannon, Simpson, Pielou, Chao1. Beta: Bray-Curtis (relative depth), Jaccard
(presence/absence), Aitchison (CLR); PCoA of Bray-Curtis for ≥ 3 samples. Limits: this is
the diversity of the **RCA-amplifiable** virome (circular DNA viruses and satellites; most
RNA viruses are missed); Chao1 treats reads as individuals, which concatemers violate (read
it as a rough lower bound); richness depends on the identity threshold; no rarefaction
is applied, so compare samples of very different depth with care.

**Negative controls** — name them in `diversity.blanks`. They go through every step, are
clustered with the samples, and are left out of alpha/beta diversity. `contamination.tsv`
lists, for each contig whose vOTU is also in a blank, the ratio of its depth to the blank's;
below `blank_ratio` (10×) it is flagged `possible_contamination`, and `genbank` then sends
the contig to review. Run `diversity` before `genbank` so the file exists. A ratio is only a
screen: index hopping from a high-titre sample into the blank looks the same as
contamination, so check the flagged vOTUs against the barcodes next to them.

**genbank** — per-contig status `submit` / `review` / `reject` with reasons in
`contigs.qc.tsv`. Checks: length, ambiguous bases, depth, ID length (≤ 25), taxonomy tier,
circularity, chimeric-read share, origin found exactly once. Circular genomes are rotated
to start at `origin_motif[origin_offset]`. With `orient_by_reference` the strand is taken from the
best reference genome (species step) and the motif is searched on that strand only; this avoids
the double hit of palindromic motifs (the begomovirus nonanucleotide `TAATATTA` is its own
reverse complement) and makes `origin_max_mismatch: 1` usable for consensus errors. If the motif
is not unique the contig starts where the reference starts (`origin_from_reference` note).
Without reference orientation both strands are searched. ORFs come from `pyrodigal-gv`;
`features.tbl` lists genes/CDS as "hypothetical protein" (`orf1`, `orf2`, …) because no
curated gene names are shipped. `submission.fsa` has table2asn-style definition lines
with `[organism=…][moltype=genomic][molecule=DNA][topology=…]` plus your `metadata`.
**You must check** the proposed organism (species name, else `Genus sp.`, else the contig
goes to review), fill country/collection date, and validate with NCBI's tool, which this
package does not run:
`table2asn -i submission.fsa -f features.tbl -t template.sbt -V vb`.
ORFs crossing the origin are missed (fine when the origin is intergenic, as in begomoviruses).

**variety** (optional) — maps all QC'd reads to `host.variety_references`, calls SNVs
(`variants.tsv`), writes `consensus.fasta` with `N` where depth < `min_depth` and
`coverage.tsv` flagged **PARTIAL** unless ≥ 90 % of each reference reaches that depth.
Coverage of host DNA in an RCA virome is partial and uneven by nature. SNVs only (no
indels), organelles assumed homoplasmic, and copies of one concatemer are not independent
evidence. Naming a cultivar further needs a panel of known cultivar haplotypes (not
included).

## Tests

```bash
pytest tests/unit          # fast logic tests
pytest tests/integration   # simulated two-sample run: reads → submission files → diversity
```

`tests/simulate.py` builds synthetic circular genomes, RCA concatemers with Nanopore-like
errors, chimeras and host reads, so no real data is needed.

## Status

Developed and tested on Linux (x86_64) with simulated data only; the Docker image builds and
passes the full test suite, including the end-to-end run. Not yet validated on real
RCA-Nanopore runs (thresholds will need tuning), on a live Slurm cluster (scripts are
generated and their command lines checked, but never submitted), or by building the
Apptainer image. macOS: every conda tool has osx-64 and osx-arm64 packages on bioconda, but
`environment.yml` has not been installed on a Mac.

## License

MIT, see `LICENSE`.
