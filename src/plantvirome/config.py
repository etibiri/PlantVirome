"""YAML configuration schema. Nothing here is virus- or host-specific."""
from pathlib import Path
from typing import Literal

import yaml
from pydantic import BaseModel, ConfigDict, Field

# Execution order; each name maps to a `plantvirome <name>` subcommand.
STEPS = ["qc", "hostdeplete", "rca", "assemble", "classify", "species", "satellite",
         "abundance", "diversity", "genbank", "variety"]
# Steps that run once for the whole cohort instead of once per sample.
COHORT_STEPS = {"diversity"}


class _Strict(BaseModel):
    model_config = ConfigDict(extra="forbid")


class Host(_Strict):
    name: str = "unspecified"
    # Any mix of nuclear (partial), organelle or transcript FASTA files; may be empty.
    references: list[Path] = []
    # Organelle/marker references used by the optional `variety` step.
    variety_references: list[Path] = []


class QC(_Strict):
    min_length: int = 200
    max_length: int | None = None
    min_quality: float = 10.0
    # Adapter trimming only; internal-adapter splitting is off (RCA reads have none).
    trim_adapters: bool = True


class HostDeplete(_Strict):
    # A read is host if alignments to the host panel cover >= min_query_cov of it
    # at >= min_identity. Partial references leave residual host reads by design.
    min_query_cov: float = 0.5
    min_identity: float = 0.8


class RCA(_Strict):
    # Tandem-repeat (concatemer) detection; the period is the genome/unit length.
    # Set min/max period to bracket the genome sizes you expect, not a specific virus.
    min_copies: int = Field(2, ge=2)
    min_period: int = 300
    max_period: int = 20000
    # A read is "clean" if one repeat region covers >= min_read_cover of it.
    min_read_cover: float = 0.8


class Assemble(_Strict):
    # units: assemble the RCA unit consensus sequences (default, fast, circular by construction).
    # flye: de novo assembly of the raw reads with Flye --meta (see steps/assemble.py).
    method: Literal["units", "flye"] = "units"
    flye_mode: str = "nano-hq"            # flye --<flye_mode>
    flye_reads: Literal["nonhost", "qc"] = "nonhost"   # host-depleted reads, or all QC'd reads
    # Random subsample of the reads to about this many bases before Flye (None = use all).
    # Concatemers are highly redundant; deep input makes Flye very slow and memory hungry.
    flye_max_bases: int | None = 20_000_000
    flye_seed: int = 1
    # Collapse contigs made of tandem copies of the genome (always the case for raw RCA reads).
    flye_collapse_multimers: bool = True
    # Drop Flye contigs below this Flye-reported depth (on the subsampled input): chimeras and
    # mis-joined concatemers sit at 3-10x, real genomes and variants at >= 17x (barcode01, 80 Mb).
    flye_min_coverage: float = 15.0
    # Units (one genome copy each, random rotation/strand) join a cluster if they align to
    # the seed unit (written twice, to absorb rotation) at >= min_identity over >= min_cover.
    min_identity: float = 0.85
    min_cover: float = 0.9
    # Clusters with fewer units are noise and are dropped.
    min_cluster_units: int = 3
    # Below this many clean units the sample is reported as empty, not an error.
    min_units: int = 20
    # Contigs are circular by construction (built from concatemers); this many
    # independent units are required before the circular flag is set.
    circ_min_support: int = 3
    # Racon polishing of the final contigs with the host-depleted reads (0 = off).
    # Circular contigs are padded with their own ends so the junction is polished too.
    polish_rounds: int = 2
    # Medaka model of the basecaller used, e.g. r1041_e82_400bps_sup_v5.2.0 (None = no Medaka).
    # Runs once after Racon; `medaka tools list_models` shows the installed models.
    polish_medaka_model: str | None = None
    polish_max_bases: int | None = 50_000_000   # subsample the reads to about this many bases


class Classify(_Strict):
    # Required in `databases`: viral_protein (DIAMOND .dmnd built with --taxonmap and
    # --taxonnodes) and taxdump (directory with NCBI nodes.dmp / names.dmp).
    evalue: float = 1e-5
    sensitivity: str = "sensitive"   # diamond --<sensitivity>
    # Hits within top_pct % of the best bitscore are combined by LCA.
    top_pct: float = 10.0
    # (min protein identity %, label), highest first. Below the last: "unclassified".
    tiers: list[tuple[float, str]] = [(85.0, "close"), (60.0, "related"), (35.0, "distant")]


class Species(_Strict):
    # Optional nucleotide-level species assignment against reference genomes made by
    # `plantvirome download-refs`. The step is skipped when `reference` is not set.
    reference: Path | None = None
    taxon: str | None = None            # genus name, used to label novel candidates ("<taxon> sp.")
    # Pairwise identity at/above which a contig belongs to the reference's species.
    # Begomovirus DNA-A: 0.91 (ICTV, Brown et al. 2015); check the current ICTV criteria.
    species_threshold: float = 0.91
    strain_threshold: float | None = None   # begomovirus: 0.94; None = no strain call
    min_cover: float = 0.8              # aligned fraction of the contig for a full comparison
    min_partial_cover: float = 0.3      # below this there is no usable alignment at all
    borderline_margin: float = 0.01     # |identity - threshold| below this is flagged borderline
    # A contig covering less of the reference genome than this is a fragment (defective or
    # subgenomic DNA, or an assembly piece), not a complete genome: caveat partial_genome.
    min_ref_cover: float = 0.8


class SatelliteGroup(_Strict):
    # FASTA made by `plantvirome download-refs` for one satellite taxon (e.g. Alphasatellitidae).
    reference: Path
    # Nucleotide identity at/above which a contig is the same species as its nearest reference.
    # ICTV (verify the current criteria): alphasatellites 0.88 (Geminialphasatellitinae) or
    # 0.80 (Nanoalphasatellitinae); betasatellites and deltasatellites 0.91. None = report the
    # identity only (satellite_related), no species call.
    species_threshold: float | None = None
    # Contig length window searched against this group (alpha/betasatellites are ~1.3-1.4 kb).
    min_length: int = 1000
    max_length: int = 2000


class Satellite(_Strict):
    # Optional: {group_name: SatelliteGroup}. The step is skipped when empty.
    groups: dict[str, SatelliteGroup] = {}
    min_cover: float = 0.8
    min_partial_cover: float = 0.3
    min_ref_cover: float = 0.8
    borderline_margin: float = 0.01


class Abundance(_Strict):
    min_identity: float = 0.8
    # A read is assigned to a contig only if >= min_read_cover of it aligns there.
    min_read_cover: float = 0.5
    # A contig counts as detected if >= min_breadth of its positions have depth >= 1.
    min_breadth: float = 0.9
    # Read pieces on two contigs, each >= chimera_min_frac of the read, flag a chimera.
    chimera_min_frac: float = 0.2


class Diversity(_Strict):
    # Contigs from all samples are clustered into vOTUs (viral operational taxonomic
    # units): >= cluster_identity nucleotide identity over >= cluster_min_cover of the
    # shorter contig. Set the identity to the species threshold of your virus group.
    cluster_identity: float = 0.90
    cluster_min_cover: float = 0.80      # aligned fraction of the LONGER contig (genome-level)
    # Shorter contigs (phi29 artifacts, fragments) are not clustered: they bridge genomes that
    # share a small region, e.g. DNA-A and DNA-B through the common region.
    min_contig_length: int = 1000
    # Leave out contigs with no protein tier and no nucleotide species call (host or
    # unidentified DNA that survived depletion would otherwise count as vOTUs).
    exclude_unclassified: bool = False
    # Negative controls (sample names). They are clustered with the samples but left out of
    # alpha/beta diversity; vOTUs also seen in a blank are reported in contamination.tsv.
    blanks: list[str] = []
    # A sample's vOTU is flagged if its depth is below blank_ratio x the depth in the blank.
    blank_ratio: float = 10.0


class GenBank(_Strict):
    min_length: int = 200                  # GenBank minimum sequence length
    max_n_fraction: float = 0.01           # ambiguous bases tolerated per contig
    min_depth: float = 20.0                # mean read depth needed to trust the consensus
    allowed_tiers: list[str] = ["close", "related"]   # classify tiers eligible for submission
    submit_partial: bool = False           # submit non-circular contigs too?
    max_chimeric_fraction: float = 0.5     # chimeric / (assigned + chimeric) reads above this -> review
    # Circular genomes are rotated so the sequence starts at origin_motif[origin_offset]
    # (either strand is searched). None = keep the assembler's coordinates.
    # begomovirus: motif TAATATTAC, offset 7 (nick site TAATATT^AC).
    origin_motif: str | None = None
    origin_offset: int = 0
    # Substitutions tolerated in the motif (consensus errors); >0 can make the match ambiguous.
    origin_max_mismatch: int = 0
    # Orient contigs on the strand of the best reference genome (species step); the motif then
    # places the origin on that strand only, which avoids palindromic double hits. Without a
    # usable motif the contig starts where the reference starts.
    orient_by_reference: bool = False
    min_orf_aa: int = 60
    # CheckV database directory: when set, CheckV runs on each sample's contigs. Contamination and
    # host genes can send a contig to review; its completeness is informational only (it compares
    # with reference genomes and underestimates small circular genomes such as begomovirus DNA-A/B).
    checkv_db: Path | None = None
    checkv_max_contamination: float = 10.0   # % above this -> review
    # Source modifiers written into FASTA deflines, e.g. {host: Manihot esculenta, country: ...}
    metadata: dict[str, str] = {}
    # Optional CSV with a `sample` column; its other columns override `metadata` per sample.
    metadata_csv: Path | None = None


class Variety(_Strict):
    min_identity: float = 0.8
    min_mapq: int = 0           # 0 keeps reads from repeats such as the chloroplast inverted repeat
    min_depth: int = 10         # positions below this are "uncovered" (N in the consensus)
    min_alt_frac: float = 0.8   # organelle genomes are treated as homoplasmic
    # Coverage is reported PARTIAL if fewer positions than this reach min_depth.
    complete_breadth: float = 0.9


class Resources(_Strict):
    cpus: int = 4
    mem_gb: int = 16
    time: str = "04:00:00"


class Slurm(_Strict):
    partition: str | None = None
    account: str | None = None
    # Shell lines run before each job, e.g. "module load conda" / "conda activate plantvirome".
    setup: list[str] = []
    default: Resources = Resources()
    steps: dict[str, Resources] = {}


class Config(_Strict):
    host: Host = Host()
    databases: dict[str, Path] = {}
    qc: QC = QC()
    hostdeplete: HostDeplete = HostDeplete()
    rca: RCA = RCA()
    assemble: Assemble = Assemble()
    classify: Classify = Classify()
    species: Species = Species()
    satellite: Satellite = Satellite()
    abundance: Abundance = Abundance()
    diversity: Diversity = Diversity()
    genbank: GenBank = GenBank()
    variety: Variety = Variety()
    slurm: Slurm = Slurm()
    threads: int = Field(4, ge=1)

    def resources(self, step: str) -> Resources:
        return self.slurm.steps.get(step, self.slurm.default)


def load_config(path: Path | str) -> Config:
    with open(path) as fh:
        data = yaml.safe_load(fh) or {}
    return Config.model_validate(data)
