"""Command line interface. Each pipeline step is registered as a subcommand."""
from pathlib import Path

import typer

from . import __version__
from .config import COHORT_STEPS, STEPS, load_config
from .logging import setup_logging
from .steps import hostdeplete as hostdeplete_step
from .steps import assemble as assemble_step
from .steps import classify as classify_step
from .steps import abundance as abundance_step
from .steps import diversity as diversity_step
from .steps import genbank as genbank_step
from .steps import variety as variety_step
from .steps import download_db as download_db_step
from .steps import download_refs as download_refs_step
from .steps import rca as rca_step
from .steps import satellite as satellite_step
from .steps import species as species_step
from .steps import qc as qc_step
from .slurm import submit as slurm_submit
from .slurm import write_scripts

app = typer.Typer(no_args_is_help=True, help="RCA-Nanopore plant virome pipeline.")


@app.callback()
def main(verbose: bool = typer.Option(False, "--verbose", "-v")) -> None:
    import logging
    setup_logging(level=logging.DEBUG if verbose else logging.INFO)


@app.command()
def version() -> None:
    """Print the version."""
    typer.echo(__version__)


@app.command("validate-config")
def validate_config(config: Path = typer.Argument(..., exists=True)) -> None:
    """Validate a YAML configuration file."""
    load_config(config)
    typer.echo("OK")


@app.command()
def submit(
    config: Path = typer.Option(..., exists=True),
    outdir: Path = typer.Option(Path("results")),
    sample: list[str] = typer.Option(..., "--sample", help="Repeat for each sample."),
    steps: str = typer.Option(",".join(STEPS), help="Comma-separated steps."),
    dry_run: bool = typer.Option(False, help="Write scripts but do not call sbatch."),
) -> None:
    """Generate Slurm scripts (and submit them unless --dry-run)."""
    cfg = load_config(config)
    scripts = write_scripts(cfg, sample, steps.split(","), config.resolve(), outdir.resolve())
    for step, smp, path in scripts:
        typer.echo(f"{step}\t{smp or '-'}\t{path}")
    if not dry_run:
        slurm_submit(scripts)



@app.command()
def qc(
    config: Path = typer.Option(..., exists=True),
    outdir: Path = typer.Option(Path("results")),
    sample: str = typer.Option(...),
    reads: list[Path] = typer.Option(None, help="FASTQ(.gz) files or directories; "
                                     "default: <outdir>/<sample>/input/*.fastq*"),
    threads: int = typer.Option(None),
    force: bool = typer.Option(False, help="Re-run even if already done."),
) -> None:
    """Adapter trimming, length/quality filtering and read statistics."""
    cfg = load_config(config)
    qc_step.run_qc(cfg, outdir / sample, reads or None, threads or cfg.threads, force)


@app.command()
def hostdeplete(
    config: Path = typer.Option(..., exists=True),
    outdir: Path = typer.Option(Path("results")),
    sample: str = typer.Option(...),
    threads: int = typer.Option(None),
    force: bool = typer.Option(False, help="Re-run even if already done."),
) -> None:
    """Split reads into host / non-host using the configured host panel."""
    cfg = load_config(config)
    hostdeplete_step.run_hostdeplete(cfg, outdir / sample, threads or cfg.threads, force)


@app.command()
def rca(
    config: Path = typer.Option(..., exists=True),
    outdir: Path = typer.Option(Path("results")),
    sample: str = typer.Option(...),
    threads: int = typer.Option(None),
    force: bool = typer.Option(False, help="Re-run even if already done."),
) -> None:
    """Split RCA concatemers into unit consensus sequences and flag chimeric reads."""
    cfg = load_config(config)
    rca_step.run_rca(cfg, outdir / sample, threads or cfg.threads, force)


@app.command()
def assemble(
    config: Path = typer.Option(..., exists=True),
    outdir: Path = typer.Option(Path("results")),
    sample: str = typer.Option(...),
    threads: int = typer.Option(None),
    force: bool = typer.Option(False, help="Re-run even if already done."),
) -> None:
    """De novo assembly of RCA unit consensus sequences, with circularity check."""
    cfg = load_config(config)
    assemble_step.run_assemble(cfg, outdir / sample, threads or cfg.threads, force)


@app.command()
def classify(
    config: Path = typer.Option(..., exists=True),
    outdir: Path = typer.Option(Path("results")),
    sample: str = typer.Option(...),
    threads: int = typer.Option(None),
    force: bool = typer.Option(False, help="Re-run even if already done."),
) -> None:
    """Taxonomic assignment of contigs (DIAMOND blastx + LCA)."""
    cfg = load_config(config)
    classify_step.run_classify(cfg, outdir / sample, threads or cfg.threads, force)


@app.command()
def abundance(
    config: Path = typer.Option(..., exists=True),
    outdir: Path = typer.Option(Path("results")),
    sample: str = typer.Option(...),
    threads: int = typer.Option(None),
    force: bool = typer.Option(False, help="Re-run even if already done."),
) -> None:
    """Map reads to contigs: depth, breadth, relative abundance, chimera check."""
    cfg = load_config(config)
    abundance_step.run_abundance(cfg, outdir / sample, threads or cfg.threads, force)


@app.command()
def diversity(
    config: Path = typer.Option(..., exists=True),
    outdir: Path = typer.Option(Path("results")),
    sample: list[str] = typer.Option(None, "--sample",
                                     help="Default: every <outdir>/*/abundance/.done"),
    threads: int = typer.Option(None),
    force: bool = typer.Option(False, help="Re-run even if already done."),
) -> None:
    """Cohort step: vOTU table, alpha and beta diversity, PCoA."""
    cfg = load_config(config)
    diversity_step.run_diversity(cfg, outdir, sample or None, threads or cfg.threads, force)


@app.command()
def genbank(
    config: Path = typer.Option(..., exists=True),
    outdir: Path = typer.Option(Path("results")),
    sample: str = typer.Option(...),
    threads: int = typer.Option(None),
    force: bool = typer.Option(False, help="Re-run even if already done."),
) -> None:
    """Contig QC, origin rotation, ORF annotation and GenBank-ready files."""
    cfg = load_config(config)
    genbank_step.run_genbank(cfg, outdir / sample, force)


@app.command()
def variety(
    config: Path = typer.Option(..., exists=True),
    outdir: Path = typer.Option(Path("results")),
    sample: str = typer.Option(...),
    threads: int = typer.Option(None),
    force: bool = typer.Option(False, help="Re-run even if already done."),
) -> None:
    """Host variety markers: organelle variants and coverage (partial by design)."""
    cfg = load_config(config)
    variety_step.run_variety(cfg, outdir / sample, threads or cfg.threads, force)


@app.command()
def run(
    config: Path = typer.Option(..., exists=True),
    outdir: Path = typer.Option(Path("results")),
    sample: str = typer.Option(...),
    reads: list[Path] = typer.Option(None, help="FASTQ(.gz) files or directories for qc."),
    steps: str = typer.Option(",".join(s for s in STEPS if s not in COHORT_STEPS),
                              help="Comma-separated; add `diversity` to also run the cohort step "
                                   "(on this sample plus nothing else: use `diversity` for a cohort)."),
    threads: int = typer.Option(None),
    force: bool = typer.Option(False),
) -> None:
    """Run steps locally, in order, for one sample."""
    cfg = load_config(config)
    threads = threads or cfg.threads
    sd = outdir / sample
    todo = {
        "qc": lambda: qc_step.run_qc(cfg, sd, reads or None, threads, force),
        "hostdeplete": lambda: hostdeplete_step.run_hostdeplete(cfg, sd, threads, force),
        "rca": lambda: rca_step.run_rca(cfg, sd, threads, force),
        "assemble": lambda: assemble_step.run_assemble(cfg, sd, threads, force),
        "classify": lambda: classify_step.run_classify(cfg, sd, threads, force),
        "species": lambda: species_step.run_species(cfg, sd, threads, force),
        "satellite": lambda: satellite_step.run_satellite(cfg, sd, threads, force),
        "abundance": lambda: abundance_step.run_abundance(cfg, sd, threads, force),
        "diversity": lambda: diversity_step.run_diversity(cfg, outdir, [sample], threads, force),
        "genbank": lambda: genbank_step.run_genbank(cfg, sd, force),
        "variety": lambda: variety_step.run_variety(cfg, sd, threads, force),
    }
    wanted = steps.split(",")
    if set(wanted) - set(STEPS):
        raise typer.BadParameter(f"unknown steps: {sorted(set(wanted) - set(STEPS))}")
    for step in sorted(wanted, key=STEPS.index):
        todo[step]()


@app.command("download-db")
def download_db(
    outdir: Path = typer.Option(Path("db"), help="Where to write the database."),
    threads: int = typer.Option(4),
    keep_raw: bool = typer.Option(False, help="Keep downloaded files (they are ~250 MB)."),
) -> None:
    """Download NCBI RefSeq viral proteins + taxdump and build the DIAMOND database."""
    dmnd = download_db_step.run_download_db(outdir, threads)
    if not keep_raw:
        import shutil
        shutil.rmtree(outdir / "raw")
    typer.echo(f"Add to your config:\ndatabases:\n  viral_protein: {dmnd.resolve()}\n"
               f"  taxdump: {(outdir / 'taxdump').resolve()}")


@app.command()
def species(
    config: Path = typer.Option(..., exists=True),
    outdir: Path = typer.Option(Path("results")),
    sample: str = typer.Option(...),
    threads: int = typer.Option(None),
    force: bool = typer.Option(False, help="Re-run even if already done."),
) -> None:
    """Species assignment of contigs by nucleotide identity to reference genomes."""
    cfg = load_config(config)
    species_step.run_species(cfg, outdir / sample, threads or cfg.threads, force)


@app.command()
def satellite(
    config: Path = typer.Option(..., exists=True),
    outdir: Path = typer.Option(Path("results")),
    sample: str = typer.Option(...),
    threads: int = typer.Option(None),
    force: bool = typer.Option(False, help="Re-run even if already done."),
) -> None:
    """Detect satellite DNAs (alpha-, betasatellites ...) among the contigs."""
    cfg = load_config(config)
    satellite_step.run_satellite(cfg, outdir / sample, threads or cfg.threads, force)


@app.command("download-refs")
def download_refs(
    taxon: str = typer.Option(..., help="Genus or other taxon name, e.g. Begomovirus."),
    outdir: Path = typer.Option(Path("db/refs")),
    extra: list[str] = typer.Option([], help="Non-RefSeq accession to add, ACC or ACC=SPECIES_TAXID "
                                    "(repeatable), e.g. AJ717531.1=62079."),
) -> None:
    """Download RefSeq reference genomes of a taxon for the species step."""
    acc = {a.partition("=")[0]: int(a.partition("=")[2] or 0) or None for a in extra}
    fasta = download_refs_step.run_download_refs(taxon, outdir, acc)
    typer.echo(f"Add to your config:\nspecies:\n  reference: {fasta.resolve()}\n  taxon: {taxon}")
