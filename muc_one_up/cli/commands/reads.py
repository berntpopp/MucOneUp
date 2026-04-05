"""Reads commands -- simulate sequencing reads from FASTA files."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import click

from .._common import require_config
from ..error_handling import cli_error_handler
from ..options import shared_read_options

# ============================================================================
# Shared batch helper
# ============================================================================


def _run_batch_simulation(
    config: dict[str, Any],
    input_fastas: tuple[str, ...],
    out_dir: str,
    out_base: str | None,
    suffix: str,
    simulator_label: str,
    track_read_source: bool,
) -> None:
    """Run read simulation over one or more FASTA files.

    Shared loop for illumina, ont, and pacbio commands. Handles output naming,
    source tracker reconstruction, and pipeline dispatch.
    """
    from ...read_simulation import simulate_reads as simulate_reads_pipeline
    from ...read_simulator.output_config import OutputConfig

    # Warn if --out-base provided for multiple files
    if len(input_fastas) > 1 and out_base:
        logging.warning(
            "--out-base '%s' will be used for all %d files. "
            "Consider omitting --out-base for auto-generated names.",
            out_base,
            len(input_fastas),
        )

    total_files = len(input_fastas)
    logging.info("Processing %d FASTA file(s) for %s read simulation", total_files, simulator_label)

    for idx, input_fasta in enumerate(input_fastas, start=1):
        # Determine output base name
        if out_base:
            actual_out_base = f"{out_base}_{idx:03d}" if total_files > 1 else out_base
        else:
            actual_out_base = None

        output_config = OutputConfig.from_input_fasta(
            input_fa=input_fasta,
            out_dir=out_dir,
            out_base=actual_out_base,
            suffix=suffix,
        )

        logging.info(
            "[%d/%d] Simulating %s reads: %s -> %s",
            idx,
            total_files,
            simulator_label,
            input_fasta,
            output_config.out_base,
        )

        # Build source tracker from companion files if requested
        source_tracker = None
        if track_read_source:
            from ...read_simulator.source_tracking import ReadSourceTracker

            # Derive stats path: input may be "base.001.simulated.fa" but stats
            # file is "base.001.simulation_stats.json" (no ".simulated" component).
            input_p = Path(input_fasta)
            if input_p.stem.endswith(".simulated"):
                stats_path = str(
                    input_p.with_name(
                        input_p.stem.removesuffix(".simulated") + ".simulation_stats.json"
                    )
                )
            else:
                stats_path = str(input_p.with_suffix(".simulation_stats.json"))
            source_tracker = ReadSourceTracker.from_companion_files(stats_path)
            if source_tracker is None:
                logging.warning("Could not reconstruct read source tracker from %s", stats_path)
            else:
                coord_map_path = output_config.derive_path_str("_repeat_coordinates.tsv")
                source_tracker.write_coordinate_map(coord_map_path)
                logging.info("Repeat coordinate map written: %s", coord_map_path)

        # Run simulation for this file
        simulate_reads_pipeline(
            config,
            input_fasta,
            source_tracker=source_tracker,
            output_config=output_config,
        )

    logging.info("%s read simulation completed for all %d file(s).", simulator_label, total_files)


def _setup_read_config(
    config: dict[str, Any],
    simulator_type: str,
    coverage: int | None,
    seed: int | None,
    seed_config_key: str = "read_simulation",
) -> None:
    """Initialize read_simulation config section with shared defaults.

    Args:
        config: Full configuration dict (mutated in place)
        simulator_type: One of "illumina", "ont", "pacbio"
        coverage: CLI --coverage value (None if not provided)
        seed: CLI --seed value (None if not provided)
        seed_config_key: Config section to store seed in
    """
    if "read_simulation" not in config:
        config["read_simulation"] = {}
    config["read_simulation"]["simulator"] = simulator_type

    if coverage is not None:
        config["read_simulation"]["coverage"] = coverage
    elif "coverage" not in config.get("read_simulation", {}):
        config["read_simulation"]["coverage"] = 30

    if seed is not None:
        if seed_config_key not in config:
            config[seed_config_key] = {}
        config[seed_config_key]["seed"] = seed
        logging.info(f"Using random seed: {seed} (results will be reproducible)")


def _apply_pacbio_params(
    config: dict[str, Any],
    model_type: str | None,
    model_file: str | None,
    seed: int | None,
    threads: int | None = None,
    pass_num: int | None = None,
    min_passes: int | None = None,
    min_rq: float | None = None,
) -> None:
    """Apply PacBio CLI overrides and validate required params.

    Used by both pacbio() and amplicon() commands.
    Mutates config in place.

    Raises:
        click.ClickException: If required params are missing after merge.
    """
    if "pacbio_params" not in config:
        config["pacbio_params"] = {}

    params = config["pacbio_params"]
    overrides = {
        "model_type": model_type,
        "model_file": model_file,
        "pass_num": pass_num,
        "min_passes": min_passes,
        "min_rq": min_rq,
    }
    if threads is not None:
        overrides["threads"] = threads

    for key, value in overrides.items():
        if value is not None:
            params[key] = value

    if seed is not None:
        params["seed"] = seed
        logging.info("Using random seed: %d (results will be reproducible)", seed)

    required = ["model_type", "model_file"]
    missing = [p for p in required if p not in params]
    if missing:
        raise click.ClickException(
            f"Missing required PacBio parameters: {', '.join(missing)}. "
            f"Provide via CLI options or config.json pacbio_params section."
        )


# ============================================================================
# READS Command Group - Pure Read Simulation Utilities
# ============================================================================


@click.group()
def reads():
    """Read simulation utilities.

    \b
    Simulate reads from ANY FASTA file.
    Works with MucOneUp outputs or external sequences.
    Requires --config at the top level.
    """


@reads.command()
@click.option(
    "--threads",
    type=int,
    default=8,
    show_default=True,
    help="Number of threads.",
)
@shared_read_options
@click.pass_context
@cli_error_handler
def illumina(ctx, input_fastas, out_dir, out_base, coverage, threads, seed, track_read_source):
    """Simulate Illumina short reads from one or more FASTA files.

    Supports batch processing following Unix philosophy:
    - Single file: muconeup reads illumina file.fa --out-base reads
    - Multiple files: muconeup reads illumina file1.fa file2.fa file3.fa
    - Glob pattern: muconeup reads illumina *.simulated.fa

    When processing multiple files, --out-base is auto-generated from input
    filenames unless explicitly provided (which applies to all files).

    \b
    Examples:
      # Single file with custom output name
      muconeup --config X reads illumina sample.001.fa --out-base my_reads

      # Multiple files (auto-generated output names)
      muconeup --config X reads illumina sample.001.fa sample.002.fa

      # Glob pattern (shell expands)
      muconeup --config X reads illumina sample.*.simulated.fa

      # Compose with shell (Unix philosophy)
      for f in *.fa; do muconeup --config X reads illumina "$f"; done
    """
    require_config(ctx)

    from ...config import load_config_raw

    config = load_config_raw(str(ctx.obj["config_path"]))
    _setup_read_config(config, "illumina", coverage, seed)
    config["read_simulation"]["threads"] = threads

    _run_batch_simulation(
        config, input_fastas, out_dir, out_base, "_reads", "Illumina", track_read_source
    )


@reads.command()
@click.option(
    "--min-read-length",
    type=int,
    default=100,
    show_default=True,
    help="Minimum read length.",
)
@shared_read_options
@click.pass_context
@cli_error_handler
def ont(ctx, input_fastas, out_dir, out_base, coverage, min_read_length, seed, track_read_source):
    """Simulate Oxford Nanopore long reads from one or more FASTA files.

    Supports batch processing following Unix philosophy:
    - Single file: muconeup reads ont file.fa --out-base reads
    - Multiple files: muconeup reads ont file1.fa file2.fa file3.fa
    - Glob pattern: muconeup reads ont *.simulated.fa

    When processing multiple files, --out-base is auto-generated from input
    filenames unless explicitly provided (which applies to all files).

    \b
    Examples:
      # Single file with custom output name
      muconeup --config X reads ont sample.001.fa --out-base my_reads

      # Multiple files (auto-generated output names)
      muconeup --config X reads ont sample.001.fa sample.002.fa

      # Glob pattern (shell expands)
      muconeup --config X reads ont sample.*.simulated.fa
    """
    require_config(ctx)

    from ...config import load_config_raw

    config = load_config_raw(str(ctx.obj["config_path"]))
    _setup_read_config(config, "ont", coverage, seed, seed_config_key="nanosim_params")

    from typing import cast

    from ...type_defs import NanosimConfig

    if "nanosim_params" not in config:
        config["nanosim_params"] = {}
    ns = cast(NanosimConfig, config["nanosim_params"])
    ns["min_read_length"] = min_read_length
    # Propagate CLI coverage to nanosim_params where the ONT backend reads it.
    # Only overwrite if CLI provided a value or config lacks ONT-specific coverage.
    if coverage is not None:
        ns["coverage"] = coverage
    elif "coverage" not in ns:
        ns["coverage"] = config["read_simulation"]["coverage"]

    _run_batch_simulation(
        config, input_fastas, out_dir, out_base, "_ont_reads", "ONT", track_read_source
    )


@reads.command()
@click.option(
    "--threads",
    type=int,
    default=4,
    show_default=True,
    help="Number of threads.",
)
@click.option(
    "--model-file",
    type=click.Path(exists=True, dir_okay=False),
    default=None,
    help="Path to pbsim3 model file (overrides config if provided).",
)
@click.option(
    "--model-type",
    type=click.Choice(["qshmm", "errhmm"]),
    default=None,
    help="pbsim3 model type (overrides config if provided).",
)
@click.option(
    "--min-rq",
    type=float,
    default=None,
    help="Minimum predicted read quality (RQ) score (0.0-1.0). 0.99=Q20 (standard HiFi), 0.999=Q30.",
)
@click.option(
    "--min-passes",
    type=int,
    default=None,
    help="Minimum passes required for CCS HiFi consensus (>=1, overrides config if provided).",
)
@click.option(
    "--pass-num",
    type=int,
    default=None,
    help="Number of passes per molecule for multi-pass CLR simulation (>=2, overrides config if provided).",
)
@shared_read_options
@click.pass_context
@cli_error_handler
def pacbio(
    ctx,
    input_fastas,
    out_dir,
    out_base,
    coverage,
    pass_num,
    min_passes,
    min_rq,
    model_type,
    model_file,
    threads,
    seed,
    track_read_source,
):
    """Simulate PacBio HiFi reads from one or more FASTA files.

    Supports batch processing following Unix philosophy:
    - Single file: muconeup reads pacbio file.fa --model-file X.model --out-base reads
    - Multiple files: muconeup reads pacbio file1.fa file2.fa --model-file X.model
    - Glob pattern: muconeup reads pacbio *.simulated.fa --model-file X.model

    When processing multiple files, --out-base is auto-generated from input
    filenames unless explicitly provided (which applies to all files).

    \b
    Workflow:
      1. Multi-pass CLR simulation (pbsim3)
      2. HiFi consensus generation (CCS)
      3. Read alignment (minimap2 with map-hifi preset)

    \b
    Examples:
      # Single file with standard HiFi settings (Q20)
      muconeup --config X reads pacbio sample.001.fa \\
        --model-file /models/QSHMM-SEQUEL.model \\
        --out-base my_hifi

      # Multiple files with high-accuracy HiFi (Q30)
      muconeup --config X reads pacbio sample.*.fa \\
        --model-file /models/QSHMM-SEQUEL.model \\
        --min-rq 0.999 --min-passes 5

      # Ultra-deep coverage simulation
      muconeup --config X reads pacbio sample.fa \\
        --model-file /models/QSHMM-SEQUEL.model \\
        --coverage 100 --pass-num 5

    \b
    Model Files:
      Download from: https://github.com/yukiteruono/pbsim3/tree/master/data
      - QSHMM-SEQUEL.model: Sequel II chemistry
      - QSHMM-RSII.model: RS II chemistry
      - ERRHMM-SEQUEL.model: Alternative error model

    \b
    Quality Control:
      - pass_num >=2 required for multi-pass (>=3 recommended)
      - min_passes controls CCS stringency (higher = better quality, lower yield)
      - min_rq=0.99 is Q20 (standard HiFi threshold)
      - min_rq=0.999 is Q30 (ultra-high accuracy)
    """
    require_config(ctx)

    from ...config import load_config_raw

    config = load_config_raw(str(ctx.obj["config_path"]))

    _setup_read_config(config, "pacbio", coverage, seed, seed_config_key="pacbio_params")
    _apply_pacbio_params(
        config,
        model_type,
        model_file,
        seed,
        threads=threads if threads != 4 else None,
        pass_num=pass_num,
        min_passes=min_passes,
        min_rq=min_rq,
    )

    _run_batch_simulation(
        config, input_fastas, out_dir, out_base, "_pacbio_hifi", "PacBio HiFi", track_read_source
    )


@reads.command()
@click.option(
    "--model-file",
    type=click.Path(exists=True, dir_okay=False),
    default=None,
    help="Path to pbsim3 model file (overrides config if provided).",
)
@click.option(
    "--model-type",
    type=click.Choice(["qshmm", "errhmm"]),
    default=None,
    help="pbsim3 model type (overrides config if provided).",
)
@click.option(
    "--pcr-preset",
    type=click.Choice(["default", "no_bias"]),
    default=None,
    help="PCR bias preset profile (default: from config or 'default').",
)
@click.option(
    "--stochastic-pcr",
    is_flag=True,
    default=False,
    help="Enable stochastic PCR bias (Galton-Watson branching process).",
)
@click.option(
    "--platform",
    type=click.Choice(["pacbio", "ont"]),
    default="pacbio",
    show_default=True,
    help="Sequencing platform for amplicon simulation.",
)
@shared_read_options
@click.pass_context
@cli_error_handler
def amplicon(
    ctx,
    input_fastas,
    out_dir,
    out_base,
    coverage,
    model_type,
    model_file,
    pcr_preset,
    stochastic_pcr,
    platform,
    seed,
    track_read_source,
):
    """Simulate amplicon reads from one or more FASTA files.

    Supports PacBio (default) and ONT platforms. PacBio uses multi-pass
    CCS consensus; ONT uses single-pass pbsim3 with map-ont alignment.

    \b
    Note: --coverage specifies the total number of template molecules
    (before CCS filtering for PacBio). Final HiFi read count may be
    lower due to CCS quality filtering (min-rq, min-passes). For diploid
    inputs the total is split between alleles by the PCR bias model.

    \b
    Pipeline (PacBio):
      1. Extract amplicon region per haplotype (primer-based)
      2. Apply PCR length bias to determine per-allele coverage
      3. Simulate full-length reads (pbsim3 --strategy templ)
      4. Generate HiFi consensus (CCS)
      5. Align to reference (minimap2 map-hifi preset)

    \b
    Pipeline (ONT):
      1-3. Same extraction and PCR bias stages
      4. Simulate single-pass reads (pbsim3 --strategy templ, pass_num=1)
      5. BAM to FASTQ (no CCS)
      6. Align to reference (minimap2 map-ont preset)

    \b
    Examples:
      # Basic PacBio amplicon simulation
      muconeup --config X reads amplicon sample.simulated.fa \\
        --model-file /models/QSHMM-SEQUEL.model

      # ONT amplicon simulation
      muconeup --config X reads amplicon --platform ont sample.fa \\
        --model-file /models/ERRHMM-ONT.model

      # High coverage with stochastic PCR bias
      muconeup --config X reads amplicon sample.fa \\
        --model-file /models/QSHMM-SEQUEL.model \\
        --coverage 1000 --stochastic-pcr --seed 42
    """
    require_config(ctx)

    # Reject --track-read-source early
    if track_read_source:
        raise click.ClickException(
            "Read source tracking is not yet supported for amplicon simulation. "
            "Remove --track-read-source to proceed."
        )

    from ...config import load_config_raw

    config = load_config_raw(str(ctx.obj["config_path"]))
    _setup_read_config(config, "amplicon", coverage, seed)

    if platform == "ont":
        config["read_simulation"]["simulator"] = "ont-amplicon"
    config["read_simulation"]["assay_type"] = "amplicon"

    # Ensure amplicon_params exists
    if "amplicon_params" not in config:
        raise click.ClickException(
            "Missing amplicon_params section in config. "
            "Add forward_primer and reverse_primer to config.json."
        )

    _apply_pacbio_params(config, model_type, model_file, seed)

    # Apply PCR bias overrides
    if "pcr_bias" not in config["amplicon_params"]:
        config["amplicon_params"]["pcr_bias"] = {}
    if pcr_preset is not None:
        config["amplicon_params"]["pcr_bias"]["preset"] = pcr_preset
    if stochastic_pcr:
        config["amplicon_params"]["pcr_bias"]["stochastic"] = True

    _run_batch_simulation(
        config,
        input_fastas,
        out_dir,
        out_base,
        "_amplicon",
        "Amplicon",
        track_read_source=False,
    )
