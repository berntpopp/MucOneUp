"""Analyze commands -- analysis utilities for FASTA and VNTR data."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

import click

from .._common import require_config
from ..error_handling import cli_error_handler
from ..outputs import generate_output_base

# ============================================================================
# ANALYZE Command Group - Pure Analysis Utilities
# ============================================================================


@click.group()
def analyze():
    """Analysis utilities.

    \b
    Single Responsibility: Analyze ANY FASTA file.
    Works with MucOneUp outputs or external sequences.
    """


@analyze.command()
@click.argument(
    "input_fastas", nargs=-1, required=True, type=click.Path(exists=True, dir_okay=False)
)
@click.option(
    "--out-dir",
    default=".",
    show_default=True,
    type=click.Path(file_okay=False),
    help="Output folder.",
)
@click.option(
    "--out-base",
    default=None,
    help="Base name for output files (auto-generated if processing multiple files).",
)
@click.option(
    "--orf-min-aa",
    type=int,
    default=100,
    show_default=True,
    help="Minimum ORF length in amino acids.",
)
@click.option(
    "--orf-aa-prefix",
    default=None,
    help="Filter ORFs by prefix (e.g., MTSSV).",
)
@click.pass_context
@cli_error_handler
def orfs(ctx, input_fastas, out_dir, out_base, orf_min_aa, orf_aa_prefix):
    """Predict ORFs and detect toxic protein features from one or more FASTA files.

    Supports batch processing following Unix philosophy:
    - Single file: muconeup analyze orfs file.fa --out-base analysis
    - Multiple files: muconeup analyze orfs file1.fa file2.fa file3.fa
    - Glob pattern: muconeup analyze orfs *.simulated.fa

    When processing multiple files, --out-base is auto-generated from input
    filenames unless explicitly provided (which applies to all files).

    \b
    Examples:
      # Single file with custom output name
      muconeup --config X analyze orfs sample.001.fa --out-base my_analysis

      # Multiple files (auto-generated output names)
      muconeup --config X analyze orfs sample.001.fa sample.002.fa

      # Glob pattern (shell expands)
      muconeup --config X analyze orfs sample.*.simulated.fa
    """
    # Validate config is provided
    require_config(ctx)

    # Load config once for toxic protein detection (DRY principle)
    from ...config import load_config_raw

    config = load_config_raw(str(ctx.obj["config_path"]))
    ref_assembly = config.get("reference_assembly", "hg38")
    assembly_constants = config.get("constants", {}).get(ref_assembly, {})
    left_const = assembly_constants.get("left")
    right_const = assembly_constants.get("right")

    # Warn if --out-base provided for multiple files
    if len(input_fastas) > 1 and out_base:
        logging.warning(
            "--out-base '%s' will be used for all %d files. "
            "Consider omitting --out-base for auto-generated names.",
            out_base,
            len(input_fastas),
        )

    # Process each file (KISS principle - simple iteration)
    total_files = len(input_fastas)
    logging.info("Processing %d FASTA file(s) for ORF prediction", total_files)

    for idx, input_fasta in enumerate(input_fastas, start=1):
        from ..analysis import run_orf_analysis_standalone

        if out_base:
            actual_out_base = f"{out_base}_{idx:03d}" if total_files > 1 else out_base
        else:
            actual_out_base = generate_output_base(Path(input_fasta), "_orfs")

        logging.info(
            "[%d/%d] Running ORF prediction: %s -> %s",
            idx,
            total_files,
            input_fasta,
            Path(out_dir) / f"{actual_out_base}.orfs.fa",
        )

        run_orf_analysis_standalone(
            input_fasta=input_fasta,
            out_dir=str(out_dir),
            out_base=actual_out_base,
            orf_min_aa=orf_min_aa,
            orf_aa_prefix=orf_aa_prefix,
            left_const=left_const,
            right_const=right_const,
        )

    logging.info("ORF prediction completed for all %d file(s).", total_files)


@analyze.command()
@click.argument(
    "input_fastas", nargs=-1, required=True, type=click.Path(exists=True, dir_okay=False)
)
@click.option(
    "--out-dir",
    default=".",
    show_default=True,
    type=click.Path(file_okay=False),
    help="Output folder.",
)
@click.option(
    "--out-base",
    default=None,
    help="Base name for output files (auto-generated if processing multiple files).",
)
@click.pass_context
@cli_error_handler
def stats(ctx, input_fastas, out_dir, out_base):
    """Generate basic sequence statistics from one or more FASTA files.

    Supports batch processing following Unix philosophy:
    - Single file: muconeup analyze stats file.fa --out-base stats
    - Multiple files: muconeup analyze stats file1.fa file2.fa file3.fa
    - Glob pattern: muconeup analyze stats *.simulated.fa

    When processing multiple files, --out-base is auto-generated from input
    filenames unless explicitly provided (which applies to all files).

    \b
    Examples:
      # Single file with custom output name
      muconeup --config X analyze stats sample.001.fa --out-base my_stats

      # Multiple files (auto-generated output names)
      muconeup --config X analyze stats sample.001.fa sample.002.fa

      # Glob pattern (shell expands)
      muconeup --config X analyze stats sample.*.simulated.fa
    """
    # Warn if --out-base provided for multiple files
    if len(input_fastas) > 1 and out_base:
        logging.warning(
            "--out-base '%s' will be used for all %d files. "
            "Consider omitting --out-base for auto-generated names.",
            out_base,
            len(input_fastas),
        )

    # Process each file (KISS principle - simple iteration)
    total_files = len(input_fastas)
    logging.info("Processing %d FASTA file(s) for statistics generation", total_files)

    for idx, input_fasta in enumerate(input_fastas, start=1):
        # Determine output base name
        if out_base:
            # User provided: use as-is (or append index for multiple files)
            actual_out_base = f"{out_base}_{idx:03d}" if total_files > 1 else out_base
        else:
            # Auto-generate from input filename
            actual_out_base = generate_output_base(Path(input_fasta), "_stats")

        logging.info("[%d/%d] Generating statistics: %s", idx, total_files, input_fasta)

        # Simple FASTA parsing
        sequences: list[tuple[str, str]] = []
        current_header: str | None = None
        current_seq: list[str] = []

        with Path(input_fasta).open() as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if current_header is not None:
                        sequences.append((current_header, "".join(current_seq)))
                    current_header = line[1:]
                    current_seq = []
                elif line:
                    current_seq.append(line)

            if current_header is not None:
                sequences.append((current_header, "".join(current_seq)))

        # Compute stats
        stats_data: dict[str, Any] = {
            "input_file": str(input_fasta),
            "num_sequences": len(sequences),
            "sequences": [],
        }

        for header, sequence in sequences:
            gc_count = sequence.upper().count("G") + sequence.upper().count("C")
            gc_content = 100.0 * gc_count / len(sequence) if sequence else 0.0

            stats_data["sequences"].append(
                {"header": header, "length": len(sequence), "gc_content": round(gc_content, 2)}
            )

        stats_file = Path(out_dir) / f"{actual_out_base}.basic_stats.json"
        with stats_file.open("w") as f:
            json.dump(stats_data, f, indent=4)

        logging.info("Statistics written: %s", stats_file)

    logging.info("Statistics generation completed for all %d file(s).", total_files)


# ============================================================================
# ANALYZE vntr-stats - VNTR Transition Probability Analysis
# ============================================================================


@analyze.command("vntr-stats")
@click.argument("input_file", type=click.Path(exists=True, dir_okay=False))
@click.option(
    "--structure-column",
    default="vntr",
    show_default=True,
    help="Column name (if header) or 0-based index containing VNTR structure.",
)
@click.option(
    "--delimiter",
    default="\t",
    show_default=True,
    help="Field delimiter for input file.",
)
@click.option(
    "--header",
    is_flag=True,
    help="Specify if input file has header row.",
)
@click.option(
    "--output",
    "-o",
    type=click.Path(),
    help="Output JSON file (default: stdout).",
)
@click.pass_context
@cli_error_handler
def vntr_stats(ctx, input_file, structure_column, delimiter, header, output):
    """Analyze VNTR structures and compute transition probabilities.

    Processes a CSV/TSV file containing VNTR structures, calculates statistics
    (min/max/mean/median repeat units), and builds a transition probability
    matrix showing the likelihood of each repeat unit following another.

    The analysis removes duplicate VNTR structures and includes an "END" state
    representing sequence termination. Unknown repeat tokens (not in config)
    trigger warnings but don't cause failure.

    \b
    Examples:
      # Analyze example VNTR database
      muconeup --config X analyze vntr-stats data/examples/vntr_database.tsv --header

      # Use custom column and save to file
      muconeup --config X analyze vntr-stats data.csv \\
        --delimiter "," --structure-column "sequence" --output stats.json

      # Column index without header
      muconeup --config X analyze vntr-stats data.tsv --structure-column 3

      # Pipe to jq for filtering
      muconeup --config X analyze vntr-stats data/examples/vntr_database.tsv \\
        --header | jq '.mean_repeats'

    \b
    Output JSON contains:
      - Statistics: min/max/mean/median repeat counts
      - Probabilities: State transition matrix (including END state)
      - Repeats: Known repeat dictionary from config
    """
    # Validate config is provided
    require_config(ctx)

    # Load config for known repeats (DRY principle - reuse existing pattern)
    from ...config import load_config_raw

    config = load_config_raw(str(ctx.obj["config_path"]))

    known_repeats = config.get("repeats", {})
    if not known_repeats:
        logging.error("Configuration file does not contain 'repeats' key")
        ctx.exit(1)

    logging.info("Analyzing VNTR structures from %s", input_file)

    # Import analysis module
    from ...analysis.vntr_statistics import analyze_vntr_sequences

    # Run analysis
    with Path(input_file).open() as f:
        results = analyze_vntr_sequences(f, structure_column, header, known_repeats, delimiter)

    # Build output (matches original helper behavior)
    output_data = {
        "min_repeats": results["min_repeats"],
        "max_repeats": results["max_repeats"],
        "mean_repeats": results["mean_repeats"],
        "median_repeats": results["median_repeats"],
        "repeats": known_repeats,
        "probabilities": results["probabilities"],
    }

    # Write output
    if output:
        output_path = Path(output)
        with output_path.open("w") as out:
            json.dump(output_data, out, indent=2)
        logging.info("VNTR statistics written to %s", output)
    else:
        # Output to stdout (pipeable)
        click.echo(json.dumps(output_data, indent=2))

    logging.info(
        "Analysis complete: %d unique structures, %d repeat types",
        len(results["probabilities"]),
        len([t for t in results["probabilities"] if t != "END"]),
    )


@analyze.command()
@click.argument("input_fasta", type=click.Path(exists=True, dir_okay=False))
@click.option(
    "--mutation",
    required=True,
    help="Mutation name to validate (e.g., 'dupC').",
)
@click.option(
    "--output",
    type=click.Path(),
    help="Output JSON file for validation results (prints to stdout if not specified).",
)
@click.pass_context
@cli_error_handler
def snapshot_validate(ctx, input_fasta, mutation, output):
    """Validate SNaPshot assay for MUC1 VNTR mutations.

    Simulates complete SNaPshot workflow: PCR -> MwoI digest -> extension -> detection.

    \b
    Examples:
      # Validate dupC mutation in a sample
      muconeup --config config.json analyze snapshot-validate sample.fa --mutation dupC

      # Save results to JSON
      muconeup --config config.json analyze snapshot-validate sample.fa --mutation dupC --output results.json
    """
    # Validate config is provided
    require_config(ctx)

    from Bio import SeqIO

    from ...analysis.snapshot_validator import SnapshotValidator
    from ...config import load_config

    # Load config
    config_path = Path(ctx.obj["config_path"])
    config = load_config(str(config_path))

    # Check if mutation is configured
    if "snapshot_validation" not in config:
        logging.error("Config missing 'snapshot_validation' section")
        ctx.exit(1)

    if mutation not in config["snapshot_validation"]:
        available = list(config["snapshot_validation"].keys())
        logging.error(
            "Mutation '%s' not found in config. Available: %s",
            mutation,
            available,
        )
        ctx.exit(1)

    # Initialize validator
    validator = SnapshotValidator(config, mutation)
    logging.info("SNaPshot validator initialized for mutation: %s", mutation)

    # Load input FASTA
    records = list(SeqIO.parse(input_fasta, "fasta"))
    if not records:
        logging.error("No sequences found in %s", input_fasta)
        ctx.exit(1)

    logging.info("Loaded %d haplotype(s) from %s", len(records), input_fasta)

    # Run validation on each haplotype
    all_results = {}
    for i, record in enumerate(records, 1):
        haplotype_name = record.id or f"haplotype_{i}"
        template = str(record.seq)

        logging.info(
            "Validating haplotype %d: %s (%d bp)",
            i,
            haplotype_name,
            len(template),
        )

        result = validator.validate_complete_workflow(template)

        all_results[haplotype_name] = {
            "mutation_detected": result["mutation_detected"],
            "pcr_products": result["pcr_results"]["product_count"],
            "digest_survivors": result["digest_results"]["survivor_count"],
            "expected_fluorescence": result.get("expected_fluorescence", ""),
            "summary": result["summary"],
        }

        # Add SNaPshot details if mutation detected
        if result["mutation_detected"] and result.get("snapshot_results"):
            snapshot_info = []
            for snap in result["snapshot_results"]:
                if snap["extension"].get("binds"):
                    snapshot_info.append(
                        {
                            "fluorophore": snap["extension"].get("fluorophore"),
                            "dye": snap["extension"].get("fluorophore_dye"),
                            "next_base": snap["extension"].get("next_base"),
                        }
                    )
            all_results[haplotype_name]["snapshot_details"] = snapshot_info

        logging.info("  Result: %s", result["summary"])

    # Output results
    output_data = {
        "mutation": mutation,
        "input_file": str(input_fasta),
        "haplotypes": all_results,
        "overall_detection": any(r["mutation_detected"] for r in all_results.values()),
    }

    if output:
        output_path = Path(output)
        with output_path.open("w") as f:
            json.dump(output_data, f, indent=2)
        logging.info("Validation results written to %s", output)
    else:
        click.echo(json.dumps(output_data, indent=2))

    # Exit code: 0 if mutation detected, 1 if not
    if output_data["overall_detection"]:
        logging.info("Mutation %s DETECTED in sample", mutation.upper())
        ctx.exit(0)
    else:
        logging.info("Mutation %s NOT detected in sample", mutation.upper())
        ctx.exit(1)
