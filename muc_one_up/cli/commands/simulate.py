"""Simulate command -- generate MUC1 VNTR diploid haplotypes."""

from __future__ import annotations

import logging
from contextlib import nullcontext

import click

from .._common import require_config
from ..config import (
    determine_simulation_mode,
    process_mutation_config,
    setup_configuration,
)
from ..error_handling import cli_error_handler
from ..options import SimulationOptions
from ..orchestration import run_single_simulation_iteration


@click.command()
@click.option(
    "--out-base",
    default="muc1_simulated",
    show_default=True,
    help="Base name for output files.",
)
@click.option(
    "--out-dir",
    default=".",
    show_default=True,
    type=click.Path(file_okay=False),
    help="Output folder.",
)
@click.option(
    "--num-haplotypes",
    default=2,
    show_default=True,
    type=int,
    help="Number of haplotypes to simulate.",
)
@click.option(
    "--seed",
    type=int,
    help="Random seed for reproducibility.",
)
@click.option(
    "--reference-assembly",
    type=click.Choice(["hg19", "hg38"]),
    help="Reference assembly (overrides config).",
)
@click.option(
    "--output-structure",
    is_flag=True,
    help="Write VNTR structure file.",
)
# Length Options
@click.option(
    "--fixed-lengths",
    multiple=True,
    type=str,
    help="Fixed VNTR lengths or ranges (e.g., '60' or '20-40').",
)
@click.option(
    "--input-structure",
    type=click.Path(exists=True, dir_okay=False),
    help="Predefined VNTR structure file.",
)
@click.option(
    "--simulate-series",
    type=int,
    help="Series step size for fixed-length ranges.",
)
# Mutation Options
@click.option(
    "--mutation-name",
    help="Mutation name. Use 'normal,mutation' for dual simulation.",
)
@click.option(
    "--mutation-targets",
    multiple=True,
    help="Mutation targets as 'hap_idx,rep_idx' pairs (1-based).",
)
# SNP Options
@click.option(
    "--snp-input-file",
    type=click.Path(exists=True, dir_okay=False),
    help="TSV file with predefined SNPs.",
)
@click.option(
    "--random-snps",
    is_flag=True,
    help="Enable random SNP generation.",
)
@click.option(
    "--random-snp-density",
    type=float,
    help="SNP density per 1000 bp.",
)
@click.option(
    "--random-snp-output-file",
    type=str,
    help="Output file for random SNPs.",
)
@click.option(
    "--random-snp-region",
    type=click.Choice(["all", "constants_only", "vntr_only"]),
    default="constants_only",
    show_default=True,
    help="Region for random SNPs.",
)
@click.option(
    "--random-snp-haplotypes",
    type=click.Choice(["all", "1", "2"]),
    default="all",
    show_default=True,
    help="Haplotypes for random SNPs.",
)
@click.option(
    "--track-read-source",
    is_flag=True,
    default=False,
    help="Generate read source tracking manifest and coordinate map alongside simulated reads.",
)
@click.pass_context
@cli_error_handler
def simulate(ctx, **kwargs):
    """Generate MUC1 VNTR diploid haplotypes.

    \b
    Generates haplotype FASTA files only. For read simulation,
    pipe output to 'reads' commands (e.g., reads illumina).

    \b
    Output:
      - {out_base}.{iteration}.simulated.fa (haplotype sequences)
      - {out_base}.{iteration}.vntr_structure.txt (if --output-structure)
      - {out_base}.{iteration}.simulation_stats.json (statistics)

    \b
    Example:
      muconeup --config config.json simulate --out-base output
    """
    # Validate config is provided
    require_config(ctx)

    # Convert Click kwargs to typed options (replaces _make_args_namespace)
    args = SimulationOptions.from_click_kwargs(ctx.obj["config_path"], kwargs)

    # IMPORTANT: Disable pipeline options (simulate is PURE)
    args.simulate_reads = None
    args.output_orfs = False

    # Delegate to existing orchestration (SOLID - Dependency Inversion)
    config, out_dir, out_base = setup_configuration(args)
    simulation_configs, predefined_chains, structure_mutation_info = determine_simulation_mode(
        args, config
    )
    dual_mutation_mode, mutation_pair, _ = process_mutation_config(args, structure_mutation_info)

    # Run haplotype generation ONLY
    total_iterations = len(simulation_configs)

    # Show progress bar only for series mode (DRY - no code duplication)
    # Use nullcontext for single iterations to avoid if/else duplication
    progress_ctx = (
        click.progressbar(
            simulation_configs,
            label=f"Simulating {total_iterations} iterations",
            show_eta=True,
            show_pos=True,
        )
        if total_iterations > 1
        else nullcontext(simulation_configs)
    )

    with progress_ctx as configs:
        for sim_index, fixed_conf in enumerate(configs, start=1):
            run_single_simulation_iteration(
                args,
                config,
                out_dir,
                out_base,
                sim_index,
                fixed_conf,
                predefined_chains,
                dual_mutation_mode,
                mutation_pair,
                structure_mutation_info,
            )

    logging.info("Haplotype generation completed successfully.")
