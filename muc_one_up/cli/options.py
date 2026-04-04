"""Typed option dataclasses for CLI commands."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import click


@dataclass
class SimulationOptions:
    """Typed options for the simulate command and orchestration pipeline.

    Replaces the SimpleNamespace created by _make_args_namespace().
    All orchestration functions access these fields by attribute name.
    """

    # Required
    config: str
    out_base: str
    out_dir: str

    # Simulation parameters
    num_haplotypes: int = 2
    seed: int | None = None
    reference_assembly: str | None = "hg38"

    # Output control
    output_structure: bool = False

    # Mode selection
    fixed_lengths: list[str] | None = None
    input_structure: str | None = None
    simulate_series: int | None = None

    # Mutation config
    mutation_name: str | None = None
    mutation_targets: list[str] | None = None

    # SNP config
    snp_input_file: str | None = None
    random_snps: bool = False
    random_snp_density: float = 0.001
    random_snp_output_file: str | None = None
    random_snp_region: str = "constants_only"
    random_snp_haplotypes: str = "all"

    # Pipeline flags (set to None/False for pure simulate)
    simulate_reads: str | None = None
    output_orfs: bool = False
    orf_min_aa: int = 100
    orf_aa_prefix: str | None = None
    track_read_source: bool = False

    @classmethod
    def from_click_kwargs(cls, config_path: str, kwargs: dict) -> SimulationOptions:
        """Construct from Click command kwargs.

        Converts Click tuples to lists and empty tuples to None.
        """
        fixed = list(kwargs["fixed_lengths"]) if kwargs["fixed_lengths"] else None
        targets = list(kwargs["mutation_targets"]) if kwargs["mutation_targets"] else None

        return cls(
            config=config_path,
            out_base=kwargs["out_base"],
            out_dir=kwargs["out_dir"],
            num_haplotypes=kwargs["num_haplotypes"],
            seed=kwargs["seed"],
            reference_assembly=kwargs["reference_assembly"],
            output_structure=kwargs["output_structure"],
            fixed_lengths=fixed,
            input_structure=kwargs["input_structure"],
            simulate_series=kwargs["simulate_series"],
            mutation_name=kwargs["mutation_name"],
            mutation_targets=targets,
            snp_input_file=kwargs["snp_input_file"],
            random_snps=kwargs["random_snps"],
            random_snp_density=kwargs["random_snp_density"],
            random_snp_output_file=kwargs["random_snp_output_file"],
            random_snp_region=kwargs["random_snp_region"],
            random_snp_haplotypes=kwargs["random_snp_haplotypes"],
            track_read_source=kwargs.get("track_read_source", False),
        )


def shared_read_options(func: Any) -> Any:
    """Apply common Click options shared across all read simulation subcommands.

    Adds: input_fastas argument, --out-dir, --out-base, --coverage, --seed,
    --track-read-source.

    Options are applied in reverse display order so that after Click's
    internal reversal the help text shows them as:
    --out-dir, --out-base, --coverage, --seed, --track-read-source.

    Command-specific options stacked above this decorator will appear
    before the shared options in ``--help`` output.
    """
    func = click.option(
        "--track-read-source",
        is_flag=True,
        default=False,
        help="Generate read source tracking manifest and coordinate map alongside simulated reads.",
    )(func)
    func = click.option(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducibility (same seed = identical reads).",
    )(func)
    func = click.option(
        "--coverage",
        type=int,
        default=None,
        help="Target sequencing coverage (overrides config if provided, defaults to config value or 30x).",
    )(func)
    func = click.option(
        "--out-base",
        default=None,
        help="Base name for output files (auto-generated if processing multiple files).",
    )(func)
    func = click.option(
        "--out-dir",
        default=".",
        show_default=True,
        type=click.Path(file_okay=False),
        help="Output folder.",
    )(func)
    func = click.argument(
        "input_fastas",
        nargs=-1,
        required=True,
        type=click.Path(exists=True, dir_okay=False),
    )(func)
    return func
