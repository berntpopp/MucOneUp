"""
Click-based CLI for MucOneUp.

Design Philosophy: Clean Separation (Unix Philosophy)
Each command does ONE thing well:
- `simulate`: ONLY generates haplotypes
- `reads`: ONLY simulates reads from FASTA
- `analyze`: ONLY analyzes FASTA
- `pipeline`: Orchestrates the above for convenience

Follows SOLID principles: Single Responsibility, Dependency Inversion.
"""

import click

from ._common import (
    CONTEXT_SETTINGS,
    configure_logging,
    print_version_callback,
    validate_config_callback,
)

# ============================================================================
# Root CLI Group
# ============================================================================


@click.group(context_settings=CONTEXT_SETTINGS)
@click.option(
    "--version",
    "-V",
    is_flag=True,
    callback=print_version_callback,
    expose_value=False,
    is_eager=True,
    help="Show the version and exit.",
)
@click.option(
    "--config",
    callback=validate_config_callback,
    is_eager=True,  # Process early to enable resilient_parsing check for help/version
    type=click.Path(exists=True, dir_okay=False),
    help="Path to JSON configuration file.",
)
@click.option(
    "--log-level",
    type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL", "NONE"]),
    default="INFO",
    help="Set logging level.",
)
@click.option(
    "--verbose",
    "-v",
    is_flag=True,
    help="Enable verbose output (sets log level to DEBUG).",
)
@click.pass_context
def cli(ctx, config, log_level, verbose):
    """MucOneUp - MUC1 VNTR diploid reference simulator.

    \b
    Philosophy: Each command does ONE thing (Unix philosophy)
      simulate  - Generate haplotypes ONLY
      reads     - Simulate reads from FASTA (supports multiple files)
      analyze   - Analyze FASTA (ORFs, stats) (supports multiple files)

    \b
    Examples:
      # Generate haplotypes
      muconeup --config X simulate --out-base Y

      # Process single file
      muconeup --config X reads illumina Y.001.simulated.fa --out-base reads

      # Process multiple files (batch)
      muconeup --config X reads illumina Y.*.simulated.fa
      muconeup --config X analyze orfs Y.*.simulated.fa

      # Shell composition (Unix philosophy)
      muconeup --config X simulate --fixed-lengths 20-40 --simulate-series 1
      for f in Y.*.fa; do muconeup --config X reads illumina "$f"; done
    """
    # KISS: Simple precedence - verbose overrides log_level
    if verbose:
        log_level = "DEBUG"

    ctx.ensure_object(dict)
    ctx.obj["config_path"] = config  # May be None during help/version
    ctx.obj["log_level"] = log_level

    # Only configure logging if not in resilient parsing mode
    if not ctx.resilient_parsing:
        configure_logging(log_level)


# Register commands from commands module
from .commands.analyze import analyze  # noqa: E402
from .commands.reads import reads  # noqa: E402
from .commands.simulate import simulate  # noqa: E402

cli.add_command(simulate)
cli.add_command(reads)
cli.add_command(analyze)


# ============================================================================
# Entry Point
# ============================================================================


def main():
    """Entry point for Click CLI."""
    cli(obj={})


if __name__ == "__main__":
    main()
