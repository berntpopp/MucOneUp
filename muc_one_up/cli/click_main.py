"""Click-based CLI for MucOneUp — command registration."""

import click

from ._common import (
    CONTEXT_SETTINGS,
    configure_logging,
    print_version_callback,
    validate_config_callback,
)


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
    is_eager=True,
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
    Commands:
      simulate  - Generate haplotypes
      reads     - Simulate reads from FASTA
      analyze   - Analyze FASTA (ORFs, stats, VNTR structure)
    """
    if verbose:
        log_level = "DEBUG"

    ctx.ensure_object(dict)
    ctx.obj["config_path"] = config
    ctx.obj["log_level"] = log_level

    if not ctx.resilient_parsing:
        configure_logging(log_level)


# Register commands
from .commands.analyze import analyze  # noqa: E402
from .commands.reads import reads  # noqa: E402
from .commands.simulate import simulate  # noqa: E402

cli.add_command(simulate)
cli.add_command(reads)
cli.add_command(analyze)


def main():
    """Entry point for Click CLI."""
    cli(obj={})


if __name__ == "__main__":
    main()
