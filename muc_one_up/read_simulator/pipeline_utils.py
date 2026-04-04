"""Shared utilities for read simulation pipelines.

Provides four helpers that are reused across Illumina, ONT, PacBio, and
Amplicon pipeline entry-points:

- resolve_pipeline_outputs: determine output directory and base name
- resolve_human_reference: locate human reference FASTA with a 3-tier fallback
- create_pipeline_metadata: thin Path-aware wrapper over write_metadata_file
- cleanup_intermediates: safely remove intermediate files without raising
"""

from __future__ import annotations

import logging
from collections.abc import Sequence
from datetime import datetime
from pathlib import Path
from typing import Any

from ..bioinformatics.reference_validation import get_reference_path_for_assembly
from ..exceptions import ConfigurationError
from .assembly_context import AssemblyContext
from .output_config import OutputConfig
from .utils.metadata_writer import write_metadata_file

logger = logging.getLogger(__name__)


def resolve_pipeline_outputs(
    input_fa: str,
    rs_config: dict[str, Any],
    output_config: OutputConfig | None,
) -> tuple[Path, str]:
    """Resolve output directory and base name for a pipeline run.

    Priority order:
    1. ``output_config`` (if not None)
    2. ``rs_config["output_dir"]`` (if present)
    3. parent directory of ``input_fa``

    The base name is always derived from the stem of ``input_fa`` unless
    ``output_config`` provides an explicit ``out_base``.

    The resolved output directory is created (including parents) if it does
    not already exist.

    Args:
        input_fa: Path to the input FASTA file.
        rs_config: Read-simulator config section (may be empty).
        output_config: Explicit output configuration, or ``None``.

    Returns:
        Tuple of ``(output_dir, output_base)``.
    """
    input_path = Path(input_fa)

    if output_config is not None:
        out_dir = output_config.out_dir
        out_base = output_config.out_base
    elif rs_config.get("output_dir"):
        out_dir = Path(rs_config["output_dir"])
        out_base = input_path.stem
    else:
        out_dir = input_path.parent
        out_base = input_path.stem

    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir, out_base


def resolve_human_reference(
    config: dict[str, Any],
    assembly_ctx: AssemblyContext,
    aligner: str = "bwa",
) -> str:
    """Resolve path to the human reference FASTA with a 3-tier fallback.

    Fallback order:
    1. ``config["reference_genomes"]`` section via
       :func:`get_reference_path_for_assembly`.
    2. ``assembly_ctx.human_reference`` (set from assembly-keyed rs_config key).
    3. ``config["read_simulation"]["human_reference"]`` (legacy flat key).

    Args:
        config: Main simulation configuration dictionary.
        assembly_ctx: Resolved assembly context for the current run.
        aligner: Aligner name passed to reference validation (default ``"bwa"``).

    Returns:
        Absolute path to the human reference FASTA as a string.

    Raises:
        ConfigurationError: If no reference can be found through any fallback.
    """
    # Tier 1: reference_genomes config section
    try:
        ref_path = get_reference_path_for_assembly(config, assembly_ctx.assembly_name)
        logger.debug(
            "Resolved human reference from reference_genomes: %s", ref_path
        )
        return str(ref_path)
    except Exception:
        logger.debug(
            "reference_genomes lookup failed for '%s', trying fallbacks",
            assembly_ctx.assembly_name,
        )

    # Tier 2: assembly_ctx.human_reference
    if assembly_ctx.human_reference:
        logger.debug(
            "Resolved human reference from assembly_ctx: %s",
            assembly_ctx.human_reference,
        )
        return assembly_ctx.human_reference

    # Tier 3: legacy config["read_simulation"]["human_reference"]
    legacy: str | None = config.get("read_simulation", {}).get("human_reference")
    if legacy:
        logger.debug("Resolved human reference from legacy read_simulation config: %s", legacy)
        return legacy

    raise ConfigurationError(
        f"human_reference not specified for assembly '{assembly_ctx.assembly_name}'. "
        "Provide it via config['reference_genomes'], assembly_ctx.human_reference, "
        "or config['read_simulation']['human_reference']."
    )


def create_pipeline_metadata(
    output_dir: Path,
    output_base: str,
    config: dict[str, Any],
    start_time: datetime,
    end_time: datetime,
    platform: str,
    tools_used: list[str],
) -> str:
    """Write pipeline metadata TSV, accepting a :class:`~pathlib.Path` for output_dir.

    This is a thin wrapper around
    :func:`~muc_one_up.read_simulator.utils.metadata_writer.write_metadata_file`
    that converts the ``Path`` argument to a string so callers do not need to
    handle the conversion themselves.

    Args:
        output_dir: Directory in which to write the metadata file.
        output_base: Base name for output files (without extension).
        config: Pipeline configuration dictionary.
        start_time: Pipeline start timestamp.
        end_time: Pipeline end timestamp.
        platform: Sequencing platform (``"Illumina"``, ``"ONT"``, or ``"PacBio"``).
        tools_used: List of tool keys actually used in this pipeline run.

    Returns:
        Path to the created metadata file as a string.
    """
    return write_metadata_file(
        str(output_dir),
        output_base,
        config,
        start_time,
        end_time,
        platform,
        tools_used,
    )


def cleanup_intermediates(file_list: Sequence[str | None]) -> None:
    """Remove intermediate files, silently skipping missing or null entries.

    Entries that are ``None`` or the empty string are skipped without logging.
    For entries that exist but cannot be removed, a warning is logged and
    execution continues — this function **never raises**.

    Args:
        file_list: Paths to remove. May contain ``None`` or empty strings.
    """
    for entry in file_list:
        if not entry:
            continue
        path = Path(entry)
        if not path.exists():
            continue
        try:
            path.unlink()
            logger.debug("Removed intermediate file: %s", path)
        except OSError as exc:
            logger.warning("Could not remove intermediate file %s: %s", path, exc)
