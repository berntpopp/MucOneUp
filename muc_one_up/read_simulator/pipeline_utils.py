"""Shared utilities for read simulation pipelines.

Provides four helpers that are reused across Illumina, ONT, PacBio, and
Amplicon pipeline entry-points:

- resolve_pipeline_outputs: determine output directory and base name
- resolve_human_reference: locate the alignment reference (one precedence for all)
- create_pipeline_metadata: thin Path-aware wrapper over write_metadata_file
- cleanup_intermediates: safely remove intermediate files without raising
"""

from __future__ import annotations

import logging
import tempfile
from collections.abc import Iterator, Mapping, Sequence
from contextlib import contextmanager
from datetime import datetime
from pathlib import Path
from typing import Any

from ..bioinformatics.reference_validation import (
    get_reference_path_for_assembly,
    validate_reference_for_assembly,
)
from ..exceptions import ConfigurationError, FileOperationError, ValidationError
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
    """Resolve the human reference FASTA used for alignment by every pipeline.

    Explicit settings win over the shared assembly table:

    1. ``read_simulation.human_reference_<assembly>`` (``assembly_ctx.human_reference``).
    2. ``read_simulation.human_reference``.
    3. ``reference_genomes[<assembly>].fasta_path`` (validated for ``aligner`` indices).

    Raises:
        ConfigurationError: If no reference is configured or the
            ``reference_genomes`` entry cannot be used (the reason is included).
    """
    explicit = assembly_ctx.human_reference or config.get("read_simulation", {}).get(
        "human_reference"
    )
    if explicit:
        logger.info("Using human reference from read_simulation config: %s", explicit)
        return str(explicit)

    try:
        ref_path = get_reference_path_for_assembly(config, assembly_ctx.assembly_name)
    except (ValidationError, FileOperationError) as exc:
        raise ConfigurationError(
            f"No human reference for assembly '{assembly_ctx.assembly_name}': set "
            f"read_simulation.human_reference or reference_genomes.{assembly_ctx.assembly_name}"
            f".fasta_path ({exc})"
        ) from exc

    logger.info("Using reference_genomes reference: %s (%s)", ref_path, assembly_ctx.assembly_name)
    validate_reference_for_assembly(config, assembly_ctx.assembly_name, aligner=aligner)  # logs
    return str(ref_path)


def create_pipeline_metadata(
    output_dir: Path,
    output_base: str,
    config: dict[str, Any],
    start_time: datetime,
    end_time: datetime,
    platform: str,
    tools_used: list[str],
    *,
    extra_rows: list[tuple[str, Any]] | None = None,
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
        extra_rows: Additional metadata rows (e.g. the read truth manifest).

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
        extra_rows=extra_rows,
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


@contextmanager
def intermediate_work_dir(
    config: Mapping[str, Any], output_dir: Path, prefix: str
) -> Iterator[Path]:
    """Scratch directory for a pipeline run.

    Kept under ``output_dir`` when ``read_simulation.keep_intermediate_files`` is
    true, otherwise a temporary directory removed on exit.
    """
    if config.get("read_simulation", {}).get("keep_intermediate_files", False):
        work_dir = Path(tempfile.mkdtemp(prefix=prefix, dir=output_dir))
        logger.info("Keeping intermediate work directory: %s", work_dir)
        yield work_dir
    else:
        with tempfile.TemporaryDirectory(prefix=prefix) as tmp:
            yield Path(tmp)


def cleanup_unless_kept(config: Mapping[str, Any], file_list: Sequence[str | None]) -> None:
    """Remove intermediates unless ``read_simulation.keep_intermediate_files`` is true."""
    if config.get("read_simulation", {}).get("keep_intermediate_files", False):
        logger.info("Keeping intermediate files (keep_intermediate_files=true)")
        return
    cleanup_intermediates(file_list)
