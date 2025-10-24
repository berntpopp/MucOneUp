#!/usr/bin/env python3
"""Metadata file generation for simulation provenance.

This module provides functionality to write TSV metadata files containing
simulation provenance information including tool versions, run parameters,
and timing information. Follows separation of concerns by creating a
separate metadata file rather than mixing with simulation statistics.

Design Principles:
    - TSV format: Human and machine readable
    - Separation of concerns: Separate file from simulation_stats.json
    - Namespace: Tool versions prefixed with "tool."
    - Based on VariantCentrifuge proven pattern

Output Example:
    Parameter                Value
    Tool                     MucOneUp
    Version                  0.15.0
    Platform                 Illumina
    tool.samtools_version    samtools 1.17
    tool.bwa_version         bwa 0.7.18-r1243
"""

import logging
from datetime import datetime
from pathlib import Path

from ...version import __version__
from .tool_version import capture_tool_versions, log_tool_versions

logger = logging.getLogger(__name__)


def write_metadata_file(
    output_dir: str,
    output_base: str,
    config: dict,
    start_time: datetime,
    end_time: datetime,
    platform: str,
) -> str:
    """
    Write TSV metadata file with provenance information.

    Args:
        output_dir: Output directory path
        output_base: Base name for output files (without extension)
        config: Pipeline configuration dictionary
        start_time: Pipeline start timestamp
        end_time: Pipeline end timestamp
        platform: Sequencing platform ("Illumina", "ONT", or "PacBio")

    Returns:
        Path to created metadata file (str)

    Example:
        >>> from datetime import datetime
        >>> config = {
        ...     "tools": {"samtools": "samtools", "bwa": "bwa"},
        ...     "seed": 42,
        ...     "read_simulation": {"coverage": 30, "fragment_size_mean": 350}
        ... }
        >>> start = datetime(2025, 10, 24, 12, 0, 0)
        >>> end = datetime(2025, 10, 24, 13, 10, 15)
        >>> path = write_metadata_file(
        ...     output_dir="/output",
        ...     output_base="sample_001",
        ...     config=config,
        ...     start_time=start,
        ...     end_time=end,
        ...     platform="Illumina"
        ... )
        # Creates: /output/sample_001_metadata.tsv

    Output Format (TSV):
        Parameter               Value
        Tool                    MucOneUp
        Version                 0.15.0
        Sequencing_platform     Illumina
        Run_start_time          2025-10-24T12:00:00
        Run_end_time            2025-10-24T13:10:15
        Run_duration_seconds    4215.0
        Seed                    42
        Coverage                30
        Fragment_size           350
        tool.samtools_version   samtools 1.17
        tool.bwa_version        bwa 0.7.18-r1243

    Design:
        - TSV format (human + machine readable)
        - Separate from simulation_stats.json (separation of concerns)
        - Tool versions namespaced with "tool." prefix
        - Follows VariantCentrifuge proven pattern
        - Platform-specific parameters conditionally included
    """
    metadata_path = Path(output_dir) / f"{output_base}_metadata.tsv"
    logger.info(f"Writing metadata: {metadata_path}")

    # Capture tool versions
    tool_versions = capture_tool_versions(config.get("tools", {}))
    log_tool_versions(tool_versions)

    # Calculate duration
    duration = (end_time - start_time).total_seconds()

    # Write TSV
    with open(metadata_path, "w", encoding="utf-8") as f:
        f.write("Parameter\tValue\n")

        # Core metadata
        f.write("Tool\tMucOneUp\n")
        f.write(f"Version\t{__version__}\n")
        f.write(f"Sequencing_platform\t{platform}\n")
        f.write(f"Run_start_time\t{start_time.isoformat()}\n")
        f.write(f"Run_end_time\t{end_time.isoformat()}\n")
        f.write(f"Run_duration_seconds\t{duration:.1f}\n")

        # Simulation parameters
        f.write(f"Output_base\t{output_base}\n")
        f.write(f"Output_dir\t{output_dir}\n")

        if "seed" in config:
            f.write(f"Seed\t{config['seed']}\n")

        # Platform-specific parameters
        if platform == "Illumina":
            read_cfg = config.get("read_simulation", {})
            f.write(f"Coverage\t{read_cfg.get('coverage', 'N/A')}\n")
            f.write(f"Fragment_size\t{read_cfg.get('fragment_size_mean', 'N/A')}\n")
        elif platform == "ONT":
            ont_cfg = config.get("nanosim_params", {})
            f.write(f"Coverage\t{ont_cfg.get('coverage', 'N/A')}\n")
            f.write(f"Min_read_length\t{ont_cfg.get('min_read_length', 'N/A')}\n")
            f.write(f"Max_read_length\t{ont_cfg.get('max_read_length', 'N/A')}\n")
        elif platform == "PacBio":
            pb_cfg = config.get("pacbio_params", {})
            f.write(f"Coverage\t{pb_cfg.get('coverage', 'N/A')}\n")
            f.write(f"Pass_num\t{pb_cfg.get('pass_num', 'N/A')}\n")

        # Tool versions (namespaced with "tool." prefix, sorted alphabetically)
        for tool_name, version in sorted(tool_versions.items()):
            f.write(f"tool.{tool_name}_version\t{version}\n")

    logger.info(f"Metadata written: {metadata_path}")
    return str(metadata_path)
