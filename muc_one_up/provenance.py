"""Provenance metadata collection for reproducibility.

Collects comprehensive provenance information for each simulation,
including software version, configuration fingerprint, random seed,
timestamps, and command-line invocation. Follows scientific software
best practices for reproducibility tracking.

References:
    - "The role of metadata in reproducible computational research"
      (Patterns, 2021): https://doi.org/10.1016/j.patter.2021.100322
    - Research Software Engineering with Python, Chapter 13:
      https://third-bit.com/py-rse/provenance.html

Example:
    >>> from muc_one_up.provenance import collect_provenance_metadata
    >>> import time
    >>>
    >>> config = {"repeats": {"X": "ACGACT"}, "seed": 42}
    >>> start = time.time()
    >>> # ... run simulation ...
    >>> end = time.time()
    >>>
    >>> provenance = collect_provenance_metadata(config, start, end)
    >>> print(provenance["software_version"])
    0.27.0
    >>> print(provenance["seed"])
    42
    >>> print(provenance["config_fingerprint"][:20])
    sha256:8f3e9a7b2c1d...

Security:
    Command-line sanitization redacts common secret patterns to prevent
    credential leakage in published datasets. However, this is NOT foolproof.
    Best practice: NEVER pass secrets via command-line arguments!
"""

from __future__ import annotations

import logging
import os
import sys
from datetime import datetime, timezone
from typing import Any, cast

from .config_fingerprint import compute_config_fingerprint
from .version import __version__

logger = logging.getLogger(__name__)

# Feature flag for emergency rollback
ENABLE_PROVENANCE = os.getenv("MUCONEUP_ENABLE_PROVENANCE", "true").lower() in (
    "true",
    "1",
    "yes",
    "on",
)


def get_software_version() -> str:
    """
    Get current software version.

    Returns:
        Version string (e.g., "0.27.0")

    Example:
        >>> from muc_one_up.provenance import get_software_version
        >>> version = get_software_version()
        >>> isinstance(version, str)
        True
        >>> version.count('.') >= 2  # Semantic versioning
        True
    """
    return __version__


def extract_seed(config: dict[str, Any]) -> int | None:
    """
    Extract random seed from configuration.

    Checks multiple possible locations in order of priority:
    1. config["seed"] (simulate command)
    2. config["read_simulation"]["seed"] (illumina reads)
    3. config["nanosim_params"]["seed"] (ONT reads)
    4. config["pacbio_params"]["seed"] (PacBio reads)

    Args:
        config: Configuration dictionary

    Returns:
        Seed value if found, None otherwise

    Example:
        >>> config = {"seed": 42}
        >>> extract_seed(config)
        42
        >>> config = {"read_simulation": {"seed": 123}}
        >>> extract_seed(config)
        123
        >>> config = {}
        >>> extract_seed(config) is None
        True
    """
    try:
        # Priority 1: Top-level seed (simulate command)
        if "seed" in config and config["seed"] is not None:
            return cast(int, config["seed"])

        # Priority 2: Read simulation seeds
        if "read_simulation" in config:
            rs_seed = config["read_simulation"].get("seed")
            if rs_seed is not None:
                return cast(int, rs_seed)

        # Priority 3: NanoSim seed
        if "nanosim_params" in config:
            ns_seed = config["nanosim_params"].get("seed")
            if ns_seed is not None:
                return cast(int, ns_seed)

        # Priority 4: PacBio seed
        if "pacbio_params" in config:
            pb_seed = config["pacbio_params"].get("seed")
            if pb_seed is not None:
                return cast(int, pb_seed)

        logger.debug("No seed found in configuration")
        return None

    except (TypeError, KeyError, AttributeError) as e:
        logger.error(f"Error extracting seed from config: {e}")
        return None


def format_timestamp_iso8601(timestamp: float) -> str:
    """
    Format Unix timestamp (float) to ISO 8601 with timezone.

    Uses RFC 3339 profile of ISO 8601 for maximum compatibility.
    Always includes microseconds and UTC timezone offset.

    Args:
        timestamp: Unix timestamp (seconds since epoch)

    Returns:
        ISO 8601 formatted string with timezone

    Raises:
        ValueError: If timestamp is invalid (negative, too large)
        OSError: If timestamp is out of range for platform

    Example:
        >>> import time
        >>> ts = 1730649330.123456  # 2024-11-03 10:15:30.123456 UTC
        >>> iso = format_timestamp_iso8601(ts)
        >>> iso
        '2024-11-03T10:15:30.123456+00:00'

    Notes:
        - Assumes timestamp is in UTC (Unix convention)
        - Always includes 6-digit microseconds
        - Always includes timezone offset (+00:00 for UTC)
    """
    try:
        # Convert Unix timestamp to timezone-aware datetime (UTC)
        dt = datetime.fromtimestamp(timestamp, tz=timezone.utc)

        # Format to ISO 8601 (RFC 3339 profile)
        return dt.isoformat()

    except (ValueError, OSError) as e:
        # Invalid timestamp (negative, too large, etc.)
        logger.error(f"Failed to format timestamp {timestamp}: {e}")
        raise


def get_command_line(sanitize: bool = True) -> str:
    """
    Get command-line invocation string.

    Args:
        sanitize: If True, redact common secret patterns (default: True)

    Returns:
        Command line string with optional secret redaction

    Security:
        Even with sanitization, this is NOT foolproof. Best practice:
        NEVER pass secrets via command-line arguments. Use environment
        variables or configuration files instead.

    Example:
        >>> import sys
        >>> # Mock sys.argv for testing
        >>> original_argv = sys.argv
        >>> sys.argv = ['muconeup', '--config', 'config.json', 'simulate']
        >>> cmd = get_command_line(sanitize=True)
        >>> 'muconeup' in cmd
        True
        >>> sys.argv = original_argv

        >>> # With secrets:
        >>> sys.argv = ['muconeup', '--api-key', 'secret123', 'run']
        >>> cmd = get_command_line(sanitize=True)
        >>> 'secret123' not in cmd
        True
        >>> '***REDACTED***' in cmd
        True
        >>> sys.argv = original_argv
    """
    try:
        cmd_parts = sys.argv[:]

        if sanitize:
            # Patterns that likely contain secrets
            secret_flags = {
                "--api-key",
                "--password",
                "--token",
                "--secret",
                "--key",
                "--auth",
                "--credential",
                "--passphrase",
            }

            sanitized = []
            skip_next = False

            for part in cmd_parts:
                if skip_next:
                    sanitized.append("***REDACTED***")
                    skip_next = False
                    continue

                # Check if this is a secret flag (--api-key VALUE)
                if part in secret_flags:
                    sanitized.append(part)
                    skip_next = True  # Redact next argument
                # Check for inline secrets (--api-key=VALUE)
                elif any(part.startswith(f"{flag}=") for flag in secret_flags):
                    flag_name = part.split("=")[0]
                    sanitized.append(f"{flag_name}=***REDACTED***")
                else:
                    sanitized.append(part)

            return " ".join(sanitized)

        return " ".join(cmd_parts)

    except Exception as e:
        logger.error(f"Failed to get command line: {e}")
        return "error:cmdline_failed"


def collect_provenance_metadata(
    config: dict[str, Any],
    start_time: float,
    end_time: float,
) -> dict[str, Any]:
    """
    Collect comprehensive provenance metadata.

    Gathers all information necessary to reproduce a simulation:
    - Software version
    - Configuration fingerprint (SHA-256)
    - Random seed (if used)
    - Timestamps (start, end, duration)
    - Command-line invocation (sanitized)

    Args:
        config: Configuration dictionary
        start_time: Simulation start time (Unix timestamp, float)
        end_time: Simulation end time (Unix timestamp, float)

    Returns:
        Provenance metadata dictionary with keys:
        - software_version: str
        - config_fingerprint: str ("sha256:..." or "error:...")
        - seed: int | None
        - start_time: str (ISO 8601)
        - end_time: str (ISO 8601)
        - duration_seconds: float
        - command_line: str (sanitized)

    Example:
        >>> import time
        >>> config = {"repeats": {"X": "ACGACT"}, "seed": 42}
        >>> start = time.time()
        >>> time.sleep(0.1)
        >>> end = time.time()
        >>> provenance = collect_provenance_metadata(config, start, end)
        >>> provenance["software_version"]
        '0.27.0'
        >>> provenance["seed"]
        42
        >>> provenance["duration_seconds"] > 0.09
        True

    Notes:
        - All timestamps converted to ISO 8601 with timezone
        - Missing seed returns None (graceful degradation)
        - Config fingerprint errors return "error:..." sentinel
        - Never raises exceptions (returns error sentinels)
        - Respects MUCONEUP_ENABLE_PROVENANCE environment variable
    """
    # Check feature flag
    if not ENABLE_PROVENANCE:
        logger.info("Provenance collection disabled via MUCONEUP_ENABLE_PROVENANCE")
        return {}

    logger.debug("Collecting provenance metadata for reproducibility tracking")

    provenance: dict[str, Any] = {}

    # Version (should never fail)
    try:
        provenance["software_version"] = get_software_version()
    except Exception as e:
        logger.error(f"Failed to get software version: {e}")
        provenance["software_version"] = "unknown"

    # Config fingerprint (may fail with complex configs)
    try:
        provenance["config_fingerprint"] = compute_config_fingerprint(config)
    except Exception as e:
        logger.error(f"Failed to compute config fingerprint: {e}")
        provenance["config_fingerprint"] = "error:fingerprint_failed"

    # Seed (should never fail)
    try:
        provenance["seed"] = extract_seed(config)
    except Exception as e:
        logger.error(f"Failed to extract seed: {e}")
        provenance["seed"] = None

    # Timestamps (may fail with invalid timestamps)
    try:
        start_iso = format_timestamp_iso8601(start_time)
        end_iso = format_timestamp_iso8601(end_time)
        provenance["start_time"] = start_iso
        provenance["end_time"] = end_iso
        provenance["duration_seconds"] = end_time - start_time
    except (ValueError, OSError) as e:
        logger.error(f"Failed to format timestamps: {e}")
        provenance["start_time"] = "error:timestamp_failed"
        provenance["end_time"] = "error:timestamp_failed"
        # Still compute duration (even if formatting failed)
        provenance["duration_seconds"] = end_time - start_time

    # Command line (should never fail)
    try:
        provenance["command_line"] = get_command_line(sanitize=True)
    except Exception as e:
        logger.error(f"Failed to get command line: {e}")
        provenance["command_line"] = "error:cmdline_failed"

    # Log summary (DEBUG level to avoid log spam in batch mode)
    if logger.isEnabledFor(logging.DEBUG):
        logger.debug(
            f"Provenance collected: version={provenance.get('software_version')}, "
            f"seed={provenance.get('seed')}, "
            f"duration={provenance.get('duration_seconds', 0):.3f}s"
        )

    return provenance
