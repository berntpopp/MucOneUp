"""Shared utilities for test data generation."""

import logging
import subprocess
from pathlib import Path


def ensure_dir(path: Path) -> Path:
    """
    Ensure directory exists, create if needed.

    Args:
        path: Directory path to ensure exists

    Returns:
        Path object (same as input)
    """
    path.mkdir(parents=True, exist_ok=True)
    return path


def run_command(
    cmd: list[str], description: str, check: bool = True
) -> subprocess.CompletedProcess:
    """
    Run a command with logging.

    Args:
        cmd: Command and arguments as list
        description: Human-readable description for logging
        check: Raise exception on non-zero exit (default: True)

    Returns:
        CompletedProcess instance

    Raises:
        subprocess.CalledProcessError: If check=True and command fails
    """
    logger = logging.getLogger(__name__)

    logger.info(f"{description}...")
    logger.debug(f"Command: {' '.join(str(c) for c in cmd)}")

    try:
        result = subprocess.run(cmd, check=check, capture_output=True, text=True)

        if result.stdout:
            logger.debug(f"stdout: {result.stdout[:500]}")
        if result.stderr:
            logger.debug(f"stderr: {result.stderr[:500]}")

        logger.info(f"✓ {description} completed")

        return result

    except subprocess.CalledProcessError as e:
        logger.error(f"✗ {description} failed")
        logger.error(f"Exit code: {e.returncode}")
        if e.stdout:
            logger.error(f"stdout: {e.stdout}")
        if e.stderr:
            logger.error(f"stderr: {e.stderr}")
        raise


def collect_all_files(root_dir: Path, exclude_patterns: list[str] | None = None) -> list[Path]:
    """
    Recursively collect all files in directory.

    Args:
        root_dir: Root directory to search
        exclude_patterns: Optional list of glob patterns to exclude

    Returns:
        Sorted list of Path objects for all files
    """
    exclude_patterns = exclude_patterns or []
    files = []

    for item in root_dir.rglob("*"):
        if item.is_file():
            # Check exclusions
            skip = False
            for pattern in exclude_patterns:
                if item.match(pattern):
                    skip = True
                    break

            if not skip:
                files.append(item)

    return sorted(files)
