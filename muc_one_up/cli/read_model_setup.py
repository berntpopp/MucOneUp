"""CLI helpers that apply read profiles, tracking and alignment choices to a config."""

from __future__ import annotations

import logging
from typing import Any

import click

from ..read_simulator.read_profiles import apply_read_profile, load_read_profile


def apply_cli_read_model(
    config: dict[str, Any], read_profile: str | None, platform: str
) -> dict[str, Any]:
    """Overlay ``--read-profile`` (or a profile named in the config) onto ``config``.

    Precedence: config file < read profile < explicit CLI flags applied afterwards.
    The profile platform must match the command platform ("ont" or "pacbio").
    """
    ref = read_profile or (config.get("read_model") or {}).get("profile")
    if not ref:
        return config
    try:
        profile = load_read_profile(str(ref))
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc
    if profile.platform != platform:
        raise click.ClickException(
            f"Read profile '{profile.name}' is for platform '{profile.platform}', not '{platform}'."
        )
    logging.info("Using read profile %s (sha256 %s)", profile.name, profile.sha256[:12])
    try:
        return apply_read_profile(config, profile)
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc


def apply_tracking_and_alignment(
    config: dict[str, Any], track_read_source: bool, no_align: bool
) -> None:
    """Record truth tracking and drop the alignment reference when requested."""
    rs = config.setdefault("read_simulation", {})
    if track_read_source:
        rs["track_read_source"] = True
    if no_align and rs.pop("human_reference", None):
        logging.info("--no-align: alignment skipped; output is FASTQ")
