"""Versioned read profiles: pbsim3/PCR config overlays plus a molecule model.

A read profile bundles everything that makes simulated long reads resemble a
specific real library type:

* ``config_overrides``: values deep-merged into existing config sections
  (``ont_amplicon_params``, ``pacbio_params``, ``amplicon_params``), validated
  against the same JSON schema as the config file;
* ``molecules``: a :class:`MoleculeModel` (strand mix, PCR artefacts, stutter);
* ``fragments``: read-length model for genomic fragment simulation;
* ``errors`` (optional): a calibrated :class:`EmpiricalErrorModel`; when present
  reads are sequenced by the empirical channel instead of pbsim3.

Profiles contain aggregate statistics only. Built-ins live in
``muc_one_up/data/read_profiles/<name>.json``; any JSON path also works.
Activating a profile records ``config["read_model"] = {"profile": ref}`` so the
config stays JSON-serialisable and pipelines can resolve the profile again.
"""

from __future__ import annotations

import copy
import hashlib
import json
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from jsonschema import ValidationError as SchemaError
from jsonschema import validate

from ..config import CONFIG_SCHEMA
from .empirical_errors import EmpiricalErrorModel
from .molecules import MoleculeModel

BUILTIN_DIR = Path(__file__).resolve().parent.parent / "data" / "read_profiles"
OVERRIDABLE_SECTIONS = ("ont_amplicon_params", "pacbio_params", "amplicon_params")
PLATFORMS = ("ont", "pacbio")
_KEYS = {
    "schema_version",
    "name",
    "platform",
    "calibration",
    "description",
    "provenance",
    "config_overrides",
    "molecules",
    "fragments",
    "errors",
}


@dataclass(frozen=True)
class FragmentModel:
    """Log-normal read-length model for genomic fragment simulation."""

    length_median: float = 5000.0
    length_sigma: float = 0.5

    def __post_init__(self) -> None:
        if self.length_median <= 0 or self.length_sigma <= 0:
            raise ValueError("fragments.length_median and fragments.length_sigma must be positive")


@dataclass(frozen=True)
class ReadProfile:
    name: str
    platform: str
    calibration: str
    description: str
    provenance: Mapping[str, Any]
    config_overrides: Mapping[str, Any]
    molecules: MoleculeModel
    fragments: FragmentModel
    sha256: str
    source: str = field(default="")
    errors: EmpiricalErrorModel | None = None


def list_builtin_profiles() -> list[str]:
    """Names of the profiles shipped with MucOneUp."""
    return sorted(p.stem for p in BUILTIN_DIR.glob("*.json"))


def _resolve(ref: str) -> Path:
    path = Path(ref)
    if path.suffix == ".json" or path.exists():
        if not path.exists():
            raise ValueError(f"Read profile file not found: {ref}")
        return path
    builtin = BUILTIN_DIR / f"{ref}.json"
    if not builtin.exists():
        raise ValueError(
            f"Unknown read profile '{ref}'. Available: {', '.join(list_builtin_profiles())}"
        )
    return builtin


def _validate_overrides(overrides: Mapping[str, Any]) -> None:
    extra = set(overrides) - set(OVERRIDABLE_SECTIONS)
    if extra:
        raise ValueError(
            f"config_overrides may only contain {OVERRIDABLE_SECTIONS}, got {sorted(extra)}"
        )


def load_read_profile(ref: str) -> ReadProfile:
    """Load and validate a built-in profile name or a profile JSON path."""
    path = _resolve(ref)
    raw_bytes = path.read_bytes()
    data = json.loads(raw_bytes)
    unknown = set(data) - _KEYS
    if unknown:
        raise ValueError(f"unknown read profile keys: {sorted(unknown)}")
    if data.get("schema_version") != 1:
        raise ValueError("read profile schema_version must be 1")
    if data.get("platform") not in PLATFORMS:
        raise ValueError(f"read profile platform must be one of {PLATFORMS}")
    overrides = data.get("config_overrides", {})
    _validate_overrides(overrides)
    return ReadProfile(
        name=str(data["name"]),
        platform=data["platform"],
        calibration=str(data.get("calibration", "generic")),
        description=str(data.get("description", "")),
        provenance=data.get("provenance", {}),
        config_overrides=overrides,
        molecules=MoleculeModel.from_dict(data.get("molecules", {})),
        fragments=FragmentModel(**data.get("fragments", {})),
        sha256=hashlib.sha256(raw_bytes).hexdigest(),
        source=str(path),
        errors=EmpiricalErrorModel.from_dict(data["errors"]) if data.get("errors") else None,
    )


def _deep_merge(base: dict[str, Any], overlay: Mapping[str, Any]) -> dict[str, Any]:
    for key, value in overlay.items():
        if isinstance(value, Mapping) and isinstance(base.get(key), dict):
            _deep_merge(base[key], value)
        else:
            base[key] = copy.deepcopy(value)
    return base


def apply_read_profile(config: Mapping[str, Any], profile: ReadProfile) -> dict[str, Any]:
    """Return a copy of ``config`` with the profile overlaid and marked active."""
    merged = copy.deepcopy(dict(config))
    for section, values in profile.config_overrides.items():
        section_config = _deep_merge(merged.setdefault(section, {}), values)
        schema = copy.deepcopy(CONFIG_SCHEMA["properties"][section])
        schema.pop("required", None)  # reads commands accept partial sections
        try:
            validate(instance=section_config, schema=schema)
        except SchemaError as exc:
            raise ValueError(
                f"read profile '{profile.name}' makes {section} invalid: {exc.message}"
            ) from exc
    merged["read_model"] = {
        "profile": profile.source or profile.name,
        "name": profile.name,
        "sha256": profile.sha256,
    }
    return merged


def active_read_profile(config: Mapping[str, Any]) -> ReadProfile | None:
    """The profile recorded in ``config['read_model']``, or None for legacy simulation."""
    read_model = config.get("read_model")
    if not read_model:
        return None
    return load_read_profile(str(read_model["profile"]))
