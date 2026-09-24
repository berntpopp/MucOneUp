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
import logging
from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, TypeVar

from jsonschema import ValidationError as SchemaError
from jsonschema import validate

from ..config import CONFIG_SCHEMA
from .empirical_errors import EmpiricalErrorModel
from .molecules import MoleculeModel

BUILTIN_DIR = Path(__file__).resolve().parent.parent / "data" / "read_profiles"
OVERRIDABLE_SECTIONS = ("ont_amplicon_params", "pacbio_params", "amplicon_params")
PLATFORMS = ("ont", "pacbio")
_T = TypeVar("_T")
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
    fragments: FragmentModel | None  # None when the profile leaves read lengths to the config
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
    if not isinstance(data.get("name"), str) or not data["name"]:
        raise ValueError(f"read profile {path}: 'name' must be a non-empty string")
    overrides = data.get("config_overrides", {})
    _validate_overrides(overrides)

    def part(name: str, build: Callable[[], _T]) -> _T:
        try:
            return build()
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"read profile {path}: invalid '{name}': {exc}") from exc

    fragments = data.get("fragments")
    errors = data.get("errors")
    molecules = part("molecules", lambda: MoleculeModel.from_dict(data.get("molecules", {})))
    error_model = part("errors", lambda: EmpiricalErrorModel.from_dict(errors)) if errors else None
    if error_model is not None:
        _check_protected_runs_are_stuttered(path, molecules, error_model)
    return ReadProfile(
        name=data["name"],
        platform=data["platform"],
        calibration=str(data.get("calibration", "generic")),
        description=str(data.get("description", "")),
        provenance=data.get("provenance", {}),
        config_overrides=overrides,
        molecules=molecules,
        fragments=part("fragments", lambda: FragmentModel(**fragments)) if fragments else None,
        sha256=hashlib.sha256(raw_bytes).hexdigest(),
        source=str(path),
        errors=error_model,
    )


def _check_protected_runs_are_stuttered(
    path: Path, molecules: MoleculeModel, errors: EmpiricalErrorModel
) -> None:
    """Reject profiles whose protected homopolymer runs could be error-free (#132).

    The empirical channel applies no base-level errors inside runs of at least
    ``errors.hp_min_len`` bases, so the stutter table is their only error
    source: it must exist, reach down to ``hp_min_len``, define a fallback for
    runs without a fitted entry, and every fitted entry for a protected run
    must have error mass (P(delta = 0) < 1).
    """
    table = molecules.stutter
    if table is None or table.fallback is None:
        raise ValueError(
            f"read profile {path}: an 'errors' model requires 'molecules.stutter' and "
            "'molecules.stutter_fallback', because homopolymer runs >= hp_min_len get no "
            "base-level errors and would otherwise be simulated error-free"
        )
    for key, pmf in table.pmfs.items():
        if int(key.partition("|")[0][1:]) >= errors.hp_min_len and dict(pmf).get(0, 0.0) >= 1.0:
            raise ValueError(
                f"read profile {path}: stutter entry '{key}' has P(delta = 0) = 1, so this "
                "protected homopolymer run would be simulated error-free; fit it from data "
                "or remove it so the stutter_fallback applies"
            )
    if errors.hp_min_len < table.min_len:
        raise ValueError(
            f"read profile {path}: errors.hp_min_len ({errors.hp_min_len}) must be >= the "
            f"stutter table's minimum run length ({table.min_len})"
        )


def _deep_merge(base: dict[str, Any], overlay: Mapping[str, Any]) -> dict[str, Any]:
    for key, value in overlay.items():
        if isinstance(value, Mapping) and isinstance(base.get(key), dict):
            _deep_merge(base[key], value)
        else:
            base[key] = copy.deepcopy(value)
    return base


def _drop_config_pcr_parameters(
    section: dict[str, Any] | None, overlay: Mapping[str, Any], profile_name: str
) -> None:
    """Keep config-file PCR parameters from silently overriding a profile's preset.

    ``PCRBiasModel.from_config`` applies every non-preset key as an override of
    the preset, so config values such as ``alpha`` would otherwise beat the
    profile. Only the ``stochastic`` mode switch is kept from the config.
    """
    pcr_overlay = overlay.get("pcr_bias") or {}
    pcr = (section or {}).get("pcr_bias")
    if "preset" not in pcr_overlay or not isinstance(pcr, dict):
        return
    dropped = sorted(k for k in pcr if k not in ("preset", "stochastic") and k not in pcr_overlay)
    for key in dropped:
        del pcr[key]
    if dropped:
        logging.info(
            "Read profile %s sets PCR preset %s; ignoring config pcr_bias %s",
            profile_name,
            pcr_overlay["preset"],
            ", ".join(dropped),
        )


def apply_read_profile(config: Mapping[str, Any], profile: ReadProfile) -> dict[str, Any]:
    """Return a copy of ``config`` with the profile overlaid and marked active."""
    merged = copy.deepcopy(dict(config))
    for section, values in profile.config_overrides.items():
        if section == "amplicon_params":
            _drop_config_pcr_parameters(merged.get(section), values, profile.name)
        section_config = _deep_merge(merged.setdefault(section, {}), values)
        schema = copy.deepcopy(CONFIG_SCHEMA["properties"][section])
        schema.pop("required", None)  # reads commands accept partial sections
        try:
            validate(instance=section_config, schema=schema)
        except SchemaError as exc:
            raise ValueError(
                f"read profile '{profile.name}' makes {section} invalid: {exc.message}"
            ) from exc
    if profile.fragments is not None:
        fragment_params = merged.setdefault("ont_fragment_params", {})
        fragment_params["length_median"] = profile.fragments.length_median
        fragment_params["length_sigma"] = profile.fragments.length_sigma
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
