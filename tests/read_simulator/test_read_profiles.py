"""Tests for read profiles: loading, validation, config overlay and built-ins."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from muc_one_up.read_simulator.read_profiles import (
    active_read_profile,
    apply_read_profile,
    list_builtin_profiles,
    load_read_profile,
)

MINIMAL = {
    "schema_version": 1,
    "name": "test_ont",
    "platform": "ont",
    "calibration": "generic",
    "description": "unit-test profile",
    "provenance": {"source": "unit test"},
    "config_overrides": {
        "ont_amplicon_params": {"accuracy_mean": 0.985, "difference_ratio": "33:31:36"},
        "amplicon_params": {"pcr_bias": {"preset": "madritsch2025_r10"}},
    },
    "molecules": {"forward_frac": 0.5, "stutter": {"C7|+": {"-1": 0.3, "0": 0.7}}},
    "fragments": {"length_median": 5000, "length_sigma": 0.45},
}


def _write(tmp_path: Path, data: dict) -> Path:
    path = tmp_path / "profile.json"
    path.write_text(json.dumps(data))
    return path


def test_load_from_path_parses_all_parts(tmp_path: Path) -> None:
    profile = load_read_profile(str(_write(tmp_path, MINIMAL)))
    assert profile.name == "test_ont" and profile.platform == "ont"
    assert profile.molecules.forward_frac == 0.5
    assert profile.fragments.length_median == 5000
    assert len(profile.sha256) == 64


@pytest.mark.parametrize(
    "mutation,match",
    [
        ({"schema_version": 2}, "schema_version"),
        ({"platform": "illumina"}, "platform"),
        ({"surprise": 1}, "unknown"),
        ({"config_overrides": {"tools": {"pbsim3": "x"}}}, "config_overrides"),
        ({"molecules": {"forward_frac": 2}}, "forward_frac"),
    ],
)
def test_invalid_profiles_rejected(tmp_path: Path, mutation: dict, match: str) -> None:
    with pytest.raises(ValueError, match=match):
        load_read_profile(str(_write(tmp_path, {**MINIMAL, **mutation})))


def test_unknown_name_lists_builtins() -> None:
    with pytest.raises(ValueError, match="Available"):
        load_read_profile("no_such_profile")


def test_apply_overlays_config_and_marks_active(tmp_path: Path) -> None:
    path = _write(tmp_path, MINIMAL)
    config = {
        "ont_amplicon_params": {"model_file": "m.model", "accuracy_mean": 0.95},
        "amplicon_params": {"forward_primer": "ACGT", "reverse_primer": "TTGG"},
    }
    merged = apply_read_profile(config, load_read_profile(str(path)))
    assert merged["ont_amplicon_params"] == {
        "model_file": "m.model",
        "accuracy_mean": 0.985,
        "difference_ratio": "33:31:36",
    }
    assert merged["amplicon_params"]["pcr_bias"] == {"preset": "madritsch2025_r10"}
    assert merged["amplicon_params"]["forward_primer"] == "ACGT"
    assert config["ont_amplicon_params"]["accuracy_mean"] == 0.95  # input not mutated
    assert active_read_profile(merged).name == "test_ont"
    assert active_read_profile(config) is None


def test_apply_rejects_overlay_that_breaks_section_schema(tmp_path: Path) -> None:
    bad = {**MINIMAL, "config_overrides": {"ont_amplicon_params": {"accuracy_mean": 3}}}
    with pytest.raises(ValueError, match="ont_amplicon_params"):
        apply_read_profile({}, load_read_profile(str(_write(tmp_path, bad))))


@pytest.mark.parametrize("name", list_builtin_profiles())
def test_builtin_profiles_load(name: str) -> None:
    profile = load_read_profile(name)
    assert profile.name == name
    assert profile.provenance.get("source")


def test_config_schema_accepts_read_model_section() -> None:
    from jsonschema import validate

    from muc_one_up.config import CONFIG_SCHEMA

    validate(instance={"profile": "x"}, schema=CONFIG_SCHEMA["properties"]["read_model"])


def test_ont_amplicon_params_reject_accuracy_sd():
    """pbsim3 has no --accuracy-sd, so the key must not validate (#115)."""
    import jsonschema

    from muc_one_up.config import CONFIG_SCHEMA

    schema = CONFIG_SCHEMA["properties"]["ont_amplicon_params"]
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate({"accuracy_sd": 0.01}, {**schema, "required": []})


def test_profile_fragments_beat_config_file(tmp_path: Path) -> None:
    """config < profile: config-file fragment lengths yield to the profile (#114)."""
    from muc_one_up.read_simulator.ont_fragment_pipeline import _fragment_settings

    config = {"ont_fragment_params": {"length_median": 1000, "length_sigma": 0.9, "n_reads": 7}}
    merged = apply_read_profile(config, load_read_profile(str(_write(tmp_path, MINIMAL))))
    assert merged["ont_fragment_params"] == {
        "length_median": 5000,
        "length_sigma": 0.45,
        "n_reads": 7,
    }
    _, fragments, _ = _fragment_settings(merged)
    assert (fragments.length_median, fragments.length_sigma) == (5000, 0.45)


def test_cli_fragment_flags_still_beat_profile(tmp_path: Path) -> None:
    from muc_one_up.read_simulator.ont_fragment_pipeline import _fragment_settings

    merged = apply_read_profile({}, load_read_profile(str(_write(tmp_path, MINIMAL))))
    merged["ont_fragment_params"]["length_median"] = 800  # what reads ont --read-length-median does
    _, fragments, _ = _fragment_settings(merged)
    assert (fragments.length_median, fragments.length_sigma) == (800, 0.45)


def test_profile_without_fragments_keeps_config_lengths(tmp_path: Path) -> None:
    data = {k: v for k, v in MINIMAL.items() if k != "fragments"}
    config = {"ont_fragment_params": {"length_median": 1000}}
    merged = apply_read_profile(config, load_read_profile(str(_write(tmp_path, data))))
    assert merged["ont_fragment_params"] == {"length_median": 1000}


def test_profile_pcr_preset_drops_config_preset_parameters(tmp_path: Path) -> None:
    """A config alpha/e_max must not silently override the profile's preset (#114)."""
    config = {
        "amplicon_params": {
            "forward_primer": "ACGT",
            "reverse_primer": "TTGG",
            "pcr_bias": {"preset": "default", "alpha": 0.5, "e_max": 0.7, "stochastic": True},
        }
    }
    merged = apply_read_profile(config, load_read_profile(str(_write(tmp_path, MINIMAL))))
    assert merged["amplicon_params"]["pcr_bias"] == {
        "preset": "madritsch2025_r10",
        "stochastic": True,
    }


@pytest.mark.parametrize(
    "mutation,match",
    [
        ({"name": None}, "name"),
        ({"fragments": {"median": 5}}, "fragments"),
        ({"errors": {"mismatch_rate": 0.01}}, "errors"),
        ({"molecules": {"forward_frac": "0.5"}}, "molecules"),
        ({"molecules": {"stutter": {"C7|+": {"-8": 1.0}}}}, "C7"),
    ],
)
def test_malformed_profiles_raise_value_error(tmp_path: Path, mutation: dict, match: str) -> None:
    """Construction errors must surface as ValueError naming the problem (#124)."""
    data = {**MINIMAL, **mutation}
    data = {k: v for k, v in data.items() if v is not None}
    with pytest.raises(ValueError, match=match):
        load_read_profile(str(_write(tmp_path, data)))
