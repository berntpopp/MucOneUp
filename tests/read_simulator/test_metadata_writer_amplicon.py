"""Tests for amplicon-specific metadata TSV fields (issue #110)."""

from datetime import datetime
from pathlib import Path
from unittest.mock import patch

import pytest

from muc_one_up.read_simulator.utils.metadata_writer import write_metadata_file


def _write(tmp_path: Path, config: dict, platform: str) -> dict[str, str]:
    """Write metadata with mocked tool versions and return it as a dict."""
    with (
        patch(
            "muc_one_up.read_simulator.utils.metadata_writer.capture_tool_versions",
            return_value={},
        ),
        patch("muc_one_up.read_simulator.utils.metadata_writer.log_tool_versions"),
    ):
        path = write_metadata_file(
            output_dir=str(tmp_path),
            output_base="amp",
            config=config,
            start_time=datetime(2025, 1, 1, 12, 0, 0),
            end_time=datetime(2025, 1, 1, 12, 1, 0),
            platform=platform,
        )
    rows = Path(path).read_text().splitlines()[1:]
    return dict(row.split("\t", 1) for row in rows)


def test_ont_amplicon_reports_amplicon_fields_not_nanosim(tmp_path):
    config = {
        "read_simulation": {"coverage": 500, "assay_type": "amplicon"},
        "nanosim_params": {"coverage": 200, "min_read_length": 2000, "max_read_length": 7500},
        "ont_amplicon_params": {"model_type": "qshmm", "model_file": "QSHMM-ONT-HQ.model"},
        "amplicon_params": {"pcr_bias": {"preset": "madritsch2025_r10"}},
    }
    meta = _write(tmp_path, config, "ONT")

    assert meta["Read_simulation_technology"] == "ONT"
    assert meta["Assay_type"] == "amplicon"
    assert meta["Template_molecules"] == "500"
    assert meta["Model_type"] == "qshmm"
    assert meta["Model_file"] == "QSHMM-ONT-HQ.model"
    assert meta["PCR_bias_preset"] == "madritsch2025_r10"
    for nanosim_key in ("Coverage", "Min_read_length", "Max_read_length"):
        assert nanosim_key not in meta


def test_pacbio_amplicon_reports_model_passes_and_default_preset(tmp_path):
    config = {
        "read_simulation": {"simulator": "amplicon", "coverage": 100},
        "pacbio_params": {
            "coverage": 30,
            "pass_num": 5,
            "model_type": "errhmm",
            "model_file": "ERRHMM-SEQUEL.model",
        },
        "amplicon_params": {},
    }
    meta = _write(tmp_path, config, "PacBio")

    assert meta["Template_molecules"] == "100"
    assert meta["Pass_num"] == "5"
    assert meta["Model_type"] == "errhmm"
    assert meta["Model_file"] == "ERRHMM-SEQUEL.model"
    assert meta["PCR_bias_preset"] == "default"
    assert "Coverage" not in meta


@pytest.mark.parametrize(
    ("pcr_bias", "expected"),
    [({"alpha": 0.1}, "custom"), ({"preset": "no_bias", "stochastic": True}, "no_bias")],
)
def test_amplicon_pcr_bias_fields(tmp_path, pcr_bias, expected):
    config = {
        "read_simulation": {"simulator": "ont-amplicon", "coverage": 50},
        "amplicon_params": {"pcr_bias": pcr_bias},
    }
    meta = _write(tmp_path, config, "ONT")

    assert meta["PCR_bias_preset"] == expected
    assert meta["PCR_stochastic"] == str(bool(pcr_bias.get("stochastic")))


def test_pacbio_amplicon_template_count_falls_back_to_pacbio_coverage(tmp_path):
    config = {
        "read_simulation": {"simulator": "amplicon"},
        "pacbio_params": {"coverage": 40},
    }
    meta = _write(tmp_path, config, "PacBio")

    assert meta["Template_molecules"] == "40"


def test_non_amplicon_ont_metadata_unchanged(tmp_path):
    config = {"nanosim_params": {"coverage": 200, "min_read_length": 2000}}
    meta = _write(tmp_path, config, "ONT")

    assert meta["Coverage"] == "200"
    assert meta["Min_read_length"] == "2000"
    assert "Template_molecules" not in meta


def _write_extra(tmp_path: Path, config: dict, platform: str, extra) -> dict[str, str]:
    with (
        patch(
            "muc_one_up.read_simulator.utils.metadata_writer.capture_tool_versions",
            return_value={},
        ),
        patch("muc_one_up.read_simulator.utils.metadata_writer.log_tool_versions"),
    ):
        path = write_metadata_file(
            output_dir=str(tmp_path),
            output_base="x",
            config=config,
            start_time=datetime(2025, 1, 1, 12, 0, 0),
            end_time=datetime(2025, 1, 1, 12, 1, 0),
            platform=platform,
            extra_rows=extra,
        )
    rows = Path(path).read_text().splitlines()[1:]
    return dict(row.split("\t", 1) for row in rows)


def test_read_profile_and_truth_are_recorded(tmp_path):
    """Truth-tracked runs record the profile identity and truth manifest (#122)."""
    config = {
        "read_simulation": {"coverage": 500, "assay_type": "amplicon"},
        "read_model": {"profile": "p.json", "name": "prof", "sha256": "ab" * 32},
    }
    meta = _write_extra(tmp_path, config, "ONT", [("Read_truth", "x_read_truth.tsv.gz")])
    assert meta["Read_profile"] == "prof"
    assert meta["Read_profile_sha256"] == "ab" * 32
    assert meta["Read_truth"] == "x_read_truth.tsv.gz"


def test_fragment_simulator_reports_fragment_not_nanosim_rows(tmp_path):
    config = {
        "read_simulation": {"simulator": "ont-fragments", "coverage": 20},
        "nanosim_params": {"coverage": 20, "min_read_length": 100},
        "ont_fragment_params": {"n_reads": 800, "length_median": 6000, "length_sigma": 0.5},
    }
    meta = _write_extra(tmp_path, config, "ONT", None)
    assert "Min_read_length" not in meta
    assert meta["Fragment_reads"] == "800"
    assert meta["Fragment_length_median"] == "6000"
    assert meta["Fragment_length_sigma"] == "0.5"


def test_legacy_metadata_has_no_new_rows(tmp_path):
    config = {"nanosim_params": {"coverage": 20, "min_read_length": 100}}
    meta = _write_extra(tmp_path, config, "ONT", None)
    assert not {"Read_profile", "Read_truth", "Fragment_reads"} & meta.keys()
    assert meta["Min_read_length"] == "100"
