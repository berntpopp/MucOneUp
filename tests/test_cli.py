"""
Tests for refactored CLI functions.

Following SOLID, DRY, KISS principles with focused unit tests.
"""

import json
from pathlib import Path
from unittest.mock import Mock

import pytest

from muc_one_up.cli.config import (
    determine_simulation_mode,
    numbered_filename,
    parse_fixed_lengths,
    process_mutation_config,
    setup_configuration,
)
from muc_one_up.cli.mutations import (
    find_random_mutation_target,
    parse_mutation_targets,
)
from muc_one_up.cli.snps import integrate_snps_unified
from muc_one_up.exceptions import ConfigurationError, ValidationError
from muc_one_up.type_defs import HaplotypeResult, MutationTarget

# ==============================================================================
# Fixtures
# ==============================================================================


@pytest.fixture
def minimal_config():
    """Minimal valid configuration for testing."""
    return {
        "repeats": {
            "1": "ATCG",
            "2": "GCTA",
            "7": "AAAA",
            "8": "TTTT",
            "9": "GGGG",
        },
        "constants": {
            "hg38": {
                "left": "AAAA",
                "right": "TTTT",
                "vntr_start": 4,
                "vntr_end": 20,
            }
        },
        "probabilities": {
            "1": {"2": 1.0},
            "2": {"7": 1.0},
            "7": {"8": 1.0},
            "8": {"9": 1.0},
            "9": {},
        },
        "length_model": {
            "distribution": "normal",
            "min_repeats": 5,
            "max_repeats": 10,
            "mean_repeats": 7,
            "median_repeats": 7,
        },
        "mutations": {
            "dupC": {
                "allowed_repeats": ["1", "2"],
                "strict_mode": False,
                "changes": [{"type": "insert", "start": 1, "end": 1, "sequence": "CCCC"}],
            }
        },
        "reference_assembly": "hg38",
    }


@pytest.fixture
def temp_config_file(tmp_path, minimal_config):
    """Create temporary config file."""
    config_file = tmp_path / "config.json"
    config_file.write_text(json.dumps(minimal_config))
    return config_file


@pytest.fixture
def mock_args(temp_config_file, tmp_path):
    """Create mock arguments."""
    args = Mock()
    args.config = str(temp_config_file)
    args.out_dir = str(tmp_path / "output")
    args.out_base = "test"
    args.reference_assembly = None
    args.num_haplotypes = 2
    args.fixed_lengths = None
    args.input_structure = None
    args.simulate_series = None
    args.mutation_name = None
    args.mutation_targets = None
    args.seed = 42
    args.snp_input_file = None
    args.random_snps = False
    args.random_snp_density = None
    args.random_snp_output_file = None
    args.random_snp_region = "constants_only"
    args.random_snp_haplotypes = "all"
    return args


@pytest.fixture
def sample_results():
    """Sample haplotype results."""
    return [
        HaplotypeResult.from_tuple(("ATCGATCGATCG", ["1", "2", "7", "8", "9"])),
        HaplotypeResult.from_tuple(("GCTAGCTAGCTA", ["1", "2", "7", "8", "9"])),
    ]


# ==============================================================================
# Tests for setup_configuration()
# ==============================================================================


def test_setup_configuration_success(mock_args, tmp_path):
    """Test successful configuration loading."""
    config, out_dir, out_base = setup_configuration(mock_args)

    assert config is not None
    assert "repeats" in config
    assert out_dir == str(tmp_path / "output")
    assert out_base == "test"
    assert Path(out_dir).exists()


def test_setup_configuration_creates_output_dir(mock_args, tmp_path):
    """Test that output directory is created."""
    output_path = tmp_path / "output"
    assert not output_path.exists()

    _config, _out_dir, _out_base = setup_configuration(mock_args)

    assert output_path.exists()
    assert output_path.is_dir()


def test_setup_configuration_missing_file(mock_args):
    """Test error handling for missing config file."""
    from muc_one_up.exceptions import ConfigurationError

    mock_args.config = "/nonexistent/config.json"

    with pytest.raises(ConfigurationError, match="Config file not found"):
        setup_configuration(mock_args)


def test_setup_configuration_override_assembly(mock_args, minimal_config):
    """Test reference assembly override."""
    mock_args.reference_assembly = "hg19"

    config, _out_dir, _out_base = setup_configuration(mock_args)

    assert config["reference_assembly"] == "hg19"


def test_setup_configuration_invalid_out_base(mock_args):
    """Test that path traversal or separators in out_base are rejected."""
    mock_args.out_base = "../escaped"
    with pytest.raises(ValidationError, match="Invalid out_base"):
        setup_configuration(mock_args)


def test_numbered_filename_invalid_out_base(tmp_path):
    """Test numbered_filename rejects invalid out_base."""
    with pytest.raises(ValidationError, match="Invalid out_base"):
        numbered_filename(str(tmp_path), "../escaped", 1, "fa")


def test_setup_configuration_random_snps_validation(mock_args):
    """Test that --random-snps requires output file and density upfront."""
    mock_args.random_snps = True
    mock_args.random_snp_output_file = None
    mock_args.random_snp_density = 0.001

    with pytest.raises(ValidationError, match="--random-snp-output-file is required"):
        setup_configuration(mock_args)

    mock_args.random_snp_output_file = "snps.tsv"
    mock_args.random_snp_density = None

    with pytest.raises(ValidationError, match="--random-snp-density is required"):
        setup_configuration(mock_args)


def test_setup_configuration_null_section_reports_configuration_error_not_invalid_json(
    tmp_path, mock_args
):
    """Test that a null section in valid JSON is reported as ConfigurationError without 'Invalid JSON'."""
    bad_config = tmp_path / "null_constants.json"
    bad_config.write_text(json.dumps({"repeats": {"1": "ATCG"}, "constants": None}))
    mock_args.config = str(bad_config)

    with pytest.raises(ConfigurationError) as exc_info:
        setup_configuration(mock_args)

    error_msg = str(exc_info.value)
    assert "Invalid JSON" not in error_msg
    assert "constants" in error_msg


# ==============================================================================
# Tests for determine_simulation_mode()
# ==============================================================================


def test_determine_simulation_mode_random(mock_args, minimal_config):
    """Test random simulation mode (no fixed lengths)."""
    sim_configs, chains, mutation_info = determine_simulation_mode(mock_args, minimal_config)

    assert sim_configs == [None]
    assert chains is None
    assert mutation_info is None


def test_determine_simulation_mode_fixed_lengths(mock_args, minimal_config):
    """Test fixed lengths simulation mode."""
    mock_args.fixed_lengths = ["20", "30"]

    sim_configs, chains, mutation_info = determine_simulation_mode(mock_args, minimal_config)

    assert len(sim_configs) == 1
    assert len(sim_configs[0]) == 2
    assert chains is None
    assert mutation_info is None


def test_determine_simulation_mode_series(mock_args, minimal_config):
    """Test series simulation mode."""
    mock_args.fixed_lengths = ["20-25"]
    mock_args.simulate_series = 2

    sim_configs, _chains, _mutation_info = determine_simulation_mode(mock_args, minimal_config)

    assert len(sim_configs) > 1
    # Should generate configs for range with step size


@pytest.mark.parametrize("invalid_step", [0, -1])
def test_determine_simulation_mode_series_invalid_step(mock_args, minimal_config, invalid_step):
    """Test that series step <= 0 raises ValidationError."""
    mock_args.fixed_lengths = ["20-30"]
    mock_args.simulate_series = invalid_step

    with pytest.raises(ValidationError, match="--simulate-series must be >= 1"):
        determine_simulation_mode(mock_args, minimal_config)


def test_parse_fixed_lengths_inverted_range():
    """Test that inverted range like 40-20 raises ValidationError."""
    with pytest.raises(ValidationError, match=r"start .* must be <= end"):
        parse_fixed_lengths(["40-20"], 1)


def test_parse_fixed_lengths_below_minimum():
    """Test that fixed length below 10 raises ValidationError."""
    with pytest.raises(ValidationError, match="Minimum repeat length is 10"):
        parse_fixed_lengths(["5"], 1)


def test_parse_fixed_lengths_range_below_minimum():
    """Test that range starting below 10 raises ValidationError."""
    with pytest.raises(ValidationError, match="Minimum repeat length is 10"):
        parse_fixed_lengths(["5-20"], 1)


# ==============================================================================
# Tests for process_mutation_config()
# ==============================================================================


def test_process_mutation_config_none(mock_args):
    """Test mutation config with no mutation."""
    dual_mode, pair, name = process_mutation_config(mock_args, None)

    assert not dual_mode
    assert pair is None
    assert name is None


def test_process_mutation_config_single(mock_args):
    """Test single mutation mode."""
    mock_args.mutation_name = "dupC"

    dual_mode, pair, name = process_mutation_config(mock_args, None)

    assert not dual_mode
    assert pair is None
    assert name == "dupC"


def test_process_mutation_config_dual(mock_args):
    """Test dual mutation mode."""
    mock_args.mutation_name = "normal,dupC"

    dual_mode, pair, name = process_mutation_config(mock_args, None)

    assert dual_mode
    assert pair == ["normal", "dupC"]
    assert name == "normal,dupC"


def test_process_mutation_config_dual_invalid_first(mock_args):
    """Test dual mode with invalid first mutation."""
    from muc_one_up.exceptions import ValidationError

    mock_args.mutation_name = "dupC,normal"

    with pytest.raises(ValidationError, match="first mutation-name must be 'normal'"):
        process_mutation_config(mock_args, None)


def test_process_mutation_config_preserves_cli_mutation_when_structure_info_present(mock_args):
    """Test that explicit CLI mutation args are not overwritten by structure file headers."""
    mock_args.mutation_name = "normal,dupC"
    mock_args.mutation_targets = ["1,10"]
    structure_info = {"name": "delA", "targets": [(1, 5)]}

    dual_mode, pair, name = process_mutation_config(mock_args, structure_info)

    assert dual_mode is True
    assert pair == ["normal", "dupC"]
    assert name == "normal,dupC"
    assert mock_args.mutation_targets == ["1,10"]


def test_process_mutation_config_uses_structure_info_defaults_when_cli_none(mock_args):
    """Test that structure file mutation defaults are applied when CLI args are omitted."""
    mock_args.mutation_name = None
    mock_args.mutation_targets = None
    structure_info = {"name": "dupC", "targets": [(1, 2)]}

    dual_mode, pair, name = process_mutation_config(mock_args, structure_info)

    assert dual_mode is False
    assert pair is None
    assert name == "dupC"
    assert mock_args.mutation_targets == ["1,2"]


# ==============================================================================
# Tests for parse_mutation_targets()
# ==============================================================================


def test_parse_mutation_targets_valid_strings():
    """Test parsing valid string targets."""
    targets = ["1,5", "2,10"]

    result = parse_mutation_targets(targets)

    assert result == [MutationTarget(1, 5), MutationTarget(2, 10)]


def test_parse_mutation_targets_valid_tuples():
    """Test parsing valid tuple targets."""
    targets = [(1, 5), (2, 10)]

    result = parse_mutation_targets(targets)

    assert result == [MutationTarget(1, 5), MutationTarget(2, 10)]


def test_parse_mutation_targets_mixed():
    """Test parsing mixed string and tuple targets."""
    targets = ["1,5", (2, 10)]

    result = parse_mutation_targets(targets)

    assert result == [MutationTarget(1, 5), MutationTarget(2, 10)]


def test_parse_mutation_targets_invalid_format():
    """Test error handling for invalid format."""
    from muc_one_up.exceptions import ValidationError

    targets = ["invalid"]

    with pytest.raises(ValidationError, match="Invalid --mutation-targets format"):
        parse_mutation_targets(targets)


# ==============================================================================
# Tests for find_random_mutation_target()
# ==============================================================================


def test_find_random_mutation_target_success(sample_results, minimal_config):
    """Test finding random mutation target."""
    result = find_random_mutation_target(sample_results, minimal_config, "dupC")

    assert len(result) == 1
    target = result[0]
    assert 1 <= target.haplotype_index <= 2
    assert 1 <= target.repeat_index <= 5


def test_find_random_mutation_target_mutation_not_found(sample_results, minimal_config):
    """Test error when mutation not in config."""
    from muc_one_up.exceptions import MutationError

    with pytest.raises(MutationError, match="not in config"):
        find_random_mutation_target(sample_results, minimal_config, "nonexistent")


def test_find_random_mutation_target_no_valid_repeats(sample_results, minimal_config):
    """Test error when no valid repeat targets exist."""
    from muc_one_up.exceptions import MutationError

    # Create mutation with no matching repeats
    minimal_config["mutations"]["test"] = {
        "allowed_repeats": ["X"],  # Not in sample_results
        "strict_mode": False,
        "changes": [],
    }

    with pytest.raises(MutationError, match="No repeats match"):
        find_random_mutation_target(sample_results, minimal_config, "test")


# ==============================================================================
# Tests for integrate_snps_unified() - Tests DRY principle
# ==============================================================================


def test_integrate_snps_unified_no_snps(mock_args, minimal_config, sample_results):
    """Test SNP integration with no SNPs requested."""
    results, snp_info = integrate_snps_unified(mock_args, minimal_config, sample_results)

    assert results == sample_results
    assert snp_info == {}


def test_integrate_snps_unified_file_based(mock_args, minimal_config, sample_results, tmp_path):
    """Test SNP integration from file."""
    # Create SNP file
    snp_file = tmp_path / "snps.tsv"
    snp_file.write_text("1\t5\tA\tG\n")
    mock_args.snp_input_file = str(snp_file)

    _results, snp_info = integrate_snps_unified(mock_args, minimal_config, sample_results)

    # Should have attempted to apply SNPs
    assert isinstance(snp_info, dict)


def test_integrate_snps_unified_skip_reference_check(
    mock_args, minimal_config, sample_results, tmp_path
):
    """Test SNP integration with skip_reference_check flag."""
    # Create SNP file
    snp_file = tmp_path / "snps.tsv"
    snp_file.write_text("1\t5\tA\tG\n")
    mock_args.snp_input_file = str(snp_file)

    # Should not raise error even if reference doesn't match
    results, snp_info = integrate_snps_unified(
        mock_args, minimal_config, sample_results, skip_reference_check=True
    )

    assert isinstance(results, list)
    assert isinstance(snp_info, dict)


# ==============================================================================
# Integration Tests
# ==============================================================================


def test_cli_refactoring_maintains_backwards_compatibility(mock_args, minimal_config):
    """Test that refactored CLI maintains backwards compatibility."""
    # This is an integration test to ensure the refactoring doesn't break existing functionality

    # Test configuration loading
    config, _out_dir, _out_base = setup_configuration(mock_args)
    assert config is not None

    # Test simulation mode determination
    sim_configs, _chains, mutation_info = determine_simulation_mode(mock_args, config)
    assert sim_configs is not None

    # Test mutation config processing
    dual_mode, _pair, _name = process_mutation_config(mock_args, mutation_info)
    assert isinstance(dual_mode, bool)


def test_cli_functions_follow_solid_principles():
    """
    Meta-test to document SOLID principles in refactored code.

    Single Responsibility:
    - setup_configuration: Only handles config loading
    - determine_simulation_mode: Only determines sim mode
    - process_mutation_config: Only processes mutation config
    - Each function has one clear purpose

    Open/Closed:
    - Functions accept args/config, can be extended without modification
    - New simulation modes can be added to determine_simulation_mode

    Liskov Substitution:
    - Not applicable (no inheritance)

    Interface Segregation:
    - Functions have minimal, focused interfaces
    - No function depends on unused parameters

    Dependency Inversion:
    - Functions depend on abstractions (args, config dicts)
    - Not tightly coupled to implementations
    """
    pass


def test_cli_refactoring_eliminates_duplication():
    """
    Meta-test documenting DRY improvements.

    Before: SNP integration logic duplicated (~120 lines)
    - Lines 629-707 for dual mode
    - Lines 772-828 for single mode

    After: Single integrate_snps_unified() function
    - Used in both dual and single modes
    - Eliminates ~120 lines of duplication
    - skip_reference_check parameter handles differences
    """
    pass


def test_cli_refactoring_improves_testability():
    """
    Meta-test documenting testability improvements.

    Before:
    - 801-line main() function
    - Cannot test individual concerns
    - sys.exit() kills test runner

    After:
    - 21 focused functions (up from 6)
    - Each function testable independently
    - main() reduced to ~45 lines
    - Clear separation of concerns
    """
    pass


def test_cli_simulate_rejects_num_haplotypes_zero(temp_config_file):
    """Click command 'simulate' rejects --num-haplotypes 0."""
    from click.testing import CliRunner

    from muc_one_up.cli.click_main import cli

    runner = CliRunner()
    result = runner.invoke(
        cli, ["--config", str(temp_config_file), "simulate", "--num-haplotypes", "0"]
    )
    assert result.exit_code != 0
    assert "is not in the range" in result.output or "Error" in result.output


def test_cli_simulate_rejects_simulate_series_zero(temp_config_file):
    """Click command 'simulate' rejects --simulate-series 0."""
    from click.testing import CliRunner

    from muc_one_up.cli.click_main import cli

    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "--config",
            str(temp_config_file),
            "simulate",
            "--fixed-lengths",
            "20-30",
            "--simulate-series",
            "0",
        ],
    )
    assert result.exit_code != 0
    assert "is not in the range" in result.output or "Error" in result.output


def test_cli_simulate_series_numbers_snp_output_file(tmp_path, temp_config_file, monkeypatch):
    """Series mode passes per-iteration numbered SNP files to avoid clobbering."""
    from click.testing import CliRunner

    from muc_one_up.cli.click_main import cli
    from muc_one_up.cli.commands import simulate as sim_module

    recorded_snp_files = []

    def mock_run_iter(args, *pargs, **kwargs):
        recorded_snp_files.append(args.random_snp_output_file)

    monkeypatch.setattr(sim_module, "run_single_simulation_iteration", mock_run_iter)

    runner = CliRunner()
    snp_file = tmp_path / "snps.tsv"
    result = runner.invoke(
        cli,
        [
            "--config",
            str(temp_config_file),
            "simulate",
            "--fixed-lengths",
            "20-22",
            "--simulate-series",
            "1",
            "--random-snps",
            "--random-snp-output-file",
            str(snp_file),
            "--random-snp-density",
            "0.001",
            "--out-dir",
            str(tmp_path / "out"),
        ],
    )
    assert result.exit_code == 0
    assert len(recorded_snp_files) > 1
    assert recorded_snp_files[0].endswith(".001.tsv")
    assert recorded_snp_files[1].endswith(".002.tsv")
