"""
Tests for exception handling across the codebase.

Verifies that:
1. Exception hierarchy is correct
2. sys.exit() has been removed from modules
3. Exceptions carry appropriate context
4. Error messages are informative
"""

from pathlib import Path
from unittest.mock import Mock

import pytest

from muc_one_up.cli.config import (
    parse_fixed_lengths,
    process_mutation_config,
    setup_configuration,
)
from muc_one_up.cli.haplotypes import generate_haplotypes
from muc_one_up.cli.mutations import (
    find_random_mutation_target,
    parse_mutation_targets,
)
from muc_one_up.exceptions import (
    ConfigurationError,
    ExternalToolError,
    FileOperationError,
    MucOneUpError,
    MutationError,
    ReadSimulationError,
    SimulationError,
    SNPIntegrationError,
    ValidationError,
)

# ==============================================================================
# Exception Hierarchy Tests
# ==============================================================================


def test_all_exceptions_inherit_from_base():
    """All custom exceptions should inherit from MucOneUpError."""
    assert issubclass(ConfigurationError, MucOneUpError)
    assert issubclass(ValidationError, MucOneUpError)
    assert issubclass(SimulationError, MucOneUpError)
    assert issubclass(MutationError, MucOneUpError)
    assert issubclass(SNPIntegrationError, MucOneUpError)
    assert issubclass(FileOperationError, MucOneUpError)
    assert issubclass(ReadSimulationError, MucOneUpError)
    assert issubclass(ExternalToolError, MucOneUpError)


def test_all_exceptions_inherit_from_exception():
    """All custom exceptions should inherit from Exception."""
    assert issubclass(MucOneUpError, Exception)
    assert issubclass(ConfigurationError, Exception)


def test_exception_instantiation():
    """Exceptions should be instantiable with message."""
    msg = "Test error message"

    exc = ConfigurationError(msg)
    assert str(exc) == msg
    assert isinstance(exc, Exception)


# ==============================================================================
# ConfigurationError Tests
# ==============================================================================


def test_setup_configuration_missing_file():
    """setup_configuration should raise ConfigurationError for missing file."""
    mock_args = Mock()
    mock_args.config = "/nonexistent/config.json"
    mock_args.out_dir = "/tmp/test"
    mock_args.out_base = "test"
    mock_args.reference_assembly = None

    with pytest.raises(ConfigurationError, match="Config file not found"):
        setup_configuration(mock_args)


def test_setup_configuration_invalid_json(tmp_path):
    """setup_configuration should raise ConfigurationError for invalid JSON."""
    invalid_json = tmp_path / "invalid.json"
    invalid_json.write_text("{ this is not valid JSON }")

    mock_args = Mock()
    mock_args.config = str(invalid_json)
    mock_args.out_dir = str(tmp_path / "output")
    mock_args.out_base = "test"
    mock_args.reference_assembly = None

    with pytest.raises(ConfigurationError, match="Invalid JSON"):
        setup_configuration(mock_args)


# ==============================================================================
# ValidationError Tests
# ==============================================================================


def test_parse_fixed_lengths_invalid_count():
    """parse_fixed_lengths should raise ValidationError for wrong count."""
    with pytest.raises(ValidationError, match="must have either 1 value or N="):
        parse_fixed_lengths(["10", "20"], num_haplotypes=3)


def test_parse_fixed_lengths_invalid_format():
    """parse_fixed_lengths should raise ValidationError for invalid format."""
    with pytest.raises(ValidationError, match="Invalid fixed-length value"):
        parse_fixed_lengths(["abc"], num_haplotypes=1)


def test_process_mutation_config_dual_invalid_first():
    """process_mutation_config should raise ValidationError for invalid dual mode."""
    mock_args = Mock()
    mock_args.mutation_name = "dupC,delA"  # Should start with 'normal'
    mock_args.mutation_targets = None

    with pytest.raises(ValidationError, match="first mutation-name must be 'normal'"):
        process_mutation_config(mock_args, None)


def test_parse_mutation_targets_invalid():
    """parse_mutation_targets should raise ValidationError for invalid format."""
    with pytest.raises(ValidationError, match="Invalid --mutation-targets format"):
        parse_mutation_targets(["invalid"])


# ==============================================================================
# SimulationError Tests
# ==============================================================================


def test_generate_haplotypes_simulation_error():
    """generate_haplotypes should raise SimulationError on failure."""
    mock_args = Mock()
    mock_args.num_haplotypes = 2
    mock_args.seed = None

    # Invalid config that will cause simulation to fail
    invalid_config = {}  # Missing required fields

    with pytest.raises(SimulationError, match="simulation failed"):
        generate_haplotypes(mock_args, invalid_config, fixed_conf=None, predefined_chains=None)


# ==============================================================================
# MutationError Tests
# ==============================================================================


def test_find_random_mutation_target_not_found():
    """find_random_mutation_target should raise MutationError when mutation doesn't exist."""
    results = [("ATCG", ["A", "B", "C"])]
    config = {"mutations": {"dupC": {"allowed_repeats": ["C"]}}}

    with pytest.raises(MutationError, match="not in config\\['mutations'\\]"):
        find_random_mutation_target(results, config, "nonexistent")


def test_find_random_mutation_target_no_valid_repeats():
    """find_random_mutation_target should raise MutationError when no valid targets."""
    results = [("ATCG", ["A", "B", "C"])]
    config = {"mutations": {"test": {"allowed_repeats": ["X", "Y", "Z"]}}}

    with pytest.raises(MutationError, match="No repeats match"):
        find_random_mutation_target(results, config, "test")


# ==============================================================================
# ExternalToolError Tests
# ==============================================================================


def test_external_tool_error_attributes():
    """ExternalToolError should preserve tool execution details."""
    error = ExternalToolError(
        tool="samtools",
        exit_code=1,
        stderr="Error: file not found",
        stdout="",
        cmd="samtools view input.bam",
    )

    assert error.tool == "samtools"
    assert error.exit_code == 1
    assert error.stderr == "Error: file not found"
    assert error.cmd == "samtools view input.bam"
    assert "samtools failed with exit code 1" in str(error)
    assert "samtools view input.bam" in str(error)


def test_external_tool_error_truncates_long_stderr():
    """ExternalToolError should truncate very long stderr output."""
    long_stderr = "X" * 1000
    error = ExternalToolError(
        tool="test_tool",
        exit_code=1,
        stderr=long_stderr,
    )

    error_str = str(error)
    # Should be truncated with "..."
    assert len(error_str) < len(long_stderr)
    assert "..." in error_str


# ==============================================================================
# sys.exit() Removal Verification
# ==============================================================================


def test_no_sys_exit_in_cli_modules():
    """CLI modules should not have sys.exit() calls (except main entry point)."""
    violations = []
    cli_dir = Path("muc_one_up/cli")

    for py_file in cli_dir.glob("*.py"):
        with py_file.open() as f:
            content = f.read()
            lines = content.split("\n")

            for line_num, line in enumerate(lines, 1):
                # Skip if it's in the __main__ block (check current and previous lines)
                if "sys.exit" in line:
                    # Check if we're in the if __main__ block
                    in_main_block = False
                    for prev_line_num in range(max(0, line_num - 10), line_num):
                        if '__name__ == "__main__"' in lines[prev_line_num]:
                            in_main_block = True
                            break

                    # Skip if in __main__ block or has OK comment
                    if not in_main_block and "# OK:" not in line:
                        violations.append(f"{py_file}:{line_num}: {line.strip()}")

    assert not violations, "Found sys.exit() in CLI modules:\\n" + "\\n".join(violations)


def test_exception_chain_preserved():
    """Exceptions should preserve the original exception chain with 'from e'."""
    # Test ConfigurationError preserves chain
    try:
        try:
            raise ValueError("Original error")
        except ValueError as e:
            raise ConfigurationError("Wrapped error") from e
    except ConfigurationError as exc:
        assert exc.__cause__ is not None
        assert isinstance(exc.__cause__, ValueError)
        assert str(exc.__cause__) == "Original error"


# ==============================================================================
# Error Message Quality Tests
# ==============================================================================


def test_configuration_error_message_quality():
    """ConfigurationError messages should be informative."""
    error = ConfigurationError("Config file not found: /path/to/config.json")
    msg = str(error)

    # Should include what went wrong
    assert "Config file not found" in msg
    # Should include the path
    assert "/path/to/config.json" in msg


def test_mutation_error_message_includes_available_mutations():
    """MutationError should list available mutations."""
    results = [("ATCG", ["A"])]
    config = {"mutations": {"dupC": {"allowed_repeats": ["C"]}, "delA": {"allowed_repeats": ["A"]}}}

    with pytest.raises(MutationError) as exc_info:
        find_random_mutation_target(results, config, "nonexistent")

    error_msg = str(exc_info.value)
    # Should mention available mutations
    assert "dupC" in error_msg or "delA" in error_msg


def test_external_tool_error_message_includes_command():
    """ExternalToolError should include the command that failed."""
    error = ExternalToolError(
        tool="bwa", exit_code=1, stderr="alignment failed", cmd="bwa mem ref.fa reads.fq"
    )

    msg = str(error)
    assert "bwa mem ref.fa reads.fq" in msg
    assert "alignment failed" in msg


# ==============================================================================
# Catch-All Tests
# ==============================================================================


def test_can_catch_all_custom_exceptions():
    """All custom exceptions should be catchable with MucOneUpError."""
    exceptions_to_test = [
        ConfigurationError("test"),
        ValidationError("test"),
        SimulationError("test"),
        MutationError("test"),
        SNPIntegrationError("test"),
        FileOperationError("test"),
        ReadSimulationError("test"),
        ExternalToolError("test", 1, "error"),
    ]

    for exc in exceptions_to_test:
        try:
            raise exc
        except MucOneUpError as e:
            assert isinstance(e, MucOneUpError)
            # Verify it's also the specific type
            assert type(e) in [
                ConfigurationError,
                ValidationError,
                SimulationError,
                MutationError,
                SNPIntegrationError,
                FileOperationError,
                ReadSimulationError,
                ExternalToolError,
            ]
