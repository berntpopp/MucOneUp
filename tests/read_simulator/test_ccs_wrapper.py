"""Tests for CCS wrapper - focuses on OUR command construction logic.

Following Phase 2 testing principles:
- Mock ONLY at system boundary (run_command)
- Test OUR code's logic (command construction, error handling, file validation)
- NOT testing CCS tool itself

These tests verify that our wrapper correctly:
1. Constructs CCS commands with proper arguments
2. Validates input BAM files
3. Detects empty/corrupted output
4. Handles quality thresholds correctly
5. Validates output files are created
"""

import pytest

from muc_one_up.exceptions import FileOperationError
from muc_one_up.read_simulator.wrappers.ccs_wrapper import (
    run_ccs_consensus,
    validate_ccs_parameters,
)


class TestRunCcsConsensus:
    """Test run_ccs_consensus command construction."""

    def test_constructs_basic_command(self, mocker, tmp_path):
        """Test that CCS command is constructed correctly with all parameters."""
        # Arrange
        input_bam = tmp_path / "input_clr.bam"
        output_bam = tmp_path / "output_hifi.bam"

        input_bam.write_text("CLR BAM data")

        # Mock run_command to create output BAM
        def mock_run_cmd(cmd, **kwargs):
            if "ccs" in str(cmd):
                # Create valid-sized output
                output_bam.write_bytes(b"HiFi BAM data" * 100)  # > 1KB

        mocker.patch(
            "muc_one_up.read_simulator.wrappers.ccs_wrapper.run_command", side_effect=mock_run_cmd
        )

        # Act
        result = run_ccs_consensus(
            ccs_cmd="ccs",
            input_bam=str(input_bam),
            output_bam=str(output_bam),
            min_passes=3,
            min_rq=0.99,
            threads=8,
            seed=42,
        )

        # Assert: Verify command arguments
        from muc_one_up.read_simulator.wrappers.ccs_wrapper import run_command

        cmd = run_command.call_args[0][0]

        assert str(input_bam) in cmd
        assert str(output_bam) in cmd
        assert "--min-passes" in cmd
        assert "3" in cmd
        assert "--min-rq" in cmd
        assert "0.99" in cmd
        assert "--num-threads" in cmd
        assert "8" in cmd
        assert "--seed" in cmd
        assert "42" in cmd

        # Verify output path returned
        assert result == str(output_bam)

    def test_validates_input_bam_exists(self, tmp_path):
        """Test that missing input BAM is detected."""
        # Arrange
        missing_bam = tmp_path / "nonexistent.bam"
        output_bam = tmp_path / "output.bam"

        # Act & Assert
        with pytest.raises(FileOperationError, match="Input multi-pass CLR BAM not found"):
            run_ccs_consensus(
                ccs_cmd="ccs",
                input_bam=str(missing_bam),  # Doesn't exist!
                output_bam=str(output_bam),
            )

    def test_validates_input_bam_not_empty(self, tmp_path):
        """Test that empty input BAM is detected."""
        # Arrange
        input_bam = tmp_path / "empty.bam"
        output_bam = tmp_path / "output.bam"

        input_bam.write_text("")  # Empty!

        # Act & Assert
        with pytest.raises(FileOperationError, match="Input multi-pass CLR BAM is empty"):
            run_ccs_consensus(
                ccs_cmd="ccs",
                input_bam=str(input_bam),
                output_bam=str(output_bam),
            )

    def test_detects_small_output_bam(self, mocker, tmp_path):
        """Test that suspiciously small output BAM is detected."""
        # Arrange
        input_bam = tmp_path / "input.bam"
        output_bam = tmp_path / "output.bam"

        input_bam.write_text("CLR data")

        # Mock run_command to create tiny output (< 1KB)
        def mock_run_cmd(cmd, **kwargs):
            if "ccs" in str(cmd):
                output_bam.write_bytes(b"tiny")  # Only 4 bytes!

        mocker.patch(
            "muc_one_up.read_simulator.wrappers.ccs_wrapper.run_command", side_effect=mock_run_cmd
        )

        # Act & Assert
        with pytest.raises(FileOperationError, match="CCS produced empty or invalid output"):
            run_ccs_consensus(
                ccs_cmd="ccs",
                input_bam=str(input_bam),
                output_bam=str(output_bam),
            )

    def test_handles_conda_command(self, mocker, tmp_path):
        """Test that conda/mamba commands are handled via build_tool_command."""
        # Arrange
        input_bam = tmp_path / "input.bam"
        output_bam = tmp_path / "output.bam"

        input_bam.write_text("CLR data")

        # Mock run_command
        def mock_run_cmd(cmd, **kwargs):
            # Verify cmd is a list (not string)
            assert isinstance(cmd, list)
            output_bam.write_bytes(b"HiFi data" * 200)

        mocker.patch(
            "muc_one_up.read_simulator.wrappers.ccs_wrapper.run_command", side_effect=mock_run_cmd
        )

        # Act
        run_ccs_consensus(
            ccs_cmd="mamba run -n env_pacbio ccs",  # Multi-word command
            input_bam=str(input_bam),
            output_bam=str(output_bam),
        )

        # Assert: Command was constructed as list
        from muc_one_up.read_simulator.wrappers.ccs_wrapper import run_command

        cmd = run_command.call_args[0][0]
        assert isinstance(cmd, list)

    def test_omits_seed_when_none(self, mocker, tmp_path):
        """Test that seed argument is omitted when None."""
        # Arrange
        input_bam = tmp_path / "input.bam"
        output_bam = tmp_path / "output.bam"

        input_bam.write_text("CLR data")

        # Mock run_command
        def mock_run_cmd(cmd, **kwargs):
            output_bam.write_bytes(b"HiFi data" * 200)

        mocker.patch(
            "muc_one_up.read_simulator.wrappers.ccs_wrapper.run_command", side_effect=mock_run_cmd
        )

        # Act
        run_ccs_consensus(
            ccs_cmd="ccs",
            input_bam=str(input_bam),
            output_bam=str(output_bam),
            seed=None,  # No seed
        )

        # Assert: No seed in command
        from muc_one_up.read_simulator.wrappers.ccs_wrapper import run_command

        cmd = run_command.call_args[0][0]
        assert "--seed" not in cmd


class TestValidateCcsParameters:
    """Test validate_ccs_parameters validation logic."""

    def test_validates_min_passes_positive(self):
        """Test min_passes must be ≥ 1."""
        # Valid
        validate_ccs_parameters(min_passes=3, min_rq=0.99)

        # Invalid - zero
        with pytest.raises(ValueError, match="Invalid min_passes"):
            validate_ccs_parameters(min_passes=0, min_rq=0.99)

        # Invalid - negative
        with pytest.raises(ValueError, match="Invalid min_passes"):
            validate_ccs_parameters(min_passes=-1, min_rq=0.99)

    def test_validates_min_rq_range(self):
        """Test min_rq must be in [0.0, 1.0]."""
        # Valid
        validate_ccs_parameters(min_passes=3, min_rq=0.99)

        # Too low
        with pytest.raises(ValueError, match="Invalid min_rq"):
            validate_ccs_parameters(min_passes=3, min_rq=-0.1)

        # Too high
        with pytest.raises(ValueError, match="Invalid min_rq"):
            validate_ccs_parameters(min_passes=3, min_rq=1.5)

    def test_warns_for_low_min_rq(self, mocker):
        """Test warning is logged for min_rq below standard HiFi threshold."""
        # Arrange
        mock_warning = mocker.patch("logging.warning")

        # Act
        validate_ccs_parameters(min_passes=3, min_rq=0.95)  # Below 0.99

        # Assert
        mock_warning.assert_called_once()
        assert "below standard HiFi threshold" in str(mock_warning.call_args)

    def test_warns_for_high_min_rq(self, mocker):
        """Test warning is logged for very high min_rq."""
        # Arrange
        mock_warning = mocker.patch("logging.warning")

        # Act
        validate_ccs_parameters(min_passes=3, min_rq=0.9999)  # Very high

        # Assert
        mock_warning.assert_called_once()
        assert "Very high min_rq value" in str(mock_warning.call_args)

    def test_warns_for_high_min_passes(self, mocker):
        """Test warning is logged for very high min_passes."""
        # Arrange
        mock_warning = mocker.patch("logging.warning")

        # Act
        validate_ccs_parameters(min_passes=51, min_rq=0.99)  # Very high

        # Assert
        mock_warning.assert_called_once()
        assert "Very high min_passes value" in str(mock_warning.call_args)
