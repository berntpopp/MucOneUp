"""Tests for pbsim3 wrapper - focuses on OUR command construction logic.

Following Phase 2 testing principles:
- Mock ONLY at system boundary (subprocess.Popen, run_command)
- Test OUR code's logic (command construction, error handling, file management)
- NOT testing pbsim3 itself

These tests verify that our wrapper correctly:
1. Constructs pbsim3 commands with proper arguments
2. Handles conda/mamba commands using build_tool_command
3. Validates model files and parameters
4. Manages SAM→BAM conversion when needed
5. Validates output files are created
"""

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from muc_one_up.exceptions import FileOperationError
from muc_one_up.read_simulator.wrappers.pbsim3_wrapper import (
    run_pbsim3_simulation,
    run_pbsim3_template_simulation,
    validate_pbsim3_parameters,
)


class TestRunPbsim3Simulation:
    """Test run_pbsim3_simulation command construction."""

    def test_constructs_basic_command(self, mocker, tmp_path):
        """Test that pbsim3 command is constructed correctly with all parameters."""
        # Arrange
        ref_fa = tmp_path / "ref.fa"
        model_file = tmp_path / "QSHMM-SEQUEL.model"
        output_prefix = tmp_path / "output"
        output_bam = tmp_path / "output.bam"

        ref_fa.write_text(">chr1\nACGT\n")
        model_file.write_text("model data")

        # Mock run_command to create output BAM
        def mock_run_cmd(cmd, **kwargs):
            if "pbsim" in str(cmd):
                output_bam.write_text("BAM data")

        mocker.patch(
            "muc_one_up.read_simulator.wrappers.pbsim3_wrapper.run_command",
            side_effect=mock_run_cmd,
        )

        # Act
        result = run_pbsim3_simulation(
            pbsim3_cmd="pbsim",
            samtools_cmd="samtools",
            reference=str(ref_fa),
            model_type="qshmm",
            model_file=str(model_file),
            coverage=30.0,
            output_prefix=str(output_prefix),
            pass_num=3,
            seed=42,
        )

        # Assert: Verify command arguments
        from muc_one_up.read_simulator.wrappers.pbsim3_wrapper import run_command

        cmd = run_command.call_args[0][0]

        assert "pbsim" in cmd
        assert "--strategy" in cmd
        assert "wgs" in cmd
        assert "--method" in cmd
        assert "qshmm" in cmd
        assert "--qshmm" in cmd
        assert str(model_file) in cmd
        assert "--depth" in cmd
        assert "--pass-num" in cmd
        assert "3" in cmd
        assert "--seed" in cmd
        assert "42" in cmd
        assert "--prefix" in cmd
        assert str(output_prefix) in cmd
        assert str(ref_fa) in cmd

        # Verify output path returned (now returns list)
        assert isinstance(result, list)
        assert len(result) == 1
        assert result[0] == str(output_bam)

    def test_validates_model_type(self, tmp_path):
        """Test that invalid model types are rejected."""
        # Arrange
        ref_fa = tmp_path / "ref.fa"
        model_file = tmp_path / "model.model"
        ref_fa.write_text(">chr1\nACGT\n")
        model_file.write_text("data")

        # Act & Assert
        with pytest.raises(FileOperationError, match="Invalid pbsim3 model type"):
            run_pbsim3_simulation(
                pbsim3_cmd="pbsim",
                samtools_cmd="samtools",
                reference=str(ref_fa),
                model_type="invalid_type",  # Invalid!
                model_file=str(model_file),
                coverage=30.0,
                output_prefix="out",
            )

    def test_validates_model_file_exists(self, tmp_path):
        """Test that missing model files are detected."""
        # Arrange
        ref_fa = tmp_path / "ref.fa"
        ref_fa.write_text(">chr1\nACGT\n")
        missing_model = tmp_path / "nonexistent.model"

        # Act & Assert
        with pytest.raises(FileOperationError, match="model file not found"):
            run_pbsim3_simulation(
                pbsim3_cmd="pbsim",
                samtools_cmd="samtools",
                reference=str(ref_fa),
                model_type="qshmm",
                model_file=str(missing_model),  # Doesn't exist!
                coverage=30.0,
                output_prefix="out",
            )

    def test_validates_reference_exists(self, tmp_path):
        """Test that missing reference files are detected."""
        # Arrange
        model_file = tmp_path / "model.model"
        model_file.write_text("data")
        missing_ref = tmp_path / "nonexistent.fa"

        # Act & Assert
        with pytest.raises(FileOperationError, match="Reference FASTA file not found"):
            run_pbsim3_simulation(
                pbsim3_cmd="pbsim",
                samtools_cmd="samtools",
                reference=str(missing_ref),  # Doesn't exist!
                model_type="qshmm",
                model_file=str(model_file),
                coverage=30.0,
                output_prefix="out",
            )

    def test_validates_pass_num_minimum(self, tmp_path):
        """Test that pass_num < 2 is rejected."""
        # Arrange
        ref_fa = tmp_path / "ref.fa"
        model_file = tmp_path / "model.model"
        ref_fa.write_text(">chr1\nACGT\n")
        model_file.write_text("data")

        # Act & Assert
        with pytest.raises(FileOperationError, match="pass_num ≥ 2"):
            run_pbsim3_simulation(
                pbsim3_cmd="pbsim",
                samtools_cmd="samtools",
                reference=str(ref_fa),
                model_type="qshmm",
                model_file=str(model_file),
                coverage=30.0,
                output_prefix="out",
                pass_num=1,  # Too low!
            )

    def test_handles_sam_to_bam_conversion(self, mocker, tmp_path):
        """Test that SAM output is automatically converted to BAM."""
        # Arrange
        ref_fa = tmp_path / "ref.fa"
        model_file = tmp_path / "model.model"
        output_prefix = tmp_path / "output"
        output_sam = tmp_path / "output.sam"
        output_bam = tmp_path / "output.bam"

        ref_fa.write_text(">chr1\nACGT\n")
        model_file.write_text("data")

        # Mock run_command to create SAM instead of BAM
        def mock_run_cmd(cmd, **kwargs):
            if "pbsim" in str(cmd):
                output_sam.write_text("SAM data")

        mocker.patch(
            "muc_one_up.read_simulator.wrappers.pbsim3_wrapper.run_command",
            side_effect=mock_run_cmd,
        )

        # Mock convert_sam_to_bam
        mock_convert = mocker.patch(
            "muc_one_up.read_simulator.wrappers.pbsim3_wrapper.convert_sam_to_bam"
        )
        mock_convert.return_value = str(output_bam)

        # Ensure BAM file is created when convert is called
        def side_effect_create_bam(*args, **kwargs):
            output_bam.write_text("BAM data")
            return str(output_bam)

        mock_convert.side_effect = side_effect_create_bam

        # Act
        result = run_pbsim3_simulation(
            pbsim3_cmd="pbsim",
            samtools_cmd="samtools",
            reference=str(ref_fa),
            model_type="qshmm",
            model_file=str(model_file),
            coverage=30.0,
            output_prefix=str(output_prefix),
        )

        # Assert: Conversion was called
        from muc_one_up.read_simulator.wrappers.pbsim3_wrapper import convert_sam_to_bam

        convert_sam_to_bam.assert_called_once()

        # Verify BAM was created (result is now a list)
        assert isinstance(result, list)
        assert len(result) == 1
        assert Path(result[0]).exists()
        assert result[0] == str(output_bam)

    def test_validates_empty_output(self, mocker, tmp_path):
        """Test that empty output BAM is detected."""
        # Arrange
        ref_fa = tmp_path / "ref.fa"
        model_file = tmp_path / "model.model"
        output_prefix = tmp_path / "output"
        output_bam = tmp_path / "output.bam"

        ref_fa.write_text(">chr1\nACGT\n")
        model_file.write_text("data")

        # Mock run_command to create empty BAM
        def mock_run_cmd(cmd, **kwargs):
            if "pbsim" in str(cmd):
                output_bam.write_text("")  # Empty!

        mocker.patch(
            "muc_one_up.read_simulator.wrappers.pbsim3_wrapper.run_command",
            side_effect=mock_run_cmd,
        )

        # Act & Assert
        with pytest.raises(FileOperationError, match="empty BAM file"):
            run_pbsim3_simulation(
                pbsim3_cmd="pbsim",
                samtools_cmd="samtools",
                reference=str(ref_fa),
                model_type="qshmm",
                model_file=str(model_file),
                coverage=30.0,
                output_prefix=str(output_prefix),
            )


class TestValidatePbsim3Parameters:
    """Test validate_pbsim3_parameters validation logic."""

    def test_validates_coverage_range(self):
        """Test coverage must be within valid range."""
        # Valid coverage
        validate_pbsim3_parameters(
            coverage=30,
            pass_num=3,
            accuracy_mean=0.85,
            accuracy_sd=0.05,
            accuracy_min=0.75,
            length_mean=15000,
            length_min=5000,
            length_max=30000,
        )

        # Coverage too low
        with pytest.raises(ValueError, match="Invalid coverage"):
            validate_pbsim3_parameters(
                coverage=0.05,  # Below MIN_COVERAGE
                pass_num=3,
                accuracy_mean=0.85,
                accuracy_sd=0.05,
                accuracy_min=0.75,
                length_mean=15000,
                length_min=5000,
                length_max=30000,
            )

    def test_validates_accuracy_range(self):
        """Test accuracy parameters must be 0.0-1.0."""
        # accuracy_mean out of range
        with pytest.raises(ValueError, match="Invalid accuracy_mean"):
            validate_pbsim3_parameters(
                coverage=30,
                pass_num=3,
                accuracy_mean=1.5,  # > 1.0
                accuracy_sd=0.05,
                accuracy_min=0.75,
                length_mean=15000,
                length_min=5000,
                length_max=30000,
            )

    def test_validates_length_constraints(self):
        """Test length_max >= length_min and length_mean in range."""
        # length_max < length_min
        with pytest.raises(ValueError, match=r"length_max.*<.*length_min"):
            validate_pbsim3_parameters(
                coverage=30,
                pass_num=3,
                accuracy_mean=0.85,
                accuracy_sd=0.05,
                accuracy_min=0.75,
                length_mean=15000,
                length_min=20000,  # min > max!
                length_max=10000,
            )


class TestRunPbsim3TemplateSimulation:
    """Tests for run_pbsim3_template_simulation()."""

    def test_builds_templ_strategy_command(self, tmp_path):
        """Verify command uses --strategy templ and --template."""
        template_fa = tmp_path / "template.fa"
        template_fa.write_text(">copy_001\nACGT\n>copy_002\nACGT\n")
        model_file = tmp_path / "test.model"
        model_file.write_text("model data")

        output_bam = tmp_path / "out.bam"
        output_bam.write_bytes(b"BAM\x01FAKE")

        with patch("muc_one_up.read_simulator.wrappers.pbsim3_wrapper.run_command") as mock_run:
            mock_run.return_value = MagicMock(returncode=0)

            run_pbsim3_template_simulation(
                pbsim3_cmd="pbsim",
                samtools_cmd="samtools",
                template_fasta=str(template_fa),
                model_type="qshmm",
                model_file=str(model_file),
                output_prefix=str(tmp_path / "out"),
            )

            cmd = mock_run.call_args[0][0]
            assert "--strategy" in cmd
            idx = cmd.index("--strategy")
            assert cmd[idx + 1] == "templ"
            assert "--template" in cmd
            assert "--genome" not in cmd
            assert "--depth" not in cmd

    def test_rejects_invalid_model_type(self, tmp_path):
        template_fa = tmp_path / "template.fa"
        template_fa.write_text(">copy\nACGT\n")
        model_file = tmp_path / "test.model"
        model_file.write_text("model")

        with pytest.raises(FileOperationError, match="Invalid pbsim3 model type"):
            run_pbsim3_template_simulation(
                pbsim3_cmd="pbsim",
                samtools_cmd="samtools",
                template_fasta=str(template_fa),
                model_type="invalid",
                model_file=str(model_file),
                output_prefix=str(tmp_path / "out"),
            )

    def test_rejects_missing_template(self, tmp_path):
        model_file = tmp_path / "test.model"
        model_file.write_text("model")

        with pytest.raises(FileOperationError, match=r"Template FASTA.*not found"):
            run_pbsim3_template_simulation(
                pbsim3_cmd="pbsim",
                samtools_cmd="samtools",
                template_fasta=str(tmp_path / "nonexistent.fa"),
                model_type="qshmm",
                model_file=str(model_file),
                output_prefix=str(tmp_path / "out"),
            )

    def test_includes_seed_when_provided(self, tmp_path):
        template_fa = tmp_path / "template.fa"
        template_fa.write_text(">copy\nACGT\n")
        model_file = tmp_path / "test.model"
        model_file.write_text("model")
        output_bam = tmp_path / "out.bam"
        output_bam.write_bytes(b"BAM\x01FAKE")

        with patch("muc_one_up.read_simulator.wrappers.pbsim3_wrapper.run_command") as mock_run:
            mock_run.return_value = MagicMock(returncode=0)

            run_pbsim3_template_simulation(
                pbsim3_cmd="pbsim",
                samtools_cmd="samtools",
                template_fasta=str(template_fa),
                model_type="qshmm",
                model_file=str(model_file),
                output_prefix=str(tmp_path / "out"),
                seed=42,
            )

            cmd = mock_run.call_args[0][0]
            assert "--seed" in cmd
            idx = cmd.index("--seed")
            assert str(cmd[idx + 1]) == "42"

    def test_pass_num_1_accepted(self, tmp_path):
        """pass_num=1 must be accepted for ONT single-pass simulation."""
        model = tmp_path / "test.model"
        model.write_text("model")
        template = tmp_path / "template.fa"
        template.write_text(">seq\nACGT\n")

        # Should not raise FileOperationError for pass_num validation.
        # It WILL raise for missing pbsim3 executable — that's expected.
        with pytest.raises(Exception) as exc_info:
            run_pbsim3_template_simulation(
                pbsim3_cmd="nonexistent_pbsim3",
                samtools_cmd="samtools",
                template_fasta=str(template),
                model_type="errhmm",
                model_file=str(model),
                output_prefix=str(tmp_path / "out"),
                pass_num=1,
            )
        # The error should NOT be the pass_num validation message
        assert "invalid pass_num" not in str(exc_info.value).lower()
        assert "requires pass_num >=" not in str(exc_info.value).lower()

    def test_pass_num_0_rejected(self, tmp_path):
        """pass_num=0 must still be rejected."""
        model = tmp_path / "test.model"
        model.write_text("model")
        template = tmp_path / "template.fa"
        template.write_text(">seq\nACGT\n")

        with pytest.raises(FileOperationError, match="pass_num"):
            run_pbsim3_template_simulation(
                pbsim3_cmd="pbsim",
                samtools_cmd="samtools",
                template_fasta=str(template),
                model_type="errhmm",
                model_file=str(model),
                output_prefix=str(tmp_path / "out"),
                pass_num=0,
            )

    def test_converts_numbered_sam_files(self, tmp_path):
        """Numbered SAM files (e.g., prefix_0001.sam) should be converted to BAM."""
        template_fa = tmp_path / "template.fa"
        template_fa.write_text(">copy\nACGT\n")
        model_file = tmp_path / "test.model"
        model_file.write_text("model")
        prefix = str(tmp_path / "out")

        # Simulate pbsim3 producing numbered SAM files
        (tmp_path / "out_0001.sam").write_text("@HD\tVN:1.6\n")

        with (
            patch("muc_one_up.read_simulator.wrappers.pbsim3_wrapper.run_command") as mock_run,
            patch(
                "muc_one_up.read_simulator.wrappers.pbsim3_wrapper.convert_sam_to_bam"
            ) as mock_convert,
        ):
            mock_run.return_value = MagicMock(returncode=0)

            # convert_sam_to_bam should create the BAM
            def create_bam(samtools_cmd, input_sam, output_bam, threads, timeout):
                Path(output_bam).write_bytes(b"BAM\x01FAKE")
                return output_bam

            mock_convert.side_effect = create_bam

            result = run_pbsim3_template_simulation(
                pbsim3_cmd="pbsim",
                samtools_cmd="samtools",
                template_fasta=str(template_fa),
                model_type="errhmm",
                model_file=str(model_file),
                output_prefix=prefix,
            )

            # Should have converted the numbered SAM
            mock_convert.assert_called_once()
            assert "out_0001.sam" in mock_convert.call_args[1]["input_sam"]
            assert len(result) == 1
            assert result[0].endswith("out_0001.bam")
