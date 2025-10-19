"""Tests for ONT pipeline - focuses on orchestration logic.

Following Phase 3 testing principles:
- Mock ONLY at wrapper boundaries (nanosim_wrapper, not NanoSim itself)
- Test OUR orchestration logic (config validation, error handling, file paths)
- NOT testing NanoSim or minimap2 themselves

These tests verify that the pipeline correctly:
1. Validates required configuration parameters
2. Constructs correct output file paths
3. Calls wrapper functions with correct parameters
4. Handles errors appropriately
5. Returns correct output BAM path
"""

import pytest

from muc_one_up.read_simulator.ont_pipeline import simulate_ont_reads_pipeline


class TestConfigValidation:
    """Test configuration parameter validation."""

    def test_raises_on_missing_nanosim_command(self, tmp_path):
        """Test that missing nanosim command raises ValueError."""
        input_fa = tmp_path / "input.fa"
        input_fa.write_text(">chr1\nACGT\n")

        config = {
            "tools": {
                "minimap2": "minimap2",
                "samtools": "samtools",
                # Missing nanosim
            },
            "nanosim_params": {
                "training_data_path": "model",
                "coverage": 30,
            },
        }

        with pytest.raises(ValueError, match="Missing 'nanosim' command"):
            simulate_ont_reads_pipeline(config, str(input_fa))

    def test_raises_on_missing_minimap2_command(self, tmp_path):
        """Test that missing minimap2 command raises ValueError."""
        input_fa = tmp_path / "input.fa"
        input_fa.write_text(">chr1\nACGT\n")

        config = {
            "tools": {
                "nanosim": "nanosim",
                "samtools": "samtools",
                # Missing minimap2
            },
            "nanosim_params": {
                "training_data_path": "model",
                "coverage": 30,
            },
        }

        with pytest.raises(ValueError, match="Missing 'minimap2' command"):
            simulate_ont_reads_pipeline(config, str(input_fa))

    def test_raises_on_missing_samtools_command(self, tmp_path):
        """Test that missing samtools command raises ValueError."""
        input_fa = tmp_path / "input.fa"
        input_fa.write_text(">chr1\nACGT\n")

        config = {
            "tools": {
                "nanosim": "nanosim",
                "minimap2": "minimap2",
                # Missing samtools
            },
            "nanosim_params": {
                "training_data_path": "model",
                "coverage": 30,
            },
        }

        with pytest.raises(ValueError, match="Missing 'samtools' command"):
            simulate_ont_reads_pipeline(config, str(input_fa))

    def test_raises_on_missing_training_model(self, tmp_path):
        """Test that missing training_data_path raises ValueError."""
        input_fa = tmp_path / "input.fa"
        input_fa.write_text(">chr1\nACGT\n")

        config = {
            "tools": {
                "nanosim": "nanosim",
                "minimap2": "minimap2",
                "samtools": "samtools",
            },
            "nanosim_params": {
                # Missing training_data_path
                "coverage": 30,
            },
        }

        with pytest.raises(ValueError, match="Missing 'training_data_path'"):
            simulate_ont_reads_pipeline(config, str(input_fa))

    def test_raises_on_missing_coverage(self, tmp_path):
        """Test that missing coverage raises ValueError."""
        input_fa = tmp_path / "input.fa"
        input_fa.write_text(">chr1\nACGT\n")

        config = {
            "tools": {
                "nanosim": "nanosim",
                "minimap2": "minimap2",
                "samtools": "samtools",
            },
            "nanosim_params": {
                "training_data_path": "model",
                # Missing coverage
            },
        }

        with pytest.raises(ValueError, match="Missing 'coverage'"):
            simulate_ont_reads_pipeline(config, str(input_fa))


class TestOutputPathGeneration:
    """Test output path generation logic."""

    def test_generates_correct_output_prefix(self, mocker, tmp_path):
        """Test that output prefix is constructed correctly."""
        input_fa = tmp_path / "test_sample.fa"
        input_fa.write_text(">chr1\nACGT\n")

        config = {
            "tools": {
                "nanosim": "nanosim",
                "minimap2": "minimap2",
                "samtools": "samtools",
            },
            "nanosim_params": {
                "training_data_path": str(tmp_path / "model"),
                "coverage": 30,
            },
            "read_simulation": {},
        }

        # Mock the wrapper functions
        mock_nanosim = mocker.patch("muc_one_up.read_simulator.ont_pipeline.run_nanosim_simulation")
        mock_nanosim.return_value = str(tmp_path / "test_sample_ont_reads.fastq")

        mocker.patch("muc_one_up.read_simulator.ont_pipeline.align_ont_reads_with_minimap2")

        # Act
        result = simulate_ont_reads_pipeline(config, str(input_fa))

        # Assert: Check that output prefix was used correctly
        assert mock_nanosim.call_count == 1
        call_kwargs = mock_nanosim.call_args[1]
        assert "test_sample_ont" in call_kwargs["output_prefix"]

        # Check BAM path returned
        assert result.endswith("test_sample_ont.bam")

    def test_uses_custom_output_dir(self, mocker, tmp_path):
        """Test that custom output_dir is respected."""
        input_fa = tmp_path / "input.fa"
        input_fa.write_text(">chr1\nACGT\n")

        custom_dir = tmp_path / "custom_output"

        config = {
            "tools": {
                "nanosim": "nanosim",
                "minimap2": "minimap2",
                "samtools": "samtools",
            },
            "nanosim_params": {
                "training_data_path": str(tmp_path / "model"),
                "coverage": 30,
            },
            "read_simulation": {
                "output_dir": str(custom_dir),
            },
        }

        # Mock wrapper functions
        mock_nanosim = mocker.patch("muc_one_up.read_simulator.ont_pipeline.run_nanosim_simulation")
        mock_nanosim.return_value = str(custom_dir / "input_ont_reads.fastq")

        mocker.patch("muc_one_up.read_simulator.ont_pipeline.align_ont_reads_with_minimap2")

        # Act
        result = simulate_ont_reads_pipeline(config, str(input_fa))

        # Assert: Output directory was created and used
        assert custom_dir.exists()
        assert str(custom_dir) in result


class TestWrapperIntegration:
    """Test integration with wrapper functions."""

    def test_calls_nanosim_with_correct_parameters(self, mocker, tmp_path):
        """Test that NanoSim wrapper is called with correct parameters."""
        input_fa = tmp_path / "input.fa"
        input_fa.write_text(">chr1\nACGT\n")

        config = {
            "tools": {
                "nanosim": "/path/to/nanosim",
                "minimap2": "minimap2",
                "samtools": "samtools",
            },
            "nanosim_params": {
                "training_data_path": "/path/to/model",
                "coverage": 50,
                "num_threads": 8,
                "min_read_length": 1000,
                "max_read_length": 50000,
                "other_options": "--dna_type linear",
            },
        }

        # Mock wrapper functions
        mock_nanosim = mocker.patch("muc_one_up.read_simulator.ont_pipeline.run_nanosim_simulation")
        mock_nanosim.return_value = str(tmp_path / "input_ont_reads.fastq")

        mocker.patch("muc_one_up.read_simulator.ont_pipeline.align_ont_reads_with_minimap2")

        # Act
        simulate_ont_reads_pipeline(config, str(input_fa))

        # Assert: Verify NanoSim was called with correct parameters
        assert mock_nanosim.call_count == 1
        call_kwargs = mock_nanosim.call_args[1]

        assert call_kwargs["nanosim_cmd"] == "/path/to/nanosim"
        assert call_kwargs["reference_fasta"] == str(input_fa)
        assert call_kwargs["training_model"] == "/path/to/model"
        assert call_kwargs["coverage"] == 50.0
        assert call_kwargs["threads"] == 8
        assert call_kwargs["min_read_length"] == 1000
        assert call_kwargs["max_read_length"] == 50000
        assert call_kwargs["other_options"] == "--dna_type linear"

    def test_calls_minimap2_with_correct_parameters(self, mocker, tmp_path):
        """Test that minimap2 alignment wrapper is called correctly."""
        input_fa = tmp_path / "input.fa"
        human_ref = tmp_path / "human.fa"
        input_fa.write_text(">chr1\nACGT\n")
        human_ref.write_text(">chr1\nACGT\n")

        config = {
            "tools": {
                "nanosim": "nanosim",
                "minimap2": "/path/to/minimap2",
                "samtools": "/path/to/samtools",
            },
            "nanosim_params": {
                "training_data_path": "model",
                "coverage": 30,
                "num_threads": 16,
            },
        }

        # Mock wrapper functions
        mock_nanosim = mocker.patch("muc_one_up.read_simulator.ont_pipeline.run_nanosim_simulation")
        fastq_output = str(tmp_path / "input_ont_reads.fastq")
        mock_nanosim.return_value = fastq_output

        mock_align = mocker.patch(
            "muc_one_up.read_simulator.ont_pipeline.align_ont_reads_with_minimap2"
        )

        # Act: Provide human_reference
        simulate_ont_reads_pipeline(config, str(input_fa), human_reference=str(human_ref))

        # Assert: Verify alignment was called with correct parameters
        assert mock_align.call_count == 1
        call_kwargs = mock_align.call_args[1]

        assert call_kwargs["minimap2_cmd"] == "/path/to/minimap2"
        assert call_kwargs["samtools_cmd"] == "/path/to/samtools"
        assert call_kwargs["human_reference"] == str(human_ref)
        assert call_kwargs["reads_fastq"] == fastq_output
        assert call_kwargs["threads"] == 16
        assert call_kwargs["output_bam"].endswith("input_ont.bam")

    def test_uses_input_fa_as_reference_when_human_ref_not_provided(self, mocker, tmp_path):
        """Test that input_fa is used for alignment if human_reference is None."""
        input_fa = tmp_path / "input.fa"
        input_fa.write_text(">chr1\nACGT\n")

        config = {
            "tools": {
                "nanosim": "nanosim",
                "minimap2": "minimap2",
                "samtools": "samtools",
            },
            "nanosim_params": {
                "training_data_path": "model",
                "coverage": 30,
            },
        }

        # Mock wrapper functions
        mock_nanosim = mocker.patch("muc_one_up.read_simulator.ont_pipeline.run_nanosim_simulation")
        mock_nanosim.return_value = str(tmp_path / "reads.fastq")

        mock_align = mocker.patch(
            "muc_one_up.read_simulator.ont_pipeline.align_ont_reads_with_minimap2"
        )

        # Act: Don't provide human_reference
        simulate_ont_reads_pipeline(config, str(input_fa))

        # Assert: input_fa should be used as reference
        call_kwargs = mock_align.call_args[1]
        assert call_kwargs["human_reference"] == str(input_fa)


class TestDefaultParameters:
    """Test default parameter handling."""

    def test_uses_default_threads_when_not_specified(self, mocker, tmp_path):
        """Test that threads defaults to 4 when not specified."""
        input_fa = tmp_path / "input.fa"
        input_fa.write_text(">chr1\nACGT\n")

        config = {
            "tools": {
                "nanosim": "nanosim",
                "minimap2": "minimap2",
                "samtools": "samtools",
            },
            "nanosim_params": {
                "training_data_path": "model",
                "coverage": 30,
                # No num_threads specified
            },
        }

        # Mock wrapper functions
        mock_nanosim = mocker.patch("muc_one_up.read_simulator.ont_pipeline.run_nanosim_simulation")
        mock_nanosim.return_value = str(tmp_path / "reads.fastq")

        mock_align = mocker.patch(
            "muc_one_up.read_simulator.ont_pipeline.align_ont_reads_with_minimap2"
        )

        # Act
        simulate_ont_reads_pipeline(config, str(input_fa))

        # Assert: Default threads = 4
        nanosim_kwargs = mock_nanosim.call_args[1]
        assert nanosim_kwargs["threads"] == 4

        align_kwargs = mock_align.call_args[1]
        assert align_kwargs["threads"] == 4

    def test_prioritizes_nanosim_threads_over_read_simulation_threads(self, mocker, tmp_path):
        """Test that nanosim_params.num_threads takes precedence."""
        input_fa = tmp_path / "input.fa"
        input_fa.write_text(">chr1\nACGT\n")

        config = {
            "tools": {
                "nanosim": "nanosim",
                "minimap2": "minimap2",
                "samtools": "samtools",
            },
            "nanosim_params": {
                "training_data_path": "model",
                "coverage": 30,
                "num_threads": 12,  # Should take precedence
            },
            "read_simulation": {
                "threads": 8,  # Should be ignored
            },
        }

        # Mock wrapper functions
        mock_nanosim = mocker.patch("muc_one_up.read_simulator.ont_pipeline.run_nanosim_simulation")
        mock_nanosim.return_value = str(tmp_path / "reads.fastq")

        mocker.patch("muc_one_up.read_simulator.ont_pipeline.align_ont_reads_with_minimap2")

        # Act
        simulate_ont_reads_pipeline(config, str(input_fa))

        # Assert: nanosim_params.num_threads (12) should be used
        nanosim_kwargs = mock_nanosim.call_args[1]
        assert nanosim_kwargs["threads"] == 12

    def test_uses_empty_string_for_missing_other_options(self, mocker, tmp_path):
        """Test that other_options defaults to empty string."""
        input_fa = tmp_path / "input.fa"
        input_fa.write_text(">chr1\nACGT\n")

        config = {
            "tools": {
                "nanosim": "nanosim",
                "minimap2": "minimap2",
                "samtools": "samtools",
            },
            "nanosim_params": {
                "training_data_path": "model",
                "coverage": 30,
                # No other_options
            },
        }

        # Mock wrapper functions
        mock_nanosim = mocker.patch("muc_one_up.read_simulator.ont_pipeline.run_nanosim_simulation")
        mock_nanosim.return_value = str(tmp_path / "reads.fastq")

        mocker.patch("muc_one_up.read_simulator.ont_pipeline.align_ont_reads_with_minimap2")

        # Act
        simulate_ont_reads_pipeline(config, str(input_fa))

        # Assert: other_options should be empty string
        nanosim_kwargs = mock_nanosim.call_args[1]
        assert nanosim_kwargs["other_options"] == ""


class TestErrorHandling:
    """Test error handling and propagation."""

    def test_raises_runtime_error_when_nanosim_fails(self, mocker, tmp_path):
        """Test that NanoSim failure raises RuntimeError."""
        input_fa = tmp_path / "input.fa"
        input_fa.write_text(">chr1\nACGT\n")

        config = {
            "tools": {
                "nanosim": "nanosim",
                "minimap2": "minimap2",
                "samtools": "samtools",
            },
            "nanosim_params": {
                "training_data_path": "model",
                "coverage": 30,
            },
        }

        # Mock NanoSim to raise exception
        mock_nanosim = mocker.patch("muc_one_up.read_simulator.ont_pipeline.run_nanosim_simulation")
        mock_nanosim.side_effect = Exception("NanoSim failed")

        # Act & Assert
        with pytest.raises(RuntimeError, match="ONT read simulation failed"):
            simulate_ont_reads_pipeline(config, str(input_fa))

    def test_raises_runtime_error_when_alignment_fails(self, mocker, tmp_path):
        """Test that alignment failure raises RuntimeError."""
        input_fa = tmp_path / "input.fa"
        input_fa.write_text(">chr1\nACGT\n")

        config = {
            "tools": {
                "nanosim": "nanosim",
                "minimap2": "minimap2",
                "samtools": "samtools",
            },
            "nanosim_params": {
                "training_data_path": "model",
                "coverage": 30,
            },
        }

        # Mock NanoSim to succeed
        mock_nanosim = mocker.patch("muc_one_up.read_simulator.ont_pipeline.run_nanosim_simulation")
        mock_nanosim.return_value = str(tmp_path / "reads.fastq")

        # Mock alignment to fail
        mock_align = mocker.patch(
            "muc_one_up.read_simulator.ont_pipeline.align_ont_reads_with_minimap2"
        )
        mock_align.side_effect = Exception("Alignment failed")

        # Act & Assert
        with pytest.raises(RuntimeError, match="ONT read alignment failed"):
            simulate_ont_reads_pipeline(config, str(input_fa))


class TestReturnValue:
    """Test that correct output BAM path is returned."""

    def test_returns_correct_bam_path(self, mocker, tmp_path):
        """Test that the function returns the correct output BAM path."""
        input_fa = tmp_path / "sample.fa"
        input_fa.write_text(">chr1\nACGT\n")

        config = {
            "tools": {
                "nanosim": "nanosim",
                "minimap2": "minimap2",
                "samtools": "samtools",
            },
            "nanosim_params": {
                "training_data_path": "model",
                "coverage": 30,
            },
        }

        # Mock wrapper functions
        mock_nanosim = mocker.patch("muc_one_up.read_simulator.ont_pipeline.run_nanosim_simulation")
        mock_nanosim.return_value = str(tmp_path / "reads.fastq")

        mocker.patch("muc_one_up.read_simulator.ont_pipeline.align_ont_reads_with_minimap2")

        # Act
        result = simulate_ont_reads_pipeline(config, str(input_fa))

        # Assert
        assert result.endswith("sample_ont.bam")
        assert str(tmp_path) in result
