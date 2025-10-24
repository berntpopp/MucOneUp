#!/usr/bin/env python3
"""
Tests for metadata file writer.

Testing Strategy:
- Mock tool version capture to avoid requiring real tools
- Test TSV format output
- Test all three platforms (Illumina, ONT, PacBio)
- Test parameter extraction from config
- Verify file creation and content
"""

from datetime import datetime
from pathlib import Path
from unittest.mock import patch

from muc_one_up.read_simulator.utils.metadata_writer import write_metadata_file


class TestMetadataFileCreation:
    """Test basic metadata file creation."""

    @patch("muc_one_up.read_simulator.utils.metadata_writer.capture_tool_versions")
    @patch("muc_one_up.read_simulator.utils.metadata_writer.log_tool_versions")
    def test_creates_tsv_file(self, mock_log_versions, mock_capture_versions, tmp_path):
        """Test that metadata file is created with correct name."""
        mock_capture_versions.return_value = {}

        config = {"tools": {}}
        start_time = datetime(2025, 10, 24, 12, 0, 0)
        end_time = datetime(2025, 10, 24, 13, 0, 0)

        result_path = write_metadata_file(
            output_dir=str(tmp_path),
            output_base="test_sample",
            config=config,
            start_time=start_time,
            end_time=end_time,
            platform="Illumina",
        )

        assert Path(result_path).exists()
        assert result_path == str(tmp_path / "test_sample_metadata.tsv")

    @patch("muc_one_up.read_simulator.utils.metadata_writer.capture_tool_versions")
    @patch("muc_one_up.read_simulator.utils.metadata_writer.log_tool_versions")
    def test_tsv_format_structure(self, mock_log_versions, mock_capture_versions, tmp_path):
        """Test TSV file has correct header and format."""
        mock_capture_versions.return_value = {}

        config = {"tools": {}}
        start_time = datetime(2025, 10, 24, 12, 0, 0)
        end_time = datetime(2025, 10, 24, 13, 0, 0)

        result_path = write_metadata_file(
            output_dir=str(tmp_path),
            output_base="test",
            config=config,
            start_time=start_time,
            end_time=end_time,
            platform="Illumina",
        )

        with open(result_path) as f:
            lines = f.readlines()

        # Check TSV header
        assert lines[0].strip() == "Parameter\tValue"
        # Check tab-separated format
        for line in lines[1:]:
            if line.strip():  # Skip empty lines
                assert "\t" in line


class TestCoreMetadata:
    """Test core metadata fields."""

    @patch("muc_one_up.read_simulator.utils.metadata_writer.capture_tool_versions")
    @patch("muc_one_up.read_simulator.utils.metadata_writer.log_tool_versions")
    @patch("muc_one_up.read_simulator.utils.metadata_writer.__version__", "0.15.0")
    def test_core_metadata_fields(self, mock_log_versions, mock_capture_versions, tmp_path):
        """Test core metadata fields are written."""
        mock_capture_versions.return_value = {}

        config = {"tools": {}}
        start_time = datetime(2025, 10, 24, 12, 30, 45)
        end_time = datetime(2025, 10, 24, 13, 40, 55)

        result_path = write_metadata_file(
            output_dir=str(tmp_path),
            output_base="sample_001",
            config=config,
            start_time=start_time,
            end_time=end_time,
            platform="Illumina",
        )

        with open(result_path) as f:
            content = f.read()

        # Check core fields
        assert "Tool\tMucOneUp" in content
        assert "Version\t0.15.0" in content
        assert "Read_simulation_technology\tIllumina" in content
        assert "Run_start_time\t2025-10-24T12:30:45" in content
        assert "Run_end_time\t2025-10-24T13:40:55" in content
        assert "Run_duration_seconds\t4210.0" in content  # 70 minutes 10 seconds
        assert "Output_base\tsample_001" in content
        assert f"Output_dir\t{tmp_path}" in content

    @patch("muc_one_up.read_simulator.utils.metadata_writer.capture_tool_versions")
    @patch("muc_one_up.read_simulator.utils.metadata_writer.log_tool_versions")
    def test_seed_parameter(self, mock_log_versions, mock_capture_versions, tmp_path):
        """Test seed parameter is written when present."""
        mock_capture_versions.return_value = {}

        config = {"tools": {}, "seed": 42}
        start_time = datetime(2025, 10, 24, 12, 0, 0)
        end_time = datetime(2025, 10, 24, 13, 0, 0)

        result_path = write_metadata_file(
            output_dir=str(tmp_path),
            output_base="test",
            config=config,
            start_time=start_time,
            end_time=end_time,
            platform="Illumina",
        )

        with open(result_path) as f:
            content = f.read()

        assert "Seed\t42" in content


class TestPlatformSpecificMetadata:
    """Test platform-specific metadata fields."""

    @patch("muc_one_up.read_simulator.utils.metadata_writer.capture_tool_versions")
    @patch("muc_one_up.read_simulator.utils.metadata_writer.log_tool_versions")
    def test_illumina_parameters(self, mock_log_versions, mock_capture_versions, tmp_path):
        """Test Illumina-specific parameters."""
        mock_capture_versions.return_value = {}

        config = {
            "tools": {},
            "read_simulation": {
                "coverage": 30,
                "fragment_size_mean": 350,
            },
        }
        start_time = datetime(2025, 10, 24, 12, 0, 0)
        end_time = datetime(2025, 10, 24, 13, 0, 0)

        result_path = write_metadata_file(
            output_dir=str(tmp_path),
            output_base="test",
            config=config,
            start_time=start_time,
            end_time=end_time,
            platform="Illumina",
        )

        with open(result_path) as f:
            content = f.read()

        assert "Coverage\t30" in content
        assert "Fragment_size\t350" in content

    @patch("muc_one_up.read_simulator.utils.metadata_writer.capture_tool_versions")
    @patch("muc_one_up.read_simulator.utils.metadata_writer.log_tool_versions")
    def test_ont_parameters(self, mock_log_versions, mock_capture_versions, tmp_path):
        """Test ONT-specific parameters."""
        mock_capture_versions.return_value = {}

        config = {
            "tools": {},
            "nanosim_params": {
                "coverage": 20,
                "min_read_length": 1500,
                "max_read_length": 5000,
            },
        }
        start_time = datetime(2025, 10, 24, 12, 0, 0)
        end_time = datetime(2025, 10, 24, 13, 0, 0)

        result_path = write_metadata_file(
            output_dir=str(tmp_path),
            output_base="test",
            config=config,
            start_time=start_time,
            end_time=end_time,
            platform="ONT",
        )

        with open(result_path) as f:
            content = f.read()

        assert "Coverage\t20" in content
        assert "Min_read_length\t1500" in content
        assert "Max_read_length\t5000" in content

    @patch("muc_one_up.read_simulator.utils.metadata_writer.capture_tool_versions")
    @patch("muc_one_up.read_simulator.utils.metadata_writer.log_tool_versions")
    def test_pacbio_parameters(self, mock_log_versions, mock_capture_versions, tmp_path):
        """Test PacBio-specific parameters."""
        mock_capture_versions.return_value = {}

        config = {
            "tools": {},
            "pacbio_params": {
                "coverage": 25,
                "pass_num": 10,
            },
        }
        start_time = datetime(2025, 10, 24, 12, 0, 0)
        end_time = datetime(2025, 10, 24, 13, 0, 0)

        result_path = write_metadata_file(
            output_dir=str(tmp_path),
            output_base="test",
            config=config,
            start_time=start_time,
            end_time=end_time,
            platform="PacBio",
        )

        with open(result_path) as f:
            content = f.read()

        assert "Coverage\t25" in content
        assert "Pass_num\t10" in content


class TestToolVersions:
    """Test tool version integration."""

    @patch("muc_one_up.read_simulator.utils.metadata_writer.capture_tool_versions")
    @patch("muc_one_up.read_simulator.utils.metadata_writer.log_tool_versions")
    def test_tool_versions_written(self, mock_log_versions, mock_capture_versions, tmp_path):
        """Test tool versions are written with correct prefix."""
        mock_capture_versions.return_value = {
            "samtools": "samtools 1.17",
            "minimap2": "minimap2 2.28-r1209",
            "bwa": "bwa 0.7.18-r1243",
        }

        config = {
            "tools": {
                "samtools": "samtools",
                "minimap2": "minimap2",
                "bwa": "bwa",
            }
        }
        start_time = datetime(2025, 10, 24, 12, 0, 0)
        end_time = datetime(2025, 10, 24, 13, 0, 0)

        result_path = write_metadata_file(
            output_dir=str(tmp_path),
            output_base="test",
            config=config,
            start_time=start_time,
            end_time=end_time,
            platform="Illumina",
        )

        with open(result_path) as f:
            content = f.read()

        # Check tool versions with "tool." prefix (sorted alphabetically)
        assert "tool.bwa_version\tbwa 0.7.18-r1243" in content
        assert "tool.minimap2_version\tminimap2 2.28-r1209" in content
        assert "tool.samtools_version\tsamtools 1.17" in content

    @patch("muc_one_up.read_simulator.utils.metadata_writer.capture_tool_versions")
    @patch("muc_one_up.read_simulator.utils.metadata_writer.log_tool_versions")
    def test_tool_version_capture_called(self, mock_log_versions, mock_capture_versions, tmp_path):
        """Test that tool version capture is called."""
        mock_capture_versions.return_value = {}

        config = {
            "tools": {
                "samtools": "samtools",
                "bwa": "bwa",
            }
        }
        start_time = datetime(2025, 10, 24, 12, 0, 0)
        end_time = datetime(2025, 10, 24, 13, 0, 0)

        write_metadata_file(
            output_dir=str(tmp_path),
            output_base="test",
            config=config,
            start_time=start_time,
            end_time=end_time,
            platform="Illumina",
        )

        # Verify capture_tool_versions was called with tools config
        mock_capture_versions.assert_called_once_with(config["tools"])
        # Verify log_tool_versions was called with captured versions
        mock_log_versions.assert_called_once_with({})

    @patch("muc_one_up.read_simulator.utils.metadata_writer.capture_tool_versions")
    @patch("muc_one_up.read_simulator.utils.metadata_writer.log_tool_versions")
    def test_tool_versions_sorted(self, mock_log_versions, mock_capture_versions, tmp_path):
        """Test tool versions are written in sorted order."""
        mock_capture_versions.return_value = {
            "samtools": "samtools 1.17",
            "bwa": "bwa 0.7.18",
            "minimap2": "minimap2 2.28",
        }

        config = {"tools": {"samtools": "samtools", "bwa": "bwa", "minimap2": "minimap2"}}
        start_time = datetime(2025, 10, 24, 12, 0, 0)
        end_time = datetime(2025, 10, 24, 13, 0, 0)

        result_path = write_metadata_file(
            output_dir=str(tmp_path),
            output_base="test",
            config=config,
            start_time=start_time,
            end_time=end_time,
            platform="Illumina",
        )

        with open(result_path) as f:
            lines = f.readlines()

        # Find tool version lines
        tool_lines = [line for line in lines if line.startswith("tool.")]

        # Verify alphabetical order
        assert tool_lines[0].startswith("tool.bwa_version")
        assert tool_lines[1].startswith("tool.minimap2_version")
        assert tool_lines[2].startswith("tool.samtools_version")


class TestEdgeCases:
    """Test edge cases and error handling."""

    @patch("muc_one_up.read_simulator.utils.metadata_writer.capture_tool_versions")
    @patch("muc_one_up.read_simulator.utils.metadata_writer.log_tool_versions")
    def test_missing_tools_config(self, mock_log_versions, mock_capture_versions, tmp_path):
        """Test handling when tools config is missing."""
        mock_capture_versions.return_value = {}

        config = {}  # No "tools" key
        start_time = datetime(2025, 10, 24, 12, 0, 0)
        end_time = datetime(2025, 10, 24, 13, 0, 0)

        result_path = write_metadata_file(
            output_dir=str(tmp_path),
            output_base="test",
            config=config,
            start_time=start_time,
            end_time=end_time,
            platform="Illumina",
        )

        # Should still create file
        assert Path(result_path).exists()
        # capture_tool_versions should be called with empty dict
        mock_capture_versions.assert_called_once_with({})

    @patch("muc_one_up.read_simulator.utils.metadata_writer.capture_tool_versions")
    @patch("muc_one_up.read_simulator.utils.metadata_writer.log_tool_versions")
    def test_missing_platform_params(self, mock_log_versions, mock_capture_versions, tmp_path):
        """Test handling when platform-specific params are missing."""
        mock_capture_versions.return_value = {}

        config = {"tools": {}}  # No read_simulation config
        start_time = datetime(2025, 10, 24, 12, 0, 0)
        end_time = datetime(2025, 10, 24, 13, 0, 0)

        result_path = write_metadata_file(
            output_dir=str(tmp_path),
            output_base="test",
            config=config,
            start_time=start_time,
            end_time=end_time,
            platform="Illumina",
        )

        with open(result_path) as f:
            content = f.read()

        # Should use "N/A" for missing values
        assert "Coverage\tN/A" in content
        assert "Fragment_size\tN/A" in content
