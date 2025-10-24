#!/usr/bin/env python3
"""
Tests for tool version detection.

Testing Strategy:
- Mock subprocess.run to avoid requiring actual tool installations
- Test all 9 tool parsers with verified real-world output
- Test graceful degradation and error handling
- Test batch operations and logging
"""

import logging
import subprocess
from unittest.mock import MagicMock, patch

from muc_one_up.read_simulator.utils.tool_version import (
    capture_tool_versions,
    get_tool_version,
    log_tool_versions,
)


class TestToolVersionParsers:
    """Test individual tool parsers with verified output."""

    @patch("muc_one_up.read_simulator.utils.tool_version.subprocess.run")
    @patch("muc_one_up.read_simulator.utils.tool_version.build_tool_command")
    def test_samtools_version(self, mock_build_cmd, mock_run):
        """Test samtools parser with real output."""
        mock_build_cmd.return_value = ["samtools", "--version"]
        mock_run.return_value = MagicMock(
            stdout="samtools 1.17\nUsing htslib 1.17\n", stderr="", returncode=0
        )

        version = get_tool_version("samtools", "samtools")
        assert version == "samtools 1.17"

    @patch("muc_one_up.read_simulator.utils.tool_version.subprocess.run")
    @patch("muc_one_up.read_simulator.utils.tool_version.build_tool_command")
    def test_minimap2_version(self, mock_build_cmd, mock_run):
        """Test minimap2 parser (version only, prepend tool name)."""
        mock_build_cmd.return_value = ["minimap2", "--version"]
        mock_run.return_value = MagicMock(stdout="2.28-r1209\n", stderr="", returncode=0)

        version = get_tool_version("minimap2", "minimap2")
        assert version == "minimap2 2.28-r1209"

    @patch("muc_one_up.read_simulator.utils.tool_version.subprocess.run")
    @patch("muc_one_up.read_simulator.utils.tool_version.build_tool_command")
    def test_nanosim_version(self, mock_build_cmd, mock_run):
        """Test NanoSim parser (binary: simulator.py)."""
        mock_build_cmd.return_value = ["simulator.py", "--version"]
        mock_run.return_value = MagicMock(stdout="NanoSim 3.2.2\n", stderr="", returncode=0)

        version = get_tool_version("simulator.py", "nanosim")
        assert version == "NanoSim 3.2.2"

    @patch("muc_one_up.read_simulator.utils.tool_version.subprocess.run")
    @patch("muc_one_up.read_simulator.utils.tool_version.build_tool_command")
    def test_ccs_version(self, mock_build_cmd, mock_run):
        """Test ccs parser."""
        mock_build_cmd.return_value = ["ccs", "--version"]
        mock_run.return_value = MagicMock(
            stdout="ccs 6.4.0 (commit v6.4.0)\n", stderr="", returncode=0
        )

        version = get_tool_version("ccs", "ccs")
        assert version == "ccs 6.4.0 (commit v6.4.0)"

    @patch("muc_one_up.read_simulator.utils.tool_version.subprocess.run")
    @patch("muc_one_up.read_simulator.utils.tool_version.build_tool_command")
    def test_reseq_version(self, mock_build_cmd, mock_run):
        """Test reseq parser."""
        mock_build_cmd.return_value = ["reseq", "--version"]
        mock_run.return_value = MagicMock(stdout="ReSeq version 1.1\n", stderr="", returncode=0)

        version = get_tool_version("reseq", "reseq")
        assert version == "ReSeq version 1.1"

    @patch("muc_one_up.read_simulator.utils.tool_version.subprocess.run")
    @patch("muc_one_up.read_simulator.utils.tool_version.build_tool_command")
    def test_bwa_version_from_stderr(self, mock_build_cmd, mock_run):
        """Test bwa parser (stderr, exit code 1)."""
        mock_build_cmd.return_value = ["bwa"]
        mock_run.return_value = MagicMock(
            stdout="",
            stderr="Program: bwa (alignment via Burrows-Wheeler transformation)\nVersion: 0.7.18-r1243-dirty\n",
            returncode=1,
        )

        version = get_tool_version("bwa", "bwa")
        assert version == "bwa 0.7.18-r1243-dirty"

    @patch("muc_one_up.read_simulator.utils.tool_version.subprocess.run")
    @patch("muc_one_up.read_simulator.utils.tool_version.build_tool_command")
    def test_pblat_version_from_usage(self, mock_build_cmd, mock_run):
        """Test pblat parser (parse from usage line)."""
        mock_build_cmd.return_value = ["pblat"]
        mock_run.return_value = MagicMock(
            stdout="pblat - BLAT with parallel supports v. 36x2 fast sequence search\n",
            stderr="",
            returncode=0,
        )

        version = get_tool_version("pblat", "pblat")
        assert version == "pblat v.36x2"

    @patch("muc_one_up.read_simulator.utils.tool_version.subprocess.run")
    @patch("muc_one_up.read_simulator.utils.tool_version.build_tool_command")
    def test_pbsim3_no_version(self, mock_build_cmd, mock_run):
        """Test pbsim3 parser (NO version in output, return generic string)."""
        mock_build_cmd.return_value = ["pbsim"]
        mock_run.return_value = MagicMock(
            stdout="USAGE: pbsim [options]\n", stderr="", returncode=255
        )

        version = get_tool_version("pbsim", "pbsim3")
        assert version == "pbsim3 (installed)"

    @patch("muc_one_up.read_simulator.utils.tool_version.subprocess.run")
    @patch("muc_one_up.read_simulator.utils.tool_version.build_tool_command")
    def test_fatotwobit_no_version(self, mock_build_cmd, mock_run):
        """Test faToTwoBit parser (NO version info anywhere)."""
        mock_build_cmd.return_value = ["faToTwoBit"]
        mock_run.return_value = MagicMock(
            stdout="faToTwoBit - Convert DNA from fasta to 2bit format\n",
            stderr="",
            returncode=255,
        )

        version = get_tool_version("faToTwoBit", "faToTwoBit")
        assert version == "faToTwoBit (version unknown)"


class TestErrorHandling:
    """Test graceful degradation and error handling."""

    @patch("muc_one_up.read_simulator.utils.tool_version.subprocess.run")
    @patch("muc_one_up.read_simulator.utils.tool_version.build_tool_command")
    def test_unknown_tool(self, mock_build_cmd, mock_run):
        """Test handling of unknown tool (not in tool_map)."""
        version = get_tool_version("unknown_tool", "unknown_tool")
        assert version == "N/A"

    @patch("muc_one_up.read_simulator.utils.tool_version.subprocess.run")
    @patch("muc_one_up.read_simulator.utils.tool_version.build_tool_command")
    def test_timeout_handling(self, mock_build_cmd, mock_run):
        """Test handling of subprocess timeout."""
        mock_build_cmd.return_value = ["samtools", "--version"]
        mock_run.side_effect = subprocess.TimeoutExpired(cmd=["samtools"], timeout=5)

        version = get_tool_version("samtools", "samtools")
        assert version == "N/A"

    @patch("muc_one_up.read_simulator.utils.tool_version.subprocess.run")
    @patch("muc_one_up.read_simulator.utils.tool_version.build_tool_command")
    def test_subprocess_exception(self, mock_build_cmd, mock_run):
        """Test handling of generic subprocess exception."""
        mock_build_cmd.return_value = ["samtools", "--version"]
        mock_run.side_effect = OSError("Tool not found")

        version = get_tool_version("samtools", "samtools")
        assert version == "N/A"

    @patch("muc_one_up.read_simulator.utils.tool_version.subprocess.run")
    @patch("muc_one_up.read_simulator.utils.tool_version.build_tool_command")
    def test_unparseable_output(self, mock_build_cmd, mock_run):
        """Test handling when parser can't extract version."""
        mock_build_cmd.return_value = ["samtools", "--version"]
        mock_run.return_value = MagicMock(
            stdout="Unexpected output format\n", stderr="", returncode=0
        )

        version = get_tool_version("samtools", "samtools")
        assert version == "N/A"


class TestBatchOperations:
    """Test batch version capture."""

    @patch("muc_one_up.read_simulator.utils.tool_version.get_tool_version")
    def test_capture_tool_versions(self, mock_get_version):
        """Test capturing versions for multiple tools."""
        mock_get_version.side_effect = lambda cmd, name: f"{name} 1.0"

        tools_config = {
            "samtools": "samtools",
            "minimap2": "minimap2",
            "bwa": "bwa",
        }

        versions = capture_tool_versions(tools_config)

        assert len(versions) == 3
        assert versions["samtools"] == "samtools 1.0"
        assert versions["minimap2"] == "minimap2 1.0"
        assert versions["bwa"] == "bwa 1.0"

    @patch("muc_one_up.read_simulator.utils.tool_version.get_tool_version")
    def test_capture_with_conda_commands(self, mock_get_version):
        """Test capturing versions with conda/mamba commands."""
        mock_get_version.side_effect = lambda cmd, name: f"{name} 1.0"

        tools_config = {
            "samtools": "samtools",
            "reseq": "mamba run -n env_wessim reseq",
            "nanosim": "conda run -n env_nanosim simulator.py",
        }

        versions = capture_tool_versions(tools_config)

        assert len(versions) == 3
        # Verify all tools called
        assert mock_get_version.call_count == 3


class TestLogging:
    """Test version logging functionality."""

    def test_log_tool_versions(self, caplog):
        """Test formatted logging output."""
        with caplog.at_level(logging.INFO):
            versions = {
                "samtools": "samtools 1.17",
                "minimap2": "minimap2 2.28-r1209",
                "bwa": "bwa 0.7.18-r1243",
            }
            log_tool_versions(versions)

        # Check that logging occurred
        assert "Tool Versions:" in caplog.text
        assert "samtools 1.17" in caplog.text
        assert "minimap2 2.28-r1209" in caplog.text
        assert "bwa 0.7.18-r1243" in caplog.text

    def test_log_empty_versions(self, caplog):
        """Test logging with empty versions dict."""
        with caplog.at_level(logging.INFO):
            log_tool_versions({})

        # Should not log anything for empty dict
        assert "Tool Versions:" not in caplog.text


class TestCondaMambaSupport:
    """Test integration with conda/mamba command building."""

    @patch("muc_one_up.read_simulator.utils.tool_version.subprocess.run")
    @patch("muc_one_up.read_simulator.utils.tool_version.build_tool_command")
    def test_conda_run_command(self, mock_build_cmd, mock_run):
        """Test conda run command is built correctly."""
        mock_build_cmd.return_value = [
            "conda",
            "run",
            "-n",
            "env_wessim",
            "reseq",
            "--version",
        ]
        mock_run.return_value = MagicMock(stdout="ReSeq version 1.1\n", stderr="", returncode=0)

        version = get_tool_version("conda run -n env_wessim reseq", "reseq")

        assert version == "ReSeq version 1.1"
        # Verify build_tool_command was called with correct args
        mock_build_cmd.assert_called_once_with("conda run -n env_wessim reseq", "--version")

    @patch("muc_one_up.read_simulator.utils.tool_version.subprocess.run")
    @patch("muc_one_up.read_simulator.utils.tool_version.build_tool_command")
    def test_mamba_run_command(self, mock_build_cmd, mock_run):
        """Test mamba run command is built correctly."""
        mock_build_cmd.return_value = [
            "mamba",
            "run",
            "-n",
            "env_nanosim",
            "simulator.py",
            "--version",
        ]
        mock_run.return_value = MagicMock(stdout="NanoSim 3.2.2\n", stderr="", returncode=0)

        version = get_tool_version("mamba run -n env_nanosim simulator.py", "nanosim")

        assert version == "NanoSim 3.2.2"
        mock_build_cmd.assert_called_once_with("mamba run -n env_nanosim simulator.py", "--version")


class TestNonZeroExitCodes:
    """Test handling of tools with non-zero exit codes."""

    @patch("muc_one_up.read_simulator.utils.tool_version.subprocess.run")
    @patch("muc_one_up.read_simulator.utils.tool_version.build_tool_command")
    def test_bwa_exit_code_1(self, mock_build_cmd, mock_run):
        """Test bwa with exit code 1 (check=False must be used)."""
        mock_build_cmd.return_value = ["bwa"]
        mock_run.return_value = MagicMock(
            stdout="",
            stderr="Version: 0.7.18-r1243\n",
            returncode=1,
        )

        version = get_tool_version("bwa", "bwa")

        assert version == "bwa 0.7.18-r1243"
        # Verify check=False was used
        _args, kwargs = mock_run.call_args
        assert kwargs.get("check") is False

    @patch("muc_one_up.read_simulator.utils.tool_version.subprocess.run")
    @patch("muc_one_up.read_simulator.utils.tool_version.build_tool_command")
    def test_pbsim_exit_code_255(self, mock_build_cmd, mock_run):
        """Test pbsim with exit code 255."""
        mock_build_cmd.return_value = ["pbsim"]
        mock_run.return_value = MagicMock(
            stdout="USAGE: pbsim [options]\n",
            stderr="",
            returncode=255,
        )

        version = get_tool_version("pbsim", "pbsim3")

        assert version == "pbsim3 (installed)"
        # Verify check=False was used
        _args, kwargs = mock_run.call_args
        assert kwargs.get("check") is False
