"""Tests for command building utilities.

These tests verify the centralized command construction logic that prevents
command injection vulnerabilities (CWE-78, Bandit B602) by safely handling
both simple commands and complex multi-word commands (conda/mamba).

Following Phase 2 testing principles:
- Test edge cases and error conditions
- Verify security properties (no shell injection)
- Test real-world usage patterns
"""

from pathlib import Path

import pytest

from muc_one_up.read_simulator.command_utils import build_tool_command, get_tool_executable


class TestBuildToolCommand:
    """Test command building utility - core security-critical functionality."""

    def test_simple_command_no_args(self):
        """Test building simple command without additional arguments."""
        cmd = build_tool_command("bwa")
        assert cmd == ["bwa"]

    def test_simple_command_with_args(self):
        """Test building simple command with arguments."""
        cmd = build_tool_command("bwa", "mem", "-t", "4", "ref.fa", "reads.fq")
        assert cmd == ["bwa", "mem", "-t", "4", "ref.fa", "reads.fq"]

    def test_conda_command_basic(self):
        """Test building conda/mamba multi-word command (basic)."""
        cmd = build_tool_command("mamba run -n env_wessim reseq", "replaceN", "-r", "input.fa")
        assert cmd == [
            "mamba",
            "run",
            "-n",
            "env_wessim",
            "reseq",
            "replaceN",
            "-r",
            "input.fa",
        ]

    def test_conda_command_with_flags(self):
        """Test building conda/mamba command with --no-capture-output flag."""
        cmd = build_tool_command(
            "mamba run --no-capture-output -n env_wessim samtools", "view", "-b", "input.bam"
        )
        assert cmd == [
            "mamba",
            "run",
            "--no-capture-output",
            "-n",
            "env_wessim",
            "samtools",
            "view",
            "-b",
            "input.bam",
        ]

    def test_quoted_arguments_with_spaces(self):
        """Test handling of quoted arguments containing spaces."""
        cmd = build_tool_command('tool --arg "value with spaces"', "input.txt")
        assert cmd == ["tool", "--arg", "value with spaces", "input.txt"]

    def test_numeric_arguments_converted_to_strings(self):
        """Test that numeric arguments are automatically converted to strings."""
        cmd = build_tool_command("samtools", "view", "-@", 4, "-b")
        assert cmd == ["samtools", "view", "-@", "4", "-b"]

    def test_float_arguments_converted_to_strings(self):
        """Test that float arguments are converted to strings."""
        cmd = build_tool_command("tool", "--coverage", 30.5, "--threshold", 0.95)
        assert cmd == ["tool", "--coverage", "30.5", "--threshold", "0.95"]

    def test_path_arguments_converted_to_strings(self):
        """Test that Path objects are converted to strings."""
        input_path = Path("/tmp/input.fa")
        output_path = Path("/tmp/output.fa")
        cmd = build_tool_command("tool", "-i", input_path, "-o", output_path)
        assert cmd == ["tool", "-i", "/tmp/input.fa", "-o", "/tmp/output.fa"]

    def test_empty_command_raises_error(self):
        """Test that empty command raises ValueError."""
        with pytest.raises(ValueError, match="cannot be empty"):
            build_tool_command("", "arg")

    def test_none_command_raises_error(self):
        """Test that None command raises ValueError."""
        with pytest.raises(ValueError, match="cannot be empty"):
            build_tool_command(None, "arg")  # type: ignore[arg-type]

    def test_command_with_equals_sign(self):
        """Test handling of arguments with equals signs."""
        cmd = build_tool_command("tool", "--option=value", "--other=123")
        assert cmd == ["tool", "--option=value", "--other=123"]

    def test_command_with_single_quotes(self):
        """Test handling of single-quoted arguments."""
        cmd = build_tool_command("tool --arg 'single quoted'", "file.txt")
        assert cmd == ["tool", "--arg", "single quoted", "file.txt"]

    def test_command_with_escaped_quotes(self):
        """Test handling of escaped quotes in arguments."""
        cmd = build_tool_command(r'tool --arg "escaped\"quote"', "file.txt")
        assert cmd == ["tool", "--arg", 'escaped"quote', "file.txt"]

    def test_real_world_reseq_command(self):
        """Test real-world reseq command from config.json."""
        cmd = build_tool_command(
            "mamba run --no-capture-output -n env_wessim reseq",
            "replaceN",
            "-r",
            "/tmp/input.fa",
            "-R",
            "/tmp/output.fa",
        )
        assert cmd == [
            "mamba",
            "run",
            "--no-capture-output",
            "-n",
            "env_wessim",
            "reseq",
            "replaceN",
            "-r",
            "/tmp/input.fa",
            "-R",
            "/tmp/output.fa",
        ]

    def test_real_world_bwa_command(self):
        """Test real-world BWA command from config.json."""
        cmd = build_tool_command(
            "mamba run --no-capture-output -n env_wessim bwa",
            "mem",
            "-t",
            "4",
            "/tmp/ref.fa",
            "/tmp/r1.fq",
            "/tmp/r2.fq",
        )
        assert cmd == [
            "mamba",
            "run",
            "--no-capture-output",
            "-n",
            "env_wessim",
            "bwa",
            "mem",
            "-t",
            "4",
            "/tmp/ref.fa",
            "/tmp/r1.fq",
            "/tmp/r2.fq",
        ]

    def test_real_world_samtools_command(self):
        """Test real-world samtools command from config.json."""
        cmd = build_tool_command(
            "mamba run --no-capture-output -n env_wessim samtools",
            "view",
            "-@",
            "4",
            "-b",
            "-o",
            "output.bam",
            "input.sam",
        )
        assert cmd == [
            "mamba",
            "run",
            "--no-capture-output",
            "-n",
            "env_wessim",
            "samtools",
            "view",
            "-@",
            "4",
            "-b",
            "-o",
            "output.bam",
            "input.sam",
        ]

    def test_security_no_shell_injection_possible(self):
        """Test that shell injection is prevented by proper list construction."""
        # Attempt shell injection with semicolon command separator
        malicious_cmd = "tool; rm -rf /"
        cmd = build_tool_command(malicious_cmd, "arg")

        # Should parse as ["tool;", "rm", "-rf", "/", "arg"]
        # When passed to subprocess with shell=False, this will fail to execute
        # because "tool;" is not a valid executable - this is CORRECT security behavior
        assert cmd == ["tool;", "rm", "-rf", "/", "arg"]

    def test_security_no_pipe_injection_possible(self):
        """Test that pipe injection is prevented."""
        # Attempt pipe injection
        malicious_cmd = "tool | cat /etc/passwd"
        cmd = build_tool_command(malicious_cmd, "arg")

        # Should parse as separate elements, not as shell pipe
        assert "|" in cmd  # Pipe is just another argument, not a shell operator
        assert "cat" in cmd

    def test_mixed_string_and_numeric_args(self):
        """Test mixing string and numeric arguments."""
        cmd = build_tool_command("tool", "-t", 4, "--input", "file.txt", "--coverage", 30.5)
        assert cmd == ["tool", "-t", "4", "--input", "file.txt", "--coverage", "30.5"]


class TestGetToolExecutable:
    """Test tool executable extraction utility."""

    def test_simple_tool(self):
        """Test extracting name from simple command."""
        assert get_tool_executable("bwa") == "bwa"

    def test_simple_tool_with_path(self):
        """Test extracting name from command with path."""
        assert get_tool_executable("/usr/bin/bwa") == "/usr/bin/bwa"

    def test_conda_command_basic(self):
        """Test extracting tool name from basic conda command."""
        assert get_tool_executable("mamba run -n env_wessim reseq") == "reseq"

    def test_conda_command_with_flags(self):
        """Test extracting tool name from conda command with flags."""
        tool_name = get_tool_executable("mamba run --no-capture-output -n env_wessim samtools")
        assert tool_name == "samtools"

    def test_conda_command_bwa(self):
        """Test extracting BWA from conda command."""
        assert get_tool_executable("mamba run -n env_wessim bwa") == "bwa"

    def test_conda_command_nanosim(self):
        """Test extracting NanoSim from conda command."""
        tool_name = get_tool_executable("mamba run -n env_nanosim nanosim")
        assert tool_name == "nanosim"

    def test_empty_command_returns_unknown(self):
        """Test that empty string returns 'unknown'."""
        assert get_tool_executable("") == "unknown"

    def test_none_command_returns_unknown(self):
        """Test that None returns 'unknown'."""
        assert get_tool_executable(None) == "unknown"  # type: ignore[arg-type]

    def test_command_only_flags(self):
        """Test command with only flags (edge case)."""
        # If all elements start with "-", should return the last one
        result = get_tool_executable("-t -n -v")
        assert result == "-v"

    def test_tool_with_version_flag(self):
        """Test command ending with version flag."""
        # Should return the tool name, not the flag
        assert get_tool_executable("bwa --version") == "bwa"

    def test_docker_wrapped_command(self):
        """Test extracting tool from Docker-wrapped command."""
        # Future extension: docker run -it container tool
        tool_name = get_tool_executable("docker run -it biocontainers/bwa bwa")
        assert tool_name == "bwa"

    def test_singularity_wrapped_command(self):
        """Test extracting tool from Singularity-wrapped command."""
        # Future extension: singularity exec image.sif tool
        tool_name = get_tool_executable("singularity exec bwa.sif bwa")
        assert tool_name == "bwa"


class TestCommandUtilsIntegration:
    """Integration tests for command_utils with real-world scenarios."""

    def test_integration_reseq_replace_ns(self):
        """Test building complete reseq replaceN command."""
        tools_dict = {"reseq": "mamba run --no-capture-output -n env_wessim reseq"}
        cmd = build_tool_command(
            tools_dict["reseq"], "replaceN", "-r", "input.fa", "-R", "output.fa"
        )

        # Verify correct structure
        assert cmd[:5] == ["mamba", "run", "--no-capture-output", "-n", "env_wessim"]
        assert cmd[5:] == ["reseq", "replaceN", "-r", "input.fa", "-R", "output.fa"]

        # Verify tool name extraction
        assert get_tool_executable(tools_dict["reseq"]) == "reseq"

    def test_integration_bwa_mem(self):
        """Test building complete BWA mem command."""
        tools_dict = {"bwa": "mamba run --no-capture-output -n env_wessim bwa"}
        threads = 4
        cmd = build_tool_command(
            tools_dict["bwa"], "mem", "-t", threads, "ref.fa", "r1.fq", "r2.fq"
        )

        # Verify correct structure
        assert "bwa" in cmd
        assert "mem" in cmd
        assert "-t" in cmd
        assert "4" in cmd  # Note: numeric arg converted to string

    def test_integration_samtools_view(self):
        """Test building complete samtools view command."""
        tools_dict = {"samtools": "mamba run --no-capture-output -n env_wessim samtools"}
        cmd = build_tool_command(tools_dict["samtools"], "view", "-b", "-o", "out.bam", "in.sam")

        # Verify correct structure
        assert cmd[:5] == ["mamba", "run", "--no-capture-output", "-n", "env_wessim"]
        assert cmd[5:] == ["samtools", "view", "-b", "-o", "out.bam", "in.sam"]

    def test_integration_simple_tool_no_conda(self):
        """Test building command for tool not wrapped in conda."""
        tools_dict = {"minimap2": "minimap2"}
        cmd = build_tool_command(tools_dict["minimap2"], "-ax", "map-ont", "ref.fa", "reads.fq")

        # Verify simple structure (no conda parsing)
        assert cmd == ["minimap2", "-ax", "map-ont", "ref.fa", "reads.fq"]
