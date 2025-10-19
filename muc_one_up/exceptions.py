"""Custom exception hierarchy for MucOneUp.

This module defines a comprehensive exception hierarchy following SOLID principles:
- Single Responsibility: Each exception type represents one category of errors
- Interface Segregation: Specific exception types for different error scenarios

All exceptions inherit from MucOneUpError for catch-all handling.
"""


class MucOneUpError(Exception):
    """Base exception for all MucOneUp errors.

    All custom exceptions inherit from this, allowing catch-all handling
    while preserving the ability to catch specific error types.
    """


class ConfigurationError(MucOneUpError):
    """Raised when configuration is invalid or cannot be loaded.

    Examples:
        - Missing config file
        - Invalid JSON syntax
        - Schema validation failure
        - Missing required fields
        - Invalid reference assembly
    """


class ValidationError(MucOneUpError):
    """Raised when input validation fails.

    Examples:
        - Invalid DNA sequence (non-ACGTN characters)
        - Invalid mutation targets
        - Invalid SNP positions
        - Invalid reference assembly
        - Invalid command-line arguments
    """


class SimulationError(MucOneUpError):
    """Raised when simulation process fails.

    Examples:
        - Cannot generate haplotype chain
        - Invalid probability distribution
        - Seed initialization failure
        - Structure file parsing failure
    """


class MutationError(MucOneUpError):
    """Raised when mutation application fails.

    Examples:
        - Mutation not defined in config
        - Invalid mutation target
        - Strict mode violation
        - Mutation position out of range
        - Cannot find random target
    """


class SNPIntegrationError(MucOneUpError):
    """Raised when SNP integration fails.

    Examples:
        - Invalid SNP file format
        - Reference base mismatch
        - SNP position out of range
        - Cannot parse SNP file
    """


class ExternalToolError(MucOneUpError):
    """Raised when external tool execution fails.

    This exception preserves detailed information about tool failures
    for debugging and user feedback.

    Attributes:
        tool: Name of the tool that failed
        exit_code: Exit code returned by tool
        stderr: Error output from tool
        stdout: Standard output from tool (optional)
        cmd: Command that was executed
    """

    def __init__(
        self,
        tool: str,
        exit_code: int,
        stderr: str = "",
        stdout: str = "",
        cmd: str = "",
    ):
        """Initialize ExternalToolError.

        Args:
            tool: Name of the external tool (e.g., "samtools", "bwa")
            exit_code: Exit code returned by the tool
            stderr: Error output from the tool
            stdout: Standard output from the tool
            cmd: Full command that was executed
        """
        self.tool = tool
        self.exit_code = exit_code
        self.stderr = stderr
        self.stdout = stdout
        self.cmd = cmd

        # Build informative error message
        message = f"{tool} failed with exit code {exit_code}"
        if cmd:
            message += f"\nCommand: {cmd}"
        if stderr:
            # Truncate very long error messages
            truncated_stderr = stderr[:500] + "..." if len(stderr) > 500 else stderr
            message += f"\nError: {truncated_stderr}"

        super().__init__(message)


class FileOperationError(MucOneUpError):
    """Raised when file operations fail.

    Examples:
        - Cannot read FASTA file
        - Cannot write output
        - Directory not found
        - Permission denied
        - I/O error
    """


class ReadSimulationError(MucOneUpError):
    """Raised when read simulation fails.

    Examples:
        - Missing reference genome
        - Invalid read simulator parameters
        - Pipeline step failure
        - Missing required tools
    """
