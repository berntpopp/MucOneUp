# Error Handling & Exceptions

**Status:** ✅ COMPLETED - October 1, 2025
**Priority:** 🔴 URGENT
**Estimated Effort:** 1-2 days
**Actual Effort:** 1.5 days
**Impact:** HIGH - Enables testing, improves error recovery

## Problem Statement

Current error handling is inconsistent and untestable:

- **60 `sys.exit()` calls** across entire codebase
  - 29 in `cli.py`
  - 31 in other modules (wrappers, etc.)
- **No custom exception hierarchy** - only generic `ValueError`, `RuntimeError`
- **Mixed strategies:** Some functions raise exceptions, others call `sys.exit()`
- **Abrupt termination:** `sys.exit()` kills entire process, including test runners
- **No error recovery:** Cannot catch and handle errors programmatically

**Result:** Code is untestable, users get poor error messages, no error recovery possible.

## Core Principles

### KISS (Keep It Simple)

**Simple Rule:** Raise exceptions everywhere, handle them once at the top level.

```python
# ❌ Complex - error handling scattered everywhere
def parse_config(path):
    if not exists(path):
        logging.error("Config not found")
        sys.exit(1)  # Handling here

def load_data(config):
    if config is None:
        logging.error("No config")
        sys.exit(1)  # And here

def process():
    if something_wrong:
        logging.error("Failed")
        sys.exit(1)  # And here

# ✅ Simple - one place handles all errors
def parse_config(path):
    if not exists(path):
        raise ConfigurationError(f"Config not found: {path}")

def load_data(config):
    if config is None:
        raise ConfigurationError("Config is None")

def process():
    if something_wrong:
        raise ProcessingError("Failed")

def main():
    try:
        # All business logic
        parse_config(...)
        load_data(...)
        process()
    except MucOneUpError as e:
        logging.error(str(e))
        return 1  # Only exit point
```

### DRY (Don't Repeat Yourself)

**Current Duplication:**
```python
# Pattern repeated 60 times:
logging.error("Something failed: %s", error)
sys.exit(1)
```

**DRY Solution:**
```python
# Single exception handler in main()
try:
    run_cli()
except MucOneUpError as e:
    logging.error(str(e))
    return 1
```

### SOLID - Single Responsibility

**Each module should handle its domain logic, NOT process termination.**

- ❌ `mutate.py` should NOT decide to exit the program
- ✅ `mutate.py` should raise `MutationError`, let caller decide

## Custom Exception Hierarchy

Create `muc_one_up/exceptions.py`:

```python
"""Custom exception hierarchy for MucOneUp."""

class MucOneUpError(Exception):
    """Base exception for all MucOneUp errors.

    All custom exceptions inherit from this, allowing catch-all handling.
    """
    pass


class ConfigurationError(MucOneUpError):
    """Raised when configuration is invalid or cannot be loaded.

    Examples:
        - Missing config file
        - Invalid JSON syntax
        - Schema validation failure
        - Missing required fields
    """
    pass


class ValidationError(MucOneUpError):
    """Raised when input validation fails.

    Examples:
        - Invalid DNA sequence (non-ACGTN characters)
        - Invalid mutation targets
        - Invalid SNP positions
        - Invalid reference assembly
    """
    pass


class SimulationError(MucOneUpError):
    """Raised when simulation process fails.

    Examples:
        - Cannot generate haplotype chain
        - Invalid probability distribution
        - Seed initialization failure
    """
    pass


class MutationError(MucOneUpError):
    """Raised when mutation application fails.

    Examples:
        - Mutation not defined in config
        - Invalid mutation target
        - Strict mode violation
        - Mutation position out of range
    """
    pass


class SNPIntegrationError(MucOneUpError):
    """Raised when SNP integration fails.

    Examples:
        - Invalid SNP file format
        - Reference base mismatch
        - SNP position out of range
    """
    pass


class ExternalToolError(MucOneUpError):
    """Raised when external tool execution fails.

    Attributes:
        tool: Name of the tool that failed
        exit_code: Exit code returned by tool
        stderr: Error output from tool
        stdout: Standard output from tool (optional)
    """

    def __init__(
        self,
        tool: str,
        exit_code: int,
        stderr: str,
        stdout: str = "",
        cmd: str = "",
    ):
        self.tool = tool
        self.exit_code = exit_code
        self.stderr = stderr
        self.stdout = stdout
        self.cmd = cmd

        message = f"{tool} failed with exit code {exit_code}"
        if cmd:
            message += f"\nCommand: {cmd}"
        if stderr:
            message += f"\nError: {stderr[:500]}"  # Truncate long errors

        super().__init__(message)


class FileOperationError(MucOneUpError):
    """Raised when file operations fail.

    Examples:
        - Cannot read FASTA file
        - Cannot write output
        - Directory not found
        - Permission denied
    """
    pass


class ReadSimulationError(MucOneUpError):
    """Raised when read simulation fails.

    Examples:
        - Missing reference genome
        - Invalid read simulator parameters
        - Pipeline step failure
    """
    pass
```

## Migration Plan

### Step 1: Create Exception Module (30 minutes)

```bash
# Create the file
touch muc_one_up/exceptions.py
# Paste the hierarchy above
```

### Step 2: Replace sys.exit() in Core Modules (3-4 hours)

**Priority order:**
1. `config.py` (configuration errors)
2. `mutate.py` (mutation errors)
3. `simulate.py` (simulation errors)
4. `snp_integrator.py` (SNP errors)
5. `read_simulator/wrappers/*.py` (external tool errors)

**Example: config.py**

```python
# ❌ Before
def load_config(config_path: str) -> Dict[str, Any]:
    if not os.path.exists(config_path):
        logging.error(f"Config file not found: {config_path}")
        sys.exit(1)

    try:
        with open(config_path) as f:
            config = json.load(f)
    except json.JSONDecodeError as e:
        logging.error(f"Invalid JSON in config: {e}")
        sys.exit(1)

    try:
        validate(config, CONFIG_SCHEMA)
    except ValidationError as e:
        logging.error(f"Config validation failed: {e}")
        sys.exit(1)

    return config

# ✅ After
from .exceptions import ConfigurationError

def load_config(config_path: str) -> Dict[str, Any]:
    """Load and validate configuration file.

    Args:
        config_path: Path to JSON config file

    Returns:
        Validated configuration dictionary

    Raises:
        ConfigurationError: If config cannot be loaded or is invalid
    """
    if not os.path.exists(config_path):
        raise ConfigurationError(f"Config file not found: {config_path}")

    try:
        with open(config_path) as f:
            config = json.load(f)
    except json.JSONDecodeError as e:
        raise ConfigurationError(f"Invalid JSON in config: {e}") from e

    try:
        validate(config, CONFIG_SCHEMA)
    except ValidationError as e:
        raise ConfigurationError(f"Config validation failed: {e}") from e

    return config
```

**Example: External Tool Wrappers**

```python
# ❌ Before (samtools_wrapper.py)
def run_samtools(cmd: List[str]) -> None:
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logging.error(f"samtools failed: {result.stderr}")
        sys.exit(1)

# ✅ After
from ..exceptions import ExternalToolError

def run_samtools(cmd: List[str]) -> None:
    """Run samtools command.

    Args:
        cmd: Command list to execute

    Raises:
        ExternalToolError: If samtools fails
    """
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise ExternalToolError(
            tool="samtools",
            exit_code=result.returncode,
            stderr=result.stderr,
            stdout=result.stdout,
            cmd=" ".join(cmd),
        )
```

### Step 3: Update CLI Main Entry Point (1 hour)

```python
# cli.py or cli/__init__.py

from .exceptions import MucOneUpError, ExternalToolError
import logging
import sys

def main() -> int:
    """Main entry point.

    Returns:
        Exit code (0 for success, 1 for expected errors, 2 for unexpected)
    """
    args = build_parser().parse_args()
    configure_logging(args.log_level)

    try:
        # All business logic here
        config, out_dir = setup_configuration(args)
        sim_configs = determine_simulation_mode(args)

        for sim_index, fixed_conf in enumerate(sim_configs, 1):
            run_single_simulation(args, config, out_dir, fixed_conf, sim_index)

        return 0  # Success

    except ConfigurationError as e:
        logging.error("Configuration error: %s", e)
        return 1

    except ValidationError as e:
        logging.error("Validation error: %s", e)
        return 1

    except MutationError as e:
        logging.error("Mutation error: %s", e)
        return 1

    except SimulationError as e:
        logging.error("Simulation error: %s", e)
        return 1

    except ExternalToolError as e:
        logging.error("External tool error: %s", e)
        if e.stderr:
            logging.debug("Tool stderr: %s", e.stderr)
        return 1

    except MucOneUpError as e:
        # Catch-all for custom exceptions
        logging.error("MucOneUp error: %s", e)
        return 1

    except KeyboardInterrupt:
        logging.warning("Interrupted by user")
        return 130  # Standard Unix code for SIGINT

    except Exception as e:
        # Unexpected errors - include traceback
        logging.exception("Unexpected error occurred: %s", e)
        return 2


# Entry point for setup.py
if __name__ == "__main__":
    sys.exit(main())
```

### Step 4: Update Tests (2 hours)

**Before:** Cannot test due to `sys.exit()`
```python
# ❌ This kills the test runner!
def test_invalid_config():
    with pytest.raises(SystemExit):
        load_config("nonexistent.json")
```

**After:** Test exception behavior
```python
# ✅ Properly tests exception
def test_invalid_config():
    with pytest.raises(ConfigurationError, match="Config file not found"):
        load_config("nonexistent.json")

def test_invalid_json():
    with pytest.raises(ConfigurationError, match="Invalid JSON"):
        load_config("tests/data/invalid.json")

def test_schema_validation_failure():
    with pytest.raises(ConfigurationError, match="validation failed"):
        load_config("tests/data/missing_required.json")
```

**Test exception attributes:**
```python
def test_external_tool_error_attributes():
    with pytest.raises(ExternalToolError) as exc_info:
        run_samtools(["samtools", "view", "nonexistent.bam"])

    error = exc_info.value
    assert error.tool == "samtools"
    assert error.exit_code != 0
    assert "nonexistent.bam" in error.stderr
    assert error.cmd is not None
```

## Error Message Quality

### Bad Error Messages (Current)
```
ERROR - Simulation failed
ERROR - Tool failed
ERROR - Invalid input
```

### Good Error Messages (Target)
```
ERROR - Configuration error: Config file not found: /path/to/config.json
ERROR - Mutation error: Mutation 'dupX' not defined in config. Available mutations: ['dupC', 'delA', 'normal']
ERROR - External tool error: samtools failed with exit code 1
  Command: samtools view -b input.sam -o output.bam
  Error: [samtools] failed to open input.sam: No such file or directory
ERROR - Validation error: Invalid DNA sequence at position 150: Found 'X' (expected ACGTN)
```

**Principles for good error messages:**
1. ✅ Include context (what failed, where, why)
2. ✅ Suggest fixes when possible
3. ✅ Include relevant values
4. ✅ Truncate very long output
5. ✅ Preserve original exception chain with `from e`

## Benefits of This Approach

| Benefit | Before | After |
|---------|--------|-------|
| **Testability** | Cannot test - `sys.exit()` kills test runner | All error paths testable with `pytest.raises()` |
| **Error Recovery** | Impossible - program terminates | Caller can catch and handle |
| **Error Messages** | Generic, unclear | Specific, actionable |
| **Debugging** | Stack trace lost at `sys.exit()` | Full exception chain preserved |
| **Programmatic Use** | Cannot use as library | Can import and handle errors |
| **Consistency** | 60 different exit points | Single exit point in `main()` |

## Checklist

### Implementation Checklist

- [x] Create `muc_one_up/exceptions.py` with hierarchy
- [x] Replace `sys.exit()` in `config.py` → `ConfigurationError`
- [x] Replace `sys.exit()` in `mutate.py` → `MutationError`
- [x] Replace `sys.exit()` in `simulate.py` → `SimulationError`
- [x] Replace `sys.exit()` in `snp_integrator.py` → `SNPIntegrationError`
- [x] Replace `sys.exit()` in `read_simulator/wrappers/*.py` → `ExternalToolError`
- [x] Replace `sys.exit()` in other modules → appropriate exception
- [x] Update `cli.py:main()` with exception handling
- [x] Verify `sys.exit()` only in `if __name__ == "__main__"` block
- [x] Add tests for all exception types
- [x] Update docstrings with `Raises:` sections
- [x] Verify all tests pass

### Verification Commands

```bash
# Find remaining sys.exit() calls (should only be in main entry point)
grep -rn "sys.exit" muc_one_up/ --include="*.py" | grep -v "# OK: top-level"

# Count should be 1-2 (only in main entry)
grep -r "sys.exit" muc_one_up/ --include="*.py" | wc -l

# Run tests - should all pass
pytest tests/ -v

# Test exception handling specifically
pytest tests/test_exceptions.py -v
```

## Testing Exception Handling

Create `tests/test_exceptions.py`:

```python
"""Tests for exception handling across the codebase."""

import pytest
from muc_one_up.exceptions import (
    MucOneUpError,
    ConfigurationError,
    MutationError,
    ExternalToolError,
)
from muc_one_up.config import load_config
from muc_one_up.mutate import apply_mutation


def test_configuration_error_hierarchy():
    """ConfigurationError should inherit from MucOneUpError."""
    error = ConfigurationError("test")
    assert isinstance(error, MucOneUpError)
    assert isinstance(error, Exception)


def test_load_config_missing_file():
    """Loading nonexistent config should raise ConfigurationError."""
    with pytest.raises(ConfigurationError, match="Config file not found"):
        load_config("nonexistent.json")


def test_external_tool_error_has_attributes():
    """ExternalToolError should preserve tool information."""
    error = ExternalToolError(
        tool="samtools",
        exit_code=1,
        stderr="Error message",
        cmd="samtools view"
    )

    assert error.tool == "samtools"
    assert error.exit_code == 1
    assert error.stderr == "Error message"
    assert "samtools failed with exit code 1" in str(error)


def test_no_sys_exit_in_modules():
    """No module except main entry point should call sys.exit()."""
    # This test scans codebase to ensure sys.exit() was removed
    import os
    import re

    violations = []
    for root, dirs, files in os.walk("muc_one_up"):
        for file in files:
            if not file.endswith(".py"):
                continue

            filepath = os.path.join(root, file)
            with open(filepath) as f:
                for line_num, line in enumerate(f, 1):
                    if "sys.exit" in line and "if __name__" not in line:
                        violations.append(f"{filepath}:{line_num}")

    # Allow sys.exit only in main entry points
    allowed = ["cli.py", "cli/__init__.py"]
    violations = [v for v in violations if not any(a in v for a in allowed)]

    assert not violations, f"Found sys.exit() in modules: {violations}"
```

## Migration Priority

**Day 1: Core modules**
1. Create `exceptions.py`
2. Update `config.py`
3. Update `mutate.py`
4. Update `simulate.py`
5. Add tests

**Day 2: External tools & CLI**
1. Update all `read_simulator/wrappers/*.py`
2. Update `snp_integrator.py`
3. Update `cli.py` main()
4. Add comprehensive exception tests
5. Verify all tests pass

## Success Criteria

- ✅ Zero `sys.exit()` calls outside of `if __name__ == "__main__"`
- ✅ All modules raise custom exceptions
- ✅ `main()` has centralized exception handling
- ✅ All exception types are tested
- ✅ All docstrings include `Raises:` sections
- ✅ Error messages are clear and actionable
- ✅ Full test suite passes

## Next Steps

After completing exception handling:
1. ➡️ Proceed to **Testing Strategy** (now possible to test error paths)
2. ➡️ Add **Type Hints** to exception hierarchy
3. ➡️ Document error handling in developer guide
