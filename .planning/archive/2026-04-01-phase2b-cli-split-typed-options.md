# Phase 2B: Split CLI God Module & Typed Options — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Split the 1536-line `click_main.py` into focused per-command modules, replace the fake `argparse.SimpleNamespace` with a typed dataclass, and centralize error handling with a decorator.

**Architecture:** Extract shared CLI utilities to `_common.py`. Create an error-handling decorator that maps domain exceptions to exit codes. Define `SimulationOptions` dataclass matching the current `SimpleNamespace` fields. Move commands into `commands/simulate.py`, `commands/reads.py`, and `commands/analyze.py`. Slim `click_main.py` to ~60 lines of group registration. All orchestration functions continue to work because `SimulationOptions` exposes the same attribute interface.

**Tech Stack:** Python 3.10+, Click, dataclasses, pytest

**Spec:** `.planning/specs/2026-04-01-codebase-refactoring-design.md` (sections 2.3, 2.4, 2.5)

---

## File Structure

### New Files
| File | Responsibility |
|------|---------------|
| `muc_one_up/cli/_common.py` | CONTEXT_SETTINGS, configure_logging, validate_config_callback, require_config, print_version_callback |
| `muc_one_up/cli/error_handling.py` | `cli_error_handler` decorator mapping exceptions to exit codes |
| `muc_one_up/cli/options.py` | `SimulationOptions` frozen dataclass |
| `muc_one_up/cli/commands/__init__.py` | Empty package init |
| `muc_one_up/cli/commands/simulate.py` | `simulate` command using SimulationOptions + error decorator |
| `muc_one_up/cli/commands/reads.py` | `reads` group + illumina, ont, pacbio subcommands |
| `muc_one_up/cli/commands/analyze.py` | `analyze` group + orfs, stats, vntr-stats, snapshot-validate |

### Modified Files
| File | Changes |
|------|---------|
| `muc_one_up/cli/click_main.py` | Slim to ~60 lines: root group + command registration |
| `muc_one_up/cli/orchestration.py` | Change `args` type annotation to `SimulationOptions` |

### Test Files
| File | Tests |
|------|-------|
| `tests/test_cli_error_handling.py` | Error decorator behavior |
| `tests/test_simulation_options.py` | SimulationOptions construction and from_click_kwargs |
| `tests/test_cli_module_structure.py` | Structural tests: click_main.py line count, no SimpleNamespace |

---

## Task 1: Extract Shared Utilities to _common.py

**Files:**
- Create: `muc_one_up/cli/_common.py`
- Modify: `muc_one_up/cli/click_main.py`

- [ ] **Step 1: Create _common.py with shared utilities**

Extract these functions from `click_main.py` (lines 36-144) into `_common.py`:

```python
# muc_one_up/cli/_common.py
"""Shared CLI utilities: settings, logging, config validation, version display."""

import logging
from pathlib import Path

import click

from ..version import __version__

CONTEXT_SETTINGS = {"help_option_names": ["-h", "--help"]}


def configure_logging(level_str: str) -> None:
    """Configure Python logging based on level string."""
    if level_str == "NONE":
        logging.disable(logging.CRITICAL)
        return

    level = getattr(logging, level_str.upper(), logging.WARNING)
    logging.basicConfig(
        level=level,
        format="%(asctime)s - %(levelname)s - %(message)s",
        force=True,
    )
    # Silence noisy third-party loggers
    for noisy in ("urllib3", "httpx", "httpcore"):
        logging.getLogger(noisy).setLevel(logging.WARNING)


def validate_config_callback(
    ctx: click.Context, param: click.Parameter, value: str | None
) -> str | None:
    """Click callback for --config option. Stores path in ctx.obj."""
    if value is not None:
        config_path = Path(value)
        if not config_path.exists():
            raise click.BadParameter(f"Config file not found: {config_path}")
        ctx.ensure_object(dict)
        ctx.obj["config_path"] = str(config_path.resolve())
    return value


def require_config(ctx: click.Context) -> None:
    """Raise UsageError if --config was not provided."""
    if not ctx.obj or "config_path" not in ctx.obj:
        raise click.UsageError(
            "Missing option '--config'. "
            "Please provide a configuration file: muconeup --config <path>"
        )


def print_version_callback(
    ctx: click.Context, param: click.Parameter, value: bool
) -> None:
    """Click callback for --version flag."""
    if not value or ctx.resilient_parsing:
        return
    click.echo(f"MucOneUp {__version__}")
    ctx.exit()
```

Note: Copy the exact implementations from `click_main.py` lines 45-144. The code above shows the structure — match the exact current implementations.

- [ ] **Step 2: Update click_main.py imports**

In `click_main.py`, replace the inline definitions with imports from `_common`:

```python
# REPLACE the inline CONTEXT_SETTINGS, configure_logging, validate_config_callback,
# require_config, print_version_callback definitions (lines 36-144)
# WITH imports:
from ._common import (
    CONTEXT_SETTINGS,
    configure_logging,
    print_version_callback,
    require_config,
    validate_config_callback,
)
```

Remove the deleted function bodies but keep the `cli()` group and everything after it.

- [ ] **Step 3: Run tests**

Run: `pytest tests/test_cli_smoke.py tests/test_click_cli.py --tb=short -q --no-cov`
Expected: All PASS (behavior unchanged, just code moved)

- [ ] **Step 4: Commit**

```bash
git add muc_one_up/cli/_common.py muc_one_up/cli/click_main.py
git commit -m "refactor: extract shared CLI utilities to _common.py"
```

---

## Task 2: Create Error Handling Decorator

**Files:**
- Create: `muc_one_up/cli/error_handling.py`
- Create: `tests/test_cli_error_handling.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/test_cli_error_handling.py
"""Tests for CLI error handling decorator."""

import logging

import click
from click.testing import CliRunner

from muc_one_up.cli.error_handling import cli_error_handler
from muc_one_up.exceptions import (
    ConfigurationError,
    ExternalToolError,
    FileOperationError,
    MucOneUpError,
    ReadSimulationError,
    SimulationError,
    ValidationError,
)


def _make_test_command(exception_to_raise):
    """Create a Click command that raises the given exception."""

    @click.command()
    @cli_error_handler
    def cmd():
        raise exception_to_raise

    return cmd


class TestCLIErrorHandler:
    def test_keyboard_interrupt_exits_130(self):
        runner = CliRunner()
        cmd = _make_test_command(KeyboardInterrupt())
        result = runner.invoke(cmd)
        assert result.exit_code == 130

    def test_muc_one_up_error_exits_1(self):
        runner = CliRunner()
        cmd = _make_test_command(MucOneUpError("test error"))
        result = runner.invoke(cmd)
        assert result.exit_code == 1

    def test_simulation_error_exits_1(self):
        runner = CliRunner()
        cmd = _make_test_command(SimulationError("sim failed"))
        result = runner.invoke(cmd)
        assert result.exit_code == 1

    def test_configuration_error_exits_1(self):
        runner = CliRunner()
        cmd = _make_test_command(ConfigurationError("bad config"))
        result = runner.invoke(cmd)
        assert result.exit_code == 1

    def test_generic_exception_exits_2(self):
        runner = CliRunner()
        cmd = _make_test_command(RuntimeError("unexpected"))
        result = runner.invoke(cmd)
        assert result.exit_code == 2

    def test_no_exception_exits_0(self):
        @click.command()
        @cli_error_handler
        def cmd():
            click.echo("ok")

        runner = CliRunner()
        result = runner.invoke(cmd)
        assert result.exit_code == 0
        assert "ok" in result.output
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_cli_error_handling.py -v`
Expected: FAIL with ImportError

- [ ] **Step 3: Implement error handling decorator**

```python
# muc_one_up/cli/error_handling.py
"""Centralized CLI error handling.

Maps domain exceptions to consistent exit codes:
- KeyboardInterrupt → 130
- MucOneUpError (and subclasses) → 1
- Unexpected exceptions → 2
"""

from __future__ import annotations

import functools
import logging

import click

from ..exceptions import MucOneUpError


def cli_error_handler(func):
    """Decorator that catches exceptions and maps them to CLI exit codes."""

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        ctx = click.get_current_context()
        try:
            return func(*args, **kwargs)
        except KeyboardInterrupt:
            logging.info("Operation cancelled by user")
            ctx.exit(130)
        except MucOneUpError as e:
            logging.error("%s: %s", type(e).__name__, e)
            ctx.exit(1)
        except Exception as e:
            logging.error("Unexpected error: %s", e, exc_info=True)
            ctx.exit(2)

    return wrapper
```

- [ ] **Step 4: Run tests**

Run: `pytest tests/test_cli_error_handling.py -v`
Expected: All PASS

- [ ] **Step 5: Commit**

```bash
git add muc_one_up/cli/error_handling.py tests/test_cli_error_handling.py
git commit -m "feat: add cli_error_handler decorator for centralized error handling"
```

---

## Task 3: Create SimulationOptions Dataclass

**Files:**
- Create: `muc_one_up/cli/options.py`
- Create: `tests/test_simulation_options.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/test_simulation_options.py
"""Tests for SimulationOptions typed dataclass."""

from muc_one_up.cli.options import SimulationOptions


class TestSimulationOptions:
    def test_basic_construction(self):
        opts = SimulationOptions(
            config="/path/to/config.json",
            out_base="muc1_simulated",
            out_dir=".",
            num_haplotypes=2,
            seed=42,
            reference_assembly="hg38",
        )
        assert opts.config == "/path/to/config.json"
        assert opts.seed == 42
        assert opts.num_haplotypes == 2

    def test_defaults(self):
        opts = SimulationOptions(
            config="/path/to/config.json",
            out_base="test",
            out_dir=".",
        )
        assert opts.num_haplotypes == 2
        assert opts.seed is None
        assert opts.reference_assembly == "hg38"
        assert opts.output_structure is False
        assert opts.fixed_lengths is None
        assert opts.input_structure is None
        assert opts.simulate_series is None
        assert opts.mutation_name is None
        assert opts.mutation_targets is None
        assert opts.simulate_reads is None
        assert opts.output_orfs is False
        assert opts.orf_min_aa == 100
        assert opts.orf_aa_prefix is None
        assert opts.track_read_source is False

    def test_from_click_kwargs(self):
        kwargs = {
            "out_base": "test",
            "out_dir": "output",
            "num_haplotypes": 4,
            "seed": 123,
            "reference_assembly": "hg19",
            "output_structure": True,
            "fixed_lengths": ("20", "30"),
            "input_structure": None,
            "simulate_series": None,
            "mutation_name": "dupC",
            "mutation_targets": (),
            "snp_input_file": None,
            "random_snps": False,
            "random_snp_density": 0.001,
            "random_snp_output_file": None,
            "random_snp_region": "constants_only",
            "random_snp_haplotypes": "all",
            "track_read_source": True,
        }
        opts = SimulationOptions.from_click_kwargs("/path/config.json", kwargs)
        assert opts.config == "/path/config.json"
        assert opts.out_base == "test"
        assert opts.num_haplotypes == 4
        assert opts.seed == 123
        assert opts.fixed_lengths == ["20", "30"]
        assert opts.mutation_targets is None  # empty tuple → None
        assert opts.track_read_source is True

    def test_from_click_kwargs_with_mutation_targets(self):
        kwargs = {
            "out_base": "test",
            "out_dir": ".",
            "num_haplotypes": 2,
            "seed": None,
            "reference_assembly": "hg38",
            "output_structure": False,
            "fixed_lengths": (),
            "input_structure": None,
            "simulate_series": None,
            "mutation_name": None,
            "mutation_targets": ("1,3", "2,5"),
            "snp_input_file": None,
            "random_snps": False,
            "random_snp_density": 0.001,
            "random_snp_output_file": None,
            "random_snp_region": "constants_only",
            "random_snp_haplotypes": "all",
            "track_read_source": False,
        }
        opts = SimulationOptions.from_click_kwargs("/config.json", kwargs)
        assert opts.mutation_targets == ["1,3", "2,5"]
        assert opts.fixed_lengths is None  # empty tuple → None
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_simulation_options.py -v`
Expected: FAIL with ImportError

- [ ] **Step 3: Implement SimulationOptions**

```python
# muc_one_up/cli/options.py
"""Typed option dataclasses for CLI commands.

Replaces the untyped argparse.SimpleNamespace with explicit, documented fields.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class SimulationOptions:
    """Typed options for the simulate command and orchestration pipeline.

    Replaces the SimpleNamespace created by _make_args_namespace().
    All orchestration functions (config.py, haplotypes.py, mutations.py,
    outputs.py, analysis.py, snps.py) access these fields by attribute name.
    """

    # Required
    config: str
    out_base: str
    out_dir: str

    # Simulation parameters
    num_haplotypes: int = 2
    seed: int | None = None
    reference_assembly: str = "hg38"

    # Output control
    output_structure: bool = False

    # Mode selection
    fixed_lengths: list[str] | None = None
    input_structure: str | None = None
    simulate_series: int | None = None

    # Mutation config
    mutation_name: str | None = None
    mutation_targets: list[str] | None = None

    # SNP config
    snp_input_file: str | None = None
    random_snps: bool = False
    random_snp_density: float = 0.001
    random_snp_output_file: str | None = None
    random_snp_region: str = "constants_only"
    random_snp_haplotypes: str = "all"

    # Pipeline flags (set to None/False for pure simulate)
    simulate_reads: str | None = None
    output_orfs: bool = False
    orf_min_aa: int = 100
    orf_aa_prefix: str | None = None
    track_read_source: bool = False

    @classmethod
    def from_click_kwargs(cls, config_path: str, kwargs: dict) -> SimulationOptions:
        """Construct from Click command kwargs.

        Converts Click tuples to lists and empty tuples to None.
        """
        fixed = list(kwargs["fixed_lengths"]) if kwargs["fixed_lengths"] else None
        targets = list(kwargs["mutation_targets"]) if kwargs["mutation_targets"] else None

        return cls(
            config=config_path,
            out_base=kwargs["out_base"],
            out_dir=kwargs["out_dir"],
            num_haplotypes=kwargs["num_haplotypes"],
            seed=kwargs["seed"],
            reference_assembly=kwargs["reference_assembly"],
            output_structure=kwargs["output_structure"],
            fixed_lengths=fixed,
            input_structure=kwargs["input_structure"],
            simulate_series=kwargs["simulate_series"],
            mutation_name=kwargs["mutation_name"],
            mutation_targets=targets,
            snp_input_file=kwargs["snp_input_file"],
            random_snps=kwargs["random_snps"],
            random_snp_density=kwargs["random_snp_density"],
            random_snp_output_file=kwargs["random_snp_output_file"],
            random_snp_region=kwargs["random_snp_region"],
            random_snp_haplotypes=kwargs["random_snp_haplotypes"],
            track_read_source=kwargs.get("track_read_source", False),
        )
```

- [ ] **Step 4: Run tests**

Run: `pytest tests/test_simulation_options.py -v`
Expected: All PASS

- [ ] **Step 5: Commit**

```bash
git add muc_one_up/cli/options.py tests/test_simulation_options.py
git commit -m "feat: add SimulationOptions typed dataclass replacing SimpleNamespace"
```

---

## Task 4: Move simulate Command to commands/simulate.py

**Files:**
- Create: `muc_one_up/cli/commands/__init__.py`
- Create: `muc_one_up/cli/commands/simulate.py`
- Modify: `muc_one_up/cli/click_main.py`

- [ ] **Step 1: Create commands package**

```python
# muc_one_up/cli/commands/__init__.py
"""CLI command modules."""
```

- [ ] **Step 2: Create commands/simulate.py**

Move the `simulate` command (lines 225-412 of click_main.py) into `commands/simulate.py`. Replace `_make_args_namespace` with `SimulationOptions.from_click_kwargs`. Apply the `@cli_error_handler` decorator.

```python
# muc_one_up/cli/commands/simulate.py
"""Simulate command: generate MUC1 VNTR diploid haplotypes."""

from __future__ import annotations

import logging
from contextlib import nullcontext

import click

from ..error_handling import cli_error_handler
from ..options import SimulationOptions
```

The command itself: copy the full `@click.command()` decorated function from `click_main.py` lines 225-412, with these changes:

1. Replace `@cli.command()` with just defining the function (it will be registered via `cli.add_command()` in click_main.py).
2. Replace the `_make_args_namespace(config_path, kwargs)` call at line 353 with:
   ```python
   args = SimulationOptions.from_click_kwargs(config_path, kwargs)
   ```
3. Replace the try/except at lines 351-412 with `@cli_error_handler` on the function.
4. Import `setup_configuration`, `determine_simulation_mode`, `process_mutation_config` from `..config`.
5. Import `run_single_simulation_iteration` from `..orchestration`.

The function should be decorated with `@click.command()` and `@cli_error_handler`.

- [ ] **Step 3: Register in click_main.py**

In `click_main.py`, replace the full simulate command definition with:

```python
from .commands.simulate import simulate
cli.add_command(simulate)
```

- [ ] **Step 4: Run tests**

Run: `pytest tests/test_cli_smoke.py tests/test_click_cli.py --tb=short -q --no-cov`
Expected: All PASS

- [ ] **Step 5: Delete _make_args_namespace from click_main.py**

Remove lines 1360-1389 (`_make_args_namespace` function) since it's replaced by `SimulationOptions.from_click_kwargs`.

- [ ] **Step 6: Run full test suite**

Run: `pytest tests/ --tb=short -q --no-cov`
Expected: All PASS

- [ ] **Step 7: Commit**

```bash
git add muc_one_up/cli/commands/ muc_one_up/cli/click_main.py
git commit -m "refactor: move simulate command to commands/simulate.py with typed options"
```

---

## Task 5: Move reads Commands to commands/reads.py

**Files:**
- Create: `muc_one_up/cli/commands/reads.py`
- Modify: `muc_one_up/cli/click_main.py`

- [ ] **Step 1: Create commands/reads.py**

Move the `reads` group and its 3 subcommands (illumina, ont, pacbio) from `click_main.py` lines 420-996. Apply `@cli_error_handler` to each subcommand.

```python
# muc_one_up/cli/commands/reads.py
"""Read simulation commands: illumina, ont, pacbio."""

from __future__ import annotations

import logging
from pathlib import Path

import click

from ..error_handling import cli_error_handler
from ..outputs import generate_output_base
from .._common import require_config
```

Define the `reads` group and all 3 subcommands. Each subcommand:
1. Keeps its `@reads.command()` decorator
2. Gets `@cli_error_handler` added (replacing the try/except block)
3. Removes the manual try/except error handling
4. Keeps all Click option decorators exactly as-is

- [ ] **Step 2: Register in click_main.py**

In `click_main.py`, replace the full reads group + 3 subcommands with:

```python
from .commands.reads import reads
cli.add_command(reads)
```

- [ ] **Step 3: Run tests**

Run: `pytest tests/test_cli_smoke.py tests/test_click_cli.py --tb=short -q --no-cov`
Expected: All PASS

- [ ] **Step 4: Commit**

```bash
git add muc_one_up/cli/commands/reads.py muc_one_up/cli/click_main.py
git commit -m "refactor: move reads commands to commands/reads.py with error decorator"
```

---

## Task 6: Move analyze Commands to commands/analyze.py

**Files:**
- Create: `muc_one_up/cli/commands/analyze.py`
- Modify: `muc_one_up/cli/click_main.py`

- [ ] **Step 1: Create commands/analyze.py**

Move the `analyze` group and its 4 subcommands (orfs, stats, vntr-stats, snapshot-validate) from `click_main.py` lines 1004-1522. Apply `@cli_error_handler` to each subcommand.

```python
# muc_one_up/cli/commands/analyze.py
"""Analysis commands: orfs, stats, vntr-stats, snapshot-validate."""

from __future__ import annotations

import json
import logging
from pathlib import Path

import click

from ..error_handling import cli_error_handler
from ..outputs import generate_output_base
from .._common import require_config
```

Each subcommand:
1. Gets `@cli_error_handler` added (replacing try/except)
2. Removes manual try/except error handling
3. Keeps all Click option decorators exactly as-is

- [ ] **Step 2: Register in click_main.py**

In `click_main.py`, replace the full analyze group + 4 subcommands with:

```python
from .commands.analyze import analyze
cli.add_command(analyze)
```

- [ ] **Step 3: Run tests**

Run: `pytest tests/test_cli_smoke.py tests/test_click_cli.py --tb=short -q --no-cov`
Expected: All PASS

- [ ] **Step 4: Commit**

```bash
git add muc_one_up/cli/commands/analyze.py muc_one_up/cli/click_main.py
git commit -m "refactor: move analyze commands to commands/analyze.py with error decorator"
```

---

## Task 7: Slim click_main.py to Registration Only

**Files:**
- Modify: `muc_one_up/cli/click_main.py`

- [ ] **Step 1: Write structural test**

```python
# tests/test_cli_module_structure.py
"""Structural tests for CLI module organization."""

import inspect


def test_click_main_under_100_lines():
    """click_main.py should be registration-only, under 100 lines."""
    import muc_one_up.cli.click_main as mod

    source = inspect.getsource(mod)
    line_count = len(source.splitlines())
    assert line_count < 100, (
        f"click_main.py is {line_count} lines. "
        f"Should be under 100 (registration only)."
    )


def test_no_simple_namespace_in_codebase():
    """No SimpleNamespace usage anywhere in CLI code."""
    import muc_one_up.cli.click_main as mod

    source = inspect.getsource(mod)
    assert "SimpleNamespace" not in source


def test_commands_registered():
    """All command groups should be registered on the root CLI."""
    from muc_one_up.cli.click_main import cli

    command_names = list(cli.commands.keys())
    assert "simulate" in command_names
    assert "reads" in command_names
    assert "analyze" in command_names
```

- [ ] **Step 2: Verify click_main.py is clean**

After Tasks 4-6, `click_main.py` should only contain:
1. Imports from `_common` and `commands`
2. The `@click.group()` `cli()` function
3. `cli.add_command()` registrations
4. The `main()` entry point

If there's leftover code (helper functions, command definitions), remove it now.

The final `click_main.py` should look approximately like:

```python
"""MucOneUp CLI — command registration."""

import click

from ..version import __version__
from ._common import (
    CONTEXT_SETTINGS,
    configure_logging,
    print_version_callback,
    validate_config_callback,
)


@click.group(context_settings=CONTEXT_SETTINGS)
@click.option(
    "--config", "-c",
    type=click.Path(exists=True),
    callback=validate_config_callback,
    expose_value=False,
    is_eager=True,
    help="Path to JSON configuration file.",
)
@click.option(
    "--log-level",
    type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR", "NONE"]),
    default="INFO",
    callback=lambda ctx, param, value: configure_logging(value),
    expose_value=False,
    is_eager=True,
    help="Set logging verbosity.",
)
@click.option(
    "--verbose", "-v",
    is_flag=True,
    callback=lambda ctx, param, value: configure_logging("DEBUG") if value else None,
    expose_value=False,
    is_eager=True,
    help="Enable verbose (DEBUG) logging.",
)
@click.option(
    "--version", "-V",
    is_flag=True,
    callback=print_version_callback,
    expose_value=False,
    is_eager=True,
    help="Show version and exit.",
)
@click.pass_context
def cli(ctx):
    """MucOneUp: Simulate MUC1 VNTR diploid references."""
    ctx.ensure_object(dict)


# Register command groups
from .commands.analyze import analyze  # noqa: E402
from .commands.reads import reads  # noqa: E402
from .commands.simulate import simulate  # noqa: E402

cli.add_command(simulate)
cli.add_command(reads)
cli.add_command(analyze)


def main():
    """Entry point for the CLI."""
    cli()
```

- [ ] **Step 3: Run structural tests**

Run: `pytest tests/test_cli_module_structure.py -v`
Expected: All PASS

- [ ] **Step 4: Run full test suite**

Run: `pytest tests/ --tb=short -q --no-cov`
Expected: All PASS

- [ ] **Step 5: Commit**

```bash
git add muc_one_up/cli/click_main.py tests/test_cli_module_structure.py
git commit -m "refactor: slim click_main.py to registration-only (~60 lines)"
```

---

## Task 8: Verify Full Suite and CI Readiness

**Files:** None (verification only)

- [ ] **Step 1: Run full test suite**

Run: `pytest --tb=short -q`
Expected: All tests pass.

- [ ] **Step 2: Run ruff linter**

Run: `ruff check muc_one_up/ tests/`
Expected: No errors.

- [ ] **Step 3: Run ruff formatter**

Run: `ruff format --check muc_one_up/ tests/`
Expected: No formatting issues. If any, run `ruff format muc_one_up/ tests/` to fix.

- [ ] **Step 4: Run mypy**

Run: `mypy muc_one_up/`
Expected: Passes with no errors.

- [ ] **Step 5: Verify success criteria**

```bash
# click_main.py under 100 lines
wc -l muc_one_up/cli/click_main.py
# Expected: < 100

# No SimpleNamespace anywhere in CLI
grep -rn "SimpleNamespace" muc_one_up/cli/
# Expected: zero hits

# No argparse import
grep -rn "argparse" muc_one_up/cli/
# Expected: zero hits

# Commands registered
python -c "from muc_one_up.cli.click_main import cli; print(list(cli.commands.keys()))"
# Expected: ['simulate', 'reads', 'analyze']

# Error handler used in all command files
grep -rn "cli_error_handler" muc_one_up/cli/commands/
# Expected: hits in simulate.py, reads.py, analyze.py
```

- [ ] **Step 6: Commit any formatting fixes if needed**

```bash
git add -u
git commit -m "style: fix formatting after Phase 2B changes"
```

- [ ] **Step 7: Run full test suite one final time**

Run: `pytest --tb=short -q`
Expected: All pass. Phase 2B is complete.
