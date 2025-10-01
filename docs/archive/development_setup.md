# Modern Development Setup - MucOneUp

**Last Updated:** 2025-09-30
**Stack:** uv + ruff + mypy + make + conda (hybrid approach)

---

## TL;DR

**Current:** Legacy packaging (setup.py/setup.cfg), no lock files, no tooling configuration
**Target:** Modern Python with uv, ruff, mypy, automated workflows, reproducible builds

**Hybrid Architecture:**
- **uv**: Python package management (fast, modern, with lock files)
- **conda**: Bioinformatics tools isolation (reseq, bwa, samtools, nanosim, minimap2)
- **ruff**: Linting + formatting (replaces black, flake8, isort, pyupgrade)
- **mypy**: Type checking
- **make**: Task automation

---

## Quick Start

### For New Developers

```bash
# Clone repository
git clone https://github.com/berntpopp/MucOneUp.git
cd MucOneUp

# Setup everything
make init          # Install uv, setup Python package
make conda-setup   # Create conda environments for bioinformatics tools
make check         # Run all quality checks

# Daily workflow
make test          # Run tests
make format        # Auto-format code
make lint          # Check code quality
```

---

## 1. Current State

### ❌ Issues
- **Legacy packaging**: 7-line `pyproject.toml`, using `setup.cfg` + `setup.py`
- **No lock files**: Non-reproducible builds
- **No tool configs**: ruff/mypy/black installed globally but not configured
- **No automation**: Manual commands for everything
- **Minimal dependencies**: Only 2 listed (orfipy, jsonschema)

### ✅ What Works
- Package installs via `pip install .`
- Conda environments properly isolate C/C++ bioinformatics tools
- Tests run with pytest

### Current File Structure
```
MucOneUp/
├── pyproject.toml        # 7 lines (minimal)
├── setup.cfg             # Legacy metadata
├── setup.py              # Stub
├── conda/
│   ├── env_wessim.yaml   # Illumina tools (reseq, bwa, pblat, samtools)
│   └── env_nanosim.yml   # ONT tools (nanosim, minimap2)
└── muc_one_up/
```

---

## 2. Modern Setup (Recommended)

### Architecture

```
┌─────────────────────────────────────────────┐
│  System: Python 3.10+ (via uv)             │
│  Package: muc_one_up (uv managed)          │
│  Tools: ruff, mypy (uv managed)            │
│  Lock: uv.lock (reproducible builds)       │
└─────────────────────────────────────────────┘
                    │
      ┌─────────────┴─────────────┐
      │                           │
┌─────▼──────────┐    ┌───────────▼─────────┐
│  env_wessim    │    │   env_nanosim       │
│  (Illumina)    │    │   (ONT)             │
│  - reseq       │    │  - nanosim          │
│  - bwa         │    │  - minimap2         │
│  - pblat       │    │  - samtools         │
│  - samtools    │    │  - bedtools         │
└────────────────┘    └─────────────────────┘
```

### Tool Responsibilities

| Tool | Purpose | Replaces |
|------|---------|----------|
| **uv** | Package manager, virtual env | pip, pip-tools, virtualenv |
| **ruff** | Linting + formatting | black, flake8, isort, pyupgrade |
| **mypy** | Type checking | - |
| **make** | Task automation | Manual commands |
| **pytest** | Testing | - |
| **conda/mamba** | Bioinformatics tools | - |

---

## 3. Complete Configuration Files

### `pyproject.toml` (Full Modern Configuration)

```toml
[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[project]
name = "muc_one_up"
version = "0.2.0"
description = "MucOneUp: a tool to simulate MUC1 VNTR diploid references"
readme = "README.md"
requires-python = ">=3.10"
license = {text = "MIT"}
authors = [
    {name = "Bernt Popp", email = "bernt.popp.md@gmail.com"},
]
keywords = ["bioinformatics", "genomics", "MUC1", "VNTR", "simulation"]
classifiers = [
    "Development Status :: 4 - Beta",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: MIT License",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.10",
    "Programming Language :: Python :: 3.11",
    "Programming Language :: Python :: 3.12",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

dependencies = [
    "orfipy>=0.0.3,<1.0",
    "jsonschema>=3.2.0,<5.0",
]

[project.optional-dependencies]
dev = [
    "pytest>=8.0.0",
    "pytest-cov>=5.0.0",
    "ruff>=0.8.0",
    "mypy>=1.8.0",
]

[project.scripts]
muconeup = "muc_one_up.cli:main"

[project.urls]
Homepage = "https://github.com/berntpopp/MucOneUp"
Repository = "https://github.com/berntpopp/MucOneUp"
Issues = "https://github.com/berntpopp/MucOneUp/issues"

# ==================== RUFF ====================
[tool.ruff]
target-version = "py310"
line-length = 100
indent-width = 4

# Exclude build artifacts and caches
exclude = [
    ".git",
    ".mypy_cache",
    ".pytest_cache",
    ".ruff_cache",
    "__pycache__",
    "build",
    "dist",
    "*.egg-info",
]

[tool.ruff.lint]
# Select comprehensive rule set
select = [
    "E",      # pycodestyle errors
    "W",      # pycodestyle warnings
    "F",      # Pyflakes
    "UP",     # pyupgrade
    "B",      # flake8-bugbear
    "SIM",    # flake8-simplify
    "I",      # isort
    "N",      # pep8-naming
    "C4",     # flake8-comprehensions
    "PTH",    # flake8-use-pathlib
    "RUF",    # Ruff-specific rules
]

ignore = [
    "E501",   # Line too long (handled by formatter)
    "B008",   # Do not perform function calls in argument defaults
    "N802",   # Function name should be lowercase (bioinformatics conventions)
]

# Allow autofix for all enabled rules
fixable = ["ALL"]
unfixable = []

# Allow unused variables when underscore-prefixed
dummy-variable-rgx = "^(_+|(_+[a-zA-Z0-9_]*[a-zA-Z0-9]+?))$"

[tool.ruff.lint.per-file-ignores]
"tests/*" = ["S101", "ARG", "PLR2004"]  # Allow asserts, unused args, magic values
"__init__.py" = ["F401"]  # Allow unused imports in __init__

[tool.ruff.lint.isort]
known-first-party = ["muc_one_up"]

[tool.ruff.format]
quote-style = "double"
indent-style = "space"
skip-magic-trailing-comma = false
line-ending = "auto"

# ==================== MYPY ====================
[tool.mypy]
python_version = "3.10"
warn_return_any = true
warn_unused_configs = true
warn_redundant_casts = true
warn_unused_ignores = true
warn_no_return = true
strict_equality = true
check_untyped_defs = true

# Start with gradual typing - can be made stricter over time
disallow_untyped_defs = false
disallow_incomplete_defs = false
no_implicit_optional = true

[[tool.mypy.overrides]]
module = "orfipy.*"
ignore_missing_imports = true

[[tool.mypy.overrides]]
module = "tests.*"
disallow_untyped_defs = false

# ==================== PYTEST ====================
[tool.pytest.ini_options]
testpaths = ["tests"]
python_files = ["test_*.py"]
python_classes = ["Test*"]
python_functions = ["test_*"]
addopts = [
    "-v",
    "--strict-markers",
    "--cov=muc_one_up",
    "--cov-report=term-missing",
    "--cov-report=html",
    "--cov-report=xml",
]
markers = [
    "slow: marks tests as slow (deselect with '-m \"not slow\"')",
    "integration: marks tests as integration tests",
]

# ==================== COVERAGE ====================
[tool.coverage.run]
source = ["muc_one_up"]
omit = [
    "*/tests/*",
    "*/__pycache__/*",
    "*/site-packages/*",
]

[tool.coverage.report]
exclude_lines = [
    "pragma: no cover",
    "def __repr__",
    "raise AssertionError",
    "raise NotImplementedError",
    "if __name__ == .__main__.:",
    "if TYPE_CHECKING:",
    "@abstractmethod",
]
```

### `Makefile` (Task Automation)

```makefile
.PHONY: help init install dev conda-setup test test-cov lint format check clean all

help:  ## Show this help message
	@echo "Usage: make [target]"
	@echo ""
	@echo "Targets:"
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2}'

# ==================== INSTALLATION ====================

init: install-uv dev  ## Initialize complete development environment

install-uv:  ## Install uv package manager (if not present)
	@command -v uv >/dev/null 2>&1 || { \
		echo "Installing uv..."; \
		curl -LsSf https://astral.sh/uv/install.sh | sh; \
	}
	@echo "✓ uv installed"

install:  ## Install package in production mode
	uv pip install .

dev:  ## Install package with development dependencies
	uv pip install -e ".[dev]"
	@echo "✓ Development environment ready"

conda-setup:  ## Create conda environments for bioinformatics tools
	@echo "Creating conda environments..."
	@command -v mamba >/dev/null 2>&1 || { \
		echo "ERROR: mamba not found. Install miniforge3 first."; \
		exit 1; \
	}
	mamba env create -f conda/env_wessim.yaml --force
	mamba env create -f conda/env_nanosim.yml --force
	@echo "✓ Conda environments created: env_wessim, env_nanosim"

# ==================== TESTING ====================

test:  ## Run tests with coverage
	uv run pytest

test-fast:  ## Run tests without coverage (faster)
	uv run pytest --no-cov -x

test-cov:  ## Run tests and open coverage report in browser
	uv run pytest
	@command -v xdg-open >/dev/null && xdg-open htmlcov/index.html || \
	 command -v open >/dev/null && open htmlcov/index.html || \
	 echo "Open htmlcov/index.html in your browser"

# ==================== CODE QUALITY ====================

lint:  ## Run ruff linter
	uv run ruff check muc_one_up/ tests/

lint-fix:  ## Run ruff linter and auto-fix issues
	uv run ruff check --fix muc_one_up/ tests/

format:  ## Format code with ruff
	uv run ruff format muc_one_up/ tests/

format-check:  ## Check if code is formatted (no changes)
	uv run ruff format --check muc_one_up/ tests/

type-check:  ## Run mypy type checker
	uv run mypy muc_one_up/

check: lint format-check type-check test  ## Run all quality checks

# ==================== CLEANUP ====================

clean:  ## Remove build artifacts and caches
	rm -rf build/ dist/ *.egg-info
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete
	find . -type d -name ".pytest_cache" -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name ".mypy_cache" -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name ".ruff_cache" -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name "htmlcov" -exec rm -rf {} + 2>/dev/null || true
	rm -f .coverage coverage.xml
	@echo "✓ Cleaned build artifacts"

clean-conda:  ## Remove conda environments
	@command -v mamba >/dev/null 2>&1 || { \
		echo "ERROR: mamba not found."; \
		exit 1; \
	}
	mamba env remove -n env_wessim -y || true
	mamba env remove -n env_nanosim -y || true
	@echo "✓ Conda environments removed"

# ==================== UTILITIES ====================

lock:  ## Update uv.lock file
	uv lock

sync:  ## Sync environment with uv.lock
	uv sync

update:  ## Update all dependencies
	uv lock --upgrade
	uv sync

all: clean init conda-setup check  ## Full setup from scratch

.DEFAULT_GOAL := help
```

### `.gitignore` Updates

Add to existing .gitignore:
```gitignore
# uv
.uv/
uv.lock

# Ruff
.ruff_cache/
```

---

## 4. Migration Steps

### Phase 1: Modern Packaging (Day 1)

```bash
# 1. Backup current files
cp pyproject.toml pyproject.toml.backup
cp setup.cfg setup.cfg.backup

# 2. Replace pyproject.toml with new config (from section 3)

# 3. Delete legacy files
rm setup.py setup.cfg

# 4. Install uv
curl -LsSf https://astral.sh/uv/install.sh | sh

# 5. Create lock file
uv lock

# 6. Test installation
uv pip install -e ".[dev]"
python -c "import muc_one_up; print(muc_one_up.__version__)"
```

### Phase 2: Add Makefile (Day 1)

```bash
# Copy Makefile from section 3
# Test commands
make help
make dev
make lint
make format
```

### Phase 3: Run Quality Checks (Day 2)

```bash
# Run ruff and auto-fix issues
make lint-fix

# Format all code
make format

# Run mypy (will show type errors - fix gradually)
make type-check

# Run tests
make test
```

### Phase 4: Verify Conda Integration (Day 2)

```bash
# Ensure conda environments work
make conda-setup

# Test that config.json still works
muconeup --config config.json --help
```

---

## 5. Daily Workflow

### Before Committing

```bash
# Format and fix linting issues
make format
make lint-fix

# Check types
make type-check

# Run tests
make test

# Or run everything at once
make check
```

### Adding Dependencies

```bash
# Add runtime dependency
uv add package-name

# Add dev dependency
uv add --dev package-name

# Update lock file
uv lock

# Sync environment
uv sync
```

### Running Commands in uv Environment

```bash
# Run any Python command in the project environment
uv run muconeup --config config.json

# Run Python script
uv run python scripts/my_script.py

# Run pytest
uv run pytest tests/test_specific.py
```

---

## 6. CI/CD Integration (Optional - Week 2)

### `.github/workflows/test.yml`

```yaml
name: Tests

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    strategy:
      matrix:
        python-version: ["3.10", "3.11", "3.12"]

    steps:
      - uses: actions/checkout@v4

      - name: Install uv
        uses: astral-sh/setup-uv@v5

      - name: Set up Python ${{ matrix.python-version }}
        run: uv python install ${{ matrix.python-version }}

      - name: Install dependencies
        run: uv sync --all-extras --dev

      - name: Run ruff
        run: uv run ruff check .

      - name: Run ruff format
        run: uv run ruff format --check .

      - name: Run mypy
        run: uv run mypy muc_one_up/

      - name: Run tests
        run: uv run pytest --cov --cov-report=xml

      - name: Upload coverage
        uses: codecov/codecov-action@v4
        with:
          file: ./coverage.xml
```

---

## 7. Troubleshooting

### uv not found
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
# Add to PATH: export PATH="$HOME/.cargo/bin:$PATH"
```

### mamba not found
```bash
# Install miniforge3 (includes mamba)
# macOS/Linux:
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
```

### uv.lock conflicts
```bash
rm uv.lock
uv lock
uv sync
```

### Conda environment issues
```bash
make clean-conda
make conda-setup
```

---

## 8. Quick Reference

### Common Commands

| Task | Command |
|------|---------|
| Install everything | `make init && make conda-setup` |
| Run tests | `make test` |
| Format code | `make format` |
| Check code quality | `make check` |
| Add dependency | `uv add package` |
| Update dependencies | `make update` |
| Clean build artifacts | `make clean` |

### File Locations

| File | Purpose |
|------|---------|
| `pyproject.toml` | All configuration (project, ruff, mypy, pytest) |
| `uv.lock` | Locked dependencies (commit to git) |
| `Makefile` | Task automation |
| `conda/env_*.yaml` | Bioinformatics tool environments |

### Tool Commands

```bash
# Ruff
uv run ruff check .                # Lint
uv run ruff check --fix .          # Lint and fix
uv run ruff format .               # Format

# Mypy
uv run mypy muc_one_up/            # Type check

# Pytest
uv run pytest                      # All tests
uv run pytest tests/test_file.py   # Specific file
uv run pytest -k test_name         # Specific test

# uv
uv add package                     # Add dependency
uv add --dev package               # Add dev dependency
uv remove package                  # Remove dependency
uv lock                            # Update lock file
uv sync                            # Sync environment
uv run command                     # Run in environment
```

---

## Summary

**Modern Stack:**
- ✅ **uv** for fast, reproducible Python package management
- ✅ **ruff** for comprehensive linting and formatting
- ✅ **mypy** for type safety
- ✅ **make** for developer convenience
- ✅ **conda** for bioinformatics tools (hybrid approach)

**Benefits:**
- 🚀 10-100x faster than pip
- 📦 Reproducible builds with uv.lock
- 🎯 Single pyproject.toml for all configuration
- ⚡ Fast linting and formatting with ruff
- 🔧 Simple automation with make
- 🔒 Type safety with mypy

**Next Steps:**
1. Replace pyproject.toml (5 minutes)
2. Add Makefile (2 minutes)
3. Run `make init` (1 minute)
4. Run `make check` (30 seconds)
5. Done! Start developing with modern tools.

---

**Questions?** Open an issue at https://github.com/berntpopp/MucOneUp/issues