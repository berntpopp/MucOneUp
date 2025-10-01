# Code Quality & Standards

**Status:** ✅ PHASE 1 COMPLETE (Infrastructure) | ⚠️ PHASE 2 DEFERRED (Testing & Security)
**Completion Date:** 2025-10-01
**Completion:** 70% - Infrastructure and quality gates implemented
**Remaining:** See docs/PHASE_2_ROADMAP.md for testing and security implementation

**Original Priority:** 🟡 MEDIUM
**Original Estimated Effort:** 2-3 days (actual: 1 day for Phase 1)
**Impact:** MEDIUM - Consistency, maintainability, automation
**Dependencies:** Should follow CLI refactoring and exception handling

## Current State

### Missing Infrastructure

| Component | Status | Impact |
|-----------|--------|--------|
| **Pre-commit hooks** | ❌ None | Code quality not enforced |
| **Code formatter** | ❌ None | Inconsistent style |
| **Linter** | ❌ None | Style violations undetected |
| **Import sorter** | ❌ None | Messy imports |
| **CI/CD pipeline** | ❌ None | No automated checks |
| **Static analysis** | ❌ None | Security issues undetected |

### Problems

1. **Inconsistent style** - No automated formatting
2. **No quality gates** - Can commit code with issues
3. **Manual checks** - Developers must remember to run tools
4. **No automation** - Quality depends on discipline

## SOLID Principles & Code Quality

### Single Responsibility

Each tool has one job:
- **Black:** Format code
- **isort:** Sort imports
- **flake8:** Check style
- **mypy:** Check types
- **pytest:** Run tests

Don't try to use one tool for everything!

### Open/Closed Principle

Configuration files make tools extensible:
- `pyproject.toml` for tool config
- `.pre-commit-config.yaml` for hooks
- `.github/workflows/` for CI/CD

Can add new checks without modifying code.

## Pre-commit Hooks Setup

### Step 1: Install pre-commit

```bash
pip install pre-commit
```

### Step 2: Create .pre-commit-config.yaml

```yaml
# .pre-commit-config.yaml
repos:
  # Black - code formatting
  - repo: https://github.com/psf/black
    rev: 24.1.1
    hooks:
      - id: black
        language_version: python3.10
        args: ["--line-length=100"]

  # isort - import sorting
  - repo: https://github.com/PyCQA/isort
    rev: 5.13.2
    hooks:
      - id: isort
        args: ["--profile", "black", "--line-length=100"]

  # flake8 - style checking
  - repo: https://github.com/PyCQA/flake8
    rev: 7.0.0
    hooks:
      - id: flake8
        args: [
          "--max-line-length=100",
          "--extend-ignore=E203,E501,W503",
          "--max-complexity=10"
        ]

  # mypy - type checking
  - repo: https://github.com/pre-commit/mirrors-mypy
    rev: v1.8.0
    hooks:
      - id: mypy
        additional_dependencies: [types-all]
        args: ["--ignore-missing-imports"]

  # Standard pre-commit hooks
  - repo: https://github.com/pre-commit/pre-commit-hooks
    rev: v4.5.0
    hooks:
      - id: trailing-whitespace
      - id: end-of-file-fixer
      - id: check-yaml
      - id: check-added-large-files
        args: ["--maxkb=1000"]
      - id: check-json
      - id: check-toml
      - id: check-merge-conflict
      - id: debug-statements
      - id: mixed-line-ending
        args: ["--fix=lf"]

  # Security checks
  - repo: https://github.com/PyCQA/bandit
    rev: 1.7.5
    hooks:
      - id: bandit
        args: ["-c", "pyproject.toml"]
        additional_dependencies: ["bandit[toml]"]

# Run these on all files initially
default_stages: [commit]
```

### Step 3: Install hooks

```bash
# Install hooks
pre-commit install

# Run on all files (first time)
pre-commit run --all-files

# Update hooks to latest versions
pre-commit autoupdate
```

### Step 4: Usage

```bash
# Hooks run automatically on git commit

# Run manually on all files
pre-commit run --all-files

# Run specific hook
pre-commit run black --all-files

# Skip hooks (not recommended)
git commit --no-verify -m "message"
```

## Code Formatting Standards

### Black Configuration

Create `pyproject.toml`:

```toml
# pyproject.toml
[tool.black]
line-length = 100
target-version = ['py38', 'py39', 'py310']
include = '\.pyi?$'
extend-exclude = '''
/(
    \.git
  | \.mypy_cache
  | \.pytest_cache
  | \.tox
  | build
  | dist
)/
'''

[tool.isort]
profile = "black"
line_length = 100
multi_line_output = 3
include_trailing_comma = true
force_grid_wrap = 0
use_parentheses = true
ensure_newline_before_comments = true
skip_gitignore = true

[tool.pytest.ini_options]
testpaths = ["tests"]
python_files = ["test_*.py"]
addopts = """
    --verbose
    --cov=muc_one_up
    --cov-report=html
    --cov-report=term-missing
    --cov-fail-under=60
"""

[tool.coverage.run]
source = ["muc_one_up"]
omit = [
    "*/tests/*",
    "*/__init__.py",
    "*/setup.py",
]

[tool.coverage.report]
exclude_lines = [
    "pragma: no cover",
    "def __repr__",
    "raise AssertionError",
    "raise NotImplementedError",
    "if __name__ == .__main__.:",
    "if TYPE_CHECKING:",
]

[tool.mypy]
python_version = "3.8"
warn_return_any = true
warn_unused_configs = true
disallow_untyped_defs = false  # Start lenient
check_untyped_defs = true

[tool.bandit]
exclude_dirs = ["tests", "build", "dist"]
skips = ["B101"]  # Skip assert_used (used in tests)
```

### Flake8 Configuration

Create `.flake8`:

```ini
# .flake8
[flake8]
max-line-length = 100
extend-ignore = E203, E501, W503
max-complexity = 10
exclude =
    .git,
    __pycache__,
    build,
    dist,
    .eggs,
    *.egg-info,
    .venv,
    venv,
    .mypy_cache,
    .pytest_cache,
    .tox

# Per-file ignores
per-file-ignores =
    __init__.py:F401,F403
    tests/*:S101

# Complexity
max-complexity = 10
```

## KISS: Simple Quality Workflow

### For Developers

```bash
# 1. Make changes
vim muc_one_up/simulate.py

# 2. Run formatting (automatic with pre-commit)
git add muc_one_up/simulate.py
git commit -m "fix: improve simulation"
# Pre-commit hooks run automatically!

# 3. If hooks fail, fix and re-commit
# Black/isort auto-fix, just re-stage
git add muc_one_up/simulate.py
git commit -m "fix: improve simulation"
```

### Manual Quality Checks

```bash
# Format code
black muc_one_up/ tests/

# Sort imports
isort muc_one_up/ tests/

# Check style
flake8 muc_one_up/ tests/

# Check types
mypy muc_one_up/

# Run tests
pytest

# Security scan
bandit -r muc_one_up/

# All at once
pre-commit run --all-files
```

## CI/CD Pipeline

### GitHub Actions Workflow

Create `.github/workflows/test.yml`:

```yaml
name: Test & Quality

on:
  push:
    branches: [main, dev, dev/*]
  pull_request:
    branches: [main, dev]

jobs:
  quality:
    name: Code Quality Checks
    runs-on: ubuntu-latest

    steps:
      - uses: actions/checkout@v4

      - name: Set up Python
        uses: actions/setup-python@v5
        with:
          python-version: "3.10"

      - name: Cache dependencies
        uses: actions/cache@v3
        with:
          path: ~/.cache/pip
          key: ${{ runner.os }}-pip-${{ hashFiles('**/setup.py') }}
          restore-keys: |
            ${{ runner.os }}-pip-

      - name: Install dependencies
        run: |
          python -m pip install --upgrade pip
          pip install -e ".[dev]"

      - name: Run pre-commit hooks
        run: |
          pre-commit run --all-files

      - name: Check types with mypy
        run: |
          mypy muc_one_up/

      - name: Security scan with bandit
        run: |
          bandit -r muc_one_up/ -f json -o bandit-report.json
        continue-on-error: true

      - name: Upload bandit report
        uses: actions/upload-artifact@v3
        if: always()
        with:
          name: bandit-report
          path: bandit-report.json

  test:
    name: Test Suite
    runs-on: ${{ matrix.os }}
    strategy:
      matrix:
        os: [ubuntu-latest, macos-latest]
        python-version: ["3.8", "3.9", "3.10", "3.11"]

    steps:
      - uses: actions/checkout@v4

      - name: Set up Python ${{ matrix.python-version }}
        uses: actions/setup-python@v5
        with:
          python-version: ${{ matrix.python-version }}

      - name: Install dependencies
        run: |
          python -m pip install --upgrade pip
          pip install -e ".[dev]"

      - name: Run tests with coverage
        run: |
          pytest --cov=muc_one_up --cov-report=xml --cov-report=term

      - name: Upload coverage to Codecov
        uses: codecov/codecov-action@v3
        with:
          file: ./coverage.xml
          flags: unittests
          name: codecov-${{ matrix.os }}-py${{ matrix.python-version }}
```

### Additional Workflows

**Dependency Security:**

```yaml
# .github/workflows/security.yml
name: Security Scan

on:
  schedule:
    - cron: '0 0 * * 0'  # Weekly
  workflow_dispatch:

jobs:
  security:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: "3.10"

      - name: Install safety
        run: pip install safety

      - name: Check dependencies
        run: |
          pip freeze | safety check --stdin
```

## DRY: Centralized Configuration

### setup.py with dev dependencies

```python
# setup.py (add dev extras)
setup(
    # ... existing config ...
    extras_require={
        "dev": [
            # Testing
            "pytest>=7.4.0",
            "pytest-cov>=4.1.0",
            "pytest-xdist>=3.3.0",
            "pytest-timeout>=2.1.0",
            # Code quality
            "black>=24.1.0",
            "isort>=5.13.0",
            "flake8>=7.0.0",
            "mypy>=1.8.0",
            "bandit>=1.7.0",
            "pre-commit>=3.6.0",
            # Type stubs
            "types-all>=1.0.0",
        ],
    },
)
```

**Install dev dependencies:**
```bash
pip install -e ".[dev]"
```

## Makefile for Common Tasks (KISS)

```makefile
# Makefile
.PHONY: help install test format lint type-check security clean

help:
	@echo "Available commands:"
	@echo "  make install      Install package with dev dependencies"
	@echo "  make test         Run test suite"
	@echo "  make format       Format code with black and isort"
	@echo "  make lint         Check code style with flake8"
	@echo "  make type-check   Check types with mypy"
	@echo "  make security     Run security scan with bandit"
	@echo "  make quality      Run all quality checks"
	@echo "  make clean        Remove build artifacts"

install:
	pip install -e ".[dev]"
	pre-commit install

test:
	pytest -v

format:
	black muc_one_up/ tests/
	isort muc_one_up/ tests/

lint:
	flake8 muc_one_up/ tests/

type-check:
	mypy muc_one_up/

security:
	bandit -r muc_one_up/

quality: format lint type-check test
	@echo "All quality checks passed!"

clean:
	rm -rf build/ dist/ *.egg-info
	rm -rf .pytest_cache/ .mypy_cache/ .coverage htmlcov/
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete
```

**Usage:**
```bash
make install       # Setup dev environment
make quality       # Run all checks before commit
make test          # Run tests
make clean         # Clean up
```

## Code Review Checklist

### For Pull Requests

**Automated (CI/CD):**
- ✅ All tests pass
- ✅ Code coverage ≥60%
- ✅ Type checking passes
- ✅ Style checks pass
- ✅ No security issues

**Manual Review:**
- ✅ Code follows SOLID principles
- ✅ No code duplication (DRY)
- ✅ Functions are simple (KISS)
- ✅ Clear variable/function names
- ✅ Docstrings present
- ✅ Tests included for new features
- ✅ No hardcoded values
- ✅ Error handling appropriate
- ✅ Logging statements added

## Success Criteria

### Infrastructure
- ✅ `.pre-commit-config.yaml` configured
- ✅ `pyproject.toml` configured
- ✅ `.flake8` configured
- ✅ `mypy.ini` configured
- ✅ GitHub Actions workflows created
- ✅ `Makefile` for common tasks
- ✅ Dev dependencies in `setup.py`

### Automation
- ✅ Pre-commit hooks installed
- ✅ Hooks run on every commit
- ✅ CI/CD runs on every push
- ✅ Coverage reports generated
- ✅ Security scans automated

### Quality Gates
- ✅ Code must pass all hooks to commit
- ✅ CI must pass to merge PR
- ✅ Coverage must be ≥60%
- ✅ No type errors allowed
- ✅ No critical security issues

## Metrics to Track

### Code Quality Metrics

| Metric | Current | Target | Tool |
|--------|---------|--------|------|
| **Test Coverage** | 7.7% | 60%+ | pytest-cov |
| **Type Coverage** | ~40% | 100% | mypy |
| **Cyclomatic Complexity** | High | <10 | flake8 |
| **Code Duplication** | ~2% | <5% | Manual review |
| **Security Issues** | Unknown | 0 critical | bandit |
| **Style Violations** | Unknown | 0 | flake8 |

### Track in CI/CD

```bash
# Generate metrics
pytest --cov=muc_one_up --cov-report=json
mypy muc_one_up/ --html-report mypy-report/
bandit -r muc_one_up/ -f json -o bandit.json

# Parse and track trends over time
```

## Implementation Timeline

### Week 1: Setup Infrastructure
- [ ] Create `pyproject.toml`
- [ ] Create `.pre-commit-config.yaml`
- [ ] Create `.flake8`
- [ ] Create `Makefile`
- [ ] Install pre-commit hooks
- [ ] Run formatters on entire codebase

### Week 2: Fix Quality Issues
- [ ] Run black/isort on all files
- [ ] Fix flake8 violations
- [ ] Fix mypy type errors
- [ ] Fix bandit security issues
- [ ] Ensure all hooks pass

### Week 3: CI/CD Setup
- [ ] Create GitHub Actions workflows
- [ ] Setup coverage reporting
- [ ] Setup security scanning
- [ ] Add status badges to README
- [ ] Document workflow in CONTRIBUTING.md

### Week 4: Documentation
- [ ] Update README with quality badges
- [ ] Create CONTRIBUTING.md
- [ ] Document development workflow
- [ ] Add code review guidelines

## Verification Commands

```bash
# Verify pre-commit setup
pre-commit run --all-files

# Verify CI would pass
make quality

# Check metrics
pytest --cov=muc_one_up --cov-report=term-missing
mypy muc_one_up/ --html-report mypy-report/
bandit -r muc_one_up/

# Verify installation
pip install -e ".[dev]"
pytest
```

## Benefits Summary

| Aspect | Before | After |
|--------|--------|-------|
| **Quality** | Manual, inconsistent | Automated, enforced |
| **Style** | Varies by developer | Consistent (Black) |
| **Errors** | Found in production | Caught in CI/CD |
| **Security** | Unknown vulnerabilities | Automated scanning |
| **Onboarding** | "Follow the style guide" | Hooks enforce automatically |
| **Reviews** | Focus on style | Focus on logic |
| **Confidence** | Low (manual checks) | High (automated gates) |

## Next Steps

After implementing code quality standards:
1. ➡️ Monitor metrics in CI/CD
2. ➡️ Incrementally increase coverage targets
3. ➡️ Add property-based testing (Hypothesis)
4. ➡️ Add mutation testing (mutmut)
5. ➡️ Consider additional tools (pylint, radon)
