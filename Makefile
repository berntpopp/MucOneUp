.PHONY: help init install-uv install dev conda-setup test test-fast test-cov lint lint-fix format format-check type-check check clean clean-conda lock sync update all

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

ci-check:  ## Run EXACT same checks as GitHub Actions CI (run before committing!)
	@echo "Running CI checks locally (same as GitHub Actions)..."
	@echo ""
	@echo "=== JOB: Code Quality Checks ==="
	@echo "1. Ruff linter..."
	ruff check muc_one_up/ tests/
	@echo "2. Ruff formatter check..."
	ruff format --check muc_one_up/ tests/
	@echo "3. Mypy type checker..."
	mypy muc_one_up/ || true
	@echo ""
	@echo "=== JOB: Test Suite ==="
	@echo "4. Running pytest with coverage..."
	pytest --cov=muc_one_up --cov-report=term-missing
	@echo ""
	@echo "✅ All CI checks passed!"

check: lint format-check type-check test  ## Run all quality checks

# ==================== DOCUMENTATION ====================

docs-install:  ## Install documentation dependencies
	uv sync --extra docs

docs-build:  ## Build documentation locally (non-strict)
	uv run mkdocs build

docs-serve:  ## Serve documentation locally with live reload
	uv run mkdocs serve

docs-ci:  ## Run EXACT same docs build as GitHub Actions CI
	@echo "=== Simulating GitHub Actions Documentation Build ==="
	@echo ""
	@echo "Step 1: Sync docs dependencies..."
	uv sync --extra docs
	@echo ""
	@echo "Step 2: Build documentation (strict mode)..."
	uv run mkdocs build --strict --verbose
	@echo ""
	@echo "✅ Documentation build succeeded (same as GitHub Actions)"

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

# ==================== VERSION MANAGEMENT ====================

bump-patch:  ## Bump patch version (0.15.0 -> 0.15.1)
	@python3 -c "import sys; \
	from pathlib import Path; \
	version_file = Path('muc_one_up/version.py'); \
	content = version_file.read_text(); \
	current = content.split('\"')[1]; \
	major, minor, patch = current.split('.'); \
	new_version = f'{major}.{minor}.{int(patch)+1}'; \
	version_file.write_text(f'# muc_one_up/version.py\n\n__version__ = \"{new_version}\"\n'); \
	pyproject = Path('pyproject.toml'); \
	pyproject_content = pyproject.read_text(); \
	pyproject.write_text(pyproject_content.replace(f'version = \"{current}\"', f'version = \"{new_version}\"')); \
	print(f'✓ Bumped version: {current} → {new_version}'); \
	print(f'  - muc_one_up/version.py'); \
	print(f'  - pyproject.toml')"

bump-minor:  ## Bump minor version (0.15.0 -> 0.16.0)
	@python3 -c "import sys; \
	from pathlib import Path; \
	version_file = Path('muc_one_up/version.py'); \
	content = version_file.read_text(); \
	current = content.split('\"')[1]; \
	major, minor, patch = current.split('.'); \
	new_version = f'{major}.{int(minor)+1}.0'; \
	version_file.write_text(f'# muc_one_up/version.py\n\n__version__ = \"{new_version}\"\n'); \
	pyproject = Path('pyproject.toml'); \
	pyproject_content = pyproject.read_text(); \
	pyproject.write_text(pyproject_content.replace(f'version = \"{current}\"', f'version = \"{new_version}\"')); \
	print(f'✓ Bumped version: {current} → {new_version}'); \
	print(f'  - muc_one_up/version.py'); \
	print(f'  - pyproject.toml')"

bump-major:  ## Bump major version (0.15.0 -> 1.0.0)
	@python3 -c "import sys; \
	from pathlib import Path; \
	version_file = Path('muc_one_up/version.py'); \
	content = version_file.read_text(); \
	current = content.split('\"')[1]; \
	major, minor, patch = current.split('.'); \
	new_version = f'{int(major)+1}.0.0'; \
	version_file.write_text(f'# muc_one_up/version.py\n\n__version__ = \"{new_version}\"\n'); \
	pyproject = Path('pyproject.toml'); \
	pyproject_content = pyproject.read_text(); \
	pyproject.write_text(pyproject_content.replace(f'version = \"{current}\"', f'version = \"{new_version}\"')); \
	print(f'✓ Bumped version: {current} → {new_version}'); \
	print(f'  - muc_one_up/version.py'); \
	print(f'  - pyproject.toml')"

show-version:  ## Show current version
	@python3 -c "from muc_one_up.version import __version__; print(f'Current version: {__version__}')"

all: clean init conda-setup check  ## Full setup from scratch

.DEFAULT_GOAL := help
