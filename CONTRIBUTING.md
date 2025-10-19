# Contributing to MucOneUp

## Before Committing

**IMPORTANT**: Always run CI checks locally before committing to avoid CI failures:

```bash
make ci-check
```

This runs the **exact same checks** as GitHub Actions CI:
1. `ruff check muc_one_up/ tests/` - Linting (RUF059, etc.)
2. `ruff format --check muc_one_up/ tests/` - Formatting check
3. `mypy muc_one_up/` - Type checking

If formatting issues are found, fix them with:
```bash
make format
```

Then re-run `make ci-check` to verify.

## Development Workflow

1. **Make changes** to code/tests
2. **Run tests**: `make test` or `pytest`
3. **Run CI checks**: `make ci-check` ⚠️ **DO NOT SKIP THIS**
4. **Fix any issues** found by CI checks
5. **Commit** your changes
6. **Push** to GitHub

## Pre-commit Hooks

The repository uses pre-commit hooks that run automatically on commit. These hooks will:
- Fix safe issues automatically (like line endings)
- Apply ruff formatting
- Run linting checks

If pre-commit modifies files, you need to re-stage and commit again.

## Common Commands

```bash
# Install dev environment
make dev

# Run tests
make test

# Run CI checks (BEFORE committing!)
make ci-check

# Format code
make format

# Run linter with auto-fix
make lint-fix

# Type check
make type-check

# Run all checks + tests
make check
```

## CI Alignment

The `make ci-check` command is designed to **exactly match** the GitHub Actions CI workflow defined in `.github/workflows/test.yml`. If `make ci-check` passes locally, CI will pass on GitHub.

## Troubleshooting

### "Would reformat" errors in CI
Run `make format` locally before committing.

### RUF059 violations (unused unpacked variables)
Prefix unused variables with `_`:
```python
# Bad
results, info = some_function()  # if info is unused

# Good
results, _info = some_function()
```

### Line ending issues
Pre-commit hooks handle this automatically. Just re-commit after the hook fixes them.
