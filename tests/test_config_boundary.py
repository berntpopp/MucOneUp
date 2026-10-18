"""Test that all CLI config reads go through the normalized load_config() loader."""

import ast
import inspect

import muc_one_up.cli.analysis as analysis_module


def test_no_raw_json_load_in_cli_commands():
    """CLI commands must not use json.load() or json.loads() to read config files.

    All config loading must go through config.load_config() or load_config_raw()
    which normalizes flat constants into the nested assembly-aware shape.
    """
    from pathlib import Path

    cli_dir = Path(__file__).resolve().parent.parent / "muc_one_up" / "cli"
    violations: list[str] = []

    for py_file in sorted(cli_dir.rglob("*.py")):
        tree = ast.parse(py_file.read_text())
        for node in ast.walk(tree):
            if isinstance(node, ast.Call):
                func = node.func
                if (
                    isinstance(func, ast.Attribute)
                    and func.attr in ("load", "loads")
                    and isinstance(func.value, ast.Name)
                    and func.value.id == "json"
                ):
                    violations.append(f"{py_file.name}:{node.lineno}")

    assert violations == [], (
        f"Found raw json.load/loads calls in CLI modules: {violations}. "
        f"Use config.load_config() or load_config_raw() instead to ensure normalization."
    )


def test_no_flat_constant_access_in_analysis():
    """analysis.py must not use config['constants']['left'] (flat format).

    After load_config() normalization, constants are nested:
    config['constants']['hg38']['left'], not config['constants']['left'].
    """
    source = inspect.getsource(analysis_module)
    assert '.get("constants", {}).get("left")' not in source, (
        "analysis.py uses flat constant access. Use assembly-keyed access instead."
    )
    assert '.get("constants", {}).get("right")' not in source, (
        "analysis.py uses flat constant access. Use assembly-keyed access instead."
    )
