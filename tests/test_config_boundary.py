"""Test that all CLI config reads go through the normalized load_config() loader."""

import ast
import inspect

import muc_one_up.cli.analysis as analysis_module
import muc_one_up.cli.click_main as cli_module


def test_no_raw_json_load_in_cli_commands():
    """CLI commands must not use json.load() to read config files.

    All config loading must go through config.load_config() which normalizes
    flat constants into the nested assembly-aware shape.
    """
    source = inspect.getsource(cli_module)
    tree = ast.parse(source)

    violations = []
    for node in ast.walk(tree):
        # Look for calls like json.load(f)
        if isinstance(node, ast.Call):
            func = node.func
            if (
                isinstance(func, ast.Attribute)
                and func.attr == "load"
                and isinstance(func.value, ast.Name)
                and func.value.id == "json"
            ):
                violations.append(node.lineno)

    assert violations == [], (
        f"Found json.load() calls at lines {violations} in click_main.py. "
        f"Use config.load_config() instead to ensure constants normalization."
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
