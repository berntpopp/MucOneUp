"""Test that config_fingerprint module imports without rfc8785 at module level."""

import ast
import inspect


def test_config_fingerprint_module_imports_without_rfc8785_at_top_level():
    """Importing config_fingerprint should not fail even if rfc8785 has issues.

    The actual rfc8785 import should only happen inside compute_config_fingerprint().
    We verify this by checking that the module-level code does not reference rfc8785.
    """
    import muc_one_up.config_fingerprint as mod

    source = inspect.getsource(mod)
    tree = ast.parse(source)

    # Check that no top-level import statement imports rfc8785
    for node in ast.iter_child_nodes(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                assert alias.name != "rfc8785", "rfc8785 should not be imported at module level"
        elif isinstance(node, ast.ImportFrom):
            assert node.module != "rfc8785", "rfc8785 should not be imported at module level"
