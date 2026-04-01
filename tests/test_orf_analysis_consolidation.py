"""Tests for consolidated ORF analysis."""

import inspect


def test_no_subprocess_orfipy_in_click_main():
    """click_main.py should not call orfipy directly via subprocess."""
    import muc_one_up.cli.click_main as cli_mod

    source = inspect.getsource(cli_mod)
    assert "orfipy" not in source, (
        "click_main.py should not reference orfipy. Use analysis.run_orf_analysis_standalone()."
    )


def test_run_orf_analysis_standalone_exists():
    """analysis.py should have run_orf_analysis_standalone."""
    from muc_one_up.cli.analysis import run_orf_analysis_standalone

    sig = inspect.signature(run_orf_analysis_standalone)
    param_names = list(sig.parameters.keys())
    assert "input_fasta" in param_names
    assert "out_dir" in param_names
    assert "out_base" in param_names
    assert "orf_min_aa" in param_names
    assert "orf_aa_prefix" in param_names
    assert "left_const" in param_names
    assert "right_const" in param_names
