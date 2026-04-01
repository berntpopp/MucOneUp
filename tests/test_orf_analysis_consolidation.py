"""Tests for consolidated ORF analysis."""

import inspect


def test_click_main_delegates_orf_analysis():
    """Analyze commands should delegate ORF analysis, not shell out directly."""
    import muc_one_up.cli.commands.analyze as analyze_mod

    source = inspect.getsource(analyze_mod)
    assert "subprocess" not in source, (
        "commands/analyze.py should not use subprocess to run ORF analysis."
    )
    assert "run_orf_analysis_standalone" in source, (
        "commands/analyze.py should delegate ORF analysis to analysis.run_orf_analysis_standalone()."
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
