"""Test generate_output_base from cli/outputs.py."""

from pathlib import Path

from muc_one_up.cli.outputs import generate_output_base


def test_generate_output_base_reads():
    assert (
        generate_output_base(Path("sample.001.simulated.fa"), "_reads")
        == "sample.001.simulated_reads"
    )


def test_generate_output_base_orfs():
    assert generate_output_base(Path("/path/sample.fa"), "_orfs") == "sample_orfs"


def test_generate_output_base_various():
    assert generate_output_base(Path("test.fasta"), "_stats") == "test_stats"
    assert generate_output_base(Path("x.fa"), "_ont_reads") == "x_ont_reads"
