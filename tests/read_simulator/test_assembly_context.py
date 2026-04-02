"""Tests for AssemblyContext."""

from __future__ import annotations

from muc_one_up.read_simulator.assembly_context import AssemblyContext


class TestAssemblyContext:
    def test_create_from_config_hg38(self):
        config = {
            "reference_assembly": "hg38",
            "constants": {"hg38": {"left": "ATGGCCCC", "right": "CAATGGTG"}},
        }
        rs_config = {
            "sample_bam_hg38": "/path/to/sample.bam",
            "vntr_region_hg38": "chr1:155000-156000",
            "human_reference_hg38": "/path/to/hg38.fa",
        }
        ctx = AssemblyContext.from_configs(config, rs_config)
        assert ctx.assembly_name == "hg38"
        assert ctx.left_constant == "ATGGCCCC"
        assert ctx.right_constant == "CAATGGTG"
        assert ctx.sample_bam == "/path/to/sample.bam"
        assert ctx.vntr_region == "chr1:155000-156000"
        assert ctx.human_reference == "/path/to/hg38.fa"

    def test_defaults_to_hg38(self):
        config = {"constants": {"hg38": {"left": "L", "right": "R"}}}
        ctx = AssemblyContext.from_configs(config, {})
        assert ctx.assembly_name == "hg38"

    def test_hg19_assembly(self):
        config = {
            "reference_assembly": "hg19",
            "constants": {"hg19": {"left": "L19", "right": "R19"}},
        }
        rs_config = {
            "sample_bam_hg19": "/bam19",
            "vntr_region_hg19": "region19",
            "human_reference_hg19": "/ref19",
        }
        ctx = AssemblyContext.from_configs(config, rs_config)
        assert ctx.assembly_name == "hg19"
        assert ctx.left_constant == "L19"
        assert ctx.sample_bam == "/bam19"

    def test_missing_keys_return_none(self):
        config = {
            "reference_assembly": "hg38",
            "constants": {"hg38": {"left": "L", "right": "R"}},
        }
        ctx = AssemblyContext.from_configs(config, {})
        assert ctx.sample_bam is None
        assert ctx.vntr_region is None
        assert ctx.human_reference is None

    def test_none_rs_config(self):
        config = {
            "reference_assembly": "hg38",
            "constants": {"hg38": {"left": "L", "right": "R"}},
        }
        ctx = AssemblyContext.from_configs(config)
        assert ctx.assembly_name == "hg38"
        assert ctx.sample_bam is None
