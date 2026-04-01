"""Tests for SimulationOptions typed dataclass."""

from muc_one_up.cli.options import SimulationOptions


class TestSimulationOptions:
    def test_basic_construction(self):
        opts = SimulationOptions(
            config="/path/config.json",
            out_base="muc1_simulated",
            out_dir=".",
            num_haplotypes=2,
            seed=42,
            reference_assembly="hg38",
        )
        assert opts.config == "/path/config.json"
        assert opts.seed == 42

    def test_defaults(self):
        opts = SimulationOptions(config="/c.json", out_base="test", out_dir=".")
        assert opts.num_haplotypes == 2
        assert opts.seed is None
        assert opts.reference_assembly == "hg38"
        assert opts.output_structure is False
        assert opts.fixed_lengths is None
        assert opts.simulate_reads is None
        assert opts.output_orfs is False
        assert opts.orf_min_aa == 100
        assert opts.track_read_source is False

    def test_from_click_kwargs(self):
        kwargs = {
            "out_base": "test",
            "out_dir": "output",
            "num_haplotypes": 4,
            "seed": 123,
            "reference_assembly": "hg19",
            "output_structure": True,
            "fixed_lengths": ("20", "30"),
            "input_structure": None,
            "simulate_series": None,
            "mutation_name": "dupC",
            "mutation_targets": (),
            "snp_input_file": None,
            "random_snps": False,
            "random_snp_density": 0.001,
            "random_snp_output_file": None,
            "random_snp_region": "constants_only",
            "random_snp_haplotypes": "all",
            "track_read_source": True,
        }
        opts = SimulationOptions.from_click_kwargs("/config.json", kwargs)
        assert opts.config == "/config.json"
        assert opts.num_haplotypes == 4
        assert opts.fixed_lengths == ["20", "30"]
        assert opts.mutation_targets is None  # empty tuple -> None
        assert opts.track_read_source is True

    def test_from_click_kwargs_with_targets(self):
        kwargs = {
            "out_base": "t",
            "out_dir": ".",
            "num_haplotypes": 2,
            "seed": None,
            "reference_assembly": "hg38",
            "output_structure": False,
            "fixed_lengths": (),
            "input_structure": None,
            "simulate_series": None,
            "mutation_name": None,
            "mutation_targets": ("1,3", "2,5"),
            "snp_input_file": None,
            "random_snps": False,
            "random_snp_density": 0.001,
            "random_snp_output_file": None,
            "random_snp_region": "constants_only",
            "random_snp_haplotypes": "all",
            "track_read_source": False,
        }
        opts = SimulationOptions.from_click_kwargs("/c.json", kwargs)
        assert opts.mutation_targets == ["1,3", "2,5"]
        assert opts.fixed_lengths is None
