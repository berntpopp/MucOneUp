"""Tests for muc_one_up.cli.orchestration module.

Covers: run_single_simulation_iteration — call ordering, data flow,
dual mutation mode, source tracker construction.
"""

from unittest.mock import patch

import pytest

from muc_one_up.type_defs import HaplotypeResult, MutationTarget, RepeatUnit

# All delegates are module-level imports in orchestration.py, so we patch there.
PATCH_PREFIX = "muc_one_up.cli.orchestration"


def _make_haplotype_results():
    """Build minimal HaplotypeResult objects for mock returns."""
    return [
        HaplotypeResult(
            sequence="AAAAATTTTT",
            chain=[RepeatUnit("1"), RepeatUnit("9")],
        ),
        HaplotypeResult(
            sequence="CCCCCGGGGG",
            chain=[RepeatUnit("1"), RepeatUnit("9")],
        ),
    ]


@pytest.fixture
def _mock_delegates():
    """Patch all 6 delegate functions called by run_single_simulation_iteration.

    Returns a dict of mock objects keyed by function name.
    """
    results = _make_haplotype_results()

    with (
        patch(f"{PATCH_PREFIX}.generate_haplotypes", return_value=results) as m_hap,
        patch(
            f"{PATCH_PREFIX}.apply_mutation_pipeline",
            return_value=(results, None, None, None),
        ) as m_mut,
        patch(
            f"{PATCH_PREFIX}.write_fasta_outputs",
            return_value=(results, None, None, None),
        ) as m_fasta,
        patch(f"{PATCH_PREFIX}.write_mutated_units") as m_units,
        patch(f"{PATCH_PREFIX}.write_structure_files") as m_struct,
        patch(f"{PATCH_PREFIX}.write_simulation_statistics") as m_stats,
    ):
        yield {
            "generate_haplotypes": m_hap,
            "apply_mutation_pipeline": m_mut,
            "write_fasta_outputs": m_fasta,
            "write_mutated_units": m_units,
            "write_structure_files": m_struct,
            "write_simulation_statistics": m_stats,
        }


class TestRunSingleSimulationIteration:
    """Tests for the main orchestration function."""

    def test_calls_all_steps_no_mutation(
        self, mock_simulation_options, minimal_config, _mock_delegates
    ):
        """All delegate functions called in correct order with no mutation."""
        from muc_one_up.cli.orchestration import run_single_simulation_iteration

        run_single_simulation_iteration(
            args=mock_simulation_options,
            config=minimal_config,
            out_dir=mock_simulation_options.out_dir,
            out_base="test",
            sim_index=1,
            fixed_conf=None,
            predefined_chains=None,
            dual_mutation_mode=False,
            mutation_pair=None,
            structure_mutation_info=None,
        )

        # All delegates called exactly once
        for name, mock in _mock_delegates.items():
            assert mock.called, f"{name} was not called"

    def test_results_flow_through_pipeline(
        self, mock_simulation_options, minimal_config, _mock_delegates
    ):
        """Return values from each step are passed as arguments to the next."""
        from muc_one_up.cli.orchestration import run_single_simulation_iteration

        results = _make_haplotype_results()
        mutated_results = _make_haplotype_results()
        mutated_units = {1: [(3, "ATCG")]}
        mutation_positions = [MutationTarget(1, 3)]
        snp_normal = [None, None]
        snp_mut = [None, None]

        _mock_delegates["generate_haplotypes"].return_value = results
        _mock_delegates["apply_mutation_pipeline"].return_value = (
            results,
            mutated_results,
            mutated_units,
            mutation_positions,
        )
        _mock_delegates["write_fasta_outputs"].return_value = (
            results,
            mutated_results,
            snp_normal,
            snp_mut,
        )

        mock_simulation_options.mutation_name = "dupC"

        run_single_simulation_iteration(
            args=mock_simulation_options,
            config=minimal_config,
            out_dir=mock_simulation_options.out_dir,
            out_base="test",
            sim_index=1,
            fixed_conf=None,
            predefined_chains=None,
            dual_mutation_mode=False,
            mutation_pair=None,
            structure_mutation_info=None,
        )

        # apply_mutation_pipeline receives results from generate_haplotypes
        mut_call = _mock_delegates["apply_mutation_pipeline"].call_args
        assert mut_call[0][2] is results  # third positional arg is results

        # write_fasta_outputs receives mutation_positions from apply_mutation_pipeline
        fasta_call = _mock_delegates["write_fasta_outputs"].call_args
        assert fasta_call[0][9] is mutation_positions  # 10th positional arg

        # write_structure_files receives mutation_positions
        struct_call = _mock_delegates["write_structure_files"].call_args
        assert struct_call[0][8] is mutation_positions

    def test_dual_mutation_mode(self, mock_simulation_options, minimal_config, _mock_delegates):
        """Dual mode passes mutation_pair through to delegates."""
        from muc_one_up.cli.orchestration import run_single_simulation_iteration

        mutation_pair = ("normal", "dupC")
        mock_simulation_options.mutation_name = "normal,dupC"

        run_single_simulation_iteration(
            args=mock_simulation_options,
            config=minimal_config,
            out_dir=mock_simulation_options.out_dir,
            out_base="test",
            sim_index=1,
            fixed_conf=None,
            predefined_chains=None,
            dual_mutation_mode=True,
            mutation_pair=mutation_pair,
            structure_mutation_info=None,
        )

        # write_fasta_outputs should receive dual_mutation_mode=True
        fasta_call = _mock_delegates["write_fasta_outputs"].call_args
        assert fasta_call[0][7] is True  # dual_mutation_mode

        # write_simulation_statistics should receive mutation_pair
        stats_call = _mock_delegates["write_simulation_statistics"].call_args
        assert stats_call[0][10] is mutation_pair

    def test_rng_threaded_to_delegates(
        self, mock_simulation_options, minimal_config, _mock_delegates
    ):
        """rng is passed to generate_haplotypes, apply_mutation_pipeline, and write_fasta_outputs."""
        import random

        from muc_one_up.cli.orchestration import run_single_simulation_iteration

        test_rng = random.Random(12345)

        run_single_simulation_iteration(
            args=mock_simulation_options,
            config=minimal_config,
            out_dir=mock_simulation_options.out_dir,
            out_base="test",
            sim_index=1,
            fixed_conf=None,
            predefined_chains=None,
            dual_mutation_mode=False,
            mutation_pair=None,
            structure_mutation_info=None,
            rng=test_rng,
        )

        assert _mock_delegates["generate_haplotypes"].call_args[1]["rng"] is test_rng
        assert _mock_delegates["apply_mutation_pipeline"].call_args[1]["rng"] is test_rng
        assert _mock_delegates["write_fasta_outputs"].call_args[1]["rng"] is test_rng

    def test_mutated_units_passed_to_statistics(
        self, mock_simulation_options, minimal_config, _mock_delegates
    ):
        """mutated_units from apply_mutation_pipeline is forwarded to write_simulation_statistics."""
        from muc_one_up.cli.orchestration import run_single_simulation_iteration

        results = _make_haplotype_results()
        mutated_units = {1: [(2, "ATCG_MUT")]}
        _mock_delegates["apply_mutation_pipeline"].return_value = (
            results,
            results,
            mutated_units,
            [MutationTarget(1, 2)],
        )

        run_single_simulation_iteration(
            args=mock_simulation_options,
            config=minimal_config,
            out_dir=mock_simulation_options.out_dir,
            out_base="test",
            sim_index=1,
            fixed_conf=None,
            predefined_chains=None,
            dual_mutation_mode=False,
            mutation_pair=None,
            structure_mutation_info=None,
        )

        stats_call = _mock_delegates["write_simulation_statistics"].call_args
        assert stats_call[1]["mutated_units"] == mutated_units
        assert stats_call[1]["mutation_positions"] == [MutationTarget(1, 2)]

    def test_statistics_receive_timing(
        self, mock_simulation_options, minimal_config, _mock_delegates
    ):
        """write_simulation_statistics receives valid start < end timing."""
        from muc_one_up.cli.orchestration import run_single_simulation_iteration

        run_single_simulation_iteration(
            args=mock_simulation_options,
            config=minimal_config,
            out_dir=mock_simulation_options.out_dir,
            out_base="test",
            sim_index=1,
            fixed_conf=None,
            predefined_chains=None,
            dual_mutation_mode=False,
            mutation_pair=None,
            structure_mutation_info=None,
        )

        stats_call = _mock_delegates["write_simulation_statistics"].call_args
        iteration_start = stats_call[0][5]
        iteration_end = stats_call[0][6]
        assert isinstance(iteration_start, float)
        assert isinstance(iteration_end, float)
        assert iteration_start <= iteration_end


class TestSeededSeriesReproducibility:
    """Tests for seed reproducibility across series iterations (#139)."""

    def test_seeded_series_iterations_reproducible_across_runs(
        self, mock_simulation_options, minimal_config
    ):
        """Sharing RNG across iterations produces different haplotypes per iteration,
        but repeating the run with the same seed reproduces the exact sequence series.
        """
        import random

        from muc_one_up.cli.haplotypes import generate_haplotypes

        mock_simulation_options.num_haplotypes = 2
        mock_simulation_options.seed = 42

        # Run 1: seed=42, 2 iterations
        rng1 = random.Random(42)
        iter1_run1 = generate_haplotypes(
            mock_simulation_options, minimal_config, None, None, rng=rng1
        )
        iter2_run1 = generate_haplotypes(
            mock_simulation_options, minimal_config, None, None, rng=rng1
        )

        # Iterations within the same run must differ (RNG advances)
        assert iter1_run1[0].sequence != iter2_run1[0].sequence

        # Run 2: same seed=42, 2 iterations
        rng2 = random.Random(42)
        iter1_run2 = generate_haplotypes(
            mock_simulation_options, minimal_config, None, None, rng=rng2
        )
        iter2_run2 = generate_haplotypes(
            mock_simulation_options, minimal_config, None, None, rng=rng2
        )

        # Matching iterations across identical-seed runs must be identical
        assert iter1_run1[0].sequence == iter1_run2[0].sequence
        assert iter2_run1[0].sequence == iter2_run2[0].sequence
        assert iter1_run1[1].sequence == iter1_run2[1].sequence
        assert iter2_run1[1].sequence == iter2_run2[1].sequence
