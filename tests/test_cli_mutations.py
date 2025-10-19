"""Comprehensive tests for muc_one_up.cli.mutations module.

Following modern testing best practices:
- Parametrized tests for multiple scenarios
- Descriptive test names
- Edge case and error condition testing
- Proper use of fixtures
- AAA pattern (Arrange-Act-Assert)
"""

from unittest.mock import Mock, patch

import pytest

from muc_one_up.cli.mutations import apply_mutation_pipeline
from muc_one_up.exceptions import MutationError


@pytest.mark.unit
class TestApplyMutationPipeline:
    """Comprehensive tests for apply_mutation_pipeline function."""

    def test_no_mutation_returns_unchanged_results(self, minimal_config: dict):
        """Given no mutation name, when applying, then returns unchanged results."""
        args = Mock()
        args.mutation_targets = None
        results = [("ATCG", ["1", "2"])]

        res, mut_res, mut_units, mut_pos = apply_mutation_pipeline(
            args=args,
            config=minimal_config,
            results=results,
            mutation_name=None,
            dual_mutation_mode=False,
            mutation_pair=None,
        )

        assert res == results
        assert mut_res is None
        assert mut_units is None
        assert mut_pos is None

    @patch("muc_one_up.cli.mutations.apply_mutations")
    def test_single_mode_with_explicit_targets(self, mock_apply: Mock, minimal_config: dict):
        """Given explicit targets in single mode, when applying, then mutates at targets."""
        args = Mock()
        args.mutation_targets = ["1,1"]

        left = minimal_config["constants"]["hg38"]["left"]
        repeat_seq = minimal_config["repeats"]["X"]
        right = minimal_config["constants"]["hg38"]["right"]
        results = [(f"{left}{repeat_seq}{right}", ["X"])]

        # Mock the apply_mutations to return expected structure
        mock_apply.return_value = (results, {"mutated": "units"})

        _res, mut_res, mut_units, mut_pos = apply_mutation_pipeline(
            args=args,
            config=minimal_config,
            results=results,
            mutation_name="dupC",
            dual_mutation_mode=False,
            mutation_pair=None,
        )

        assert mut_pos == [(1, 1)]
        assert mut_res is None  # Single mode doesn't create separate mutated results
        assert mut_units is not None
        mock_apply.assert_called_once()

    @patch("muc_one_up.cli.mutations.apply_mutations")
    def test_single_mode_with_random_target(self, mock_apply: Mock, minimal_config: dict):
        """Given no explicit targets in single mode, when applying, then finds random target."""
        args = Mock()
        args.mutation_targets = None

        left = minimal_config["constants"]["hg38"]["left"]
        repeat_seq = minimal_config["repeats"]["X"]
        right = minimal_config["constants"]["hg38"]["right"]
        results = [(f"{left}{repeat_seq}{right}", ["X"])]

        mock_apply.return_value = (results, {"mutated": "units"})

        _res, _mut_res, mut_units, mut_pos = apply_mutation_pipeline(
            args=args,
            config=minimal_config,
            results=results,
            mutation_name="dupC",
            dual_mutation_mode=False,
            mutation_pair=None,
        )

        assert mut_pos is not None
        assert len(mut_pos) == 1
        assert mut_units is not None

    @patch("muc_one_up.cli.mutations.apply_mutations")
    def test_dual_mode_with_explicit_targets(self, mock_apply: Mock, minimal_config: dict):
        """Given explicit targets in dual mode, when applying, then creates mutated version."""
        args = Mock()
        args.mutation_targets = ["1,1"]

        left = minimal_config["constants"]["hg38"]["left"]
        repeat_seq = minimal_config["repeats"]["X"]
        right = minimal_config["constants"]["hg38"]["right"]
        results = [(f"{left}{repeat_seq}{right}", ["X"])]

        mock_apply.return_value = (results, {"mutated": "units"})

        res, mut_res, mut_units, mut_pos = apply_mutation_pipeline(
            args=args,
            config=minimal_config,
            results=results,
            mutation_name="dupC",
            dual_mutation_mode=True,
            mutation_pair=("normal", "dupC"),
        )

        assert res == results  # Normal results unchanged in dual mode
        assert mut_res is not None  # Mutated results created
        assert mut_units is not None
        assert mut_pos == [(1, 1)]

    @patch("muc_one_up.cli.mutations.apply_mutations")
    def test_dual_mode_with_random_target(self, mock_apply: Mock, minimal_config: dict):
        """Given no explicit targets in dual mode, when applying, then finds random target."""
        args = Mock()
        args.mutation_targets = None

        left = minimal_config["constants"]["hg38"]["left"]
        repeat_seq = minimal_config["repeats"]["X"]
        right = minimal_config["constants"]["hg38"]["right"]
        results = [(f"{left}{repeat_seq}{right}", ["X"])]

        mock_apply.return_value = (results, {"mutated": "units"})

        res, mut_res, _mut_units, mut_pos = apply_mutation_pipeline(
            args=args,
            config=minimal_config,
            results=results,
            mutation_name="dupC",
            dual_mutation_mode=True,
            mutation_pair=("normal", "dupC"),
        )

        assert res == results
        assert mut_res is not None
        assert mut_pos is not None
        assert len(mut_pos) == 1

    def test_single_mode_mutation_error_propagates(self, minimal_config: dict):
        """Given invalid mutation in single mode, when applying, then raises MutationError."""
        args = Mock()
        args.mutation_targets = ["1,1"]

        left = minimal_config["constants"]["hg38"]["left"]
        repeat_seq = minimal_config["repeats"]["X"]
        right = minimal_config["constants"]["hg38"]["right"]
        results = [(f"{left}{repeat_seq}{right}", ["X"])]

        # Try to apply non-existent mutation
        with pytest.raises(MutationError, match="Mutation application failed"):
            apply_mutation_pipeline(
                args=args,
                config=minimal_config,
                results=results,
                mutation_name="nonExistentMutation",
                dual_mutation_mode=False,
                mutation_pair=None,
            )

    def test_dual_mode_mutation_error_propagates(self, minimal_config: dict):
        """Given invalid mutation in dual mode, when applying, then raises MutationError."""
        args = Mock()
        args.mutation_targets = ["1,1"]

        left = minimal_config["constants"]["hg38"]["left"]
        repeat_seq = minimal_config["repeats"]["X"]
        right = minimal_config["constants"]["hg38"]["right"]
        results = [(f"{left}{repeat_seq}{right}", ["X"])]

        with pytest.raises(MutationError, match="Dual mutation application failed"):
            apply_mutation_pipeline(
                args=args,
                config=minimal_config,
                results=results,
                mutation_name="dupC",
                dual_mutation_mode=True,
                mutation_pair=("normal", "nonExistentMutation"),
            )

    @patch("muc_one_up.cli.mutations.apply_mutations")
    def test_preserves_chain_in_dual_mode(self, mock_apply: Mock, minimal_config: dict):
        """Given dual mode, when applying mutation, then preserves original chains."""
        args = Mock()
        args.mutation_targets = ["1,1"]

        left = minimal_config["constants"]["hg38"]["left"]
        repeat_seq = minimal_config["repeats"]["X"]
        right = minimal_config["constants"]["hg38"]["right"]
        original_chain = ["X", "B"]
        results = [(f"{left}{repeat_seq}{right}", original_chain)]

        # Mock returns modified chain for mutated results
        mutated_chain = ["Xm", "B"]
        mock_apply.return_value = (
            [(f"{left}{repeat_seq}{right}", mutated_chain)],
            {"mutated": "units"},
        )

        res, mut_res, _mut_units, _mut_pos = apply_mutation_pipeline(
            args=args,
            config=minimal_config,
            results=results,
            mutation_name="dupC",
            dual_mutation_mode=True,
            mutation_pair=("normal", "dupC"),
        )

        # Original chain should be unchanged
        assert res[0][1] == original_chain
        # Mutated chain should be different
        assert mut_res[0][1] != original_chain

    @patch("muc_one_up.cli.mutations.apply_mutations")
    def test_multiple_haplotypes_with_random_target(self, mock_apply: Mock, minimal_config: dict):
        """Given multiple haplotypes, when finding random target, then selects valid repeat."""
        args = Mock()
        args.mutation_targets = None

        left = minimal_config["constants"]["hg38"]["left"]
        repeat_x = minimal_config["repeats"]["X"]
        repeat_b = minimal_config["repeats"]["B"]
        right = minimal_config["constants"]["hg38"]["right"]
        results = [
            (f"{left}{repeat_x}{right}", ["X"]),
            (f"{left}{repeat_b}{right}", ["B"]),
        ]

        mock_apply.return_value = (results, {"mutated": "units"})

        _res, _mut_res, _mut_units, mut_pos = apply_mutation_pipeline(
            args=args,
            config=minimal_config,
            results=results,
            mutation_name="dupC",
            dual_mutation_mode=False,
            mutation_pair=None,
        )

        assert mut_pos is not None
        assert len(mut_pos) == 1
        # Should target one of the valid repeats
        hap_idx, rep_idx = mut_pos[0]
        assert 1 <= hap_idx <= 2
        assert rep_idx == 1
