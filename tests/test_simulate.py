"""Comprehensive tests for muc_one_up.simulate module.

Following modern testing best practices:
- Parametrized tests for multiple scenarios
- Descriptive test names
- Edge case and error condition testing
- Proper use of fixtures
- AAA pattern (Arrange-Act-Assert)
"""

import random

import pytest

from muc_one_up.simulate import (
    assemble_haplotype_from_chain,
    pick_next_symbol_no_end,
    simulate_diploid,
    simulate_from_chains,
    simulate_single_haplotype,
)


@pytest.fixture
def simple_config():
    """Returns a minimal config dictionary for testing."""
    return {
        "repeats": {
            "1": "AAA",
            "2": "CCC",
            "9": "GGG",  # final block usage
        },
        "constants": {"hg38": {"left": "TTTT", "right": "AAAA"}},  # Use nested format
        "probabilities": {"1": {"2": 1.0}, "2": {"9": 1.0}, "9": {"END": 1.0}},
        "length_model": {
            "distribution": "normal",
            "min_repeats": 2,
            "max_repeats": 5,
            "mean_repeats": 3,
        },
    }


def test_simulate_single_haplotype_min_len(simple_config):
    """If the target length is below min_length=10 in the function, it raises ValueError."""
    with pytest.raises(ValueError):
        simulate_single_haplotype(simple_config, target_length=5)  # below default min_length=10


def test_simulate_single_haplotype_override_min_len(simple_config):
    """
    If we override the function to allow min_length=2,
    we can test building a short chain without error.
    """

    def patched(cfg, target_length):
        return simulate_single_haplotype(cfg, target_length, min_length=2)

    seq, chain = patched(simple_config, target_length=3)
    # We have a small chain. Possibly 2 or 3 repeats, depending on logic.
    # Just confirm it produced something valid:
    assert seq.startswith("TTTT")  # left const
    assert "AAA" in seq or "CCC" in seq  # we used '1' -> '2'
    assert seq.endswith(
        "AAAA"
    )  # right const if final repeat was '9'? Maybe not if we never forced it.


def test_simulate_diploid_fixed_length(simple_config):
    """
    Here we override the min_length so that fixed_length=3 won't raise an error.
    This is just to ensure the code can produce 2 haplotypes of length=3
    in a synthetic scenario.
    """

    def patched_sim_single(cfg, target_length):
        return simulate_single_haplotype(cfg, target_length, min_length=2)

    # We'll monkeypatch simulate_single_haplotype inside simulate_diploid
    original_func = simulate_single_haplotype
    try:
        from muc_one_up import simulate

        simulate.simulate_single_haplotype = patched_sim_single

        results = simulate_diploid(
            config=simple_config, num_haplotypes=2, fixed_lengths=[3, 3], seed=42
        )
        assert len(results) == 2
        for seq, chain in results:
            # chain might have 2 or 3 repeats. Just confirm the sequence is non-empty
            assert len(chain) > 0
            assert seq.startswith("TTTT")
    finally:
        # Restore original function so other tests remain unaffected
        simulate.simulate_single_haplotype = original_func


# ====================================================================================
# Modern comprehensive tests following best practices
# ====================================================================================


@pytest.mark.unit
class TestPickNextSymbolNoEnd:
    """Comprehensive tests for pick_next_symbol_no_end function."""

    def test_returns_none_when_symbol_not_in_probabilities(self):
        """Given unknown symbol, when picking next, then returns None."""
        probabilities = {"1": {"2": 1.0}}
        result = pick_next_symbol_no_end(probabilities, "unknown")
        assert result is None

    def test_returns_none_when_all_options_are_forbidden(self):
        """Given only forbidden symbols, when picking, then returns None."""
        probabilities = {"test": {"6": 0.25, "6p": 0.25, "9": 0.25, "END": 0.25}}
        result = pick_next_symbol_no_end(probabilities, "test")
        assert result is None

    def test_picks_from_valid_symbols_excluding_forbidden(self, minimal_config: dict):
        """Given mix of valid and forbidden symbols, when picking, then selects only valid."""
        random.seed(42)
        probabilities = minimal_config["probabilities"]

        # Run multiple times to verify randomness
        results = [pick_next_symbol_no_end(probabilities, "1") for _ in range(10)]

        forbidden = {"6", "6p", "9", "END"}
        assert all(r not in forbidden for r in results if r is not None)


@pytest.mark.unit
class TestAssembleHaplotypeFromChain:
    """Comprehensive tests for assemble_haplotype_from_chain function."""

    def test_assembles_with_left_and_right_constants(self, minimal_config: dict):
        """Given valid chain, when assembling, then includes both constants."""
        chain = ["1", "2", "X"]
        result = assemble_haplotype_from_chain(chain, minimal_config)

        left = minimal_config["constants"]["hg38"]["left"]
        right = minimal_config["constants"]["hg38"]["right"]
        assert result.startswith(left)
        assert result.endswith(right)

    def test_handles_mutation_markers(self, minimal_config: dict):
        """Given chain with 'm' markers, when assembling, then strips marker for lookup."""
        chain = ["1", "2m", "Xm"]
        result = assemble_haplotype_from_chain(chain, minimal_config)

        # Should successfully use base symbols
        assert minimal_config["repeats"]["2"] in result
        assert minimal_config["repeats"]["X"] in result

    def test_raises_error_for_unknown_symbol(self, minimal_config: dict):
        """Given unknown symbol, when assembling, then raises ValueError."""
        chain = ["1", "UNKNOWN"]

        with pytest.raises(ValueError, match="Symbol 'UNKNOWN' not found"):
            assemble_haplotype_from_chain(chain, minimal_config)

    def test_uses_hg19_when_specified(self, minimal_config: dict):
        """Given hg19 assembly, when assembling, then uses hg19 constants."""
        minimal_config["constants"]["hg19"] = {
            "left": "GGGG" * 20,
            "right": "CCCC" * 20,
        }
        minimal_config["reference_assembly"] = "hg19"
        chain = ["1"]

        result = assemble_haplotype_from_chain(chain, minimal_config)

        assert result.startswith("GGGG" * 20)
        assert result.endswith("CCCC" * 20)

    def test_assembles_empty_chain(self, minimal_config: dict):
        """Given empty chain, when assembling, then returns only constants."""
        chain: list[str] = []
        result = assemble_haplotype_from_chain(chain, minimal_config)

        left = minimal_config["constants"]["hg38"]["left"]
        right = minimal_config["constants"]["hg38"]["right"]
        assert result == left + right


@pytest.mark.unit
class TestSimulateFromChains:
    """Comprehensive tests for simulate_from_chains function."""

    def test_creates_haplotypes_from_chains(self, minimal_config: dict):
        """Given predefined chains, when simulating, then creates haplotypes."""
        chains = [
            ["1", "2", "X", "B", "6", "7", "8", "9"],
            ["1", "2", "A", "B", "6p", "7", "8", "9"],
        ]

        results = simulate_from_chains(chains, minimal_config)

        assert len(results) == 2
        assert all(isinstance(r, tuple) and len(r) == 2 for r in results)

    def test_preserves_original_chains(self, minimal_config: dict):
        """Given chains, when simulating, then preserves chain structure."""
        chains = [["1", "2", "X"], ["1", "A", "B"]]

        results = simulate_from_chains(chains, minimal_config)

        _, chain1 = results[0]
        _, chain2 = results[1]
        assert chain1 == ["1", "2", "X"]
        assert chain2 == ["1", "A", "B"]

    def test_applies_mutation_to_specific_haplotype(self, minimal_config: dict):
        """Given mutation target, when simulating, then applies mutation correctly."""
        chains = [["1", "2", "X"], ["1", "A", "B"]]
        targets = [(1, 2)]  # Haplotype 1, position 2

        results = simulate_from_chains(chains, minimal_config, "dupC", targets)

        _, chain1 = results[0]
        _, chain2 = results[1]
        assert chain1[1].endswith("m")  # Position 2 (0-indexed: 1)
        assert not any(r.endswith("m") for r in chain2)

    def test_applies_multiple_mutations(self, minimal_config: dict):
        """Given multiple targets, when simulating, then applies all mutations."""
        chains = [["1", "2", "X"], ["1", "A", "B"]]
        targets = [(1, 2), (2, 3)]

        results = simulate_from_chains(chains, minimal_config, "dupC", targets)

        _, chain1 = results[0]
        _, chain2 = results[1]
        assert chain1[1].endswith("m")
        assert chain2[2].endswith("m")

    def test_skips_out_of_range_mutations(self, minimal_config: dict):
        """Given out-of-range target, when simulating, then skips mutation."""
        chains = [["1", "2"]]  # Only 2 repeats
        targets = [(1, 10)]  # Position 10 doesn't exist

        results = simulate_from_chains(chains, minimal_config, "dupC", targets)

        _, chain = results[0]
        assert not any(r.endswith("m") for r in chain)

    def test_does_not_duplicate_mutation_markers(self, minimal_config: dict):
        """Given already-marked repeat, when mutating, then doesn't add duplicate marker."""
        chains = [["1", "2m"]]  # Already mutated
        targets = [(1, 2)]  # Target same position

        results = simulate_from_chains(chains, minimal_config, "dupC", targets)

        _, chain = results[0]
        assert chain[1] == "2m"  # Not "2mm"


@pytest.mark.unit
class TestSimulateSingleHaplotypeAdvanced:
    """Advanced tests for simulate_single_haplotype function."""

    def test_produces_non_empty_haplotype(self, minimal_config: dict):
        """Given target length, when simulating, then produces valid haplotype."""
        random.seed(42)
        target = 20

        seq, chain = simulate_single_haplotype(minimal_config, target)

        assert len(chain) > 0
        assert len(seq) > 0
        # Verify has constants
        left = minimal_config["constants"]["hg38"]["left"]
        right = minimal_config["constants"]["hg38"]["right"]
        assert seq.startswith(left)
        assert seq.endswith(right)

    def test_uses_hg38_constants_by_default(self, minimal_config: dict):
        """Given no assembly specified, when simulating, then uses hg38."""
        random.seed(42)

        seq, _ = simulate_single_haplotype(minimal_config, 15)

        left = minimal_config["constants"]["hg38"]["left"]
        assert seq.startswith(left)

    def test_chain_always_starts_with_1(self, minimal_config: dict):
        """Given any simulation, when complete, then chain starts with '1'."""
        random.seed(42)

        _, chain = simulate_single_haplotype(minimal_config, 15)

        assert chain[0] == "1"


@pytest.mark.unit
class TestSimulateDiploidAdvanced:
    """Advanced tests for simulate_diploid function."""

    def test_creates_multiple_haplotypes(self, minimal_config: dict):
        """Given num_haplotypes, when simulating, then creates correct count."""
        random.seed(42)

        results = simulate_diploid(minimal_config, num_haplotypes=3)

        assert len(results) == 3
        assert all(isinstance(r, tuple) and len(r) == 2 for r in results)

    def test_produces_valid_sequences(self, minimal_config: dict):
        """Given simulation, when complete, then sequences are valid."""
        random.seed(42)

        results = simulate_diploid(minimal_config, num_haplotypes=2)

        left = minimal_config["constants"]["hg38"]["left"]
        right = minimal_config["constants"]["hg38"]["right"]

        for seq, chain in results:
            assert len(seq) > 0
            assert len(chain) > 0
            assert seq.startswith(left)
            assert seq.endswith(right)

    def test_seed_ensures_reproducibility(self, minimal_config: dict):
        """Given same seed, when simulating, then produces identical results."""
        r1 = simulate_diploid(minimal_config, num_haplotypes=2, seed=42)
        r2 = simulate_diploid(minimal_config, num_haplotypes=2, seed=42)

        assert r1[0][1] == r2[0][1]  # Same chains
        assert r1[1][1] == r2[1][1]


@pytest.mark.integration
class TestSimulationIntegration:
    """Integration tests for complete workflows."""

    def test_complete_diploid_workflow(self, minimal_config: dict):
        """Test end-to-end diploid simulation."""
        random.seed(42)

        results = simulate_diploid(minimal_config, num_haplotypes=2, seed=42)

        # Verify structure
        assert len(results) == 2
        seq1, chain1 = results[0]
        seq2, chain2 = results[1]

        # Verify basic properties
        assert len(chain1) > 0
        assert len(chain2) > 0
        assert chain1[0] == "1"
        assert chain2[0] == "1"

        # Verify sequences have correct structure
        left = minimal_config["constants"]["hg38"]["left"]
        right = minimal_config["constants"]["hg38"]["right"]
        assert seq1.startswith(left) and seq1.endswith(right)
        assert seq2.startswith(left) and seq2.endswith(right)

    def test_predefined_chains_with_mutations(self, minimal_config: dict):
        """Test predefined chains simulation with mutations."""
        chains = [
            ["1", "2", "X", "B", "6", "7", "8", "9"],
            ["1", "2", "A", "B", "6p", "7", "8", "9"],
        ]
        targets = [(1, 3), (2, 4)]

        results = simulate_from_chains(chains, minimal_config, "dupC", targets)

        _, chain1 = results[0]
        _, chain2 = results[1]

        # Verify mutations
        assert chain1[2] == "Xm"
        assert chain2[3] == "Bm"
