import pytest

from muc_one_up.simulate import simulate_diploid, simulate_single_haplotype


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
