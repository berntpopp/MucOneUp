import pytest

from muc_one_up.mutate import apply_mutations, validate_allowed_repeats


@pytest.fixture
def mutation_config():
    """A minimal config that includes a mutation definition."""
    return {
        "repeats": {"X": "XXXXX", "C": "CCCCC"},  # 5-base sequence
        "constants": {"hg38": {"left": "TTTT", "right": "GGGG"}},  # Use nested format
        "probabilities": {},
        "length_model": {},
        "mutations": {
            "testMut": {
                "allowed_repeats": ["X"],
                "changes": [{"type": "replace", "start": 2, "end": 2, "sequence": "Z"}],
            },
            "strictMut": {
                "allowed_repeats": ["X"],
                "strict_mode": True,
                "changes": [{"type": "replace", "start": 2, "end": 2, "sequence": "Z"}],
            },
            "invalidAllowedMut": {
                "allowed_repeats": ["Y"],  # Y is not a valid repeat
                "changes": [{"type": "replace", "start": 2, "end": 2, "sequence": "Z"}],
            },
        },
    }


def test_apply_mutations_forced_change(mutation_config):
    """
    Test that if the current symbol is not allowed then it is forcibly changed
    and the mutation applied. In this case, the chain has "C" but only "X" is allowed.
    After replacement, the mutated unit should be "XZXXX" and the chain marker becomes "Xm".
    """
    # Original haplotype: left constant "TTTT", repeat "CCCCC", right constant not appended because repeat != "9".
    results = [
        ("TTTTCCCCC", ["C"]),
    ]
    updated, mutated_units = apply_mutations(
        config=mutation_config,
        results=results,
        mutation_name="testMut",
        targets=[(1, 1)],
    )
    # Verify updated results
    assert len(updated) == 1
    new_seq, new_chain = updated[0]
    # Forced change: "C" should be replaced with "X" and then mutated marker appended ("Xm")
    assert new_chain[0] == "Xm"
    # Verify that the mutated unit "XZXXX" appears in the new sequence (i.e. left constant + mutated unit)
    assert "XZXXX" in new_seq
    # Ensure we no longer see the unmodified "XXXXX"
    assert "XXXXX" not in new_seq

    # Verify that the mutated_units dict has a record for haplotype 1, repeat 1.
    assert 1 in mutated_units
    found = False
    for rep_idx, unit_seq in mutated_units[1]:
        if rep_idx == 1:
            assert unit_seq == "XZXXX"
            found = True
    assert found


def test_apply_mutations_replace_ok(mutation_config):
    """Test a normal replacement mutation in an allowed repeat (X)."""
    results = [
        ("TTTTXXXXXGGGG", ["X"]),  # chain has allowed symbol "X"
    ]
    updated, mutated_units = apply_mutations(
        config=mutation_config,
        results=results,
        mutation_name="testMut",
        targets=[(1, 1)],
    )
    new_seq, new_chain = updated[0]
    # After mutation, chain symbol becomes "Xm"
    assert new_chain[0] == "Xm"
    # With the replacement at base position 2, "XXXXX" should become "XZXXX"
    assert "XZXXX" in new_seq
    # The right constant ("GGGG") should remain unchanged.
    assert new_seq.endswith("GGGG")
    # Verify mutated_units contains the expected mutated repeat.
    assert 1 in mutated_units
    found = False
    for rep_idx, unit_seq in mutated_units[1]:
        if rep_idx == 1:
            assert unit_seq == "XZXXX"
            found = True
    assert found


def test_apply_mutations_out_of_range(mutation_config):
    """
    If the change parameters are out of bounds for the repeat, a ValueError should be raised.
    """
    # Set change to be out-of-bounds (for a 5-base repeat)
    mutation_config["mutations"]["testMut"]["changes"] = [
        {"type": "replace", "start": 100, "end": 101, "sequence": "ABC"}
    ]
    results = [
        ("TTTTXXXXXGGGG", ["X"]),
    ]
    with pytest.raises(ValueError) as exc:
        apply_mutations(
            config=mutation_config,
            results=results,
            mutation_name="testMut",
            targets=[(1, 1)],
        )
    assert "out of bounds" in str(exc.value)


def test_validate_allowed_repeats():
    """
    Test that validate_allowed_repeats correctly identifies invalid repeat symbols.
    """
    config = {"repeats": {"X": "XXXXX", "C": "CCCCC"}}

    # Valid repeats should not raise an error
    mutation_def = {"allowed_repeats": ["X", "C"]}
    valid_repeats = validate_allowed_repeats(mutation_def, config)
    assert valid_repeats == {"X", "C"}

    # Invalid repeats should raise an error
    mutation_def = {"allowed_repeats": ["X", "Y", "Z"]}
    with pytest.raises(ValueError) as exc:
        validate_allowed_repeats(mutation_def, config)
    assert "Invalid repeats in allowed_repeats" in str(exc.value)
    assert "Y, Z" in str(exc.value)
    assert "Valid repeats are: C, X" in str(exc.value)


def test_strict_mode_enforcement(mutation_config):
    """
    Test that strict mode raises an error when encountering a disallowed repeat.
    """
    results = [
        (
            "TTTTCCCCC",
            ["C"],
        ),  # Chain has symbol "C" but strict mutation only allows "X"
    ]

    # In strict mode, an error should be raised
    with pytest.raises(ValueError) as exc:
        apply_mutations(
            config=mutation_config,
            results=results,
            mutation_name="strictMut",
            targets=[(1, 1)],
        )
    assert "Cannot apply mutation 'strictMut'" in str(exc.value)
    assert "Repeat symbol 'C' is not in allowed_repeats" in str(exc.value)

    # With non-strict mutation, it should force a change (tested in test_apply_mutations_forced_change)


def test_invalid_allowed_repeats_in_mutations(mutation_config):
    """
    Test that using a mutation with invalid allowed_repeats raises an error.
    """
    results = [
        ("TTTTCCCCC", ["C"]),
    ]

    with pytest.raises(ValueError) as exc:
        apply_mutations(
            config=mutation_config,
            results=results,
            mutation_name="invalidAllowedMut",
            targets=[(1, 1)],
        )
    assert "Invalid repeats in allowed_repeats" in str(exc.value)
    assert "Y" in str(exc.value)
    assert "Valid repeats are: C, X" in str(exc.value)
