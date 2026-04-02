import pytest

from muc_one_up.mutate import apply_mutations, validate_allowed_repeats
from muc_one_up.type_defs import HaplotypeResult, MutationTarget, RepeatUnit


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
    results = [
        HaplotypeResult(sequence="TTTTCCCCC", chain=[RepeatUnit("C")]),
    ]
    updated, mutated_units = apply_mutations(
        config=mutation_config,
        results=results,
        mutation_name="testMut",
        targets=[MutationTarget(1, 1)],
    )
    # Verify updated results
    assert len(updated) == 1
    hr = updated[0]
    # Forced change: "C" should be replaced with "X" and then mutated
    assert hr.chain[0].symbol == "X"
    assert hr.chain[0].mutated is True
    # Verify that the mutated unit "XZXXX" appears in the new sequence
    assert "XZXXX" in hr.sequence
    # Ensure we no longer see the unmodified "XXXXX"
    assert "XXXXX" not in hr.sequence

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
        HaplotypeResult(
            sequence="TTTTXXXXXGGGG", chain=[RepeatUnit("X")]
        ),  # chain has allowed symbol "X"
    ]
    updated, mutated_units = apply_mutations(
        config=mutation_config,
        results=results,
        mutation_name="testMut",
        targets=[MutationTarget(1, 1)],
    )
    hr = updated[0]
    # After mutation, chain symbol is "X" with mutated=True
    assert hr.chain[0].symbol == "X"
    assert hr.chain[0].mutated is True
    # With the replacement at base position 2, "XXXXX" should become "XZXXX"
    assert "XZXXX" in hr.sequence
    # The right constant ("GGGG") should remain unchanged.
    assert hr.sequence.endswith("GGGG")
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
        HaplotypeResult(sequence="TTTTXXXXXGGGG", chain=[RepeatUnit("X")]),
    ]
    with pytest.raises(ValueError) as exc:
        apply_mutations(
            config=mutation_config,
            results=results,
            mutation_name="testMut",
            targets=[MutationTarget(1, 1)],
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
        HaplotypeResult(
            sequence="TTTTCCCCC",
            chain=[RepeatUnit("C")],
        ),  # Chain has symbol "C" but strict mutation only allows "X"
    ]

    # In strict mode, an error should be raised
    with pytest.raises(ValueError) as exc:
        apply_mutations(
            config=mutation_config,
            results=results,
            mutation_name="strictMut",
            targets=[MutationTarget(1, 1)],
        )
    assert "Cannot apply mutation 'strictMut'" in str(exc.value)
    assert "Repeat symbol 'C' is not in allowed_repeats" in str(exc.value)

    # With non-strict mutation, it should force a change (tested in test_apply_mutations_forced_change)


def test_invalid_allowed_repeats_in_mutations(mutation_config):
    """
    Test that using a mutation with invalid allowed_repeats raises an error.
    """
    results = [
        HaplotypeResult(sequence="TTTTCCCCC", chain=[RepeatUnit("C")]),
    ]

    with pytest.raises(ValueError) as exc:
        apply_mutations(
            config=mutation_config,
            results=results,
            mutation_name="invalidAllowedMut",
            targets=[MutationTarget(1, 1)],
        )
    assert "Invalid repeats in allowed_repeats" in str(exc.value)
    assert "Y" in str(exc.value)
    assert "Valid repeats are: C, X" in str(exc.value)


@pytest.mark.unit
class TestMutateErrorConditions:
    """Comprehensive tests for error conditions in mutate module."""

    def test_missing_mutations_section_in_config(self):
        """Test error when config lacks mutations section."""
        config = {
            "repeats": {"X": "XXX"},
            "constants": {"hg38": {"left": "TT", "right": "GG"}},
        }
        results = [HaplotypeResult(sequence="TTXXXGG", chain=[RepeatUnit("X")])]

        with pytest.raises(ValueError, match="No 'mutations' section in config"):
            apply_mutations(config, results, "someMut", [MutationTarget(1, 1)])

    def test_mutation_name_not_found(self, mutation_config):
        """Test error when mutation name doesn't exist."""
        results = [HaplotypeResult(sequence="TTTTXXXXXGGGG", chain=[RepeatUnit("X")])]

        with pytest.raises(ValueError, match="Mutation 'nonExistent' not found"):
            apply_mutations(mutation_config, results, "nonExistent", [MutationTarget(1, 1)])

    def test_haplotype_index_out_of_range(self, mutation_config):
        """Test error when haplotype index exceeds number of haplotypes."""
        results = [HaplotypeResult(sequence="TTTTXXXXXGGGG", chain=[RepeatUnit("X")])]

        with pytest.raises(ValueError, match="Haplotype index 5 out of range"):
            apply_mutations(mutation_config, results, "testMut", [MutationTarget(5, 1)])

    def test_repeat_index_out_of_range(self, mutation_config):
        """Test error when repeat index exceeds chain length."""
        results = [HaplotypeResult(sequence="TTTTXXXXXGGGG", chain=[RepeatUnit("X")])]

        with pytest.raises(ValueError, match="Repeat index 10 out of range"):
            apply_mutations(mutation_config, results, "testMut", [MutationTarget(1, 10)])

    def test_empty_allowed_repeats_raises_error(self, mutation_config):
        """Test error when mutation has no allowed_repeats and symbol doesn't match."""
        mutation_config["mutations"]["emptyAllowed"] = {
            "allowed_repeats": [],
            "changes": [{"type": "replace", "start": 1, "end": 1, "sequence": "Z"}],
        }
        results = [HaplotypeResult(sequence="TTTTCCCCC", chain=[RepeatUnit("C")])]

        with pytest.raises(ValueError, match="has no allowed_repeats, cannot fix symbol"):
            apply_mutations(mutation_config, results, "emptyAllowed", [MutationTarget(1, 1)])

    def test_malformed_change_missing_fields(self, mutation_config):
        """Test error when change definition is missing required fields."""
        mutation_config["mutations"]["badChange"] = {
            "allowed_repeats": ["X"],
            "changes": [{"type": "replace", "start": 1}],  # Missing 'end'
        }
        results = [HaplotypeResult(sequence="TTTTXXXXXGGGG", chain=[RepeatUnit("X")])]

        with pytest.raises(ValueError, match=r"Malformed change.*missing fields"):
            apply_mutations(mutation_config, results, "badChange", [MutationTarget(1, 1)])

    def test_insert_out_of_bounds(self, mutation_config):
        """Test error when insert position is out of bounds."""
        mutation_config["mutations"]["badInsert"] = {
            "allowed_repeats": ["X"],
            "changes": [{"type": "insert", "start": 100, "end": 100, "sequence": "ZZZ"}],
        }
        results = [HaplotypeResult(sequence="TTTTXXXXXGGGG", chain=[RepeatUnit("X")])]

        with pytest.raises(ValueError, match="Insert out of bounds"):
            apply_mutations(mutation_config, results, "badInsert", [MutationTarget(1, 1)])

    def test_delete_out_of_bounds(self, mutation_config):
        """Test error when delete range is out of bounds."""
        mutation_config["mutations"]["badDelete"] = {
            "allowed_repeats": ["X"],
            "changes": [{"type": "delete", "start": 1, "end": 100, "sequence": ""}],
        }
        results = [HaplotypeResult(sequence="TTTTXXXXXGGGG", chain=[RepeatUnit("X")])]

        with pytest.raises(ValueError, match="Delete out of bounds"):
            apply_mutations(mutation_config, results, "badDelete", [MutationTarget(1, 1)])

    def test_delete_insert_operation(self, mutation_config):
        """Test delete_insert mutation type."""
        mutation_config["mutations"]["delIns"] = {
            "allowed_repeats": ["X"],
            "changes": [
                {"type": "delete_insert", "start": 1, "end": 4, "sequence": "ZZZ"}
            ],  # Delete between positions 1 and 4, insert ZZZ
        }
        results = [HaplotypeResult(sequence="TTTTXXXXXGGGG", chain=[RepeatUnit("X")])]

        updated, _ = apply_mutations(mutation_config, results, "delIns", [MutationTarget(1, 1)])

        # Should have applied delete_insert successfully
        assert len(updated) == 1
        hr = updated[0]
        assert hr.chain[0].symbol == "X"
        assert hr.chain[0].mutated is True

    def test_delete_insert_out_of_bounds(self, mutation_config):
        """Test error when delete_insert range is invalid."""
        mutation_config["mutations"]["badDelIns"] = {
            "allowed_repeats": ["X"],
            "changes": [{"type": "delete_insert", "start": 1, "end": 100, "sequence": "Z"}],
        }
        results = [HaplotypeResult(sequence="TTTTXXXXXGGGG", chain=[RepeatUnit("X")])]

        with pytest.raises(ValueError, match="delete_insert out of bounds"):
            apply_mutations(mutation_config, results, "badDelIns", [MutationTarget(1, 1)])

    def test_unknown_mutation_type(self, mutation_config):
        """Test error for unknown mutation type."""
        mutation_config["mutations"]["unknownType"] = {
            "allowed_repeats": ["X"],
            "changes": [{"type": "unknown_op", "start": 1, "end": 1, "sequence": "Z"}],
        }
        results = [HaplotypeResult(sequence="TTTTXXXXXGGGG", chain=[RepeatUnit("X")])]

        with pytest.raises(ValueError, match="Unknown mutation type 'unknown_op'"):
            apply_mutations(mutation_config, results, "unknownType", [MutationTarget(1, 1)])

    def test_chain_not_ending_with_9(self, mutation_config):
        """Test sequence assembly when chain doesn't end with '9'."""
        from muc_one_up.mutate import rebuild_haplotype_sequence

        mutation_config["repeats"]["A"] = "AAAAA"
        chain = ["X", "A"]  # Doesn't end with "9"

        seq = rebuild_haplotype_sequence(chain, mutation_config)

        # Should not include right constant
        left = mutation_config["constants"]["hg38"]["left"]
        right = mutation_config["constants"]["hg38"]["right"]
        assert seq.startswith(left)
        assert not seq.endswith(right)  # Right constant not added

    def test_multi_repeat_offset_calculation(self, mutation_config):
        """Test mutation targeting later repeats in chain (offset calculation)."""
        mutation_config["repeats"]["A"] = "AAAAA"
        mutation_config["repeats"]["B"] = "BBBBB"
        mutation_config["mutations"]["testMut"]["allowed_repeats"] = ["X", "A", "B"]

        # Chain with multiple repeats
        results = [
            HaplotypeResult(
                sequence="TTTTXXXXXAAAAABBBBB",
                chain=[RepeatUnit("X"), RepeatUnit("A"), RepeatUnit("B")],
            )
        ]

        # Target the 3rd repeat (B)
        updated, _ = apply_mutations(mutation_config, results, "testMut", [MutationTarget(1, 3)])

        hr = updated[0]
        # Third repeat should be marked
        assert hr.chain[2].symbol == "B"
        assert hr.chain[2].mutated is True
