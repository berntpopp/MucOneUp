import pytest

from muc_one_up.mutate import (
    apply_changes_to_repeat,
    apply_mutations,
    validate_allowed_repeats,
)
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
    valid_symbols = {"X", "C"}

    # Valid repeats should not raise an error
    mutation_def = {"allowed_repeats": ["X", "C"]}
    valid_repeats = validate_allowed_repeats(mutation_def, valid_symbols)
    assert valid_repeats == {"X", "C"}

    # Invalid repeats should raise an error
    mutation_def = {"allowed_repeats": ["X", "Y", "Z"]}
    with pytest.raises(ValueError) as exc:
        validate_allowed_repeats(mutation_def, valid_symbols)
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


class TestMutationAliasing:
    """Prove that apply_mutations does not mutate caller-owned data."""

    def test_original_results_not_mutated(self, mutation_config):
        """Caller's results list and HaplotypeResult objects must be untouched."""
        original_chain = [RepeatUnit("X")]
        original_seq = "TTTTXXXXXGGGG"
        results = [HaplotypeResult(sequence=original_seq, chain=original_chain)]

        # Snapshot before mutation
        chain_before = [RepeatUnit(ru.symbol, ru.mutated) for ru in original_chain]
        seq_before = original_seq

        _updated, _ = apply_mutations(
            config=mutation_config,
            results=results,
            mutation_name="testMut",
            targets=[MutationTarget(1, 1)],
        )

        # Caller's objects must be completely unchanged
        assert results[0].sequence == seq_before
        assert len(results[0].chain) == len(chain_before)
        for orig, before in zip(results[0].chain, chain_before, strict=True):
            assert orig.symbol == before.symbol
            assert orig.mutated == before.mutated

    def test_chain_objects_not_shared(self, mutation_config):
        """Updated chain list must be a different object from the original."""
        original_chain = [RepeatUnit("X")]
        results = [HaplotypeResult(sequence="TTTTXXXXXGGGG", chain=original_chain)]

        updated, _ = apply_mutations(
            config=mutation_config,
            results=results,
            mutation_name="testMut",
            targets=[MutationTarget(1, 1)],
        )

        # The chain list in updated must be a different object
        assert updated[0].chain is not results[0].chain

    def test_forced_change_does_not_alias(self, mutation_config):
        """Even the forced-change path (non-strict, wrong symbol) must not alias."""
        original_chain = [RepeatUnit("C")]
        results = [HaplotypeResult(sequence="TTTTCCCCC", chain=original_chain)]

        _updated, _ = apply_mutations(
            config=mutation_config,
            results=results,
            mutation_name="testMut",
            targets=[MutationTarget(1, 1)],
        )

        # Original chain must still be "C", not mutated
        assert results[0].chain[0].symbol == "C"
        assert results[0].chain[0].mutated is False


def test_multiple_length_changing_mutations_order_independence():
    """Verify that multiple length-changing mutations on the same haplotype
    in both ascending and descending target order produce the exact expected
    DNA sequence without coordinate drift (#134).
    """
    config = {
        "repeats": {"1": "AAAA", "2": "CCCC", "3": "GGGG"},
        "constants": {"hg38": {"left": "LL", "right": "RR"}},
        "mutations": {
            "dup": {
                "allowed_repeats": ["1", "2", "3"],
                "strict_mode": True,
                "changes": [{"type": "insert", "start": 2, "end": 2, "sequence": "IN"}],
            },
            "del": {
                "allowed_repeats": ["1", "2", "3"],
                "strict_mode": True,
                "changes": [{"type": "delete", "start": 2, "end": 3, "sequence": ""}],
            },
        },
    }
    init_seq = "LLAAAACCCCGGGGRR"

    # Testing insertions on repeat 1 and repeat 3:
    # Repeat 1: "AAAA" with insert at pos 2 -> "A" + "IN" + "AAA" = "AINAAA"
    # Repeat 2: "CCCC"
    # Repeat 3: "GGGG" with insert at pos 2 -> "G" + "IN" + "GGG" = "GINGGG"
    expected_dup_seq = "LLAINAAACCCCGINGGGRR"

    # Ascending target order: repeat 1 then repeat 3
    results_asc = [
        HaplotypeResult(
            sequence=init_seq,
            chain=[RepeatUnit("1"), RepeatUnit("2"), RepeatUnit("3")],
        )
    ]
    updated_asc, _ = apply_mutations(
        config=config,
        results=results_asc,
        mutation_name="dup",
        targets=[MutationTarget(1, 1), MutationTarget(1, 3)],
    )
    assert updated_asc[0].sequence == expected_dup_seq

    # Descending target order: repeat 3 then repeat 1
    results_desc = [
        HaplotypeResult(
            sequence=init_seq,
            chain=[RepeatUnit("1"), RepeatUnit("2"), RepeatUnit("3")],
        )
    ]
    updated_desc, _ = apply_mutations(
        config=config,
        results=results_desc,
        mutation_name="dup",
        targets=[MutationTarget(1, 3), MutationTarget(1, 1)],
    )
    assert updated_desc[0].sequence == expected_dup_seq
    assert updated_asc[0].sequence == updated_desc[0].sequence

    # Testing mixed length changes: deletion on adjacent repeats in ascending vs descending
    # Repeat 1: "AAAA" -> "AA" (delete bases 2..3)
    # Repeat 2: "CCCC" -> "CC" (delete bases 2..3)
    expected_del_seq = "LLAACCGGGGRR"

    res_del_asc = [
        HaplotypeResult(
            sequence=init_seq,
            chain=[RepeatUnit("1"), RepeatUnit("2"), RepeatUnit("3")],
        )
    ]
    updated_del_asc, _ = apply_mutations(
        config=config,
        results=res_del_asc,
        mutation_name="del",
        targets=[MutationTarget(1, 1), MutationTarget(1, 2)],
    )
    assert updated_del_asc[0].sequence == expected_del_seq

    res_del_desc = [
        HaplotypeResult(
            sequence=init_seq,
            chain=[RepeatUnit("1"), RepeatUnit("2"), RepeatUnit("3")],
        )
    ]
    updated_del_desc, _ = apply_mutations(
        config=config,
        results=res_del_desc,
        mutation_name="del",
        targets=[MutationTarget(1, 2), MutationTarget(1, 1)],
    )
    assert updated_del_desc[0].sequence == expected_del_seq
    assert updated_del_asc[0].sequence == updated_del_desc[0].sequence


def test_non_strict_symbol_conversion_preserves_earlier_indel():
    """Verify that non-strict repeat symbol conversion at a second target
    preserves an indel mutation already applied at an earlier target on the same haplotype (#134).
    """
    config = {
        "repeats": {"X": "XXXX", "C": "CCCC"},
        "constants": {"hg38": {"left": "LL", "right": "RR"}},
        "mutations": {
            "dupNonStrict": {
                "allowed_repeats": ["X"],
                "strict_mode": False,
                "changes": [{"type": "insert", "start": 2, "end": 2, "sequence": "IN"}],
            },
        },
    }
    init_seq = "LLXXXXCCCCRR"
    results = [
        HaplotypeResult(
            sequence=init_seq,
            chain=[RepeatUnit("X"), RepeatUnit("C")],
        )
    ]

    # Target 1: repeat 1 has allowed symbol 'X' -> mutates to 'XINXXX' (len 6, +2 bases)
    # Target 2: repeat 2 has disallowed symbol 'C' -> auto-converts to 'X' ('XXXX'), then mutates to 'XINXXX'
    # Auto-conversion must NOT reassemble and wipe repeat 1's insertion!
    updated, _ = apply_mutations(
        config=config,
        results=results,
        mutation_name="dupNonStrict",
        targets=[MutationTarget(1, 1), MutationTarget(1, 2)],
    )
    assert updated[0].sequence == "LLXINXXXXINXXXRR"
    assert updated[0].chain[0].symbol == "X"
    assert updated[0].chain[0].mutated is True
    assert updated[0].chain[1].symbol == "X"
    assert updated[0].chain[1].mutated is True


def test_two_mutations_on_same_repeat_preserves_edits():
    """Verify that applying two mutations to the same (haplotype_index, repeat_index)
    preserves the first edit and applies the second on top of it (#134).
    """
    config = {
        "repeats": {"1": "ABCDE"},
        "constants": {"hg38": {"left": "LL", "right": "RR"}},
        "mutations": {
            "ins": {
                "allowed_repeats": ["1"],
                "strict_mode": True,
                "changes": [{"type": "insert", "start": 2, "end": 2, "sequence": "Z"}],
            },
            "rep": {
                "allowed_repeats": ["1"],
                "strict_mode": True,
                "changes": [{"type": "replace", "start": 1, "end": 1, "sequence": "W"}],
            },
        },
    }
    init_seq = "LLABCDERR"
    results = [HaplotypeResult(sequence=init_seq, chain=[RepeatUnit("1")])]

    # Target 1 inserts "Z" at 2: "ABCDE" -> "AZBCDE"
    # Target 2 inserts "Z" at 2 again on top of the modified repeat: "AZBCDE" -> "AZZBCDE"
    updated, mutated_units = apply_mutations(
        config=config,
        results=results,
        mutation_name="ins",
        targets=[MutationTarget(1, 1), MutationTarget(1, 1)],
    )
    assert updated[0].sequence == "LLAZZBCDERR"
    assert mutated_units[1] == [(1, "AZZBCDE")]

    # Sequential calls with different mutations on the same repeat:
    # First: replace at 1 with 'W' ("ABCDE" -> "WBCDE")
    # Second: replace at 5 with 'Z' ("WBCDE" -> "WBCDZ")
    res1, _ = apply_mutations(config, results, "rep", [MutationTarget(1, 1)])
    assert res1[0].sequence == "LLWBCDERR"

    config["mutations"]["rep5"] = {
        "allowed_repeats": ["1"],
        "strict_mode": True,
        "changes": [{"type": "replace", "start": 5, "end": 5, "sequence": "Z"}],
    }
    res2, _ = apply_mutations(config, res1, "rep5", [MutationTarget(1, 1)])
    assert res2[0].sequence == "LLWBCDZRR"


def test_apply_mutations_sequential_calls_tracks_repeat_lengths():
    """Calling apply_mutations twice sequentially where call 1 inserts into repeat 1
    and call 2 mutates repeat 2 correctly offsets repeat 2 bounds ([AAAA][CCCC] -> [AGAAA][TCCC]).
    """
    config = {
        "repeats": {"1": "AAAA", "2": "CCCC"},
        "constants": {"hg38": {"left": "LL", "right": "RR"}},
        "mutations": {
            "ins": {
                "allowed_repeats": ["1"],
                "strict_mode": True,
                "changes": [{"type": "insert", "start": 2, "end": 2, "sequence": "G"}],
            },
            "sub": {
                "allowed_repeats": ["2"],
                "strict_mode": True,
                "changes": [{"type": "replace", "start": 1, "end": 1, "sequence": "T"}],
            },
        },
    }
    init_seq = "LLAAAACCCCRR"
    results = [HaplotypeResult(sequence=init_seq, chain=[RepeatUnit("1"), RepeatUnit("2")])]

    # Call 1: insert into repeat 1: "AAAA" -> "AGAAA" (length grows from 4 to 5)
    res1, units1 = apply_mutations(config, results, "ins", [MutationTarget(1, 1)])
    assert res1[0].sequence == "LLAGAAACCCCRR"
    assert res1[0]._repeat_lengths == [5, 4]
    assert units1[1] == [(1, "AGAAA")]

    # Call 2: mutate repeat 2 on res1: "CCCC" -> "TCCC"
    res2, units2 = apply_mutations(config, res1, "sub", [MutationTarget(1, 2)])
    assert res2[0].sequence == "LLAGAAATCCCRR"
    assert res2[0]._repeat_lengths == [5, 4]
    assert units2[1] == [(2, "TCCC")]

    # Also test fallback delta inference when _repeat_lengths is stripped/None
    res1_stripped = [HaplotypeResult(sequence=res1[0].sequence, chain=res1[0].chain)]
    assert res1_stripped[0]._repeat_lengths is None
    res2_fallback, units2_fb = apply_mutations(config, res1_stripped, "sub", [MutationTarget(1, 2)])
    assert res2_fallback[0].sequence == "LLAGAAATCCCRR"
    assert res2_fallback[0]._repeat_lengths == [5, 4]
    assert units2_fb[1] == [(2, "TCCC")]


def test_apply_changes_to_repeat_multi_operation_updates_bounds():
    """Multi-operation mutation in apply_changes_to_repeat validates bounds
    against the updated length after prior insertions/deletions.
    """
    config = {
        "repeats": {"1": "AAAA"},
        "constants": {"hg38": {"left": "LL", "right": "RR"}},
    }
    chain = [RepeatUnit("1")]
    seq = "LLAAAARR"

    # Growth: insert 4 chars so repeat length grows from 4 to 8,
    # then modify at position 7 (which would fail if length was still 4).
    changes_grow = [
        {"type": "insert", "start": 2, "end": 2, "sequence": "XYZW"},
        {"type": "replace", "start": 7, "end": 8, "sequence": "TT"},
    ]
    new_seq, _, mutated = apply_changes_to_repeat(
        seq=seq,
        chain=chain,
        repeat_index=0,
        changes=changes_grow,
        config=config,
        mutation_name="grow_test",
        start=2,
        stop=6,
    )
    assert mutated == "AXYZWATT"
    assert new_seq == "LLAXYZWATTRR"

    # Shrink: delete 3 chars so repeat length shrinks from 4 to 1,
    # then an operation at position 3 should fail bounds validation against length 1.
    changes_shrink = [
        {"type": "delete", "start": 2, "end": 4, "sequence": ""},
        {"type": "replace", "start": 3, "end": 3, "sequence": "C"},
    ]
    with pytest.raises(ValueError, match=r"Replace out of bounds.*repeat length=1"):
        apply_changes_to_repeat(
            seq=seq,
            chain=chain,
            repeat_index=0,
            changes=changes_shrink,
            config=config,
            mutation_name="shrink_test",
            start=2,
            stop=6,
        )
