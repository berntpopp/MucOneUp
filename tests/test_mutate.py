import pytest
from muc_one_up.mutate import apply_mutations

@pytest.fixture
def mutation_config():
    """A minimal config that includes a mutation definition."""
    return {
        "repeats": {
            "X": "XXXXX",   # 5-base sequence
            "C": "CCCCC"
        },
        "constants": {
            "left": "TTTT",
            "right": "GGGG"
        },
        "probabilities": {},
        "length_model": {},
        "mutations": {
            "testMut": {
                "allowed_repeats": ["X"],
                "changes": [
                    {
                        "type": "replace",
                        "start": 2,
                        "end": 2,
                        "sequence": "Z"
                    }
                ]
            }
        }
    }

def test_apply_mutations_forced_change(mutation_config):
    """
    If current_symbol is not in allowed_repeats, we forcibly change it to a random allowed symbol.
    Then we also apply the mutation. 'X' is 5 bases => 'XXXXX'.
    The replacement at start=2 => 'XZXXX'. Then the final 'm' => 'Xm'.
    """
    # "testMut" only allows "X", but our chain has "C". We'll see if it forces "X".
    results = [
        ("TTTTCCCCC", ["C"]),  # single haplotype with chain=[C]
    ]
    updated = apply_mutations(
        config=mutation_config,
        results=results,
        mutation_name="testMut",
        targets=[(1, 1)]
    )
    assert len(updated) == 1
    new_seq, new_chain = updated[0]
    # The chain symbol should now be "Xm" (forced from 'C' -> 'X' -> plus 'm')
    assert new_chain[0] == "Xm"

    # Because the code does not append the right constant unless final symbol == '9', 
    # new_seq will be:
    #  left const (TTTT) + mutated 'X' => 'XZXXX' + no right const
    # => 'TTTTXZXXX'
    #
    # Check that 'XZXXX' substring is present:
    assert "XZXXX" in new_seq
    # Ensure we do NOT see "XXXXX" unmutated
    assert "XXXXX" not in new_seq

def test_apply_mutations_replace_ok(mutation_config):
    """Test a normal replace in an allowed repeat (X)."""
    results = [
        ("TTTTXXXXXGGGG", ["X"]),  # chain has X => no forced change
    ]
    updated = apply_mutations(
        config=mutation_config,
        results=results,
        mutation_name="testMut",
        targets=[(1, 1)]
    )
    new_seq, new_chain = updated[0]
    # The chain symbol becomes "Xm" after mutation
    assert new_chain[0] == "Xm"
    # The mutation replaces base 2 with "Z", so 'XXXXX' => 'XZXXX'
    assert "XZXXX" in new_seq
    # Right const is "GGGG" so final seq is TTTT + XZXXX + GGGG = TTTTXZXXXGGGG
    assert new_seq.endswith("GGGG")

def test_apply_mutations_out_of_range(mutation_config):
    """If start/end are out of bounds for the repeat, we should get a ValueError."""
    mutation_config["mutations"]["testMut"]["changes"] = [
        {
            "type": "replace",
            "start": 100,  # out of range for 5-base repeat
            "end": 101,
            "sequence": "ABC"
        }
    ]
    results = [
        ("TTTTXXXXXGGGG", ["X"]),
    ]
    with pytest.raises(ValueError) as exc:
        apply_mutations(
            config=mutation_config,
            results=results,
            mutation_name="testMut",
            targets=[(1, 1)]
        )
    assert "out of bounds" in str(exc.value)
