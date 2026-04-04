"""Tests for centralized sequence assembly."""

from __future__ import annotations

import pytest

from muc_one_up.assembly import assemble_sequence
from muc_one_up.type_defs import RepeatUnit


class TestAssembleSequence:
    """Tests for assemble_sequence()."""

    @pytest.fixture()
    def config(self):
        return {
            "repeats": {
                "1": "GCCACCACCTCCAACTCCT",
                "2": "GCCACCACC",
                "X": "GCCACCACCTCCAACTCCTGCCACC",
                "9": "TCTAGGACCTAGCTCCT",
            },
            "constants": {
                "hg38": {
                    "left": "ATGGCCCC",
                    "right": "CAATGGTG",
                }
            },
            "reference_assembly": "hg38",
        }

    def test_simple_chain(self, config):
        chain = [RepeatUnit("1"), RepeatUnit("2"), RepeatUnit("9")]
        repeats = config["repeats"]
        constants = config["constants"]["hg38"]
        seq = assemble_sequence(chain, repeats, constants)
        left = constants["left"]
        right = constants["right"]
        r1 = repeats["1"]
        r2 = repeats["2"]
        r9 = repeats["9"]
        assert seq == left + r1 + r2 + r9 + right

    def test_mutation_marker_stripped(self, config):
        """Chains with mutated RepeatUnit should look up the base symbol."""
        chain = [RepeatUnit("1"), RepeatUnit("X", True), RepeatUnit("9")]
        repeats = config["repeats"]
        constants = config["constants"]["hg38"]
        seq = assemble_sequence(chain, repeats, constants)
        left = constants["left"]
        right = constants["right"]
        r1 = repeats["1"]
        r_x = repeats["X"]
        r9 = repeats["9"]
        assert seq == left + r1 + r_x + r9 + right

    def test_empty_chain(self, config):
        """Empty chain produces only the left constant (no terminal '9' means no right)."""
        repeats = config["repeats"]
        constants = config["constants"]["hg38"]
        seq = assemble_sequence([], repeats, constants)
        left = constants["left"]
        assert seq == left

    def test_chain_not_ending_with_9_omits_right_constant(self, config):
        """If chain doesn't end with 9, right constant is omitted."""
        chain = [RepeatUnit("1"), RepeatUnit("2")]
        repeats = config["repeats"]
        constants = config["constants"]["hg38"]
        seq = assemble_sequence(chain, repeats, constants)
        left = constants["left"]
        r1 = repeats["1"]
        r2 = repeats["2"]
        assert seq == left + r1 + r2

    def test_unknown_symbol_raises(self, config):
        repeats = config["repeats"]
        constants = config["constants"]["hg38"]
        with pytest.raises(KeyError, match="UNKNOWN"):
            assemble_sequence([RepeatUnit("UNKNOWN")], repeats, constants)

    def test_hg19_assembly(self):
        repeats = {"1": "ACGT"}
        constants = {"left": "LEFT", "right": "RIGHT"}
        # Chain doesn't end with "9", so no right constant
        seq = assemble_sequence([RepeatUnit("1")], repeats, constants)
        assert seq == "LEFT" + "ACGT"

    def test_single_repeat_ending_with_9(self, config):
        chain = [RepeatUnit("9")]
        repeats = config["repeats"]
        constants = config["constants"]["hg38"]
        seq = assemble_sequence(chain, repeats, constants)
        left = constants["left"]
        right = constants["right"]
        r9 = repeats["9"]
        assert seq == left + r9 + right


def test_simulate_delegates_to_centralized_assembly():
    """simulate.assemble_haplotype_from_chain should delegate to assemble_sequence."""
    from unittest.mock import patch

    from muc_one_up.simulate import assemble_haplotype_from_chain

    config = {
        "repeats": {"1": "ACGT", "9": "TGCA"},
        "constants": {"hg38": {"left": "L", "right": "R"}},
        "reference_assembly": "hg38",
    }
    with patch("muc_one_up.simulate.assemble_sequence", return_value="MOCKED") as mock:
        result = assemble_haplotype_from_chain(["1", "9"], config)
    mock.assert_called_once_with(
        [RepeatUnit("1"), RepeatUnit("9")],
        config["repeats"],
        config["constants"]["hg38"],
    )
    assert result == "MOCKED"


def test_mutate_delegates_to_centralized_assembly():
    """mutate.rebuild_haplotype_sequence should delegate to assemble_sequence."""
    from unittest.mock import patch

    from muc_one_up.mutate import rebuild_haplotype_sequence

    config = {
        "repeats": {"1": "ACGT", "9": "TGCA"},
        "constants": {"hg38": {"left": "L", "right": "R"}},
        "reference_assembly": "hg38",
    }
    with patch("muc_one_up.mutate.assemble_sequence", return_value="MOCKED") as mock:
        result = rebuild_haplotype_sequence(["1", "9"], config)
    mock.assert_called_once_with(
        [RepeatUnit("1"), RepeatUnit("9")],
        config["repeats"],
        config["constants"]["hg38"],
    )
    assert result == "MOCKED"
