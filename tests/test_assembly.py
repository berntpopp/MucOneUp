"""Tests for centralized sequence assembly."""

from __future__ import annotations

import pytest

from muc_one_up.assembly import assemble_sequence


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
        chain = ["1", "2", "9"]
        seq = assemble_sequence(chain, config)
        left = config["constants"]["hg38"]["left"]
        right = config["constants"]["hg38"]["right"]
        r1 = config["repeats"]["1"]
        r2 = config["repeats"]["2"]
        r9 = config["repeats"]["9"]
        assert seq == left + r1 + r2 + r9 + right

    def test_mutation_marker_stripped(self, config):
        """Chains with 'm' suffix should look up the base symbol."""
        chain = ["1", "Xm", "9"]
        seq = assemble_sequence(chain, config)
        left = config["constants"]["hg38"]["left"]
        right = config["constants"]["hg38"]["right"]
        r1 = config["repeats"]["1"]
        r_x = config["repeats"]["X"]
        r9 = config["repeats"]["9"]
        assert seq == left + r1 + r_x + r9 + right

    def test_empty_chain(self, config):
        """Empty chain produces only the left constant (no terminal '9' means no right)."""
        seq = assemble_sequence([], config)
        left = config["constants"]["hg38"]["left"]
        assert seq == left

    def test_chain_not_ending_with_9_omits_right_constant(self, config):
        """If chain doesn't end with 9, right constant is omitted."""
        chain = ["1", "2"]
        seq = assemble_sequence(chain, config)
        left = config["constants"]["hg38"]["left"]
        r1 = config["repeats"]["1"]
        r2 = config["repeats"]["2"]
        assert seq == left + r1 + r2

    def test_unknown_symbol_raises(self, config):
        with pytest.raises(KeyError, match="UNKNOWN"):
            assemble_sequence(["UNKNOWN"], config)

    def test_hg19_assembly(self):
        config = {
            "repeats": {"1": "ACGT"},
            "constants": {
                "hg19": {"left": "LEFT", "right": "RIGHT"},
                "hg38": {"left": "OTHER", "right": "OTHER"},
            },
            "reference_assembly": "hg19",
        }
        # Chain doesn't end with "9", so no right constant
        seq = assemble_sequence(["1"], config)
        assert seq == "LEFT" + "ACGT"

    def test_single_repeat_ending_with_9(self, config):
        chain = ["9"]
        seq = assemble_sequence(chain, config)
        left = config["constants"]["hg38"]["left"]
        right = config["constants"]["hg38"]["right"]
        r9 = config["repeats"]["9"]
        assert seq == left + r9 + right


def test_simulate_uses_centralized_assembly():
    """simulate.py should delegate to assembly.assemble_sequence."""
    import inspect

    import muc_one_up.simulate as mod

    source = inspect.getsource(mod)
    assert (
        "from .assembly import assemble_sequence" in source
        or "from muc_one_up.assembly import assemble_sequence" in source
    ), "simulate.py should import assemble_sequence from assembly module"


def test_mutate_uses_centralized_assembly():
    """mutate.py should delegate to assembly.assemble_sequence."""
    import inspect

    import muc_one_up.mutate as mod

    source = inspect.getsource(mod)
    assert (
        "from .assembly import assemble_sequence" in source
        or "from muc_one_up.assembly import assemble_sequence" in source
    ), "mutate.py should import assemble_sequence from assembly module"
