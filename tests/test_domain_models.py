"""Tests for typed domain models."""

from __future__ import annotations

from muc_one_up.type_defs import RepeatUnit


class TestRepeatUnit:
    """Tests for RepeatUnit dataclass."""

    def test_create_normal_repeat(self):
        ru = RepeatUnit(symbol="1", mutated=False)
        assert ru.symbol == "1"
        assert ru.mutated is False

    def test_create_mutated_repeat(self):
        ru = RepeatUnit(symbol="X", mutated=True)
        assert ru.symbol == "X"
        assert ru.mutated is True

    def test_str_normal(self):
        ru = RepeatUnit(symbol="1", mutated=False)
        assert str(ru) == "1"

    def test_str_mutated(self):
        ru = RepeatUnit(symbol="X", mutated=True)
        assert str(ru) == "Xm"

    def test_from_str_normal(self):
        ru = RepeatUnit.from_str("1")
        assert ru.symbol == "1"
        assert ru.mutated is False

    def test_from_str_mutated(self):
        ru = RepeatUnit.from_str("Xm")
        assert ru.symbol == "X"
        assert ru.mutated is True

    def test_from_str_multi_char_symbol(self):
        ru = RepeatUnit.from_str("6p")
        assert ru.symbol == "6p"
        assert ru.mutated is False

    def test_from_str_multi_char_mutated(self):
        ru = RepeatUnit.from_str("6pm")
        assert ru.symbol == "6p"
        assert ru.mutated is True

    def test_equality(self):
        assert RepeatUnit("1", False) == RepeatUnit("1", False)
        assert RepeatUnit("1", True) != RepeatUnit("1", False)

    def test_chain_to_str_list(self):
        chain = [RepeatUnit("1", False), RepeatUnit("X", True), RepeatUnit("9", False)]
        assert [str(ru) for ru in chain] == ["1", "Xm", "9"]

    def test_chain_from_str_list(self):
        raw = ["1", "Xm", "9"]
        chain = [RepeatUnit.from_str(s) for s in raw]
        assert chain[0] == RepeatUnit("1", False)
        assert chain[1] == RepeatUnit("X", True)
        assert chain[2] == RepeatUnit("9", False)
