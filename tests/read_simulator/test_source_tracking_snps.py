"""Regression tests for SNP info handed to ReadSourceTracker (#111)."""

from __future__ import annotations

from types import SimpleNamespace
from unittest.mock import patch

import pytest

from muc_one_up.read_simulator.source_tracking import ReadSourceTracker

SNP_A = {"position": 10, "ref_base": "A", "alt_base": "G"}
SNP_B = {"position": 20, "ref_base": "C", "alt_base": "T"}
RESULTS = [SimpleNamespace(chain=[]), SimpleNamespace(chain=[])]
CONFIG = {"constants": {"hg38": {"left": "ACGT"}}, "repeats": {}}


@pytest.mark.parametrize(
    "applied,expected",
    [
        ({0: [SNP_A], 1: [SNP_B]}, {0: [SNP_A], 1: [SNP_B]}),  # producer returns a dict
        ([[SNP_A], [SNP_B]], {0: [SNP_A], 1: [SNP_B]}),  # list form stays supported
        ({0: [], 1: [SNP_B]}, {1: [SNP_B]}),
        (None, None),  # no SNPs: constructor receives None (existing contract)
    ],
)
def test_snp_info_keyed_by_haplotype(applied, expected) -> None:
    with patch.object(ReadSourceTracker, "__init__", return_value=None) as init:
        ReadSourceTracker.from_simulation_results(RESULTS, CONFIG, applied_snp_info=applied)
    assert init.call_args.kwargs["snp_info"] == expected
