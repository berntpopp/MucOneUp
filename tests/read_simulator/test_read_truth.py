"""Tests for pbsim3 MAF parsing, read relabelling and the read truth manifest."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from muc_one_up.exceptions import ReadSimulationError
from muc_one_up.read_simulator.molecules import HpEdit, Molecule
from muc_one_up.read_simulator.read_truth import (
    TRUTH_COLUMNS,
    molecule_index,
    parse_maf_read_templates,
    relabel_reads,
    template_id,
)

MAF = """a
s m0000001 0 1500 + 1500 ACGT
s ont_1  0 1490 + 1490 ACGT

a
s m0000002 0 1500 + 1500 ACGT
s ont_2  0 1488 + 1488 ACGT

"""


def _write_gz(path: Path, text: str) -> Path:
    with gzip.open(path, "wt") as handle:
        handle.write(text)
    return path


def _molecules() -> dict[int, Molecule]:
    return {
        1: Molecule(1, 1, "full", "+", "ACGT", 0, 4),
        2: Molecule(2, 2, "smear", "-", "AC", 0, 4, (HpEdit(10, "C", 7, 8),), "deletion:1-3"),
    }


def test_template_id_roundtrip() -> None:
    assert template_id(42) == "m0000042"
    assert molecule_index("m0000042") == 42


def test_parse_maf_pairs_reads_with_templates(tmp_path: Path) -> None:
    maf = _write_gz(tmp_path / "x.maf.gz", MAF)
    assert parse_maf_read_templates(maf) == {"ont_1": "m0000001", "ont_2": "m0000002"}


def test_parse_maf_rejects_unpaired_block(tmp_path: Path) -> None:
    maf = _write_gz(tmp_path / "x.maf.gz", "a\ns m0000001 0 4 + 4 ACGT\n\n")
    with pytest.raises(ReadSimulationError, match="MAF"):
        parse_maf_read_templates(maf)


def test_relabel_reads_writes_unique_names_and_truth(tmp_path: Path) -> None:
    fq = _write_gz(tmp_path / "in.fq.gz", "@ont_2\nAC\n+\nII\n@ont_1\nACGT\n+\nIIII\n")
    mapping = {"ont_1": "m0000001", "ont_2": "m0000002"}
    out_fq, truth = tmp_path / "out.fastq", tmp_path / "truth.tsv.gz"
    n = relabel_reads([fq], mapping, _molecules(), "sample", out_fq, truth)
    assert n == 2
    names = [line.strip() for line in out_fq.read_text().splitlines()[::4]]
    assert names == ["@sample_h2_m0000002", "@sample_h1_m0000001"]
    with gzip.open(truth, "rt") as handle:
        rows = [line.rstrip("\n").split("\t") for line in handle]
    assert rows[0] == list(TRUTH_COLUMNS)
    row = dict(zip(TRUTH_COLUMNS, rows[1], strict=True))
    assert row["read_id"] == "sample_h2_m0000002" and row["kind"] == "smear"
    assert row["strand"] == "-" and row["n_hp_edits"] == "1" and row["detail"] == "deletion:1-3"
    assert row["hp_edits"] == "10:C:7>8"


def test_relabel_resolves_ccs_names_via_zmw(tmp_path: Path) -> None:
    fq = tmp_path / "hifi.fastq"
    fq.write_text("@hifi/1/ccs\nACGT\n+\nIIII\n")
    mapping = {"hifi/1/0": "m0000001", "hifi/1/1": "m0000001"}
    n = relabel_reads([fq], mapping, _molecules(), "s", tmp_path / "o.fastq", tmp_path / "t.tsv.gz")
    assert n == 1


def test_relabel_fails_loudly_on_unknown_read(tmp_path: Path) -> None:
    fq = tmp_path / "x.fastq"
    fq.write_text("@stranger\nA\n+\nI\n")
    with pytest.raises(ReadSimulationError, match="stranger"):
        relabel_reads([fq], {}, _molecules(), "s", tmp_path / "o.fastq", tmp_path / "t.tsv.gz")
