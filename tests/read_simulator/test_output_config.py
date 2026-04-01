"""Tests for OutputConfig dataclass."""

from __future__ import annotations

from pathlib import Path

import pytest

from muc_one_up.read_simulator.output_config import OutputConfig


class TestOutputConfigConstruction:
    """Test basic construction of OutputConfig."""

    def test_basic_construction(self) -> None:
        cfg = OutputConfig(out_dir=Path("/tmp/out"), out_base="sample")
        assert cfg.out_dir == Path("/tmp/out")
        assert cfg.out_base == "sample"

    def test_frozen(self) -> None:
        cfg = OutputConfig(out_dir=Path("/tmp/out"), out_base="sample")
        with pytest.raises(AttributeError):
            cfg.out_dir = Path("/other")  # type: ignore[misc]


class TestDerivePath:
    """Test derive_path and derive_path_str."""

    def test_derive_path_bam(self) -> None:
        cfg = OutputConfig(out_dir=Path("/data/output"), out_base="reads")
        assert cfg.derive_path(".bam") == Path("/data/output/reads.bam")

    def test_derive_path_fastq_gz(self) -> None:
        cfg = OutputConfig(out_dir=Path("/data/output"), out_base="reads")
        assert cfg.derive_path("_R1.fastq.gz") == Path("/data/output/reads_R1.fastq.gz")

    def test_derive_path_str(self) -> None:
        cfg = OutputConfig(out_dir=Path("/data/output"), out_base="reads")
        assert cfg.derive_path_str(".bam") == "/data/output/reads.bam"


class TestFromInputFasta:
    """Test from_input_fasta classmethod."""

    def test_no_overrides(self) -> None:
        cfg = OutputConfig.from_input_fasta(
            input_fa="/data/refs/sample.fasta",
            out_dir=None,
            out_base=None,
            suffix="_sim",
        )
        assert cfg.out_dir == Path("/data/refs")
        assert cfg.out_base == "sample_sim"

    def test_both_overrides(self) -> None:
        cfg = OutputConfig.from_input_fasta(
            input_fa="/data/refs/sample.fasta",
            out_dir="/custom/dir",
            out_base="custom_base",
            suffix="_sim",
        )
        assert cfg.out_dir == Path("/custom/dir")
        assert cfg.out_base == "custom_base"

    def test_only_out_dir(self) -> None:
        cfg = OutputConfig.from_input_fasta(
            input_fa="/data/refs/sample.fasta",
            out_dir="/custom/dir",
            out_base=None,
            suffix="_sim",
        )
        assert cfg.out_dir == Path("/custom/dir")
        assert cfg.out_base == "sample_sim"

    def test_only_out_base(self) -> None:
        cfg = OutputConfig.from_input_fasta(
            input_fa="/data/refs/sample.fasta",
            out_dir=None,
            out_base="my_output",
            suffix="_sim",
        )
        assert cfg.out_dir == Path("/data/refs")
        assert cfg.out_base == "my_output"


class TestEnsureDir:
    """Test ensure_dir creates directories."""

    def test_ensure_dir_creates_nested(self, tmp_path: Path) -> None:
        nested = tmp_path / "a" / "b" / "c"
        cfg = OutputConfig(out_dir=nested, out_base="test")
        assert not nested.exists()
        cfg.ensure_dir()
        assert nested.is_dir()

    def test_ensure_dir_existing(self, tmp_path: Path) -> None:
        cfg = OutputConfig(out_dir=tmp_path, out_base="test")
        cfg.ensure_dir()  # should not raise
        assert tmp_path.is_dir()
