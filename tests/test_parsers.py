"""Tests for platform-specific read origin parsers."""

from __future__ import annotations

import csv
import gzip
from pathlib import Path

from muc_one_up.read_simulator.parsers.illumina_parser import (
    parse_illumina_reads,
    write_fragment_origins,
)
from muc_one_up.read_simulator.parsers.ont_parser import parse_nanosim_reads
from muc_one_up.read_simulator.parsers.pacbio_parser import parse_pacbio_reads


class TestONTParser:
    """Tests for NanoSim read name parser."""

    def test_parse_aligned_read(self, tmp_path):
        fastq = tmp_path / "reads.fastq"
        fastq.write_text(
            "@haplotype_1_500_aligned_0_F_10_100_5\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n"
        )
        origins = parse_nanosim_reads(str(fastq), haplotype_map=None)
        assert len(origins) == 1
        assert origins[0].read_id == "haplotype_1_500_aligned_0_F_10_100_5"
        assert origins[0].haplotype == 1
        assert origins[0].ref_start == 500
        assert origins[0].strand == "+"
        assert origins[0].ref_end == 500 + 16

    def test_parse_reverse_strand(self, tmp_path):
        fastq = tmp_path / "reads.fastq"
        fastq.write_text("@haplotype_2_1000_aligned_1_R_5_200_10\nACGTACGT\n+\nIIIIIIII\n")
        origins = parse_nanosim_reads(str(fastq), haplotype_map=None)
        assert origins[0].strand == "-"
        assert origins[0].haplotype == 2

    def test_parse_unaligned_read(self, tmp_path):
        fastq = tmp_path / "reads.fastq"
        fastq.write_text("@haplotype_1_0_unaligned_0_F_0_100_0\nACGT\n+\nIIII\n")
        origins = parse_nanosim_reads(str(fastq), haplotype_map=None)
        assert len(origins) == 1
        assert origins[0].ref_start == 0

    def test_diploid_split_sim_haplotype_map(self, tmp_path):
        fastq = tmp_path / "reads.fastq"
        fastq.write_text("@chr1_500_aligned_0_F_10_100_5\nACGTACGT\n+\nIIIIIIII\n")
        origins = parse_nanosim_reads(str(fastq), haplotype_map=2)
        assert origins[0].haplotype == 2

    def test_multiple_reads(self, tmp_path):
        fastq = tmp_path / "reads.fastq"
        fastq.write_text(
            "@haplotype_1_100_aligned_0_F_5_50_3\n"
            "ACGT\n"
            "+\n"
            "IIII\n"
            "@haplotype_1_200_aligned_1_R_3_80_2\n"
            "ACGTACGT\n"
            "+\n"
            "IIIIIIII\n"
        )
        origins = parse_nanosim_reads(str(fastq), haplotype_map=None)
        assert len(origins) == 2
        assert origins[0].ref_start == 100
        assert origins[1].ref_start == 200

    def test_empty_fastq(self, tmp_path):
        fastq = tmp_path / "reads.fastq"
        fastq.write_text("")
        origins = parse_nanosim_reads(str(fastq), haplotype_map=None)
        assert len(origins) == 0


class TestPacBioParser:
    """Tests for pbsim3 MAF parser."""

    def test_parse_maf_file(self, tmp_path):
        maf_content = (
            "a\n"
            "s haplotype_1 100 60 + 5000 " + "A" * 60 + "\n"
            "s S0_0 0 60 + 60 " + "A" * 60 + "\n"
            "\n"
            "a\n"
            "s haplotype_1 300 40 + 5000 " + "A" * 40 + "\n"
            "s S0_1 0 40 + 40 " + "A" * 40 + "\n"
            "\n"
        )
        maf_path = tmp_path / "test.maf"
        maf_path.write_text(maf_content)
        origins = parse_pacbio_reads(
            maf_paths=[str(maf_path)],
            haplotype_index=1,
            aligned_bam=None,
        )
        assert len(origins) == 2
        assert origins[0].haplotype == 1
        assert origins[0].ref_start == 100
        assert origins[0].ref_end == 160
        assert origins[0].read_id == "S0_0"
        assert origins[1].ref_start == 300
        assert origins[1].ref_end == 340

    def test_parse_multiple_maf_files(self, tmp_path):
        maf1 = tmp_path / "sd_0001.maf"
        maf1.write_text(
            "a\ns haplotype_1 100 50 + 5000 " + "A" * 50 + "\ns S0_0 0 50 + 50 " + "A" * 50 + "\n\n"
        )
        maf2 = tmp_path / "sd_0002.maf"
        maf2.write_text(
            "a\ns haplotype_2 200 30 + 4000 " + "A" * 30 + "\ns S1_0 0 30 + 30 " + "A" * 30 + "\n\n"
        )
        origins_hap1 = parse_pacbio_reads(
            maf_paths=[str(maf1)], haplotype_index=1, aligned_bam=None
        )
        origins_hap2 = parse_pacbio_reads(
            maf_paths=[str(maf2)], haplotype_index=2, aligned_bam=None
        )
        assert all(o.haplotype == 1 for o in origins_hap1)
        assert all(o.haplotype == 2 for o in origins_hap2)

    def test_empty_maf(self, tmp_path):
        maf_path = tmp_path / "empty.maf"
        maf_path.write_text("")
        origins = parse_pacbio_reads(
            maf_paths=[str(maf_path)],
            haplotype_index=1,
            aligned_bam=None,
        )
        assert len(origins) == 0

    def test_missing_maf(self, tmp_path):
        origins = parse_pacbio_reads(
            maf_paths=[str(tmp_path / "nonexistent.maf")],
            haplotype_index=1,
            aligned_bam=None,
        )
        assert len(origins) == 0

    def test_parse_gzipped_maf(self, tmp_path):
        maf_content = (
            "a\ns haplotype_1 100 60 + 5000 " + "A" * 60 + "\ns S0_0 0 60 + 60 " + "A" * 60 + "\n\n"
        )
        maf_gz_path = tmp_path / "test.maf.gz"
        with gzip.open(maf_gz_path, "wt") as f:
            f.write(maf_content)
        origins = parse_pacbio_reads(
            maf_paths=[str(maf_gz_path)],
            haplotype_index=1,
        )
        assert len(origins) == 1
        assert origins[0].ref_start == 100

    def test_deduplicate_maf_and_maf_gz(self, tmp_path):
        """When both .maf and .maf.gz exist, only parse once."""
        maf_content = (
            "a\ns ref 50 20 + 200 " + "A" * 20 + "\ns read_0 0 20 + 20 " + "A" * 20 + "\n\n"
        )
        maf_path = tmp_path / "sd_0001.maf"
        maf_path.write_text(maf_content)
        maf_gz_path = tmp_path / "sd_0001.maf.gz"
        with gzip.open(maf_gz_path, "wt") as f:
            f.write(maf_content)
        # Pass both — should only get 1 read, not 2
        origins = parse_pacbio_reads(
            maf_paths=[str(maf_path), str(maf_gz_path)],
            haplotype_index=1,
        )
        assert len(origins) == 1


class TestIlluminaParser:
    """Tests for Illumina fragment sidecar and read parser."""

    def test_write_fragment_origins(self, tmp_path):
        origins_path = str(tmp_path / "fragment_origins.tsv")
        fragments = [
            {
                "fragment_index": 1,
                "chrom": "haplotype_1",
                "fstart": 100,
                "fend": 250,
                "strand": "+",
            },
            {
                "fragment_index": 2,
                "chrom": "haplotype_2",
                "fstart": 300,
                "fend": 450,
                "strand": "-",
            },
        ]
        write_fragment_origins(fragments, origins_path)
        assert Path(origins_path).exists()
        with open(origins_path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = list(reader)
        assert len(rows) == 2
        assert rows[0]["fragment_index"] == "1"
        assert rows[0]["chrom"] == "haplotype_1"
        assert rows[0]["fstart"] == "100"

    def test_parse_illumina_reads_from_sidecar(self, tmp_path):
        origins_path = str(tmp_path / "fragment_origins.tsv")
        with open(origins_path, "w") as f:
            f.write("fragment_index\tchrom\tfstart\tfend\tstrand\n")
            f.write("1\thaplotype_1\t100\t250\t+\n")
            f.write("2\thaplotype_2\t300\t450\t-\n")

        fq1 = tmp_path / "reads_R1.fastq.gz"
        with gzip.open(str(fq1), "wt") as f:
            f.write("@read_0001/1\nACGT\n+\nIIII\n@read_0002/1\nACGT\n+\nIIII\n")

        seq_names = {"haplotype_1": 1, "haplotype_2": 2}
        origins = parse_illumina_reads(origins_path, str(fq1), seq_names)
        assert len(origins) == 2
        assert origins[0].haplotype == 1
        assert origins[0].ref_start == 100
        assert origins[0].ref_end == 250
        assert origins[0].read_id == "read_0001/1"
        assert origins[1].haplotype == 2
        assert origins[1].ref_start == 300

    def test_parse_without_fastq(self, tmp_path):
        origins_path = str(tmp_path / "fragment_origins.tsv")
        with open(origins_path, "w") as f:
            f.write("fragment_index\tchrom\tfstart\tfend\tstrand\n")
            f.write("1\thaplotype_1\t100\t250\t+\n")

        origins = parse_illumina_reads(origins_path, None, {"haplotype_1": 1})
        assert len(origins) == 1
        assert origins[0].read_id == "fragment_1"

    def test_empty_sidecar(self, tmp_path):
        origins_path = str(tmp_path / "fragment_origins.tsv")
        with open(origins_path, "w") as f:
            f.write("fragment_index\tchrom\tfstart\tfend\tstrand\n")
        origins = parse_illumina_reads(origins_path, None, {})
        assert len(origins) == 0

    def test_missing_sidecar(self, tmp_path):
        origins = parse_illumina_reads(str(tmp_path / "nonexistent.tsv"), None, {})
        assert len(origins) == 0
