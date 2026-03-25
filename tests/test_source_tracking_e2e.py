"""End-to-end integration tests for read source tracking.

These tests exercise the full flow: create realistic simulator output files,
run the platform-specific parsers, feed results through ReadSourceTracker
annotation, write manifest/coordinate map files, and verify the output
contents are correct and consistent.
"""

from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path

import pytest

from muc_one_up.read_simulator.parsers.illumina_parser import (
    parse_illumina_reads,
    write_fragment_origins,
)
from muc_one_up.read_simulator.parsers.ont_parser import parse_nanosim_reads
from muc_one_up.read_simulator.parsers.pacbio_parser import parse_pacbio_reads
from muc_one_up.read_simulator.source_tracking import (
    ReadOrigin,
    ReadSourceTracker,
)

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

REPEATS_DICT = {
    "1": "AACCCTCCC",  # 9 bp
    "2": "AACCCTCCC",  # 9 bp
    "X": "AACCCTCCCAACCCTCCC",  # 18 bp
    "7": "AACCCTCCC",  # 9 bp
    "8": "AACCCTCCC",  # 9 bp
    "9": "AACCCTCCC",  # 9 bp
}

LEFT_CONST_LEN = 50
# VNTR starts at 50, regions: 1(50-59), 2(59-68), X(68-86), 7(86-95), 8(95-104), 9(104-113)
# Total haplotype reference length ~ 113 + right_const


@pytest.fixture()
def tracker_with_mutation():
    """Tracker with mutation on haplotype 1, repeat 3 (Xm), SNP at pos 55."""
    return ReadSourceTracker(
        repeat_chains={
            1: ["1", "2", "Xm", "7", "8", "9"],
            2: ["1", "2", "X", "7", "8", "9"],
        },
        repeats_dict=REPEATS_DICT,
        left_const_len=LEFT_CONST_LEN,
        mutation_positions=[(1, 3)],
        mutation_name="dupC",
        snp_info={0: [{"position": 55, "ref_base": "A", "alt_base": "T"}]},
    )


@pytest.fixture()
def tracker_no_mutation():
    """Tracker with no mutations or SNPs."""
    return ReadSourceTracker(
        repeat_chains={
            1: ["1", "2", "X", "7", "8", "9"],
            2: ["1", "2", "X", "7", "8", "9"],
        },
        repeats_dict=REPEATS_DICT,
        left_const_len=LEFT_CONST_LEN,
    )


# ---------------------------------------------------------------------------
# ONT Parser E2E
# ---------------------------------------------------------------------------


class TestONTEndToEnd:
    """End-to-end test: NanoSim FASTQ -> parser -> tracker -> manifest."""

    def _write_nanosim_fastq(self, path: Path, reads: list[tuple[str, str]]):
        """Write a minimal NanoSim-style FASTQ file.

        Args:
            reads: List of (read_name, sequence) tuples.
        """
        with open(path, "w") as f:
            for name, seq in reads:
                f.write(f"@{name}\n")
                f.write(f"{seq}\n")
                f.write("+\n")
                f.write("I" * len(seq) + "\n")

    def test_full_ont_pipeline(self, tracker_with_mutation, tmp_path):
        """Parse NanoSim reads, annotate, write manifest, verify contents."""
        # Create realistic NanoSim FASTQ with reads at known positions
        reads = [
            # Read from haplotype_1 at pos 70 (inside Xm repeat 68-86), forward, 20bp
            ("haplotype_1_70_aligned_0_F_5_10_5", "A" * 20),
            # Read from haplotype_1 at pos 50 (start of VNTR, repeat 1), reverse, 15bp
            ("haplotype_1_50_aligned_1_R_3_9_3", "C" * 15),
            # Read from haplotype_2 at pos 70 (inside X repeat, no mutation), forward, 25bp
            ("haplotype_2_70_aligned_2_F_5_15_5", "G" * 25),
            # Read from haplotype_1 at pos 10 (outside VNTR, in left constant), forward, 30bp
            ("haplotype_1_10_aligned_3_F_5_20_5", "T" * 30),
        ]
        fastq_path = tmp_path / "nanosim_reads.fastq"
        self._write_nanosim_fastq(fastq_path, reads)

        # Parse
        origins = parse_nanosim_reads(str(fastq_path))
        assert len(origins) == 4

        # Verify parsed origins
        assert origins[0].haplotype == 1
        assert origins[0].ref_start == 70
        assert origins[0].ref_end == 90  # 70 + 20
        assert origins[0].strand == "+"

        assert origins[1].haplotype == 1
        assert origins[1].ref_start == 50
        assert origins[1].strand == "-"

        assert origins[2].haplotype == 2
        assert origins[2].ref_start == 70

        assert origins[3].haplotype == 1
        assert origins[3].ref_start == 10
        assert origins[3].ref_end == 40  # 10 + 30

        # Annotate
        annotated = list(tracker_with_mutation.annotate_reads(origins))
        assert len(annotated) == 4

        # Read 0: haplotype 1, pos 70-90, overlaps Xm(68-86) and 7(86-95)
        r0 = annotated[0]
        assert r0.overlaps_vntr is True
        assert r0.overlaps_mutation is True
        assert r0.mutation_name == "dupC"
        assert "3" in r0.repeat_units  # Xm is index 3
        assert "4" in r0.repeat_units  # 7 is index 4

        # Read 1: haplotype 1, pos 50-65, overlaps repeat 1(50-59) and 2(59-68)
        r1 = annotated[1]
        assert r1.overlaps_vntr is True
        assert r1.overlaps_mutation is False
        assert r1.overlaps_snp is True  # SNP at pos 55
        assert "1" in r1.repeat_units
        assert "2" in r1.repeat_units

        # Read 2: haplotype 2, pos 70-95, overlaps X(68-86) - no mutation on hap 2
        r2 = annotated[2]
        assert r2.overlaps_vntr is True
        assert r2.overlaps_mutation is False
        assert r2.mutation_name == ""

        # Read 3: haplotype 1, pos 10-40, entirely in left constant
        r3 = annotated[3]
        assert r3.overlaps_vntr is False
        assert r3.overlaps_mutation is False
        assert r3.overlaps_snp is False
        assert r3.repeat_units == ""

        # Write manifest and verify
        manifest_path = str(tmp_path / "ont_read_manifest.tsv.gz")
        tracker_with_mutation.write_manifest(annotated, manifest_path)

        assert Path(manifest_path).exists()
        with gzip.open(manifest_path, "rt") as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = list(reader)

        assert len(rows) == 4
        # Verify manifest columns match annotation
        assert rows[0]["overlaps_mutation"] == "true"
        assert rows[0]["mutation_name"] == "dupC"
        assert rows[1]["overlaps_snp"] == "true"
        assert rows[2]["overlaps_mutation"] == "false"
        assert rows[3]["overlaps_vntr"] == "false"
        assert rows[3]["repeat_units"] == "."

    def test_ont_split_sim_haplotype_override(self, tracker_no_mutation, tmp_path):
        """Test split-sim mode where haplotype_map overrides parsed haplotype."""
        reads = [
            ("haplotype_1_60_aligned_0_F_5_10_5", "A" * 15),
        ]
        fastq_path = tmp_path / "split_hap2.fastq"
        self._write_nanosim_fastq(fastq_path, reads)

        # Parse with haplotype_map=2 (split-sim: this file is haplotype 2)
        origins = parse_nanosim_reads(str(fastq_path), haplotype_map=2)
        assert len(origins) == 1
        assert origins[0].haplotype == 2  # overridden, not 1

    def test_ont_empty_fastq(self, tmp_path):
        """Parser handles empty FASTQ gracefully."""
        fastq_path = tmp_path / "empty.fastq"
        fastq_path.write_text("")
        origins = parse_nanosim_reads(str(fastq_path))
        assert origins == []


# ---------------------------------------------------------------------------
# PacBio Parser E2E
# ---------------------------------------------------------------------------


class TestPacBioEndToEnd:
    """End-to-end test: pbsim3 MAF -> parser -> tracker -> manifest."""

    def _write_maf(self, path: Path, blocks: list[tuple[str, str]]):
        """Write a minimal MAF file.

        Args:
            blocks: List of (ref_s_line, read_s_line) tuples.
        """
        with open(path, "w") as f:
            f.write("##maf version=1 scoring=pbsim3\n\n")
            for ref_line, read_line in blocks:
                f.write("a\n")
                f.write(f"{ref_line}\n")
                f.write(f"{read_line}\n")
                f.write("\n")

    def test_full_pacbio_pipeline(self, tracker_with_mutation, tmp_path):
        """Parse MAF files, annotate, write manifest, verify contents."""
        # Create MAF files for haplotype 1 and 2
        # MAF s-line format: s name start size strand srcSize sequence
        hap1_maf = tmp_path / "sd_0001.maf"
        self._write_maf(
            hap1_maf,
            [
                # Read at pos 65, spanning 30bp (crosses Xm repeat)
                (
                    "s ref 65 30 + 200 AACCCTCCCAACCCTCCCAACCCTCCCAAC",
                    "s S1_1 0 30 + 30 AACCCTCCCAACCCTCCCAACCCTCCCAAC",
                ),
                # Read at pos 0, spanning 45bp (in left constant only)
                (
                    "s ref 0 45 + 200 " + "A" * 45,
                    "s S1_2 0 45 + 45 " + "A" * 45,
                ),
            ],
        )

        hap2_maf = tmp_path / "sd_0002.maf"
        self._write_maf(
            hap2_maf,
            [
                # Read at pos 86, spanning 20bp (repeat 7 region)
                (
                    "s ref 86 20 + 200 " + "G" * 20,
                    "s S2_1 0 20 + 20 " + "G" * 20,
                ),
            ],
        )

        # Parse haplotype 1
        origins_hap1 = parse_pacbio_reads([str(hap1_maf)], haplotype_index=1)
        assert len(origins_hap1) == 2
        assert origins_hap1[0].ref_start == 65
        assert origins_hap1[0].ref_end == 95
        assert origins_hap1[0].haplotype == 1
        assert origins_hap1[1].ref_start == 0
        assert origins_hap1[1].ref_end == 45

        # Parse haplotype 2
        origins_hap2 = parse_pacbio_reads([str(hap2_maf)], haplotype_index=2)
        assert len(origins_hap2) == 1
        assert origins_hap2[0].haplotype == 2
        assert origins_hap2[0].ref_start == 86

        # Combine and annotate
        all_origins = origins_hap1 + origins_hap2
        annotated = list(tracker_with_mutation.annotate_reads(all_origins))
        assert len(annotated) == 3

        # Read 0: hap1, pos 65-95, overlaps repeats 2(59-68), Xm(68-86), 7(86-95)
        r0 = annotated[0]
        assert r0.overlaps_vntr is True
        assert r0.overlaps_mutation is True
        assert r0.mutation_name == "dupC"
        assert "2" in r0.repeat_units
        assert "3" in r0.repeat_units  # Xm
        assert "4" in r0.repeat_units  # 7

        # Read 1: hap1, pos 0-45, entirely in left constant
        r1 = annotated[1]
        assert r1.overlaps_vntr is False
        assert r1.overlaps_mutation is False

        # Read 2: hap2, pos 86-106, overlaps repeat 7(86-95) and 8(95-104)
        r2 = annotated[2]
        assert r2.overlaps_vntr is True
        assert r2.overlaps_mutation is False  # no mutation on hap2

        # Write and verify manifest
        manifest_path = str(tmp_path / "pacbio_read_manifest.tsv.gz")
        tracker_with_mutation.write_manifest(annotated, manifest_path)
        with gzip.open(manifest_path, "rt") as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = list(reader)
        assert len(rows) == 3
        assert rows[0]["mutation_name"] == "dupC"

    def test_pacbio_missing_maf(self, tmp_path):
        """Parser returns empty list for missing MAF file."""
        origins = parse_pacbio_reads([str(tmp_path / "nonexistent.maf")], haplotype_index=1)
        assert origins == []

    def test_pacbio_empty_maf(self, tmp_path):
        """Parser handles empty MAF gracefully."""
        maf = tmp_path / "empty.maf"
        maf.write_text("##maf version=1\n\n")
        origins = parse_pacbio_reads([str(maf)], haplotype_index=1)
        assert origins == []


# ---------------------------------------------------------------------------
# Illumina Parser E2E
# ---------------------------------------------------------------------------


class TestIlluminaEndToEnd:
    """End-to-end test: fragment sidecar -> parser -> tracker -> manifest."""

    def _write_fastq_gz(self, path: Path, read_names: list[str]):
        """Write a gzipped FASTQ with given read names."""
        with gzip.open(path, "wt") as f:
            for name in read_names:
                f.write(f"@{name}\n")
                f.write("ACGT" * 10 + "\n")
                f.write("+\n")
                f.write("I" * 40 + "\n")

    def test_full_illumina_pipeline(self, tracker_with_mutation, tmp_path):
        """Write sidecar, parse with FASTQ, annotate, write manifest."""
        # Write fragment origins sidecar
        fragments = [
            {
                "fragment_index": 0,
                "chrom": "haplotype_1",
                "fstart": 55,
                "fend": 90,
                "strand": "+",
            },
            {
                "fragment_index": 1,
                "chrom": "haplotype_2",
                "fstart": 86,
                "fend": 110,
                "strand": "-",
            },
            {
                "fragment_index": 2,
                "chrom": "haplotype_1",
                "fstart": 5,
                "fend": 35,
                "strand": "+",
            },
        ]
        sidecar_path = str(tmp_path / "fragment_origins.tsv")
        write_fragment_origins(fragments, sidecar_path)

        # Verify sidecar was written correctly
        with open(sidecar_path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = list(reader)
        assert len(rows) == 3
        assert rows[0]["chrom"] == "haplotype_1"
        assert rows[0]["fstart"] == "55"

        # Write a gzipped FASTQ with matching read names
        fastq_path = tmp_path / "reads_R1.fastq.gz"
        self._write_fastq_gz(
            fastq_path,
            ["read_pair_0/1", "read_pair_1/1", "read_pair_2/1"],
        )

        # Parse
        seq_name_to_haplotype = {"haplotype_1": 1, "haplotype_2": 2}
        origins = parse_illumina_reads(
            sidecar_path,
            str(fastq_path),
            seq_name_to_haplotype,
        )
        assert len(origins) == 3

        # Verify parsed origins
        assert origins[0].haplotype == 1
        assert origins[0].ref_start == 55
        assert origins[0].ref_end == 90
        assert origins[0].read_id == "read_pair_0/1"  # full ID preserved

        assert origins[1].haplotype == 2
        assert origins[1].ref_start == 86
        assert origins[1].strand == "-"

        assert origins[2].haplotype == 1
        assert origins[2].ref_start == 5
        assert origins[2].ref_end == 35

        # Annotate
        annotated = list(tracker_with_mutation.annotate_reads(origins))

        # Fragment 0: hap1, pos 55-90, overlaps repeats 1(50-59), 2(59-68), Xm(68-86), 7(86-95)
        r0 = annotated[0]
        assert r0.overlaps_vntr is True
        assert r0.overlaps_mutation is True
        assert r0.mutation_name == "dupC"
        assert r0.overlaps_snp is True  # SNP at 55

        # Fragment 1: hap2, pos 86-110, overlaps 7(86-95), 8(95-104), 9(104-113)
        r1 = annotated[1]
        assert r1.overlaps_vntr is True
        assert r1.overlaps_mutation is False  # no mutation on hap 2

        # Fragment 2: hap1, pos 5-35, in left constant
        r2 = annotated[2]
        assert r2.overlaps_vntr is False

        # Write manifest
        manifest_path = str(tmp_path / "illumina_read_manifest.tsv.gz")
        tracker_with_mutation.write_manifest(annotated, manifest_path)

        with gzip.open(manifest_path, "rt") as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = list(reader)
        assert len(rows) == 3

    def test_illumina_without_fastq(self, tmp_path):
        """Parser works without FASTQ, using fragment_index as read_id."""
        fragments = [
            {
                "fragment_index": 0,
                "chrom": "haplotype_1",
                "fstart": 60,
                "fend": 80,
                "strand": "+",
            },
        ]
        sidecar_path = str(tmp_path / "fragment_origins.tsv")
        write_fragment_origins(fragments, sidecar_path)

        origins = parse_illumina_reads(
            sidecar_path,
            None,  # no FASTQ
            {"haplotype_1": 1},
        )
        assert len(origins) == 1
        assert origins[0].read_id == "fragment_0"

    def test_illumina_missing_sidecar(self, tmp_path):
        """Parser returns empty list when sidecar file is missing."""
        origins = parse_illumina_reads(
            str(tmp_path / "nonexistent.tsv"),
            None,
            {},
        )
        assert origins == []


# ---------------------------------------------------------------------------
# Coordinate Map Persistence E2E
# ---------------------------------------------------------------------------


class TestCoordinateMapEndToEnd:
    """Test writing and verifying coordinate map TSV files."""

    def test_coordinate_map_with_mutations_and_snps(self, tracker_with_mutation, tmp_path):
        """Write coordinate map and verify all columns are correct."""
        coord_path = str(tmp_path / "repeat_coordinates.tsv")
        tracker_with_mutation.write_coordinate_map(coord_path)

        with open(coord_path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = list(reader)

        # 6 repeats x 2 haplotypes = 12 rows
        assert len(rows) == 12

        # Verify haplotype 1 entries
        hap1_rows = [r for r in rows if r["haplotype"] == "1"]
        assert len(hap1_rows) == 6

        # Xm repeat (index 3, haplotype 1) should be mutated
        xm_row = next(r for r in hap1_rows if r["index"] == "3")
        assert xm_row["is_mutated"] == "true"
        assert xm_row["mutation_name"] == "dupC"
        assert xm_row["repeat_type"] == "Xm"

        # First repeat (index 1, haplotype 1) should have SNP
        r1_row = next(r for r in hap1_rows if r["index"] == "1")
        assert r1_row["snp_count"] == "1"
        assert r1_row["snp_positions"] == "55"

        # Haplotype 2 X repeat should NOT be mutated
        hap2_rows = [r for r in rows if r["haplotype"] == "2"]
        x_hap2 = next(r for r in hap2_rows if r["index"] == "3")
        assert x_hap2["is_mutated"] == "false"
        assert x_hap2["mutation_name"] == "."


# ---------------------------------------------------------------------------
# from_companion_files E2E
# ---------------------------------------------------------------------------


class TestCompanionFilesEndToEnd:
    """Test reconstructing tracker from simulation_stats.json and using it."""

    def test_full_roundtrip(self, tmp_path):
        """Create stats JSON, reconstruct tracker, annotate reads, write manifest."""
        stats = {
            "haplotypes": {
                "haplotype_1": {"repeat_chain": "1-2-Xm-7-8-9"},
                "haplotype_2": {"repeat_chain": "1-2-X-7-8-9"},
            },
            "mutation_details": {
                "mutation_name": "dupC",
                "targets": [[1, 3]],
            },
            "snp_info": {
                "1": [{"position": 55, "ref_base": "A", "alt_base": "T"}],
            },
            "config": {
                "repeats": REPEATS_DICT,
                "constants": {"hg38": {"left": "A" * LEFT_CONST_LEN}},
            },
            "reference_assembly": "hg38",
        }
        stats_path = str(tmp_path / "simulation_stats.json")
        with open(stats_path, "w") as f:
            json.dump(stats, f)

        # Reconstruct tracker
        tracker = ReadSourceTracker.from_companion_files(stats_path)
        assert tracker is not None
        assert 1 in tracker.coordinate_maps
        assert 2 in tracker.coordinate_maps

        # Verify coordinate maps match expectations
        hap1_map = tracker.coordinate_maps[1]
        assert len(hap1_map.regions) == 6
        assert hap1_map.vntr_start == LEFT_CONST_LEN
        assert hap1_map.regions[2].is_mutated is True  # Xm

        # Annotate a read and verify
        origins = [
            ReadOrigin(
                read_id="companion_read",
                haplotype=1,
                ref_start=70,
                ref_end=90,
                strand="+",
            ),
        ]
        annotated = list(tracker.annotate_reads(origins))
        assert len(annotated) == 1
        assert annotated[0].overlaps_mutation is True
        assert annotated[0].mutation_name == "dupC"

        # Write manifest and verify it's valid
        manifest_path = str(tmp_path / "companion_manifest.tsv.gz")
        tracker.write_manifest(annotated, manifest_path)
        with gzip.open(manifest_path, "rt") as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = list(reader)
        assert len(rows) == 1
        assert rows[0]["mutation_name"] == "dupC"

    def test_companion_missing_mutation_info(self, tmp_path):
        """Tracker works when mutation_details is absent (normal simulation)."""
        stats = {
            "haplotypes": {
                "haplotype_1": {"repeat_chain": "1-2-X-7-8-9"},
            },
            "mutation_details": {},
            "snp_info": {},
            "config": {
                "repeats": REPEATS_DICT,
                "constants": {"hg38": {"left": "A" * LEFT_CONST_LEN}},
            },
            "reference_assembly": "hg38",
        }
        stats_path = str(tmp_path / "normal_stats.json")
        with open(stats_path, "w") as f:
            json.dump(stats, f)

        tracker = ReadSourceTracker.from_companion_files(stats_path)
        assert tracker is not None

        origins = [
            ReadOrigin(read_id="r1", haplotype=1, ref_start=70, ref_end=85, strand="+"),
        ]
        annotated = list(tracker.annotate_reads(origins))
        assert annotated[0].overlaps_vntr is True
        assert annotated[0].overlaps_mutation is False


# ---------------------------------------------------------------------------
# Cross-Platform Consistency
# ---------------------------------------------------------------------------


class TestCrossPlatformConsistency:
    """Verify that the same read position produces identical annotations
    regardless of which parser produced the ReadOrigin."""

    def test_same_position_same_annotation(self, tracker_with_mutation, tmp_path):
        """Three parsers producing origins at the same position should yield
        identical annotation results."""
        # Create ONT read at pos 70, haplotype 1, 20bp
        ont_fastq = tmp_path / "ont.fastq"
        with open(ont_fastq, "w") as f:
            f.write("@haplotype_1_70_aligned_0_F_5_10_5\n")
            f.write("A" * 20 + "\n")
            f.write("+\n")
            f.write("I" * 20 + "\n")
        ont_origins = parse_nanosim_reads(str(ont_fastq))

        # Create PacBio MAF read at same position
        pb_maf = tmp_path / "pb.maf"
        with open(pb_maf, "w") as f:
            f.write("##maf version=1\n\n")
            f.write("a\n")
            f.write(f"s ref 70 20 + 200 {'A' * 20}\n")
            f.write(f"s pb_read 0 20 + 20 {'A' * 20}\n")
            f.write("\n")
        pb_origins = parse_pacbio_reads([str(pb_maf)], haplotype_index=1)

        # Create Illumina fragment at same position
        sidecar = tmp_path / "frags.tsv"
        write_fragment_origins(
            [
                {
                    "fragment_index": 0,
                    "chrom": "haplotype_1",
                    "fstart": 70,
                    "fend": 90,
                    "strand": "+",
                }
            ],
            str(sidecar),
        )
        ill_origins = parse_illumina_reads(str(sidecar), None, {"haplotype_1": 1})

        # All three should produce same region: hap1, 70-90
        assert ont_origins[0].ref_start == pb_origins[0].ref_start == ill_origins[0].ref_start
        assert ont_origins[0].ref_end == pb_origins[0].ref_end == ill_origins[0].ref_end

        # Annotate all three
        ont_ann = next(iter(tracker_with_mutation.annotate_reads(ont_origins)))
        pb_ann = next(iter(tracker_with_mutation.annotate_reads(pb_origins)))
        ill_ann = next(iter(tracker_with_mutation.annotate_reads(ill_origins)))

        # All three should have identical annotation fields
        assert ont_ann.overlaps_vntr == pb_ann.overlaps_vntr == ill_ann.overlaps_vntr == True  # noqa: E712
        assert (
            ont_ann.overlaps_mutation
            == pb_ann.overlaps_mutation
            == ill_ann.overlaps_mutation
            is True
        )
        assert ont_ann.mutation_name == pb_ann.mutation_name == ill_ann.mutation_name == "dupC"
        assert ont_ann.repeat_units == pb_ann.repeat_units == ill_ann.repeat_units
        assert ont_ann.overlaps_snp == pb_ann.overlaps_snp == ill_ann.overlaps_snp == False  # noqa: E712
