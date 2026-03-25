"""Tests for read source tracking data models, coordinate map, and tracker."""

from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path

import pytest

from muc_one_up.read_simulator.source_tracking import (
    ReadOrigin,
    ReadSourceTracker,
    RepeatRegion,
    build_coordinate_map,
)


class TestRepeatRegion:
    """Tests for RepeatRegion dataclass."""

    def test_basic_creation(self):
        region = RepeatRegion(
            index=1,
            repeat_type="X",
            start=100,
            end=160,
            is_mutated=False,
            mutation_name=None,
        )
        assert region.index == 1
        assert region.repeat_type == "X"
        assert region.start == 100
        assert region.end == 160
        assert region.length == 60

    def test_mutated_region(self):
        region = RepeatRegion(
            index=3,
            repeat_type="Xm",
            start=200,
            end=261,
            is_mutated=True,
            mutation_name="dupC",
        )
        assert region.is_mutated is True
        assert region.mutation_name == "dupC"
        assert region.length == 61


class TestBuildCoordinateMap:
    """Tests for building coordinate maps from repeat chains."""

    @pytest.fixture()
    def repeats_dict(self):
        return {
            "1": "AACCCTCCC",
            "2": "AACCCTCCC",
            "X": "AACCCTCCCAACCCTCCC",
            "7": "AACCCTCCC",
            "8": "AACCCTCCC",
            "9": "AACCCTCCC",
        }

    def test_simple_chain(self, repeats_dict):
        chain = ["1", "2", "X", "7", "8", "9"]
        coord_map = build_coordinate_map(
            haplotype=1,
            chain=chain,
            repeats_dict=repeats_dict,
            left_const_len=50,
            mutation_positions=[],
            mutation_name=None,
            snp_info=[],
        )
        assert coord_map.haplotype == 1
        assert len(coord_map.regions) == 6
        assert coord_map.regions[0].start == 50
        assert coord_map.regions[0].end == 59
        assert coord_map.regions[1].start == 59
        assert coord_map.regions[1].end == 68
        assert coord_map.regions[2].start == 68
        assert coord_map.regions[2].end == 86
        assert coord_map.vntr_start == 50
        assert coord_map.vntr_end == 86 + 9 + 9 + 9

    def test_with_mutation_markers(self, repeats_dict):
        chain = ["1", "2", "Xm", "7", "8", "9"]
        coord_map = build_coordinate_map(
            haplotype=1,
            chain=chain,
            repeats_dict=repeats_dict,
            left_const_len=50,
            mutation_positions=[(1, 3)],
            mutation_name="dupC",
            snp_info=[],
        )
        assert coord_map.regions[2].repeat_type == "Xm"
        assert coord_map.regions[2].is_mutated is True
        assert coord_map.regions[2].mutation_name == "dupC"
        # Non-mutated regions
        assert coord_map.regions[0].is_mutated is False
        assert coord_map.regions[0].mutation_name is None

    def test_with_snps(self, repeats_dict):
        chain = ["1", "2", "X"]
        snps = [{"position": 55, "ref_base": "A", "alt_base": "T"}]
        coord_map = build_coordinate_map(
            haplotype=1,
            chain=chain,
            repeats_dict=repeats_dict,
            left_const_len=50,
            mutation_positions=[],
            mutation_name=None,
            snp_info=snps,
        )
        assert len(coord_map.snp_positions) == 1
        assert coord_map.snp_positions[0].position == 55
        # Position 55 falls in first repeat (50-59)
        assert coord_map.snp_positions[0].repeat_index == 1

    def test_snp_outside_vntr(self, repeats_dict):
        chain = ["1"]
        snps = [{"position": 10, "ref_base": "G", "alt_base": "C"}]
        coord_map = build_coordinate_map(
            haplotype=1,
            chain=chain,
            repeats_dict=repeats_dict,
            left_const_len=50,
            mutation_positions=[],
            mutation_name=None,
            snp_info=snps,
        )
        assert coord_map.snp_positions[0].repeat_index is None

    def test_empty_chain(self, repeats_dict):
        coord_map = build_coordinate_map(
            haplotype=1,
            chain=[],
            repeats_dict=repeats_dict,
            left_const_len=50,
            mutation_positions=[],
            mutation_name=None,
            snp_info=[],
        )
        assert len(coord_map.regions) == 0
        assert coord_map.vntr_start == 50
        assert coord_map.vntr_end == 50

    def test_single_repeat(self, repeats_dict):
        chain = ["X"]
        coord_map = build_coordinate_map(
            haplotype=1,
            chain=chain,
            repeats_dict=repeats_dict,
            left_const_len=0,
            mutation_positions=[],
            mutation_name=None,
            snp_info=[],
        )
        assert len(coord_map.regions) == 1
        assert coord_map.regions[0].start == 0
        assert coord_map.regions[0].end == 18
        assert coord_map.vntr_start == 0
        assert coord_map.vntr_end == 18

    def test_unknown_repeat_symbol_raises(self, repeats_dict):
        chain = ["1", "Z"]  # Z is not in repeats_dict
        with pytest.raises(ValueError, match="Unknown repeat symbol 'Z'"):
            build_coordinate_map(
                haplotype=1,
                chain=chain,
                repeats_dict=repeats_dict,
                left_const_len=50,
                mutation_positions=[],
                mutation_name=None,
                snp_info=[],
            )


class TestReadSourceTracker:
    """Tests for ReadSourceTracker annotation and manifest writing."""

    @pytest.fixture()
    def repeats_dict(self):
        return {
            "1": "AACCCTCCC",
            "2": "AACCCTCCC",
            "X": "AACCCTCCCAACCCTCCC",
            "7": "AACCCTCCC",
            "8": "AACCCTCCC",
            "9": "AACCCTCCC",
        }

    @pytest.fixture()
    def tracker(self, repeats_dict):
        return ReadSourceTracker(
            repeat_chains={
                1: ["1", "2", "Xm", "7", "8", "9"],
                2: ["1", "2", "X", "7", "8", "9"],
            },
            repeats_dict=repeats_dict,
            left_const_len=50,
            mutation_positions=[(1, 3)],
            mutation_name="dupC",
            snp_info={0: [{"position": 55, "ref_base": "A", "alt_base": "T"}]},
        )

    def test_coordinate_maps_built(self, tracker):
        assert 1 in tracker.coordinate_maps
        assert 2 in tracker.coordinate_maps
        assert len(tracker.coordinate_maps[1].regions) == 6
        assert len(tracker.coordinate_maps[2].regions) == 6

    def test_annotation_overlaps_vntr(self, tracker):
        origins = [
            ReadOrigin(
                read_id="read1",
                haplotype=1,
                ref_start=60,
                ref_end=90,
                strand="+",
            ),
        ]
        annotated = list(tracker.annotate_reads(origins))
        assert len(annotated) == 1
        assert annotated[0].overlaps_vntr is True

    def test_annotation_outside_vntr(self, tracker):
        origins = [
            ReadOrigin(
                read_id="read2",
                haplotype=1,
                ref_start=0,
                ref_end=40,
                strand="+",
            ),
        ]
        annotated = list(tracker.annotate_reads(origins))
        assert annotated[0].overlaps_vntr is False
        assert annotated[0].repeat_units == ""

    def test_annotation_overlaps_mutation(self, tracker):
        # Xm repeat on haplotype 1 is at positions 68-86
        origins = [
            ReadOrigin(
                read_id="read3",
                haplotype=1,
                ref_start=70,
                ref_end=80,
                strand="+",
            ),
        ]
        annotated = list(tracker.annotate_reads(origins))
        assert annotated[0].overlaps_mutation is True
        assert annotated[0].mutation_name == "dupC"

    def test_annotation_no_mutation_on_haplotype2(self, tracker):
        # Haplotype 2 has no mutation
        origins = [
            ReadOrigin(
                read_id="read4",
                haplotype=2,
                ref_start=70,
                ref_end=80,
                strand="+",
            ),
        ]
        annotated = list(tracker.annotate_reads(origins))
        assert annotated[0].overlaps_mutation is False

    def test_annotation_overlaps_snp(self, tracker):
        # SNP at position 55 on haplotype 1
        origins = [
            ReadOrigin(
                read_id="read5",
                haplotype=1,
                ref_start=50,
                ref_end=60,
                strand="+",
            ),
        ]
        annotated = list(tracker.annotate_reads(origins))
        assert annotated[0].overlaps_snp is True

    def test_annotation_no_snp_overlap(self, tracker):
        origins = [
            ReadOrigin(
                read_id="read6",
                haplotype=1,
                ref_start=70,
                ref_end=80,
                strand="+",
            ),
        ]
        annotated = list(tracker.annotate_reads(origins))
        assert annotated[0].overlaps_snp is False

    def test_annotation_repeat_units(self, tracker):
        # Read spanning repeats 2 (59-68) and 3 (68-86) on haplotype 1
        origins = [
            ReadOrigin(
                read_id="read7",
                haplotype=1,
                ref_start=59,
                ref_end=86,
                strand="+",
            ),
        ]
        annotated = list(tracker.annotate_reads(origins))
        assert annotated[0].repeat_units == "2,3"

    def test_annotation_unknown_haplotype(self, tracker):
        origins = [
            ReadOrigin(
                read_id="read8",
                haplotype=99,
                ref_start=0,
                ref_end=100,
                strand="+",
            ),
        ]
        annotated = list(tracker.annotate_reads(origins))
        assert annotated[0].overlaps_vntr is False

    def test_write_manifest(self, tracker, tmp_path):
        origins = [
            ReadOrigin(
                read_id="read1",
                haplotype=1,
                ref_start=60,
                ref_end=90,
                strand="+",
            ),
            ReadOrigin(
                read_id="read2",
                haplotype=2,
                ref_start=0,
                ref_end=40,
                strand="-",
            ),
        ]
        annotated = list(tracker.annotate_reads(origins))
        manifest_path = str(tmp_path / "test_read_manifest.tsv.gz")
        tracker.write_manifest(annotated, manifest_path)

        assert Path(manifest_path).exists()
        with gzip.open(manifest_path, "rt") as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = list(reader)
        assert len(rows) == 2
        assert rows[0]["read_id"] == "read1"
        assert rows[0]["overlaps_vntr"] == "true"
        assert rows[0]["haplotype"] == "1"
        assert rows[1]["read_id"] == "read2"
        assert rows[1]["overlaps_vntr"] == "false"

    def test_write_coordinate_map(self, tracker, tmp_path):
        coord_path = str(tmp_path / "test_repeat_coordinates.tsv")
        tracker.write_coordinate_map(coord_path)
        assert Path(coord_path).exists()
        with open(coord_path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = list(reader)
        # 6 repeats x 2 haplotypes = 12 rows
        assert len(rows) == 12
        assert rows[0]["haplotype"] == "1"
        assert rows[0]["repeat_type"] == "1"
        # Check SNP columns exist
        assert "snp_count" in rows[0]
        assert "snp_positions" in rows[0]

    def test_manifest_boolean_format(self, tracker, tmp_path):
        origins = [
            ReadOrigin(
                read_id="r1",
                haplotype=1,
                ref_start=70,
                ref_end=80,
                strand="+",
            ),
        ]
        annotated = list(tracker.annotate_reads(origins))
        path = str(tmp_path / "bool_test.tsv.gz")
        tracker.write_manifest(annotated, path)

        with gzip.open(path, "rt") as f:
            reader = csv.DictReader(f, delimiter="\t")
            row = next(reader)
        # Booleans must be lowercase
        assert row["overlaps_vntr"] in ("true", "false")
        assert row["overlaps_mutation"] in ("true", "false")
        assert row["overlaps_snp"] in ("true", "false")

    def test_manifest_empty_repeat_units_as_dot(self, tracker, tmp_path):
        origins = [
            ReadOrigin(
                read_id="r1",
                haplotype=1,
                ref_start=0,
                ref_end=10,
                strand="+",
            ),
        ]
        annotated = list(tracker.annotate_reads(origins))
        path = str(tmp_path / "dot_test.tsv.gz")
        tracker.write_manifest(annotated, path)

        with gzip.open(path, "rt") as f:
            reader = csv.DictReader(f, delimiter="\t")
            row = next(reader)
        assert row["repeat_units"] == "."

    def test_manifest_empty_mutation_name_as_dot(self, tracker, tmp_path):
        # Read outside mutation region on haplotype 2 (no mutation)
        origins = [
            ReadOrigin(
                read_id="r1",
                haplotype=2,
                ref_start=70,
                ref_end=80,
                strand="+",
            ),
        ]
        annotated = list(tracker.annotate_reads(origins))
        path = str(tmp_path / "mut_dot_test.tsv.gz")
        tracker.write_manifest(annotated, path)

        with gzip.open(path, "rt") as f:
            reader = csv.DictReader(f, delimiter="\t")
            row = next(reader)
        assert row["mutation_name"] == "."


class TestTrackerFromCompanionFiles:
    """Tests for ReadSourceTracker.from_companion_files()."""

    def test_reconstruct_from_stats(self, tmp_path):
        stats = {
            "haplotypes": {
                "haplotype_1": {"repeat_chain": "1-2-X-7-8-9"},
                "haplotype_2": {"repeat_chain": "1-2-X-7-8-9"},
            },
            "mutation_details": {
                "mutation_name": "dupC",
                "targets": [[1, 3]],
            },
            "snp_info": {},
            "config": {
                "repeats": {
                    "1": "AACCCTCCC",
                    "2": "AACCCTCCC",
                    "X": "AACCCTCCCAACCCTCCC",
                    "7": "AACCCTCCC",
                    "8": "AACCCTCCC",
                    "9": "AACCCTCCC",
                },
                "constants": {"hg38": {"left": "A" * 50}},
            },
            "reference_assembly": "hg38",
        }
        stats_path = str(tmp_path / "test.simulation_stats.json")
        with open(stats_path, "w") as f:
            json.dump(stats, f)

        tracker = ReadSourceTracker.from_companion_files(stats_path)
        assert tracker is not None
        assert 1 in tracker.coordinate_maps
        assert 2 in tracker.coordinate_maps
        assert len(tracker.coordinate_maps[1].regions) == 6

    def test_missing_file_returns_none(self):
        result = ReadSourceTracker.from_companion_files("/nonexistent/path.json")
        assert result is None

    def test_missing_haplotypes_returns_none(self, tmp_path):
        stats = {"config": {"repeats": {}}}
        stats_path = str(tmp_path / "bad_stats.json")
        with open(stats_path, "w") as f:
            json.dump(stats, f)
        result = ReadSourceTracker.from_companion_files(stats_path)
        assert result is None
