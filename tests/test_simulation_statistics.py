"""Tests for muc_one_up.simulation_statistics module.

Tests cover:
- GC content calculation
- VNTR region extraction
- Repeat length calculation
- Repeat type counting
- Mutation detection
- Haplotype statistics generation
- Overall statistics aggregation
- Full simulation report generation
"""

import json
from pathlib import Path

import pytest

from muc_one_up.simulation_statistics import (
    compute_gc_content,
    count_repeat_types,
    extract_vntr_region,
    generate_haplotype_stats,
    generate_overall_stats,
    generate_simulation_statistics,
    get_mutation_details,
    get_repeat_lengths,
    write_statistics_report,
)


@pytest.mark.unit
class TestComputeGcContent:
    """Tests for compute_gc_content function."""

    def test_gc_content_50_percent(self):
        """Test sequence with 50% GC content."""
        seq = "ATCGATCG"  # 4 G/C out of 8 = 50%
        gc = compute_gc_content(seq)
        assert gc == 50.0

    def test_gc_content_100_percent(self):
        """Test sequence with 100% GC content."""
        seq = "GCGCGCGC"
        gc = compute_gc_content(seq)
        assert gc == 100.0

    def test_gc_content_0_percent(self):
        """Test sequence with 0% GC content."""
        seq = "ATATATAT"
        gc = compute_gc_content(seq)
        assert gc == 0.0

    def test_gc_content_case_insensitive(self):
        """Test that GC calculation is case-insensitive."""
        seq_upper = "ATCG"
        seq_lower = "atcg"
        assert compute_gc_content(seq_upper) == compute_gc_content(seq_lower)

    def test_gc_content_empty_sequence(self):
        """Test GC content of empty sequence."""
        gc = compute_gc_content("")
        assert gc == 0.0


@pytest.mark.unit
class TestExtractVntrRegion:
    """Tests for extract_vntr_region function."""

    def test_extract_vntr_with_flanks(self, minimal_config: dict):
        """Test extracting VNTR region with known flanks."""
        left = minimal_config["constants"]["hg38"]["left"]
        right = minimal_config["constants"]["hg38"]["right"]
        vntr = "ATCGATCGATCG"
        seq = left + vntr + right

        extracted = extract_vntr_region(seq, minimal_config)

        assert extracted == vntr

    def test_extract_vntr_without_flanks(self, minimal_config: dict):
        """Test extracting VNTR when flanks are not present."""
        seq = "ATCGATCGATCG"

        extracted = extract_vntr_region(seq, minimal_config)

        # Should return the sequence as-is if flanks not found
        assert extracted == seq

    def test_extract_vntr_empty_sequence(self, minimal_config: dict):
        """Test extracting VNTR from empty sequence."""
        extracted = extract_vntr_region("", minimal_config)
        assert extracted == ""


@pytest.mark.unit
class TestGetRepeatLengths:
    """Tests for get_repeat_lengths function."""

    def test_repeat_lengths_basic(self, minimal_config: dict):
        """Test getting repeat lengths for basic chain."""
        chain = ["1", "2", "X"]

        lengths = get_repeat_lengths(chain, minimal_config)

        expected = [
            len(minimal_config["repeats"]["1"]),
            len(minimal_config["repeats"]["2"]),
            len(minimal_config["repeats"]["X"]),
        ]
        assert lengths == expected

    def test_repeat_lengths_with_mutation_markers(self, minimal_config: dict):
        """Test getting repeat lengths with mutation markers."""
        chain = ["1", "2m", "X"]  # "2m" should be treated as "2"

        lengths = get_repeat_lengths(chain, minimal_config)

        expected = [
            len(minimal_config["repeats"]["1"]),
            len(minimal_config["repeats"]["2"]),
            len(minimal_config["repeats"]["X"]),
        ]
        assert lengths == expected

    def test_repeat_lengths_empty_chain(self, minimal_config: dict):
        """Test getting repeat lengths for empty chain."""
        lengths = get_repeat_lengths([], minimal_config)
        assert lengths == []


@pytest.mark.unit
class TestCountRepeatTypes:
    """Tests for count_repeat_types function."""

    def test_count_repeat_types_basic(self):
        """Test counting repeat types in chain."""
        chain = ["1", "2", "X", "X", "A", "2"]

        counts = count_repeat_types(chain)

        assert counts == {"1": 1, "2": 2, "X": 2, "A": 1}

    def test_count_repeat_types_with_mutation_markers(self):
        """Test counting ignores mutation markers."""
        chain = ["1", "2m", "X", "Xm", "A", "2"]

        counts = count_repeat_types(chain)

        # Should count "2m" as "2" and "Xm" as "X"
        assert counts == {"1": 1, "2": 2, "X": 2, "A": 1}

    def test_count_repeat_types_empty_chain(self):
        """Test counting repeat types in empty chain."""
        counts = count_repeat_types([])
        assert counts == {}

    def test_count_repeat_types_single_repeat(self):
        """Test counting single repeat type."""
        chain = ["X", "X", "X", "X"]

        counts = count_repeat_types(chain)

        assert counts == {"X": 4}


@pytest.mark.unit
class TestGetMutationDetails:
    """Tests for get_mutation_details function."""

    def test_mutation_details_basic(self):
        """Test getting mutation details."""
        chain = ["1", "2m", "X", "Am", "B"]

        mutations = get_mutation_details(chain)

        assert len(mutations) == 2
        assert mutations[0] == {"position": 2, "repeat": "2"}
        assert mutations[1] == {"position": 4, "repeat": "A"}

    def test_mutation_details_no_mutations(self):
        """Test getting mutation details when no mutations present."""
        chain = ["1", "2", "X", "A", "B"]

        mutations = get_mutation_details(chain)

        assert mutations == []

    def test_mutation_details_all_mutated(self):
        """Test getting mutation details when all repeats mutated."""
        chain = ["1m", "2m", "Xm"]

        mutations = get_mutation_details(chain)

        assert len(mutations) == 3
        assert mutations[0] == {"position": 1, "repeat": "1"}
        assert mutations[1] == {"position": 2, "repeat": "2"}
        assert mutations[2] == {"position": 3, "repeat": "X"}

    def test_mutation_details_position_is_one_based(self):
        """Test that positions are 1-based."""
        chain = ["Xm"]

        mutations = get_mutation_details(chain)

        assert mutations[0]["position"] == 1  # Not 0


@pytest.mark.unit
class TestGenerateHaplotypeStats:
    """Tests for generate_haplotype_stats function."""

    def test_haplotype_stats_single_haplotype(
        self, minimal_config: dict, sample_haplotype_sequences: list
    ):
        """Test generating stats for single haplotype."""
        # Use first haplotype from fixture
        seq, _ = sample_haplotype_sequences[0]
        chain = ["1", "2", "X", "B", "6", "7", "8", "9"]
        simulation_results = [(seq, chain)]

        stats = generate_haplotype_stats(simulation_results, minimal_config)

        assert len(stats) == 1
        hap_stat = stats[0]
        assert hap_stat["repeat_count"] == 8
        assert hap_stat["vntr_length"] > 0
        assert 0 <= hap_stat["gc_content"] <= 100
        assert len(hap_stat["repeat_lengths"]) == 8
        assert hap_stat["repeat_type_counts"]["1"] == 1
        assert hap_stat["mutant_repeat_count"] == 0
        assert hap_stat["mutation_details"] == []

    def test_haplotype_stats_with_mutations(self, minimal_config: dict):
        """Test generating stats for haplotype with mutations."""
        left = minimal_config["constants"]["hg38"]["left"]
        right = minimal_config["constants"]["hg38"]["right"]
        vntr = "ATCGATCGATCG"
        seq = left + vntr + right
        chain = ["1", "2m", "X", "Am", "6", "7", "8", "9"]
        simulation_results = [(seq, chain)]

        stats = generate_haplotype_stats(simulation_results, minimal_config)

        hap_stat = stats[0]
        assert hap_stat["mutant_repeat_count"] == 2
        assert len(hap_stat["mutation_details"]) == 2

    def test_haplotype_stats_multiple_haplotypes(
        self, minimal_config: dict, sample_haplotype_sequences: list
    ):
        """Test generating stats for multiple haplotypes."""
        seq1, _ = sample_haplotype_sequences[0]
        seq2, _ = sample_haplotype_sequences[1]
        chain1 = ["1", "2", "X", "B", "6", "7", "8", "9"]
        chain2 = ["1", "2", "A", "B", "6p", "7", "8", "9"]
        simulation_results = [(seq1, chain1), (seq2, chain2)]

        stats = generate_haplotype_stats(simulation_results, minimal_config)

        assert len(stats) == 2
        assert stats[0]["repeat_count"] == 8
        assert stats[1]["repeat_count"] == 8

    def test_haplotype_stats_repeat_lengths_summary(self, minimal_config: dict):
        """Test that repeat_lengths_summary is calculated correctly."""
        left = minimal_config["constants"]["hg38"]["left"]
        right = minimal_config["constants"]["hg38"]["right"]
        vntr = "ATCGATCGATCG"
        seq = left + vntr + right
        chain = ["1", "2", "X"]
        simulation_results = [(seq, chain)]

        stats = generate_haplotype_stats(simulation_results, minimal_config)

        summary = stats[0]["repeat_lengths_summary"]
        lengths = stats[0]["repeat_lengths"]
        assert summary["min"] == min(lengths)
        assert summary["max"] == max(lengths)
        assert summary["average"] == sum(lengths) / len(lengths)


@pytest.mark.unit
class TestGenerateOverallStats:
    """Tests for generate_overall_stats function."""

    def test_overall_stats_basic(self):
        """Test generating overall stats from haplotype stats."""
        hap_stats = [
            {
                "repeat_count": 10,
                "vntr_length": 200,
                "gc_content": 45.0,
                "mutant_repeat_count": 2,
                "snp_count": 5,
            },
            {
                "repeat_count": 15,
                "vntr_length": 300,
                "gc_content": 55.0,
                "mutant_repeat_count": 3,
                "snp_count": 7,
            },
        ]

        overall = generate_overall_stats(hap_stats)

        # Check repeat_count structure
        assert overall["repeat_count"]["average"] == 12.5
        assert overall["repeat_count"]["min"] == 10
        assert overall["repeat_count"]["max"] == 15

        # Check vntr_length structure
        assert overall["vntr_length"]["average"] == 250.0
        assert overall["vntr_length"]["min"] == 200
        assert overall["vntr_length"]["max"] == 300

        # Check gc_content structure
        assert overall["gc_content"]["average"] == 50.0
        assert overall["gc_content"]["min"] == 45.0
        assert overall["gc_content"]["max"] == 55.0

    def test_overall_stats_single_haplotype(self):
        """Test overall stats with single haplotype."""
        hap_stats = [
            {
                "repeat_count": 10,
                "vntr_length": 200,
                "gc_content": 42.0,
                "mutant_repeat_count": 2,
                "snp_count": 5,
            }
        ]

        overall = generate_overall_stats(hap_stats)

        assert overall["repeat_count"]["average"] == 10.0
        assert overall["repeat_count"]["min"] == 10
        assert overall["repeat_count"]["max"] == 10


@pytest.mark.integration
class TestGenerateSimulationStatistics:
    """Integration tests for generate_simulation_statistics function."""

    def test_full_statistics_generation(
        self, minimal_config: dict, sample_haplotype_sequences: list
    ):
        """Test full statistics report generation."""
        seq1, _ = sample_haplotype_sequences[0]
        seq2, _ = sample_haplotype_sequences[1]
        chain1 = ["1", "2", "X", "B", "6", "7", "8", "9"]
        chain2 = ["1", "2", "A", "B", "6p", "7", "8", "9"]
        simulation_results = [(seq1, chain1), (seq2, chain2)]

        start_time = 100.0
        end_time = 110.5
        mutation_info = {"name": "dupC", "targets": [(1, 25)]}

        report = generate_simulation_statistics(
            start_time=start_time,
            end_time=end_time,
            simulation_results=simulation_results,
            config=minimal_config,
            mutation_info=mutation_info,
        )

        assert "runtime_seconds" in report
        assert report["runtime_seconds"] == 10.5
        assert "mutation_info" in report
        assert report["mutation_info"]["name"] == "dupC"
        assert report["mutation_info"]["targets"] == [(1, 25)]
        assert "haplotype_statistics" in report
        assert len(report["haplotype_statistics"]) == 2
        assert "overall_statistics" in report
        assert "reference_assembly" in report

    def test_statistics_generation_no_mutation(
        self, minimal_config: dict, sample_haplotype_sequences: list
    ):
        """Test statistics generation with no mutation."""
        seq1, _ = sample_haplotype_sequences[0]
        chain1 = ["1", "2", "X", "B", "6", "7", "8", "9"]
        simulation_results = [(seq1, chain1)]

        start_time = 100.0
        end_time = 105.0

        report = generate_simulation_statistics(
            start_time=start_time,
            end_time=end_time,
            simulation_results=simulation_results,
            config=minimal_config,
            mutation_info=None,
        )

        assert report["mutation_info"] == {}


@pytest.mark.unit
class TestWriteStatisticsReport:
    """Tests for write_statistics_report function."""

    def test_write_statistics_report(self, tmp_path: Path):
        """Test writing statistics report to file."""
        output_file = tmp_path / "stats.json"
        report = {
            "runtime_seconds": 10.5,
            "haplotypes": [{"repeat_count": 10}],
            "overall": {"total_repeat_count": 10},
        }

        write_statistics_report(report, str(output_file))

        assert output_file.exists()

        # Verify JSON is valid and content matches
        with output_file.open() as f:
            loaded = json.load(f)
        assert loaded == report

    def test_write_statistics_report_requires_existing_directory(self, tmp_path: Path):
        """Test that writing report requires parent directories to exist."""
        output_file = tmp_path / "subdir" / "stats.json"
        report = {"runtime_seconds": 10.5}

        # Should fail because parent directory doesn't exist
        with pytest.raises(FileNotFoundError):
            write_statistics_report(report, str(output_file))


@pytest.mark.bioinformatics
class TestStatisticsBioinformatics:
    """Bioinformatics-specific tests for statistics."""

    def test_gc_content_realistic_sequences(self):
        """Test GC content calculation on realistic sequences."""
        # Human genome average GC content is ~40-42%
        # MUC1 VNTR GC content varies by repeat type
        gc_rich_seq = "GCGCGCGC" * 10  # High GC
        at_rich_seq = "ATATATAT" * 10  # Low GC

        gc_high = compute_gc_content(gc_rich_seq)
        gc_low = compute_gc_content(at_rich_seq)

        assert gc_high > 90.0
        assert gc_low < 10.0

    def test_repeat_statistics_realistic_vntr(self, minimal_config: dict):
        """Test repeat statistics for realistic VNTR."""
        # Typical MUC1 VNTR: 20-120 repeats
        chain = ["1", "2", "X", "A", "B", "C", "6", "7", "8", "9"] * 5  # 50 repeats
        left = minimal_config["constants"]["hg38"]["left"]
        right = minimal_config["constants"]["hg38"]["right"]
        vntr = "ATCGATCG" * 50
        seq = left + vntr + right

        simulation_results = [(seq, chain)]
        stats = generate_haplotype_stats(simulation_results, minimal_config)

        assert stats[0]["repeat_count"] == 50
        assert stats[0]["vntr_length"] > 0
        assert len(stats[0]["repeat_type_counts"]) > 1


@pytest.mark.unit
class TestProvenanceIntegration:
    """Tests for provenance metadata integration in simulation statistics."""

    def test_provenance_in_report(self, minimal_config: dict):
        """Test provenance section is included when provided."""
        results = [("ACGTACGT", ["X", "A"])]
        start = 1730649330.0
        end = 1730649331.0

        provenance_info = {
            "software_version": "0.27.0",
            "config_fingerprint": "sha256:abc123...",
            "seed": 42,
            "start_time": "2025-11-03T10:00:00.000000+00:00",
            "end_time": "2025-11-03T10:00:01.000000+00:00",
            "duration_seconds": 1.0,
            "command_line": "muconeup --config config.json simulate",
        }

        report = generate_simulation_statistics(
            start_time=start,
            end_time=end,
            simulation_results=results,
            config=minimal_config,
            provenance_info=provenance_info,
        )

        assert "provenance" in report
        assert report["provenance"] == provenance_info

    def test_provenance_optional_backward_compatible(self, minimal_config: dict):
        """Test backward compatibility when provenance not provided."""
        results = [("ACGTACGT", ["X", "A"])]
        start = 1730000000.0
        end = start + 1.0

        # OLD CALL: no provenance_info parameter
        report = generate_simulation_statistics(
            start_time=start,
            end_time=end,
            simulation_results=results,
            config=minimal_config,
        )

        assert "provenance" in report
        assert report["provenance"] == {}  # Empty dict (graceful degradation)

    def test_runtime_computed_from_provenance(self, minimal_config: dict):
        """Test runtime_seconds uses provenance.duration_seconds if available."""
        results = [("ACGTACGT", ["X", "A"])]
        start = 1000.0
        end = 1001.5

        provenance_info = {
            "duration_seconds": 1.5,
            "software_version": "0.27.0",
        }

        report = generate_simulation_statistics(
            start_time=start,
            end_time=end,
            simulation_results=results,
            config=minimal_config,
            provenance_info=provenance_info,
        )

        # runtime_seconds should match provenance (single source of truth)
        assert report["runtime_seconds"] == 1.5
        assert report["runtime_seconds"] == provenance_info["duration_seconds"]

    def test_runtime_falls_back_to_calculation(self, minimal_config: dict):
        """Test runtime_seconds falls back to calculation if provenance missing."""
        results = [("ACGTACGT", ["X", "A"])]
        start = 1000.0
        end = 1002.5

        # No provenance provided
        report = generate_simulation_statistics(
            start_time=start,
            end_time=end,
            simulation_results=results,
            config=minimal_config,
        )

        # Should fall back to direct calculation
        assert report["runtime_seconds"] == 2.5

    def test_provenance_with_all_fields(self, minimal_config: dict):
        """Test report with complete provenance metadata."""
        results = [("ACGTACGT", ["X", "A"])]
        start = 1730649330.0
        end = 1730649331.5

        provenance_info = {
            "software_version": "0.28.0",
            "config_fingerprint": "sha256:8f3e9a7b2c1d4e5f6a7b8c9d0e1f2a3b4c5d6e7f8a9b0c1d2e3f4a5b6c7d8e9f",
            "seed": 42,
            "start_time": "2024-11-03T10:15:30.000000+00:00",
            "end_time": "2024-11-03T10:15:31.500000+00:00",
            "duration_seconds": 1.5,
            "command_line": "muconeup --config config.json simulate --seed 42",
        }

        report = generate_simulation_statistics(
            start_time=start,
            end_time=end,
            simulation_results=results,
            config=minimal_config,
            provenance_info=provenance_info,
        )

        assert report["provenance"]["software_version"] == "0.28.0"
        assert report["provenance"]["config_fingerprint"].startswith("sha256:")
        assert report["provenance"]["seed"] == 42
        assert "T" in report["provenance"]["start_time"]  # ISO 8601
        assert "T" in report["provenance"]["end_time"]  # ISO 8601

    def test_provenance_with_error_sentinels(self, minimal_config: dict):
        """Test report handles provenance with error sentinels."""
        results = [("ACGTACGT", ["X", "A"])]
        start = 1000.0
        end = 1001.0

        provenance_info = {
            "software_version": "0.28.0",
            "config_fingerprint": "error:fingerprint_failed:TypeError",
            "seed": None,
            "start_time": "error:timestamp_failed",
            "end_time": "error:timestamp_failed",
            "duration_seconds": 1.0,
            "command_line": "error:cmdline_failed",
        }

        report = generate_simulation_statistics(
            start_time=start,
            end_time=end,
            simulation_results=results,
            config=minimal_config,
            provenance_info=provenance_info,
        )

        # Should not crash, includes error sentinels
        assert report["provenance"]["config_fingerprint"].startswith("error:")
        assert report["provenance"]["seed"] is None
        assert report["runtime_seconds"] == 1.0  # Duration still works
