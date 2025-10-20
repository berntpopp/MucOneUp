"""
Integration tests for SNaPshot validator with real simulated data.

Tests complete workflow on real dupC mutant and normal samples.
"""

from pathlib import Path

import pytest
from Bio import SeqIO

from muc_one_up.analysis.snapshot_validator import SnapshotValidator
from muc_one_up.config import load_config


@pytest.fixture
def config():
    """Load actual config.json with snapshot_validation section."""
    config_path = Path("config.json")
    if not config_path.exists():
        pytest.skip("config.json not found")
    return load_config(str(config_path))


@pytest.fixture
def test_data_dir():
    """Test data directory."""
    return Path("output/dupC")


@pytest.fixture
def mutant_sample_path(test_data_dir):
    """Path to mutant sample FASTA."""
    # Try different possible naming patterns
    possible_paths = [
        test_data_dir / "dupC.001.mut.simulated.fa",
        test_data_dir / "dupC.mut.simulated.fa",
        test_data_dir / "dupC.001.simulated.fa",
    ]

    for path in possible_paths:
        if path.exists():
            return path

    pytest.skip(
        f"Mutant sample not found in {test_data_dir}. Generate with: muconeup --config config.json simulate --out-base dupC --mutation-name dupC --fixed-lengths 27 --out-dir output/dupC/"
    )


@pytest.fixture
def normal_sample_path(test_data_dir):
    """Path to normal sample FASTA."""
    possible_paths = [
        test_data_dir / "dupC.001.normal.simulated.fa",
        test_data_dir / "dupC.normal.simulated.fa",
    ]

    for path in possible_paths:
        if path.exists():
            return path

    pytest.skip(
        f"Normal sample not found in {test_data_dir}. Generate with: muconeup --config config.json simulate --out-base dupC --mutation-name normal,dupC --fixed-lengths 27 --out-dir output/dupC/"
    )


class TestIntegrationWithRealData:
    """Integration tests using real simulated dupC data."""

    def test_mutant_sample_complete_workflow(self, config, mutant_sample_path):
        """Test complete workflow on dupC mutant sample (positive control)."""
        # Initialize validator
        validator = SnapshotValidator(config, "dupC")

        # Load mutant sample (haplotype 1 or combined)
        records = list(SeqIO.parse(str(mutant_sample_path), "fasta"))
        assert len(records) >= 1, "No sequences found in mutant sample"

        # Test on first haplotype
        template = str(records[0].seq)

        # Run complete workflow
        result = validator.validate_complete_workflow(template)

        # Assertions for mutant sample
        assert result["pcr_results"]["success"] is True, "PCR should succeed"
        assert result["pcr_results"]["product_count"] > 0, "Should have PCR products"

        # Should have survivors (8C amplicons lack MwoI site)
        assert result["digest_results"]["survivor_count"] >= 1, (
            "Should have digest survivors (8C amplicons)"
        )

        # Mutation should be detected
        assert result["mutation_detected"] is True, "dupC mutation should be detected"

        # Should have SNaPshot results
        assert len(result["snapshot_results"]) >= 1, "Should have SNaPshot results"

        # Check for Black (dTAMRA) fluorescence
        black_peaks = [
            r for r in result["snapshot_results"] if r["extension"].get("fluorophore") == "Black"
        ]
        assert len(black_peaks) >= 1, "Should have Black (dTAMRA) fluorescence peaks"

        # Check summary
        assert "DETECTED" in result["summary"], "Summary should indicate detection"

        print("\n✓ Mutant Sample Results:")
        print(f"  PCR products: {result['pcr_results']['product_count']}")
        print(f"  Digest survivors: {result['digest_results']['survivor_count']}")
        print(f"  Mutation detected: {result['mutation_detected']}")
        print(f"  Summary: {result['summary']}")

    def test_normal_sample_complete_workflow(self, config, normal_sample_path):
        """Test complete workflow on normal sample (negative control)."""
        # Initialize validator
        validator = SnapshotValidator(config, "dupC")

        # Load normal sample
        records = list(SeqIO.parse(str(normal_sample_path), "fasta"))
        assert len(records) >= 1, "No sequences found in normal sample"

        # Test on first haplotype
        template = str(records[0].seq)

        # Run complete workflow
        result = validator.validate_complete_workflow(template)

        # Assertions for normal sample
        assert result["pcr_results"]["success"] is True, "PCR should succeed"
        assert result["pcr_results"]["product_count"] > 0, "Should have PCR products"

        # Normal sample: all 7C amplicons should be digested
        # May have 0 survivors (all digested) or very few
        # Mutation should NOT be detected
        assert result["mutation_detected"] is False, (
            "dupC mutation should NOT be detected in normal sample"
        )

        # Check summary indicates no mutation
        assert "No dupC" in result["summary"] or "normal" in result["summary"].lower(), (
            "Summary should indicate no mutation"
        )

        print("\n✓ Normal Sample Results:")
        print(f"  PCR products: {result['pcr_results']['product_count']}")
        print(f"  Digest survivors: {result['digest_results']['survivor_count']}")
        print(f"  Mutation detected: {result['mutation_detected']}")
        print(f"  Summary: {result['summary']}")

    def test_both_haplotypes_mutant(self, config, mutant_sample_path):
        """Test workflow on both haplotypes of mutant sample."""
        validator = SnapshotValidator(config, "dupC")

        records = list(SeqIO.parse(str(mutant_sample_path), "fasta"))

        results = []
        for i, record in enumerate(records, 1):
            template = str(record.seq)
            result = validator.validate_complete_workflow(template)
            results.append((i, result))

            print(f"\n✓ Haplotype {i}:")
            print(f"  PCR products: {result['pcr_results']['product_count']}")
            print(f"  Digest survivors: {result['digest_results']['survivor_count']}")
            print(f"  Mutation detected: {result['mutation_detected']}")

        # At least one haplotype should show mutation
        mutations_detected = sum(1 for _, r in results if r["mutation_detected"])
        assert mutations_detected >= 1, "At least one haplotype should show mutation"

    def test_amplicon_validation_mutant(self, config, mutant_sample_path):
        """Test amplicon validation on mutant sample."""
        validator = SnapshotValidator(config, "dupC")

        # Load sample
        records = list(SeqIO.parse(str(mutant_sample_path), "fasta"))
        template = str(records[0].seq)

        # Run PCR to get amplicons
        pcr_result = validator.pcr.amplify(
            template=template,
            mutation_patterns={
                "mutant": validator.mutation_patterns["mutant_pattern"],
                "normal": validator.mutation_patterns["normal_pattern"],
            },
        )

        # Test validate_amplicon on each product
        mutant_amplicons = 0
        for product in pcr_result["products"][:10]:  # Test first 10
            validation = validator.validate_amplicon(product["sequence"])

            if validation["mutation_present"]:
                mutant_amplicons += 1
                # Mutant amplicons should survive digest
                assert validation["survives_digest"] is True, (
                    "Mutant amplicon should survive digest"
                )
                assert validation["has_restriction_site"] is False, (
                    "Mutant amplicon should not have MwoI site"
                )
                assert validation["detectable"] is True, "Mutant amplicon should be detectable"

        print(f"\n✓ Amplicon Validation: {mutant_amplicons} mutant amplicons found")
        assert mutant_amplicons >= 1, "Should find mutant amplicons"

    def test_performance_benchmark(self, config, mutant_sample_path):
        """Benchmark performance of complete workflow."""
        import time

        validator = SnapshotValidator(config, "dupC")

        records = list(SeqIO.parse(str(mutant_sample_path), "fasta"))
        template = str(records[0].seq)

        # Benchmark
        start_time = time.time()
        result = validator.validate_complete_workflow(template)
        elapsed_time = time.time() - start_time

        print("\n✓ Performance Benchmark:")
        print(f"  Elapsed time: {elapsed_time:.3f} seconds")
        print(f"  PCR products: {result['pcr_results']['product_count']}")
        print(f"  Digest survivors: {result['digest_results']['survivor_count']}")

        # Should complete in reasonable time (< 5 seconds for single haplotype)
        assert elapsed_time < 5.0, f"Workflow should complete in < 5s, took {elapsed_time:.3f}s"

    def test_false_positive_rate(self, config, normal_sample_path):
        """Test false positive rate on normal sample."""
        validator = SnapshotValidator(config, "dupC")

        # Load normal sample
        records = list(SeqIO.parse(str(normal_sample_path), "fasta"))

        false_positives = 0
        total_haplotypes = 0

        for record in records:
            template = str(record.seq)
            result = validator.validate_complete_workflow(template)

            total_haplotypes += 1
            if result["mutation_detected"]:
                false_positives += 1

        false_positive_rate = false_positives / total_haplotypes if total_haplotypes > 0 else 0

        print(f"\n✓ False Positive Rate: {false_positive_rate:.1%}")
        print(f"  False positives: {false_positives}/{total_haplotypes}")

        # Should have 0% false positive rate
        assert false_positive_rate == 0.0, "Should have 0% false positive rate on normal samples"


class TestConfigBasedValidation:
    """Test that implementation is truly config-based (no hardcoded values)."""

    def test_config_mutation_patterns_used(self, config):
        """Verify mutation patterns come from config."""
        validator = SnapshotValidator(config, "dupC")

        # Check patterns come from config
        assert validator.mutation_patterns["mutant_pattern"] == "GCCCCCCCCAGC"
        assert validator.mutation_patterns["normal_pattern"] == "GCCCCCCCAGC"
        assert validator.mutation_patterns["expected_mutant_fluorescence"] == "Black"

    def test_config_primers_used(self, config):
        """Verify primers come from config."""
        validator = SnapshotValidator(config, "dupC")

        # Check PCR primers from config
        assert validator.pcr.forward_primer == "GGCCGGCCCCGGGCTCCACC"
        assert validator.pcr.reverse_primer == "TGTCACCTCGGCCCCGGA"  # As-is from config

        # Check SNaPshot primers from config
        assert "primer_7c" in validator.extension.primers
        assert validator.extension.primers["primer_7c"] == "CGGGCTCCACCGCCCCCCC"

    def test_config_enzyme_used(self, config):
        """Verify enzyme comes from config."""
        validator = SnapshotValidator(config, "dupC")

        assert validator.digest.enzyme_name == "MwoI"
        assert validator.digest.recognition_site == "GCNNNNNNNGC"

    def test_config_size_range_used(self, config):
        """Verify size range comes from config."""
        validator = SnapshotValidator(config, "dupC")

        assert validator.pcr.size_min == 50
        assert validator.pcr.size_max == 65

    def test_config_fluorophore_map_used(self, config):
        """Verify fluorophore mapping comes from config."""
        validator = SnapshotValidator(config, "dupC")

        fluor_map = validator.extension.fluorophore_map

        assert fluor_map["C"]["color"] == "Black"
        assert fluor_map["C"]["dye"] == "dTAMRA"
        assert fluor_map["A"]["color"] == "Green"
        assert fluor_map["A"]["dye"] == "dR6G"
