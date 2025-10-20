"""
Unit tests for SNaPshot validator module.

Tests each component independently with mocked config following
SOLID principles and modular design.
"""

import pytest

from muc_one_up.analysis.snapshot_validator import (
    DigestSimulator,
    PCRSimulator,
    SnapshotExtensionSimulator,
    SnapshotValidator,
)


@pytest.fixture
def pcr_config():
    """Mock PCR configuration."""
    return {
        "pcr": {
            "forward_primer": "GGCCGGCCCCGGGCTCCACC",
            "reverse_primer": "TCCGGGGCCGAGGTGACA",  # Already RC'd for testing
            "reverse_needs_rc": False,  # False since we're providing it already RC'd
            "max_products": 100,
            "size_range": {"min": 50, "max": 65},
        }
    }


@pytest.fixture
def digest_config():
    """Mock digest configuration."""
    return {"digest": {"enzyme": "MwoI", "recognition_site": "GCNNNNNNNGC"}}


@pytest.fixture
def snapshot_config():
    """Mock SNaPshot configuration."""
    return {
        "snapshot": {
            "primers": {
                "primer_7c": "CGGGCTCCACCGCCCCCCC",
                "primer_repeat_r": "TGTCACCTCGGCCCCGGA",
            },
            "fluorophore_map": {
                "A": {"color": "Green", "dye": "dR6G"},
                "C": {"color": "Black", "dye": "dTAMRA"},
                "G": {"color": "Blue", "dye": "dR110"},
                "T": {"color": "Red", "dye": "dROX"},
            },
        }
    }


@pytest.fixture
def complete_config(pcr_config, digest_config, snapshot_config):
    """Complete configuration for SnapshotValidator."""
    return {
        "snapshot_validation": {
            "dupC": {
                "description": "8C mutation validation",
                **pcr_config,
                **digest_config,
                **snapshot_config,
                "validation": {
                    "mutant_pattern": "GCCCCCCCCAGC",
                    "normal_pattern": "GCCCCCCCAGC",
                    "expected_mutant_fluorescence": "Black",
                },
            }
        }
    }


# PCRSimulator Tests
class TestPCRSimulator:
    """Test PCRSimulator class."""

    def test_initialization(self, pcr_config):
        """Test PCR simulator initializes correctly from config."""
        pcr_sim = PCRSimulator(pcr_config)

        assert pcr_sim.forward_primer == "GGCCGGCCCCGGGCTCCACC"
        assert pcr_sim.max_products == 100
        assert pcr_sim.size_min == 50
        assert pcr_sim.size_max == 65
        # Reverse primer should stay as-is (already provided in correct form)
        assert pcr_sim.reverse_primer == "TCCGGGGCCGAGGTGACA"

    def test_amplify_no_binding(self, pcr_config):
        """Test PCR when primers don't bind."""
        pcr_sim = PCRSimulator(pcr_config)

        # Template without primer binding sites
        template = "ATATATATATATATATATAT"

        result = pcr_sim.amplify(
            template=template,
            mutation_patterns={"mutant": "GCCCCCCCCAGC", "normal": "GCCCCCCCAGC"},
        )

        assert result["success"] is False
        assert result["product_count"] == 0
        assert len(result["products"]) == 0

    def test_amplify_with_products(self, pcr_config):
        """Test PCR amplification produces products."""
        pcr_sim = PCRSimulator(pcr_config)

        # Create template with primer binding sites and flanked region
        # [Forward primer] - [Flanked] - [Reverse primer RC]
        template = (
            "GGCCGGCCCCGGGCTCCACC"  # Forward primer
            + "GCCCCCCCAGCCCACGG"  # Flanked region (7C)
            + "TCCGGGGCCGAGGTGACA"  # Reverse primer (RC)
        )

        result = pcr_sim.amplify(
            template=template,
            mutation_patterns={"mutant": "GCCCCCCCCAGC", "normal": "GCCCCCCCAGC"},
        )

        assert result["success"] is True
        assert result["product_count"] >= 1
        assert len(result["products"]) >= 1

        # Check first product
        product = result["products"][0]
        assert "sequence" in product
        assert "length" in product
        assert "flanked_region" in product
        assert "mutation_type" in product

    def test_amplify_mutation_categorization(self, pcr_config):
        """Test products are categorized by mutation type."""
        pcr_sim = PCRSimulator(pcr_config)

        # Template with 8C mutant pattern
        template_8c = (
            "GGCCGGCCCCGGGCTCCACC"
            + "GCCCCCCCCAGCCCACGG"  # 8C mutant
            + "TCCGGGGCCGAGGTGACA"
        )

        result = pcr_sim.amplify(
            template=template_8c,
            mutation_patterns={"mutant": "GCCCCCCCCAGC", "normal": "GCCCCCCCAGC"},
        )

        assert result["success"] is True
        # Should have at least one product categorized as mutant
        mutant_products = [p for p in result["products"] if p["mutation_type"] == "mutant"]
        assert len(mutant_products) >= 1

    def test_amplify_size_filtering(self, pcr_config):
        """Test size filtering removes too small/large products."""
        pcr_sim = PCRSimulator(pcr_config)

        # Template that would produce product too small (< 50bp)
        template_small = "GGCCGGCCCCGGGCTCCACC" + "ATCG" + "TCCGGGGCCGAGGTGACA"

        result = pcr_sim.amplify(
            template=template_small,
            mutation_patterns={"mutant": "GCCCCCCCCAGC", "normal": "GCCCCCCCAGC"},
        )

        # Should have no products (filtered by size)
        assert result["product_count"] == 0


# DigestSimulator Tests
class TestDigestSimulator:
    """Test DigestSimulator class."""

    def test_initialization(self, digest_config):
        """Test digest simulator initializes correctly."""
        digest_sim = DigestSimulator(digest_config)

        assert digest_sim.enzyme_name == "MwoI"
        assert digest_sim.recognition_site == "GCNNNNNNNGC"
        assert digest_sim.enzyme is not None

    def test_initialization_invalid_enzyme(self):
        """Test initialization with invalid enzyme name."""
        config = {"digest": {"enzyme": "InvalidEnzyme123", "recognition_site": "XXX"}}

        with pytest.raises(ValueError, match="Unknown restriction enzyme"):
            DigestSimulator(config)

    def test_digest_7c_has_site(self, digest_config):
        """Test that 7C amplicon has MwoI site and gets digested."""
        digest_sim = DigestSimulator(digest_config)

        # 7C amplicon: GC[CCCCCCC]AGC - exactly 7 bases = HAS SITE
        amplicon_7c = {
            "sequence": "GGCCGGCCCCGGGCTCCACCGCCCCCCCAGCCCACGGTCCGGGGCCGAGGTGACA",
            "length": 55,
            "mutation_type": "normal",
        }

        result = digest_sim.digest([amplicon_7c])

        assert result["total_products"] == 1
        assert result["survivor_count"] == 0  # 7C gets digested
        assert result["digested_count"] == 1
        assert len(result["digested"]) == 1
        assert result["digested"][0]["has_restriction_site"] is True
        assert result["digested"][0]["site_count"] >= 1

    def test_digest_8c_no_site(self, digest_config):
        """Test that 8C amplicon lacks MwoI site and survives."""
        digest_sim = DigestSimulator(digest_config)

        # 8C amplicon: GC[CCCCCCCC]AGC - 8 bases = NO SITE
        amplicon_8c = {
            "sequence": "GGCCGGCCCCGGGCTCCACCGCCCCCCCCAGCCCACGGTCCGGGGCCGAGGTGACA",
            "length": 56,
            "mutation_type": "mutant",
        }

        result = digest_sim.digest([amplicon_8c])

        assert result["total_products"] == 1
        assert result["survivor_count"] == 1  # 8C survives
        assert result["digested_count"] == 0
        assert len(result["survivors"]) == 1
        assert result["survivors"][0]["has_restriction_site"] is False
        assert result["survivors"][0]["site_count"] == 0

    def test_digest_mixed_amplicons(self, digest_config):
        """Test digest with mix of 7C and 8C amplicons."""
        digest_sim = DigestSimulator(digest_config)

        amplicons = [
            {
                "sequence": "GGCCGGCCCCGGGCTCCACCGCCCCCCCAGCCCACGGTCCGGGGCCGAGGTGACA",
                "mutation_type": "normal",
            },  # 7C
            {
                "sequence": "GGCCGGCCCCGGGCTCCACCGCCCCCCCCAGCCCACGGTCCGGGGCCGAGGTGACA",
                "mutation_type": "mutant",
            },  # 8C
            {
                "sequence": "GGCCGGCCCCGGGCTCCACCGCCCCCCCAGCCCACGGTCCGGGGCCGAGGTGACA",
                "mutation_type": "normal",
            },  # 7C
        ]

        result = digest_sim.digest(amplicons)

        assert result["total_products"] == 3
        assert result["survivor_count"] == 1  # Only 8C survives
        assert result["digested_count"] == 2  # Both 7C digested


# SnapshotExtensionSimulator Tests
class TestSnapshotExtensionSimulator:
    """Test SnapshotExtensionSimulator class."""

    def test_initialization(self, snapshot_config):
        """Test extension simulator initializes correctly."""
        ext_sim = SnapshotExtensionSimulator(snapshot_config)

        assert "primer_7c" in ext_sim.primers
        assert ext_sim.primers["primer_7c"] == "CGGGCTCCACCGCCCCCCC"
        assert "C" in ext_sim.fluorophore_map
        assert ext_sim.fluorophore_map["C"]["color"] == "Black"
        assert ext_sim.fluorophore_map["C"]["dye"] == "dTAMRA"

    def test_extend_no_binding(self, snapshot_config):
        """Test extension when primer doesn't bind."""
        ext_sim = SnapshotExtensionSimulator(snapshot_config)

        # Template without primer binding site
        template = "ATATATATATATATAT"

        result = ext_sim.extend(template=template, primer_name="primer_7c")

        assert result["binds"] is False
        assert result["position"] == -1

    def test_extend_8c_black_peak(self, snapshot_config):
        """Test that 8C mutation produces Black (dTAMRA) fluorescence."""
        ext_sim = SnapshotExtensionSimulator(snapshot_config)

        # Template: Primer + 8C (next base is C)
        # CGGGCTCCACCGCCCCCCC (primer) + C (8C extension) + rest
        template = "CGGGCTCCACCGCCCCCCCCAGCCCACGG"

        result = ext_sim.extend(template=template, primer_name="primer_7c")

        assert result["binds"] is True
        assert result["next_base"] == "C"
        assert result["fluorophore"] == "Black"
        assert result["fluorophore_dye"] == "dTAMRA"
        assert result["ddNTP"] == "ddCTP"
        assert result["peak_size"] == 20  # 19bp primer + 1
        assert "8C mutation detected" in result["interpretation"]

    def test_extend_7c_green_peak(self, snapshot_config):
        """Test that 7C produces Green (dR6G) fluorescence."""
        ext_sim = SnapshotExtensionSimulator(snapshot_config)

        # Template: Primer + 7C (next base is A)
        # CGGGCTCCACCGCCCCCCC (primer) + A (7C continues to AGC) + rest
        template = "CGGGCTCCACCGCCCCCCCAGCCCACGG"

        result = ext_sim.extend(template=template, primer_name="primer_7c")

        assert result["binds"] is True
        assert result["next_base"] == "A"
        assert result["fluorophore"] == "Green"
        assert result["fluorophore_dye"] == "dR6G"
        assert result["ddNTP"] == "ddATP"

    def test_extend_invalid_primer_name(self, snapshot_config):
        """Test extension with invalid primer name raises error."""
        ext_sim = SnapshotExtensionSimulator(snapshot_config)

        with pytest.raises(ValueError, match="Unknown primer"):
            ext_sim.extend(template="ATCGATCG", primer_name="invalid_primer")


# SnapshotValidator Tests
class TestSnapshotValidator:
    """Test SnapshotValidator orchestrator class."""

    def test_initialization(self, complete_config):
        """Test validator initializes correctly."""
        validator = SnapshotValidator(complete_config, "dupC")

        assert validator.mutation_name == "dupC"
        assert validator.pcr is not None
        assert validator.digest is not None
        assert validator.extension is not None
        assert validator.mutation_patterns["mutant_pattern"] == "GCCCCCCCCAGC"

    def test_initialization_missing_config(self):
        """Test initialization with missing config sections."""
        with pytest.raises(ValueError, match="missing 'snapshot_validation'"):
            SnapshotValidator({}, "dupC")

    def test_initialization_missing_mutation(self, complete_config):
        """Test initialization with non-existent mutation."""
        with pytest.raises(ValueError, match="not found in config"):
            SnapshotValidator(complete_config, "nonexistent_mutation")

    def test_validate_amplicon_8c(self, complete_config):
        """Test validation of 8C amplicon."""
        validator = SnapshotValidator(complete_config, "dupC")

        # 8C amplicon (no MwoI site, has mutant pattern)
        amplicon_8c = "GGCCGGCCCCGGGCTCCACCGCCCCCCCCAGCCCACGGTCCGGGGCCGAGGTGACA"

        result = validator.validate_amplicon(amplicon_8c)

        assert result["mutation_present"] is True
        assert result["has_restriction_site"] is False
        assert result["survives_digest"] is True
        assert result["site_count"] == 0
        assert result["detectable"] is True  # Mutant + survives = detectable

    def test_validate_amplicon_7c(self, complete_config):
        """Test validation of 7C amplicon."""
        validator = SnapshotValidator(complete_config, "dupC")

        # 7C amplicon (has MwoI site, has normal pattern)
        amplicon_7c = "GGCCGGCCCCGGGCTCCACCGCCCCCCCAGCCCACGGTCCGGGGCCGAGGTGACA"

        result = validator.validate_amplicon(amplicon_7c)

        assert result["has_normal"] is True
        assert result["mutation_present"] is False
        assert result["has_restriction_site"] is True
        assert result["survives_digest"] is False
        assert result["site_count"] >= 1
        assert result["detectable"] is False  # Normal + digested = not detectable

    def test_complete_workflow_mutant(self, complete_config):
        """Test complete workflow on mutant sample."""
        validator = SnapshotValidator(complete_config, "dupC")

        # Simple template with 8C amplicon
        template = (
            "ATCGATCG" * 50  # Padding
            + "GGCCGGCCCCGGGCTCCACC"  # Forward primer
            + "GCCCCCCCCAGCCCACGG"  # 8C flanked
            + "TCCGGGGCCGAGGTGACA"  # Reverse primer RC
            + "ATCGATCG" * 50  # Padding
        )

        result = validator.validate_complete_workflow(template)

        assert result["pcr_results"]["success"] is True
        assert result["pcr_results"]["product_count"] >= 1
        assert result["digest_results"]["survivor_count"] >= 1
        assert result["mutation_detected"] is True
        assert "DETECTED" in result["summary"]

    def test_complete_workflow_normal(self, complete_config):
        """Test complete workflow on normal sample."""
        validator = SnapshotValidator(complete_config, "dupC")

        # Template with only 7C (normal)
        template = (
            "ATCGATCG" * 50
            + "GGCCGGCCCCGGGCTCCACC"
            + "GCCCCCCCAGCCCACGG"  # 7C flanked
            + "TCCGGGGCCGAGGTGACA"
            + "ATCGATCG" * 50
        )

        result = validator.validate_complete_workflow(template)

        assert result["pcr_results"]["success"] is True
        # 7C should be digested
        assert result["digest_results"]["survivor_count"] == 0
        assert result["mutation_detected"] is False
        assert "No dupC mutation detected" in result["summary"]

    def test_complete_workflow_no_pcr_products(self, complete_config):
        """Test workflow when PCR fails."""
        validator = SnapshotValidator(complete_config, "dupC")

        # Template without primer binding sites
        template = "ATATATATATAT" * 100

        result = validator.validate_complete_workflow(template)

        assert result["pcr_results"]["success"] is False
        assert result["mutation_detected"] is False
        assert "PCR failed" in result["summary"]


# Edge Cases and Error Handling
class TestEdgeCases:
    """Test edge cases and error conditions."""

    def test_pcr_max_products_limit(self, pcr_config):
        """Test that max_products limit is enforced."""
        config = pcr_config.copy()
        config["pcr"]["max_products"] = 5  # Limit to 5

        pcr_sim = PCRSimulator(config)

        # Template with many binding sites
        template = ("GGCCGGCCCCGGGCTCCACC" + "ATCG" * 20) * 10 + "TCCGGGGCCGAGGTGACA"

        result = pcr_sim.amplify(
            template=template,
            mutation_patterns={"mutant": "GCCCCCCCCAGC", "normal": "GCCCCCCCAGC"},
        )

        # Should be limited to max_products
        assert result["product_count"] <= 5

    def test_digest_empty_amplicons(self, digest_config):
        """Test digest with empty amplicon list."""
        digest_sim = DigestSimulator(digest_config)

        result = digest_sim.digest([])

        assert result["total_products"] == 0
        assert result["survivor_count"] == 0
        assert result["digested_count"] == 0

    def test_extension_primer_at_end(self, snapshot_config):
        """Test extension when primer binds at template end."""
        ext_sim = SnapshotExtensionSimulator(snapshot_config)

        # Template ends exactly at primer
        template = "CGGGCTCCACCGCCCCCCC"  # Just the primer, no extension possible

        result = ext_sim.extend(template=template, primer_name="primer_7c")

        assert result["binds"] is True
        assert result["next_base"] is None
        assert "template end" in result["interpretation"]
