"""Tests for muc_one_up/validation.py - general validation functions."""

import pytest

from muc_one_up.exceptions import FileOperationError, ValidationError
from muc_one_up.validation import (
    validate_directory_exists,
    validate_file_exists,
    validate_haplotype_index,
    validate_mutation_exists,
    validate_mutation_targets,
    validate_non_negative_integer,
    validate_positive_integer,
    validate_probability,
    validate_reference_assembly,
    validate_repeat_index,
    validate_repeat_symbol,
)


class TestValidateHaplotypeIndex:
    """Test haplotype index validation."""

    def test_valid_index(self):
        """Valid haplotype indices should not raise errors."""
        validate_haplotype_index(0, 2)
        validate_haplotype_index(1, 2)

    def test_negative_index(self):
        """Negative indices should raise ValidationError."""
        with pytest.raises(ValidationError, match="out of range"):
            validate_haplotype_index(-1, 2)

    def test_index_too_large(self):
        """Indices >= num_haplotypes should raise ValidationError."""
        with pytest.raises(ValidationError, match="out of range"):
            validate_haplotype_index(2, 2)

    def test_boundary_conditions(self):
        """Test boundary conditions."""
        validate_haplotype_index(0, 1)  # Should pass
        with pytest.raises(ValidationError):
            validate_haplotype_index(1, 1)  # Should fail


class TestValidateRepeatIndex:
    """Test repeat index validation."""

    def test_valid_index(self):
        """Valid repeat indices should not raise errors."""
        validate_repeat_index(0, 10)
        validate_repeat_index(9, 10)

    def test_with_context(self):
        """Error message should include context."""
        with pytest.raises(ValidationError, match="mutation application"):
            validate_repeat_index(10, 10, context="mutation application")

    def test_negative_index(self):
        """Negative indices should raise ValidationError."""
        with pytest.raises(ValidationError):
            validate_repeat_index(-1, 10)


class TestValidateFileExists:
    """Test file existence validation."""

    def test_existing_file(self, tmp_path):
        """Should not raise error for existing file."""
        test_file = tmp_path / "test.txt"
        test_file.write_text("test")
        validate_file_exists(test_file)

    def test_missing_file(self):
        """Should raise FileOperationError for missing file."""
        with pytest.raises(FileOperationError, match="not found"):
            validate_file_exists("/nonexistent/file.txt")

    def test_custom_description(self):
        """Error message should include custom description."""
        with pytest.raises(FileOperationError, match="Config file not found"):
            validate_file_exists("/nonexistent/config.json", description="Config file")


class TestValidateDirectoryExists:
    """Test directory existence validation."""

    def test_existing_directory(self, tmp_path):
        """Should not raise error for existing directory."""
        validate_directory_exists(tmp_path)

    def test_missing_directory(self):
        """Should raise FileOperationError for missing directory."""
        with pytest.raises(FileOperationError, match="not found"):
            validate_directory_exists("/nonexistent/dir")

    def test_file_instead_of_directory(self, tmp_path):
        """Should raise error if path is file not directory."""
        test_file = tmp_path / "test.txt"
        test_file.write_text("test")
        with pytest.raises(FileOperationError, match="not a directory"):
            validate_directory_exists(test_file)


class TestValidateMutationExists:
    """Test mutation existence validation."""

    def test_existing_mutation(self):
        """Should not raise error for existing mutation."""
        config = {"mutations": {"dupC": {"changes": []}}}
        validate_mutation_exists("dupC", config)

    def test_missing_mutation(self):
        """Should raise ValidationError with available mutations."""
        config = {"mutations": {"dupC": {}, "delA": {}}}
        with pytest.raises(ValidationError, match=r"dupC.*delA"):
            validate_mutation_exists("missing", config)

    def test_no_mutations_section(self):
        """Should raise ValidationError if no mutations section."""
        config = {}
        with pytest.raises(ValidationError):
            validate_mutation_exists("dupC", config)


class TestValidateMutationTargets:
    """Test mutation targets validation."""

    def test_valid_targets(self):
        """Valid targets should not raise errors."""
        targets = [(1, 5), (2, 10)]
        chains = [["1", "2", "7", "8", "9"] * 2, ["1", "2", "7", "8", "9"] * 3]
        validate_mutation_targets(targets, 2, chains)

    def test_invalid_haplotype_index(self):
        """Should raise ValidationError for invalid haplotype."""
        targets = [(3, 5)]
        chains = [["1", "2"], ["1", "2"]]
        with pytest.raises(ValidationError, match=r"haplotype index.*out of range"):
            validate_mutation_targets(targets, 2, chains)

    def test_invalid_repeat_index(self):
        """Should raise ValidationError for invalid repeat position."""
        targets = [(1, 10)]
        chains = [["1", "2", "7"]]
        with pytest.raises(ValidationError, match=r"repeat index.*out of range"):
            validate_mutation_targets(targets, 1, chains)


class TestValidateRepeatSymbol:
    """Test repeat symbol validation."""

    def test_valid_symbol(self):
        """Valid symbols should not raise errors."""
        valid = {"1", "2", "7", "8", "9", "X"}
        validate_repeat_symbol("1", valid)
        validate_repeat_symbol("X", valid)

    def test_symbol_with_mutation_marker(self):
        """Symbols with 'm' suffix should be valid."""
        valid = {"1", "2", "7"}
        validate_repeat_symbol("1m", valid)

    def test_invalid_symbol(self):
        """Invalid symbols should raise ValidationError."""
        valid = {"1", "2", "7"}
        with pytest.raises(ValidationError, match="Invalid repeat symbol"):
            validate_repeat_symbol("Z", valid)


class TestValidatePositiveInteger:
    """Test positive integer validation."""

    def test_valid_positive(self):
        """Positive integers should not raise errors."""
        validate_positive_integer(1, "count")
        validate_positive_integer(100, "count")

    def test_zero(self):
        """Zero should raise ValidationError."""
        with pytest.raises(ValidationError, match="positive"):
            validate_positive_integer(0, "count")

    def test_negative(self):
        """Negative values should raise ValidationError."""
        with pytest.raises(ValidationError, match="positive"):
            validate_positive_integer(-1, "count")


class TestValidateNonNegativeInteger:
    """Test non-negative integer validation."""

    def test_valid_values(self):
        """Non-negative values should not raise errors."""
        validate_non_negative_integer(0, "count")
        validate_non_negative_integer(100, "count")

    def test_negative(self):
        """Negative values should raise ValidationError."""
        with pytest.raises(ValidationError, match="non-negative"):
            validate_non_negative_integer(-1, "count")


class TestValidateProbability:
    """Test probability validation."""

    def test_valid_probabilities(self):
        """Values in [0.0, 1.0] should not raise errors."""
        validate_probability(0.0, "prob")
        validate_probability(0.5, "prob")
        validate_probability(1.0, "prob")

    def test_below_range(self):
        """Values < 0.0 should raise ValidationError."""
        with pytest.raises(ValidationError, match=r"between 0\.0 and 1\.0"):
            validate_probability(-0.1, "prob")

    def test_above_range(self):
        """Values > 1.0 should raise ValidationError."""
        with pytest.raises(ValidationError, match=r"between 0\.0 and 1\.0"):
            validate_probability(1.1, "prob")


class TestValidateReferenceAssembly:
    """Test reference assembly validation."""

    def test_valid_assemblies(self):
        """hg19 and hg38 should be valid."""
        validate_reference_assembly("hg19")
        validate_reference_assembly("hg38")

    def test_invalid_assembly(self):
        """Other assemblies should raise ValidationError."""
        with pytest.raises(ValidationError, match="Invalid reference assembly"):
            validate_reference_assembly("hg37")

        with pytest.raises(ValidationError, match="Invalid reference assembly"):
            validate_reference_assembly("grch38")
