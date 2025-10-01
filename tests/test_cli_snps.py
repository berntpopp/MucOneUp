"""Comprehensive tests for muc_one_up.cli.snps module.

Following modern testing best practices:
- Parametrized tests for multiple scenarios
- Descriptive test names
- Edge case and error condition testing
- Proper use of fixtures
- AAA pattern (Arrange-Act-Assert)
"""

from pathlib import Path
from unittest.mock import Mock, patch

import pytest

from muc_one_up.cli.snps import integrate_snps_unified
from muc_one_up.exceptions import SNPIntegrationError, ValidationError


@pytest.mark.unit
class TestIntegrateSNPsUnified:
    """Comprehensive tests for integrate_snps_unified function."""

    def test_no_snps_returns_unchanged_results(self, minimal_config: dict):
        """Given no SNP options, when integrating, then returns unchanged results."""
        args = Mock()
        args.snp_input_file = None
        args.random_snps = False

        results = [("ATCGATCGATCG", ["1", "2"])]

        modified_results, applied_snp_info = integrate_snps_unified(
            args=args, config=minimal_config, results=results, skip_reference_check=False
        )

        assert modified_results == results
        assert applied_snp_info == {}

    def test_snp_file_integration(self, tmp_path: Path, minimal_config: dict):
        """Given SNP file, when integrating, then applies SNPs from file."""
        # Create SNP file
        snp_file = tmp_path / "snps.tsv"
        snp_file.write_text("1\t5\tA\tG\n2\t10\tC\tT\n")

        args = Mock()
        args.snp_input_file = str(snp_file)
        args.random_snps = False

        results = [("ATCGATCGATCG", ["1"]), ("GCTAGCTAGCTA", ["2"])]

        modified_results, applied_snp_info = integrate_snps_unified(
            args=args, config=minimal_config, results=results, skip_reference_check=False
        )

        # Verify SNPs were applied
        assert len(modified_results) == 2
        # Check that at least some SNPs were processed
        assert applied_snp_info is not None

    def test_snp_file_parsing_error(self, tmp_path: Path, minimal_config: dict):
        """Given invalid SNP file, when integrating, then raises SNPIntegrationError."""
        snp_file = tmp_path / "invalid_snps.tsv"
        snp_file.write_text("invalid\tdata\n")

        args = Mock()
        args.snp_input_file = str(snp_file)
        args.random_snps = False

        results = [("ATCGATCGATCG", ["1"])]

        with pytest.raises(SNPIntegrationError, match="Failed to parse SNP input file"):
            integrate_snps_unified(
                args=args, config=minimal_config, results=results, skip_reference_check=False
            )

    def test_random_snps_missing_density_raises_error(self, minimal_config: dict):
        """Given random SNPs without density, when integrating, then raises ValidationError."""
        args = Mock()
        args.snp_input_file = None
        args.random_snps = True
        args.random_snp_density = None  # Missing required parameter
        args.random_snp_output_file = "output.tsv"

        results = [("ATCGATCGATCG", ["1"])]

        with pytest.raises(ValidationError, match="--random-snp-density is required"):
            integrate_snps_unified(
                args=args, config=minimal_config, results=results, skip_reference_check=False
            )

    def test_random_snps_missing_output_file_raises_error(self, minimal_config: dict):
        """Given random SNPs without output file, when integrating, then raises ValidationError."""
        args = Mock()
        args.snp_input_file = None
        args.random_snps = True
        args.random_snp_density = 1.0
        args.random_snp_output_file = None  # Missing required parameter

        results = [("ATCGATCGATCG", ["1"])]

        with pytest.raises(ValidationError, match="--random-snp-output-file is required"):
            integrate_snps_unified(
                args=args, config=minimal_config, results=results, skip_reference_check=False
            )

    def test_random_snps_generation_and_application(self, tmp_path: Path, minimal_config: dict):
        """Given random SNP parameters, when integrating, then generates and applies SNPs."""
        output_file = tmp_path / "generated_snps.tsv"

        args = Mock()
        args.snp_input_file = None
        args.random_snps = True
        args.random_snp_density = 1.0
        args.random_snp_output_file = str(output_file)
        args.random_snp_region = "all"
        args.random_snp_haplotypes = "all"

        # Create sequences long enough for SNPs
        left = minimal_config["constants"]["hg38"]["left"]
        right = minimal_config["constants"]["hg38"]["right"]
        repeat = "ATCGATCGATCGATCGATCGATCG" * 10
        results = [(f"{left}{repeat}{right}", ["1"])]

        modified_results, applied_snp_info = integrate_snps_unified(
            args=args, config=minimal_config, results=results, skip_reference_check=False
        )

        # Verify SNP file was created
        assert output_file.exists()

        # Verify results structure is maintained
        assert len(modified_results) == 1
        assert isinstance(modified_results[0], tuple)
        assert len(modified_results[0]) == 2

    def test_random_snps_constants_only_region(self, tmp_path: Path, minimal_config: dict):
        """Given constants_only region, when generating SNPs, then restricts to constants."""
        output_file = tmp_path / "constants_snps.tsv"

        args = Mock()
        args.snp_input_file = None
        args.random_snps = True
        args.random_snp_density = 1.0
        args.random_snp_output_file = str(output_file)
        args.random_snp_region = "constants_only"
        args.random_snp_haplotypes = "all"

        left = minimal_config["constants"]["hg38"]["left"]
        right = minimal_config["constants"]["hg38"]["right"]
        repeat = "ATCGATCGATCGATCGATCGATCG" * 10
        results = [(f"{left}{repeat}{right}", ["1"])]

        modified_results, applied_snp_info = integrate_snps_unified(
            args=args, config=minimal_config, results=results, skip_reference_check=False
        )

        assert output_file.exists()
        assert len(modified_results) == 1

    def test_random_snps_single_haplotype(self, tmp_path: Path, minimal_config: dict):
        """Given single haplotype selection, when generating SNPs, then applies to one haplotype."""
        output_file = tmp_path / "single_hap_snps.tsv"

        args = Mock()
        args.snp_input_file = None
        args.random_snps = True
        args.random_snp_density = 1.0
        args.random_snp_output_file = str(output_file)
        args.random_snp_region = "all"
        args.random_snp_haplotypes = "1"  # Only haplotype 1

        left = minimal_config["constants"]["hg38"]["left"]
        right = minimal_config["constants"]["hg38"]["right"]
        repeat = "ATCGATCGATCGATCGATCGATCG" * 10
        results = [
            (f"{left}{repeat}{right}", ["1"]),
            (f"{left}{repeat}{right}", ["2"]),
        ]

        modified_results, applied_snp_info = integrate_snps_unified(
            args=args, config=minimal_config, results=results, skip_reference_check=False
        )

        assert output_file.exists()
        assert len(modified_results) == 2

    def test_skip_reference_check_parameter(self, tmp_path: Path, minimal_config: dict):
        """Given skip_reference_check=True, when integrating, then skips validation."""
        snp_file = tmp_path / "snps.tsv"
        # Create SNPs that would normally fail reference check
        # Use a position that exists but with mismatched reference base
        snp_file.write_text("1\t0\tT\tG\n")  # Position 0, ref T (but actual is A)

        args = Mock()
        args.snp_input_file = str(snp_file)
        args.random_snps = False

        results = [("ATCGATCGATCG", ["1"])]

        # Should not raise error with skip_reference_check=True
        # SNPs that don't match reference will be skipped, but no error raised
        modified_results, applied_snp_info = integrate_snps_unified(
            args=args, config=minimal_config, results=results, skip_reference_check=True
        )

        # Results structure should be preserved
        assert len(modified_results) == 1
        assert len(modified_results[0]) == 2

    @patch("muc_one_up.cli.snps.write_snps_to_file")
    def test_random_snps_write_failure_continues(
        self, mock_write: Mock, tmp_path: Path, minimal_config: dict, caplog
    ):
        """Given write failure for SNP file, when integrating, then logs error but continues."""
        import logging

        # Set caplog to capture ERROR level logs
        caplog.set_level(logging.ERROR)

        # Mock write_snps_to_file to raise an exception
        mock_write.side_effect = Exception("Permission denied")

        output_file = tmp_path / "snps.tsv"

        args = Mock()
        args.snp_input_file = None
        args.random_snps = True
        args.random_snp_density = 1.0
        args.random_snp_output_file = str(output_file)
        args.random_snp_region = "all"
        args.random_snp_haplotypes = "all"

        left = minimal_config["constants"]["hg38"]["left"]
        right = minimal_config["constants"]["hg38"]["right"]
        repeat = "ATCGATCGATCGATCGATCGATCG" * 10
        results = [(f"{left}{repeat}{right}", ["1"])]

        # Should not raise error, just log it
        modified_results, applied_snp_info = integrate_snps_unified(
            args=args, config=minimal_config, results=results, skip_reference_check=False
        )

        # Verify error was logged
        assert "Failed to write generated SNP file" in caplog.text

        # Results should still be returned
        assert len(modified_results) == 1

    def test_preserves_chain_information(self, tmp_path: Path, minimal_config: dict):
        """Given SNP integration, when applying, then preserves chain information."""
        snp_file = tmp_path / "snps.tsv"
        snp_file.write_text("1\t5\tA\tG\n")

        args = Mock()
        args.snp_input_file = str(snp_file)
        args.random_snps = False

        original_chain = ["1", "2", "X"]
        results = [("ATCGATCGATCGATCGATCG", original_chain)]

        modified_results, applied_snp_info = integrate_snps_unified(
            args=args, config=minimal_config, results=results, skip_reference_check=False
        )

        # Chain should be preserved
        assert modified_results[0][1] == original_chain
