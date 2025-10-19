#!/usr/bin/env python3
"""
Allelic Balance Analyzer for NanoSim read simulations.

This module analyzes BAM files from diploid read simulations to detect
allelic imbalance between haplotypes. It calculates per-haplotype coverage,
read counts, and performs statistical tests to determine if coverage deviates
from the expected 50:50 ratio.

Key features:
- Extract coverage per haplotype from aligned BAM files
- Calculate allelic ratios and deviations
- Perform statistical tests (Chi-square, binomial)
- Generate detailed reports for scientific analysis

Author: Claude Code (for Issue #31 empirical testing)
Date: 2025-10-19
"""

import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pysam
from scipy.stats import chisquare

# Handle scipy version compatibility
try:
    from scipy.stats import binomtest
except ImportError:
    # Fallback for older scipy versions
    from scipy.stats import binom_test as binomtest


@dataclass
class HaplotypeCoverage:
    """Container for haplotype coverage metrics."""

    name: str
    read_count: int
    total_bases: int
    reference_length: int
    mean_coverage: float

    @property
    def coverage(self) -> float:
        """Alias for mean_coverage for backwards compatibility."""
        return self.mean_coverage


@dataclass
class AllelicBalanceReport:
    """Container for allelic balance analysis results."""

    test_id: str
    haplotype_1: HaplotypeCoverage
    haplotype_2: HaplotypeCoverage
    total_reads: int
    allelic_ratio: float  # ratio of hap1 / (hap1 + hap2)
    expected_ratio: float  # expected ratio (default 0.5)
    deviation_percent: float  # absolute deviation from expected
    chi_square_statistic: float
    chi_square_pvalue: float
    binomial_pvalue: float
    is_balanced: bool  # True if within tolerance
    tolerance: float

    def to_dict(self) -> dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return {
            "test_id": self.test_id,
            "haplotype_1": {
                "name": self.haplotype_1.name,
                "read_count": self.haplotype_1.read_count,
                "total_bases": self.haplotype_1.total_bases,
                "reference_length": self.haplotype_1.reference_length,
                "mean_coverage": self.haplotype_1.mean_coverage,
            },
            "haplotype_2": {
                "name": self.haplotype_2.name,
                "read_count": self.haplotype_2.read_count,
                "total_bases": self.haplotype_2.total_bases,
                "reference_length": self.haplotype_2.reference_length,
                "mean_coverage": self.haplotype_2.mean_coverage,
            },
            "total_reads": self.total_reads,
            "allelic_ratio": self.allelic_ratio,
            "expected_ratio": self.expected_ratio,
            "deviation_percent": self.deviation_percent,
            "chi_square_statistic": self.chi_square_statistic,
            "chi_square_pvalue": self.chi_square_pvalue,
            "binomial_pvalue": self.binomial_pvalue,
            "is_balanced": self.is_balanced,
            "tolerance": self.tolerance,
        }


def calculate_haplotype_coverage(
    bam_file: str,
    haplotype_name: str,
    reference_length: int,
) -> HaplotypeCoverage:
    """
    Calculate coverage for a specific haplotype from aligned BAM.

    Args:
        bam_file: Path to indexed BAM file
        haplotype_name: Name of the haplotype reference (contig name)
        reference_length: Length of the haplotype reference sequence

    Returns:
        HaplotypeCoverage object with read count and coverage metrics

    Raises:
        FileNotFoundError: If BAM file or index doesn't exist
        ValueError: If haplotype not found in BAM
    """
    bam_path = Path(bam_file)
    if not bam_path.exists():
        raise FileNotFoundError(f"BAM file not found: {bam_file}")

    bai_path = Path(str(bam_file) + ".bai")
    if not bai_path.exists():
        raise FileNotFoundError(f"BAM index not found: {bai_path}")

    logging.info(f"Analyzing coverage for haplotype: {haplotype_name}")

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        # Check if haplotype exists in BAM
        if haplotype_name not in bam.references:
            available = ", ".join(bam.references)
            raise ValueError(
                f"Haplotype '{haplotype_name}' not found in BAM. "
                f"Available references: {available}"
            )

        # Count reads mapping to this haplotype
        read_count = 0
        total_bases = 0

        for read in bam.fetch(haplotype_name):
            if not read.is_unmapped and not read.is_secondary and not read.is_supplementary:
                read_count += 1
                # Count aligned bases (excludes soft/hard clipped)
                total_bases += read.reference_length

        # Calculate mean coverage
        mean_coverage = total_bases / reference_length if reference_length > 0 else 0.0

        logging.info(
            f"  Haplotype: {haplotype_name} - "
            f"Reads: {read_count}, "
            f"Total bases: {total_bases:,}, "
            f"Mean coverage: {mean_coverage:.2f}x"
        )

        return HaplotypeCoverage(
            name=haplotype_name,
            read_count=read_count,
            total_bases=total_bases,
            reference_length=reference_length,
            mean_coverage=mean_coverage,
        )


def analyze_allelic_balance(
    bam_file: str,
    haplotype_1_name: str,
    haplotype_2_name: str,
    haplotype_1_length: int,
    haplotype_2_length: int,
    test_id: str = "unknown",
    tolerance: float = 0.05,
    expected_ratio: float = 0.5,
) -> AllelicBalanceReport:
    """
    Analyze allelic balance between two haplotypes.

    Args:
        bam_file: Path to aligned BAM file containing reads from both haplotypes
        haplotype_1_name: Name of first haplotype (contig name in BAM)
        haplotype_2_name: Name of second haplotype (contig name in BAM)
        haplotype_1_length: Length of first haplotype reference
        haplotype_2_length: Length of second haplotype reference
        test_id: Identifier for this test case
        tolerance: Acceptable deviation from expected ratio (default 0.05 = 5%)
        expected_ratio: Expected allelic ratio (default 0.5 for diploid)

    Returns:
        AllelicBalanceReport with detailed metrics and statistical tests

    Raises:
        FileNotFoundError: If BAM file not found
        ValueError: If haplotypes not found in BAM
    """
    logging.info(f"=== Analyzing Allelic Balance: {test_id} ===")

    # Calculate coverage for each haplotype
    hap1_cov = calculate_haplotype_coverage(bam_file, haplotype_1_name, haplotype_1_length)
    hap2_cov = calculate_haplotype_coverage(bam_file, haplotype_2_name, haplotype_2_length)

    # Calculate metrics
    total_reads = hap1_cov.read_count + hap2_cov.read_count

    if total_reads == 0:
        logging.warning("No reads found for either haplotype!")
        allelic_ratio = 0.5  # Default if no reads
    else:
        allelic_ratio = hap1_cov.read_count / total_reads

    deviation_percent = abs(allelic_ratio - expected_ratio) * 100
    is_balanced = deviation_percent <= (tolerance * 100)

    # Statistical tests
    # Chi-square test: Are observed counts consistent with expected 50:50?
    observed = [hap1_cov.read_count, hap2_cov.read_count]
    expected = [total_reads * expected_ratio, total_reads * (1 - expected_ratio)]

    if total_reads > 0:
        chi2_stat, chi2_pval = chisquare(observed, expected)
        # Binomial test: Is the ratio consistent with p=0.5?
        # Handle scipy version compatibility
        if callable(binomtest) and not hasattr(binomtest, "pvalue"):
            # Old scipy: binom_test returns p-value directly
            binom_pval = binomtest(
                hap1_cov.read_count, total_reads, expected_ratio, alternative="two-sided"
            )
        else:
            # New scipy: binomtest returns object with pvalue attribute
            result = binomtest(
                hap1_cov.read_count, n=total_reads, p=expected_ratio, alternative="two-sided"
            )
            binom_pval = result.pvalue
    else:
        chi2_stat, chi2_pval = 0.0, 1.0
        binom_pval = 1.0

    report = AllelicBalanceReport(
        test_id=test_id,
        haplotype_1=hap1_cov,
        haplotype_2=hap2_cov,
        total_reads=total_reads,
        allelic_ratio=allelic_ratio,
        expected_ratio=expected_ratio,
        deviation_percent=deviation_percent,
        chi_square_statistic=chi2_stat,
        chi_square_pvalue=chi2_pval,
        binomial_pvalue=binom_pval,
        is_balanced=is_balanced,
        tolerance=tolerance,
    )

    # Log summary
    logging.info(f"Total reads: {total_reads}")
    logging.info(f"Allelic ratio (Hap1): {allelic_ratio:.4f} (expected: {expected_ratio:.4f})")
    logging.info(f"Deviation: {deviation_percent:.2f}%")
    logging.info(f"Chi-square p-value: {chi2_pval:.6f}")
    logging.info(f"Binomial p-value: {binom_pval:.6f}")
    logging.info(f"Balanced: {is_balanced} (tolerance: {tolerance*100}%)")

    return report


def save_report(report: AllelicBalanceReport, output_file: str) -> None:
    """
    Save allelic balance report to JSON file.

    Args:
        report: AllelicBalanceReport object
        output_file: Path to output JSON file
    """
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w") as f:
        json.dump(report.to_dict(), f, indent=2)

    logging.info(f"Report saved to: {output_file}")


def generate_summary_markdown(reports: list[AllelicBalanceReport], output_file: str) -> None:
    """
    Generate human-readable markdown summary of multiple test cases.

    Args:
        reports: List of AllelicBalanceReport objects
        output_file: Path to output markdown file
    """
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w") as f:
        f.write("# Allelic Balance Analysis Summary\n\n")
        f.write("## Test Results\n\n")

        f.write(
            "| Test ID | Hap1 Reads | Hap2 Reads | Allelic Ratio | Deviation | Balanced? | p-value |\n"
        )
        f.write(
            "|---------|-----------|-----------|---------------|-----------|-----------|----------|\n"
        )

        for report in reports:
            status = "✅" if report.is_balanced else "❌"
            f.write(
                f"| {report.test_id} | "
                f"{report.haplotype_1.read_count} | "
                f"{report.haplotype_2.read_count} | "
                f"{report.allelic_ratio:.4f} | "
                f"{report.deviation_percent:.2f}% | "
                f"{status} | "
                f"{report.binomial_pvalue:.4e} |\n"
            )

        f.write("\n## Detailed Metrics\n\n")
        for report in reports:
            f.write(f"### {report.test_id}\n\n")
            f.write(f"**Haplotype 1: {report.haplotype_1.name}**\n")
            f.write(f"- Reads: {report.haplotype_1.read_count}\n")
            f.write(f"- Coverage: {report.haplotype_1.mean_coverage:.2f}x\n")
            f.write(f"- Reference length: {report.haplotype_1.reference_length:,} bp\n\n")

            f.write(f"**Haplotype 2: {report.haplotype_2.name}**\n")
            f.write(f"- Reads: {report.haplotype_2.read_count}\n")
            f.write(f"- Coverage: {report.haplotype_2.mean_coverage:.2f}x\n")
            f.write(f"- Reference length: {report.haplotype_2.reference_length:,} bp\n\n")

            f.write("**Statistics**\n")
            f.write(f"- Chi-square statistic: {report.chi_square_statistic:.4f}\n")
            f.write(f"- Chi-square p-value: {report.chi_square_pvalue:.6f}\n")
            f.write(f"- Binomial p-value: {report.binomial_pvalue:.6f}\n\n")

    logging.info(f"Summary saved to: {output_file}")
