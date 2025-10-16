#!/usr/bin/env python3
"""
Coverage Ratchet Checker for Phase 2 Testing & Security.

This script ensures test coverage doesn't decrease during Phase 2 implementation.
It compares current coverage against a baseline and fails if any module's coverage
drops below the baseline (with a small tolerance for floating-point variations).

Usage:
    python scripts/check_coverage_ratchet.py
    python scripts/check_coverage_ratchet.py --baseline coverage.baseline.json --current coverage.json
    python scripts/check_coverage_ratchet.py --tolerance 1.0
"""

import argparse
import json
import sys
from pathlib import Path


def check_ratchet(baseline_path: str, current_path: str, tolerance: float = 0.5) -> bool:
    """
    Check if coverage has regressed compared to baseline.

    Args:
        baseline_path: Path to baseline coverage JSON file
        current_path: Path to current coverage JSON file
        tolerance: Allowed decrease in coverage percentage (default: 0.5%)

    Returns:
        True if no regressions detected, False otherwise

    Raises:
        FileNotFoundError: If baseline or current coverage file not found
        json.JSONDecodeError: If coverage files are invalid JSON
    """
    # Load coverage files
    baseline_file = Path(baseline_path)
    current_file = Path(current_path)

    if not baseline_file.exists():
        print(f"❌ ERROR: Baseline file not found: {baseline_path}")
        print("   Run 'bash scripts/capture_coverage_baseline.sh' first")
        sys.exit(1)

    if not current_file.exists():
        print(f"❌ ERROR: Current coverage file not found: {current_path}")
        print("   Run 'pytest --cov=muc_one_up --cov-report=json:coverage.json' first")
        sys.exit(1)

    try:
        baseline = json.loads(baseline_file.read_text())
        current = json.loads(current_file.read_text())
    except json.JSONDecodeError as e:
        print(f"❌ ERROR: Invalid JSON in coverage file: {e}")
        sys.exit(1)

    # Check overall coverage
    baseline_total = baseline["totals"]["percent_covered"]
    current_total = current["totals"]["percent_covered"]

    # Check per-file coverage
    regressions = []
    improvements = []

    for file_path, current_data in current["files"].items():
        baseline_data = baseline["files"].get(file_path, {})
        baseline_cov = baseline_data.get("summary", {}).get("percent_covered", 0)
        current_cov = current_data["summary"]["percent_covered"]

        # Calculate difference
        diff = current_cov - baseline_cov

        if current_cov < baseline_cov - tolerance:
            regressions.append(
                {"file": file_path, "baseline": baseline_cov, "current": current_cov, "diff": diff}
            )
        elif current_cov > baseline_cov + tolerance:
            improvements.append(
                {"file": file_path, "baseline": baseline_cov, "current": current_cov, "diff": diff}
            )

    # Report results
    print("=" * 70)
    print("📊 Coverage Ratchet Check")
    print("=" * 70)
    print(f"Baseline:  {baseline_total:.1f}%")
    print(f"Current:   {current_total:.1f}%")
    print(f"Change:    {current_total - baseline_total:+.1f}%")
    print(f"Tolerance: ±{tolerance:.1f}%")
    print("=" * 70)

    if regressions:
        print("")
        print("❌ COVERAGE REGRESSIONS DETECTED:")
        print("")
        for reg in regressions:
            print(f"  {reg['file']}")
            print(f"    {reg['baseline']:.1f}% → {reg['current']:.1f}% ({reg['diff']:.1f}%)")
        print("")
        print("Coverage must not decrease during Phase 2 implementation.")
        print("Please add tests to restore coverage to baseline levels.")
        return False

    print("")
    if improvements:
        print("✅ COVERAGE IMPROVEMENTS:")
        print("")
        for imp in improvements[:10]:  # Show top 10
            print(f"  {imp['file']}")
            print(f"    {imp['baseline']:.1f}% → {imp['current']:.1f}% ({imp['diff']:+.1f}%)")
        if len(improvements) > 10:
            print(f"  ... and {len(improvements) - 10} more files improved")
        print("")

    print("✅ No coverage regressions detected!")
    print(f"   Overall coverage: {current_total:.1f}%")
    print("")

    return True


def main():
    """Main entry point for coverage ratchet checker."""
    parser = argparse.ArgumentParser(
        description="Check that test coverage hasn't regressed",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Use default paths
  python scripts/check_coverage_ratchet.py

  # Specify custom paths
  python scripts/check_coverage_ratchet.py --baseline old.json --current new.json

  # Allow 1% tolerance
  python scripts/check_coverage_ratchet.py --tolerance 1.0
        """,
    )

    parser.add_argument(
        "--baseline",
        default="coverage.baseline.json",
        help="Path to baseline coverage JSON file (default: coverage.baseline.json)",
    )

    parser.add_argument(
        "--current",
        default="coverage.json",
        help="Path to current coverage JSON file (default: coverage.json)",
    )

    parser.add_argument(
        "--tolerance",
        type=float,
        default=0.5,
        help="Allowed decrease in coverage percentage (default: 0.5)",
    )

    args = parser.parse_args()

    # Run ratchet check
    success = check_ratchet(
        baseline_path=args.baseline, current_path=args.current, tolerance=args.tolerance
    )

    # Exit with appropriate code
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
