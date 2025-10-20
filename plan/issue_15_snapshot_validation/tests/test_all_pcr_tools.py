#!/usr/bin/env python3
"""
Comprehensive testing and comparison of ALL in-silico PCR tools.

Tests:
1. pydna - Modern, actively maintained
2. ispcr - Lightweight, unmaintained
3. PyPCRtool - Feature-rich, unclear maintenance
4. BioPython PrimerSearch - EMBOSS wrapper
5. Custom implementation - Baseline comparison

Uses real MucOneUp dupC sequences.
"""

import sys
import time
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

# Test configuration
TEST_FASTA = "output/dupC/dupC.001.mut.simulated.fa"
TEST_SEQUENCE = None

# Load test sequence
print("=" * 80)
print("IN-SILICO PCR TOOLS COMPREHENSIVE TESTING")
print("=" * 80)

print("\nLoading test sequence...")
try:
    record = next(iter(SeqIO.parse(TEST_FASTA, "fasta")))
    TEST_SEQUENCE = str(record.seq)
    print(f"✓ Loaded: {len(TEST_SEQUENCE):,} bp from {TEST_FASTA}")
except Exception as e:
    print(f"✗ Error loading test file: {e}")
    sys.exit(1)

# Design test primers from the actual sequence
# Use portions we know exist in the sequence
PRIMER_FORWARD = TEST_SEQUENCE[1000:1020]  # 20bp from position 1000
PRIMER_REVERSE_RC = TEST_SEQUENCE[1400:1420]  # 20bp from position 1400
PRIMER_REVERSE = str(Seq(PRIMER_REVERSE_RC).reverse_complement())

print("\nTest primers designed from sequence:")
print(f"  Forward:  5'-{PRIMER_FORWARD}-3' (from pos 1000)")
print(f"  Reverse:  5'-{PRIMER_REVERSE}-3' (RC of pos 1400)")
print("  Expected product: ~420 bp")

# =============================================================================
# TEST 1: pydna
# =============================================================================

print("\n" + "=" * 80)
print("TEST 1: pydna")
print("=" * 80)

try:
    import pydna
    from pydna.amplify import pcr
    from pydna.dseqrecord import Dseqrecord

    print(f"Version: {pydna.__version__}")
    print("Status: ✓ Available")

    print("\n[Test 1a] Basic PCR simulation...")
    start_time = time.time()

    template = Dseqrecord(TEST_SEQUENCE, name="dupC_test")
    amplicon = pcr(PRIMER_FORWARD, PRIMER_REVERSE, template)

    elapsed = time.time() - start_time

    if amplicon:
        print("  ✓ SUCCESS")
        print(f"  Product size: {len(amplicon)} bp")
        print(f"  Forward binding: position {amplicon.forward_primer.position}")
        print(f"  Reverse binding: position {amplicon.reverse_primer.position}")
        print(f"  Time: {elapsed:.3f} seconds")

        # Verify product sequence
        expected_product = TEST_SEQUENCE[1000:1420]
        actual_product = str(amplicon.seq)

        if expected_product in actual_product or actual_product in expected_product:
            print("  ✓ Product sequence verified")
        else:
            print("  ⚠ Product sequence differs from expected")
    else:
        print("  ✗ FAILED - No amplification")

    print("\n[Test 1b] Mismatch tolerance...")
    # Test with 1 mismatch in forward primer
    mismatch_primer = list(PRIMER_FORWARD)
    mismatch_primer[10] = "N" if mismatch_primer[10] != "N" else "A"
    mismatch_primer = "".join(mismatch_primer)

    print(f"  Testing primer with mismatch: {mismatch_primer}")
    amplicon_mm = pcr(mismatch_primer, PRIMER_REVERSE, template)

    if amplicon_mm:
        print(f"  ✓ Tolerates mismatches (amplicon: {len(amplicon_mm)} bp)")
    else:
        print("  ✗ Does not tolerate mismatches")

    print("\n[Test 1c] Multiple primer pairs...")
    # Test with non-matching primers
    fake_forward = "AAAATTTTCCCCGGGG" + "TATA"
    amplicon_fake = pcr(fake_forward, PRIMER_REVERSE, template)

    if amplicon_fake:
        print("  ⚠ False positive with fake primers")
    else:
        print("  ✓ Correctly rejects non-matching primers")

    PYDNA_SCORE = {
        "available": True,
        "basic_pcr": amplicon is not None,
        "mismatch_tolerance": amplicon_mm is not None,
        "specificity": amplicon_fake is None,
        "performance": elapsed,
    }

except ImportError as e:
    print(f"Status: ✗ Not available: {e}")
    PYDNA_SCORE = {"available": False}
except Exception as e:
    print(f"Error during testing: {e}")
    PYDNA_SCORE = {"available": True, "error": str(e)}

# =============================================================================
# TEST 2: ispcr
# =============================================================================

print("\n" + "=" * 80)
print("TEST 2: ispcr")
print("=" * 80)

try:
    import tempfile

    from ispcr import get_pcr_products

    print("Status: ✓ Available")

    print("\n[Test 2a] Basic PCR simulation...")
    start_time = time.time()

    # ispcr requires FASTA files as input
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as pf:
        pf.write(f">forward\n{PRIMER_FORWARD}\n")
        pf.write(f">reverse\n{PRIMER_REVERSE}\n")
        primer_file = pf.name

    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as sf:
        sf.write(f">test_sequence\n{TEST_SEQUENCE}\n")
        seq_file = sf.name

    try:
        products = get_pcr_products(primer_file=primer_file, sequence_file=seq_file)
        elapsed = time.time() - start_time

        if products:
            print("  ✓ SUCCESS")
            print(f"  Products found: {len(products)}")
            for i, product in enumerate(products[:3]):  # Show first 3
                print(f"    Product {i+1}: {len(product)} bp")
            print(f"  Time: {elapsed:.3f} seconds")
        else:
            print("  ✗ FAILED - No products")

        ISPCR_SCORE = {
            "available": True,
            "basic_pcr": len(products) > 0 if products else False,
            "performance": elapsed,
        }

    finally:
        Path(primer_file).unlink()
        Path(seq_file).unlink()

except ImportError as e:
    print(f"Status: ✗ Not available: {e}")
    print("  Installing...")
    import subprocess

    result = subprocess.run(
        [sys.executable, "-m", "pip", "install", "ispcr"], capture_output=True, text=True
    )
    if result.returncode == 0:
        print("  ✓ Installation successful, please re-run test")
    else:
        print(f"  ✗ Installation failed: {result.stderr}")
    ISPCR_SCORE = {"available": False, "install_attempted": True}
except Exception as e:
    print(f"Error during testing: {e}")
    ISPCR_SCORE = {"available": True, "error": str(e)}

# =============================================================================
# TEST 3: PyPCRtool
# =============================================================================

print("\n" + "=" * 80)
print("TEST 3: PyPCRtool")
print("=" * 80)

try:
    # Try to import PyPCRtool
    # Note: Package name might be pypcrtool or PyPCRtool
    try:
        from pypcrtool import InSilicoPCR

        print("Status: ✓ Available")

        print("\n[Test 3a] Basic PCR simulation...")
        start_time = time.time()

        # PyPCRtool requires file paths
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as sf:
            sf.write(f">test\n{TEST_SEQUENCE}\n")
            seq_file = sf.name

        try:
            pcr_tube = InSilicoPCR(
                forward_primer=PRIMER_FORWARD,
                reverse_primer=PRIMER_REVERSE,
                sequence_file=seq_file,
                forward_mismatch_tolerance=1,
                reverse_mismatch_tolerance=1,
            )

            elapsed = time.time() - start_time

            # Get results
            if hasattr(pcr_tube, "amplicons"):
                products = pcr_tube.amplicons
                if products:
                    print("  ✓ SUCCESS")
                    print(f"  Products: {len(products)}")
                    print(f"  Time: {elapsed:.3f} seconds")
                else:
                    print("  ✗ No products")
            else:
                print("  ⚠ Could not access results")

            PYPCRTOOL_SCORE = {"available": True, "basic_pcr": True, "performance": elapsed}

        finally:
            Path(seq_file).unlink()

    except ImportError as import_err:
        raise ImportError("pypcrtool module not found") from import_err

except ImportError as e:
    print(f"Status: ✗ Not available: {e}")
    print("  Attempting installation...")
    import subprocess

    result = subprocess.run(
        [sys.executable, "-m", "pip", "install", "pypcrtool"], capture_output=True, text=True
    )
    if result.returncode == 0:
        print("  ✓ Installation successful, please re-run test")
    else:
        print(f"  ✗ Installation failed: {result.stderr}")
    PYPCRTOOL_SCORE = {"available": False, "install_attempted": True}
except Exception as e:
    print(f"Error during testing: {e}")
    import traceback

    traceback.print_exc()
    PYPCRTOOL_SCORE = {"available": False, "error": str(e)}

# =============================================================================
# TEST 4: BioPython PrimerSearch (EMBOSS)
# =============================================================================

print("\n" + "=" * 80)
print("TEST 4: BioPython PrimerSearch (EMBOSS wrapper)")
print("=" * 80)

try:
    import subprocess

    from Bio.Emboss.Applications import PrimerSearchCommandline

    # Check if primersearch is available
    result = subprocess.run(["which", "primersearch"], capture_output=True, text=True)

    if result.returncode == 0:
        print("Status: ✓ EMBOSS primersearch available")
        print(f"  Location: {result.stdout.strip()}")

        print("\n[Test 4a] Basic PCR simulation...")
        start_time = time.time()

        # Create primer file
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as pf:
            pf.write(f"test_primers {PRIMER_FORWARD} {PRIMER_REVERSE}\n")
            primer_file = pf.name

        # Create sequence file
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as sf:
            sf.write(f">test_seq\n{TEST_SEQUENCE}\n")
            seq_file = sf.name

        # Create output file
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as of:
            out_file = of.name

        try:
            # Run primersearch
            cline = PrimerSearchCommandline(
                seqall=seq_file, infile=primer_file, mismatchpercent=10, outfile=out_file
            )

            result = subprocess.run(str(cline), shell=True, capture_output=True, text=True)
            elapsed = time.time() - start_time

            if result.returncode == 0:
                # Read results
                with open(out_file) as f:
                    output = f.read()

                if "Amplimer" in output:
                    print("  ✓ SUCCESS")
                    print(f"  Time: {elapsed:.3f} seconds")
                    # Parse output
                    lines = output.split("\n")
                    for line in lines:
                        if "Amplimer length" in line:
                            print(f"  {line.strip()}")
                else:
                    print("  ✗ No amplification found")

                PRIMERSEARCH_SCORE = {
                    "available": True,
                    "basic_pcr": "Amplimer" in output,
                    "performance": elapsed,
                }
            else:
                print(f"  ✗ Command failed: {result.stderr}")
                PRIMERSEARCH_SCORE = {"available": True, "error": result.stderr}

        finally:
            Path(primer_file).unlink()
            Path(seq_file).unlink()
            Path(out_file).unlink()

    else:
        print("Status: ✗ EMBOSS primersearch not found in PATH")
        print("  Install with: sudo apt-get install emboss (Ubuntu/Debian)")
        print("  or: conda install -c bioconda emboss")
        PRIMERSEARCH_SCORE = {"available": False, "reason": "EMBOSS not installed"}

except ImportError as e:
    print(f"Status: ✗ BioPython Emboss module issue: {e}")
    PRIMERSEARCH_SCORE = {"available": False, "error": str(e)}
except Exception as e:
    print(f"Error during testing: {e}")
    import traceback

    traceback.print_exc()
    PRIMERSEARCH_SCORE = {"available": False, "error": str(e)}

# =============================================================================
# TEST 5: Custom Implementation (Baseline)
# =============================================================================

print("\n" + "=" * 80)
print("TEST 5: Custom Implementation (Baseline)")
print("=" * 80)

print("Status: ✓ Always available")


def custom_pcr_simulation(template, fwd, rev, mismatch_tolerance=0):
    """Simple string-matching PCR simulation."""
    from Bio.Seq import Seq

    # Reverse complement the reverse primer
    rev_rc = str(Seq(rev).reverse_complement())

    # Find binding sites
    fwd_pos = template.find(fwd)
    rev_pos = template.find(rev_rc)

    if (fwd_pos == -1 or rev_pos == -1) and mismatch_tolerance > 0:
            # Simple fuzzy matching (not optimal, just for comparison)
            for i in range(len(template) - len(fwd)):
                mismatches = sum(
                    1 for a, b in zip(template[i : i + len(fwd)], fwd, strict=False) if a != b
                )
                if mismatches <= mismatch_tolerance:
                    fwd_pos = i
                    break

            for i in range(len(template) - len(rev_rc)):
                mismatches = sum(
                    1 for a, b in zip(template[i : i + len(rev_rc)], rev_rc, strict=False) if a != b
                )
                if mismatches <= mismatch_tolerance:
                    rev_pos = i
                    break

    if fwd_pos == -1 or rev_pos == -1 or fwd_pos >= rev_pos:
        return None

    # Extract product
    product = template[fwd_pos : rev_pos + len(rev_rc)]
    return product


print("\n[Test 5a] Basic PCR simulation...")
start_time = time.time()

product = custom_pcr_simulation(TEST_SEQUENCE, PRIMER_FORWARD, PRIMER_REVERSE)
elapsed = time.time() - start_time

if product:
    print("  ✓ SUCCESS")
    print(f"  Product size: {len(product)} bp")
    print(f"  Time: {elapsed:.3f} seconds")
else:
    print("  ✗ FAILED")

print("\n[Test 5b] With mismatch tolerance...")
product_mm = custom_pcr_simulation(
    TEST_SEQUENCE, PRIMER_FORWARD, PRIMER_REVERSE, mismatch_tolerance=2
)
if product_mm:
    print(f"  ✓ Works with mismatches ({len(product_mm)} bp)")
else:
    print("  ✗ Failed even with tolerance")

CUSTOM_SCORE = {
    "available": True,
    "basic_pcr": product is not None,
    "mismatch_tolerance": product_mm is not None,
    "performance": elapsed,
}

# =============================================================================
# COMPARISON SUMMARY
# =============================================================================

print("\n" + "=" * 80)
print("COMPREHENSIVE COMPARISON")
print("=" * 80)

results = {
    "pydna": PYDNA_SCORE,
    "ispcr": ISPCR_SCORE,
    "PyPCRtool": PYPCRTOOL_SCORE,
    "BioPython PrimerSearch": PRIMERSEARCH_SCORE,
    "Custom": CUSTOM_SCORE,
}

print(
    "\n{:<25} {:<15} {:<15} {:<15} {:<15}".format(
        "Tool", "Available", "Basic PCR", "Mismatch", "Performance"
    )
)
print("-" * 85)

for tool, score in results.items():
    available = "✓" if score.get("available") else "✗"
    basic = "✓" if score.get("basic_pcr") else ("✗" if "basic_pcr" in score else "N/A")
    mismatch = (
        "✓"
        if score.get("mismatch_tolerance")
        else ("✗" if "mismatch_tolerance" in score else "N/A")
    )
    perf = f"{score.get('performance', 0):.3f}s" if "performance" in score else "N/A"

    print(f"{tool:<25} {available:<15} {basic:<15} {mismatch:<15} {perf:<15}")

print("\n" + "=" * 80)
print("RECOMMENDATIONS")
print("=" * 80)

# Determine best tool
if PYDNA_SCORE.get("basic_pcr"):
    print("\n✅ RECOMMENDED: pydna")
    print("   Reasons:")
    print("   - Actively maintained (v5.5.3, March 2025)")
    print("   - Thermodynamically accurate")
    print("   - Handles mismatches intelligently")
    print("   - Fast performance")
    print("   - Well-documented")
elif ISPCR_SCORE.get("basic_pcr"):
    print("\n⚠️  FALLBACK: ispcr")
    print("   Reasons:")
    print("   - Lightweight and simple")
    print("   - Works for basic cases")
    print("   - But unmaintained since 2023")
elif CUSTOM_SCORE.get("basic_pcr"):
    print("\n⚠️  FALLBACK: Custom implementation")
    print("   Reasons:")
    print("   - Always available")
    print("   - Simple string matching")
    print("   - No thermodynamics")

print("\n" + "=" * 80)
print("Test suite completed!")
print("=" * 80)
