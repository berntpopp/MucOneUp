#!/usr/bin/env python3
"""
Comprehensive test suite for SNaPshot validation using REAL MucOneUp output.

Tests the complete workflow:
1. PCR simulation (pydna + fallback)
2. Restriction digest (BioPython)
3. Extension primer validation (primer3-py)
4. End-to-end integration

Uses actual dupC mutation sequences from MucOneUp.
"""

import sys
from pathlib import Path

# Test configuration
TEST_FILES = {
    "dupC_mutant": Path("output/dupC/dupC.001.mut.simulated.fa"),
    "dupC_normal": Path("output/dupC/dupC.001.normal.simulated.fa"),
}

# Expected mutation position for dupC (from config: insert at position 60)
# This is within the X repeat, position 60-61
DUPC_MUTATION_REGION = (4000, 5000)  # Approximate VNTR region

# PCR primers (designed to flank MUC1 VNTR)
PCR_PRIMERS = {
    "forward": "GGTGAGCTCAAGGGGGGCTG",  # From left constant region
    "reverse": "TGTGAGCCACCGTGCCACGC",  # From right constant region
}

# Extension primer (designed to end just before mutation)
EXTENSION_PRIMER = "CCCGCCCGCCCGGGGCCCGG"

print("=" * 80)
print("SNAPSHOT VALIDATION TEST SUITE - REAL DATA")
print("=" * 80)

# ============================================================================
# TEST 1: Check Dependencies
# ============================================================================

print("\n[TEST 1] Checking Dependencies...")
print("-" * 80)

dependencies_ok = True

try:
    import Bio
    from Bio import Seq
    from Bio.Restriction import Restriction

    print(f"✓ BioPython {Bio.__version__}")
except ImportError as e:
    print(f"✗ BioPython not available: {e}")
    dependencies_ok = False

try:
    import primer3

    print("✓ primer3-py available")
except ImportError as e:
    print(f"✗ primer3-py not available: {e}")
    dependencies_ok = False

try:
    import pydna
    from pydna.amplify import pcr
    from pydna.dseqrecord import Dseqrecord

    print(f"✓ pydna {pydna.__version__}")
    PYDNA_AVAILABLE = True
except ImportError as e:
    print(f"⚠ pydna not available (will use fallback): {e}")
    PYDNA_AVAILABLE = False

if not dependencies_ok:
    print("\n✗ CRITICAL: Missing required dependencies")
    print("Install with: pip install biopython primer3-py pydna")
    sys.exit(1)

print("\n✓ All critical dependencies available")

# ============================================================================
# TEST 2: Load Real FASTA Files
# ============================================================================

print("\n[TEST 2] Loading Real MucOneUp FASTA Files...")
print("-" * 80)


def load_fasta(filepath):
    """Load first sequence from FASTA file."""
    from Bio import SeqIO

    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"File not found: {filepath}")

    records = list(SeqIO.parse(filepath, "fasta"))
    if not records:
        raise ValueError(f"No sequences in {filepath}")

    return str(records[0].seq)


try:
    mutant_seq = load_fasta(TEST_FILES["dupC_mutant"])
    normal_seq = load_fasta(TEST_FILES["dupC_normal"])

    print(f"✓ Loaded dupC mutant:  {len(mutant_seq):,} bp")
    print(f"✓ Loaded dupC normal:  {len(normal_seq):,} bp")

    # Find the difference
    diff_positions = [
        i for i in range(min(len(mutant_seq), len(normal_seq))) if mutant_seq[i] != normal_seq[i]
    ]

    if diff_positions:
        print(f"\n✓ Found {len(diff_positions)} difference(s) between sequences")
        print(f"  First difference at position: {diff_positions[0]}")

        # Show the mutation context
        pos = diff_positions[0]
        context_start = max(0, pos - 20)
        context_end = min(len(mutant_seq), pos + 20)

        print(f"\n  Normal:  ...{normal_seq[context_start:context_end]}...")
        print(f"  Mutant:  ...{mutant_seq[context_start:context_end]}...")
    else:
        print("\n⚠  WARNING: No differences found between normal and mutant!")

except Exception as e:
    print(f"\n✗ ERROR loading files: {e}")
    print("\nExpected files:")
    for name, path in TEST_FILES.items():
        print(f"  {name}: {path}")
    sys.exit(1)

# ============================================================================
# TEST 3: PCR Simulation with pydna
# ============================================================================

print("\n[TEST 3] PCR Simulation with pydna...")
print("-" * 80)

if PYDNA_AVAILABLE:
    from pydna.amplify import pcr
    from pydna.dseqrecord import Dseqrecord

    def test_pcr_pydna(sequence, fwd, rev, name):
        """Test PCR with pydna."""
        print(f"\nTesting {name}:")
        print(f"  Template: {len(sequence):,} bp")
        print(f"  Forward:  5'-{fwd}-3'")
        print(f"  Reverse:  5'-{rev}-3'")

        try:
            template = Dseqrecord(sequence)
            amplicon = pcr(fwd, rev, template)

            if amplicon:
                print(f"  ✓ Amplified: {len(amplicon)} bp")
                print(f"    Fwd binding: {amplicon.forward_primer.position}")
                print(f"    Rev binding: {amplicon.reverse_primer.position}")
                return str(amplicon.seq)
            else:
                print("  ✗ No amplification")
                return None
        except Exception as e:
            print(f"  ✗ Error: {e}")
            return None

    mutant_amplicon = test_pcr_pydna(
        mutant_seq, PCR_PRIMERS["forward"], PCR_PRIMERS["reverse"], "Mutant"
    )

    normal_amplicon = test_pcr_pydna(
        normal_seq, PCR_PRIMERS["forward"], PCR_PRIMERS["reverse"], "Normal"
    )

    if mutant_amplicon and normal_amplicon:
        print("\n✓ Both sequences amplified successfully")
        print(f"  Mutant amplicon: {len(mutant_amplicon):,} bp")
        print(f"  Normal amplicon: {len(normal_amplicon):,} bp")
        print(f"  Size difference: {abs(len(mutant_amplicon) - len(normal_amplicon))} bp")
    else:
        print("\n⚠  PCR amplification incomplete")

else:
    print("⚠  Skipping pydna test (not installed)")
    mutant_amplicon = None
    normal_amplicon = None

# ============================================================================
# TEST 4: Restriction Digest Simulation
# ============================================================================

print("\n[TEST 4] Restriction Digest with BioPython...")
print("-" * 80)

from Bio.Restriction import Restriction
from Bio.Seq import Seq


def test_restriction_digest(sequence, enzyme_name, seq_name):
    """Test restriction enzyme digest."""
    print(f"\nTesting {seq_name} with {enzyme_name}:")
    print(f"  Sequence length: {len(sequence):,} bp")

    try:
        # Get enzyme
        enzyme = Restriction.MwoI if enzyme_name == "MwoI" else getattr(Restriction, enzyme_name)

        seq_obj = Seq(sequence)

        # Search for sites
        sites = enzyme.search(seq_obj)

        print(f"  {enzyme_name} recognition: {enzyme.site}")
        print(f"  Sites found: {len(sites)}")

        if sites:
            print(f"  Site positions: {sites[:5]}{'...' if len(sites) > 5 else ''}")

            # Perform digest
            fragments = enzyme.catalyse(seq_obj)
            frag_sizes = [len(f) for f in fragments]
            frag_sizes.sort(reverse=True)

            print(f"  Fragments: {len(fragments)}")
            print(f"  Fragment sizes: {frag_sizes[:5]}{'...' if len(frag_sizes) > 5 else ''}")

            return {
                "cuts": True,
                "num_sites": len(sites),
                "num_fragments": len(fragments),
                "largest_fragment": max(frag_sizes),
            }
        else:
            print(f"  ✓ No {enzyme_name} sites - sequence intact")
            return {
                "cuts": False,
                "num_sites": 0,
                "num_fragments": 1,
                "largest_fragment": len(sequence),
            }

    except Exception as e:
        print(f"  ✗ Error: {e}")
        return None


# Test on full sequences first
print("\n--- Testing Full Sequences ---")
mutant_digest = test_restriction_digest(mutant_seq, "MwoI", "Mutant (full)")
normal_digest = test_restriction_digest(normal_seq, "MwoI", "Normal (full)")

# Test on amplicons if available
if mutant_amplicon and normal_amplicon:
    print("\n--- Testing PCR Amplicons ---")
    mutant_amplicon_digest = test_restriction_digest(mutant_amplicon, "MwoI", "Mutant (amplicon)")
    normal_amplicon_digest = test_restriction_digest(normal_amplicon, "MwoI", "Normal (amplicon)")

    print("\n✓ Comparison:")
    print(f"  Mutant: {mutant_amplicon_digest['num_sites']} sites")
    print(f"  Normal: {normal_amplicon_digest['num_sites']} sites")

    if mutant_amplicon_digest["num_sites"] < normal_amplicon_digest["num_sites"]:
        print("  ✓ Mutant has fewer MwoI sites (mutation disrupts restriction)")
    elif mutant_amplicon_digest["num_sites"] > normal_amplicon_digest["num_sites"]:
        print("  ⚠  Mutant has MORE MwoI sites (unexpected)")
    else:
        print("  ⚠  Same number of sites (mutation may not affect MwoI)")

# ============================================================================
# TEST 5: Extension Primer Validation
# ============================================================================

print("\n[TEST 5] Extension Primer Validation with primer3-py...")
print("-" * 80)

import primer3


def test_extension_primer(primer_seq):
    """Validate extension primer thermodynamics."""
    print(f"\nExtension Primer: 5'-{primer_seq}-3'")
    print(f"Length: {len(primer_seq)} bp")

    # Calculate Tm
    tm = primer3.calc_tm(primer_seq, mv_conc=50.0, dv_conc=1.5, dntp_conc=0.6, dna_conc=50.0)

    print(f"Melting Temp (Tm): {tm:.2f}°C")

    # Check hairpins
    hairpin = primer3.calc_hairpin(primer_seq)

    if hairpin.structure_found:
        print("Hairpin detected:")
        print(f"  Tm: {hairpin.tm:.2f}°C")
        print(f"  ΔG: {hairpin.dg:.2f} kcal/mol")

        if hairpin.tm > 45:
            print("  ⚠  High hairpin Tm may interfere")
            hairpin_ok = False
        else:
            print("  ✓  Hairpin Tm acceptable")
            hairpin_ok = True
    else:
        print("✓ No hairpin detected")
        hairpin_ok = True

    # GC content
    gc_count = primer_seq.count("G") + primer_seq.count("C")
    gc_percent = (gc_count / len(primer_seq)) * 100
    print(f"GC content: {gc_percent:.1f}%")

    # Validate
    tm_ok = 50.0 <= tm <= 65.0

    print("\nValidation:")
    print(f"  Tm in range (50-65°C): {'✓' if tm_ok else '✗'}")
    print(f"  No problematic hairpin: {'✓' if hairpin_ok else '✗'}")
    print(f"  Overall: {'✓ PASS' if (tm_ok and hairpin_ok) else '✗ FAIL'}")

    return tm_ok and hairpin_ok


extension_valid = test_extension_primer(EXTENSION_PRIMER)

# ============================================================================
# TEST 6: End-to-End Integration Test
# ============================================================================

print("\n[TEST 6] End-to-End Integration Test...")
print("-" * 80)


def end_to_end_validation(sequence, seq_name, amplicon=None):
    """Complete SNaPshot validation."""
    print(f"\n=== Validating {seq_name} ===")

    # Use amplicon if provided, otherwise full sequence
    test_seq = amplicon if amplicon else sequence

    # Step 1: Check for MwoI sites
    from Bio.Restriction import MwoI
    from Bio.Seq import Seq

    seq_obj = Seq(test_seq)
    sites = MwoI.search(seq_obj)

    print(f"[1/3] MwoI digest: {len(sites)} sites")

    if sites:
        print("      → Cuts (wild-type phenotype)")
        accessible = False
        reason = "Wild-type MwoI site present - would be digested"
    else:
        print("      → No cut (mutant phenotype)")
        accessible = True
        reason = "Mutation disrupts MwoI site"

    # Step 2: Extension primer
    print(f"[2/3] Extension primer: {extension_valid}")

    # Step 3: Signal prediction
    # For dupC, expect additional C, so complement would be G
    expected_ddNTP = "ddGTP"
    expected_fluorophore = "Blue (dR110)"

    print(f"[3/3] Expected signal: {expected_ddNTP} ({expected_fluorophore})")

    # Final verdict
    print("\nResult:")
    print(f"  Accessible by SNaPshot: {accessible}")
    print(f"  Reason: {reason}")
    print(f"  Extension primer valid: {extension_valid}")
    print(f"  Expected detection: {expected_ddNTP}")

    return {
        "accessible": accessible,
        "reason": reason,
        "extension_valid": extension_valid,
        "expected_signal": expected_ddNTP,
    }


# Run for both sequences
mutant_result = end_to_end_validation(mutant_seq, "Mutant", mutant_amplicon)

normal_result = end_to_end_validation(normal_seq, "Normal", normal_amplicon)

# ============================================================================
# FINAL SUMMARY
# ============================================================================

print("\n" + "=" * 80)
print("FINAL SUMMARY")
print("=" * 80)

print("\nDependencies:")
print("  ✓ BioPython: Available")
print("  ✓ primer3-py: Available")
print(
    f"  {'✓' if PYDNA_AVAILABLE else '⚠'} pydna: {'Available' if PYDNA_AVAILABLE else 'Not available (fallback mode)'}"
)

print("\nTest Files:")
print(f"  ✓ Mutant sequence: {len(mutant_seq):,} bp")
print(f"  ✓ Normal sequence: {len(normal_seq):,} bp")

if PYDNA_AVAILABLE and mutant_amplicon and normal_amplicon:
    print("\nPCR Simulation:")
    print(f"  ✓ Mutant amplicon: {len(mutant_amplicon):,} bp")
    print(f"  ✓ Normal amplicon: {len(normal_amplicon):,} bp")

print("\nSNaPshot Validation Results:")
print(f"  Mutant accessible: {mutant_result['accessible']}")
print(f"  Normal accessible: {normal_result['accessible']}")

if mutant_result["accessible"] and not normal_result["accessible"]:
    print("\n✅ EXPECTED RESULT: Mutant detectable, normal not detectable")
    print("   This is the correct behavior for dupC mutation validation")
elif not mutant_result["accessible"] and not normal_result["accessible"]:
    print("\n⚠️  UNEXPECTED: Neither sequence would be detectable")
    print("   May need different primers or enzyme")
else:
    print("\n⚠️  REVIEW NEEDED: Results may indicate issues with validation approach")

print("\n" + "=" * 80)
print("Test suite completed!")
print("=" * 80)
