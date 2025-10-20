#!/usr/bin/env python3
"""
Validate SNaPshot protocol findings with actual MucOneUp dupC data.

Based on Protocol MUC1-Englisch_Arif.docx
Tests the complete workflow:
1. MwoI recognition of 7C vs 8C
2. Primer design from X repeat
3. PCR simulation with pydna
4. Complete workflow validation
"""

import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Restriction import MwoI
from Bio.Seq import Seq

print("=" * 80)
print("PROTOCOL VALIDATION TEST - Real MUC1 SNaPshot Workflow")
print("=" * 80)

# =============================================================================
# TEST 1: Verify MwoI Recognition Pattern
# =============================================================================

print("\n[TEST 1] MwoI Recognition: 7C vs 8C")
print("-" * 80)

# MwoI recognition: GCNNNNNNNGC (GC...7 bases...GC)

# Wild-type (7 Cs)
wt_7c_sequence = "CGGGCTCCACCGCCCCCCCA"  # 7 Cs in the middle

# Mutant (8 Cs) - dupC mutation
mut_8c_sequence = "CGGGCTCCACCGCCCCCCCCA"  # 8 Cs in the middle

print("\nSequences:")
print(f"  WT (7C):  {wt_7c_sequence}")
print(f"  Mut (8C): {mut_8c_sequence}")
print(f"  Difference: Position {wt_7c_sequence.index('GCCCCCCCA')} - extra C inserted")

# Test with longer context including MwoI sites
# MwoI needs: GC.......GC (7 bases between the GC pairs)

# Create test sequences with clear MwoI sites
test_wt = "AGCTAGCTAGCGGGCTCCACCGCCCCCCACAGCTAGCTAGC"
test_mut = "AGCTAGCTAGCGGGCTCCACCGCCCCCCCACAGCTAGCTAGC"

print("\nTest with flanking sequences:")
print(f"  WT:  {test_wt}")
print(f"  Mut: {test_mut}")

wt_seq = Seq(test_wt)
mut_seq = Seq(test_mut)

wt_sites = MwoI.search(wt_seq)
mut_sites = MwoI.search(mut_seq)

print("\nMwoI site search:")
print(f"  WT sites found: {len(wt_sites)} at positions {wt_sites if wt_sites else 'None'}")
print(f"  Mut sites found: {len(mut_sites)} at positions {mut_sites if mut_sites else 'None'}")

if len(wt_sites) != len(mut_sites):
    print("  ✓ Mutation changes MwoI recognition!")
else:
    print("  ⚠ Warning: Same number of sites in WT and mutant")

# =============================================================================
# TEST 2: Analyze X Repeat Structure
# =============================================================================

print("\n[TEST 2] X Repeat Structure Analysis")
print("-" * 80)

# X repeat from config
X_REPEAT = "GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCCCCCCCA"

print("\nX repeat (60 bp):")
print(f"  {X_REPEAT}")
print(f"  Length: {len(X_REPEAT)} bp")

# Find the 7C region
print("\nSearching for C-rich regions:")
for i in range(len(X_REPEAT) - 7):
    window = X_REPEAT[i : i + 10]
    c_count = window.count("C")
    if c_count >= 7:
        print(f"  Position {i:2d}: {window} ({c_count} Cs)")

# The critical region according to protocol
critical_region = "GGGCTCCACCGCCCCCCCA"  # Contains the 7C
print("\nCritical region (7C):")
print(f"  {critical_region}")

if critical_region in X_REPEAT:
    pos = X_REPEAT.index(critical_region)
    print(f"  Found at position {pos} in X repeat")
else:
    print("  ⚠ Not found in X repeat - need to check full repeat context")

# Check MwoI sites in X repeat
x_seq = Seq(X_REPEAT)
x_sites = MwoI.search(x_seq)
print(f"\nMwoI sites in X repeat: {len(x_sites)} at positions {x_sites if x_sites else 'None'}")

# =============================================================================
# TEST 3: Design Primers from X Repeat
# =============================================================================

print("\n[TEST 3] Primer Design from X Repeat")
print("-" * 80)

# According to protocol: "Primers located in 2 contiguous repeats flanking 7C/8C"

# Strategy: Design primers that:
# 1. Forward: Early in repeat (positions 0-20)
# 2. Reverse: Late in repeat (positions 40-60)
# 3. Product should be ~60bp (one repeat) or ~120bp (two repeats)

# Forward primer candidates
print("\nForward primer candidates (early in X repeat):")
for start in range(0, 20, 5):
    fwd = X_REPEAT[start : start + 22]  # 22bp primer
    gc = (fwd.count("G") + fwd.count("C")) / len(fwd) * 100
    print(f"  Pos {start:2d}-{start+22:2d}: {fwd} (GC: {gc:.1f}%)")

# Select forward primer (avoid 7C region which is around pos 40-55)
FORWARD_PRIMER = X_REPEAT[0:22]  # Position 0-22
print(f"\n✓ Selected Forward: {FORWARD_PRIMER}")

# Reverse primer candidates
print("\nReverse primer candidates (late in X repeat):")
for start in range(35, 45, 2):
    rev_rc = X_REPEAT[start : start + 22]  # Take sequence
    rev = str(Seq(rev_rc).reverse_complement())  # Reverse complement
    gc = (rev.count("G") + fwd.count("C")) / len(rev) * 100
    print(f"  Pos {start:2d}-{start+22:2d} (RC): {rev[:30]}... (GC: {gc:.1f}%)")

# Select reverse primer (should cover region after 7C)
REVERSE_PRIMER_RC = X_REPEAT[38:60]  # Position 38-60 (end of repeat)
REVERSE_PRIMER = str(Seq(REVERSE_PRIMER_RC).reverse_complement())
print(f"\n✓ Selected Reverse (RC): {REVERSE_PRIMER}")

print("\nPrimer Summary:")
print(f"  Forward:  5'-{FORWARD_PRIMER}-3'")
print(f"  Reverse:  5'-{REVERSE_PRIMER}-3'")
print("  Expected product: ~60 bp (single repeat)")

# =============================================================================
# TEST 4: Load Real dupC Data
# =============================================================================

print("\n[TEST 4] Load Real dupC Mutation Data")
print("-" * 80)

dupC_mutant_file = Path("output/dupC/dupC.001.mut.simulated.fa")
dupC_normal_file = Path("output/dupC/dupC.001.normal.simulated.fa")

if not dupC_mutant_file.exists():
    print(f"✗ ERROR: {dupC_mutant_file} not found")
    print("  Please run MucOneUp to generate dupC mutation files first")
    sys.exit(1)

# Load sequences
mutant_record = next(iter(SeqIO.parse(dupC_mutant_file, "fasta")))
normal_record = next(iter(SeqIO.parse(dupC_normal_file, "fasta")))

mutant_seq = str(mutant_record.seq)
normal_seq = str(normal_record.seq)

print(f"✓ Loaded mutant:  {len(mutant_seq):,} bp")
print(f"✓ Loaded normal:  {len(normal_seq):,} bp")

# Search for 8C pattern in mutant
pattern_8c = "GCCCCCCCCAGC"  # 8 Cs
pattern_7c = "GCCCCCCCAGC"  # 7 Cs

count_8c_mut = mutant_seq.count(pattern_8c)
count_7c_mut = mutant_seq.count(pattern_7c)
count_8c_norm = normal_seq.count(pattern_8c)
count_7c_norm = normal_seq.count(pattern_7c)

print("\nPattern occurrence:")
print(f"  Mutant - 8C pattern: {count_8c_mut} occurrences")
print(f"  Mutant - 7C pattern: {count_7c_mut} occurrences")
print(f"  Normal - 8C pattern: {count_8c_norm} occurrences")
print(f"  Normal - 7C pattern: {count_7c_norm} occurrences")

if count_8c_mut > count_8c_norm:
    print("  ✓ Mutant has more 8C patterns (dupC mutation confirmed)")
else:
    print("  ⚠ Warning: Expected mutant to have more 8C patterns")

# =============================================================================
# TEST 5: Search for X Repeats in Sequences
# =============================================================================

print("\n[TEST 5] X Repeat Locations in dupC Sequences")
print("-" * 80)


def find_all_occurrences(sequence, pattern):
    """Find all occurrences of pattern in sequence."""
    positions = []
    start = 0
    while True:
        pos = sequence.find(pattern, start)
        if pos == -1:
            break
        positions.append(pos)
        start = pos + 1
    return positions


# Search for exact X repeat
x_positions_mut = find_all_occurrences(mutant_seq, X_REPEAT)
x_positions_norm = find_all_occurrences(normal_seq, X_REPEAT)

print("\nExact X repeat occurrences:")
print(
    f"  Mutant: {len(x_positions_mut)} at positions {x_positions_mut[:5]}{'...' if len(x_positions_mut) > 5 else ''}"
)
print(
    f"  Normal: {len(x_positions_norm)} at positions {x_positions_norm[:5]}{'...' if len(x_positions_norm) > 5 else ''}"
)

# Search for forward primer
fwd_positions_mut = find_all_occurrences(mutant_seq, FORWARD_PRIMER)
fwd_positions_norm = find_all_occurrences(normal_seq, FORWARD_PRIMER)

print("\nForward primer binding sites:")
print(f"  Mutant: {len(fwd_positions_mut)} occurrences")
print(f"  Normal: {len(fwd_positions_norm)} occurrences")

if len(fwd_positions_mut) > 1:
    print("  ⚠ Multiple binding sites! Primers may not be specific enough")
    print(f"     Positions in mutant: {fwd_positions_mut[:10]}")

# Search for reverse primer (reverse complement)
rev_rc = str(Seq(REVERSE_PRIMER).reverse_complement())
rev_positions_mut = find_all_occurrences(mutant_seq, rev_rc)
rev_positions_norm = find_all_occurrences(normal_seq, rev_rc)

print("\nReverse primer binding sites (RC):")
print(f"  Mutant: {len(rev_positions_mut)} occurrences")
print(f"  Normal: {len(rev_positions_norm)} occurrences")

# =============================================================================
# TEST 6: PCR Simulation with pydna
# =============================================================================

print("\n[TEST 6] PCR Simulation with pydna")
print("-" * 80)

try:
    import pydna
    from pydna.amplify import pcr
    from pydna.dseqrecord import Dseqrecord

    print(f"✓ pydna {pydna.__version__} loaded")
except ImportError as e:
    print(f"✗ pydna not available: {e}")
    print("  Install with: pip install pydna")
    sys.exit(1)

print("\nRunning PCR simulation...")
print(f"  Template: dupC mutant ({len(mutant_seq):,} bp)")
print(f"  Forward:  {FORWARD_PRIMER}")
print(f"  Reverse:  {REVERSE_PRIMER}")

template_mut = Dseqrecord(mutant_seq, name="dupC_mutant")

import time

start_time = time.time()
amplicon_mut = pcr(FORWARD_PRIMER, REVERSE_PRIMER, template_mut)
elapsed = time.time() - start_time

if amplicon_mut:
    print("\n✓ PCR SUCCESS")
    print(f"  Product size: {len(amplicon_mut)} bp")
    print(f"  Forward binding: position {amplicon_mut.forward_primer.position}")
    print(f"  Reverse binding: position {amplicon_mut.reverse_primer.position}")
    print(f"  Execution time: {elapsed:.3f} seconds")

    # Check if product contains 8C mutation
    product_seq = str(amplicon_mut.seq)
    if pattern_8c in product_seq:
        print(f"  ✓ Product contains 8C pattern: {product_seq.count(pattern_8c)} occurrences")
    else:
        print("  ⚠ Product does not contain 8C pattern")

    # Show product sequence context
    print("\n  Product sequence (first 100 bp):")
    print(f"  {product_seq[:100]}...")

else:
    print("\n✗ PCR FAILED - No amplification")
    print("  Primers may not bind to template")

# Test with normal sequence
print("\n\nRunning PCR on normal sequence...")
template_norm = Dseqrecord(normal_seq, name="dupC_normal")
amplicon_norm = pcr(FORWARD_PRIMER, REVERSE_PRIMER, template_norm)

if amplicon_norm:
    print("✓ PCR SUCCESS (normal)")
    print(f"  Product size: {len(amplicon_norm)} bp")
    print(f"  Contains 7C: {pattern_7c in str(amplicon_norm.seq)}")
    print(f"  Contains 8C: {pattern_8c in str(amplicon_norm.seq)}")
else:
    print("✗ PCR FAILED (normal)")

# =============================================================================
# TEST 7: MwoI Digest of PCR Products
# =============================================================================

print("\n[TEST 7] MwoI Digest of PCR Products")
print("-" * 80)

if amplicon_mut and amplicon_norm:
    # Test MwoI sites
    mut_product_seq = Seq(str(amplicon_mut.seq))
    norm_product_seq = Seq(str(amplicon_norm.seq))

    mut_sites = MwoI.search(mut_product_seq)
    norm_sites = MwoI.search(norm_product_seq)

    print("\nMwoI sites in PCR products:")
    print(f"  Mutant product: {len(mut_sites)} sites at {mut_sites if mut_sites else 'None'}")
    print(f"  Normal product: {len(norm_sites)} sites at {norm_sites if norm_sites else 'None'}")

    if len(mut_sites) < len(norm_sites):
        print("  ✓ Expected result: Mutant has fewer MwoI sites")
        print("    → Protocol prediction: Mutant should NOT be cut (8C disrupts MwoI)")
        print("    → Protocol prediction: Normal should be cut (7C has MwoI site)")
    elif len(mut_sites) == len(norm_sites):
        print("  ⚠ Same number of sites - need to investigate further")
    else:
        print("  ✗ Unexpected: Mutant has MORE MwoI sites than normal")

    # Simulate digest
    if norm_sites:
        norm_fragments = MwoI.catalyse(norm_product_seq)
        print("\n  Normal digest simulation:")
        print(f"    Fragments: {len(norm_fragments)}")
        print(f"    Fragment sizes: {[len(f) for f in norm_fragments]}")

    if mut_sites:
        mut_fragments = MwoI.catalyse(mut_product_seq)
        print("\n  Mutant digest simulation:")
        print(f"    Fragments: {len(mut_fragments)}")
        print(f"    Fragment sizes: {[len(f) for f in mut_fragments]}")
    else:
        print("\n  Mutant digest simulation:")
        print("    ✓ No MwoI sites - product remains intact (expected for 8C mutation)")

# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)

print("\n✓ Tests Completed:")
print("  [1] MwoI recognition pattern analyzed")
print("  [2] X repeat structure examined")
print("  [3] Primers designed from X repeat")
print("  [4] Real dupC data loaded")
print("  [5] X repeat locations identified")
print("  [6] PCR simulation executed")
print("  [7] MwoI digest simulated")

print("\n⚠ Critical Findings:")
print("  - Primers from X repeat may bind to multiple locations")
print("  - Need 21bp tags to make primers specific (not yet included)")
print("  - Protocol uses REPEAT-based primers, not constant-region primers")
print("  - Complete workflow requires: Digest → PCR → Digest → SNaPshot")

print("\n📋 Next Steps:")
print("  1. Find/design the 21bp primer tags")
print("  2. Test with tagged primers")
print("  3. Design SNaPshot extension primers (Primer 7C and Primer repeat R)")
print("  4. Implement complete workflow simulation")

print("\n" + "=" * 80)
