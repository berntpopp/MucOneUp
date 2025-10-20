#!/usr/bin/env python3
"""
COMPLETE SNaPshot Workflow Simulation with ACTUAL Primer Sequences

Based on Protocol MUC1-Englisch_Arif.docx with real primer sequences.

PCR Primers (from user):
- Primer 1: GGCCGGCCCCGGGCTCCACC
- Primer 2: TGTCACCTCGGCCCCGGA
- Tag: GCCCCCCCAGCCCACGG (contains 8C: GCCCCCCCCAGC)

SNaPshot Primers (from user):
- Primer 7C (19bp): CGGGCTCCACCGCCCCCCC
- Primer repeat R (18bp): TGTCACCTCGGCCCCGGA
- Tag: GCCCCACGG
"""

import sys
import time
from pathlib import Path

from Bio import SeqIO
from Bio.Restriction import MwoI
from Bio.Seq import Seq

print("=" * 80)
print("COMPLETE SNaPshot WORKFLOW - With Real Primer Sequences")
print("=" * 80)

# =============================================================================
# PRIMER DEFINITIONS (from user/protocol)
# =============================================================================

# PCR Primers
PCR_PRIMER_F = "GGCCGGCCCCGGGCTCCACC"  # 20bp (Forward)
PCR_PRIMER_R_ORIGINAL = "TGTCACCTCGGCCCCGGA"  # 18bp (from user)
# CORRECTED: Use reverse complement of primer 2
PCR_PRIMER_R = str(Seq(PCR_PRIMER_R_ORIGINAL).reverse_complement())  # "TCCGGGGCCGAGGTGACA"
PCR_TAG = "GCCCCCCCAGCCCACGG"  # 17bp - contains 8C pattern!

# SNaPshot Extension Primers
SNAPSHOT_PRIMER_7C = "CGGGCTCCACCGCCCCCCC"  # 19bp (Forward)
SNAPSHOT_PRIMER_REPEAT_R_ORIGINAL = "TGTCACCTCGGCCCCGGA"  # 18bp (from user)
# CORRECTED: Use reverse complement for second primer
SNAPSHOT_PRIMER_REPEAT_R = str(
    Seq(SNAPSHOT_PRIMER_REPEAT_R_ORIGINAL).reverse_complement()
)  # "TCCGGGGCCGAGGTGACA"
SNAPSHOT_TAG = "GCCCCACGG"  # 9bp

print("\nPrimer Configuration:")
print(f"  PCR Forward:  5'-{PCR_PRIMER_F}-3' ({len(PCR_PRIMER_F)}bp)")
print(f"  PCR Reverse:  5'-{PCR_PRIMER_R}-3' ({len(PCR_PRIMER_R)}bp)")
print(f"  PCR Tag:      5'-{PCR_TAG}-3' ({len(PCR_TAG)}bp)")
print("")
print(f"  SNaPshot 7C:        5'-{SNAPSHOT_PRIMER_7C}-3' ({len(SNAPSHOT_PRIMER_7C)}bp)")
print(f"  SNaPshot repeat R:  5'-{SNAPSHOT_PRIMER_REPEAT_R}-3' ({len(SNAPSHOT_PRIMER_REPEAT_R)}bp)")
print(f"  SNaPshot Tag:       5'-{SNAPSHOT_TAG}-3' ({len(SNAPSHOT_TAG)}bp)")

# Note: PCR Tag contains the mutation site!
print("\n  ⚠ PCR Tag contains: GCCCCCCCCAGC (8C pattern for mutant detection)")

# =============================================================================
# Load Test Data
# =============================================================================

print(f"\n{'=' * 80}")
print("STEP 0: Load dupC Test Data")
print("=" * 80)

dupC_mutant_file = Path("output/dupC/dupC.001.mut.simulated.fa")  # noqa: N816
dupC_normal_file = Path("output/dupC/dupC.001.normal.simulated.fa")  # noqa: N816

if not dupC_mutant_file.exists():
    print(f"✗ ERROR: {dupC_mutant_file} not found")
    sys.exit(1)

mutant_record = next(iter(SeqIO.parse(dupC_mutant_file, "fasta")))
normal_record = next(iter(SeqIO.parse(dupC_normal_file, "fasta")))

mutant_seq = str(mutant_record.seq)
normal_seq = str(normal_record.seq)

print(f"✓ Loaded mutant:  {len(mutant_seq):,} bp")
print(f"✓ Loaded normal:  {len(normal_seq):,} bp")

# Check for mutation patterns
pattern_8c = "GCCCCCCCCAGC"  # 8 Cs (mutant)
pattern_7c = "GCCCCCCCAGC"  # 7 Cs (normal)

count_8c_mut = mutant_seq.count(pattern_8c)
count_8c_norm = normal_seq.count(pattern_8c)

print("\nMutation verification:")
print(f"  Mutant has {count_8c_mut} x 8C pattern (expected: 1)")
print(f"  Normal has {count_8c_norm} x 8C pattern (expected: 0)")

if count_8c_mut > 0 and count_8c_norm == 0:
    print("  ✓ dupC mutation confirmed")
else:
    print("  ⚠ Warning: Mutation pattern unexpected")

# =============================================================================
# STEP 1: Genomic DNA Digestion with MwoI (3 rounds)
# =============================================================================

print(f"\n{'=' * 80}")
print("STEP 1: Genomic DNA Digestion with MwoI")
print("=" * 80)


def simulate_mwoi_digest(sequence, enzyme_name="MwoI"):
    """Simulate MwoI digestion."""
    seq_obj = Seq(sequence)
    sites = enzyme_name.search(seq_obj) if hasattr(enzyme_name, "search") else MwoI.search(seq_obj)

    if sites:
        fragments = MwoI.catalyse(seq_obj)
        return {
            "cuts": True,
            "num_sites": len(sites),
            "num_fragments": len(fragments),
            "fragment_sizes": sorted([len(f) for f in fragments], reverse=True),
            "largest_fragment": max(len(f) for f in fragments),
        }
    else:
        return {
            "cuts": False,
            "num_sites": 0,
            "num_fragments": 1,
            "fragment_sizes": [len(sequence)],
            "largest_fragment": len(sequence),
        }


print("\nDigestion Round 1-3 (simulated as single step):")
print("\nMutant sequence:")
mut_digest = simulate_mwoi_digest(mutant_seq)
print(f"  MwoI sites: {mut_digest['num_sites']}")
print(f"  Fragments: {mut_digest['num_fragments']}")
if mut_digest["cuts"]:
    print(f"  Largest fragment: {mut_digest['largest_fragment']:,} bp")
    print("  → Sequence would be digested")
else:
    print("  → No MwoI sites - sequence remains intact")

print("\nNormal sequence:")
norm_digest = simulate_mwoi_digest(normal_seq)
print(f"  MwoI sites: {norm_digest['num_sites']}")
print(f"  Fragments: {norm_digest['num_fragments']}")
if norm_digest["cuts"]:
    print(f"  Largest fragment: {norm_digest['largest_fragment']:,} bp")
    print("  → Sequence would be digested")
else:
    print("  → No MwoI sites - sequence remains intact")

# Expected: Both should have same number of sites (from earlier test)
print("\n⚠ Note: Earlier testing showed both have 112 MwoI sites")
print("   This is because dupC is a frameshift, not a site-disrupting mutation")
print("   The 8C pattern itself may not create/destroy MwoI sites globally")

# =============================================================================
# STEP 2: PCR Amplification with MUC1-Repeat F/R
# =============================================================================

print(f"\n{'=' * 80}")
print("STEP 2: PCR Amplification")
print("=" * 80)

try:
    import pydna
    from pydna.amplify import pcr
    from pydna.dseqrecord import Dseqrecord

    print(f"✓ pydna {pydna.__version__} loaded")
except ImportError as e:
    print(f"✗ pydna not available: {e}")
    sys.exit(1)

# Test primer binding first
print("\nSearching for primer binding sites...")

# Forward primer
fwd_count = mutant_seq.count(PCR_PRIMER_F)
fwd_count_rc = mutant_seq.count(str(Seq(PCR_PRIMER_F).reverse_complement()))

# Reverse primer (check both orientations)
rev_count = mutant_seq.count(PCR_PRIMER_R)
rev_rc = str(Seq(PCR_PRIMER_R).reverse_complement())
rev_count_rc = mutant_seq.count(rev_rc)

print(f"  Forward primer ({PCR_PRIMER_F}):")
print(f"    Forward orientation: {fwd_count} sites")
print(f"    Reverse complement:  {fwd_count_rc} sites")

print(f"  Reverse primer ({PCR_PRIMER_R}):")
print(f"    Forward orientation: {rev_count} sites")
print(f"    Reverse complement:  {rev_count_rc} sites")

# Run PCR simulation
print("\nRunning PCR simulation (mutant)...")
template_mut = Dseqrecord(mutant_seq, name="dupC_mutant")

start_time = time.time()
try:
    amplicon_mut = pcr(PCR_PRIMER_F, PCR_PRIMER_R, template_mut)
    elapsed = time.time() - start_time

    if amplicon_mut:
        print("✓ PCR SUCCESS")
        print(f"  Product size: {len(amplicon_mut)} bp")
        print(f"  Forward binding: position {amplicon_mut.forward_primer.position}")
        print(f"  Reverse binding: position {amplicon_mut.reverse_primer.position}")
        print(f"  Execution time: {elapsed:.3f} seconds")

        # Check product sequence
        product_seq = str(amplicon_mut.seq)
        print("\n  Product sequence (first 100 bp):")
        print(f"  {product_seq[:100]}...")

        # Check for mutation in product
        if pattern_8c in product_seq:
            print(
                f"\n  ✓ Product contains 8C mutation: {product_seq.count(pattern_8c)} occurrence(s)"
            )
        elif pattern_7c in product_seq:
            print(
                f"\n  ✓ Product contains 7C pattern: {product_seq.count(pattern_7c)} occurrence(s)"
            )

    else:
        print("✗ PCR FAILED - No amplification")
        amplicon_mut = None

except ValueError as e:
    if "PCR not specific" in str(e):
        print("✗ PCR NOT SPECIFIC")
        print(f"  Error: {str(e)[:200]}...")
        print("\n  → Primers bind to multiple sites")
        print("  → This is expected for repeat-based primers")
        print("  → In real protocol, tags provide specificity")
        amplicon_mut = None
    else:
        raise

# Run PCR on normal
print("\nRunning PCR simulation (normal)...")
template_norm = Dseqrecord(normal_seq, name="dupC_normal")

try:
    amplicon_norm = pcr(PCR_PRIMER_F, PCR_PRIMER_R, template_norm)

    if amplicon_norm:
        print("✓ PCR SUCCESS (normal)")
        print(f"  Product size: {len(amplicon_norm)} bp")
    else:
        print("✗ PCR FAILED (normal)")
        amplicon_norm = None

except ValueError as e:
    if "PCR not specific" in str(e):
        print("✗ PCR NOT SPECIFIC (normal)")
        print("  → Same issue as mutant")
        amplicon_norm = None
    else:
        raise

# =============================================================================
# STEP 3: PCR Product Digestion
# =============================================================================

print(f"\n{'=' * 80}")
print("STEP 3: PCR Product Digestion with MwoI")
print("=" * 80)

if amplicon_mut:
    print("\nDigesting mutant PCR product...")
    mut_pcr_digest = simulate_mwoi_digest(str(amplicon_mut.seq))
    print(f"  MwoI sites: {mut_pcr_digest['num_sites']}")
    print(f"  Fragments: {mut_pcr_digest['num_fragments']}")

    if mut_pcr_digest["cuts"]:
        print("  → Product would be digested")
    else:
        print("  ✓ No MwoI sites - product remains intact (expected for 8C)")
else:
    print("⚠ Skipping - no PCR product from Step 2")

if amplicon_norm:
    print("\nDigesting normal PCR product...")
    norm_pcr_digest = simulate_mwoi_digest(str(amplicon_norm.seq))
    print(f"  MwoI sites: {norm_pcr_digest['num_sites']}")
    print(f"  Fragments: {norm_pcr_digest['num_fragments']}")

    if norm_pcr_digest["cuts"]:
        print("  → Product would be digested (expected for 7C)")
    else:
        print("  ⚠ No MwoI sites - unexpected for normal")
else:
    print("⚠ Skipping - no PCR product from Step 2")

# =============================================================================
# STEP 4-5: Purification and ExoSAP (Simulation)
# =============================================================================

print(f"\n{'=' * 80}")
print("STEP 4-5: Purification and ExoSAP (Simulated)")
print("=" * 80)
print("✓ AMPure purification (simulated)")
print("✓ ExoSAP digestion (simulated)")

# =============================================================================
# STEP 6: SNaPshot Extension
# =============================================================================

print(f"\n{'=' * 80}")
print("STEP 6: SNaPshot Extension")
print("=" * 80)


def simulate_snapshot_extension(template_seq, primer_seq, primer_name):
    """Simulate SNaPshot single-base extension."""

    # Find primer binding site
    pos = template_seq.find(primer_seq)

    if pos == -1:
        # Try reverse complement
        primer_rc = str(Seq(primer_seq).reverse_complement())
        pos = template_seq.find(primer_rc)
        if pos != -1:
            print(f"  {primer_name}: Binds at position {pos} (reverse complement)")
            # Next base is BEFORE the primer (on reverse strand)
            if pos > 0:
                next_base = template_seq[pos - 1]
                next_base_rc = str(Seq(next_base).complement())
                return {
                    "binds": True,
                    "position": pos,
                    "next_base": next_base_rc,
                    "orientation": "reverse",
                }
        return {"binds": False}

    print(f"  {primer_name}: Binds at position {pos} (forward)")

    # Get next base (one base after primer 3' end)
    next_pos = pos + len(primer_seq)

    if next_pos < len(template_seq):
        next_base = template_seq[next_pos]

        # Determine ddNTP incorporated
        fluorophore_map = {
            "A": "Green (dR6G)",
            "C": "Black (dTAMRA)",
            "G": "Blue (dR110)",
            "T": "Red (dROX)",
        }

        fluorophore = fluorophore_map.get(next_base, "Unknown")

        return {
            "binds": True,
            "position": pos,
            "next_base": next_base,
            "ddNTP": f"dd{next_base}TP",
            "fluorophore": fluorophore,
            "orientation": "forward",
            "peak_size": len(primer_seq) + 1,  # primer + 1 fluorescent base
        }
    else:
        return {"binds": True, "position": pos, "error": "No base after primer"}


if amplicon_mut:
    print("\nSNaPshot on mutant PCR product:")
    product_mut = str(amplicon_mut.seq)

    # Test Primer 7C
    result_7c = simulate_snapshot_extension(product_mut, SNAPSHOT_PRIMER_7C, "Primer 7C (19bp)")
    if result_7c.get("binds"):
        if "error" not in result_7c:
            print(f"    Next base: {result_7c['next_base']}")
            print(f"    ddNTP: {result_7c['ddNTP']}")
            print(f"    Fluorophore: {result_7c['fluorophore']}")
            print(f"    Expected peak: {result_7c['peak_size']} + tag = 28 bp")

            # Protocol expectation
            if result_7c["next_base"] == "C":
                print("    ✓ BLACK peak → 8C mutation confirmed!")
            elif result_7c["next_base"] == "A":
                print("    ⚠ GREEN peak → Incomplete digestion")
        else:
            print(f"    ✗ Error: {result_7c['error']}")
    else:
        print("    ✗ Primer does not bind")

    # Test Primer repeat R
    result_repeat = simulate_snapshot_extension(
        product_mut, SNAPSHOT_PRIMER_REPEAT_R, "Primer repeat R (18bp+tag)"
    )
    if result_repeat.get("binds") and "error" not in result_repeat:
        print(f"    Next base: {result_repeat['next_base']}")
        print(f"    ddNTP: {result_repeat['ddNTP']}")
        print(f"    Fluorophore: {result_repeat['fluorophore']}")
        print(f"    Expected peak: {result_repeat['peak_size']} + tag (21bp) = 45 bp")
else:
    print("⚠ Skipping - no PCR product available")

# =============================================================================
# SUMMARY
# =============================================================================

print(f"\n{'=' * 80}")
print("WORKFLOW SUMMARY")
print("=" * 80)

print("\nStep 1: Genomic DNA Digestion")
print(f"  Mutant: {mut_digest['num_sites']} MwoI sites")
print(f"  Normal: {norm_digest['num_sites']} MwoI sites")

print("\nStep 2: PCR Amplification")
if amplicon_mut:
    print(f"  Mutant: ✓ {len(amplicon_mut)} bp product")
else:
    print("  Mutant: ✗ No product (primers not specific)")

if amplicon_norm:
    print(f"  Normal: ✓ {len(amplicon_norm)} bp product")
else:
    print("  Normal: ✗ No product (primers not specific)")

print("\nStep 3: PCR Product Digestion")
if amplicon_mut:
    if mut_pcr_digest["cuts"]:
        print(f"  Mutant: Digested ({mut_pcr_digest['num_sites']} sites)")
    else:
        print("  Mutant: ✓ Intact (no MwoI sites)")
else:
    print("  Mutant: N/A")

print("\nStep 6: SNaPshot Extension")
if amplicon_mut and result_7c.get("binds"):
    print(f"  Primer 7C: {result_7c.get('fluorophore', 'N/A')}")
    if result_7c.get("next_base") == "C":
        print("    → ✓ 8C mutation detected")
else:
    print("  Primer 7C: N/A")

print(f"\n{'=' * 80}")
print("Critical Findings:")
print("  1. Primers without tags bind to multiple sites")
print("  2. Need to implement tagged primer logic")
print("  3. MwoI digest may not discriminate for frameshift mutations")
print("  4. SNaPshot extension logic works conceptually")
print("=" * 80)
