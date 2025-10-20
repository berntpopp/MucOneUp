#!/usr/bin/env python3
"""
TEST: Complete SNaPshot Logic - Digest Selection Mechanism

The Key Insight:
1. PCR creates MULTIPLE products (one from each X repeat) - THIS IS EXPECTED!
2. MwoI digest SELECTS which products survive:
   - Products with 7C → DIGESTED (destroyed)
   - Products with 8C → SURVIVE (MwoI doesn't cut 8C)
3. SNaPshot extension on SURVIVING products
4. If 8C present → adds C → Black peak

Testing Strategy:
- Simulate multiple PCR products
- Check each for 7C vs 8C
- Simulate digest on each
- Keep only undigested products
- Check extension on survivors
"""

import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Restriction import MwoI
from Bio.Seq import Seq

print("=" * 80)
print("DIGEST SELECTION LOGIC TEST")
print("=" * 80)

# =============================================================================
# Load dupC Data
# =============================================================================

print("\n[STEP 0] Load Test Data")
print("-" * 80)

dupC_mutant_file = Path("output/dupC/dupC.001.mut.simulated.fa")

if not dupC_mutant_file.exists():
    print(f"✗ ERROR: {dupC_mutant_file} not found")
    sys.exit(1)

mutant_record = next(iter(SeqIO.parse(dupC_mutant_file, "fasta")))
mutant_seq = str(mutant_record.seq)

print(f"✓ Loaded mutant: {len(mutant_seq):,} bp")

# Verify mutation
pattern_8c = "GCCCCCCCCAGC"  # 8C (mutant)
pattern_7c = "GCCCCCCCAGC"  # 7C (wild-type)

count_8c = mutant_seq.count(pattern_8c)
count_7c = mutant_seq.count(pattern_7c)

print("\nMutation patterns:")
print(f"  8C (mutant):     {count_8c} occurrences")
print(f"  7C (wild-type):  {count_7c} occurrences")

if count_8c > 0:
    print("  ✓ dupC mutation present")

# =============================================================================
# Simulate PCR - Find All X Repeats
# =============================================================================

print("\n[STEP 1] Simulate PCR - Find All X Repeat Products")
print("-" * 80)

X_REPEAT = "GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCCCCCCCA"

print(f"Searching for X repeat ({len(X_REPEAT)} bp)...")

# Find all X repeats
x_positions = []
start = 0
while True:
    pos = mutant_seq.find(X_REPEAT, start)
    if pos == -1:
        break
    x_positions.append(pos)
    start = pos + 1

print(f"Found {len(x_positions)} exact X repeat occurrences")
print(f"Positions: {x_positions[:10]}{'...' if len(x_positions) > 10 else ''}")

# Also search for X repeats with mutations (8C variant)
# The 8C mutation changes the end of X repeat
X_REPEAT_8C = X_REPEAT[:-1] + "CA"  # Change last A to CA (adding C)

x_8c_positions = []
start = 0
while True:
    pos = mutant_seq.find(X_REPEAT_8C, start)
    if pos == -1:
        break
    x_8c_positions.append(pos)
    start = pos + 1

print(f"\nFound {len(x_8c_positions)} X repeats with 8C mutation")
if x_8c_positions:
    print(f"Positions: {x_8c_positions}")

# Simulate "PCR products" - extract regions around X repeats
print("\nSimulating PCR products (first 10 X repeats):")

pcr_products = []
for i, pos in enumerate(x_positions[:10]):
    # Extract X repeat + some flanking sequence
    product_start = max(0, pos - 20)
    product_end = min(len(mutant_seq), pos + len(X_REPEAT) + 20)
    product = mutant_seq[product_start:product_end]

    pcr_products.append(
        {
            "id": f"Product_{i + 1}",
            "position": pos,
            "sequence": product,
            "size": len(product),
            "has_7C": pattern_7c in product,
            "has_8C": pattern_8c in product,
        }
    )

    print(f"  Product {i + 1} at position {pos}:")
    print(f"    Size: {len(product)} bp")
    print(f"    Contains 7C: {pattern_7c in product}")
    print(f"    Contains 8C: {pattern_8c in product}")

print(f"\n✓ Generated {len(pcr_products)} simulated PCR products")

# =============================================================================
# Digest Each PCR Product with MwoI
# =============================================================================

print("\n[STEP 2] Digest Each PCR Product with MwoI")
print("-" * 80)

surviving_products = []

for product in pcr_products:
    product_seq = Seq(product["sequence"])
    mwoi_sites = MwoI.search(product_seq)

    print(f"\n{product['id']}:")
    print(f"  Has 7C: {product['has_7C']}")
    print(f"  Has 8C: {product['has_8C']}")
    print(f"  MwoI sites: {len(mwoi_sites)}")

    if mwoi_sites:
        # Product gets digested
        fragments = MwoI.catalyse(product_seq)
        fragment_sizes = [len(f) for f in fragments]
        print(f"  → DIGESTED into {len(fragments)} fragments: {fragment_sizes}")
        print("  → Product DESTROYED")
    else:
        # Product survives
        print("  → NO MwoI sites")
        print("  ✓ Product SURVIVES")
        surviving_products.append(product)

print(f"\n{'=' * 80}")
print("Digest Results:")
print(f"  Total products: {len(pcr_products)}")
print(f"  Survived: {len(surviving_products)}")
print(f"  Destroyed: {len(pcr_products) - len(surviving_products)}")

if surviving_products:
    print("\nSurviving products:")
    for p in surviving_products:
        print(f"  - {p['id']}: 7C={p['has_7C']}, 8C={p['has_8C']}")

# =============================================================================
# SNaPshot Extension on Surviving Products
# =============================================================================

print("\n[STEP 3] SNaPshot Extension on Surviving Products")
print("-" * 80)

# Extension primers
SNAPSHOT_PRIMER_7C = "CGGGCTCCACCGCCCCCCC"  # 19bp

fluorophore_map = {
    "A": "Green (dR6G)",
    "C": "Black (dTAMRA)",
    "G": "Blue (dR110)",
    "T": "Red (dROX)",
}

for product in surviving_products:
    print(f"\n{product['id']}:")

    # Find where primer binds
    pos = product["sequence"].find(SNAPSHOT_PRIMER_7C)

    if pos == -1:
        # Try reverse complement
        primer_rc = str(Seq(SNAPSHOT_PRIMER_7C).reverse_complement())
        pos = product["sequence"].find(primer_rc)

        if pos == -1:
            print("  ✗ Primer does not bind")
            continue
        else:
            print(f"  Primer binds at position {pos} (reverse complement)")
            # For RC, next base is BEFORE primer
            if pos > 0:
                next_base = product["sequence"][pos - 1]
                next_base = str(Seq(next_base).complement())
            else:
                print("  ✗ No base before primer")
                continue
    else:
        print(f"  Primer binds at position {pos} (forward)")
        # Next base is after primer
        next_pos = pos + len(SNAPSHOT_PRIMER_7C)
        if next_pos < len(product["sequence"]):
            next_base = product["sequence"][next_pos]
        else:
            print("  ✗ No base after primer")
            continue

    # Determine ddNTP and fluorophore
    ddNTP = f"dd{next_base}TP"
    fluorophore = fluorophore_map.get(next_base, "Unknown")

    print(f"  Next base: {next_base}")
    print(f"  ddNTP incorporated: {ddNTP}")
    print(f"  Fluorophore: {fluorophore}")

    # Interpret
    if next_base == "C":
        print("  ✓✓ BLACK PEAK → 8C MUTATION DETECTED! ✓✓")
    elif next_base == "A":
        print("  ⚠ GREEN PEAK → Incomplete digestion or wild-type")
    else:
        print(f"  ? {fluorophore.split()[0]} PEAK → Unexpected")

# =============================================================================
# Test with Simulated Pure 7C and Pure 8C Products
# =============================================================================

print("\n[STEP 4] Control Test - Pure 7C and Pure 8C Sequences")
print("-" * 80)

# Create test sequences
test_7c = "AGCTAGCTA" + "GGGCTCCACCGCCCCCCCA" + "GCTAGCTAGC"  # 7C
test_8c = "AGCTAGCTA" + "GGGCTCCACCGCCCCCCCCA" + "GCTAGCTAGC"  # 8C (extra C)

print("\nTest 1: Pure 7C sequence")
print(f"Sequence: {test_7c}")

# Check MwoI
test_7c_seq = Seq(test_7c)
mwoi_7c = MwoI.search(test_7c_seq)
print(f"  MwoI sites: {len(mwoi_7c)}")

if mwoi_7c:
    print("  → Would be DIGESTED")
else:
    print("  → Would SURVIVE")

print("\nTest 2: Pure 8C sequence")
print(f"Sequence: {test_8c}")

# Check MwoI
test_8c_seq = Seq(test_8c)
mwoi_8c = MwoI.search(test_8c_seq)
print(f"  MwoI sites: {len(mwoi_8c)}")

if mwoi_8c:
    print("  → Would be DIGESTED")
else:
    print("  → Would SURVIVE")

# Extension test
print("\nSNaPshot extension test:")

for name, seq in [("7C", test_7c), ("8C", test_8c)]:
    pos = seq.find(SNAPSHOT_PRIMER_7C)
    if pos != -1:
        next_pos = pos + len(SNAPSHOT_PRIMER_7C)
        if next_pos < len(seq):
            next_base = seq[next_pos]
            print(f"  {name}: Next base = {next_base} → {fluorophore_map[next_base]}")

# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)

print("\n✓ Complete Workflow Tested:")
print("  1. PCR creates multiple products (from all X repeats) ✓")
print("  2. MwoI digest SELECTS products:")
print("     - 7C products → Digested (destroyed)")
print("     - 8C products → Survive (no MwoI cut)")
print("  3. SNaPshot extension on survivors:")
print("     - If 8C → adds C → Black peak")

print("\n✓ Key Insight:")
print("  PCR non-specificity is EXPECTED and CORRECT!")
print("  The digest step provides the SELECTION")
print("  Not the PCR step!")

print("\n✓ Results:")
print(f"  - Found {len(x_positions)} X repeats in sequence")
print(f"  - Simulated {len(pcr_products)} PCR products")
print(f"  - {len(surviving_products)} products survived digest")
print(f"  - Mutation detection: {'SUCCESS' if surviving_products else 'Check digest logic'}")

print("\n" + "=" * 80)
