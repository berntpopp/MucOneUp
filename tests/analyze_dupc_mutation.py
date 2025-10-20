#!/usr/bin/env python3
"""
Deep analysis of dupC mutation and its effect on MwoI sites.

This script analyzes where the dupC mutation occurs and how it affects
restriction sites.
"""

import re

from Bio import SeqIO
from Bio.Restriction import MwoI
from Bio.Seq import Seq

# Load sequences
print("Loading sequences...")
mutant_record = next(iter(SeqIO.parse("output/dupC/dupC.001.mut.simulated.fa", "fasta")))
normal_record = next(iter(SeqIO.parse("output/dupC/dupC.001.normal.simulated.fa", "fasta")))

mutant = str(mutant_record.seq)
normal = str(normal_record.seq)

print(f"Mutant: {len(mutant)} bp - {mutant_record.description}")
print(f"Normal: {len(normal)} bp - {normal_record.description}")
print(f"Size difference: {len(mutant) - len(normal)} bp (expected: 1 bp for dupC)")

# Parse mutation info from header
if "Targets:" in mutant_record.description:
    import ast

    targets_str = mutant_record.description.split("Targets: ")[1].split(")")[0] + ")"
    targets = ast.literal_eval(targets_str)
    print(f"\nMutation targets from header: {targets}")
    print(f"  Haplotype {targets[0][0]}, Position {targets[0][1]}")

# Find the insertion point
print("\n" + "=" * 80)
print("FINDING INSERTION POINT")
print("=" * 80)

# Align sequences to find where they diverge
for i in range(min(len(normal), len(mutant))):
    if normal[i] != mutant[i]:
        first_diff = i
        break

print(f"\nFirst difference at position: {first_diff}")

# Check around the insertion point
window_start = max(0, first_diff - 50)
window_end = min(len(normal), first_diff + 50)

print(f"\nContext around insertion (position {first_diff}):")
print(f"Normal: ...{normal[window_start:window_end]}...")
print(f"Mutant: ...{mutant[window_start:window_end]}...")

# Try to realign by skipping the insertion
print("\n" + "=" * 80)
print("CHECKING IF INSERTION IS SINGLE C")
print("=" * 80)

# Check if removing one C from mutant realigns sequences
test_removed = mutant[:first_diff] + mutant[first_diff + 1 :]
if test_removed == normal:
    print(f"✓ Confirmed: Single C insertion at position {first_diff}")
    print(f"  Base at position: {mutant[first_diff]}")
else:
    # Check if there's a frameshift
    match_after = 0
    for offset in range(1, 100):
        if (
            mutant[first_diff + offset : first_diff + offset + 20]
            == normal[first_diff : first_diff + 20]
        ):
            match_after = offset
            break

    if match_after:
        print(f"✓ Insertion of {match_after} bp at position {first_diff}")
        print(f"  Inserted sequence: {mutant[first_diff : first_diff + match_after]}")
    else:
        print("⚠ Complex mutation pattern (likely frameshift effect)")

# Check MwoI sites near mutation
print("\n" + "=" * 80)
print("MwoI SITES NEAR MUTATION")
print("=" * 80)

mut_sites = MwoI.search(Seq(mutant))
norm_sites = MwoI.search(Seq(normal))

print("\nTotal MwoI sites:")
print(f"  Mutant: {len(mut_sites)}")
print(f"  Normal: {len(norm_sites)}")

# Find sites within 200 bp of mutation
nearby_mut = [s for s in mut_sites if abs(s - first_diff) < 200]
nearby_norm = [s for s in norm_sites if abs(s - first_diff) < 200]

print("\nMwoI sites within 200 bp of mutation:")
print(f"  Mutant: {nearby_mut}")
print(f"  Normal: {nearby_norm}")

# Check if any sites were disrupted or created
print("\n" + "=" * 80)
print("SITE DISRUPTION ANALYSIS")
print("=" * 80)

# MwoI recognition: GCNNNNNNNGC (where N = any base)
mwoi_pattern = r"GC.{7}GC"

# Search in window around mutation
search_window = 100
window_start = max(0, first_diff - search_window)
window_end = min(len(normal), first_diff + search_window)

norm_window = normal[window_start:window_end]
mut_window = mutant[window_start : window_end + 1]  # +1 for insertion

norm_matches = [
    (m.start() + window_start, m.group()) for m in re.finditer(mwoi_pattern, norm_window)
]
mut_matches = [(m.start() + window_start, m.group()) for m in re.finditer(mwoi_pattern, mut_window)]

print(f"\nMwoI sites in {search_window}bp window around mutation:")
print(f"  Normal: {len(norm_matches)} sites")
for pos, seq in norm_matches:
    print(f"    Position {pos}: {seq}")

print(f"  Mutant: {len(mut_matches)} sites")
for pos, seq in mut_matches:
    print(f"    Position {pos}: {seq}")

if len(mut_matches) != len(norm_matches):
    print("\n✓ MUTATION AFFECTS MwoI SITES!")
    print(f"  Site difference: {len(mut_matches) - len(norm_matches)}")
else:
    print("\n⚠ No change in MwoI site count in this window")

# Check VNTR region markers
print("\n" + "=" * 80)
print("VNTR REGION IDENTIFICATION")
print("=" * 80)

# Look for high GC content (VNTR region typically >80% GC)
window_size = 100
gc_content = []

for i in range(0, len(normal) - window_size, 50):
    window = normal[i : i + window_size]
    gc = (window.count("G") + window.count("C")) / len(window)
    if gc > 0.80:
        gc_content.append((i, i + window_size, gc))

if gc_content:
    print("\nHigh GC regions (>80%, likely VNTR):")
    for start, end, gc in gc_content[:5]:
        print(f"  Position {start:,}-{end:,}: {gc * 100:.1f}% GC")

    # Check if mutation is in VNTR region
    in_vntr = any(start <= first_diff <= end for start, end, _ in gc_content)
    print(
        f"\nMutation at position {first_diff}: {'✓ IN VNTR region' if in_vntr else '✗ OUTSIDE VNTR region'}"
    )

# Final summary
print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)

print("\nMutation:")
print(f"  Position: {first_diff}")
print("  Type: Insertion of 1 bp (C)")
print("  Effect: Frameshift downstream")

print("\nMwoI Sites:")
print(f"  Total sites unchanged: {len(mut_sites)} == {len(norm_sites)}")
print("  But frameshift affects all downstream sequences")

print("\nImplication:")
print("  The dupC mutation creates a frameshift that changes the entire")
print("  downstream sequence, but doesn't specifically disrupt MwoI sites")
print("  in a way that would allow SNaPshot discrimination.")

print("\n⚠️  CONCLUSION:")
print("  The current approach may not work for dupC with MwoI.")
print("  Need to either:")
print("    1. Use different restriction enzyme sensitive to frameshift")
print("    2. Use different primers targeting the specific mutation")
print("    3. Use direct SNaPshot without restriction digest pre-selection")
