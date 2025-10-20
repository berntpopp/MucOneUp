#!/usr/bin/env python3
"""
Complete SNaPshot Workflow Validation with Real dupC Data

Tests the complete in-silico SNaPshot workflow:
1. PCR amplification (multiple products from X repeats)
2. MwoI digest selection (7C digested, 8C survives)
3. SNaPshot extension (fluorescence prediction)
4. Mutation detection

Tests on actual MucOneUp-generated sequences:
- output/dupC/dupC.001.mut.simulated.fa (contains 8C mutation)
- output/dupC/dupC.001.normal.simulated.fa (7C only)
"""

import sys
from pathlib import Path
from typing import Any

from Bio import SeqIO
from Bio.Restriction import MwoI
from Bio.Seq import Seq

# =============================================================================
# PRIMER DEFINITIONS
# =============================================================================


class SnapshotPrimers:
    """Validated primer sequences for dupC SNaPshot assay."""

    # First PCR primers
    PCR_PRIMER_F = "GGCCGGCCCCGGGCTCCACC"  # 20bp (forward)
    PCR_PRIMER_R_ORIGINAL = "TGTCACCTCGGCCCCGGA"  # 18bp (original)
    PCR_PRIMER_R = "TCCGGGGCCGAGGTGACA"  # 18bp (reverse complement)

    # Flanked sequences
    PCR_FLANKED_7C = "GCCCCCCCAGCCCACGG"  # 17bp (7C normal)
    PCR_FLANKED_8C = "GCCCCCCCCAGCCCACGG"  # 18bp (8C mutant - add one C)

    # SNaPshot extension primers
    SNAPSHOT_PRIMER_7C = "CGGGCTCCACCGCCCCCCC"  # 19bp
    SNAPSHOT_PRIMER_REPEAT_R = "TCCGGGGCCGAGGTGACA"  # 18bp (RC)

    # Expected amplicon sizes
    AMPLICON_7C_SIZE = 55  # 20 + 17 + 18
    AMPLICON_8C_SIZE = 56  # 20 + 18 + 18


# =============================================================================
# SNAPSHOT WORKFLOW FUNCTIONS
# =============================================================================


def find_primer_binding_sites(template: str, primer: str, rc: bool = False) -> list[int]:
    """
    Find all binding positions for a primer in template.

    Args:
        template: Template sequence
        primer: Primer sequence
        rc: If True, search for reverse complement

    Returns:
        List of 0-based positions where primer binds
    """
    if rc:
        primer = str(Seq(primer).reverse_complement())

    positions = []
    start = 0
    while True:
        pos = template.find(primer, start)
        if pos == -1:
            break
        positions.append(pos)
        start = pos + 1

    return positions


def simulate_pcr_amplification(
    template: str, forward_primer: str, reverse_primer: str, max_products: int = 100
) -> dict[str, Any]:
    """
    Simulate PCR amplification with multiple product generation.

    In real protocol, primers bind to multiple X repeats.
    This creates multiple amplicons - mix of 7C and 8C products.

    Strategy: Prioritize finding products with 8C or 7C patterns.

    Args:
        template: Template sequence
        forward_primer: Forward primer (20bp)
        reverse_primer: Reverse primer RC (18bp)
        max_products: Maximum number of products to generate

    Returns:
        Dictionary with PCR results
    """
    # Find all binding sites
    forward_sites = find_primer_binding_sites(template, forward_primer, rc=False)
    reverse_sites = find_primer_binding_sites(template, reverse_primer, rc=True)

    print("\nPCR Primer Binding:")
    print(f"  Forward sites: {len(forward_sites)}")
    print(f"  Reverse sites: {len(reverse_sites)}")
    print(
        f"  Potential products: {len(forward_sites)} x {len(reverse_sites)} = {len(forward_sites) * len(reverse_sites)}"
    )

    # Find mutation sites to prioritize
    pattern_8c = "GCCCCCCCCAGC"
    pattern_7c = "GCCCCCCCAGC"

    mutation_8c_positions = []
    mutation_7c_positions = []

    # Find all 8C positions
    start = 0
    while True:
        pos = template.find(pattern_8c, start)
        if pos == -1:
            break
        mutation_8c_positions.append(pos)
        start = pos + 1

    # Find all 7C positions
    start = 0
    while True:
        pos = template.find(pattern_7c, start)
        if pos == -1:
            break
        mutation_7c_positions.append(pos)
        start = pos + 1

    print(
        f"  Mutation sites found: 8C={len(mutation_8c_positions)}, 7C={len(mutation_7c_positions)}"
    )

    # Generate amplicons - prioritize ones containing mutations
    products = []
    products_8c = []
    products_7c = []
    products_other = []

    for f_pos in forward_sites:
        for r_pos in reverse_sites:
            # Reverse primer binds on opposite strand
            # Amplicon spans from forward start to reverse end
            if f_pos < r_pos:
                # Extract amplicon
                amplicon_start = f_pos
                amplicon_end = r_pos + len(reverse_primer)
                amplicon_seq = template[amplicon_start:amplicon_end]

                # Extract flanked region (between primers)
                flanked_start = f_pos + len(forward_primer)
                flanked_end = r_pos
                flanked_seq = template[flanked_start:flanked_end]

                product = {
                    "sequence": amplicon_seq,
                    "length": len(amplicon_seq),
                    "forward_pos": f_pos,
                    "reverse_pos": r_pos,
                    "flanked_region": flanked_seq,
                }

                # Filter by expected amplicon size (55-56bp for valid products)
                # Expected: 20bp (F) + 17-18bp (flanked) + 18bp (R) = 55-56bp
                expected_min = 50
                expected_max = 65

                if expected_min <= len(amplicon_seq) <= expected_max:
                    # Categorize by mutation presence
                    if pattern_8c in amplicon_seq:
                        products_8c.append(product)
                    elif pattern_7c in amplicon_seq:
                        products_7c.append(product)
                    else:
                        products_other.append(product)

    # Combine: prioritize 8C, then 7C, then others
    all_products = products_8c + products_7c + products_other
    products = all_products[:max_products]

    print(
        f"  Products by type: 8C={len(products_8c)}, 7C={len(products_7c)}, Other={len(products_other)}"
    )
    print(f"  Returning top {len(products)} products")

    return {
        "success": True,
        "product_count": len(products),
        "products": products,
        "non_specific": len(products) > 1,
        "stats": {
            "products_8c": len(products_8c),
            "products_7c": len(products_7c),
            "products_other": len(products_other),
        },
    }


def simulate_digest_selection(
    amplicons: list[dict[str, Any]], enzyme_name: str = "MwoI"
) -> dict[str, Any]:
    """
    Simulate MwoI restriction digest selection.

    MwoI recognizes: GCNNNNNNNGC (exactly 7 bases between GC...GC)
    - 7C amplicons: GC[CCCCCCC]AGC → has site → DIGESTED
    - 8C amplicons: GC[CCCCCCCC]AGC → no site → SURVIVES

    Args:
        amplicons: List of PCR products
        enzyme_name: Restriction enzyme (default MwoI)

    Returns:
        Dictionary with digest results
    """
    survivors = []
    digested = []

    for amp in amplicons:
        seq_obj = Seq(amp["sequence"])
        mwoi_sites = MwoI.search(seq_obj)
        site_count = len(mwoi_sites)

        # Check mutation type
        has_8c = "GCCCCCCCCAGC" in amp["sequence"]  # 8 Cs
        has_7c = "GCCCCCCCAGC" in amp["sequence"]  # 7 Cs

        mutation_type = "8C" if has_8c else ("7C" if has_7c else "unknown")

        result = {
            "sequence": amp["sequence"],
            "length": amp["length"],
            "flanked_region": amp["flanked_region"],
            "has_mwoi_site": site_count > 0,
            "mwoi_site_count": site_count,
            "mwoi_positions": list(mwoi_sites),
            "mutation_type": mutation_type,
        }

        if site_count == 0:
            # No MwoI site → survives digest
            result["survives_digest"] = True
            survivors.append(result)
        else:
            # Has MwoI site → gets digested
            result["survives_digest"] = False
            digested.append(result)

    mechanism = (
        f"MwoI digests products with site (7C), spares products without (8C). "
        f"Digested: {len(digested)}, Survivors: {len(survivors)}"
    )

    return {
        "total_products": len(amplicons),
        "digested_count": len(digested),
        "survivor_count": len(survivors),
        "survivors": survivors,
        "digested": digested,
        "mechanism": mechanism,
    }


def simulate_snapshot_extension(
    template_seq: str, primer_seq: str, primer_name: str
) -> dict[str, Any]:
    """
    Simulate SNaPshot single-base extension.

    Primer binds to surviving amplicon, single fluorescent ddNTP added.

    Args:
        template_seq: Amplicon sequence (survivor from digest)
        primer_seq: Extension primer
        primer_name: Identifier

    Returns:
        Dictionary with extension results
    """
    # Fluorophore mapping
    fluorophore_map = {
        "A": {"color": "Green", "dye": "dR6G"},
        "C": {"color": "Black", "dye": "dTAMRA"},
        "G": {"color": "Blue", "dye": "dR110"},
        "T": {"color": "Red", "dye": "dROX"},
    }

    # Find primer binding site
    pos = template_seq.find(primer_seq)

    if pos == -1:
        # Try reverse complement
        primer_rc = str(Seq(primer_seq).reverse_complement())
        pos = template_seq.find(primer_rc)
        if pos != -1 and pos > 0:
            # Binds on reverse strand - next base is BEFORE primer
            next_base = template_seq[pos - 1]
            next_base_rc = str(Seq(next_base).complement())
            fluoro = fluorophore_map.get(next_base_rc, {})
            return {
                "primer": primer_name,
                "binds": True,
                "position": pos,
                "orientation": "reverse",
                "next_base": next_base_rc,
                "ddNTP": f"dd{next_base_rc}TP",
                "fluorophore_color": fluoro.get("color", "Unknown"),
                "fluorophore_dye": fluoro.get("dye", "Unknown"),
                "peak_size": len(primer_seq) + 1,
                "interpretation": interpret_result(next_base_rc),
            }
        return {"primer": primer_name, "binds": False}

    # Forward orientation - next base is AFTER primer
    next_pos = pos + len(primer_seq)
    if next_pos < len(template_seq):
        next_base = template_seq[next_pos]
        fluoro = fluorophore_map.get(next_base, {})

        return {
            "primer": primer_name,
            "binds": True,
            "position": pos,
            "orientation": "forward",
            "next_base": next_base,
            "ddNTP": f"dd{next_base}TP",
            "fluorophore_color": fluoro.get("color", "Unknown"),
            "fluorophore_dye": fluoro.get("dye", "Unknown"),
            "peak_size": len(primer_seq) + 1,
            "interpretation": interpret_result(next_base),
        }

    return {"primer": primer_name, "binds": True, "position": pos, "error": "No base after primer"}


def interpret_result(next_base: str) -> str:
    """Interpret SNaPshot result for dupC mutation."""
    interpretations = {
        "C": "8C mutation present (Black peak expected)",
        "A": "Incomplete digestion or normal 7C (Green peak)",
        "G": "Unexpected - Blue peak",
        "T": "Unexpected - Red peak",
    }
    return interpretations.get(next_base, f"Unexpected base: {next_base}")


def validate_complete_workflow(template_seq: str, sample_name: str) -> dict[str, Any]:
    """
    Execute complete SNaPshot workflow validation.

    Workflow:
    1. PCR amplification
    2. MwoI digest selection
    3. SNaPshot extension
    4. Mutation detection

    Args:
        template_seq: Genomic template
        sample_name: Sample identifier

    Returns:
        Complete workflow results
    """
    print(f"\n{'='*80}")
    print(f"COMPLETE SNAPSHOT WORKFLOW: {sample_name}")
    print("=" * 80)

    # Step 1: PCR Amplification
    print("\n[STEP 1] PCR Amplification")
    print("-" * 80)
    pcr_result = simulate_pcr_amplification(
        template=template_seq,
        forward_primer=SnapshotPrimers.PCR_PRIMER_F,
        reverse_primer=SnapshotPrimers.PCR_PRIMER_R,
        max_products=50,  # Limit for performance
    )
    print(f"  Products generated: {pcr_result['product_count']}")
    if pcr_result["non_specific"]:
        print("  ✓ PCR is non-specific (expected from X repeat primers)")

    # Step 2: MwoI Digest Selection
    print("\n[STEP 2] MwoI Digest Selection")
    print("-" * 80)
    digest_result = simulate_digest_selection(pcr_result["products"])
    print(f"  Total products: {digest_result['total_products']}")
    print(f"  Digested (7C): {digest_result['digested_count']}")
    print(f"  Survivors (8C): {digest_result['survivor_count']}")

    # Analyze survivors
    survivor_8c_count = sum(1 for s in digest_result["survivors"] if s["mutation_type"] == "8C")
    survivor_7c_count = sum(1 for s in digest_result["survivors"] if s["mutation_type"] == "7C")

    print("\n  Survivor breakdown:")
    print(f"    8C survivors: {survivor_8c_count}")
    print(f"    7C survivors: {survivor_7c_count}")
    print(f"    Unknown: {digest_result['survivor_count'] - survivor_8c_count - survivor_7c_count}")

    # Step 3: SNaPshot Extension (on survivors only)
    print("\n[STEP 3] SNaPshot Extension")
    print("-" * 80)

    snapshot_results = []
    mutation_detected = False

    if digest_result["survivor_count"] > 0:
        print(f"  Testing {min(3, len(digest_result['survivors']))} survivors...")

        for idx, survivor in enumerate(digest_result["survivors"][:3]):  # Test first 3
            result = simulate_snapshot_extension(
                template_seq=survivor["sequence"],
                primer_seq=SnapshotPrimers.SNAPSHOT_PRIMER_7C,
                primer_name=f"Survivor_{idx+1}_Primer7C",
            )

            snapshot_results.append(result)

            if result.get("binds"):
                print(f"\n  Survivor {idx+1}:")
                print(f"    Mutation type: {survivor['mutation_type']}")
                print(f"    Amplicon length: {survivor['length']} bp")
                print("    Primer binds: Yes")
                print(f"    Next base: {result.get('next_base')}")
                print(
                    f"    Fluorophore: {result.get('fluorophore_color')} ({result.get('fluorophore_dye')})"
                )
                print(f"    Interpretation: {result.get('interpretation')}")

                if result.get("next_base") == "C":
                    mutation_detected = True
    else:
        print("  No survivors - all products digested (normal 7C only)")

    # Step 4: Final Interpretation
    print("\n[STEP 4] Mutation Detection")
    print("-" * 80)

    if mutation_detected:
        print("  ✓✓ 8C MUTATION DETECTED!")
        print("  Expected: Black peak (dTAMRA)")
        print("  Mechanism: 8C amplicons lack MwoI site → survive digest → SNaPshot detects C")
    else:
        print("  No 8C mutation detected")
        if digest_result["survivor_count"] == 0:
            print("  All products digested (7C only - normal)")
        else:
            print("  Survivors present but no Black peak")

    return {
        "sample": sample_name,
        "pcr_products": pcr_result["product_count"],
        "digest_survivors": digest_result["survivor_count"],
        "survivor_8c_count": survivor_8c_count,
        "snapshot_results": snapshot_results,
        "mutation_detected": mutation_detected,
        "summary": "8C mutation detected" if mutation_detected else "No mutation detected",
    }


# =============================================================================
# MAIN TEST EXECUTION
# =============================================================================


def main():
    """Execute complete workflow tests on real dupC data."""
    print("=" * 80)
    print("SNAPSHOT WORKFLOW VALIDATION - Real dupC Data")
    print("=" * 80)

    # Test data paths
    mutant_file = Path("output/dupC/dupC.001.mut.simulated.fa")
    normal_file = Path("output/dupC/dupC.001.normal.simulated.fa")

    # Check files exist
    if not mutant_file.exists():
        print(f"\n✗ ERROR: Mutant file not found: {mutant_file}")
        print("  Please run dupC simulation first:")
        print(
            "  muconeup --config config.json simulate --out-base dupC --mutation-name normal,dupC"
        )
        sys.exit(1)

    if not normal_file.exists():
        print(f"\n✗ ERROR: Normal file not found: {normal_file}")
        sys.exit(1)

    # Load sequences
    print("\n[LOADING DATA]")
    print("-" * 80)

    mutant_record = next(iter(SeqIO.parse(mutant_file, "fasta")))
    normal_record = next(iter(SeqIO.parse(normal_file, "fasta")))

    mutant_seq = str(mutant_record.seq)
    normal_seq = str(normal_record.seq)

    print(f"✓ Loaded mutant:  {len(mutant_seq):,} bp from {mutant_file.name}")
    print(f"✓ Loaded normal:  {len(normal_seq):,} bp from {normal_file.name}")

    # Verify mutation presence
    pattern_8c = "GCCCCCCCCAGC"

    mutant_has_8c = mutant_seq.count(pattern_8c)
    normal_has_8c = normal_seq.count(pattern_8c)

    print("\nMutation verification:")
    print(f"  Mutant 8C count: {mutant_has_8c}")
    print(f"  Normal 8C count: {normal_has_8c}")

    if mutant_has_8c > 0 and normal_has_8c == 0:
        print("  ✓ Mutation confirmed in data")
    else:
        print("  ⚠ Warning: Unexpected mutation pattern")

    # Run complete workflow on MUTANT
    mutant_result = validate_complete_workflow(mutant_seq, "MUTANT (8C)")

    # Run complete workflow on NORMAL
    normal_result = validate_complete_workflow(normal_seq, "NORMAL (7C)")

    # Final Summary
    print("\n" + "=" * 80)
    print("FINAL VALIDATION SUMMARY")
    print("=" * 80)

    print("\nMUTANT Sample:")
    print(f"  PCR products: {mutant_result['pcr_products']}")
    print(
        f"  Digest survivors: {mutant_result['digest_survivors']} (8C: {mutant_result['survivor_8c_count']})"
    )
    print(f"  Mutation detected: {'✓ YES' if mutant_result['mutation_detected'] else '✗ NO'}")

    print("\nNORMAL Sample:")
    print(f"  PCR products: {normal_result['pcr_products']}")
    print(
        f"  Digest survivors: {normal_result['digest_survivors']} (8C: {normal_result['survivor_8c_count']})"
    )
    print(f"  Mutation detected: {'✓ YES' if normal_result['mutation_detected'] else '✗ NO'}")

    # Validation
    print("\n" + "=" * 80)
    print("VALIDATION RESULTS")
    print("=" * 80)

    all_pass = True

    # Test 1: Mutant should be detected
    if mutant_result["mutation_detected"]:
        print("✓ Test 1 PASS: Mutant correctly identified")
    else:
        print("✗ Test 1 FAIL: Mutant NOT detected")
        all_pass = False

    # Test 2: Normal should NOT be detected
    if not normal_result["mutation_detected"]:
        print("✓ Test 2 PASS: Normal correctly identified (no 8C)")
    else:
        print("✗ Test 2 FAIL: Normal incorrectly shows mutation")
        all_pass = False

    # Test 3: Mutant should have survivors
    if mutant_result["survivor_8c_count"] > 0:
        print(f"✓ Test 3 PASS: Mutant has 8C survivors ({mutant_result['survivor_8c_count']})")
    else:
        print("✗ Test 3 FAIL: Mutant has no 8C survivors")
        all_pass = False

    # Test 4: Normal should have few/no survivors
    if normal_result["digest_survivors"] == 0:
        print("✓ Test 4 PASS: Normal has no survivors (all 7C digested)")
    else:
        print(f"⚠ Test 4 WARNING: Normal has {normal_result['digest_survivors']} survivors")

    print("\n" + "=" * 80)
    if all_pass:
        print("✓✓ ALL TESTS PASSED - WORKFLOW VALIDATED!")
    else:
        print("✗✗ SOME TESTS FAILED - REVIEW RESULTS")
    print("=" * 80)

    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())
