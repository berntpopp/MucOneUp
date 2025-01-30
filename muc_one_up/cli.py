# muc_one_up/cli.py

import argparse
import json
import sys
from .simulate import simulate_diploid

def build_parser():
    parser = argparse.ArgumentParser(
        prog="MucOneUp",
        description="Tool to simulate MUC1 VNTR diploid references."
    )
    parser.add_argument(
        "--config",
        required=True,
        help="Path to JSON configuration file."
    )
    parser.add_argument(
        "--output",
        default="muc1_simulated.fa",
        help="Output FASTA filename."
    )
    parser.add_argument(
        "--structure-output",
        default=None,
        help="Optional: Output text file describing the VNTR structure of each haplotype."
    )
    parser.add_argument(
        "--num-haplotypes",
        default=2,
        type=int,
        help="Number of haplotypes to simulate. Typically 2 for diploid."
    )
    # NEW: multiple fixed lengths
    parser.add_argument(
        "--fixed-lengths",
        nargs="+",         # allow multiple arguments
        type=int,
        default=None,
        help=(
            "One or more fixed lengths (in VNTR repeats). "
            "If one value is provided, it will apply to all haplotypes. "
            "If multiple, must match --num-haplotypes."
        )
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducible simulations."
    )
    return parser

def main():
    parser = build_parser()
    args = parser.parse_args()

    # Load config
    try:
        with open(args.config) as fh:
            config = json.load(fh)
    except (FileNotFoundError, json.JSONDecodeError) as e:
        print(f"[ERROR] Could not load config: {e}", file=sys.stderr)
        sys.exit(1)

    # Process --fixed-lengths logic
    fixed_lengths = None
    if args.fixed_lengths is not None:
        if len(args.fixed_lengths) == 1:
            # Replicate one length for all haplotypes
            fixed_lengths = [args.fixed_lengths[0]] * args.num_haplotypes
        elif len(args.fixed_lengths) == args.num_haplotypes:
            # Match each haplotype
            fixed_lengths = args.fixed_lengths
        else:
            print(
                "[ERROR] Number of --fixed-lengths values must be either 1 or exactly equal to --num-haplotypes.",
                file=sys.stderr
            )
            sys.exit(1)

    # Simulate
    try:
        results = simulate_diploid(
            config=config,
            num_haplotypes=args.num_haplotypes,
            fixed_lengths=fixed_lengths,
            seed=args.seed
        )
    except Exception as e:
        print(f"[ERROR] Simulation failed: {e}", file=sys.stderr)
        sys.exit(1)

    # results is a list of (sequence, chain) for each haplotype
    # Write FASTA
    try:
        with open(args.output, "w") as out_fh:
            for i, (sequence, chain) in enumerate(results, start=1):
                out_fh.write(f">haplotype_{i}\n")
                out_fh.write(sequence + "\n")
    except Exception as e:
        print(f"[ERROR] Failed to write FASTA: {e}", file=sys.stderr)
        sys.exit(1)

    # Optionally write VNTR structure
    if args.structure_output:
        try:
            with open(args.structure_output, "w") as struct_fh:
                for i, (sequence, chain) in enumerate(results, start=1):
                    chain_str = "-".join(chain)
                    struct_fh.write(f"haplotype_{i}\t{chain_str}\n")
        except Exception as e:
            print(f"[ERROR] Failed to write VNTR structure file: {e}", file=sys.stderr)
            sys.exit(1)

if __name__ == "__main__":
    main()
