# muc_one_up/cli.py

import argparse
import json
import sys
import random
from .simulate import simulate_diploid
from .mutate import apply_mutations

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
    # multiple fixed lengths
    parser.add_argument(
        "--fixed-lengths",
        nargs="+",
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
    # NEW: mutation arguments
    parser.add_argument(
        "--mutation-name",
        default=None,
        help="Name of the mutation to apply (must match an entry in config['mutations'])."
    )
    parser.add_argument(
        "--mutation-targets",
        nargs="+",
        default=None,
        help=(
            "One or more 'haplotype_index,repeat_index' pairs, e.g. '1,5' '2,7'. "
            "If provided, each pair indicates the haplotype and repeat number "
            "at which to apply the named mutation once. "
            "If omitted, the mutation (if any) is applied once at a random allowed repeat."
        )
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

    # -------------------------------------------------------------------------
    # MUTATION LOGIC
    # -------------------------------------------------------------------------
    if args.mutation_name:
        # 1) If user provided specific targets => apply once per target
        if args.mutation_targets:
            mutation_positions = []
            for t in args.mutation_targets:
                try:
                    hap_i_str, rep_i_str = t.split(",")
                    hap_i = int(hap_i_str)
                    rep_i = int(rep_i_str)
                    mutation_positions.append((hap_i, rep_i))
                except ValueError:
                    print(f"[ERROR] Invalid --mutation-targets format: '{t}'", file=sys.stderr)
                    sys.exit(1)

            try:
                results = apply_mutations(
                    config=config,
                    results=results,
                    mutation_name=args.mutation_name,
                    targets=mutation_positions
                )
            except Exception as e:
                print(f"[ERROR] Failed to apply mutation '{args.mutation_name}': {e}", file=sys.stderr)
                sys.exit(1)

        # 2) Otherwise, pick a random haplotype & repeat that is allowed for this mutation
        else:
            try:
                # Quick validation
                if "mutations" not in config or args.mutation_name not in config["mutations"]:
                    raise ValueError(f"Mutation '{args.mutation_name}' not in config['mutations']")

                mutation_def = config["mutations"][args.mutation_name]
                allowed_repeats = set(mutation_def["allowed_repeats"])  # for fast 'in' checks

                possible_targets = []  # collect all (hap_idx, repeat_idx) that are allowed
                for hap_idx, (seq, chain) in enumerate(results, start=1):
                    for rep_idx, sym in enumerate(chain, start=1):
                        # strip off 'm' if present in the chain
                        pure_sym = sym.replace("m", "")
                        if pure_sym in allowed_repeats:
                            possible_targets.append((hap_idx, rep_idx))

                if not possible_targets:
                    raise ValueError(
                        f"No repeats match 'allowed_repeats' for mutation '{args.mutation_name}'"
                    )

                # pick one random target
                random_target = random.choice(possible_targets)
                # apply mutation once
                results = apply_mutations(
                    config=config,
                    results=results,
                    mutation_name=args.mutation_name,
                    targets=[random_target]
                )
            except Exception as e:
                print(f"[ERROR] Failed to apply mutation '{args.mutation_name}': {e}", file=sys.stderr)
                sys.exit(1)
    # -------------------------------------------------------------------------
    # end MUTATION LOGIC

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
