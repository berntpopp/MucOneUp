#!/usr/bin/env python3
"""Calibrate a read profile against real-data aggregate targets and validate it.

Two engines:

* ``--engine empirical`` (recommended for ONT R10): the profile's ``errors``
  section is derived directly from context-split error rates (non-homopolymer
  mismatch/insertion/deletion per exposed base, indel length pmfs, per-read error
  spread) and the stutter table equals the measured homopolymer spectrum,
  because the empirical channel never touches homopolymer runs. Every per-strand
  key with ``n >= --min-target-n`` is fitted (``|delta| <= --max-delta``), and
  ``molecules.stutter_fallback`` is derived for runs without a fitted entry
  (see ``muc_one_up/read_simulator/stutter_fit.py``); the fit is recorded in
  ``provenance.stutter_fit``.
* ``--engine pbsim3``: pbsim3 adds its own context-free homopolymer errors, so
  the injected table is *deconvolved* (Richardson-Lucy against a stutter-free
  pbsim3 run, then multiplicative fixed-point corrections). The output profile
  uses pbsim3 as the sequencer, so an ``errors`` model, ``stutter_fallback``
  and ``provenance.stutter_fit`` inherited from the base profile are dropped.

Either way the profile is then simulated through MucOneUp's own truth-tracked
path and re-measured. ``provenance.calibration_report`` stores, per
base/length/strand of the targets, the simulated vs target observed-minus-true
pmf (``per_key``: p_correct, max absolute residual, pass/fail against
``--tolerance``, or ``not_in_template``). Use a template that contains the runs
to validate (for MUC1, a dupC haplotype for C8).

Inputs are aggregate statistics only (a targets JSON with ``hp_P_obs_given_true``
and ``error_rates_per_ref_base``). ``--engine pbsim3`` requires pbsim3 on PATH
(or --pbsim3) and ``--model-file``; both engines need the ``edlib`` Python
package (``pip install edlib``).

Example:
    python helpers/calibrate_read_profile.py --targets targets.json \\
        --target-key ont_amplicon_PRJEB92208 --template amplicon_dupc.fa \\
        --base-profile muc_one_up/data/read_profiles/ont_r10_sup_amplicon_v1.json \\
        --n-molecules 4000 --out muc_one_up/data/read_profiles/ont_r10_sup_amplicon_v1.json
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import random
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

import edlib
from Bio import SeqIO

from muc_one_up.read_simulator.empirical_errors import EmpiricalErrorModel
from muc_one_up.read_simulator.molecule_pipeline import (
    EmpiricalSequencer,
    PbsimRun,
    Sequencer,
    simulate_molecule_reads,
)
from muc_one_up.read_simulator.molecules import (
    MoleculeModel,
    StutterFallback,
    StutterTable,
    build_amplicon_molecules,
    homopolymer_runs,
    reverse_complement,
)
from muc_one_up.read_simulator.stutter import MIN_RUN_LEN
from muc_one_up.read_simulator.stutter_fit import (
    compare_to_targets,
    fit_stutter_pmfs,
    round_pmf,
    stutter_profile_section,
)


def _t2q(read: str, template: str) -> list[int]:
    """Template position -> read position map from a global alignment."""
    cigar = edlib.align(read, template, mode="NW", task="path")["cigar"]
    t2q, ti, qi, num = [0] * (len(template) + 1), 0, 0, ""
    for ch in cigar:
        if ch.isdigit():
            num += ch
            continue
        n = int(num)
        num = ""
        if ch in "=X":
            for _ in range(n):
                t2q[ti] = qi
                ti += 1
                qi += 1
        elif ch == "I":
            qi += n
        else:
            for _ in range(n):
                t2q[ti] = qi
                ti += 1
    t2q[len(template)] = qi
    return t2q


def _observed_run(read: str, t2q: list[int], start: int, end: int, base: str) -> int:
    a = t2q[start]
    while a > 0 and read[a - 1] == base:
        a -= 1
    b = max(t2q[end], t2q[start])
    while b < len(read) and read[b] == base:
        b += 1
    return max(
        (len(m) for m in "".join(c if c == base else " " for c in read[a:b]).split()), default=0
    )


def measure(
    fastq: Path, truth: Path, template: str, max_reads: int, max_len: int
) -> tuple[dict[str, Counter], float]:
    """Observed homopolymer deltas per key and the total edit rate, from truth-tracked reads."""
    strands = {}
    with gzip.open(truth, "rt") as handle:
        next(handle)
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields[3] == "full":
                strands[fields[0]] = fields[4]
    runs = [(s, e, b) for s, e, b in homopolymer_runs(template, 2) if e - s <= max_len]
    spectra: dict[str, Counter] = defaultdict(Counter)
    edits = bases = 0
    lines = fastq.read_text().splitlines()
    for i in range(0, min(len(lines), 4 * max_reads), 4):
        name, seq = lines[i][1:], lines[i + 1]
        if name not in strands:
            continue
        strand = strands[name]
        forward = seq if strand == "+" else reverse_complement(seq)
        edits += edlib.align(forward, template, mode="NW")["editDistance"]
        bases += len(template)
        t2q = _t2q(forward, template)
        for start, end, base in runs:
            delta = _observed_run(forward, t2q, start, end, base) - (end - start)
            spectra[f"{base}{end - start}|{strand}"][delta] += 1
    return spectra, edits / max(bases, 1)


def _pmf(counter: Counter) -> dict[int, float]:
    total = sum(counter.values()) or 1
    return {d: counter.get(d, 0) / total for d in range(-6, 7)}


def _convolve(q: dict[int, float], kernel: dict[int, float]) -> dict[int, float]:
    out: dict[int, float] = defaultdict(float)
    for d, pq in q.items():
        for e, pk in kernel.items():
            out[d + e] += pq * pk
    return out


def richardson_lucy(
    target: dict[int, float], kernel: dict[int, float], deltas: range, iters: int = 300
) -> dict[int, float]:
    q = {d: 1.0 / len(deltas) for d in deltas}
    for _ in range(iters):
        model = _convolve(q, kernel)
        new = {}
        for d in deltas:
            new[d] = q[d] * sum(
                pk * target.get(d + e, 0.0) / max(model.get(d + e, 0.0), 1e-12)
                for e, pk in kernel.items()
            )
        total = sum(new.values()) or 1.0
        q = {d: v / total for d, v in new.items()}
    return q


def simulate(
    template: str, n: int, model: MoleculeModel, run: Sequencer, seed: int, max_len: int
) -> tuple[dict[str, Counter], float]:
    molecules = build_amplicon_molecules([template], [n], model, random.Random(seed))
    with tempfile.TemporaryDirectory(prefix="calibrate_") as tmp:
        work = Path(tmp)
        fastq, truth = work / "r.fastq", work / "t.tsv.gz"
        simulate_molecule_reads(molecules, run, work / "sim", fastq, truth, "cal", seed)
        return measure(fastq, truth, template, n, max_len)


def empirical_errors_from_targets(targets: dict, sigma: float) -> dict:
    """Per-exposed-base rates (median over libraries) from context-split measurements."""
    rates: dict[str, list[float]] = defaultdict(list)
    for lib in targets["error_context_per_ref_base"].values():
        exposed = 1.0 - lib["hp_base_frac"]
        for key, name in (
            ("X_nonhp", "mismatch_rate"),
            ("I_nonhp", "insertion_rate"),
            ("D_nonhp", "deletion_rate"),
        ):
            rates[name].append(lib[key] / exposed)
    median = {k: sorted(v)[len(v) // 2] for k, v in rates.items()}
    pmfs = {
        f"{kind}_len_pmf": {k: round(p / sum(v.values()), 5) for k, p in v.items() if p > 0}
        for kind, v in (
            ("insertion", targets["indel_len_pmf"]["ins"]),
            ("deletion", targets["indel_len_pmf"]["del"]),
        )
    }
    for pmf in pmfs.values():  # absorb rounding drift into length 1
        pmf["1"] = round(pmf["1"] + 1.0 - sum(pmf.values()), 5)
    return {**{k: round(v, 6) for k, v in median.items()}, **pmfs, "read_error_sigma": sigma}


def read_error_sigma(targets: dict) -> float:
    """Log-normal sigma of per-read error from the median library's p5/p95 spread."""
    spreads = sorted(
        math.log(q["p95"] / q["p5"]) / (2 * 1.645)
        for q in targets["per_read_err_quantiles"].values()
    )
    return round(spreads[len(spreads) // 2], 3)


def _key_deltas(key: str, deltas: range) -> range:
    """Deltas allowed for a key: a run cannot lose more bases than it has."""
    return range(max(deltas.start, -int(key.partition("|")[0][1:])), deltas.stop)


def _parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--targets", type=Path, required=True)
    ap.add_argument("--target-key", required=True)
    ap.add_argument(
        "--template", type=Path, required=True, help="FASTA; first record, forward orientation"
    )
    ap.add_argument(
        "--base-profile", type=Path, required=True, help="profile JSON to fill with stutter"
    )
    ap.add_argument("--engine", choices=["empirical", "pbsim3"], default="empirical")
    ap.add_argument("--model-file", help="pbsim3 model (required for --engine pbsim3)")
    ap.add_argument("--model-type", default="qshmm")
    ap.add_argument("--pbsim3", default="pbsim")
    ap.add_argument("--n-molecules", type=int, default=600)
    ap.add_argument(
        "--min-target-n",
        type=int,
        default=500,
        help="fit a base/length/strand key only with at least this many target observations",
    )
    ap.add_argument(
        "--max-delta", type=int, default=6, help="largest |observed - true| length kept"
    )
    ap.add_argument(
        "--min-lengths-for-slope",
        type=int,
        default=3,
        help="fitted lengths a base/strand needs to contribute to the extrapolation slope",
    )
    ap.add_argument(
        "--generic-ref-len",
        type=int,
        default=3,
        help="run length whose pooled fitted pmfs form the generic fallback pmf",
    )
    ap.add_argument(
        "--max-extrapolation-bases",
        type=int,
        default=3,
        help="cap on how many bases beyond the fitted range the log-odds slope is extrapolated",
    )
    ap.add_argument(
        "--tolerance",
        type=float,
        default=0.03,
        help="per-key pass threshold on the max absolute pmf residual",
    )
    ap.add_argument("--rounds", type=int, default=2)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    if args.engine == "pbsim3" and not args.model_file:
        ap.error("--engine pbsim3 requires --model-file")
    return args


def main() -> None:
    args = _parse_args()
    targets = json.loads(args.targets.read_text())[args.target_key]
    hp_all = targets["hp_P_obs_given_true"]
    max_len = max(int(k.partition("|")[0][1:]) for k in hp_all)
    deltas = range(-args.max_delta, args.max_delta + 1)
    profile = json.loads(args.base_profile.read_text())
    template = str(next(SeqIO.parse(args.template, "fasta")).seq).upper()
    report: dict = {
        "engine": args.engine,
        "target_total_error": targets["error_rates_per_ref_base"]["all:total"]["median"],
    }
    molecules = profile.setdefault("molecules", {})
    provenance = profile.setdefault("provenance", {})
    sequencer: Sequencer
    if args.engine == "empirical":
        derived = empirical_errors_from_targets(targets, read_error_sigma(targets))
        profile["errors"] = {**profile.get("errors", {}), **derived}
        error_model = EmpiricalErrorModel.from_dict(profile["errors"])
        if error_model.hp_min_len != MIN_RUN_LEN:
            # Runs >= hp_min_len are protected, runs >= MIN_RUN_LEN are stuttered;
            # a mismatch would leave runs protected but unstuttered or doubly noised.
            raise SystemExit(
                f"errors.hp_min_len ({error_model.hp_min_len}) must equal the stutter "
                f"table's minimum run length ({MIN_RUN_LEN}) for the empirical engine"
            )
        sequencer = EmpiricalSequencer(error_model)
        section, provenance["stutter_fit"] = stutter_profile_section(
            hp_all,
            min_n=args.min_target_n,
            min_len=error_model.hp_min_len,
            max_delta=args.max_delta,
            min_lengths_for_slope=args.min_lengths_for_slope,
            generic_ref_len=args.generic_ref_len,
            max_extrapolation_bases=args.max_extrapolation_bases,
        )
        provenance["stutter_fit"]["target_key"] = args.target_key
        molecules.update(section)
        provenance.pop("uncovered_runs", None)  # superseded by stutter_fit (#132)
        stutter = {k: {int(d): p for d, p in v.items()} for k, v in section["stutter"].items()}
        fallback = StutterFallback.from_dict(section["stutter_fallback"])
        rounds = 1
    else:
        ont = profile.get("config_overrides", {}).get("ont_amplicon_params", {})
        sequencer = PbsimRun(
            args.pbsim3,
            "samtools",
            args.model_type,
            args.model_file,
            accuracy_mean=ont.get("accuracy_mean", 0.95),
            difference_ratio=ont.get("difference_ratio"),
        )
        hp_targets = fit_stutter_pmfs(
            hp_all,
            min_n=args.min_target_n,
            min_len=MIN_RUN_LEN,
            max_delta=args.max_delta,
        )
        pbsim_spectra, pbsim_err = simulate(
            template,
            args.n_molecules,
            MoleculeModel(forward_frac=0.5),
            sequencer,
            args.seed,
            max_len,
        )
        report["pbsim_only_total_error"] = round(pbsim_err, 5)
        stutter = {
            key: richardson_lucy(
                target, _pmf(pbsim_spectra.get(key, Counter({0: 1}))), _key_deltas(key, deltas)
            )
            for key, target in hp_targets.items()
        }
        fallback = None
        rounds = args.rounds
    for rnd in range(1, rounds + 1):
        table = StutterTable.from_dict(
            {k: {str(d): p for d, p in q.items()} for k, q in stutter.items()}, fallback=fallback
        )
        model = MoleculeModel(forward_frac=0.5, stutter=table)
        observed, total_err = simulate(
            template, args.n_molecules, model, sequencer, args.seed + rnd, max_len
        )
        per_key = compare_to_targets(
            observed, hp_all, max_delta=args.max_delta, tolerance=args.tolerance
        )
        compared = [row for row in per_key.values() if row["status"] == "compared"]
        report[f"validation_round{rnd}"] = {
            "total_error": round(total_err, 5),
            "template": args.template.name,
            "n_molecules": args.n_molecules,
            "tolerance": args.tolerance,
            "keys_within_tolerance": f"{sum(r['within_tolerance'] for r in compared)}/{len(compared)}",
            "per_key": per_key,
        }
        if args.engine == "pbsim3":
            for key, target in hp_targets.items():
                obs = _pmf(observed.get(key, Counter({0: 1})))
                corrected = {
                    d: stutter[key][d] * (target.get(d, 0) + 1e-4) / (obs.get(d, 0) + 1e-4)
                    for d in _key_deltas(key, deltas)
                }
                total = sum(corrected.values())
                stutter[key] = {d: v / total for d, v in corrected.items()}
    if args.engine == "pbsim3":
        molecules["stutter"] = {k: round_pmf(q) for k, q in sorted(stutter.items())}
        # The profile is calibrated for pbsim3 as the sequencer. An empirical
        # errors model, its stutter_fallback and fit record inherited from the
        # base profile describe a different table: drop them. pbsim3 adds its own
        # homopolymer errors, so unfitted runs are not error-free on this path.
        molecules.pop("stutter_fallback", None)
        profile.pop("errors", None)
        provenance.pop("stutter_fit", None)
    provenance["calibration_report"] = report
    args.out.write_text(json.dumps(profile, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
