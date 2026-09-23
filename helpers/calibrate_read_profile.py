#!/usr/bin/env python3
"""Calibrate a read profile against real-data aggregate targets and validate it.

Two engines:

* ``--engine empirical`` (recommended for ONT R10): the profile's ``errors``
  section is derived directly from context-split error rates (non-homopolymer
  mismatch/insertion/deletion per exposed base, indel length pmfs, per-read error
  spread) and the stutter table equals the measured homopolymer spectrum,
  because the empirical channel never touches homopolymer runs.
* ``--engine pbsim3``: pbsim3 adds its own context-free homopolymer errors, so
  the injected table is *deconvolved* (Richardson-Lucy against a stutter-free
  pbsim3 run, then multiplicative fixed-point corrections).

Either way the profile is then simulated through MucOneUp's own truth-tracked
path and re-measured; residuals are stored in ``provenance.calibration_report``.

Inputs are aggregate statistics only (a targets JSON with ``hp_P_obs_given_true``
and ``error_rates_per_ref_base``). Requires pbsim3 on PATH (or --pbsim3) and the
``edlib`` Python package (``pip install edlib``).

Example:
    python helpers/calibrate_read_profile.py --targets targets.json \\
        --target-key ont_amplicon_PRJEB92208 --template amplicon.fa \\
        --base-profile base.json --model-file reference/pbsim3/QSHMM-ONT-HQ.model \\
        --out muc_one_up/data/read_profiles/ont_r10_sup_amplicon_v1.json
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
    StutterTable,
    build_amplicon_molecules,
    homopolymer_runs,
    reverse_complement,
)

DELTAS = range(-3, 4)
MAX_LEN = 8


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
    fastq: Path, truth: Path, template: str, max_reads: int
) -> tuple[dict[str, Counter], float]:
    """Observed homopolymer deltas per key and the total edit rate, from truth-tracked reads."""
    strands = {}
    with gzip.open(truth, "rt") as handle:
        next(handle)
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields[3] == "full":
                strands[fields[0]] = fields[4]
    runs = [(s, e, b) for s, e, b in homopolymer_runs(template, 2) if e - s <= MAX_LEN]
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
    target: dict[int, float], kernel: dict[int, float], iters: int = 300
) -> dict[int, float]:
    q = {d: 1.0 / len(DELTAS) for d in DELTAS}
    for _ in range(iters):
        model = _convolve(q, kernel)
        new = {}
        for d in DELTAS:
            new[d] = q[d] * sum(
                pk * target.get(d + e, 0.0) / max(model.get(d + e, 0.0), 1e-12)
                for e, pk in kernel.items()
            )
        total = sum(new.values()) or 1.0
        q = {d: v / total for d, v in new.items()}
    return q


def simulate(
    template: str, n: int, model: MoleculeModel, run: Sequencer, seed: int
) -> tuple[dict[str, Counter], float]:
    molecules = build_amplicon_molecules([template], [n], model, random.Random(seed))
    with tempfile.TemporaryDirectory(prefix="calibrate_") as tmp:
        work = Path(tmp)
        fastq, truth = work / "r.fastq", work / "t.tsv.gz"
        simulate_molecule_reads(molecules, run, work / "sim", fastq, truth, "cal", seed)
        return measure(fastq, truth, template, n)


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


def main() -> None:
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
    ap.add_argument("--model-file", required=True)
    ap.add_argument("--model-type", default="qshmm")
    ap.add_argument("--pbsim3", default="pbsim")
    ap.add_argument("--n-molecules", type=int, default=600)
    ap.add_argument("--min-target-n", type=int, default=2000)
    ap.add_argument("--rounds", type=int, default=2)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    targets = json.loads(args.targets.read_text())[args.target_key]
    hp_targets = {
        k: {int(d): p for d, p in v["p_obs_minus_true"].items()}
        for k, v in targets["hp_P_obs_given_true"].items()
        if not k.endswith("|both") and v["n"] >= args.min_target_n and int(k[1:].split("|")[0]) >= 3
    }
    profile = json.loads(args.base_profile.read_text())
    ont = profile.get("config_overrides", {}).get("ont_amplicon_params", {})
    run = PbsimRun(
        args.pbsim3,
        "samtools",
        args.model_type,
        args.model_file,
        accuracy_mean=ont.get("accuracy_mean", 0.95),
        difference_ratio=ont.get("difference_ratio"),
    )
    template = str(next(SeqIO.parse(args.template, "fasta")).seq).upper()
    base_model = MoleculeModel(forward_frac=0.5)
    report: dict = {
        "engine": args.engine,
        "target_total_error": targets["error_rates_per_ref_base"]["all:total"]["median"],
    }
    sequencer: Sequencer = run
    if args.engine == "empirical":
        profile["errors"] = empirical_errors_from_targets(targets, read_error_sigma(targets))
        sequencer = EmpiricalSequencer(EmpiricalErrorModel.from_dict(profile["errors"]))
        stutter = {k: {d: t.get(d, 0.0) for d in DELTAS} for k, t in hp_targets.items()}
        for q in stutter.values():  # restrict to DELTAS and renormalise
            total = sum(q.values())
            for d in q:
                q[d] /= total
        rounds = 1
    else:
        pbsim_spectra, pbsim_err = simulate(template, args.n_molecules, base_model, run, args.seed)
        report["pbsim_only_total_error"] = round(pbsim_err, 5)
        stutter = {
            key: richardson_lucy(target, _pmf(pbsim_spectra.get(key, Counter({0: 1}))))
            for key, target in hp_targets.items()
        }
        rounds = args.rounds
    for rnd in range(1, rounds + 1):
        table = StutterTable.from_dict(
            {k: {str(d): p for d, p in q.items()} for k, q in stutter.items()}
        )
        model = MoleculeModel(forward_frac=0.5, stutter=table)
        observed, total_err = simulate(
            template, args.n_molecules, model, sequencer, args.seed + rnd
        )
        residual = {}
        for key, target in hp_targets.items():
            obs = _pmf(observed.get(key, Counter({0: 1})))
            residual[key] = round(max(abs(obs.get(d, 0) - target.get(d, 0)) for d in DELTAS), 4)
            if args.engine == "pbsim3":
                corrected = {
                    d: stutter[key][d] * (target.get(d, 0) + 1e-4) / (obs.get(d, 0) + 1e-4)
                    for d in DELTAS
                }
                total = sum(corrected.values())
                stutter[key] = {d: v / total for d, v in corrected.items()}
        report[f"validation_round{rnd}"] = {
            "total_error": round(total_err, 5),
            "max_abs_residual_per_key": residual,
        }

    profile.setdefault("molecules", {})["stutter"] = {
        k: {str(d): round(p, 5) for d, p in q.items() if p >= 1e-5}
        for k, q in sorted(stutter.items())
    }
    for key, pmf in profile["molecules"]["stutter"].items():  # renormalise after rounding
        total = sum(pmf.values())
        profile["molecules"]["stutter"][key] = {d: round(p / total, 5) for d, p in pmf.items()}
        drift = round(1.0 - sum(profile["molecules"]["stutter"][key].values()), 5)
        profile["molecules"]["stutter"][key]["0"] = round(
            profile["molecules"]["stutter"][key].get("0", 0) + drift, 5
        )
    profile.setdefault("provenance", {})["calibration_report"] = report
    args.out.write_text(json.dumps(profile, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
