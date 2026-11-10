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
  ``provenance.stutter_fit``. With ``--fit-strand-rates`` the channel's
  mismatch/insertion/deletion rates are then fitted per read strand
  (``errors.strand_rates``) so the simulated per-strand totals, stutter
  included, match ``error_rates_per_ref_base["<strand>:<type>"]`` (#131).
  With ``--fit-unit-stutter`` (either engine) ``molecules.unit_stutter_rate`` and
  ``unit_stutter_pmf`` (whole-repeat-unit PCR slippage) are fitted from the
  per-allele ``per_allele_stutter_and_length_noise`` ratios and recorded in
  ``provenance.unit_stutter_fit`` (#131). With ``--fit-unit-stutter ROUNDS --vntr START END``
  they are then refitted by simulation, so that reads of the template measured
  with the real profiler's anchor-span rule reproduce those ratios
  (``provenance.unit_stutter_fit.refit``). ``strand_fit.runs_per_60bp`` records
  the template's homopolymer runs against the target's: per-strand rates fitted
  on a template with fewer long runs than real alleles come out too high.
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

_TYPES = {"X": "mismatch", "I": "ins", "D": "del"}


def _t2q(cigar: str, template: str, types: Counter | None = None) -> list[int]:
    """Template position -> read position map from a global alignment's CIGAR."""
    t2q, ti, qi, num = [0] * (len(template) + 1), 0, 0, ""
    for ch in cigar:
        if ch.isdigit():
            num += ch
            continue
        n = int(num)
        num = ""
        if types is not None and ch in _TYPES:
            types[_TYPES[ch]] += n
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


def _strands(truth: Path, kinds: tuple[str, ...]) -> dict[str, str]:
    """Read id -> strand of the truth rows of the given kinds."""
    strands = {}
    with gzip.open(truth, "rt") as handle:
        next(handle)
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields[3] in kinds:
                strands[fields[0]] = fields[4]
    return strands


def measure(
    fastq: Path, truth: Path, template: str, max_reads: int, max_len: int
) -> tuple[dict[str, Counter], float, dict[str, float]]:
    """Homopolymer deltas per key, total edit rate and ``"<strand>:<type>"`` rates per base."""
    strands = _strands(truth, ("full",))
    runs = [(s, e, b) for s, e, b in homopolymer_runs(template, 2) if e - s <= max_len]
    spectra: dict[str, Counter] = defaultdict(Counter)
    edits = bases = 0
    types: dict[str, Counter] = {"+": Counter(), "-": Counter()}
    lines = fastq.read_text().splitlines()
    for i in range(0, min(len(lines), 4 * max_reads), 4):
        name, seq = lines[i][1:], lines[i + 1]
        if name not in strands:
            continue
        strand = strands[name]
        forward = seq if strand == "+" else reverse_complement(seq)
        aln = edlib.align(forward, template, mode="NW", task="path")
        edits += aln["editDistance"]
        bases += len(template)
        types[strand]["bases"] += len(template)
        t2q = _t2q(aln["cigar"], template, types[strand])
        for start, end, base in runs:
            delta = _observed_run(forward, t2q, start, end, base) - (end - start)
            spectra[f"{base}{end - start}|{strand}"][delta] += 1
    rates = {
        f"{strand}:{kind}": c[kind] / max(c["bases"], 1)
        for strand, c in types.items()
        for kind in _TYPES.values()
    }
    return spectra, edits / max(bases, 1), rates


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


def _simulated(template, n, model, run, seed, measure_fn, unit_bounds=()):
    """Simulate ``n`` molecules of ``template`` and return ``measure_fn(fastq, truth)``."""
    molecules = build_amplicon_molecules(
        [template],
        [n],
        model,
        random.Random(seed),
        unit_bounds=[list(unit_bounds)] if unit_bounds else (),
    )
    with tempfile.TemporaryDirectory(prefix="calibrate_") as tmp:
        work = Path(tmp)
        fastq, truth = work / "r.fastq", work / "t.tsv.gz"
        simulate_molecule_reads(molecules, run, work / "sim", fastq, truth, "cal", seed)
        return measure_fn(fastq, truth)


def simulate(
    template: str, n: int, model: MoleculeModel, run: Sequencer, seed: int, max_len: int
) -> tuple[dict[str, Counter], float, dict[str, float]]:
    return _simulated(
        template, n, model, run, seed, lambda f, t: measure(f, t, template, n, max_len)
    )


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


def fit_strand_rates(
    errors: dict,
    targets: dict,
    template: str,
    model: MoleculeModel,
    n: int,
    seed: int,
    max_len: int,
    rounds: int,
) -> dict:
    """Fit ``errors.strand_rates`` so simulated per-strand type rates hit the targets.

    Stutter adds homopolymer indels on top of the channel, and the channel never
    touches protected runs, so each round moves a channel rate by the remaining
    per-base gap divided by the exposed (non-homopolymer) fraction of the template.
    """
    min_len = EmpiricalErrorModel.from_dict(errors).hp_min_len
    exposed = 1.0 - sum(e - s for s, e, _ in homopolymer_runs(template, min_len)) / len(template)
    target = {k: v["median"] for k, v in targets["error_rates_per_ref_base"].items()}
    names = {"mismatch": "mismatch_rate", "ins": "insertion_rate", "del": "deletion_rate"}
    rates = {s: {name: errors[name] for name in names.values()} for s in ("+", "-")}
    history = []
    for rnd in range(rounds + 1):
        run = EmpiricalSequencer(EmpiricalErrorModel.from_dict({**errors, "strand_rates": rates}))
        _, _, sim = simulate(template, n, model, run, seed + 100 + rnd, max_len)
        history.append({k: round(v, 5) for k, v in sorted(sim.items())})
        if rnd == rounds:
            break
        for key, value in sim.items():
            strand, kind = key.split(":")
            name = names[kind]
            step = (target[key] - value) / exposed
            rates[strand][name] = round(min(0.2, max(0.0, rates[strand][name] + step)), 6)
    report = {
        "exposed_frac": round(exposed, 4),
        "target": {k: target[k] for k in sorted(history[-1])},
        "simulated_per_round": history,
    }
    return {"strand_rates": rates, "report": report}


UNIT_STUTTER_KEYS = {-2: "minus2_per_onpeak", -1: "minus1_per_onpeak", 1: "plus1_per_onpeak"}


def fit_unit_stutter(targets: dict) -> tuple[dict, dict]:
    """Whole-unit PCR stutter (``molecules`` keys and a fit record) from per-allele ratios.

    A full product slips with probability ``rate * units``. Per allele the target
    gives slipped products per on-peak read ``r_k`` for k = -2, -1, +1 units, so
    ``n_on_peak * sum(r_k)`` products slipped out of ``n_on_peak * (1 + sum(r_k))``
    molecules. ``rate`` pools them over alleles per molecule unit; the pmf is the
    pooled share of each k.
    """
    rows = targets["per_allele_stutter_and_length_noise"]
    mass = {k: sum(r["n_on_peak"] * r[key] for r in rows) for k, key in UNIT_STUTTER_KEYS.items()}
    slipped = sum(mass.values())
    exposure = sum(r["n_on_peak"] * r["units"] for r in rows) + sum(
        r["n_on_peak"] * r[key] * r["units"] for r in rows for key in UNIT_STUTTER_KEYS.values()
    )
    pmf = {str(k): round(v / slipped, 5) for k, v in mass.items()}
    pmf["-1"] = round(pmf["-1"] + 1.0 - sum(pmf.values()), 5)  # absorb rounding drift
    section = {"unit_stutter_rate": round(slipped / exposure, 6), "unit_stutter_pmf": pmf}
    record = {
        "method": "P(slip) of a full product = unit_stutter_rate * VNTR units; rate = "
        "pooled slipped products / pooled molecule units, pmf = pooled share of each k, "
        "over per_allele_stutter_and_length_noise (k = -2, -1, +1 units)",
        "n_alleles": len(rows),
        "n_libraries": len({r["lib"] for r in rows}),
        "n_on_peak": sum(r["n_on_peak"] for r in rows),
        "units_range": [min(r["units"] for r in rows), max(r["units"] for r in rows)],
    }
    return section, record


UNIT = 60  # MUC1 repeat unit (bp); the real profiler's anchors are one unit long
ANCHOR_MAX_EDITS, ANCHOR_MAX_HITS = 12, 6  # real profiler anchor search (edlib infix)


def _anchor_hits(anchor: str, read: str) -> list[tuple[int, int]]:
    """Sorted non-overlapping ``[start, end)`` hits of an anchor, as the real profiler finds them."""
    hits = []
    for _ in range(ANCHOR_MAX_HITS):
        found = edlib.align(anchor, read, mode="HW", task="locations", k=ANCHOR_MAX_EDITS)
        if found["editDistance"] < 0:
            break
        s, e = found["locations"][0]
        hits.append((s, e + 1))
        read = read[:s] + "N" * (e + 1 - s) + read[e + 1 :]
    return sorted(hits)


def unit_offsets(fastq: Path, truth: Path, template: str, vntr: tuple[int, int]) -> dict:
    """The real profiler's per-allele stutter statistic on simulated reads of one template.

    Span = first left anchor (the VNTR's first unit) to the first right anchor (its
    last unit) after it; ``dist`` = span - VNTR length. On-peak: ``|dist| <=
    max(30, 0.012 * length)``; a read is k units off if ``round(dist / UNIT) == k``.
    """
    lo, hi = vntr
    left, right, allele = template[lo : lo + UNIT], template[hi - UNIT : hi], hi - lo
    strands = _strands(truth, ("full", "smear", "concatemer"))
    on: list[int] = []
    off: Counter = Counter()
    lines = fastq.read_text().splitlines()
    for i in range(0, len(lines), 4):
        name, seq = lines[i][1:], lines[i + 1]
        if name not in strands:
            continue
        read = seq if strands[name] == "+" else reverse_complement(seq)
        first_left = _anchor_hits(left, read)[:1]
        ends = [e for s, e in _anchor_hits(right, read) if first_left and s >= first_left[0][1]]
        if not ends:
            continue
        dist = ends[0] - first_left[0][0] - allele
        if abs(dist) <= max(30.0, 0.012 * allele):
            on.append(dist)
        off[round(dist / UNIT)] += 1
    n = max(len(on), 1)
    return {
        "n_on_peak": len(on),
        "per_onpeak": {k: off[k] / n for k in UNIT_STUTTER_KEYS},
        "on_peak_median_bp": sorted(on)[len(on) // 2] if on else None,
    }


def refit_unit_stutter(
    section: dict,
    molecules: dict,
    targets: dict,
    template: str,
    vntr: tuple[int, int],
    run: Sequencer,
    n: int,
    seed: int,
    rounds: int,
) -> tuple[dict, dict]:
    """Refit ``unit_stutter_rate``/``pmf`` so simulated reads reproduce the target statistic.

    The profiler counts a read -1 unit off when ``dist`` lies in [-90, -30] bp, which
    overlaps the on-peak window, so the target ratios also count the peak's own
    indel tail and short smears; `fit_unit_stutter` reads them as pure slippage.
    Here the target is the pooled ratio per allele unit times the template's units,
    measured on reads simulated with the profile's full molecule model. Each round
    moves the per-unit slip rate of each k by the remaining gap per unit.
    """
    lo, hi = vntr
    if (hi - lo) % UNIT or hi - lo < 2 * UNIT:
        raise SystemExit(f"--vntr {lo} {hi}: the VNTR must be whole {UNIT} bp units")
    units = (hi - lo) // UNIT
    rows = targets["per_allele_stutter_and_length_noise"]
    exposure = sum(r["n_on_peak"] * r["units"] for r in rows)
    target = {
        k: units * sum(r["n_on_peak"] * r[key] for r in rows) / exposure
        for k, key in UNIT_STUTTER_KEYS.items()
    }
    slip = {
        int(k): section["unit_stutter_rate"] * p for k, p in section["unit_stutter_pmf"].items()
    }
    history = []
    for rnd in range(rounds + 1):
        rate = sum(slip.values())
        if not rate:
            raise SystemExit(
                "unit-stutter refit: the template's reads already exceed every target ratio"
            )
        pmf = {str(k): round(v / rate, 5) for k, v in sorted(slip.items())}
        pmf["-1"] = round(pmf["-1"] + 1.0 - sum(pmf.values()), 5)  # absorb rounding drift
        section = {"unit_stutter_rate": round(rate, 6), "unit_stutter_pmf": pmf}
        model = MoleculeModel.from_dict({**molecules, **section, "offtarget_frac": 0.0})
        stats = _simulated(
            template,
            n,
            model,
            run,
            seed + 200 + rnd,
            lambda f, t: unit_offsets(f, t, template, vntr),
            range(lo, hi + 1, UNIT),
        )
        history.append(
            {
                **section,
                "per_onpeak": {str(k): round(v, 4) for k, v in stats["per_onpeak"].items()},
                "on_peak_median_bp": stats["on_peak_median_bp"],
            }
        )
        if rnd == rounds:
            break
        for k in slip:
            slip[k] = max(0.0, slip[k] + (target[k] - stats["per_onpeak"][k]) / units)
    record = {
        "method": "simulation refit: per-unit slip rate of each k moved by (target - "
        "simulated) / units per round, where target = pooled per_allele ratio per unit "
        "x template units and simulated = the real profiler's anchor-span rule "
        "(round(dist / 60) == k, per on-peak read) on reads of the template",
        "template_vntr_units": units,
        "target_per_onpeak": {str(k): round(v, 4) for k, v in target.items()},
        "rounds": history,
    }
    return section, record


def run_density(template: str, targets: dict, min_len: int) -> dict[str, list[float]]:
    """Homopolymer runs per 60 bp, ``[template, target]``, for every '+' target key.

    Stutter indels scale with these counts: a template with fewer long runs than the
    real alleles makes `fit_strand_rates` raise the channel rates to compensate.
    """
    bases = sum(lib["+"]["n_ref_bases"] for lib in targets["error_rates_per_lib"].values())
    found = Counter(f"{b}{e - s}" for s, e, b in homopolymer_runs(template, min_len))
    return {
        run: [round(UNIT * found[run] / len(template), 3), round(UNIT * row["n"] / bases, 3)]
        for key, row in sorted(targets["hp_P_obs_given_true"].items())
        for run in [key.partition("|")[0]]
        if key.endswith("|+") and int(run[1:]) >= min_len
    }


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
    ap.add_argument(
        "--fit-strand-rates",
        type=int,
        default=0,
        metavar="ROUNDS",
        help="empirical engine: fit errors.strand_rates in this many rounds (0 = off)",
    )
    ap.add_argument(
        "--fit-unit-stutter",
        type=int,
        nargs="?",
        const=0,
        default=None,
        metavar="ROUNDS",
        help="fit molecules.unit_stutter_rate/pmf from per_allele_stutter_and_length_noise; "
        "with ROUNDS > 0, then refit it in that many simulation rounds (needs --vntr)",
    )
    ap.add_argument(
        "--vntr",
        type=int,
        nargs=2,
        metavar=("START", "END"),
        help="the template's VNTR [START, END) (whole repeat units), for unit-stutter rounds",
    )
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    if args.engine == "pbsim3" and not args.model_file:
        ap.error("--engine pbsim3 requires --model-file")
    if args.engine == "pbsim3" and args.fit_strand_rates:
        ap.error("--fit-strand-rates needs --engine empirical")
    if args.fit_unit_stutter and not args.vntr:
        ap.error("--fit-unit-stutter ROUNDS needs --vntr")
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
        profile["errors"].pop("strand_rates", None)  # fitted for the old rates; refit below
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
        if args.fit_strand_rates:
            table = StutterTable.from_dict(section["stutter"], fallback=fallback)
            fit = fit_strand_rates(
                profile["errors"],
                targets,
                template,
                MoleculeModel(forward_frac=0.5, stutter=table),
                args.n_molecules,
                args.seed,
                max_len,
                args.fit_strand_rates,
            )
            profile["errors"]["strand_rates"] = fit["strand_rates"]
            report["strand_fit"] = fit["report"]
            report["strand_fit"]["runs_per_60bp"] = run_density(
                template, targets, error_model.hp_min_len
            )
            sequencer = EmpiricalSequencer(EmpiricalErrorModel.from_dict(profile["errors"]))
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
        pbsim_spectra, pbsim_err, _ = simulate(
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
        observed, total_err, _ = simulate(
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
    if args.fit_unit_stutter is not None:
        section, provenance["unit_stutter_fit"] = fit_unit_stutter(targets)
        provenance["unit_stutter_fit"]["target_key"] = args.target_key
        if args.fit_unit_stutter:
            section, provenance["unit_stutter_fit"]["refit"] = refit_unit_stutter(
                section,
                molecules,
                targets,
                template,
                tuple(args.vntr),
                sequencer,
                args.n_molecules,
                args.seed,
                args.fit_unit_stutter,
            )
        molecules.update(section)
    provenance["calibration_report"] = report
    args.out.write_text(json.dumps(profile, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
