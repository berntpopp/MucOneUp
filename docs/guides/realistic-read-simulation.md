# Realistic, truth-tracked long reads

MucOneUp can simulate long reads whose library and error properties match
real MUC1 data, and record the exact origin of every read. This is intended
for benchmarking callers such as allele reconstruction, phasing and
frameshift detection.

## Quick start

```bash
# Realistic ONT R10.4.1 PCR amplicons with per-read truth, FASTQ only
muconeup --config config.json reads amplicon --platform ont sample.001.simulated.fa \
  --read-profile ont_r10_sup_amplicon_v4 --coverage 3000 --seed 7 --no-align

# Genomic / targeted ONT reads from sampled fragments
muconeup --config config.json reads ont --simulator pbsim3-fragments sample.001.simulated.fa \
  --read-profile ont_r10_genomic_v2 --n-reads 800 --seed 7 --no-align

# List built-in profiles
muconeup reads profiles
```

Every run writes `{base}_read_truth.tsv.gz` next to the reads.

With a read profile, `--platform` (amplicon) defaults to the profile's platform
and `reads ont --simulator` defaults to `pbsim3-fragments`, so both flags are
optional in the examples above.

## What a read profile contains

| Part | Effect |
| --- | --- |
| `config_overrides` | Values merged into `ont_amplicon_params`, `pacbio_params`, `amplicon_params` (for example the PCR preset), validated against the config schema |
| `molecules` | Molecule model: strand mix, PCR smear (single-junction deletion products), whole-repeat-unit PCR stutter, between-allele chimeras, concatemers, off-target short products, a strand-aware homopolymer stutter table and its `stutter_fallback` for runs without a fitted entry |
| `errors` (optional) | Empirical error channel used instead of pbsim3 |
| `fragments` | Log-normal read-length model for `--simulator pbsim3-fragments` |
| `provenance` | Data source, derivation method, the stutter fit (`stutter_fit`: fitted and excluded keys, minimum n, fallback rule) and the validation report |

Precedence is **config file < read profile < explicit CLI flags**. For example,
`--pcr-preset no_bias` overrides the profile's PCR preset. When a profile sets
a PCR preset, preset parameters from the config file (`alpha`, `e_max`,
`cycles`, `denaturation_time`) are ignored; `stochastic` is kept. A profile's
`fragments` lengths replace `ont_fragment_params.length_median`/`length_sigma`
from the config file, and `--read-length-median`/`--read-length-sigma` replace
the profile's values. A profile can also
be named in the config file as `"read_model": {"profile": "<name or path>"}`.

## Built-in profiles

| Profile | Engine | Calibrated against | Validation |
| --- | --- | --- | --- |
| `ont_r10_sup_amplicon_v1` | empirical | 9 public R10.4.1 MUC1 amplicon libraries (ENA PRJEB92208) | total error 2.27% vs 2.14%; fitted keys C3–C8 and G3–G4 on both strands within 0.015 (dupC C8: p_correct 0.36/0.73 on +/− vs 0.36/0.75) |
| `ont_r10_genomic_v1` | empirical | 2 public R10.4.1 WGS runs (PRJEB92208) | total error 0.46% vs 0.40%; fitted keys C3, C4, C7 and G3 within 0.003; C5, C6 and C8 extrapolated |
| `ont_r10_sup_amplicon_v2` | empirical | as `_v1`, plus per-strand error rates | `_v1` plus `errors.strand_rates`: every per-strand mismatch, insertion and deletion rate within 11% of the target on the MucSim-Bench standard set (30 cases; `_v1`: −30% to +39%) |
| `ont_r10_sup_amplicon_v3` | empirical | as `_v2`, plus whole-unit PCR stutter | `_v2` plus `molecules.unit_stutter_rate`/`unit_stutter_pmf`: span-offset Jensen-Shannon distance 0.162 / 0.181 (< 55 / ≥ 55 units; `_v2`: 0.186 / 0.187; target ≤ 0.1); per-strand error rates unchanged (30 cases) |
| `ont_r10_sup_amplicon_v4` | empirical | as `_v3`, refitted | `_v3` with `errors.strand_rates` and the unit stutter refitted: span-offset Jensen-Shannon distance 0.075 / 0.076 (`_v3`: 0.086 / 0.107; symmetric scoring, see [Calibration status](#calibration-status)); every per-strand rate within 2% (30 cases) |
| `ont_r10_genomic_v2` | empirical | as `_v1`, plus per-strand error rates | `_v1` plus `errors.strand_rates`: every per-strand rate within 10% (30 cases; `_v1`: −27% to +40%). Read lengths are not fitted, so the spanning fraction is not calibrated |
| `hifi_amplicon_v1` | pbsim3 + ccs | none (generic) | strand mix only; not calibrated to real HiFi data |

Use `ont_r10_sup_amplicon_v4` and `ont_r10_genomic_v2` for new work. Older
profiles are kept unchanged so that runs made with them reproduce: `_v2`
differs from `_v1` only in the error channel's strand symmetry, the amplicon
`_v3` differs from `_v2` only in the unit stutter, and `_v4` differs from
`_v3` only in the refitted values of both (#131). There is no
genomic `_v3`: WGS libraries have no PCR, and the public WGS target has a
single allele with 27 on-peak reads, too few to fit a stutter rate.

The genomic profiles' read lengths are a generic default (median 6 kb,
sigma 0.5), not a fit: the two public WGS runs have read-length medians of
about 107 kb and 0.9 kb, which no single log-normal fits. Set
`--read-length-median`/`--read-length-sigma` for your library.

The amplicon profile reproduces the measured library composition:

- smear products in about 24% of amplicon slots;
- chimeras 2.3% and concatemers 2.4%;
- whole-unit PCR stutter: a full product loses or gains whole repeat units
  with probability 0.0016 per VNTR unit in `_v4` (about 9% of a 60-unit
  allele; −1 / −2 units in 77 / 23% of cases), 0.0024 in `_v3` (−1 / −2 / +1
  units in 78 / 16 / 6%);
- about 30% off-target reads (median ~367 bp);
- the measured PCR length bias (`madritsch2025_r10`: ln(long/short) ≈ −0.056 per repeat unit).

Profiles contain aggregate statistics only, never reads or sequences.

## Calibration status

Measured on the MucSim-Bench v4 `standard` designs (30 ONT amplicon or 30 ONT
genomic cases; same haplotypes, seeds and amounts for every profile) and
scored with MucOneSpan `benchsim.realism` against the PRJEB92208 targets.
Error rates are the median over cases, as a deviation from the per-strand
target. Span JS is the Jensen-Shannon distance of the span-offset histogram
(< 55 / ≥ 55 units, limit 0.1) under two scorings: *old* (MucOneSpan `main`;
MucOneUp #208, #210) and *symmetric* (MucOneSpan #93: the span is measured with
the real profiler's first-anchor rule, and only alleles inside the real 40–83
unit support are scored).

| Profile | Per-strand error rates vs target | Span JS, old scoring | Span JS, symmetric scoring | Between alleles (target median 0.023) |
| --- | --- | --- | --- | --- |
| `ont_r10_sup_amplicon_v1` | −30% to +39% (#208) | 0.222 / 0.214 (#208) | not measured | 0.0 median (#208) |
| `ont_r10_sup_amplicon_v2` | −3% to +7%; all-strand deletion +11% (#208) | 0.186 / 0.187 (#208) | 0.133 / 0.163 (MucOneSpan #93) | 0.0 median (#208) |
| `ont_r10_sup_amplicon_v3` | as `_v2`: `+` deletion +7% (#210) | 0.162 / 0.181 (#210) | 0.086 / 0.107 (MucOneSpan #93) | 0.053 median of the 8 cases with alleles ≥ 6 units apart (this release) |
| `ont_r10_sup_amplicon_v4` | −1% to +1%; all-strand deletion +7% (this release) | 0.168 / 0.141 (this release) | **0.075 / 0.076** (this release) | 0.050 median of the 8 cases with alleles ≥ 6 units apart (this release) |
| `ont_r10_genomic_v1` | −27% to +40% (#208) | no span target | no span target | no span target |
| `ont_r10_genomic_v2` | −7% to +10% (#208) | no span target | no span target | no span target |

What the `_v4` refit changed, in plain words:

- **Per-strand error rates.** Most homopolymer indels come from the stutter
  table, so the error channel's own rates only fill the gap to the target.
  `_v2` fitted that gap on a template with fewer C7 runs (about 0.54 per
  60 bp) than real MUC1 alleles (0.62), so on real-like alleles the `+` strand
  lost 7% too many bases. `_v4` fits on a template with real-like runs
  (0.626 per 60 bp); its `+` channel deletion rate drops from 0.0047 to
  0.0037. On simulated alleles inside the real support, the span bias (median
  span minus allele length) now grows by −0.060 bp per unit, against −0.067
  real (`_v3`: −0.076).
- **Whole-unit stutter.** The real profiler counts a read as one unit short
  when its span is 30–90 bp short, so the real ratio also counts reads from the
  left edge of the main peak. `_v3` read that ratio as pure slippage: measured
  the profiler's way on a 73-unit template, its −1 unit ratio was 0.20 against
  0.16 in the target, and the −97 to −52 bp bins (≥ 55 units) held 0.112 of the reads
  against 0.075 real. `_v4` fits the rate so that simulated reads, measured
  the same way, give the real ratio: 0.0016 per unit instead of 0.0024, and
  those bins now hold 0.077. The +1 unit ratio is already reached by the
  peak's right edge, so its fitted share is 0.

The span and per-strand error criteria of #131 are met with symmetric scoring.
Between alleles is not: in the 8 standard cases whose alleles differ by
6 units or more (as in all 9 real libraries), the simulated fraction is about
twice the real one (0.050 against 0.016–0.026 for real libraries with
12–22 units difference). The other 22 cases differ by less than 6 units
(16 by 2 units or less), so the between-alleles window, more than 1.5 units
from both alleles, is narrow or empty, and the median over all 30 cases is 0.

## Why an empirical error channel?

Measured with MucOneUp's truth-tracked path, pbsim3's ONT models
(QSHMM-ONT-HQ, QSHMM-ONT, ERRHMM-ONT) have an error floor of about 3.5–4% in
template mode, with C7 read correctly only 78–84% of the time. Real R10.4.1 sup
MUC1 amplicons show about 2.1% error, with C7 read correctly 52% (+ strand) and
89% (− strand) of the time.

Injected stutter can add noise but cannot remove pbsim3's own, so realistic ONT
profiles use the empirical channel:

- mismatch, insertion and deletion rates outside homopolymer runs, optionally
  per read strand (`errors.strand_rates`, keys `+`/`-`), because R10.4.1
  errors are strand-asymmetric (for example deletions 1.14% on `+` and 0.49%
  on `−` in the amplicon target);
- measured indel lengths;
- a log-normal per-read error spread;
- Phred qualities scaled for basecaller over-confidence.

Homopolymer runs are left to the stutter table, so that table equals the
measured spectrum. pbsim3 remains the default when no profile is used.

## Homopolymer stutter coverage

The empirical channel applies no base-level errors inside runs of at least
`errors.hp_min_len` bases (default 3). The stutter table is their only error
source, so every such run must resolve to a stutter pmf. A profile with an
`errors` section is rejected unless it defines `molecules.stutter` and
`molecules.stutter_fallback`.

Stutter keys are `{base}{length}|{strand}` in haplotype orientation, e.g.
`C8|+`. Case is ignored: a soft-masked `cccC` is a C4 run, and stutter keeps
the run's own case. Runs of `N` (or any non-ACGT character) are unknown
sequence, not homopolymers: they get no stutter and are not protected from
base-level errors. The built-in ONT profiles fit every per-strand key of the target with
enough observations (`provenance.stutter_fit.min_target_n`: 500 for the
amplicon profile, 300 for the genomic profile), including the dupC C8 run in
the amplicon profile. Runs without a fitted entry use the fallback:

```json
"stutter_fallback": {
  "rule": "log_odds_interpolate",
  "log_odds_slope_per_base": 0.57825,
  "max_extrapolation_bases": 3,
  "generic": {
    "ref_len": 3,
    "pmf": {"-3": 0.00072, "-2": 0.0034, "-1": 0.02715, "0": 0.94838,
            "1": 0.01744, "2": 0.00266, "3": 0.00025}
  }
}
```

(values from `ont_r10_sup_amplicon_v1`)

The error log-odds of a pmf are log(P(Δ≠0)/P(Δ=0)); its error shape is the
pmf conditional on an error. Deltas that would remove more bases than the run
has are dropped.

1. **Fitted lengths.** The fitted entries of the same base and strand (or
   strand `both` if the base has no strand-specific entries).
2. **Interpolation.** A run between two fitted lengths takes the log-odds
   interpolated linearly between the nearest fitted lengths on both sides,
   with the error shapes mixed by the same weights. For example the genomic
   profile's C5 and C6 lie between its fitted C4 and C7.
3. **Extrapolation.** A run longer than the longest fitted length (or shorter
   than the shortest) multiplies that entry's error odds by
   `exp(log_odds_slope_per_base × d)` and keeps its error shape. `d` is the
   length difference capped at `max_extrapolation_bases`, so a T19 flank run
   is treated like a run three bases beyond the data instead of being
   misread almost always.
4. **Generic.** A base without any fitted entry (for example A and T runs,
   which the PRJEB92208 targets do not cover at length ≥ 3) uses
   `generic.pmf` at `generic.ref_len`, scaled as in step 3. MucOneUp logs a
   warning once per base. The generic pmf is also used when a reference keeps
   no error delta for a short run, so a run resolved by the fallback is never
   error-free.

The calibration helper derives both values from the fitted table: the slope is
the median least-squares slope of error log-odds against run length over
base/strand groups with at least three fitted lengths (clamped at 0), and the
generic pmf is the n-weighted pool of the fitted length-3 entries.
`max_extrapolation_bases` is configured (`--max-extrapolation-bases`, default
3). Profiles
without an `errors` section (pbsim3 path) may omit the fallback; missing keys
then get no extra stutter, because pbsim3 adds its own homopolymer errors.

## Truth manifest

`{base}_read_truth.tsv.gz` has one row per read, in FASTQ order:

| Column | Meaning |
| --- | --- |
| `read_id` | `{base}_h{hap}_m{molecule}` (unique) |
| `hap` | 1-based haplotype of origin (off-target products: the haplotype whose flank they were cut from) |
| `molecule` | Molecule id |
| `kind` | `full`, `smear`, `chimera`, `concatemer`, `offtarget` or `fragment`. Unit-stutter products are `smear` (an off-peak product of one parent, so they never count as correct-length `full` reads); `detail` tells them apart |
| `strand` | `+`: read in haplotype orientation (MUC1 C-runs read as C); `-`: reverse complement |
| `src_start`, `src_end` | 0-based, end-exclusive source interval: amplicon products use amplicon coordinates, off-target products haplotype coordinates, fragments coordinates in the fragment source (left flank + haplotype + right flank, so haplotype position = `src_start` − left flank length when `--flank-fasta` is used) |
| `n_hp_edits`, `hp_edits` | Injected homopolymer changes as `pos:base:true>new`; `pos` is 0-based in the product before strand orientation (after a smear deletion, chimera junction or concatemer join), not in amplicon coordinates |
| `detail` | Smear deletion interval (`deletion:a-b`), unit stutter (`unit_stutter:k:a-b`: whole units `[a, b)` deleted for k < 0 or duplicated for k > 0, amplicon coordinates) or chimera junction (`junction:a|b`) |

The run's `*_metadata.tsv` records `Read_profile`, `Read_profile_sha256` and
`Read_truth` (the manifest file name). Fragment runs record `Fragment_reads`
(or `Coverage`), `Fragment_length_median`, `Fragment_length_sigma` and
`Flank_fasta` instead of the NanoSim read-length rows.

`--track-read-source` without a profile uses a no-op molecule model. Reads then
get truth but keep the legacy error model.

## Calibrating your own profile

`helpers/calibrate_read_profile.py` derives a profile from aggregate targets
(homopolymer spectra, context-split error rates, indel lengths, per-read error
quantiles). It simulates through MucOneUp, re-measures, and stores the
validation in the profile's provenance. It needs the `edlib` Python package.

With `--engine empirical` it fits every per-strand homopolymer key with at
least `--min-target-n` observations (default 500), keeping deltas up to
`--max-delta` (default 6), and derives `stutter_fallback`
(`--min-lengths-for-slope`, `--generic-ref-len`). The validation report
(`calibration_report.validation_round1.per_key`) compares the simulated and
target observed-minus-true pmf for every base/length/strand of the target:
p_correct, the largest absolute residual and a pass/fail against
`--tolerance` (default 0.03). Keys whose run is not in the template are marked
`not_in_template`, so validate with a template that contains the runs you care
about (for MUC1, the amplicon of a dupC haplotype for C8).

`--fit-strand-rates ROUNDS` (empirical engine) then fits `errors.strand_rates`:
each round simulates the template, measures the per-strand mismatch, insertion
and deletion rates (stutter included) and moves each channel rate by the gap to
`error_rates_per_ref_base["<strand>:<type>"]` divided by the template's
non-homopolymer fraction. The rounds are stored in
`calibration_report.strand_fit`. Stutter adds most of the homopolymer indels,
so the fitted channel rates depend on how many long runs the template has.
`strand_fit.runs_per_60bp` records the template's runs per 60 bp against the
target's (`[template, target]` per run, e.g. `"C7": [0.626, 0.622]`); pick a
template whose long runs match the target.

`--fit-unit-stutter` (either engine) fits whole-repeat-unit PCR stutter from
`per_allele_stutter_and_length_noise`, the per-allele −2, −1 and +1 unit
products per on-peak read. A full product slips with probability
`unit_stutter_rate` × VNTR units, so the rate is the pooled number of slipped
products over the pooled molecule units, and `unit_stutter_pmf` is the pooled
share of each k. The fit is stored in `provenance.unit_stutter_fit`. Slips
delete or duplicate whole units at boundaries read from the haplotype: the
VNTR lies between the config's left and right constants, and units are matched
against the config repeats from both ends, so a mutated unit is never cut.
The real alleles span 40–83 units; longer alleles use the same per-unit rate.

`--fit-unit-stutter ROUNDS --vntr START END` then refits the stutter by
simulation. The real profiler counts a read as −1 unit off when its span is
30–90 bp short, which overlaps the on-peak window, so the target ratios also
contain the peak's own indel tail. Each round simulates the template (VNTR at
`[START, END)`, whole 60 bp units) with the full molecule model, measures the
same ratios with the profiler's anchor rule, and moves each k's per-unit rate
by the remaining gap. The rounds are stored in
`provenance.unit_stutter_fit.refit`.

For `--engine pbsim3`, it deconvolves the stutter table against pbsim3's own
homopolymer errors instead.

## Performance

- **Alignment:** use `--no-align` when FASTQ is enough. When alignment is
  needed, build the index once so it is reused
  (`minimap2 -x map-ont -d GRCh38.fa.map-ont.mmi GRCh38.fa`). Only an index
  named `{reference}.{preset}.mmi` and newer than the FASTA is used.
- **Speed:** without alignment, the empirical ONT profile produces about 3,000
  amplicon reads in a few seconds.
- **Genomic flanks:** fragments are sampled as from a window of a longer
  genome: coverage is uniform along the source, and reads that overlap a source
  end are clipped there (flanks are 10 kb by default). `--coverage` is the mean
  depth. Use `--flank-fasta` with `left`/`right` records to extend the flanks
  so that fewer reads are clipped.
