# Realistic, truth-tracked long reads

MucOneUp can simulate long reads whose library and error properties match
real MUC1 data, and record the exact origin of every read. This is intended
for benchmarking callers such as allele reconstruction, phasing and
frameshift detection.

## Quick start

```bash
# Realistic ONT R10.4.1 PCR amplicons with per-read truth, FASTQ only
muconeup --config config.json reads amplicon --platform ont sample.001.simulated.fa \
  --read-profile ont_r10_sup_amplicon_v1 --coverage 3000 --seed 7 --no-align

# Genomic / targeted ONT reads from sampled fragments
muconeup --config config.json reads ont --simulator pbsim3-fragments sample.001.simulated.fa \
  --read-profile ont_r10_genomic_v1 --n-reads 800 --seed 7 --no-align

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
| `molecules` | Molecule model: strand mix, PCR smear (single-junction deletion products), between-allele chimeras, concatemers, off-target short products, a strand-aware homopolymer stutter table and its `stutter_fallback` for runs without a fitted entry |
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
| `hifi_amplicon_v1` | pbsim3 + ccs | none (generic) | strand mix only; not calibrated to real HiFi data |

The amplicon profile reproduces the measured library composition:

- smear products in about 24% of amplicon slots;
- chimeras 2.3% and concatemers 2.4%;
- about 30% off-target reads (median ~367 bp);
- the measured PCR length bias (`madritsch2025_r10`: ln(long/short) ≈ −0.056 per repeat unit).

Profiles contain aggregate statistics only, never reads or sequences.

## Why an empirical error channel?

Measured with MucOneUp's truth-tracked path, pbsim3's ONT models
(QSHMM-ONT-HQ, QSHMM-ONT, ERRHMM-ONT) have an error floor of about 3.5–4% in
template mode, with C7 read correctly only 78–84% of the time. Real R10.4.1 sup
MUC1 amplicons show about 2.1% error, with C7 read correctly 52% (+ strand) and
89% (− strand) of the time.

Injected stutter can add noise but cannot remove pbsim3's own, so realistic ONT
profiles use the empirical channel:

- mismatch, insertion and deletion rates outside homopolymer runs;
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
`C8|+`. The built-in ONT profiles fit every per-strand key of the target with
enough observations (`provenance.stutter_fit.min_target_n`: 500 for the
amplicon profile, 300 for the genomic profile), including the dupC C8 run in
the amplicon profile. Runs without a fitted entry use the fallback:

```json
"stutter_fallback": {
  "rule": "log_odds_linear",
  "log_odds_slope_per_base": 0.57825,
  "generic": {
    "ref_len": 3,
    "pmf": {"-3": 0.00072, "-2": 0.0034, "-1": 0.02715, "0": 0.94838,
            "1": 0.01744, "2": 0.00266, "3": 0.00025}
  }
}
```

(values from `ont_r10_sup_amplicon_v1`)

1. **Reference.** The fitted entry of the same base and strand (or strand
   `both` if the base has no strand-specific entries) at the longest fitted
   length that is not longer than the run; if every fitted length is longer,
   the shortest.
2. **Generic.** A base without any fitted entry (for example A and T runs,
   which the PRJEB92208 targets do not cover at length ≥ 3) uses
   `generic.pmf` at `generic.ref_len`. MucOneUp logs a warning once per base.
3. **Scaling.** The reference's error odds P(Δ≠0)/P(Δ=0) are multiplied by
   `exp(log_odds_slope_per_base × (run length − reference length))`. The shape
   of the error deltas is kept, and deltas that would remove more bases than
   the run has are dropped.

The calibration helper derives both values from the fitted table: the slope is
the median least-squares slope of error log-odds against run length over
base/strand groups with at least three fitted lengths (clamped at 0), and the
generic pmf is the n-weighted pool of the fitted length-3 entries. Profiles
without an `errors` section (pbsim3 path) may omit the fallback; missing keys
then get no extra stutter, because pbsim3 adds its own homopolymer errors.

## Truth manifest

`{base}_read_truth.tsv.gz` has one row per read, in FASTQ order:

| Column | Meaning |
| --- | --- |
| `read_id` | `{base}_h{hap}_m{molecule}` (unique) |
| `hap` | 1-based haplotype of origin (off-target products: the haplotype whose flank they were cut from) |
| `molecule` | Molecule id |
| `kind` | `full`, `smear`, `chimera`, `concatemer`, `offtarget` or `fragment` |
| `strand` | `+`: read in haplotype orientation (MUC1 C-runs read as C); `-`: reverse complement |
| `src_start`, `src_end` | 0-based, end-exclusive source interval: amplicon products use amplicon coordinates, off-target products haplotype coordinates, fragments coordinates in the fragment source (left flank + haplotype + right flank, so haplotype position = `src_start` − left flank length when `--flank-fasta` is used) |
| `n_hp_edits`, `hp_edits` | Injected homopolymer changes as `pos:base:true>new`; `pos` is 0-based in the product before strand orientation (after a smear deletion, chimera junction or concatemer join), not in amplicon coordinates |
| `detail` | Smear deletion interval or chimera junction |

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
