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

## What a read profile contains

| Part | Effect |
| --- | --- |
| `config_overrides` | Values merged into `ont_amplicon_params`, `pacbio_params`, `amplicon_params` (for example the PCR preset), validated against the config schema |
| `molecules` | Molecule model: strand mix, PCR smear (single-junction deletion products), between-allele chimeras, concatemers, off-target short products, and a strand-aware homopolymer stutter table |
| `errors` (optional) | Empirical error channel used instead of pbsim3 |
| `fragments` | Log-normal read-length model for `--simulator pbsim3-fragments` |
| `provenance` | Data source, derivation method, uncovered homopolymer keys and the validation report |

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
| `ont_r10_sup_amplicon_v1` | empirical | 9 public R10.4.1 MUC1 amplicon libraries (ENA PRJEB92208) | total error 2.26% vs 2.14%; homopolymer keys C3–C7 and G3–G4 on both strands within 0.028 |
| `ont_r10_genomic_v1` | empirical | 2 public R10.4.1 WGS runs (PRJEB92208) | total error 0.43% vs 0.40%; residuals ≤ 0.0035 |
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

## Truth manifest

`{base}_read_truth.tsv.gz` has one row per read, in FASTQ order:

| Column | Meaning |
| --- | --- |
| `read_id` | `{base}_h{hap}_m{molecule}` (unique) |
| `hap` | 1-based haplotype of origin (off-target products: the haplotype whose flank they were cut from) |
| `molecule` | Molecule id |
| `kind` | `full`, `smear`, `chimera`, `concatemer`, `offtarget` or `fragment` |
| `strand` | `+`: read in haplotype orientation (MUC1 C-runs read as C); `-`: reverse complement |
| `src_start`, `src_end` | 0-based, end-exclusive source interval: amplicon products use amplicon coordinates, off-target products haplotype coordinates |
| `n_hp_edits`, `hp_edits` | Injected homopolymer changes as `pos:base:true>new` (source orientation) |
| `detail` | Smear deletion interval or chimera junction |

`--track-read-source` without a profile uses a no-op molecule model. Reads then
get truth but keep the legacy error model.

## Calibrating your own profile

`helpers/calibrate_read_profile.py` derives a profile from aggregate targets
(homopolymer spectra, context-split error rates, indel lengths, per-read error
quantiles). It simulates through MucOneUp, re-measures, and stores the
residuals in the profile's provenance. It needs the `edlib` Python package.

For `--engine pbsim3`, it deconvolves the stutter table against pbsim3's own
homopolymer errors instead.

## Performance

- **Alignment:** use `--no-align` when FASTQ is enough. When alignment is
  needed, build the index once so it is reused
  (`minimap2 -x map-ont -d GRCh38.fa.map-ont.mmi GRCh38.fa`).
- **Speed:** without alignment, the empirical ONT profile produces about 3,000
  amplicon reads in a few seconds.
- **Genomic flanks:** fragment simulation clips reads at the ends of the
  simulated flanks (10 kb by default). Use `--flank-fasta` with `left`/`right`
  records to extend them.
