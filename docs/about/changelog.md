# Changelog

All notable changes to MucOneUp will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

## [0.46.0] - 2026-09-24

### Breaking change: read profile format
Custom read profiles with an `errors` section (empirical error channel) must
be updated (#132). Without these, homopolymer runs could be simulated
error-free, so such profiles are now rejected with a message that points to
the calibration helper and the
[Homopolymer stutter coverage](../guides/realistic-read-simulation.md#homopolymer-stutter-coverage)
section:
- `molecules.stutter` and `molecules.stutter_fallback` are required, with
  `rule` (`log_odds_interpolate`), `log_odds_slope_per_base`,
  `max_extrapolation_bases` and `generic` (`ref_len`, `pmf`);
- a fitted stutter entry for a protected run (length >= `errors.hp_min_len`)
  must have P(delta = 0) < 1, and the generic pmf must keep error mass at the
  table's minimum run length;
- `errors.hp_min_len` must not be below the stutter table's minimum run length.

Regenerate with `helpers/calibrate_read_profile.py --engine empirical`, or
copy `stutter_fallback` from a built-in profile. Profiles without `errors`
(pbsim3 path) are unaffected.

### Fixed
- **Homopolymer runs without a stutter entry were simulated error-free**
  (#132). The empirical error channel protects every run >= `hp_min_len` from
  base-level errors, and a missing stutter key drew no length change, so the
  dupC C8 run, all A/T runs and G runs >= 5 of the ONT profiles were read
  perfectly. Now:
  - `ont_r10_sup_amplicon_v1` fits every per-strand key of the PRJEB92208
    amplicon target with n >= 500, adding `C8|+` and `C8|-` (p_correct 0.362
    and 0.747); `ont_r10_genomic_v1` is refit from the WGS target (n >= 300),
    which has no C8 key;
  - new `molecules.stutter_fallback` (rule `log_odds_interpolate`) resolves
    runs without a fitted entry: it interpolates the error log-odds between
    the nearest fitted lengths of the same base and strand, extrapolates
    beyond the fitted range (at most `max_extrapolation_bases`, 3 in the
    built-ins, so long flank runs such as T19 are not almost always
    misread), or uses a generic pmf with a warning (once per base) when the
    base has no fitted entries;
  - profiles with an `errors` model are validated so no protected run is
    error-free by construction (see **Breaking change** above);
  - homopolymer runs ignore case (a soft-masked `cccC` is C4 and keeps its
    case when stuttered); runs of `N` get no stutter and are not protected
    from base-level errors;
  - `provenance.stutter_fit` records the fitted and excluded keys, the minimum
    n and the fallback derivation; runs resolved by interpolation or
    extrapolation are logged at INFO once per profile load.
- `helpers/calibrate_read_profile.py` validates the simulated observed-minus-
  true pmf against the target for every base/length/strand
  (`calibration_report.validation_round*.per_key`) and gains `--min-target-n`
  (default now 500), `--max-delta` (default now 6, was 3),
  `--min-lengths-for-slope`, `--generic-ref-len`, `--max-extrapolation-bases`
  and `--tolerance`;
  `--model-file` is only required for `--engine pbsim3`. With
  `--engine pbsim3` it drops an inherited `errors` model and
  `stutter_fallback`; with `--engine empirical` it requires
  `errors.hp_min_len` to equal the stutter table's minimum run length (3).
  Fitting fails with a message instead of dividing by zero when a target key
  has no mass within `--max-delta`, and rounding never makes P(0) negative.

### Changed (intentional output differences)
- Reads simulated with `ont_r10_sup_amplicon_v1` or `ont_r10_genomic_v1`
  change for the same seed: all homopolymer runs are now stuttered and the
  refit tables keep deltas up to ±6. Simulations without a read profile are
  unchanged (seeded legacy ONT and HiFi amplicon FASTQs are byte-identical to
  0.45.0). User profiles with an `errors` section but no
  `stutter_fallback` are now rejected (see **Breaking change** above).

### Tests
- Legacy byte-identity guard: a unit test pins the template FASTA md5s and
  tool argument lists of seeded legacy ONT and HiFi amplicon runs (tools
  mocked), and an integration test checks the golden FASTQ md5s with real
  pbsim3/ccs/samtools (skips when they are missing).

---

## [0.45.0] - 2026-09-24

Realistic, truth-tracked long-read simulation. See the
[Realistic Truth-Tracked Reads](../guides/realistic-read-simulation.md) guide.
Simulations without a read profile or `--track-read-source` produce the same
reads (a seeded legacy ONT amplicon FASTQ is byte-identical to 0.44.5). The bug
fixes listed under **Changed** alter outputs of some existing configurations on
purpose. Requires samtools >= 1.13 (`view --subsample`).

### Added
- **Read profiles** (`--read-profile NAME|PATH`, `muconeup reads profiles`):
  versioned JSON bundles of config overlays, a molecule model and an optional
  empirical error model. Built-ins: `ont_r10_sup_amplicon_v1` and
  `ont_r10_genomic_v1` (calibrated to public PRJEB92208 aggregates) and
  `hifi_amplicon_v1` (generic). Refs #103
- **Molecule model**: strand mix, PCR smear, chimeras, concatemers, off-target
  products and strand-aware homopolymer stutter (#101, #103)
- **Per-read truth**: unique read names `{base}_h{hap}_m{molecule}` and
  `{base}_read_truth.tsv.gz` for amplicon and fragment simulation.
  `--track-read-source` is now supported in amplicon mode (#100)
- **Empirical error channel**, calibrated to R10.4.1 sup, as an alternative to
  pbsim3; pbsim3's ONT models floor at ~3.5–4% error in template mode (#103)
- **Genomic ONT from fragments**: `reads ont --simulator pbsim3-fragments`
  (`--n-reads`, `--read-length-median`, `--read-length-sigma`,
  `--flank-fasta`) (#107)
- **PCR preset `madritsch2025_r10`** reproducing ln(long/short) ≈ −0.056 per
  repeat unit (#104)
- **pbsim3 options**: `difference_ratio` (ONT amplicon) and read id prefix,
  passed only when set (#105). pbsim3 has no `--accuracy-sd` option, so
  `ont_amplicon_params` does not accept `accuracy_sd` (#115)
- **`--no-align`** for `reads amplicon`, `reads ont` and `reads pacbio` (#106)
- `helpers/calibrate_read_profile.py` to derive and validate profiles

### Fixed
- **minimap2**: a prebuilt `{reference}.{preset}.mmi` that is not older than
  the FASTA is now reused instead of re-indexing the whole genome on every run
  (#106). Generic `{reference}.mmi` files are ignored because their preset is
  unknown (#117)
- **Seed provenance**: `simulate --seed N` is recorded in `provenance.seed` (#109)
- **NanoSim read-source tracking**: `haplotype-N` names are parsed and the
  required companion keys are written (#102)
- **samtools downsampling**: now uses `--subsample`/`--subsample-seed`; the old
  `-s SEED.FRAC` form misread fractions such as 0.05 (#97)
- **`reads ont --min-read-length`**: no longer silently overrides the config (#108)
- **SNP info** in `ReadSourceTracker.from_simulation_results` is keyed by
  haplotype (#111)
- **Amplicon housekeeping** (#110):
  - metadata TSV reports amplicon settings;
  - PacBio amplicon mode honours `keep_intermediate_files`;
  - docs no longer reference the unshipped `QSHMM-SEQUEL.model`;
  - the pbsim3 version probe no longer logs a spurious ERROR.

### Changed (intentional output differences for existing configurations)
- samtools downsampling uses the exact fraction: `-s SEED.FRAC` turned 0.05
  into ~50%, and the wrapper truncated fractions to 4 decimals (#97)
- `reads ont` keeps `nanosim_params.min_read_length` from the config instead of
  overwriting it with 100 when `--min-read-length` is omitted (#108)
- `simulate --seed` is recorded in provenance, which changes
  `config_fingerprint` and adds a `Seed` row (#109)
- Amplicon metadata TSV reports `Template_molecules`, model and PCR rows
  instead of NanoSim/PacBio WGS rows; PacBio amplicon mode honours
  `keep_intermediate_files` (#110)
- The NanoSim read-name parser raises on names without a haplotype instead of
  assigning haplotype 1 (#102)
- A generic `{reference}.mmi` next to the reference is no longer used (#117)

### Pre-release review fixes
An independent review of the release branch found and fixed:
- off-target molecules sampled across artificial flank/haplotype junctions;
  they now come from real flank intervals with haplotype coordinates (#112)
- truth-tracked runs yielding zero reads now fail; zero fragment requests are
  rejected (#113)
- config-file fragment lengths and PCR parameters no longer override the read
  profile (#114)
- `--accuracy-sd` removed: pbsim3 rejects it (#115)
- `reads ont --no-align` now skips NanoSim alignment (#116)
- empirical deletions no longer remove protected homopolymer bases (#118)
- smear retention bounds (`smear_min_keep` in [0, 0.95], at least one base) (#119)
- `--platform`/`--simulator` default to the read profile's platform and to
  `pbsim3-fragments` with a profile (#120)
- the empirical error channel uses its own seed stream (#121)
- metadata records `Read_profile`, `Read_profile_sha256`, `Read_truth` and
  fragment settings (#122)
- fragment sampling gives uniform coverage along the source; truth coordinate
  systems are documented (#123)
- malformed read profiles raise validation errors; full-run homopolymer
  deletions in the stutter table are applied instead of clamped (#124)

---

## [0.44.3] - 2026-04-07

### Fixed
- ONT amplicon pipeline now reads from dedicated `ont_amplicon_params` config section instead of `pacbio_params` -- was using PacBio Sequel model (ERRHMM-SEQUEL) for ONT reads (#79)
- CLI routes `--model-file` overrides to `ont_amplicon_params` for ONT and `pacbio_params` for PacBio

### Added
- `ont_amplicon_params` config section with ONT-specific defaults (QSHMM-ONT-HQ.model)
- `OntAmpliconConfig` TypedDict for type-safe ONT parameter access
- Config schema validation and path resolution for `ont_amplicon_params`
- Test for ONT platform using `ont_amplicon_params` (not `pacbio_params`)
- Documentation: ONT configuration section in amplicon simulation guide

---

## [0.43.4] - 2026-04-06

### Fixed
- Race condition in VNTR efficiency temp directory — replaced deterministic path with `tempfile.TemporaryDirectory` for process-safe isolation in parallel runs (#82)

---

## [0.43.3] - 2026-04-06

### Fixed
- ONT amplicon: warn when model file does not appear to be an ONT model (#79)
- pbsim3 wrapper: detect and convert numbered SAM files (`prefix_0001.sam`) produced with template mode (#79)
- Config paths (model_file, reseq_model, human_reference) now resolved relative to config file directory, not working directory (#80)
- Default `expected_product_range` widened from `[1500, 6000]` to `[1500, 15000]` to support VNTRs with >99 repeats (#81)

### Added
- 5 new tests for config path resolution, ONT model warning, and numbered SAM detection

---

## [0.43.2] - 2026-04-06

### Added
- `--read-number` CLI option for `reads illumina` to control fragment generation count
- 3 new CLI tests for `--read-number` option

### Changed
- Default `read_number` raised from 10,000 to 100,000 in config.json (was limited by original ReSeq deadlock bug, now resolved by ReSeq2)

---

## [0.43.1] - 2026-04-05

### Changed
- Replaced ReSeq with [ReSeq2](https://github.com/berntpopp/ReSeq2) as the recommended Illumina error modeling tool
- Removed `seqToIllumina` timeout workaround — ReSeq2 v2.0.3 fixes the upstream deadlock bug that limited the original ReSeq to ~10,000 reads
- Default `seqToIllumina` timeout raised from 120s to 600s (safety net only)

### Removed
- `seqtoillumina_timeout` config parameter (no longer needed)

---

## [0.43.0] - 2026-04-05

### Added
- ONT amplicon simulation via `reads amplicon --platform ont`
- `--platform` option for `reads amplicon` command (choices: `pacbio`, `ont`)
- Shared amplicon preparation module (`amplicon_common.py`) used by both platforms
- `Assay_type` field in simulation metadata TSV
- `assay_type` and `ont-amplicon` to config schema
- `keep_intermediate_files` support in ONT amplicon pipeline
- 13 new tests for ONT amplicon pipeline, shared helpers, CLI routing, and metadata

### Changed
- Relaxed pbsim3 `pass_num` validation from `>= 2` to `>= 1` (ONT single-pass)
- Refactored PacBio amplicon pipeline to use shared extraction stages
- Metadata platform field changed from `PacBio-Amplicon` to `PacBio` with separate `assay_type`

---

## [0.40.0] - 2026-04-04

### Added
- PacBio amplicon read simulation using PBSIM3 template mode (`reads amplicon` command)
- PCR length bias model with exponential decay, calibrated to Madritsch et al. 2026 empirical data
- Deterministic and stochastic (Galton-Watson) PCR bias modes
- Preset profiles (`default`, `no_bias`) for PCR bias configuration
- Primer-based amplicon extraction from diploid VNTR references
- Shared primer binding site utility (refactored from snapshot validator)
- `amplicon_params` configuration section with primer sequences and PCR bias settings
- 69 new tests covering all amplicon simulation components

### Changed
- Updated GitHub Actions to Node.js 24 compatible versions (checkout v6, setup-python v6, setup-uv v7)
- Extended config schema to accept `"amplicon"` and `"pacbio"` as simulator types

---

## [0.19.0] - 2025-10-20

### Added
- MkDocs Material documentation system with GitHub Actions deployment
- Comprehensive user guides (simulation, toxic protein detection, SNaPshot validation)
- Auto-generated CLI documentation via mkdocs-click
- Dark/light mode toggle, instant search, mobile responsive design
- Professional documentation at https://berntpopp.github.io/MucOneUp/

---

## [0.15.0] - 2025-10-18

### Added
- In silico SNaPshot assay validation for MUC1 dupC mutation

---

## [0.14.0] - 2025-10-15

### Added
- Diploid split-simulation for ONT reads

### Fixed
- ONT read simulation bias in diploid references

---

## [0.13.0] - 2025-10-10

### Added
- Toxic protein detection algorithm
- ORF prediction with orfipy integration

---

## Earlier Versions

See [GitHub Releases](https://github.com/berntpopp/MucOneUp/releases) for earlier version history.

---

[Unreleased]: https://github.com/berntpopp/MucOneUp/compare/v0.40.0...HEAD
[0.40.0]: https://github.com/berntpopp/MucOneUp/releases/tag/v0.40.0
[0.19.0]: https://github.com/berntpopp/MucOneUp/releases/tag/v0.19.0
[0.15.0]: https://github.com/berntpopp/MucOneUp/releases/tag/v0.15.0
[0.14.0]: https://github.com/berntpopp/MucOneUp/releases/tag/v0.14.0
[0.13.0]: https://github.com/berntpopp/MucOneUp/releases/tag/v0.13.0
