# Changelog

All notable changes to MucOneUp will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

## [0.48.0] - 2026-10-19

### Fixed
- **Amplicon simulation pipelines (#147)**:
  - Prevent duplicate read IDs in diploid PacBio/ONT amplicon runs when `--seed` is set by deriving distinct per-haplotype seeds (`seed + 1`, `seed + 2`) and prefixing diploid pbsim3/empirical amplicon reads (`hap1_`, `hap2_`).
  - Reject polyploid FASTAs (>2 records) in diploid amplicon pipelines instead of silently dropping haplotype 3+; reject `--coverage` < 2 in diploid amplicon runs before extraction.
  - Propagate original `ExternalToolError` exceptions without double-wrapping, honour `OutputConfig.out_dir` / `out_base`, and clean up all `pbsim3` multi-template sidecars (`*.fq.gz`, `*.maf.gz`, `*.ref`).
- **External tool wrappers & VNTR efficiency (#148)**:
  - Pass `--id-prefix` to `pbsim3` to avoid read-ID collisions across diploid splits; remove duplicate coverage logging; validate `min_Bind` (`0.0 <= min_Bind <= 1.0`) and `min_fragment_size <= max_fragment_size` in `prepare_fragments`; validate `threads >= 1` in `reseq_wrapper`.
  - Stream FASTA split output in `nanosim_wrapper.split_fasta_for_diploid` via `SeqIO.parse`; surface actual `samtools` stderr on empty BAM/FASTQ outputs instead of masking with generic errors; remove redundant `samtools index` subprocess invocation after `samtools sort --write-index`.
  - Include temp directory in `VNTREfficiencyModel` cache path to prevent multi-user `/tmp` collisions and add retry backoff when removing transient temp directories.
- **CLI read simulation options & pipeline routing (#149)**:
  - Propagate `reads pacbio --coverage` into `config["pacbio_params"]["coverage"]`; inherit root `config["seed"]` when `--seed` and `read_simulation.seed` are omitted; accept float `--coverage` across `reads illumina`, `reads ont`, and `reads pacbio`.
  - Wire `--threads` (`-t`) and `--no-align` on `reads amplicon` through `pacbio_params` and `read_simulation`; remove dead legacy `simulate_reads` dispatcher and simplify `simulate_reads_pipeline` signature.
- **Long-read pipelines, parsers & truth tracking (#127, #150)**:
  - Honour `keep_intermediate_files` in `pacbio_pipeline` and `ont_pipeline`; remove CCS report/metrics sidecars (`*ccs_report.txt`, `*ccs_metrics.json`, `*.zmw_metrics.json.gz`) during intermediate cleanup.
  - Reconcile PacBio HiFi and ONT source-tracking manifests with the final surviving FASTQ reads (preventing inflated manifest rows after CCS filtering or NanoSim split merging) and support post-filtering in `ReadSourceTracker.write_read_manifest`.
  - Stream FASTQ reading/writing in `renumber_fastq_reads` and `parse_pacbio_reads` without buffering whole files in memory; reject polyploid (>2 record) FASTAs in `split_diploid_fasta`; guard empty reference dicts in `annotate_reads_from_umi`.
- **SNP integration & coordinate translation (#151)**:
  - Eliminate unconditional `set(range(seq_length))` memory spike in `generate_random_snps` by passing `Sequence[int]` ranges directly (`_SampledPositions` exclusions + binary search over VNTR bounds).
  - Preserve `HaplotypeResult._repeat_lengths` through `integrate_snps_unified`.
  - Translate `--snp-input-file` coordinates across indel mutations in single-mutation mode and support multi-repeat indels in fallback repeat-length calculation via shared `_translate_snps_to_mutated`; accept haplotype indices >= 3 in `parse_snp_file` for polyploid simulations.
- **Haplotype simulation, assembly & configuration (#152)**:
  - Align `simulate_single_haplotype` dead-end retry semantics with `from_predefined_chains`; validate `mutation_targets` syntax (`haplotype_index >= 1`, `repeat_index >= 1`) and `fixed_lengths` positivity before simulation.
  - Reject path separators and traversal in `--out-base` (`validate_out_base`); validate `--reference-assembly` against configured reference assemblies before simulation; replace fragile `sys.path.insert` hacks in `tests/conftest.py`.
- **ORF analysis, translation & SNaPshot validation (#153)**:
  - Use `orfipy_core.orfs()` in-memory in `run_orf_finder` instead of spawning `orfipy` CLI subprocesses and temp directories; forward loaded config into `analyze orfs`; strip trailing `*` stop codons in `translate_dna` before regex matching.
  - Propagate non-zero exit codes when `analyze orfs` or `analyze stats` fails across all inputs; clean up `SnapshotValidator.validate()` temporary directories using `tempfile.TemporaryDirectory`.
- **Illumina fragment simulation & empirical error model (#154)**:
  - Slice `Wessim2` fragment start/end tendency arrays by actual nucleotide coordinates (`frag_start:frag_end`) rather than repeat-symbol indices (`range(len(hr.chain))`).
  - Raise `ValueError` on out-of-range nucleotide offsets in `EmpiricalErrorModel.sample_indel`; surface `OSError` failures during `cleanup_intermediates`; reduce unit test subprocess overhead.
- **Reference helper security hardening (#155)**:
  - Validate download URLs (`https://` scheme and trusted domain allowlist), enforce SHA-256 verification, stream downloads with a 500 MB size cap, and apply `filter="data"` on `tarfile.extractall()` in `helpers/download_references.py`.
  - Add `--preset` (`map-ont`, `map-hifi`, `asm20`) to `helpers/install_references.py` so pre-built minimap2 `.mmi` indexes match runtime preset filenames (`{stem}_{preset}.mmi`).
- **Packaging, CI & documentation (#156)**:
  - Change default `reseq` command in `config.json` from a hardcoded conda binary path to `"reseq"`; fix `Makefile` targets (`format`, `security`, `docs-serve`) to use `uv run` and `muc_one_up`; run `pytest -m "not slow"` in the CI `coverage` job so codecov uploads complete on schedule; update `.trivyignore` and installation/quickstart documentation.

## [0.47.0] - 2026-10-18

### Fixed
- **Mutation & sequence assembly (#134, #135)**:
  - Preserve prior mutations and live per-repeat lengths across multiple targets and sequential `apply_mutations` calls via `HaplotypeResult._repeat_lengths`; update repeat length after each step in multi-operation mutations; replace only the target repeat slice on non-strict symbol conversion.
  - Bound normal repeat-count rejection sampling (`1,000,000` draws) and reject `min_rep > max_rep`; validate `fixed_lengths` count in `simulate_diploid`; assemble haplotype sequences in $O(N)$.
- **ORF & toxic protein detection (#136, #137)**:
  - Clamp `composition_similarity` to `[0.0, 1.0]`; enforce `min_nt_length = max(min_len, orf_min_aa * 3)` in `run_orf_finder`; strip constant flanks per ORF record so multi-ORF FASTAs in different frames or with SNP-modified flanks score identically regardless of record order.
  - Support configurable `snapshot.primary_primer` and verify expected fluorophores in `SnapshotValidator`.
- **VNTR boundaries, repeat lengths & source tracking (#126, #138)**:
  - Propagate exact mutated repeat lengths through `generate_simulation_statistics`, `write_simulation_statistics`, `ReadSourceTracker.from_simulation_results`, and `companion_stats_fields`; unify VNTR boundary calculation (`get_haplotype_vntr_bounds`) across `simulation_statistics` and `snp_integrator`, handling early-`9` termination without right flank and 0-bp repeat deletions.
- **Configuration, CLI RNG threading & dual-mode SNPs (#139, #140)**:
  - Enforce single-CLI-entry config loading boundary; resolve relative NanoSim `training_data_path` prefixes; expand `CONFIG_SCHEMA` with `read_simulation`, `nanosim_params`, `pacbio_params`, `snapshot.primary_primer`, and `validation.mutant_extension`.
  - Thread seeded `rng` across CLI simulation, mutation, and SNP integration (honoring both `--seed` and `config["seed"]`); implement exact per-repeat coordinate translation for dual-mode SNPs across indels and enforce 0-based SNP input coordinates (with 1-based haplotype indices).
- **Read simulator pipelines, truth tracking & wrappers (#129, #142, #143, #144, #145)**:
  - Guard negative start coordinates in `amplicon_Clip`, preserve diploid amplicon templates when only one haplotype extracts, make `StutterTable`, `MoleculeModel`, and `ReadProfile` hashable, and validate 2-record diploid FASTAs.
  - Compute minimum fragment probe overlap via `math.ceil(plen * bind)`, guard invalid `bind`/insert bounds, deduplicate `reverse_complement` and IUPAC primer matching (including ambiguous template bases) into `sequence_utils.py`, and fix `EmpiricalErrorModel` weighted sampling boundary.
  - Split VNTR and non-VNTR reads in a single `samtools view -L -U` pass, honor `seed` across `VNTREfficiencyModel` and Stage 10/11 downsampling, and stream `samtools depth` in `check_vntr_coverage`.
  - Route ONT alignment through `minimap2_wrapper.align_reads_with_minimap2(..., exclude_flags=4)`, clamp extreme diploid PCR bias counts to `[1, total - 1]`, and stream FASTA parsing in `amplicon_extractor` and `diploid_handler`.
  - Escalate subprocess and pipeline timeouts (`SIGTERM` -> `SIGKILL`) across the entire process group (`os.killpg`), preserve upstream pipeline stage error messages, and validate tool parameters and non-empty outputs across `bwa`, `ccs`, `pbsim3`, `ucsc`, and `samtools` wrappers.

### Removed
- Deleted dead modules and duplicate helpers (#134, #135, #137, #141, #142, #144, #145): `muc_one_up/validation.py`, `muc_one_up/bioinformatics/validation.py`, `muc_one_up/probabilities.py`, `muc_one_up/cli/source_tracking_setup.py`, `rebuild_haplotype_sequence`, `estimate_hifi_yield`, `_count_fastq_reads`, `align_ont_reads_with_minimap2`, and duplicate `reverse_complement`/`comp` implementations.

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
