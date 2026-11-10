# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed
- A read profile is now loaded once per run and checked against the
  `read_model.sha256` recorded when it was applied; a profile edited on disk
  mid-run raises an error instead of being used silently (#128).
- A `read_model` in a config passed to `simulate_reads()` (library API) now
  gets the same platform check and `config_overrides` as the CLI (#128).
- `reads ont --simulator pbsim3-fragments` rejects a profile with amplicon-only
  molecule rates (smear, chimera, concatemer, off-target), which fragment mode
  would ignore (#128).

## [0.50.1] - 2026-11-10

**Seeded outputs change** for haplotypes that contain L or R repeats.

### Fixed
- Repeat units L, R and S now use the sequences from Vrbacká et al. 2025
  (pristanna/muc1repeats). They came from the VNtyper motif set and each
  differed by 1 bp (#201).

### Changed
- `dupA` (60dupA) is also allowed on repeat C, where Vrbacká et al. 2025
  report it (C:60dupA). Random `dupA` targets can now land on C (#201).

### Tests
- The real-tool e2e golden checksums are re-measured for the new seed-42
  haplotype (#179, #201). Reads from the same FASTA are byte-identical to
  0.49.0. The read-tracking e2e test sets its molecule count in
  `ont_amplicon_params` (#181).

## [0.50.0] - 2026-11-10

Resolves the 2026-11-09 follow-up audit (issues #159–#201).

**Seeded outputs change.**
- The VNTR end block is now drawn from the transition table: V always
  precedes 6/6p, and 6p can skip to 8 (#179).
- Amplicon and fragment long-read modes default to 30 molecules, using the
  platform's own `*_params.coverage` instead of the Illumina target of 150
  (#181).
- The Vrbacká variants `insG_pos58`, `insG_pos54`, `insA_pos54` and
  `ins16bp` now insert at the published position. `insG_pos54` and
  `insA_pos54` target the new units `bJ` and `aH` (#190).
- The benign controls `delCCC_benign` and `del3bp_pos30_benign` delete
  3 bp, not 4 (#159).

### Breaking
- Duplicate mutation targets are now an error. The CLI used to drop them
  silently, and the Python API applied them twice (#166).
- Config validation is stricter: only `distribution: "normal"` is accepted
  (`"uniform"` never worked; it gave fixed lengths), mutation positions
  must be ≥ 1, transition weights must be > 0, and transition targets must
  be known repeats or `END` (#180).
- `--coverage` must be > 0, and `--threads`, `--read-number`, `--n-reads`
  and `--min-read-length` must be ≥ 1 (#169).
- ONT without a human reference returns FASTQ instead of aligning to the
  simulated haplotypes (#183).

### Fixed
- The benign deletion controls were frameshifts and were flagged as toxic
  (#159).
- hg19 runs use the hg19 VNTR region for capture bias instead of the
  hard-coded hg38 region, which dropped the reads (#160).
- Dual mode with `--input-structure` reuses the file's targets and writes
  correct structure files (#161).
- Ctrl-C and SIGTERM stop the external tools instead of leaving them
  running (#162). A failed pipe reports the stage that really failed (#165).
- A dead end in the repeat chain, or an end block that cannot reach 9,
  raises `SimulationError` instead of returning a truncated haplotype
  (#167, #179).
- Kept ONT haplotype FASTQs are prefixed with the sample name, so samples
  no longer overwrite each other (#168).
- Docker: NanoSim models load again (scikit-learn 0.23.2), pbsim3 model
  paths resolve, and the image installs from `uv.lock` (#175, #177, #191).
- A missing explicit `human_reference` fails before simulation, and PacBio
  no longer deletes the HiFi FASTQ when alignment fails (#176).
- An unknown ReSeq version is treated as ReSeq2, so input is not chunked
  (#182). The metadata `Seed` row records the seed actually used (#185).
- `snapshot-validate` no longer crashes on PCR failure, and 0 products are
  reported as a PCR failure instead of normal (#184).
- pbsim3 accepts `accuracy_mean` < 0.75 and validates `length_sd` (#186).
- The Vrbacká citation DOI is corrected (#178).
- PacBio and ONT amplicon defaults agree (threads 8, a clear error for a
  missing PacBio `model_type`), and `analyze stats` reads `.gz` (#187).

### Performance
- Read-source tracking streams its manifest: peak memory at 1M read pairs
  drops from 743 MB to 46 MB, with byte-identical output (#189).
- The homopolymer scan is about 4x faster (#188).

### CI and tests
- The e2e tests detect the real `simulator.py` and `pbsim` binaries, and
  all 7 pass with the real tools (#163, #192).
- The Docker image is tested with a real `simulate` + `reads` run (#163).
- Version tags on the image are set only on releases (#164).
- The scikit-learn CVEs that NanoSim cannot avoid are documented in
  `.trivyignore`, with expiry dates (#175).

## [0.49.0] - 2026-11-09

Resolves the 2026-11-09 codebase audit
(`.planning/reports/2026-11-09-codebase-audit.md`) in six workstreams.

**Seeded outputs change.** Repeat counts are now rounded instead of
truncated. Random SNP counts use stochastic rounding, so low densities no
longer give 0 SNPs every time. ONT `pbsim3-fragments` now treats `--coverage`
as the total over both haplotypes (it used to give about 2.3x the request)
and weights haplotypes by length. The empirical error model no longer
lengthens protected homopolymers without recording it. Unseeded runs now
draw a seed and record it in provenance.

### Fixed
- Illumina simulation works again with real tools:
  - Gapped pblat matches no longer make insert-length sampling impossible.
  - The conda env, Docker image and `config.json` ran the legacy bioconda
    ReSeq 1.1 (since #156), which deadlocks above 10,000 fragments. They
    now use bioconda `reseq2=2.0.3`. If `tools.reseq` still points at
    ReSeq 1.x, `seqToIllumina` input is fed in batches of at most 10,000
    fragments. The old 600 s cap on that step is gone. Existing `env_wessim`
    environments need `mamba install -n env_wessim -c bioconda reseq2=2.0.3`.
- All pipelines now resolve the alignment reference once, in one order:
  `human_reference_<asm>`, then `human_reference`, then
  `reference_genomes.<asm>.fasta_path`.
- VNTR downsampling and the VNTR capture model now keep read pairs together
  (no orphan mates). The VNTR BED is now correctly 0-based.
- Truth manifests:
  - ONT read end coordinates come from the aligned length.
  - Illumina reads are matched to fragments by name.
  - A HiFi read takes the strand of its first subread.
- CLI:
  - A duplicate `--mutation-targets` is no longer applied twice.
  - Conflicting or ignored option combinations now exit 1.
  - Batch output-name collisions are rejected.
  - `--random-snp-output-file` is written under `--out-dir`.
  - `reads --track-read-source` no longer crashes when `--out-dir` does not
    exist yet.
- `keep_intermediate_files` now works in the amplicon and ONT pipelines.
  Illumina `--no-align --track-read-source` now writes a manifest.
- Config:
  - The fingerprint no longer depends on the working directory and is
    computed once per run.
  - The schema rejects repeat counts below 10 and non-integer mutation
    positions.
- Subprocesses:
  - Timeouts now kill the whole process group.
  - The failing pipeline stage's stderr is kept.
  - The bwa and reseq timeouts default to none; the new
    `read_simulation.alignment_timeout` sets one.
- The minimap2 index check now matches the index the aligner actually reuses.
- Reference helpers:
  - Downloads and extractions are atomic and MD5-checked.
  - Existing verified downloads are skipped.
  - The GATK `.dict` file name is correct.

### Changed
- CI installs from `uv.lock`, installs samtools, and pins mypy and setup-uv.
- The sdist shrank from 41 MB to 1.6 MB.
- The snapshot integration tests now run on simulated data instead of
  skipping.
- pysam is a new optional `tracking` extra.
- Test-data tarballs moved from git to GitHub Releases (`testdata-v0.25.0`);
  `data/testdata_releases/` is now gitignored. `make publish-testdata` uploads a
  new dataset with a `.sha256` file, and the generator no longer writes an
  in-repo `_latest/` copy. Old blobs remain in git history.

### Removed
- About 1,800 lines of dead or duplicated code: unused type aliases, a
  test-only ORF pipeline, duplicate samtools and `reverse_complement`
  helpers, unused constants and helpers, and no-op tests.

## [0.48.0] - 2026-10-19

Resolves 11 codebase audit issues (#127 and #147-#156) across PacBio/ONT
amplicon and long-read simulation pipelines, external tool wrappers, CLI read
simulation options, SNP integration and coordinate translation, haplotype
assembly and config validation, ORF and SNaPshot analysis, Illumina fragment
simulation, reference helper security hardening, and packaging/CI/documentation.
See `docs/about/changelog.md`.

## [0.47.0] - 2026-10-18

Resolves 14 codebase audit issues (#126, #129, #134-#145) across haplotype
mutation/assembly, VNTR boundary and repeat-length tracking, dual-mode SNP
coordinate translation, read simulation pipelines, external tool wrappers,
subprocess process-group timeout cleanup, and dead-code removal (#146). See
`docs/about/changelog.md`.

## [0.46.0] - 2026-09-24

Homopolymer runs without a stutter-table entry (e.g. the dupC C8 run) are no
longer simulated error-free by the ONT read profiles (#132). **Breaking change
for custom read profiles with an `errors` section:** they now need
`molecules.stutter_fallback`. See `docs/about/changelog.md`.

## [0.45.0] - 2026-09-24

Realistic, truth-tracked long-read simulation (read profiles, molecule model,
empirical error channel, per-read truth, fragment-based genomic ONT) and fixes
#97, #102, #104-#124. The maintained changelog is `docs/about/changelog.md`.

## [0.28.1] - 2026-03-16

### Fixed

- **Add biopython to core dependencies** — `biopython` was only listed in optional `dev` and `docs` dependency groups, but is imported at module level in `read_simulator/utils/reference_utils.py` and used in CLI commands, causing `ModuleNotFoundError: No module named 'Bio'` on fresh installs

## [0.20.0] - 2025-10-23

### Added
- Tool version tracking and metadata improvements

### Changed
- Update gitignore configuration

### Fixed
- Configure Trivy to ignore unfixable base image vulnerabilities
- Remove bulk dismissal script (no longer needed)

---

*For complete version history, see git log.*
