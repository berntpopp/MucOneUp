# MucOneUp Codebase Audit — 2026-11-09

Base: `main` @ `1ab736d` (v0.48.0). Read-only audit. Five parallel reviewers covered: core generation, read-simulation pipelines, tool wrappers/subprocess utils, CLI, and helpers/packaging/tests/CI. The top findings were re-checked by hand against the source.

**What this repo does:** MucOneUp is a Python CLI that builds synthetic MUC1 VNTR diploid haplotypes, with chosen repeat counts, mutations and SNPs. It then simulates Illumina, ONT and PacBio reads from them, driving external tools (reseq, NanoSim, pbsim3, ccs, bwa, minimap2, samtools), and writes truth files for benchmarking variant callers.

**Load assumed:** one researcher on one machine. Runs can be large batches (series, many iterations) on big FASTA/BAM/FASTQ files against GRCh38. Nothing runs at the same time inside the code.

**Test status:** 1832 passed, 13 skipped, 0 failed (≈35 s). Coverage 92%. CI enforces `--cov-fail-under=80`, so the "30%" in `CLAUDE.md` is out of date.

**Main theme:** the code is well tested, but some **results are silently wrong**: wrong sequences, wrong counts, wrong truth files, broken read pairs. For a benchmarking tool this is the costliest kind of bug, because nothing crashes and nobody notices.

---

## Must fix

### 1. The same mutation target given twice is applied twice
`muc_one_up/cli/mutations.py:12-41`, applied in `muc_one_up/mutate.py:175-240`
- **What this is:** the code that turns `--mutation-targets` into a list of (haplotype, repeat) positions.
- **Problem:** duplicates are not removed. `--mutation-name dupC --mutation-targets 1,25 --mutation-targets 1,25` exits 0. Haplotype 1 then has *two* inserted Cs (22402 bp instead of 22401), but the structure file reports one mutation.
- **Fix:** end `parse_mutation_targets` with `return list(dict.fromkeys(mutation_positions))`. `MutationTarget` is a frozen dataclass, so it can be hashed. Add one test.
- **If we skip it:** the FASTA no longer matches the truth file, which corrupts any benchmark built from it.

### 2. Pipelines align to different reference genomes for the same config
`muc_one_up/read_simulation.py:188-213`, `pipeline.py:128`, `ont_pipeline.py:293-303`, `pacbio_pipeline.py:361-379`, `amplicon_pipeline.py:240`, `ont_amplicon_pipeline.py:183`, `ont_fragment_pipeline.py:131`, `pipeline_utils.py:98-125`
- **What this is:** each pipeline picks the human reference for alignment.
- **Problem:** each pipeline looks for it in a different place.
  - Illumina prefers `reference_genomes[asm].fasta_path`.
  - ONT prefers `read_simulation.human_reference`, then the other key.
  - PacBio, amplicon and ONT-fragments only read `read_simulation.human_reference`.
  - With the shipped `config.json`, Illumina aligns to `hg38.fa` and PacBio/ONT align to `GCA_…no_alt_analysis_set.fna`.
  - With only `reference_genomes` set, PacBio silently writes FASTQ only.
  - The first lookup step catches all exceptions at debug level, so a typo in the path silently switches to the other reference.
- **Fix:** resolve the reference once in `simulate_reads` in one fixed order (skip it when `skip_alignment` is set) and pass it into every pipeline. Log a failed lookup as a warning, not a debug message.
- **If we skip it:** BAMs from different platforms cannot be compared, and nothing tells you.

### 3. VNTR downsampling breaks read pairs (two places)
`muc_one_up/read_simulator/wrappers/samtools_coverage.py:94-160` and `muc_one_up/read_simulator/vntr_efficiency.py:155-189`
- **What this is:** reads in the VNTR are subsampled to mimic capture bias, and flank reads are kept.
- **Problem:** `samtools view -L` assigns each *read* to a set, not each *pair*.
  - Take a pair with R1 in the flank and R2 in the VNTR. R2 can be dropped while R1 always stays.
  - The result is orphan mates whose flags still say "mate mapped". That is about 8% of VNTR pairs at a 350 bp insert size.
  - `convert_bam_to_paired_fastq` then quietly drops both mates, so flank depth is also lower than intended.
- **Fix:** split by read name.
  1. Collect the VNTR read names: `samtools view -L bed | cut -f1 | sort -u`.
  2. Use `samtools view -N names.txt -U flank.bam` to split whole templates.
  3. Subsample the VNTR template BAM. `samtools -s` hashes the read name, so both mates share the same fate.
- **Notes:** make one shared helper for both call sites. No test covers this today; add one that asserts there are no orphan reads.
- **If we skip it:** the final BAMs fail Picard `ValidateSamFile`, and the VNTR-to-flank depth ratio is wrong.

### 4. Repeat counts are biased low
`muc_one_up/distribution.py:52`
- **What this is:** samples the repeat count from a normal distribution.
- **Problem:** `int(rng.gauss(...))` cuts off decimals toward zero instead of rounding.
  - The mean shifts down by about 0.5.
  - Measured with mean 2 and range 0–4: 0 came up 16% of the time, 4 only 2.2%.
- **Fix:** `val = round(_rng.gauss(...))`. This changes seeded outputs, so update the pinned test values and note it in the CHANGELOG.
- **If we skip it:** the simulated length distribution does not match the configured one.

### 5. Random SNP counts are rounded down, often to zero
`muc_one_up/snp_integrator.py:226`
- **What this is:** decides how many random SNPs to place.
- **Problem:** `int(len(valid_positions) * rate_per_base)` cuts off decimals. A 1.5 kb VNTR at 0.5 SNP/kb gives 0 SNPs on *every* iteration, when the expected average is 0.75.
- **Fix:** draw the count from the seeded RNG: `n = int(x) + (rng.random() < x % 1)`, or use a binomial draw.
- **If we skip it:** low-density SNP simulations silently contain no SNPs.

### 6. The empirical error model creates fake dupC-like homopolymers with no truth record
`muc_one_up/read_simulator/empirical_errors.py:157-176`
- **What this is:** injects sequencing errors and protects homopolymer runs from changes that are not tracked.
- **Problem:** deletions respect the protected runs, but insertions and mismatches *next to* a run do not. For `GGGTGCCCCCCCAGTGG` at 5% insertion and 5% mismatch, 4.6% of reads came out with C8 or longer while `hp_edits` stayed empty.
- **Fix:** when the next or previous base is protected, redraw any inserted or substituted base that equals the run's base.
- **If we skip it:** false dupC-like signals show up in benchmarks, and the truth file says they are not there.

### 7. The download helpers accept broken installs as finished
`helpers/download_references.py:283-317`, `helpers/install_references.py:123-141, 293-308`
- **What this is:** scripts that download and unpack the NanoSim model and the GRCh38 reference.
- **Problem:**
  - **NanoSim:** the target folder is created before extracting. After a failed download, the next run sees the folder and returns it as "already exists", and config.json then points at an empty model (reproduced). `requests.get` also has no timeout, and the tarball's MD5 is never checked, although `install_references_config.json` already lists it.
  - **GRCh38:** gzip extraction writes straight to the final `.fna`. A Ctrl-C or a full disk during the ~3 GB unpack leaves a truncated genome that later runs skip over and index.
- **Fix:**
  - **NanoSim:** delete `extract_path` in the `except` branch, add `timeout=60`, and compare the MD5.
  - **GRCh38:** extract to `path.with_suffix(".part")`, then `os.replace` it to the final name.
- **If we skip it:** confusing NanoSim or aligner errors far from their cause, or reads from the missing chromosomes that silently never align.

---

## Should fix

### 8. ONT-fragments: about 2.3× the coverage asked for, and the shorter allele gets more depth
`muc_one_up/read_simulator/ont_fragment_pipeline.py:72-83`, `molecules.py:324-329`
- **Problem:**
  - **Coverage:** `_n_reads` multiplies by `len(sources)`, so each haplotype gets the full coverage, and it uses the median read length where the mean belongs. `--coverage 10` gave 23× in total. NanoSim and PacBio treat `--coverage` as the total over both haplotypes.
  - **Allele balance:** the haplotype is picked uniformly regardless of its length. With 21 kb vs 26 kb alleles, the short allele got 1.19× the depth.
- **Fix:**
  - **Coverage:** drop the `* len(sources)` factor and use the mean length `median*exp(sigma**2/2)`.
  - **Allele balance:** draw the read length first, then pick the haplotype with `weights=[len(s)+length-1 for s in sources]`.
  - Both change seeded output.
- **If we skip it:** the same `--coverage` gives different depths per simulator, and there is built-in bias against the expanded allele.

### 9. Truth files report wrong read coordinates or origin
`parsers/ont_parser.py:50-60`, `parsers/illumina_parser.py:73-84`, `parsers/pacbio_parser.py:33-38`
- **What this is:** parses simulator read names into the per-read truth manifest.
- **Problem:**
  - **ONT:** `ref_end = ref_start + len(seq)` counts the unaligned head and tail of the read. It should use the aligned length, which is field 7 of the NanoSim read name.
  - **Illumina:** fragment *i* is matched to read *i* by position. If there are fewer reads than fragments, it silently falls back to `fragment_N` IDs, so one dropped fragment shifts every read after it.
  - **HiFi:** each subread overwrites the strand, so a CCS read takes the strand of the *last* subread.
- **Fix:**
  - **ONT:** `ref_start + int(match.group(7))`.
  - **Illumina:** raise when `len(read_names) != len(fragments)`.
  - **HiFi:** `zmw_map.setdefault(...)`.
  - Update `tests/test_parsers.py`.
- **If we skip it:** overlap flags for the VNTR, mutations and SNPs are wrong in the truth manifest.

### 10. CLI options that are silently ignored or overridden
`muc_one_up/cli/config.py:156-260`, `cli/mutations.py:125`, `cli/snps.py:52-70`
- **Problem:** each of these exits 0 and quietly does less than asked:
  - `--simulate-series` without a `--fixed-lengths` *range* writes one file instead of a series.
  - `--mutation-targets` without `--mutation-name` writes unmutated output.
  - `--snp-input-file` together with `--random-snps` drops the random SNPs, and `--random-snp-output-file` is never written.
  - `--fixed-lengths` together with `--input-structure` drops the lengths without a warning.
- **Fix:** raise `ValidationError` for each conflicting combination, in `setup_configuration` and `determine_simulation_mode`, where the other checks already live. Four `if` statements and four tests.
- **If we skip it:** samples labelled "mutated", "with SNPs" or "series" that are not.

### 11. Batch outputs overwrite each other
`cli/commands/analyze.py:117-120, 203-209`, `cli/commands/reads.py:56-68`, `read_simulator/output_config.py:64-67`, `cli/snps.py:72-80`
- **Problem:**
  - The automatic output name is just the input's file stem, and `--out-dir` defaults to `.`. So `analyze stats r1/x.fa r2/x.fa --out-dir out` writes one JSON that holds only r2's data. `reads` and `analyze orfs` behave the same way.
  - A relative `--random-snp-output-file` is written relative to the current directory, not to `--out-dir`.
- **Fix:**
  - Before each batch loop, raise when the set of resolved output bases has duplicates, and tell the user to pass `--out-base`.
  - Join a bare `--random-snp-output-file` name onto `out_dir`.
- **If we skip it:** results lost without any message.

### 12. The `--simulate-series` docs describe 5 runs, while the command really runs 81
`docs/guides/simulation.md:159-180, 761-765`, `docs/getting-started/quickstart.md:343-350`
- **Problem:**
  - The docs say `--fixed-lengths 40-80 --simulate-series 5` gives 5 samples. In the code, 5 is a *step size*, and the range is applied to both haplotypes and combined, so the command runs 9 × 9 = 81 iterations.
  - The documented `--out-base output/series` is now rejected.
  - The documented `--mutation-targets "1,25 2,30"` fails to parse.
- **Fix:** correct the docs: step semantics, the per-haplotype product, `--out-dir` instead of a path in `--out-base`, and repeated `--mutation-targets` flags.
- **If we skip it:** users start 16× more work than they expect, and the copy-paste examples fail.

### 13. bwa alignment against GRCh38 is hard-capped at 300 s
`read_simulator/wrappers/bwa_wrapper.py:24,59`, caller `stages/alignment.py:86`
- **Problem:** `timeout or 300` turns "no timeout" into 300 s. That budget has to cover loading the ~5 GB index plus all alignment. A cold disk cache or high coverage kills the run and can leave a partial BAM behind.
- **Fix:** default `timeout=None`, pass it through as is, and let `alignment.py` pass `rs_config.get("alignment_timeout")`.

### 14. The error output of the failing tool is lost in piped commands
`read_simulator/utils/common_utils.py:335, 403-414`
- **Problem:** upstream stages such as bwa and minimap2 write stderr to the terminal, and the raised error only shows the *last* stage's stderr. If the bwa index is missing, the error says "samtools sort: failed to read header", and the real bwa message is in neither the log nor the exception.
- **Fix:** give each upstream stage `stderr=tempfile.TemporaryFile()` (not PIPE, so it cannot deadlock), and read the failing stage's file into the error.

### 15. The minimap2 index check and the aligner disagree, so GRCh38 is re-indexed on every run
`bioinformatics/reference_validation.py:83-108` vs `wrappers/minimap2_wrapper.py:58-82`
- **Problem:** the validator accepts `ref.fa.mmi` and `{stem}_ont.mmi`, but the aligner only reuses `{ref}.{preset}.mmi`. A user with `ref.fa.mmi` sees no warning, yet every long-read run builds the index in memory again: about 12 GB RAM and several minutes.
- **Fix:** validate `Path(f"{ref}.{preset}.mmi")` for the preset actually in use, and fix the docstring at `minimap2_wrapper.py:155`.

### 16. Reproducibility gaps: unseeded runs, and a fingerprint that depends on the current directory
`simulate.py:181-186`, `provenance.py:112-136, 333`, `config_fingerprint.py:41-56`, `cli/analysis.py:66-71`
- **What this is:** provenance metadata that lets you check and rebuild a run.
- **Problem:**
  - Without `--seed`, runs use the global `random` and record `seed: null`, so a surprising result cannot be rebuilt.
  - `PATH_KEYS` leaves out `reseq_model`, `sample_bam_*`, `sample_target_bed`, `capture_bed`, `flank_fasta` and `reference_genomes.*.fasta_path`. Because `load_config` makes some paths absolute only when needed, the same config gives different fingerprints from different directories (`7286ca…` vs `1c4953…`).
  - The fingerprint is recomputed every iteration. That was 57% of the iteration time in an 81-run series.
- **Fix:**
  - When no seed is given, draw one with `random.SystemRandom().randrange(2**32)` and store it.
  - Add the missing keys to `PATH_KEYS`.
  - Compute the fingerprint once in `simulate` and pass it into `collect_provenance_metadata`.

### 17. Output options that work in some pipelines but not others
`amplicon_pipeline.py:170`, `ont_amplicon_pipeline.py:121`, `ont_pipeline.py:239`, `pipeline.py:155-168`
- **Problem:**
  - `keep_intermediate_files` is ignored in the amplicon, ONT-amplicon and ONT-split pipelines; one has `keep_intermediate=False` hard-coded. The #150 fix only covered PacBio WGS.
  - Illumina with `--no-align --track-read-source` writes no read manifest and leaves a `_fragment_origins.tsv` behind. ONT and PacBio do write a manifest without alignment.
- **Fix:**
  - When the option is on, create the work directory with `mkdtemp` under `output_dir`.
  - Always call `generate_read_manifest`; it already falls back to all origins when there is no BAM.

### 18. CI does not run some tests or use the lockfile, and releases ship large files
`.github/workflows/test.yml`, `pyproject.toml`, `docker/Dockerfile:62`, `data/testdata_releases/`, `tests/test_bam_to_fastq_conversion.py:34-45`, `tests/test_snapshot_validator_integration.py:25-63`
- **Problem:**
  - **Skipped tests:**
    - 12 of the 14 BAM-to-FASTQ regression tests skip in CI because samtools is not installed there.
    - 6 snapshot-validator tests *always* skip, because they need `output/dupC/`, which is gitignored and never generated.
  - **Lockfile and pins:**
    - CI installs with `uv pip install -e .` and the Dockerfile uses `pip3 install`, so `uv.lock` is used by neither.
    - mypy is unpinned.
    - `make ci-check` runs pytest without the 80% coverage threshold.
  - **Package and repo size:**
    - The sdist is **40 MB**, because it includes `.planning/`, the test-data tarballs, the BAMs in `data/` and the docs images.
    - The 17 MB test-data tarball is committed twice (once as `_latest`), and each new release adds about 17 MB to git history.
- **Fix:**
  - **CI tests:** add `apt-get install -y samtools` to the CI job, and give the snapshot tests a `tmp_path` fixture that runs `simulate --fixed-lengths 27 --mutation-name normal,dupC`.
  - **Lockfile:** use `uv sync --locked --extra dev` in CI, pin mypy, and add `--cov-fail-under=80` to `ci-check`.
  - **Package size:** add `[tool.hatch.build.targets.sdist] exclude = [".planning", "data/testdata_releases", "data/*.bam*", "docs/assets"]`.
  - **Tarballs:** move them to GitHub Releases, gitignore the folder, and update the link in `docs/reference/test-datasets.md:7`.

### 19. Duplicated code that has already drifted apart
- **reverse_complement exists three times, with three behaviours:**
  - `translate.py:146` turns `R` into `N`.
  - `read_simulator/utils/sequence_utils.py:32` keeps IUPAC codes (`R` to `Y`).
  - `helpers/download_references.py:145` does not complement IUPAC codes at all.
  - Keep `sequence_utils`. Check the `cli/analysis.py:308` caller before switching `translate.py` to it.
- **Two ORF pipelines with different start codons:** `translate.py:268-396` (`predict_orfs_in_haplotypes`, `run_orf_finder_in_memory`, `write_peptides_to_fasta`) is used only by tests. It treats ATG, TTG and CTG as starts and translates TTG/CTG as `L`, while the real CLI path uses ATG only. Delete it, and set `find_orfs_in_memory`'s default start codons to `["ATG"]`.
- **The "how much did a mutation change the length" rule is written twice:** `mutate.py:151-170` and `simulation_statistics.py:162-167`, with different flank checks. Both are wrong when 2 or more repeats are mutated and no stored lengths are given. Use `get_repeat_lengths` in both places and raise an error in the ambiguous case.
- **Samtools helpers overlap:** two `downsample_bam` functions (`wrappers/samtools_coverage.py:94` vs `utils/samtools.py:183`), `merge_bam_files` vs `merge_bams`, and `calculate_target_coverage` vs `calculate_mean_coverage`. Keep one of each.
- **Near-copy pipelines:** `amplicon_pipeline.py:165-285` and `ont_amplicon_pipeline.py:117-227` share the same flow. Factor out the shared body so fixes like #17 only need to be made once.

### 20. Dead code
Each item below was grepped across the whole tree, including tests, config and docs.
- **`muc_one_up/type_defs.py:139-155, 242-297`:** 17 type aliases and protocols that nothing uses, for example `ToolWrapper`, `ConfigLoader`, `SNPRecord`, `MutatedUnits`; `MutatedUnits` also has the wrong shape.
- **`common_utils.py`:** `fix_field` and `cleanup_files` (`:441, :526`) are unused.
- **`RunResult.__bool__` (`common_utils.py:50`):** inverted, so a successful result is falsy. Delete it.
- **Used only by tests:**
  - `sort_and_index_bam`, `validate_preset`, `get_tool_executable`, `check_bedtools_available`, `validate_bam_file`, `validate_bed_file`, `get_muc1_region_for_assembly`.
  - `ReadSourceTracker.from_simulation_results`, `PCRBiasModel.from_params`, `OutputConfig.ensure_dir`, `assemble_haplotype_from_chain`.
- **`read_simulator/constants.py`:** 8 unused constants (`:138, 153, 157, 165, 169, 188, 191, 198`).
- **`cli/_common.py:48-87`:** `validate_config_callback` returns its input on both branches. Delete it along with its `callback=` and `is_eager=`; `click.Path(exists=True)` already checks the file.
- **`simulate.py:184-186`:** `elif rng is None: rng = None` does nothing.
- **Tests:**
  - `tests/analyze_dupc_mutation.py` (190 lines) is never collected and runs against output that does not exist.
  - `tests/test_cli.py:493-552` has three "tests" that are a docstring plus `pass`.

---

## Nice to have (smaller items)

| # | Item | Where | Fix |
|---|------|-------|-----|
| a | `split_reads` can silently get R1 and R2 out of step on an incomplete file; it also uses gzip level 9 and no `with` block | `wrappers/reseq_wrapper.py:144-199` | Check that each R2 header matches its R1 and raise if not; `gzip.open(..., compresslevel=1)` inside `with` |
| b | VNTR BED start is 1-based where BED expects 0-based | `config.json:576`, `vntr_efficiency.py:19-24` | Use 155188486, as `samtools_coverage.py:107` already does |
| c | Sorting an already sorted merged BAM | `vntr_efficiency.py:261-281` | Merge straight to the output and drop the sort |
| d | Tool versions probed twice per run, which takes seconds under `mamba run` | `utils/tool_version.py:30` | `@functools.lru_cache` |
| e | `check_external_tools` breaks on multi-word commands such as `singularity exec … samtools` | `common_utils.py:505,514` | Build the command with `build_tool_command(...)` |
| f | Timeouts in capture and to-file modes do not kill grandchild processes | `common_utils.py:105-188` | Use `Popen(start_new_session=True)` and `killpg`, as streaming mode already does |
| g | `.replace(".fa","")` turns `x.fasta` into `xsta` | `wrappers/samtools_core.py:38` | `removesuffix` / `Path.stem` |
| h | `analyze stats` crashes with a traceback when `--out-dir` does not exist | `cli/commands/analyze.py:247` | `mkdir(parents=True, exist_ok=True)` |
| i | A bad `reads --out-base` exits 2 with a traceback | `read_simulator/output_config.py` | Call `validate_output_base` first, as `analyze` does |
| j | The config schema accepts `min_repeats: 5` and float mutation positions, which fail mid-batch | `config.py:243-281` | Use `"integer"` with `"minimum": 10` |
| k | The reference installer downloads ~900 MB again on every run | `helpers/install_references.py:237` | Skip when the file exists and its MD5 matches |
| l | Wrong `.dict` filename check, so GATK reruns on every install | `helpers/install_references.py:389-406` | `installed_path.with_suffix(".dict")` |
| m | `create_tarball` crashes on a relative path | `helpers/generate_test_data.py:741` | `str(tarball_path)` |
| n | `random_snp_density = 0.001` default is never used (Click passes None), so the docs claim is false | `cli/options.py:41` | Set a Click default or remove the claim |
| o | `pysam` is optional but not declared as an extra; without it the read manifest quietly includes unaligned reads | `stages/source_manifest.py` | Declare it as an extra, or warn when it is missing |
| p | `CLAUDE.md` says the coverage threshold is 30%; CI enforces 80% | `CLAUDE.md` | Change it to 80% |

---

## Suggested order of work

1. **Wrong results, cheap fixes (#1, #4, #5, #10):** about 10 lines plus tests. Note in the CHANGELOG that seeded outputs change.
2. **Read-pair integrity and truth accuracy (#3, #6, #9):** these decide whether the benchmark data can be trusted.
3. **One reference resolver (#2):** a single function, called from every pipeline.
4. **Installer robustness (#7):** write to a temp name and rename; check MD5s.
5. **CI hygiene (#18):** turns the 18 skipped tests into real coverage, so later fixes are actually checked.
6. **Delete dead code and duplicates (#19, #20)** in one cleanup PR.

**Verdict:** the code is in good shape structurally and well tested (92% coverage), but several bugs silently produce wrong data. Fix #1–#6 first.

**Lean:** about 700 fewer lines possible (type_defs ≈75, the ORF pipeline used only by tests ≈130, dead helpers and constants ≈200, `analyze_dupc_mutation.py` 190, no-op tests and callback ≈100). No dependencies to remove; `biopython` is listed twice, once in the Dockerfile.

**Not checked:**
- Real runs of external tools (reseq, pbsim3, NanoSim, ccs, GATK). Whether reseq keeps read order, pbsim3's strand convention and NanoSim's `-x` flag were judged from the code, not run.
- Only parts of `stages/alignment.py`, `pipeline_utils.py`, `generate_test_data.py:280-700` and `calibrate_read_profile.py`.
- `docs.yml`, `docker.yml` and the bodies of `conda/*.yml`.
- The per-base Python loop in `apply_errors` is slow at high coverage. It was left out because speeding it up would change seeded output.
