# GENOMICON-Seq evaluation (issue #52, Phase 1)

Date: 2026-11-10. Desk review only: repo read at commit of 2026-11-03, nothing installed or run.

## Question

Should MucOneUp use GENOMICON-Seq to get realistic PCR duplicate families
(one source molecule giving many reads) in the Illumina path?

## Findings

### The tool

- Repo: https://github.com/Rounge-lab/GENOMICON-Seq (the Gabaldonlab URL in the
  issue is wrong). MIT license, compatible with MucOneUp's MIT.
- Maintained: latest release v1.2.4 (2026-11-03). Docker images v1.2 to v1.2.4
  from 2025-08 to 2026-11 (https://hub.docker.com/r/mimsto86/genomicon-seq).
  2 stars, 0 open issues: one lab maintains it.
- Publication: Stosic et al., *Sci Rep* 15, 23003 (2025),
  https://doi.org/10.1038/s41598-025-05267-8. The DOI in the issue
  (10.1093/bioinformatics/btac179) points to a different Bioinformatics 2022
  article, not to GENOMICON-Seq.
- Install: Docker only (README: "Currently, it is only available as a Docker
  container", about 1 GB image). No conda or PyPI. It is a Snakemake pipeline of R, C++,
  Python and shell scripts and needs 16 GB RAM or more.

### Fit with MucOneUp

| Need | GENOMICON-Seq (WES mode) | Fit |
|------|--------------------------|-----|
| Custom FASTA (2 × ~46 kb haplotypes) | Any FASTA or multi-FASTA, with per-sequence copy numbers | Yes |
| Custom probes/targets | `PROBES_FASTA` or `PROBES_BED` in input-FASTA coordinates, plus a targets BED | Yes, but we would have to project probe coordinates onto each haplotype |
| Capture model | Probe/fragment match length ≥ `MATCHING_LENGTH` | Coarser than our pBLAT + bind-fraction model |
| PCR duplicates | Per-fragment binomial growth over cycles with efficiency `1/(1+exp(k·(cycle-midpoint)))` (`scripts/wes/PCR_cpp.cpp`) | This is the idea #52 wants |
| PCR errors shared within a family | Not in WES mode: the README says the WES PCR runs "excluding the polymerase error rate". Only ampliseq mode has `POL_ERROR_RATE` | No. The issue's main selling point does not hold for capture |
| Sequencing errors | Modified InSilicoSeq 1.6.0 (HiSeq/NextSeq/NovaSeq/MiSeq KDE models) | Would replace ReSeq, or give a second error model to maintain |
| Outputs | FASTQ.gz plus mutation CSVs. No BAM, no per-read truth in our sidecar format | We would still have to align and rebuild source tracking ourselves |
| Runtime | Builds an SQL database per input sequence. Its reported runs are about 1 h for chr1 on 60 CPUs | Too heavy for 92 kb of input and for CI (16 GB, Docker) |

### What MucOneUp already does

- Illumina path: `pipeline.py` → `stages/fragment_preparation.py` →
  `fragment_simulation.simulate_fragments`. It is Wessim2-derived (no Wessim
  binary any more), and ReSeq applies the errors. Each read pair is an
  independent draw (`pick_on_match` → `get_insert_length` → `pick_fragment`),
  so duplicates only appear when two draws hit the same coordinates by chance.
  This matches the hypothesis of #51 (closed).
- `PCRBiasModel` (`read_simulator/pcr_bias.py`) models allele length bias for
  long-range amplicons (PacBio/ONT), with presets and stochastic mode. It
  splits coverage between alleles and does not model duplicates.
- `molecules.py` and `read_profiles.py` model long-read amplicon artefacts
  (strand, chimeras, smear, stutter). There is no duplicate or UMI concept in
  any pipeline, and there is no duplicate marking.

## Recommendation

**Do not integrate GENOMICON-Seq as a backend. Copying the duplicate-family idea
into `simulate_fragments` is the right approach, but only once a consumer needs it.**

Reasons:

1. The only feature we would gain is per-fragment copy numbers from a sigmoid
   PCR. That is about 6 lines of C++ upstream and about 15 lines of Python here.
2. It does not deliver what the issue wanted most: in capture/WES mode,
   duplicates do not share polymerase errors.
3. A backend would add Docker as a hard dependency, a second error model
   (InSilicoSeq next to ReSeq), probe coordinate projection, an alignment step,
   and lost per-read truth. It is also far slower than the current pipeline,
   which already fails the issue's "<20% runtime overhead" criterion.
4. No current consumer in this repo needs it: there is no duplicate marking or
   UMI handling anywhere in the pipelines. Build it when a benchmark depends on
   duplicate structure, for example consensus or duplicate-aware callers.

## Scoped follow-up (only when needed)

One small PR in `muc_one_up/read_simulator/fragment_simulation.py`:

- New optional config `read_simulation.pcr_duplicates`:
  `{unique_molecules: int, cycles: 12, k: 1.0, midpoint: null}`. When it is
  absent, behaviour and golden md5s stay exactly as they are.
- When set: draw `unique_molecules` fragments with the existing
  `pick_on_match`/`pick_fragment`. Give each a copy count by running
  `copies += binomial(copies, 1/(1+exp(k*(c-mid))))` for `c` in `1..cycles`
  (numpy is already a dependency; `pcr_bias.py` uses `np.random.default_rng`). Then draw `read_number` pairs from the molecules,
  weighted by copy count.
- ReSeq still adds independent sequencing errors per read, which is correct.
  Shared PCR errors stay out of scope; GENOMICON-Seq WES mode leaves them out too.
- Add a `molecule_id` column to the fragment-origins sidecar so duplicate
  families are part of the truth.
- Test: with a fixed seed, check that the duplicate rate rises as
  `unique_molecules` falls, that family sizes are skewed rather than all 1 or 2,
  and that the setting being absent leaves output unchanged.

Out of scope: the GENOMICON-Seq Docker backend, UMIs, polymerase errors, and
probe BED kits.
