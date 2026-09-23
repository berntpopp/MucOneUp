# Design: realistic, truth-tracked long-read simulation

Date: 2026-09-24. Consumer: MucOneSpan benchmark
(`MucOneSpan/.planning/2026-09-23-realistic-benchmark-spec.md`, option 1).

## 1. Motivation (measured)

A MucOneSpan audit of MucOneUp 0.44.5, plus calibration against real R10.4.1
data (PRJEB92208: 9 ONT MUC1 amplicon libraries; 5 in-house ONT genomic
samples), found:

| Property | Real data | MucOneUp today |
| --- | --- | --- |
| Amplicon read strands | both | one strand (600/600 reads flag 16) |
| Read → allele truth in amplicon mode | needed | absent; `--track-read-source` rejected; read names collide across haplotypes (162/600 duplicated) |
| Genomic (NanoSim) read-source tracking | needed | broken from the CLI: stats JSON lacks `haplotypes`/`config` keys; NanoSim renames `haplotype_1` → `haplotype-1`, so the parser regex assigns every read to haplotype 1 |
| VNTR error rate (ONT amplicon) | 2.1% (mismatch/ins/del 0.70/0.67/0.77) | 3.75%, insertion-heavy (44/30/26; the PacBio RS II difference ratio is used) |
| C7 homopolymer read correctly, + / − strand | 0.52 / 0.89 | ~0.85, symmetric |
| Single-junction short PCR products ("smear") | 8–52% of spanning reads (median 24%) | none |
| Between-allele chimeras / concatemers | ~2.3% / ~2.4% | none |
| Short off-target reads | 6–77% of reads (median 30%, ~370 bp) | none |
| PCR length bias, ln(long/short) per repeat unit | −0.056 (R² 0.85) | default preset ≈ −0.037 |
| Runtime / RAM per case | — | 34–66 s / ~12 GB, dominated by rebuilding the full-hg38 minimap2 index |
| Genomic read lengths (targeted) | median 3.5–6.9 kb; 1.8–9.1% of reads span the VNTR | NanoSim `-med/-sd` only via `other_options`; `--min-read-length` silently overrides config |

## 2. Goals

1. **Per-read truth in every long-read mode**: haplotype, molecule id,
   strand, product kind, source interval, and injected edits. Unique read
   names.
2. **A molecule model** for amplicon and fragment templates: strand mixing,
   PCR smear, chimeras, concatemers, off-target short products, and a
   strand-aware homopolymer stutter table. It is applied before pbsim3, so
   pbsim's own errors add on top.
3. **Read profiles**: versioned JSON files bundling molecule model + pbsim
   parameters + PCR preset. Ship `ont_r10_sup_amplicon_v1`,
   `ont_r10_genomic_v1`, and `hifi_amplicon_v1`. Select with
   `--read-profile NAME|PATH`.
4. **Genomic ONT via pbsim3 fragments**: `reads ont --simulator
   pbsim3-fragments` samples fragments (log-normal lengths, uniform starts,
   optional extended flanks) through the same molecule model and pbsim3
   template path. NanoSim remains the default, and its source tracking is
   fixed.
5. **Speed**: a documented `--no-align` path; amplicon output is FASTQ with
   no hg38 index rebuild. When alignment is requested, reuse the `.mmi`.
6. **Backward compatibility**: without `--read-profile` (or a `read_model`
   config section), outputs are **byte-identical** to 0.44.5 for the same
   seed. Every new behaviour is opt-in.

Non-goals: raw-signal simulation; Illumina changes; changing haplotype
generation. A separate small fix records the provenance seed and warns on
non-strict mutation rewrites.

## 3. Design

### 3.1 Molecule model (`read_simulator/molecules.py`)

```
Molecule(id, hap, kind, strand, seq, src_start, src_end, edits: list[Edit])
kind ∈ {full, smear, chimera, concatemer, offtarget, fragment}
Edit(pos_in_source, kind="hp", base, true_len, new_len)
```

`build_amplicon_molecules(amplicons: list[str], counts: list[int],
model: MoleculeModel, rng) -> list[Molecule]`:

1. **Full molecules**: `counts[h]` per haplotype, from the existing PCR split.
2. **Smear**: a fraction `smear_rate` of full molecules becomes
   single-junction deletion products. The junction start follows
   `Beta(a, b)` × length (real median relative position ~0.2–0.3). The
   retained fraction is uniform over `[smear_min_keep, 0.95]`. The deletion
   spans from the junction start to the resumption point.
3. **Chimeras**: a fraction `chimera_rate` joins a hap-A prefix to a hap-B
   suffix at a homologous position (same offset in repeat units). This only
   happens when there are two haplotypes.
4. **Concatemers**: a fraction `concatemer_rate` joins two full molecules
   head-to-tail.
5. **Off-target**: add `offtarget_frac × total` short molecules from the
   flanks with log-normal length (median `offtarget_median_bp`).
6. **Strand**: each molecule is reverse-complemented with probability
   `1 - forward_frac`.
7. **Stutter**: for every homopolymer run ≥3 in the molecule's source
   orientation, draw Δ from `stutter[base+len|strand]` (a pmf over −3..+3)
   and apply it to the template. Each change is recorded as an Edit.
   "+ strand" means the molecule reads the VNTR motif1→motif9, matching the
   MucOneSpan calibration convention.

`build_fragment_molecules(haplotypes, n_reads, length_model, flank_seq,
model, rng)` samples `(hap ~ 1:1, start ~ U, length ~ LogNormal(median,
sigma))` over `flank_left + haplotype + flank_right`. It clips to bounds and
applies strand and stutter as above.

### 3.2 Truth-preserving pbsim3 run

Templates are written with ids `m{molecule_id:07d}`. pbsim3 `templ` writes
reads plus a MAF; the read → template id pairs are parsed from the MAF
(`*.maf.gz`). Reads are renamed to `{base}_h{hap}_m{molecule_id}` (unique)
and merged. `{base}_read_truth.tsv.gz` gets one row per read:
`read_id, hap, molecule, kind, strand, src_start, src_end, n_hp_edits,
hp_edits` (compact JSON list).

pbsim3 wrapper additions: `difference_ratio` (e.g. "39:24:36") and the
existing `accuracy_mean`, passed only when set, so legacy commands are
unchanged. (pbsim3 3.0.x has no `--accuracy-sd`; it was removed in #115.)

### 3.3 Read profiles

`muc_one_up/data/read_profiles/<name>.json`:

```json
{
  "schema_version": 1,
  "name": "ont_r10_sup_amplicon_v1",
  "provenance": {"source": "PRJEB92208 aggregate statistics (public)", "derived_by": "MucOneSpan realprofile", "date": "2026-09-23"},
  "platform": "ont",
  "pbsim": {"model_type": "qshmm", "model_file": "QSHMM-ONT-HQ.model", "accuracy_mean": 0.985, "difference_ratio": "39:24:36"},
  "pcr_bias": {"preset": "madritsch2025_r10"},
  "molecules": {"forward_frac": 0.5, "smear_rate": 0.24, "smear_junction_beta": [2.0, 5.0], "smear_min_keep": 0.15,
                "chimera_rate": 0.023, "concatemer_rate": 0.024, "offtarget_frac": 0.30, "offtarget_median_bp": 370},
  "stutter": {"C7|+": {"-2": 0.08, "-1": 0.22, "0": 0.64, "1": 0.06}, "...": {}}
}
```

- The stutter pmfs are **injected** tables. MucOneSpan fits them so that
  simulated *observed* homopolymer spectra match real data (deconvolving
  pbsim's own errors), and then publishes the fitted profile.
- Only aggregate statistics are shipped, never reads or sequences.
- `ont_r10_genomic_v1` is shipped only if the owner approves releasing the
  in-house aggregate statistics. Otherwise it ships with public-literature
  defaults and is marked `calibration: generic`.
- The profile SHA-256 is recorded in the simulation metadata.

### 3.4 PCR preset `madritsch2025_r10`

Real slope: ln(long/short) = −0.056 per 60-bp unit.
For the model `E = e_max·exp(−αL)` and `ratio = ((1+E1)/(1+E2))^cycles`
with e_max 0.95 and 25 cycles, solving the exact model numerically at
2.7 vs 4.7 kb gives **α ≈ 9.3e-5** (slope −0.052 to −0.060 over 2–6 kb).
A test asserts that the preset reproduces −0.056 ± 0.006 per unit
over 2–6 kb.

### 3.5 CLI and config

- `reads amplicon ... --read-profile NAME|PATH`: new. Also
  `--strand-mix/--no-strand-mix` as a convenience override.
- `reads ont ... --simulator {nanosim,pbsim3-fragments}` (default nanosim),
  `--read-profile`, `--read-length-median`, `--read-length-sigma`,
  `--n-reads`, `--flank-fasta PATH` (extra flank sequence around the
  haplotypes for genomic fragments).
- `--no-align`: skip alignment even when `human_reference` is configured.
- Truth manifest: written whenever a read profile is active, or when
  `--track-read-source` is given. It is no longer rejected for amplicon mode.
- Config: an optional `read_model` section mirrors the profile keys;
  `--read-profile` overrides it.

### 3.6 NanoSim source-tracking fixes

- `simulate --track-read-source` writes the `haplotypes` and `config` keys
  required by `ReadSourceTracker.from_companion_files`.
- The ONT read-name parser accepts `haplotype[-_](\d+)`.
- A regression test uses a real NanoSim-style name:
  `haplotype-2_1234_aligned_5_R_0_4000_0`.

## 4. Testing

- Unit (no binaries): molecule model statistics with a seeded RNG (smear
  rate, junction distribution, strand fraction, stutter pmf within binomial
  tolerance), MAF parsing, renaming, truth TSV, profile loading and
  validation, PCR preset slope, and name-parser fixes.
- Byte-identity regression: with no profile, template FASTA and pbsim
  argument lists equal 0.44.5 (golden files).
- Integration (marked; needs pbsim3): 200-molecule ONT amplicon with a
  profile → reads exist, names unique, truth rows == reads, observed strand
  fraction 0.5 ± 0.1.

## 5. Risks

- **pbsim3 MAF format differences across versions**: parse defensively and
  pin the tested pbsim3 version in docs.
- **Stutter injected before a context-free error model** only approximates
  R10. The MucOneSpan realism report quantifies the residual gap.
- **Profile defaults could silently shift users' results**: the opt-in
  design avoids that. The CHANGELOG states it clearly.
