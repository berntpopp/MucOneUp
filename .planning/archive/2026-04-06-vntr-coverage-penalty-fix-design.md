# VNTR Coverage & Penalty Factor Fix — Design Spec

**Date:** 2026-04-06
**Issues:** #85 (coverage calculation bug), #86 (re-derive penalty factor)
**Scope:** Bug fix + empirical recalibration + documentation

---

## Problem Statement

Two related issues affect the accuracy of simulated Illumina VNTR coverage:

1. **Bug (#85):** `calculate_mean_coverage()` in `utils/samtools.py` omits the `-a` flag when calling `samtools depth`, excluding zero-coverage positions from the mean. This inflates reported VNTR coverage in `vntr_efficiency_stats.json`. The main downsampling pipeline (`wrappers/samtools_coverage.py`) correctly uses `-a`.

2. **Stale calibration (#86):** The penalty factor 0.375 was derived from 277 real + 59 simulated samples (October 2024). Analysis of 1,043 CerKiD Berlin Twist v2 exomes yields a more accurate value of **0.39** (median implied penalty from observed VNTR/flanking ratio of 3.66 vs simulated base ratio of 9.32).

---

## Changes

### 1. Fix `-a` flag in `calculate_mean_coverage()`

**File:** `muc_one_up/read_simulator/utils/samtools.py`

**Current (line 241):**
```python
["samtools", "depth", "-b", str(region_bed), str(bam_file)]
```

**Fixed:**
```python
["samtools", "depth", "-a", "-b", str(region_bed), str(bam_file)]
```

This makes the function include all positions (including zero-depth) in the mean, consistent with `calculate_vntr_coverage()` in `wrappers/samtools_coverage.py`.

### 2. Add post-downsampling validation

**File:** `muc_one_up/read_simulator/stages/alignment.py`

After `downsample_bam()` completes (around line 252), re-calculate VNTR coverage on the output BAM and log the result:

```
Post-downsampling VNTR coverage: 152.3x (target: 150.0x, deviation: +1.5%)
```

If deviation exceeds 20%, emit a warning:

```
WARNING: Post-downsampling coverage 195.0x deviates >20% from target 150.0x
```

No failure — informational only. Applies to both `vntr` and `non_vntr` downsampling modes, re-using whichever coverage function (`calculate_vntr_coverage` or `calculate_target_coverage`) was used for the initial measurement.

### 3. Update penalty factor 0.375 → 0.39

**Files:**
- `muc_one_up/read_simulator/vntr_efficiency.py` — default parameter in `__init__`
- `config.json` — `vntr_capture_efficiency.penalty_factor`

The config field remains a single float. No structural change.

**Empirical basis (n=1,043 CerKiD Berlin Twist v2 exomes):**

| Metric | Value |
|--------|-------|
| Samples | 1,043 (all Twist v2) |
| VNTR mean coverage | 191.4 ± 96.3 (median 174.3) |
| Flanking mean coverage | 47.9 ± 10.3 (median 46.4) |
| VNTR/Flanking ratio | 4.03 ± 1.97 (median 3.66) |
| VNTR % uncovered | 10.3 ± 2.7 (median 10.1) |
| Simulated base ratio (no bias) | 9.32 |
| **Implied penalty factor** | **mean 0.43, median 0.39** |
| 95% CI | [0.20, 1.28] |

Median (0.39) chosen over mean (0.43) because the distribution is right-skewed (a few outlier samples with very high ratios).

### 4. Documentation update

**File:** `docs/guides/vntr-capture-efficiency.md`

Rewrite the Empirical Validation section:
- Replace "277 real + 59 simulated" with "1,043 CerKiD Berlin Twist v2 exomes"
- Update penalty factor from 0.375 to 0.39
- Add the cohort statistics table above
- Update coverage ratio tables (expected with penalty 0.39)
- Update interpretation guide ranges
- Add note: all samples in the calibration cohort use Twist Bioscience v2 capture panels
- Update version history

### 5. Tests

**Existing tests to update:**
- `tests/test_vntr_efficiency.py` — update hardcoded 0.375 references to 0.39
- `tests/read_simulator/test_alignment_stage.py` — update if penalty_factor is referenced
- `tests/read_simulator/test_pipeline_contracts.py` — update if penalty_factor is referenced

**New tests:**
- `tests/read_simulator/test_samtools_utils_coverage.py` — test `calculate_mean_coverage()` with a mock that verifies `-a` is in the command
- `tests/read_simulator/test_alignment_stage.py` — test post-downsampling validation logging (mock `calculate_vntr_coverage` to return a value, verify log output)

---

## Out of Scope

- Per-kit penalty profiles (not needed — single Twist v2 value covers all current kits)
- Config schema changes (penalty_factor remains a float)
- Fixing the 23.5% vs 10.3% uncovered bases discrepancy (separate issue related to simulation fidelity, not the penalty model)
- Re-running simulations to validate the new factor end-to-end

---

## Risk

Low. The `-a` flag fix only affects reported statistics, not the downsampling decision (which already uses `-a`). The penalty factor change (0.375 → 0.39) is a small adjustment within the original confidence interval.
