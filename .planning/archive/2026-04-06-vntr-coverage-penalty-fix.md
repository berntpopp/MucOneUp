# VNTR Coverage & Penalty Factor Fix — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix the `-a` flag bug in coverage reporting, add post-downsampling validation, update the penalty factor from 0.375 to 0.39, and rewrite the documentation with 1,043-sample cohort data.

**Architecture:** Four independent changes touching different layers: utils fix (samtools), pipeline enhancement (alignment stage), constant update (vntr_efficiency + config), and documentation rewrite. All existing tests updated to reflect 0.39.

**Tech Stack:** Python 3.10+, samtools, pytest, ruff, mypy

---

### Task 1: Fix `-a` flag in `calculate_mean_coverage()`

**Files:**
- Modify: `muc_one_up/read_simulator/utils/samtools.py:241`
- Test: `tests/test_samtools.py`

- [ ] **Step 1: Write failing test for `-a` flag**

Add a test that verifies the `-a` flag is passed to `samtools depth`:

```python
# In tests/test_samtools.py, add to the appropriate test class:

class TestCalculateMeanCoverage:
    """Tests for calculate_mean_coverage."""

    def test_includes_all_positions_flag(self, mocker, tmp_path):
        """Verify samtools depth is called with -a to include zero-coverage positions."""
        bam_file = tmp_path / "input.bam"
        bam_file.touch()
        region_bed = tmp_path / "region.bed"
        region_bed.write_text("chr1\t100\t200\n")

        mock_result = mocker.MagicMock()
        mock_result.stdout = "chr1\t100\t10\nchr1\t101\t0\nchr1\t102\t5\n"
        mock_run = mocker.patch(
            "muc_one_up.read_simulator.utils.samtools.run_command",
            return_value=mock_result,
        )

        from muc_one_up.read_simulator.utils.samtools import calculate_mean_coverage

        result = calculate_mean_coverage(bam_file, region_bed)

        # Verify -a flag is present in the command
        cmd = mock_run.call_args[0][0]
        assert "-a" in cmd, "samtools depth must be called with -a to include zero-coverage positions"

        # Verify mean includes the zero-coverage position: (10 + 0 + 5) / 3 = 5.0
        assert result == 5.0
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `python -m pytest tests/test_samtools.py::TestCalculateMeanCoverage::test_includes_all_positions_flag -v`
Expected: FAIL — `-a` not in the command list

- [ ] **Step 3: Add `-a` flag to `calculate_mean_coverage()`**

In `muc_one_up/read_simulator/utils/samtools.py`, change line 241 from:

```python
            ["samtools", "depth", "-b", str(region_bed), str(bam_file)],
```

to:

```python
            ["samtools", "depth", "-a", "-b", str(region_bed), str(bam_file)],
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `python -m pytest tests/test_samtools.py::TestCalculateMeanCoverage -v`
Expected: PASS

- [ ] **Step 5: Run full test suite to check for regressions**

Run: `python -m pytest tests/test_samtools.py -v --tb=short`
Expected: All tests pass

- [ ] **Step 6: Commit**

```bash
git add muc_one_up/read_simulator/utils/samtools.py tests/test_samtools.py
git commit -m "fix: add -a flag to calculate_mean_coverage() for zero-coverage positions (#85)"
```

---

### Task 2: Add post-downsampling validation logging

**Files:**
- Modify: `muc_one_up/read_simulator/stages/alignment.py:262`
- Test: `tests/read_simulator/test_alignment_stage.py`

- [ ] **Step 1: Write failing test for post-downsampling validation**

Add a test that verifies coverage is re-measured after downsampling and logged:

```python
# In tests/read_simulator/test_alignment_stage.py, add to TestAlignAndRefine:

    def test_post_downsampling_validation_logged(self, mocker, tmp_path, tools, assembly_ctx):
        """After downsampling, VNTR coverage is re-measured and logged."""
        mocker.patch(f"{MODULE}.align_reads")

        # Return coverages: first call = pre-downsample (200x), second call = post-downsample (155x)
        mock_calc_vntr = mocker.patch(
            f"{MODULE}.calculate_vntr_coverage",
            side_effect=[(200.0, "/tmp/depth1.txt"), (155.0, "/tmp/depth2.txt")],
        )
        mocker.patch(f"{MODULE}.downsample_bam")

        rs_config = {
            "threads": 4,
            "coverage": 150,
            "downsample_mode": "vntr",
            "downsample_seed": 42,
            "vntr_capture_efficiency": {"enabled": False},
        }
        assembly_ctx_with_vntr = mocker.MagicMock()
        assembly_ctx_with_vntr.vntr_region = "chr1:155188487-155192239"
        assembly_ctx_with_vntr.assembly_name = "hg38"

        result = align_and_refine(
            tools=tools,
            rs_config=rs_config,
            r1=str(tmp_path / "r1.fastq.gz"),
            r2=str(tmp_path / "r2.fastq.gz"),
            human_ref="/ref/hg38.fa",
            output_dir=tmp_path,
            output_base="test_sample",
            assembly_ctx=assembly_ctx_with_vntr,
        )

        # calculate_vntr_coverage called twice: pre-downsample + post-validation
        assert mock_calc_vntr.call_count == 2

    def test_post_downsampling_validation_warns_on_deviation(
        self, mocker, tmp_path, tools, assembly_ctx, caplog
    ):
        """Warning logged when post-downsampling coverage deviates >20% from target."""
        import logging

        mocker.patch(f"{MODULE}.align_reads")
        # Pre-downsample: 500x, Post-downsample: 250x (67% deviation from target 150)
        mocker.patch(
            f"{MODULE}.calculate_vntr_coverage",
            side_effect=[(500.0, "/tmp/d1.txt"), (250.0, "/tmp/d2.txt")],
        )
        mocker.patch(f"{MODULE}.downsample_bam")

        rs_config = {
            "threads": 4,
            "coverage": 150,
            "downsample_mode": "vntr",
            "downsample_seed": 42,
            "vntr_capture_efficiency": {"enabled": False},
        }
        assembly_ctx_with_vntr = mocker.MagicMock()
        assembly_ctx_with_vntr.vntr_region = "chr1:155188487-155192239"
        assembly_ctx_with_vntr.assembly_name = "hg38"

        with caplog.at_level(logging.WARNING):
            align_and_refine(
                tools=tools,
                rs_config=rs_config,
                r1=str(tmp_path / "r1.fastq.gz"),
                r2=str(tmp_path / "r2.fastq.gz"),
                human_ref="/ref/hg38.fa",
                output_dir=tmp_path,
                output_base="test_sample",
                assembly_ctx=assembly_ctx_with_vntr,
            )

        assert any("deviates" in record.message for record in caplog.records)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/read_simulator/test_alignment_stage.py::TestAlignAndRefine::test_post_downsampling_validation_logged tests/read_simulator/test_alignment_stage.py::TestAlignAndRefine::test_post_downsampling_validation_warns_on_deviation -v`
Expected: FAIL — `calculate_vntr_coverage` called only once

- [ ] **Step 3: Add post-downsampling validation to `align_and_refine()`**

In `muc_one_up/read_simulator/stages/alignment.py`, after line 262 (`current_bam = downsampled_bam`), add:

```python
            # Post-downsampling validation: re-measure and log actual coverage
            actual_cov, _ = calculate_vntr_coverage(
                tools["samtools"],
                downsampled_bam,
                vntr_region if mode == "vntr" else "",
                threads,
                str(output_dir),
                f"{output_base}_post_downsample",
            ) if mode == "vntr" else calculate_target_coverage(
                tools["samtools"],
                downsampled_bam,
                rs_config.get("sample_target_bed", ""),
                threads,
                str(output_dir),
                f"{output_base}_post_downsample",
            )
            deviation_pct = ((actual_cov - target_coverage) / target_coverage) * 100
            logger.info(
                "Post-downsampling coverage: %.1fx (target: %.1fx, deviation: %+.1f%%)",
                actual_cov,
                target_coverage,
                deviation_pct,
            )
            if abs(deviation_pct) > 20:
                logger.warning(
                    "Post-downsampling coverage %.1fx deviates >20%% from target %.1fx",
                    actual_cov,
                    target_coverage,
                )
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `python -m pytest tests/read_simulator/test_alignment_stage.py -v --tb=short`
Expected: All tests pass (including existing ones — they may need the mock side_effect adjusted to return two values)

- [ ] **Step 5: Fix any existing tests that break**

Existing tests that mock `calculate_vntr_coverage` with a single return value will now fail because the function is called twice. Update those mocks to use `side_effect` with two return tuples. For tests where downsampling doesn't trigger (coverage below target), only one call happens — no change needed.

Specifically, check `test_vntr_downsampling_applied` and similar tests in `test_alignment_stage.py`. If they mock `calculate_vntr_coverage` returning a single tuple, change to:

```python
mocker.patch(f"{MODULE}.calculate_vntr_coverage", side_effect=[
    (200.0, "/tmp/depth.txt"),  # Pre-downsampling measurement
    (150.0, "/tmp/depth_post.txt"),  # Post-downsampling validation
])
```

- [ ] **Step 6: Commit**

```bash
git add muc_one_up/read_simulator/stages/alignment.py tests/read_simulator/test_alignment_stage.py
git commit -m "feat: add post-downsampling coverage validation logging (#85)"
```

---

### Task 3: Update penalty factor from 0.375 to 0.39

**Files:**
- Modify: `muc_one_up/read_simulator/vntr_efficiency.py:8,43,53`
- Modify: `muc_one_up/read_simulator/stages/alignment.py:95`
- Modify: `muc_one_up/config.py:318`
- Modify: `config.json:563`
- Modify: `config_experiment.json:566`
- Test: `tests/test_vntr_efficiency.py`
- Test: `tests/read_simulator/test_alignment_stage.py`

- [ ] **Step 1: Update tests to expect 0.39**

In `tests/test_vntr_efficiency.py`:

Line 63 — change:
```python
        assert model.penalty_factor == 0.375
```
to:
```python
        assert model.penalty_factor == 0.39
```

Line 208 — change:
```python
        model = VNTREfficiencyModel(penalty_factor=0.375, seed=42)
```
to:
```python
        model = VNTREfficiencyModel(penalty_factor=0.39, seed=42)
```

Line 248 — change:
```python
        model = VNTREfficiencyModel(penalty_factor=0.375, seed=42)
```
to:
```python
        model = VNTREfficiencyModel(penalty_factor=0.39, seed=42)
```

Line 254 — change:
```python
        assert stats["penalty_factor"] == 0.375
```
to:
```python
        assert stats["penalty_factor"] == 0.39
```

Line 383 — change:
```python
        model = VNTREfficiencyModel(penalty_factor=0.375, seed=42)
```
to:
```python
        model = VNTREfficiencyModel(penalty_factor=0.39, seed=42)
```

Line 400 — change:
```python
        assert stats["penalty_factor"] == 0.375
```
to:
```python
        assert stats["penalty_factor"] == 0.39
```

Line 454 — change:
```python
        for penalty in [0.357, 0.375, 0.395, 0.4, 0.5, 1.0]:
```
to:
```python
        for penalty in [0.357, 0.39, 0.395, 0.4, 0.5, 1.0]:
```

In `tests/read_simulator/test_alignment_stage.py`, line 122 — change:
```python
                "penalty_factor": 0.375,
```
to:
```python
                "penalty_factor": 0.39,
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/test_vntr_efficiency.py::TestVNTREfficiencyModelInit::test_default_initialization -v`
Expected: FAIL — `assert model.penalty_factor == 0.39` fails because default is still 0.375

- [ ] **Step 3: Update penalty factor in source code**

In `muc_one_up/read_simulator/vntr_efficiency.py`:

Line 8 — change docstring:
```python
using penalty factor 0.375 derived from comparison of 277 real Twist v2
```
to:
```python
using penalty factor 0.39 derived from 1,043 CerKiD Berlin Twist v2 exomes
```

Line 43 — change docstring:
```python
        penalty_factor: Fraction of VNTR reads to retain (default: 0.375)
```
to:
```python
        penalty_factor: Fraction of VNTR reads to retain (default: 0.39)
```

Line 53 — change default:
```python
        penalty_factor: float = 0.375,
```
to:
```python
        penalty_factor: float = 0.39,
```

In `muc_one_up/read_simulator/stages/alignment.py`, line 95 — change:
```python
            penalty_factor: float = vntr_config.get("penalty_factor", 0.375)
```
to:
```python
            penalty_factor: float = vntr_config.get("penalty_factor", 0.39)
```

In `muc_one_up/config.py`, line 318 — change:
```python
                            "default": 0.375,
```
to:
```python
                            "default": 0.39,
```

- [ ] **Step 4: Update config files**

In `config.json`, line 563 — change:
```python
      "penalty_factor": 0.375,
```
to:
```python
      "penalty_factor": 0.39,
```

In `config_experiment.json`, line 566 — change the same way.

- [ ] **Step 5: Run tests to verify they pass**

Run: `python -m pytest tests/test_vntr_efficiency.py tests/read_simulator/test_alignment_stage.py -v --tb=short`
Expected: All tests pass

- [ ] **Step 6: Lint and type-check**

Run: `ruff check muc_one_up/ tests/ && python -m mypy muc_one_up/`
Expected: Clean

- [ ] **Step 7: Commit**

```bash
git add muc_one_up/read_simulator/vntr_efficiency.py muc_one_up/read_simulator/stages/alignment.py muc_one_up/config.py config.json config_experiment.json tests/test_vntr_efficiency.py tests/read_simulator/test_alignment_stage.py
git commit -m "fix: update VNTR penalty factor 0.375 → 0.39 based on n=1,043 cohort (#86)"
```

---

### Task 4: Update documentation

**Files:**
- Modify: `docs/guides/vntr-capture-efficiency.md`

- [ ] **Step 1: Rewrite the Empirical Validation section**

Replace the current Dataset table (lines 36-39) with:

```markdown
### Dataset

| Category | Count | Source |
|----------|-------|--------|
| **Real samples** | 1,043 | CerKiD Berlin Twist Bioscience v2 exomes |
| **Validation method** | 200-sample subset | Random sample with per-base coverage (samtools depth -a) |
| **Region** | chr1:155,188,487-155,192,239 | MUC1 VNTR (hg38, 3,753 bp) |
```

Replace the Statistical Results table (lines 44-48) with:

```markdown
### Statistical Results

| Metric | Value | 95% CI |
|--------|-------|--------|
| **Penalty factor** | **0.39** | [0.20, 1.28] |
| **Cohort size** | 1,043 samples | All Twist Bioscience v2 |
| **Implied penalty (median)** | 0.392 | From observed ratio / simulated base ratio |
| **Effect size** | ~2.6x coverage reduction | VNTR reads retained |

### Cohort Coverage Profile

| Metric | Mean | Median | SD |
|--------|------|--------|-----|
| VNTR coverage (x) | 191.4 | 174.3 | 96.3 |
| Flanking coverage (x) | 47.9 | 46.4 | 10.3 |
| VNTR/Flanking ratio | 4.03 | 3.66 | 1.97 |
| VNTR % uncovered | 10.3% | 10.1% | 2.7% |
```

Replace the Coverage Ratios table (lines 52-57) with:

```markdown
### Coverage Ratios (VNTR:Flanking)

| Dataset | Mean Ratio | Median | Std Dev |
|---------|------------|--------|---------|
| **Simulated (no bias)** | 9.32 | — | ±1.25 |
| **Real (CerKiD cohort)** | 4.03 | 3.66 | ±1.97 |
| **Expected (penalty 0.39)** | ~3.6 | — | Theoretical |
```

- [ ] **Step 2: Update all 0.375 references to 0.39 throughout the doc**

Replace every occurrence of `0.375` with `0.39` in the file. Update log output examples, config examples, interpretation tables, and parameter defaults.

Key sections to update:
- Line 28: "penalty factor of 0.375" → "penalty factor of 0.39"
- Line 111: log example "Downsampling VNTR reads by 0.375" → "0.39"
- Line 124, 174, 244: config examples
- Line 152: parameter table default
- Line 200: algorithm step
- Line 260: interpretation table header
- Line 271: expected range "2.0 - 3.5" → "2.5 - 5.0"
- Line 284: troubleshooting "should be 0.375" → "should be 0.39"

Update interpretation guide:

```markdown
| Ratio | Interpretation |
|-------|---------------|
| **2.5 - 5.0** | Expected with penalty 0.39 |
| **8.0 - 10.0** | No bias applied (unrealistic) |
| **< 2.0** | Over-penalized |
| **> 6.0** | Under-penalized |
```

Update the expected behavior table:

```markdown
| Metric | No Bias | With Bias (0.39) |
|--------|---------|-------------------|
| VNTR coverage | 200x | 78x |
| Flanking coverage | 200x | 200x |
| Coverage ratio | 1.0 | ~2.6 |
| Reads retained | 100% | ~52% |
```

- [ ] **Step 3: Update version history**

Add entry:

```markdown
| **v0.44.0** | 2026-04-06 | Re-derived penalty 0.39 from n=1,043 CerKiD cohort; fixed coverage reporting bug |
```

- [ ] **Step 4: Update references section**

Replace the Data Sources subsection:

```markdown
### Data Sources
- **Calibration cohort:** 1,043 CerKiD Berlin Twist Bioscience v2 exomes (chr1 MUC1 region extracts)
- **Coverage analysis:** 200-sample random subset, samtools depth -a with MucOneUp VNTR coordinates
- **MUC1 VNTR region:** chr1:155,188,487-155,192,239 (hg38, 3,753 bp)
- **Previous calibration:** 277 real + 59 simulated samples (October 2024, penalty 0.375)
```

- [ ] **Step 5: Commit**

```bash
git add docs/guides/vntr-capture-efficiency.md
git commit -m "docs: update VNTR efficiency docs with n=1,043 cohort validation (#86)"
```

---

### Task 5: Run full test suite, lint, and type-check

**Files:** None (validation only)

- [ ] **Step 1: Run full test suite**

Run: `python -m pytest --tb=short -q`
Expected: All tests pass, coverage ≥ 30%

- [ ] **Step 2: Lint and format**

Run: `ruff check muc_one_up/ tests/ && ruff format --check muc_one_up/ tests/`
Expected: Clean

- [ ] **Step 3: Type-check**

Run: `python -m mypy muc_one_up/`
Expected: Clean

- [ ] **Step 4: Fix any issues found in steps 1-3**

If any failures, fix and re-run until clean.

- [ ] **Step 5: Final commit if any fixes were needed**

```bash
git add -A
git commit -m "chore: fix lint/type issues from coverage/penalty update"
```
