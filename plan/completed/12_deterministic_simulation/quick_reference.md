# Quick Reference: Issue #30 Implementation

**⏱️ 4-6 hours total** | **📝 6 files** | **🧪 10+ tests** | **🎯 Risk: LOW**

---

## 🚀 Implementation Checklist

### Phase 1: Config Schema (30 min)

**File**: `muc_one_up/config.py`

```python
# Line 232 (in read_simulation properties)
"seed": {"type": ["number", "null"]},

# Line 117 (in nanosim_params properties)
"seed": {"type": ["number", "null"]},
```

**Test**: Run `pytest tests/test_config.py -v`

---

### Phase 2: CLI Commands (40 min)

**File**: `muc_one_up/cli/click_main.py`

#### Illumina Command (after line 345)
```python
@click.option(
    "--seed",
    type=int,
    default=None,
    help="Random seed for reproducibility.",
)
def illumina(ctx, input_fastas, out_dir, out_base, coverage, threads, seed):
    # ... existing code ...

    # ADD after line 385
    if seed is not None:
        config["read_simulation"]["seed"] = seed
        logging.info(f"Using random seed: {seed}")
```

#### ONT Command (after line 457)
```python
@click.option(
    "--seed",
    type=int,
    default=None,
    help="Random seed for reproducibility.",
)
def ont(ctx, input_fastas, out_dir, out_base, coverage, min_read_length, seed):
    # ... existing code ...

    # ADD after line 496
    if seed is not None:
        config["nanosim_params"]["seed"] = seed
        logging.info(f"Using random seed: {seed}")
```

**Test**: Run `pytest tests/test_click_cli.py -v`

---

### Phase 3: Fragment Simulation (30 min)

**File**: `muc_one_up/read_simulator/fragment_simulation.py`

```python
# Line 306 - Update signature
def simulate_fragments(
    ref_fa: str,
    syser_file: str,
    psl_file: str,
    read_number: int,
    fragment_size: int,
    fragment_sd: int,
    min_fragment: int,
    bind: float,
    output_fragments: str,
    seed: int | None = None,  # ← ADD
) -> None:
    """..."""

    # Line 325 - ADD seed initialization (before logging)
    if seed is not None:
        random.seed(seed)
        logging.info(f"Fragment simulation using random seed: {seed}")

    # ... rest unchanged ...
```

**Test**: Run `pytest tests/read_simulator/test_fragment_simulation.py -v`

---

### Phase 4: Illumina Pipeline (15 min)

**File**: `muc_one_up/read_simulator/pipeline.py`

```python
# Line 202 - Extract seed
seed = rs_config.get("seed")

# Line 213 - Pass to simulate_fragments
simulate_fragments(
    no_ns_fa,
    syser_fq,
    psl_file,
    read_number,
    fragment_size,
    fragment_sd,
    min_fragment,
    bind,
    fragments_fa,
    seed=seed,  # ← ADD
)
```

**Test**: Manual run with `--seed 42`

---

### Phase 5: NanoSim Wrapper (30 min)

**File**: `muc_one_up/read_simulator/wrappers/nanosim_wrapper.py`

```python
# Line 30 - Update signature
def run_nanosim_simulation(
    nanosim_cmd: str,
    reference_fasta: str,
    output_prefix: str,
    training_model: str,
    coverage: float,
    threads: int = 4,
    min_read_length: int | None = None,
    max_read_length: int | None = None,
    other_options: str = "",
    timeout: int = 3600,
    seed: int | None = None,  # ← ADD
) -> str:
    """..."""

    # Line 74 - Add to command (after coverage)
    if seed is not None:
        cmd_list.extend(["--seed", str(seed)])
        logging.info(f"[NanoSim] Using random seed: {seed}")

    # ... rest unchanged ...
```

**Test**: Run `pytest tests/read_simulator/test_nanosim_wrapper.py -v`

---

### Phase 6: ONT Pipeline (10 min)

**File**: `muc_one_up/read_simulator/ont_pipeline.py`

```python
# Line 98 - Extract seed
seed = ns_params.get("seed")

# Line 122 - Pass to run_nanosim_simulation
fastq_file = run_nanosim_simulation(
    nanosim_cmd=nanosim_cmd,
    reference_fasta=input_fa,
    output_prefix=output_prefix,
    training_model=training_model,
    coverage=float(coverage),
    threads=int(threads),
    min_read_length=min_read_length,
    max_read_length=max_read_length,
    other_options=other_options,
    seed=seed,  # ← ADD
)
```

**Test**: Manual run with `--seed 42`

---

## 🧪 Testing Commands

```bash
# Run all tests
pytest tests/ -v

# Run specific test files
pytest tests/test_config.py -v
pytest tests/test_click_cli.py -v
pytest tests/read_simulator/test_fragment_simulation.py -v
pytest tests/read_simulator/test_nanosim_wrapper.py -v

# Manual determinism test (Illumina)
muconeup --config config.json reads illumina test.fa --seed 42 --out-base run1
muconeup --config config.json reads illumina test.fa --seed 42 --out-base run2
diff run1_R1.fastq.gz run2_R1.fastq.gz  # Should be identical

# Manual determinism test (ONT)
muconeup --config config.json reads ont test.fa --seed 42 --out-base run1
muconeup --config config.json reads ont test.fa --seed 42 --out-base run2
diff run1_ont_aligned_reads.fastq run2_ont_aligned_reads.fastq  # Should be identical
```

---

## 📋 Test Files to Create

### 1. Fragment Simulation Tests
**File**: `tests/read_simulator/test_fragment_simulation.py` (append)

```python
class TestFragmentSimulationSeeding:
    def test_same_seed_produces_identical_fragments(self, tmp_path):
        # ... see full plan for implementation ...

    def test_different_seeds_produce_different_fragments(self, tmp_path):
        # ... see full plan for implementation ...

    def test_no_seed_produces_random_output(self, tmp_path):
        # ... see full plan for implementation ...
```

### 2. NanoSim Wrapper Tests
**File**: `tests/read_simulator/test_nanosim_wrapper.py` (append)

```python
def test_run_nanosim_with_seed(mocker, tmp_path):
    # ... see full plan for implementation ...

def test_run_nanosim_without_seed(mocker, tmp_path):
    # ... see full plan for implementation ...
```

### 3. Config Tests
**File**: `tests/test_config.py` (append)

```python
def test_read_simulation_accepts_seed(self, tmp_path):
    # ... see full plan for implementation ...

def test_nanosim_params_accepts_seed(self, tmp_path):
    # ... see full plan for implementation ...

def test_seed_can_be_null(self, tmp_path):
    # ... see full plan for implementation ...
```

### 4. CLI Tests
**File**: `tests/test_click_cli.py` (create or append)

```python
def test_illumina_reads_accepts_seed_parameter(tmp_path, mocker):
    # ... see full plan for implementation ...

def test_ont_reads_accepts_seed_parameter(tmp_path, mocker):
    # ... see full plan for implementation ...
```

---

## 📚 Documentation Updates

### Update CLAUDE.md

**Add after "Read Simulation" section**:

```markdown
### Deterministic Read Simulation

Generate reproducible reads using the `--seed` parameter:

```bash
# Illumina reads with seed
muconeup --config config.json reads illumina sample.fa --seed 42

# ONT reads with seed
muconeup --config config.json reads ont sample.fa --seed 42
```

**Important**: Same seed guarantees identical output ONLY when:
- Using identical input files
- Running on same platform/architecture
- Using same tool versions (NanoSim, Python)
```

### Update config.json Example

**Add seed to examples**:

```json
{
  "read_simulation": {
    "seed": 42
  },
  "nanosim_params": {
    "seed": 42
  }
}
```

---

## ✅ Pre-Commit Checklist

- [ ] All 6 files modified
- [ ] No hardcoded seeds
- [ ] All new tests pass
- [ ] All existing tests pass
- [ ] Manual verification complete
- [ ] Documentation updated
- [ ] Type hints correct
- [ ] Docstrings updated
- [ ] Logging added
- [ ] Backward compatible

---

## 🐛 Common Issues & Fixes

### Issue: Tests fail with "seed not found"
**Fix**: Make sure seed parameter has default `None` in all function signatures

### Issue: NanoSim doesn't accept --seed
**Fix**: Check NanoSim version (≥ 3.0.0 required)

### Issue: Global random state pollution
**Fix**: Ensure `random.seed()` only called at top level, not in helper functions

### Issue: Config validation fails
**Fix**: Use `["number", "null"]` in schema, not just `"number"`

---

## 📞 Help & References

- **Full Implementation Plan**: `plan/issue_30_deterministic_simulation_implementation_plan.md`
- **Summary**: `plan/IMPLEMENTATION_SUMMARY_ISSUE_30.md`
- **GitHub Issue**: https://github.com/berntpopp/MucOneUp/issues/30
- **NanoSim Docs**: Check `nanosim-simulator --help` for seed syntax

---

**Status**: ⏳ Ready to Implement
**Priority**: Phase 1 (Week 1)
**Last Updated**: 2025-10-19
