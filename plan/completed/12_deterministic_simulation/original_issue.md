# Issue #30: Deterministic Simulation with Random Seeds

## Problem
Read simulators (NanoSim, w-Wessim2) use system time as default seed, making results non-reproducible. Scientific reproducibility requires deterministic simulations with explicit seeding.

## Research Background
- **NanoSim**: Supports `--seed` parameter (NanoSim-H fork documentation)
- **pbsim3**: Uses `--seed` (default: Unix time)
- **Reproducibility standards**: 2025 bioinformatics best practices require explicit seeding for all stochastic processes

## Implementation

### 1. Add Seed Parameters to Config
```json
"read_simulation": {
  "random_seed": 42,  // Master seed for reproducibility
  "use_deterministic_seed": true  // Default: false for backward compat
},
"nanosim_params": {
  "seed": null  // Auto-generated from master seed if null
}
```

### 2. Modify NanoSim Wrapper (`nanosim_wrapper.py:20`)
```python
def run_nanosim_simulation(
    ...
    seed: int | None = None,  # NEW parameter
):
    cmd_list = build_tool_command(...)

    if seed is not None:
        cmd_list.extend(["--seed", str(seed)])

    logging.info(f"NanoSim seed: {seed or 'system time (non-deterministic)'}")
```

### 3. Modify ONT Pipeline (`ont_pipeline.py:112`)
```python
# Extract master seed
master_seed = rs_config.get("random_seed")
use_deterministic = rs_config.get("use_deterministic_seed", False)

# Generate simulation seed
sim_seed = None
if use_deterministic and master_seed is not None:
    sim_seed = master_seed
    logging.info(f"Using deterministic seed: {sim_seed}")

fastq_file = run_nanosim_simulation(
    ...,
    seed=sim_seed
)
```

### 4. Illumina Pipeline Update (`pipeline.py`)
w-Wessim2 port already uses Python's `random` module:
```python
import random

def simulate_reads_pipeline(..., seed: int | None = None):
    if seed is not None:
        random.seed(seed)
        logging.info(f"Seeded random generator: {seed}")
```

### 5. Simulation Statistics Tracking
Add to `simulation_stats.json`:
```json
{
  "read_simulation": {
    "seed_used": 42,
    "deterministic": true,
    "simulator": "nanosim"
  }
}
```

## Testing
- **Reproducibility test**: Run same simulation 3x with seed=123, verify identical FASTQ
- **Null seed test**: Without seed, verify 3 runs produce different results
- **Backward compatibility**: Existing configs without seed should work unchanged

## CLI Enhancement
```bash
# New optional flag
muconeup --config config.json simulate --out-base test \
  --simulate-reads --read-seed 42

# Or use config.json
muconeup --config config.json simulate --out-base test --simulate-reads
```

## Documentation
Update README.md:
```markdown
### Reproducible Simulations

For reproducible read simulation, set a random seed:

**Config method:**
```json
"read_simulation": {"random_seed": 42, "use_deterministic_seed": true}
```

**CLI method:**
```bash
muconeup --config config.json simulate --out-base test \
  --simulate-reads --read-seed 42
```

**Note**: Seeds only affect read simulation (NanoSim/w-Wessim2), not VNTR structure
generation which uses system randomness.
```

## Files Modified
- `muc_one_up/read_simulator/nanosim_wrapper.py` (add seed param)
- `muc_one_up/read_simulator/ont_pipeline.py` (pass seed)
- `muc_one_up/read_simulator/pipeline.py` (Illumina seed)
- `muc_one_up/simulation_statistics.py` (track seed)
- `muc_one_up/cli/click_main.py` (add --read-seed flag)
- `tests/read_simulator/test_deterministic_simulation.py` (NEW)
