# Phase 2 Roadmap: Testing & Security

**Status:** 📋 PLANNED - Not Yet Implemented
**Prerequisite:** Phase 1 (Infrastructure) - ✅ Complete
**Estimated Effort:** 3-5 days

## What Was NOT Completed in Phase 1

This document outlines the remaining requirements from `05_code_quality_standards.md` that were deferred to Phase 2.

## Gap Analysis

### 1. ❌ Test Coverage: 30% → 60%+ (CRITICAL)

**Current State:**
- Total coverage: 30%
- Existing tests: 66 tests (all passing)
- No new tests added in Phase 1

**Severely Under-Tested Modules:**

| Module | Current Coverage | Target | Priority |
|--------|-----------------|--------|----------|
| `cli/main.py` | 9% | 80%+ | 🔴 CRITICAL |
| `cli/outputs.py` | 10% | 80%+ | 🔴 CRITICAL |
| `cli/analysis.py` | 11% | 80%+ | 🔴 CRITICAL |
| `read_simulator/fragment_simulation.py` | 7% | 60%+ | 🟡 HIGH |
| `read_simulator/pipeline.py` | 10% | 60%+ | 🟡 HIGH |
| `read_simulator/ont_pipeline.py` | 10% | 60%+ | 🟡 HIGH |
| `io.py` | 10% | 70%+ | 🟡 HIGH |
| `fasta_writer.py` | 17% | 70%+ | 🟡 HIGH |
| `toxic_protein_detector.py` | 0% | 60%+ | 🟢 MEDIUM |
| `simulation_statistics.py` | 17% | 60%+ | 🟢 MEDIUM |
| `snp_integrator.py` | 38% | 70%+ | 🟢 MEDIUM |

**Required Actions:**

#### Critical Priority (Week 1)
- [ ] Create `tests/test_cli_main.py`
  - Test main orchestration logic
  - Test argument parsing edge cases
  - Test error handling paths
  - Target: Bring cli/main.py to 80%+

- [ ] Create `tests/test_cli_outputs.py`
  - Test FASTA writing (single & dual mode)
  - Test structure file generation
  - Test mutation comments
  - Target: Bring cli/outputs.py to 80%+

- [ ] Create `tests/test_cli_analysis.py`
  - Test ORF prediction workflow
  - Test read simulation triggering
  - Test statistics generation
  - Target: Bring cli/analysis.py to 80%+

- [ ] Expand `tests/test_io.py`
  - Test structure file parsing
  - Test mutation info extraction
  - Test malformed file handling
  - Target: Bring io.py to 70%+

#### High Priority (Week 2)
- [ ] Create `tests/test_read_simulator_pipeline.py`
  - Test Illumina pipeline stages
  - Test file cleanup
  - Test error recovery
  - Mock external tools (reseq, bwa, samtools)
  - Target: Bring read_simulator/pipeline.py to 60%+

- [ ] Create `tests/test_read_simulator_ont.py`
  - Test ONT pipeline stages
  - Test NanoSim integration
  - Test alignment workflow
  - Mock external tools
  - Target: Bring read_simulator/ont_pipeline.py to 60%+

- [ ] Create `tests/test_fasta_writer.py`
  - Test FASTA writing with comments
  - Test multi-sequence handling
  - Test edge cases
  - Target: Bring fasta_writer.py to 70%+

#### Medium Priority (Week 3)
- [ ] Create `tests/test_snp_integrator.py`
  - Test SNP parsing from TSV
  - Test random SNP generation
  - Test SNP application
  - Test validation logic
  - Target: Bring snp_integrator.py to 70%+

- [ ] Create `tests/test_toxic_protein.py`
  - Test protein sequence analysis
  - Test toxic feature detection
  - Test ORF scanning
  - Target: Bring toxic_protein_detector.py to 60%+

- [ ] Create `tests/test_simulation_stats.py`
  - Test statistics generation
  - Test JSON output
  - Test mutation tracking
  - Target: Bring simulation_statistics.py to 60%+

**Implementation Strategy:**

```python
# Example test structure for cli/main.py
# tests/test_cli_main.py

import pytest
from pathlib import Path
from unittest.mock import Mock, patch
from muc_one_up.cli.main import main, run_simulation
from muc_one_up.exceptions import ConfigurationError

class TestMainOrchestration:
    """Test main CLI orchestration logic."""

    def test_run_simulation_single_iteration(self, tmp_path, sample_config):
        """Test single simulation iteration completes successfully."""
        # Setup test fixtures
        # Call run_simulation
        # Assert files created
        # Assert logs correct
        pass

    def test_run_simulation_dual_mutation_mode(self, tmp_path, sample_config):
        """Test dual mutation mode creates both normal and mutated outputs."""
        pass

    def test_run_simulation_handles_missing_config(self, tmp_path):
        """Test graceful error handling for missing config."""
        with pytest.raises(ConfigurationError):
            # Test error path
            pass

    def test_run_simulation_with_snps(self, tmp_path, sample_config):
        """Test SNP integration workflow."""
        pass

    # Add 20+ more tests to reach 80% coverage
```

### 2. ❌ Security Scanning (Bandit) - NOT IMPLEMENTED

**Current State:**
- No security scanning configured
- Unknown security vulnerabilities
- No automated security checks

**Required Actions:**

#### Add Bandit to Pre-commit
```yaml
# Update .pre-commit-config.yaml
- repo: https://github.com/PyCQA/bandit
  rev: 1.7.9
  hooks:
    - id: bandit
      args: ["-c", "pyproject.toml"]
      additional_dependencies: ["bandit[toml]"]
```

#### Add Bandit Configuration
```toml
# Add to pyproject.toml
[tool.bandit]
exclude_dirs = ["tests", "build", "dist"]
skips = ["B101"]  # Skip assert_used (used in tests)
```

#### Add Security Workflow
```yaml
# Create .github/workflows/security.yml
name: Security Scan

on:
  push:
    branches: [main, dev, dev/*]
  schedule:
    - cron: '0 0 * * 0'  # Weekly
  workflow_dispatch:

jobs:
  bandit:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: "3.10"
      - name: Install bandit
        run: pip install bandit[toml]
      - name: Run bandit
        run: bandit -r muc_one_up/ -f json -o bandit-report.json
      - name: Upload results
        uses: actions/upload-artifact@v3
        with:
          name: bandit-report
          path: bandit-report.json
```

#### Address Findings
- [ ] Run bandit locally: `bandit -r muc_one_up/`
- [ ] Fix any HIGH or MEDIUM severity issues
- [ ] Document any accepted risks
- [ ] Ensure 0 critical security issues

### 3. ⚠️ Progressive Coverage Thresholds

**Current State:**
- Coverage threshold: 30% (matches current coverage)
- No progression plan

**Required Actions:**

#### Update CI/CD for Progressive Thresholds
```yaml
# .github/workflows/test.yml
- name: Check coverage threshold
  run: |
    pytest --cov=muc_one_up --cov-fail-under=30  # Current
    # After Phase 2 Week 1: --cov-fail-under=40
    # After Phase 2 Week 2: --cov-fail-under=50
    # After Phase 2 Week 3: --cov-fail-under=60
```

#### Create Coverage Tracking Issue
- [ ] Create GitHub issue to track coverage milestones
- [ ] Set up coverage badge in README
- [ ] Monitor coverage trends over time

## Implementation Timeline

### Week 1: Critical CLI Testing
**Goal:** 30% → 40% coverage
- Day 1-2: `test_cli_main.py` (bring cli/main.py to 80%+)
- Day 3: `test_cli_outputs.py` (bring cli/outputs.py to 80%+)
- Day 4: `test_cli_analysis.py` (bring cli/analysis.py to 80%+)
- Day 5: `test_io.py` expansion (bring io.py to 70%+)

**Deliverable:** Critical CLI paths fully tested, coverage at 40%

### Week 2: Read Simulator Testing
**Goal:** 40% → 50% coverage
- Day 1-2: `test_read_simulator_pipeline.py` (mock external tools)
- Day 3: `test_read_simulator_ont.py` (mock NanoSim)
- Day 4: `test_fasta_writer.py` (bring to 70%+)
- Day 5: Security scanning setup (Bandit)

**Deliverable:** Read simulation tested, security scanning active, coverage at 50%

### Week 3: Complete Remaining Modules
**Goal:** 50% → 60%+ coverage
- Day 1: `test_snp_integrator.py` (bring to 70%+)
- Day 2: `test_toxic_protein.py` (bring to 60%+)
- Day 3: `test_simulation_stats.py` (bring to 60%+)
- Day 4: Fix any security findings from Bandit
- Day 5: Final validation, documentation update

**Deliverable:** 60%+ coverage achieved, 0 critical security issues, Phase 2 complete

## Success Criteria

### Test Coverage
- ✅ Overall coverage ≥60%
- ✅ All CLI modules ≥80%
- ✅ Core modules ≥70%
- ✅ Read simulator ≥60%
- ✅ All tests passing
- ✅ No flaky tests

### Security
- ✅ Bandit integrated in pre-commit
- ✅ Bandit in CI/CD
- ✅ Security workflow running weekly
- ✅ 0 HIGH or CRITICAL security issues
- ✅ All MEDIUM issues documented/accepted

### Automation
- ✅ Progressive coverage thresholds enforced
- ✅ Coverage badge in README
- ✅ Security scan reports uploaded
- ✅ Coverage tracking issue created

## Verification Commands

```bash
# After Phase 2 completion
pytest --cov=muc_one_up --cov-report=term-missing --cov-fail-under=60
bandit -r muc_one_up/
pre-commit run --all-files  # Should include Bandit

# Check all quality gates pass
make check  # Should enforce 60% coverage
```

## Why Phase 2 Was Deferred

**Decision:** Focus Phase 1 on infrastructure (linting, typing, CI/CD) to establish quality gates before investing in test development.

**Rationale:**
1. Infrastructure provides immediate value (catches regressions)
2. Test development is time-intensive (3-5 days for 60% coverage)
3. Security issues are low-priority given defensive use case
4. Existing 30% coverage provides baseline protection

**Trade-offs:**
- ✅ Fast delivery of quality infrastructure
- ✅ No regressions introduced
- ⚠️ Coverage gap remains (30% vs 60% target)
- ⚠️ Unknown security vulnerabilities

## When to Implement Phase 2

**Recommended Timing:**
- Before production deployment
- Before major feature additions
- When bandwidth allows 3-5 day investment
- If security concerns arise

**Can Be Deferred If:**
- Tool is used internally only
- Risk tolerance is high
- Time constraints are critical
- Coverage is "good enough" for now

## Documentation Updates Needed After Phase 2

- [ ] Update REFACTOR_SUMMARY.md with Phase 2 completion
- [ ] Add coverage badge to README
- [ ] Document security findings and resolutions
- [ ] Update CONTRIBUTING.md with testing guidelines
- [ ] Archive this roadmap document

---

**Note:** This roadmap represents the remaining ~30% of `05_code_quality_standards.md` requirements. Phase 1 (infrastructure) successfully completed the foundational 70%.
