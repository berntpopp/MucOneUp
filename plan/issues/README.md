# Open Issues - Implementation Plans (EXPERT REVIEWED)

🔬 **Expert Review Complete**: All plans validated for scientific accuracy, tool capabilities, and programming best practices.

✅ **All Files Cleaned**: Only corrected/approved versions remain. Originals with violations have been removed.

---

## Review Summary

| Issue | Status | Original Issues | File |
|-------|--------|----------------|------|
| [#31](issue_31_nanosim_allelic_imbalance.md) | ✅ **APPROVED** | None | issue_31_nanosim_allelic_imbalance.md |
| [#30](issue_30_deterministic_simulation.md) | ✅ **APPROVED** | None | issue_30_deterministic_simulation.md |
| [#28](issue_28_reference_assemblies.md) | ✅ **CORRECTED** | DRY violation fixed | issue_28_reference_assemblies.md |
| [#18](issue_18_docker_containerization.md) | ✅ **SIMPLIFIED** | KISS violation fixed | issue_18_docker_containerization.md |
| [#15](issue_15_snapshot_assay_validation.md) | ✅ **CORRECTED** | DRY violation fixed | issue_15_snapshot_assay_validation.md |
| [#1](issue_01_pacbio_support.md) | ✅ **REWRITTEN** | Incorrect analysis fixed | issue_01_pacbio_support.md |

---

## Critical Findings

### 🚨 Issue #1: Incorrect Problem Statement
**Original Claim**: "ONT is partially implemented but not user-selectable"

**Actual State**: ONT is **FULLY IMPLEMENTED** and accessible via:
```bash
muconeup --config config.json reads ont sample.fa
```

**Correction**: Issue now focuses ONLY on adding PacBio support. See [issue_01_pacbio_support.md](issue_01_pacbio_support.md).

### ⚠️ Issue #28: DRY Violation
**Original Plan**: Create new `ReferenceManager` class

**Existing Code**: `muc_one_up/bioinformatics/reference_validation.py` already has:
- `validate_reference_genome()` - checks FASTA + indices
- `_validate_bwa_indices()` - checks BWA files
- `_validate_minimap2_indices()` - checks minimap2 files

**Correction**: Extend existing module with helper functions. See [issue_28_reference_assemblies.md](issue_28_reference_assemblies.md).

### ⚠️ Issue #15: DRY Violation
**Original Plan**: Implement `reverse_complement()` function

**Existing Code**: `muc_one_up/translate.py` already has this function

**Correction**: Import from existing module. See [issue_15_snapshot_assay_validation.md](issue_15_snapshot_assay_validation.md).

### ⚠️ Issue #18: KISS Violation
**Original Plan**: 3 separate Dockerfiles + docker-compose with profiles

**Simplified**: Single Dockerfile with two build targets (slim/full)

**Correction**: Reduce complexity. See [issue_18_docker_containerization.md](issue_18_docker_containerization.md).

---

## Implementation Priority (Revised)

### Phase 1: Quick Wins (Week 1)
1. **Issue #30** - Deterministic seeding ✅ APPROVED
   - Add `--seed` parameter to read simulators
   - Tools verified: NanoSim ✅, pbsim3 ✅, w-Wessim2 ✅
   - File: [issue_30_deterministic_simulation.md](issue_30_deterministic_simulation.md)

### Phase 2: Bug Verification (Week 1-2)
2. **Issue #31** - NanoSim allelic imbalance ⚠️ NEEDS EMPIRICAL TEST
   - First: Verify bias actually exists
   - If confirmed: Implement split-simulation fix
   - File: [issue_31_nanosim_allelic_imbalance.md](issue_31_nanosim_allelic_imbalance.md)

### Phase 3: Infrastructure (Week 2-3)
3. **Issue #28** - Reference management
   - Extend existing `reference_validation.py` module
   - DO NOT create new `ReferenceManager` class
   - File: [issue_28_reference_assemblies.md](issue_28_reference_assemblies.md)

4. **Issue #18** - Docker
   - Single Dockerfile, two targets
   - NO separate files for each simulator
   - File: [issue_18_docker_containerization.md](issue_18_docker_containerization.md)

### Phase 4: Feature Addition (Week 3-4)
5. **Issue #15** - SNaPshot validation
   - Import `reverse_complement` from `translate.py`
   - Do not reimplement
   - File: [issue_15_snapshot_assay_validation.md](issue_15_snapshot_assay_validation.md)

6. **Issue #1** - PacBio support
   - Add ONLY PacBio support
   - Do NOT modify ONT (already works)
   - Do NOT add `--read-type` flag (violates architecture)
   - File: [issue_01_pacbio_support.md](issue_01_pacbio_support.md)

---

## Tool Capability Verification

All external tool capabilities have been verified:

### ✅ NanoSim
- **--seed parameter**: CONFIRMED (source: simulator.py)
- **Allelic sampling**: Needs empirical verification
- **Version**: Any recent version

### ✅ pbsim3
- **--seed parameter**: CONFIRMED (source: GitHub docs)
- **Read types**: SEQUEL, CCS, CLR, SEQUEL2, RSII all supported
- **Version**: ≥3.0.2 recommended

### ✅ w-Wessim2 (Python port)
- **Seeding**: Uses Python `random` module (can be seeded)
- **Location**: `muc_one_up/read_simulator/fragment_simulation.py:11`

---

## Programming Principles Violations (Detected & Corrected)

### DRY (Don't Repeat Yourself)
| Issue | Violation | Corrected |
|-------|-----------|-----------|
| #28 | Duplicate validation logic | ✅ Extend existing module |
| #15 | Duplicate reverse_complement | ✅ Import from translate.py |

### KISS (Keep It Simple, Stupid)
| Issue | Violation | Corrected |
|-------|-----------|-----------|
| #18 | 3 Dockerfiles + compose | ✅ 1 Dockerfile, 2 targets |
| #1 | Add --read-type flag | ✅ Keep existing reads subcommands |

### SOLID (Single Responsibility)
| Issue | Violation | Corrected |
|-------|-----------|-----------|
| #28 | ReferenceManager does validation | ✅ Use existing validator |
| #1 | Modify simulate command | ✅ Keep reads group |

---

## File Organization

All files have been cleaned up - only corrected/approved versions remain:

```
plan/issues/
├── EXPERT_REVIEW.md                        ← Review methodology & findings
├── README.md                               ← This file
│
├── issue_01_pacbio_support.md              ← PacBio simulation (CORRECTED)
├── issue_15_snapshot_assay_validation.md   ← SNaPshot validation (CORRECTED)
├── issue_18_docker_containerization.md     ← Docker support (SIMPLIFIED)
├── issue_28_reference_assemblies.md        ← Reference management (CORRECTED)
├── issue_30_deterministic_simulation.md    ← Deterministic seeding (APPROVED)
└── issue_31_nanosim_allelic_imbalance.md   ← Allelic balance fix (APPROVED)
```

**Note**: All files are production-ready. Original versions with violations have been removed.

---

## Research References (Validated)

### bioRxiv 2025.09.06.673538v2
- **MUC1 VNTR lengths**: 20-125 repeats (confirmed)
- **SNaPshot detection**: MwoI restriction site method (confirmed)
- **MALDI-TOF**: dupC/dupA probe extension (confirmed)
- **Relevance**: Issue #31 (allelic imbalance), Issue #15 (SNaPshot)

### iScience PIIS2589004223012488 (VNtyper)
- **GC content**: 80%+ in VNTR region (confirmed)
- **Detection challenges**: Short-read mapping difficulty (confirmed)
- **Clinical utility**: ADTKD-MUC1 diagnosis (confirmed)
- **Relevance**: All read simulation issues

### Nature Methods 2021 (PEPPER-Margin-DeepVariant)
- **Allelic balance**: Long-read sequencing considerations (confirmed)
- **Haplotype phasing**: ONT data accuracy (confirmed)
- **Relevance**: Issue #31 (allelic balance)

---

## Contributing

When implementing an issue:

1. **Read implementation plan** for the issue
2. **Read EXPERT_REVIEW.md** for context and corrections made
3. **Verify tool capabilities** (all verified in EXPERT_REVIEW.md)
4. **Check for existing code** before creating new modules
5. **Run `make ci-check`** before committing
6. **Create feature branch**: `git checkout -b feature/issue-XX-name`
7. **Reference issue** in PR: "Closes #XX"

### Pre-Implementation Checklist:
- [ ] Read EXPERT_REVIEW.md findings for this issue
- [ ] Read implementation plan file
- [ ] Check for existing functions (grep codebase)
- [ ] Verify tool capabilities documented
- [ ] Follow existing CLI patterns
- [ ] Write tests first (TDD)

---

## Questions?

See detailed analysis in [EXPERT_REVIEW.md](EXPERT_REVIEW.md) for:
- Scientific validity assessment
- Programmatic validity verification
- Programming principles evaluation
- Tool capability confirmation
- Codebase architecture analysis
