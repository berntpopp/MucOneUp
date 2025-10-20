# SNaPshot Assay Validation

In-silico validation of SNaPshot assay for detecting MUC1 VNTR mutations. Simulates complete laboratory workflow: PCR amplification, restriction enzyme digest, single-base extension, and fluorescence detection.

---

## Scientific Background

### SNaPshot Technology

**SNaPshot** is a multiplexed single-nucleotide polymorphism (SNP) detection method using:

- **Single-base extension** - Primer anneals adjacent to mutation site
- **Fluorescent ddNTPs** - Each base (A, C, G, T) labeled with different fluorophore
- **Capillary electrophoresis** - Detects fluorescent peak at specific size

### Clinical Application

**ADTKD-MUC1 Detection:**
- **Challenge:** MUC1 VNTR region (GC-rich repetitive structure) difficult to sequence
- **Solution:** SNaPshot assay selectively detects frameshift mutations
- **Mechanism:** Restriction digest eliminates wild-type, enriches mutant products

---

## Workflow Overview

```
┌─────────────────────────────────────────────────────────────┐
│               SNaPshot Validation Workflow                  │
└─────────────────────────────────────────────────────────────┘

Input: Genomic DNA template (simulated FASTA)
  │
  ▼
┌──────────────────────────┐
│ Step 1: PCR Amplification│
│ • Forward primer binds   │
│ • Reverse primer binds   │
│ • Generates amplicon(s)  │
│ • Size filter applied    │
└────────┬─────────────────┘
         │
         ▼
┌──────────────────────────┐
│ Step 2: Restriction      │
│         Digest           │
│ • MwoI enzyme recognizes │
│   GCNNNNNNNGC site       │
│ • Cuts wild-type products│
│ • Mutant products survive│
└────────┬─────────────────┘
         │
         ▼
┌──────────────────────────┐
│ Step 3: SNaPshot         │
│         Extension        │
│ • Primer anneals to      │
│   surviving amplicons    │
│ • ddNTP incorporates     │
│ • Fluorophore detected   │
└────────┬─────────────────┘
         │
         ▼
┌──────────────────────────┐
│ Output: Mutation Status  │
│ • Fluorescence detected  │
│ • Mutation type inferred │
│ • Peak size calculated   │
└──────────────────────────┘
```

---

## Step 1: PCR Amplification

### Purpose
Amplify VNTR region containing mutation site using specific primers.

### Simulator: `PCRSimulator`

**Input:**
- Template DNA sequence
- Forward primer: 5' → 3'
- Reverse primer: 5' → 3' (auto reverse-complemented)
- Size range filter (min-max bp)

**Process:**

1. **Find binding sites** for both primers (handles multiple sites)
2. **Generate amplicons** for all valid forward-reverse pairs
3. **Filter by size** (e.g., 300-800 bp)
4. **Categorize** by mutation pattern (mutant/normal/unknown)
5. **Prioritize** mutant products first

**Output:**
```json
{
  "success": true,
  "product_count": 3,
  "non_specific": true,
  "products": [
    {
      "sequence": "ACGT...TGCA",
      "length": 450,
      "forward_pos": 1523,
      "reverse_pos": 1973,
      "flanked_region": "CGT...TGC",
      "mutation_type": "mutant"
    }
  ]
}
```

**Implementation:** `muc_one_up/analysis/snapshot_validator.py:27-181`

---

## Step 2: Restriction Digest Selection

### Purpose
Selectively digest products based on presence/absence of restriction site.

### Enzyme: MwoI

**Recognition site:** `GCNNN|NNNGC` (7 bp between GC pairs, cuts in middle)

**Mechanism:**
```
Wild-type:  GCNNNNNNNGC → Digested (cut into fragments)
Mutant:     GC───────GC → Survives (no site due to mutation)
            (site disrupted by frameshift)
```

### Simulator: `DigestSimulator`

**Input:**
- List of PCR amplicons
- Enzyme specification (e.g., "MwoI")

**Process:**

1. **Search for recognition sites** in each amplicon (using Bio.Restriction)
2. **Classify amplicons:**
   - **Has site(s)** → Gets digested (destroyed)
   - **No sites** → Survives digest

**Output:**
```json
{
  "total_products": 3,
  "digested_count": 2,
  "survivor_count": 1,
  "survivors": [
    {
      "sequence": "ACGT...TGCA",
      "has_restriction_site": false,
      "site_count": 0,
      "survives_digest": true,
      "mutation_type": "mutant"
    }
  ],
  "digested": [...]
}
```

**Key Insight:** Frameshift mutations disrupt MwoI recognition site, allowing selective enrichment.

**Implementation:** `muc_one_up/analysis/snapshot_validator.py:183-276`

---

## Step 3: SNaPshot Extension

### Purpose
Single-base extension to identify specific mutation nucleotide.

### Primer Design

**Extension primer** anneals immediately 5' of mutation site:

```
Template:  5'-...ACGT[C]GATC...-3'  (mutant, C insertion)
                     ↑
Primer:    3'-...TGCA-5'
                     ↓
Extension:         ddCTP incorporated
                   (Black fluorophore)
```

### Simulator: `SnapshotExtensionSimulator`

**Input:**
- Amplicon sequence (from digest survivors)
- Primer name (e.g., "primer_7c")

**Process:**

1. **Find primer binding site** in amplicon
2. **Identify next base** after primer 3' end
3. **Determine ddNTP** that incorporates (complementary to next base)
4. **Map to fluorophore** based on base identity

**Fluorophore Mapping:**

| Base | ddNTP | Fluorophore | Dye Name | Interpretation |
|------|-------|-------------|----------|----------------|
| **C** | ddCTP | **Black** | ddCTP | Mutant pattern (dupC insertion) |
| A | ddATP | Green | ddATP | Incomplete extension or alternate |
| G | ddGTP | Blue | ddGTP | Unexpected (verify sequence) |
| T | ddTTP | Red | ddTTP | Unexpected (verify sequence) |

**Output:**
```json
{
  "binds": true,
  "position": 234,
  "next_base": "C",
  "ddNTP": "ddCTP",
  "fluorophore": "Black",
  "fluorophore_dye": "ddCTP",
  "peak_size": 23,
  "interpretation": "Mutant pattern detected (C signal)"
}
```

**Peak size** = primer length + 1 (for capillary electrophoresis)

**Implementation:** `muc_one_up/analysis/snapshot_validator.py:279-374`

---

## Complete Workflow Orchestration

### Validator: `SnapshotValidator`

Coordinates all three steps for complete assay simulation.

**Input:**
- Template genomic sequence (FASTA)
- Mutation name (e.g., "dupC")
- Configuration (primers, enzyme, patterns)

**Configuration Structure:**

```json
{
  "snapshot_validation": {
    "dupC": {
      "pcr": {
        "forward_primer": "ACGTACGTACGT",
        "reverse_primer": "TGCATGCATGCA",
        "reverse_needs_rc": true,
        "max_products": 10,
        "size_range": {"min": 300, "max": 800}
      },
      "digest": {
        "enzyme": "MwoI",
        "recognition_site": "GCNNNNNNNGC"
      },
      "snapshot": {
        "primers": {
          "primer_7c": "GCCACCCCTTCTCCC"
        },
        "fluorophore_map": {
          "C": {"color": "Black", "dye": "ddCTP"},
          "A": {"color": "Green", "dye": "ddATP"},
          "G": {"color": "Blue", "dye": "ddGTP"},
          "T": {"color": "Red", "dye": "ddTTP"}
        }
      },
      "validation": {
        "mutant_pattern": "GCCCACGATGTCACCTCAGCC",
        "normal_pattern": "GCCCACGGTGTCACCTCGGCC",
        "expected_mutant_fluorescence": "Black"
      }
    }
  }
}
```

**Output:**

```json
{
  "pcr_results": {
    "success": true,
    "product_count": 3,
    "non_specific": true,
    "products": [...]
  },
  "digest_results": {
    "survivor_count": 1,
    "digested_count": 2,
    "survivors": [...]
  },
  "snapshot_results": [
    {
      "amplicon": {...},
      "extension": {
        "next_base": "C",
        "fluorophore": "Black",
        "peak_size": 23
      }
    }
  ],
  "mutation_detected": true,
  "expected_fluorescence": "Black",
  "summary": "dupC mutation DETECTED: 1 survivor(s), fluorescence: Black (ddCTP)"
}
```

**Implementation:** `muc_one_up/analysis/snapshot_validator.py:377-591`

---

## Detailed Flow Diagram

```
┌────────────────────────────────────────────────────────────────────┐
│                      SNaPshot Validation                           │
│                    (Component Interaction)                         │
└────────────────────────────────────────────────────────────────────┘

                    ┌──────────────────────┐
                    │  Input: Template DNA │
                    │  (simulated FASTA)   │
                    └───────────┬──────────┘
                                │
                                ▼
                    ┌──────────────────────┐
                    │   PCRSimulator       │
                    │ ─────────────────    │
                    │ • _find_binding_sites│
                    │   (forward primer)   │
                    │ • _find_binding_sites│
                    │   (reverse primer)   │
                    │ • Generate amplicons │
                    │   for all F×R pairs  │
                    │ • Filter by size     │
                    │   (min-max bp)       │
                    │ • Categorize by      │
                    │   mutation pattern   │
                    │ • Prioritize mutant  │
                    │   products first     │
                    └───────────┬──────────┘
                                │
                  ┌─────────────┴─────────────┐
                  │  PCR Products (List)      │
                  │  [amplicon1, amplicon2,   │
                  │   amplicon3, ...]         │
                  └─────────────┬─────────────┘
                                │
                                ▼
                    ┌──────────────────────┐
                    │  DigestSimulator     │
                    │ ─────────────────    │
                    │ • Load enzyme (MwoI) │
                    │   from Bio.Restriction│
                    │ • For each amplicon: │
                    │   - Search for sites │
                    │   - Has site? →      │
                    │     Digested list    │
                    │   - No site? →       │
                    │     Survivor list    │
                    └───────────┬──────────┘
                                │
                  ┌─────────────┴─────────────┐
                  │  Digest Results           │
                  │  • Survivors: [amp1]      │
                  │  • Digested: [amp2, amp3] │
                  └─────────────┬─────────────┘
                                │
                                ▼
                    ┌──────────────────────────┐
                    │ SnapshotExtensionSimulator│
                    │ ─────────────────────     │
                    │ For each survivor:        │
                    │ • Find primer binding     │
                    │ • Get next base after     │
                    │   primer 3' end           │
                    │ • Determine ddNTP         │
                    │ • Map to fluorophore      │
                    │   (Black/Green/Blue/Red)  │
                    │ • Calculate peak size     │
                    └───────────┬───────────────┘
                                │
                  ┌─────────────┴─────────────┐
                  │  Extension Results        │
                  │  [{next_base: "C",        │
                  │    fluorophore: "Black",  │
                  │    peak_size: 23}]        │
                  └─────────────┬─────────────┘
                                │
                                ▼
                    ┌──────────────────────────┐
                    │  Mutation Detection Logic│
                    │ ──────────────────────   │
                    │ • Any survivor shows     │
                    │   "C" (Black)?           │
                    │   → mutation_detected=True│
                    │ • Else:                  │
                    │   → mutation_detected=False│
                    └───────────┬──────────────┘
                                │
                                ▼
                    ┌──────────────────────────┐
                    │  Output: Validation      │
                    │  Report (JSON)           │
                    │  • mutation_detected     │
                    │  • fluorescence colors   │
                    │  • summary string        │
                    └──────────────────────────┘
```

---

## Usage Examples

### Command Line

```bash
# Validate dupC mutation in simulated sample
muconeup --config config.json analyze snapshot-validate \
  sample.001.mut.fa \
  --mutation dupC \
  --output validation_results.json

# Validate on batch of samples
for file in *.mut.fa; do
  muconeup --config config.json analyze snapshot-validate \
    "$file" --mutation dupC --output "${file%.fa}.validation.json"
done
```

### Python API

```python
from Bio import SeqIO
from muc_one_up.analysis.snapshot_validator import SnapshotValidator
from muc_one_up.config import load_config

# Load configuration
config = load_config("config.json")

# Initialize validator for dupC mutation
validator = SnapshotValidator(config, mutation_name="dupC")

# Load template sequence
template_seq = str(SeqIO.read("sample.001.mut.fa", "fasta").seq)

# Run complete validation workflow
results = validator.validate_complete_workflow(template_seq)

# Check results
if results["mutation_detected"]:
    print(f"✓ Mutation detected: {results['expected_fluorescence']}")
    print(f"  Survivors: {results['digest_results']['survivor_count']}")
else:
    print(f"✗ No mutation detected")

print(f"\nSummary: {results['summary']}")
```

### Output Interpretation

**Mutation Detected:**
```
dupC mutation DETECTED: 1 survivor(s), fluorescence: Black (ddCTP)
```

**No Mutation:**
```
No dupC mutation detected: all products digested (normal phenotype)
```

---

## Biological Mechanism

### Why Restriction Digest Works

**Wild-type MUC1 VNTR:**
```
5'-GC[CCCACGGT]GC-3'  ← MwoI recognition site intact
      │││││││││
      NNNNNNNN  (7 bases between GC pairs)
```
→ **Digested** by MwoI

**Mutant MUC1 VNTR (dupC frameshift):**
```
5'-GCC[CACGATG]TC-3'  ← MwoI site disrupted (extra C shifts frame)
       │││││││││
       NNNNNNNN  (site no longer recognized)
```
→ **Survives** digest

### Mutation Types Detected

| Mutation | Mechanism | Next Base | Fluorophore |
|----------|-----------|-----------|-------------|
| **dupC** | C insertion | **C** | **Black** |
| dupT | T insertion | T | Red |
| delC | C deletion | A | Green |
| Complex | Multiple changes | Variable | Multiple |

---

## Validation and Quality Control

### Test Cases

Implemented in `tests/test_snapshot_validator.py`:

**Test 1: Positive control (dupC mutant)**
```python
# Template with dupC mutation → expect detection
template = "...GCCCACGATGTCACCTC..."  # mutant pattern
results = validator.validate_complete_workflow(template)
assert results["mutation_detected"] == True
assert results["expected_fluorescence"] == "Black"
```

**Test 2: Negative control (wild-type)**
```python
# Template without mutation → expect no detection
template = "...GCCCACGGTGTCACCTC..."  # normal pattern
results = validator.validate_complete_workflow(template)
assert results["mutation_detected"] == False
```

**Test 3: Non-specific PCR**
```python
# Multiple primer binding sites → prioritize mutant
# Should still detect mutation if any survivor is mutant
```

### Performance Metrics

- **Sensitivity:** Detects mutation if ≥1 mutant survivor
- **Specificity:** No false positives (wild-type all digested)
- **Runtime:** < 100 ms per sample

---

## Configuration Guidelines

### Primer Design

**Extension primer requirements:**
- Anneals immediately 5' of mutation site
- Length: 15-25 bp (typical)
- Tm: 50-60°C
- No self-complementarity

**PCR primer requirements:**
- Flank VNTR region
- Amplicon size: 300-800 bp (typical for SNaPshot)
- Must include mutation site and restriction site

### Enzyme Selection

**MwoI characteristics:**
- Recognition: `GCNNNNNNNGC`
- Cuts in middle: `GCNNN|NNNGC`
- Thermostable (good for PCR cleanup)
- Common in MUC1 VNTR wild-type

**Alternative enzymes:**
- Must have recognition site disrupted by mutation
- Use `Bio.Restriction` to validate

---

## Troubleshooting

### PCR Failure

**Symptom:** `"success": false`, `"product_count": 0`

**Causes:**
- Primers not in template sequence
- Reverse primer not reverse-complemented (check `reverse_needs_rc`)
- Size filter too restrictive

**Solutions:**
```bash
# Check primer presence
grep "ACGTACGTACGT" template.fa

# Adjust size range in config
"size_range": {"min": 200, "max": 1000}
```

### All Products Digested

**Symptom:** `"survivor_count": 0`, `"digested_count": 3`

**Causes:**
- Mutation not present (wild-type template)
- Mutation doesn't disrupt restriction site
- Wrong enzyme specified

**Solutions:**
```bash
# Verify mutation pattern in template
grep "GCCCACGATGTCACCTC" template.fa

# Check restriction sites
python -c "from Bio.Restriction import MwoI; from Bio.Seq import Seq; print(MwoI.search(Seq('YOUR_SEQUENCE')))"
```

### No Fluorescence Signal

**Symptom:** `"binds": false`

**Causes:**
- Extension primer not in amplicon sequence
- Primer at very end of amplicon (no base to extend)

**Solutions:**
- Verify primer design
- Check amplicon sequence in digest survivors

---

## Implementation Details

### Module Structure

**Source:** `muc_one_up/analysis/snapshot_validator.py`

**Classes:**
- `PCRSimulator` (lines 27-181) - Amplification with multi-product handling
- `DigestSimulator` (lines 183-276) - Enzyme digest selection
- `SnapshotExtensionSimulator` (lines 279-374) - Single-base extension
- `SnapshotValidator` (lines 377-591) - Workflow orchestrator

**Design Principles:**
- **SOLID** - Single Responsibility, Dependency Injection
- **DRY** - Configuration-based (no hardcoded values)
- **Modular** - Each component independently testable
- **Type-safe** - Full type hints

### Dependencies

```python
from Bio.Restriction import Restriction  # Enzyme database
from Bio.Seq import Seq                 # Sequence manipulation
```

**Install:**
```bash
pip install biopython
```

---

## See Also

- **[Mutation Guide](mutations.md)** - Creating mutations for validation
- **[Configuration Reference](../reference/configuration.md)** - SNaPshot config structure
- **[API Reference](../reference/api/snapshot_validator.md)** - Python API docs

---

**Algorithm Source:** `muc_one_up/analysis/snapshot_validator.py`
**Implementation:** Lines 1-591
**Tests:** `tests/test_snapshot_validator.py`
**CLI Integration:** `muc_one_up/cli/click_main.py:1235-1285`
