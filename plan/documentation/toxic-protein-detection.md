# Toxic Protein Detection

Quantitative algorithm to detect toxic protein features in ORF translations caused by MUC1 VNTR frameshift mutations.

---

## Overview

**Purpose:** Identify ORFs with toxic repeat patterns characteristic of ADTKD-MUC1 frameshift mutations.

**Method:** Combines repeat structure detection with amino acid composition analysis to generate a toxicity score.

**Output:** JSON report with scores and toxic flag (1.0 = toxic, 0.0 = normal)

---

## Quick Example

### Input
```bash
muconeup --config config.json analyze orfs \
  sample.001.mut.fa \
  --out-base toxic_analysis
```

### Output (Toxic ORF Detected)

```json
{
  "haplotype_1_ORF.9 [5063-8693](+) type:complete length:3630 frame:3 start:ATG stop:TGA": {
    "repeat_count": 19,
    "avg_repeat_identity": 0.98,
    "repeat_score": 0.98,
    "composition_similarity": 0.44,
    "overall_score": 0.77,
    "toxic_flag": 1.0
  }
}
```

**Interpretation:** `overall_score = 0.77 > 0.5` → **TOXIC** (frameshift mutation detected)

### Output (Normal ORF)

```json
{
  "haplotype_2_ORF.7 [5063-8957](+) type:complete length:3894 frame:3 start:ATG stop:TGA": {
    "repeat_count": 0,
    "avg_repeat_identity": 0.0,
    "repeat_score": 0.0,
    "composition_similarity": 0.25,
    "overall_score": 0.10,
    "toxic_flag": 0.0
  }
}
```

**Interpretation:** `overall_score = 0.10 < 0.5` → **NOT TOXIC** (wild-type or non-repetitive)

---

## How It Works

### Algorithm Overview

The algorithm combines **pattern matching** and **composition analysis** to detect toxic proteins:

```
Input Protein Sequence
         ↓
    [STEP 1: Sliding Window Scan]
    Find repeating toxic motifs
         ↓
    repeat_count, avg_identity
         ↓
    [STEP 2: Calculate Repeat Score]
    Normalize by expected count
         ↓
    repeat_score
         ↓
    [STEP 3: Composition Analysis]
    Compare R/C/H frequencies
         ↓
    composition_similarity
         ↓
    [STEP 4: Weighted Combination]
    60% repeat + 40% composition
         ↓
    overall_score
         ↓
    [DECISION: overall_score > 0.5?]
         ↓
    toxic_flag = 1.0 (TOXIC)
         or
    toxic_flag = 0.0 (NOT TOXIC)
```

---

### Step 1: Sliding Window Repeat Detection

**Concept:** Scan the protein sequence with a moving window to find repetitive toxic patterns.

**Toxic consensus motif:** `RCHLGPGHQAGPGLHR` (16 amino acids)

This is the **pathogenic repeat pattern** produced by MUC1 frameshift mutations. When a frameshift occurs (e.g., dupC insertion), the altered reading frame translates to this repetitive toxic sequence instead of the normal MUC1 protein.

**Visual representation:**

```
Protein: MTSSV...RCHLGPGHQAGPGLHR...RCHLGPGHQAGPGLHR...RCHLGPGAQAGPGLHR...
         [----------------]                                     Window 1
                [----------------]                              Window 2
                     [----------------]                         Window 3
                          ...scan entire sequence...
```

**Identity calculation** (Hamming distance):

For each window, compare to consensus:

```
Window:    R C H L G P G H Q A G P G L H R
Consensus: R C H L G P G H Q A G P G L H R
Match:     ✓ ✓ ✓ ✓ ✓ ✓ ✓ ✓ ✓ ✓ ✓ ✓ ✓ ✓ ✓ ✓  → 16/16 = 1.00 (100%)

Window:    R C H L G P G A Q A G P G L H R
Consensus: R C H L G P G H Q A G P G L H R
Match:     ✓ ✓ ✓ ✓ ✓ ✓ ✓ ✗ ✓ ✓ ✓ ✓ ✓ ✓ ✓ ✓  → 15/16 = 0.94 (94%)
```

**Formula:**

```
Identity(window) = (Number of matching amino acids) / (Window length)

If Identity ≥ 0.8 (80% threshold):
    ✓ Count this window as a toxic repeat
```

**Output:**
- `N` = repeat_count (number of windows with identity ≥ 0.8)
- `I_avg` = avg_repeat_identity (average identity across all matches)

---

### Step 2: Repeat Score Calculation

**Concept:** Normalize repeat count by expected value to get a 0-1 score.

**Formula:**

```
S_repeat = I_avg × min(N, E) / E

Where:
  N = repeat_count (detected repeats)
  E = expected_repeat_count (default: 10)
  I_avg = avg_repeat_identity
```

**Interpretation:**

- If N = 0 → S_repeat = 0.0 (no toxic repeats)
- If N ≥ E → S_repeat ≈ I_avg (saturates at expected count)
- Higher identity → higher score

**Example calculation:**

```
Given: N = 19 repeats, I_avg = 0.98, E = 10

S_repeat = 0.98 × min(19, 10) / 10
         = 0.98 × 10 / 10
         = 0.98
```

---

### Step 3: Composition Similarity

**Concept:** Measure how similar amino acid frequencies are to wild-type MUC1.

**Key residues:** R (Arginine), C (Cysteine), H (Histidine)

**Formula:**

```
S_comp = 1 - [ Σ |f_mut(aa) - f_wt(aa)| ] / [ Σ f_wt(aa) ]

Where:
  f_mut(aa) = frequency of amino acid 'aa' in mutant protein
  f_wt(aa) = frequency of amino acid 'aa' in wild-type protein
  Summation over aa ∈ {R, C, H}
```

**Interpretation:**

- S_comp = 1.0 → identical to wild-type
- S_comp = 0.0 → maximally different
- Frameshift mutations alter reading frame → different amino acid frequencies

**Visual representation:**

```
Wild-type:  R=5%, C=2%, H=3%    (normal MUC1 composition)
            [R][R][R][R][R]  [C][C]  [H][H][H]

Frameshift: R=12%, C=8%, H=6%   (toxic protein composition)
            [R][R][R][R][R][R][R][R][R][R][R][R]  [C][C][C][C][C][C][C][C]  [H][H][H][H][H][H]
            ↑ More R, C, H due to altered reading frame
```

**Example calculation:**

```
Given:
  f_mut(R) = 0.12, f_mut(C) = 0.08, f_mut(H) = 0.06
  f_wt(R) = 0.05, f_wt(C) = 0.02, f_wt(H) = 0.03

Differences:
  |0.12 - 0.05| = 0.07  (R)
  |0.08 - 0.02| = 0.06  (C)
  |0.06 - 0.03| = 0.03  (H)
  Sum = 0.16

Wildtype sum:
  0.05 + 0.02 + 0.03 = 0.10

S_comp = 1 - (0.16 / 0.10) = 1 - 1.6 = -0.6 → clipped to 0.0
```

(Note: Negative values are clipped to 0.0 minimum)

---

### Step 4: Combined Scoring

**Formula:**

```
S_overall = w_r × S_repeat + w_c × S_comp

Where:
  w_r = 0.6  (weight for repeat score)
  w_c = 0.4  (weight for composition score)
```

**Decision rule:**

```
If S_overall > 0.5:
    toxic_flag = 1.0  (TOXIC)
Else:
    toxic_flag = 0.0  (NOT TOXIC)
```

**Rationale:** Repeat pattern is more diagnostic (60% weight) than composition (40% weight).

---

### Worked Example (Toxic Protein)

**Input:** ORF with frameshift mutation

**Step 1: Detect repeats**
- Scanned 3630 amino acids
- Found N = 19 windows with identity ≥ 0.8
- Average identity I_avg = 0.98

**Step 2: Calculate repeat score**
```
S_repeat = 0.98 × min(19, 10) / 10 = 0.98
```

**Step 3: Calculate composition similarity**
```
S_comp = 0.44  (differs from wild-type)
```

**Step 4: Combine scores**
```
S_overall = 0.6 × 0.98 + 0.4 × 0.44
          = 0.588 + 0.176
          = 0.764
          ≈ 0.77
```

**Decision:**
```
0.77 > 0.5 → toxic_flag = 1.0 (TOXIC)
```

**Output:**
```json
{
  "repeat_count": 19,
  "avg_repeat_identity": 0.98,
  "repeat_score": 0.98,
  "composition_similarity": 0.44,
  "overall_score": 0.77,
  "toxic_flag": 1.0
}
```

---

### Worked Example (Normal Protein)

**Input:** Wild-type ORF

**Step 1: Detect repeats**
- Scanned 3894 amino acids
- Found N = 0 windows with identity ≥ 0.8
- Average identity I_avg = 0.0

**Step 2: Calculate repeat score**
```
S_repeat = 0.0 × min(0, 10) / 10 = 0.0
```

**Step 3: Calculate composition similarity**
```
S_comp = 0.25  (different frame, but no toxic pattern)
```

**Step 4: Combine scores**
```
S_overall = 0.6 × 0.0 + 0.4 × 0.25
          = 0.0 + 0.10
          = 0.10
```

**Decision:**
```
0.10 ≤ 0.5 → toxic_flag = 0.0 (NOT TOXIC)
```

**Output:**
```json
{
  "repeat_count": 0,
  "avg_repeat_identity": 0.0,
  "repeat_score": 0.0,
  "composition_similarity": 0.25,
  "overall_score": 0.10,
  "toxic_flag": 0.0
}
```

---

## Interpreting Scores

### Toxic ORF (Frameshift)

```json
{
  "repeat_count": 19,           // Many toxic repeats found
  "avg_repeat_identity": 0.98,  // High match to toxic pattern
  "repeat_score": 0.98,         // Strong toxic repeat signal
  "composition_similarity": 0.44, // Moderately altered composition
  "overall_score": 0.77,        // Above threshold (0.5)
  "toxic_flag": 1.0             // TOXIC
}
```

**Characteristics:**
- High `repeat_count` (>10)
- High `repeat_score` (>0.8)
- `overall_score` > 0.5
- Indicates frameshift created toxic repeat pattern

### Normal ORF (Wild-type)

```json
{
  "repeat_count": 0,            // No toxic repeats
  "avg_repeat_identity": 0.0,   // No pattern match
  "repeat_score": 0.0,          // No toxic repeat signal
  "composition_similarity": 0.25, // Different from wild-type
  "overall_score": 0.10,        // Below threshold (0.5)
  "toxic_flag": 0.0             // NOT TOXIC
}
```

**Characteristics:**
- Low `repeat_count` (0-5)
- Low `repeat_score` (<0.3)
- `overall_score` < 0.5
- Normal or non-repetitive protein

### Borderline Cases

```json
{
  "overall_score": 0.48,  // Just below threshold
  "toxic_flag": 0.0       // NOT TOXIC (but close)
}
```

Scores near 0.5 may need manual review.

---

## Algorithm Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `consensus` | `RCHLGPGHQAGPGLHR` | Toxic repeat motif pattern |
| `identity_threshold` | 0.8 | Minimum match for repeat (80%) |
| `key_residues` | `["R", "C", "H"]` | Amino acids for composition |
| `expected_repeat_count` | 10 | Expected repeats in toxic protein |
| `w_repeat` | 0.6 | Weight for repeat score |
| `w_composition` | 0.4 | Weight for composition score |
| `toxic_cutoff` | 0.5 | Threshold for toxic flag |

---

## Usage

### Command Line

```bash
# Analyze ORFs and detect toxic proteins
muconeup --config config.json analyze orfs \
  sample.001.simulated.fa \
  --out-base analysis \
  --orf-min-aa 100

# View results
cat analysis.orf_stats.txt
```

### Output Files

```
analysis.pep.fa           # Predicted peptides (FASTA)
analysis.orf_stats.txt    # Summary statistics
analysis_toxic_scores.json # Detailed scores (if saved)
```

### Text Summary Format

```
Haplotype Statistics:

haplotype_1:
  Total ORFs: 12
  Toxic ORFs: 3
  Longest ORF: 1210 aa
  Toxic Features:
    - ORF.9: overall_score=0.77 (TOXIC)
    - ORF.11: overall_score=0.65 (TOXIC)
    - ORF.15: overall_score=0.52 (TOXIC)

haplotype_2:
  Total ORFs: 10
  Toxic ORFs: 0
```

### Python API

```python
from muc_one_up.toxic_protein_detector import detect_toxic_protein_in_sequence

# Analyze single ORF
result = detect_toxic_protein_in_sequence(
    protein_seq="MTSSV...RCHLGPGHQAGPGLHR...",
    consensus="RCHLGPGHQAGPGLHR",
    identity_threshold=0.8,
    expected_repeat_count=10,
    toxic_detection_cutoff=0.5
)

# Check results
if result['toxic_flag'] == 1.0:
    print(f"TOXIC: score={result['overall_score']:.2f}")
    print(f"  Repeats: {result['repeat_count']}")
    print(f"  Identity: {result['avg_repeat_identity']:.2f}")
else:
    print(f"NOT TOXIC: score={result['overall_score']:.2f}")
```

---

## Biological Context

### ADTKD-MUC1

**Disease:** Autosomal Dominant Tubulointerstitial Kidney Disease caused by MUC1 mutations

**Mechanism:**
- Frameshift mutations in MUC1 VNTR region
- Altered reading frame produces toxic protein
- Toxic protein accumulates in kidney cells
- Leads to kidney dysfunction

### Why This Algorithm Works

**Frameshift creates repetitive toxic pattern:**

The dupC mutation (and similar frameshifts) shifts the reading frame by +1 nucleotide, creating a completely different amino acid sequence:

- **Wild-type translation:** `...STAPPA HGVTSA PDTRPA...` (normal MUC1 VNTR repeats)
- **Frameshift translation:** `...RCHLGP GHQAGP GLHRRC HLGPGH QAGPGL HR...` (toxic repeating motif)

The algorithm detects this pathogenic pattern:
- Scans for `RCHLGPGHQAGPGLHR` repeating motifs
- High repeat count (>10) → high repeat_score
- Overall_score > 0.5 → **flagged as TOXIC**

**Wild-type or truncated proteins:**
- No `RCHLGPGHQAGPGLHR` pattern present
- `repeat_count = 0`
- `overall_score < 0.5`
- **NOT flagged as toxic**

---

## Troubleshooting

### All ORFs Scoring 0.0

**Cause:** ORF sequences don't contain toxic consensus pattern

**Solution:**
- Normal if analyzing wild-type samples
- Check if ORF sequences are correct (use `--orf-min-aa` to filter)

### Unexpected Toxic Flags

**Cause:** Consensus motif may not match your mutation type

**Solution:**
- Adjust threshold via Python API for custom detection parameters
- Contact developers for additional consensus motif patterns

### Low Composition Similarity

**Cause:** Amino acid composition differs from wild-type model

**Note:** This is expected for frameshifts - focus on `overall_score` and `toxic_flag`

---

## Limitations

1. **Single consensus pattern** - Default pattern may not detect all toxic mechanisms
2. **Threshold sensitivity** - Fixed at 0.5, may need calibration for specific cases
3. **No structural analysis** - Only sequence-based, doesn't predict protein folding

---

## Implementation

**Source:** `muc_one_up/toxic_protein_detector.py`
- `detect_toxic_protein_in_sequence()` - Main analysis function (lines 145-249)
- `sliding_window_repeat_analysis()` - Repeat detection (lines 63-96)
- `composition_similarity()` - Composition scoring (lines 120-142)

**Tests:** `tests/test_toxic_protein_detector.py`

**CLI:** `muc_one_up/cli/click_main.py:790-952` (analyze orfs command)

---

## See Also

- **[ORF Prediction](orfs.md)** - Complete ORF analysis workflow
- **[Mutation Guide](mutations.md)** - Applying mutations to VNTR
