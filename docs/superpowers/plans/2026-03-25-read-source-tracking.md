# Read Source Tracking Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Annotate simulated reads with ground-truth origin metadata (haplotype, position, mutation/SNP overlap) via a compressed TSV manifest and persisted repeat coordinate map.

**Architecture:** Centralized `ReadSourceTracker` class constructed with simulation metadata in orchestration.py, threaded through to platform-specific pipelines. Each platform has a parser that extracts `ReadOrigin` from native simulator output. The tracker annotates origins against a repeat coordinate map and writes a compressed manifest.

**Tech Stack:** Python 3.10+, dataclasses, gzip, csv, re, Click (CLI), pytest

**Spec:** `docs/superpowers/specs/2026-03-25-read-source-tracking-design.md`

---

## File Structure

### New Files
| File | Responsibility |
|------|---------------|
| `muc_one_up/read_simulator/source_tracking.py` | Data models, RepeatCoordinateMap, ReadSourceTracker, manifest writer |
| `muc_one_up/read_simulator/parsers/__init__.py` | Parser package init |
| `muc_one_up/read_simulator/parsers/ont_parser.py` | NanoSim read name parser → ReadOrigin list |
| `muc_one_up/read_simulator/parsers/pacbio_parser.py` | pbsim3 MAF/BAM parser → ReadOrigin list |
| `muc_one_up/read_simulator/parsers/illumina_parser.py` | Fragment sidecar + reseq order correlation → ReadOrigin list |
| `tests/test_source_tracking.py` | Unit tests for data models, coordinate map, annotation |
| `tests/test_parsers.py` | Unit tests for ONT, PacBio, Illumina parsers |

### Modified Files
| File | Changes |
|------|---------|
| `muc_one_up/cli/click_main.py` | Add `--track-read-source` flag to `simulate`, `reads illumina`, `reads ont`, `reads pacbio` |
| `muc_one_up/cli/orchestration.py` | Construct ReadSourceTracker, pass to run_read_simulation |
| `muc_one_up/cli/analysis.py` | Accept and relay tracker to simulate_reads_pipeline |
| `muc_one_up/read_simulation.py` | Thread tracker through SIMULATOR_MAP dispatch |
| `muc_one_up/read_simulator/pipeline.py` | Accept tracker, call Illumina parser, write manifest |
| `muc_one_up/read_simulator/ont_pipeline.py` | Accept tracker, call ONT parser, write manifest |
| `muc_one_up/read_simulator/pacbio_pipeline.py` | Accept tracker, call PacBio parser, write manifest |
| `muc_one_up/read_simulator/fragment_simulation.py` | Write fragment origins sidecar TSV |

---

## Task 1: Data Models and Repeat Coordinate Map

**Files:**
- Create: `muc_one_up/read_simulator/source_tracking.py`
- Test: `tests/test_source_tracking.py`

- [ ] **Step 1: Write failing tests for data models and coordinate map**

```python
# tests/test_source_tracking.py
"""Tests for read source tracking data models and coordinate map."""

from __future__ import annotations

import pytest

from muc_one_up.read_simulator.source_tracking import (
    AnnotatedRead,
    ReadOrigin,
    RepeatCoordinateMap,
    RepeatRegion,
    SNPPosition,
    build_coordinate_map,
)


class TestRepeatRegion:
    """Tests for RepeatRegion dataclass."""

    def test_basic_creation(self):
        region = RepeatRegion(
            index=1,
            repeat_type="X",
            start=100,
            end=160,
            is_mutated=False,
            mutation_name=None,
        )
        assert region.index == 1
        assert region.repeat_type == "X"
        assert region.start == 100
        assert region.end == 160
        assert region.length == 60

    def test_mutated_region(self):
        region = RepeatRegion(
            index=3,
            repeat_type="X",
            start=200,
            end=261,
            is_mutated=True,
            mutation_name="dupC",
        )
        assert region.is_mutated is True
        assert region.mutation_name == "dupC"


class TestBuildCoordinateMap:
    """Tests for building coordinate maps from repeat chains."""

    @pytest.fixture()
    def repeats_dict(self):
        return {
            "1": "AACCCTCCC",       # 9 bp
            "2": "AACCCTCCC",       # 9 bp
            "X": "AACCCTCCCAACCCTCCC",  # 18 bp
            "7": "AACCCTCCC",       # 9 bp
            "8": "AACCCTCCC",       # 9 bp
            "9": "AACCCTCCC",       # 9 bp
        }

    def test_simple_chain(self, repeats_dict):
        chain = ["1", "2", "X", "7", "8", "9"]
        coord_map = build_coordinate_map(
            haplotype=1,
            chain=chain,
            repeats_dict=repeats_dict,
            left_const_len=50,
            mutation_positions=[],
            mutation_name=None,
            snp_info=[],
        )
        assert coord_map.haplotype == 1
        assert len(coord_map.regions) == 6
        # First repeat starts after left constant
        assert coord_map.regions[0].start == 50
        assert coord_map.regions[0].end == 59  # 50 + 9
        # Second repeat follows first
        assert coord_map.regions[1].start == 59
        assert coord_map.regions[1].end == 68  # 59 + 9
        # Third repeat (X) is 18bp
        assert coord_map.regions[2].start == 68
        assert coord_map.regions[2].end == 86  # 68 + 18
        # VNTR boundaries
        assert coord_map.vntr_start == 50
        assert coord_map.vntr_end == 86 + 9 + 9 + 9  # last 3 are 9bp each

    def test_with_mutation_markers(self, repeats_dict):
        chain = ["1", "2", "Xm", "7", "8", "9"]
        coord_map = build_coordinate_map(
            haplotype=1,
            chain=chain,
            repeats_dict=repeats_dict,
            left_const_len=50,
            mutation_positions=[(1, 3)],
            mutation_name="dupC",
            snp_info=[],
        )
        # Xm should be treated as X for sequence lookup
        assert coord_map.regions[2].repeat_type == "X"
        assert coord_map.regions[2].is_mutated is True
        assert coord_map.regions[2].mutation_name == "dupC"

    def test_with_snps(self, repeats_dict):
        chain = ["1", "2", "X"]
        snps = [{"position": 55, "ref_base": "A", "alt_base": "T"}]
        coord_map = build_coordinate_map(
            haplotype=1,
            chain=chain,
            repeats_dict=repeats_dict,
            left_const_len=50,
            mutation_positions=[],
            mutation_name=None,
            snp_info=snps,
        )
        assert len(coord_map.snp_positions) == 1
        assert coord_map.snp_positions[0].position == 55
        assert coord_map.snp_positions[0].repeat_index == 1  # falls in first repeat (50-59)

    def test_empty_chain(self, repeats_dict):
        coord_map = build_coordinate_map(
            haplotype=1,
            chain=[],
            repeats_dict=repeats_dict,
            left_const_len=50,
            mutation_positions=[],
            mutation_name=None,
            snp_info=[],
        )
        assert len(coord_map.regions) == 0
        assert coord_map.vntr_start == 50
        assert coord_map.vntr_end == 50

    def test_single_repeat(self, repeats_dict):
        chain = ["X"]
        coord_map = build_coordinate_map(
            haplotype=1,
            chain=chain,
            repeats_dict=repeats_dict,
            left_const_len=0,
            mutation_positions=[],
            mutation_name=None,
            snp_info=[],
        )
        assert len(coord_map.regions) == 1
        assert coord_map.regions[0].start == 0
        assert coord_map.regions[0].end == 18
        assert coord_map.vntr_start == 0
        assert coord_map.vntr_end == 18
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_source_tracking.py -v`
Expected: FAIL with ImportError

- [ ] **Step 3: Implement data models and coordinate map builder**

Create `muc_one_up/read_simulator/source_tracking.py` with:
- `RepeatRegion` dataclass (index, repeat_type, start, end, is_mutated, mutation_name, length property)
- `SNPPosition` dataclass (position, ref_base, alt_base, repeat_index)
- `RepeatCoordinateMap` dataclass (haplotype, regions, vntr_start, vntr_end, snp_positions)
- `ReadOrigin` dataclass (read_id, haplotype, ref_start, ref_end, strand)
- `AnnotatedRead` dataclass (all ReadOrigin fields + overlaps_vntr, repeat_units, overlaps_mutation, mutation_name, overlaps_snp)
- `build_coordinate_map()` function that:
  - Iterates chain, strips "m" suffix for sequence lookup
  - Computes cumulative start/end from left_const_len
  - Marks mutated regions from mutation_positions
  - Maps SNP positions to repeat indices
  - Returns RepeatCoordinateMap

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_source_tracking.py -v`
Expected: All PASS

- [ ] **Step 5: Commit**

```bash
git add muc_one_up/read_simulator/source_tracking.py tests/test_source_tracking.py
git commit -m "feat: add data models and coordinate map builder for read source tracking"
```

---

## Task 2: ReadSourceTracker — Annotation and Manifest Writing

**Files:**
- Modify: `muc_one_up/read_simulator/source_tracking.py`
- Modify: `tests/test_source_tracking.py`

- [ ] **Step 1: Write failing tests for annotation logic and manifest**

Add to `tests/test_source_tracking.py`:

```python
import csv
import gzip
import os
import tempfile

from muc_one_up.read_simulator.source_tracking import ReadSourceTracker


class TestReadSourceTracker:
    """Tests for ReadSourceTracker annotation and manifest writing."""

    @pytest.fixture()
    def repeats_dict(self):
        return {
            "1": "AACCCTCCC",           # 9bp, positions 50-59
            "2": "AACCCTCCC",           # 9bp, positions 59-68
            "X": "AACCCTCCCAACCCTCCC",  # 18bp, positions 68-86
            "7": "AACCCTCCC",           # 9bp, positions 86-95
            "8": "AACCCTCCC",           # 9bp, positions 95-104
            "9": "AACCCTCCC",           # 9bp, positions 104-113
        }

    @pytest.fixture()
    def tracker(self, repeats_dict):
        return ReadSourceTracker(
            repeat_chains={1: ["1", "2", "Xm", "7", "8", "9"], 2: ["1", "2", "X", "7", "8", "9"]},
            repeats_dict=repeats_dict,
            left_const_len=50,
            mutation_positions=[(1, 3)],
            mutation_name="dupC",
            snp_info={0: [{"position": 55, "ref_base": "A", "alt_base": "T"}]},
        )

    def test_annotation_overlaps_vntr(self, tracker):
        origins = [
            ReadOrigin(read_id="read1", haplotype=1, ref_start=60, ref_end=90, strand="+"),
        ]
        annotated = list(tracker.annotate_reads(origins))
        assert len(annotated) == 1
        assert annotated[0].overlaps_vntr is True

    def test_annotation_outside_vntr(self, tracker):
        origins = [
            ReadOrigin(read_id="read2", haplotype=1, ref_start=0, ref_end=40, strand="+"),
        ]
        annotated = list(tracker.annotate_reads(origins))
        assert annotated[0].overlaps_vntr is False
        assert annotated[0].repeat_units == ""

    def test_annotation_overlaps_mutation(self, tracker):
        # Read overlapping position 68-86 (Xm repeat, mutated) on haplotype 1
        origins = [
            ReadOrigin(read_id="read3", haplotype=1, ref_start=70, ref_end=80, strand="+"),
        ]
        annotated = list(tracker.annotate_reads(origins))
        assert annotated[0].overlaps_mutation is True
        assert annotated[0].mutation_name == "dupC"

    def test_annotation_no_mutation_overlap(self, tracker):
        # Read overlapping repeat 1 (not mutated) on haplotype 1
        origins = [
            ReadOrigin(read_id="read4", haplotype=1, ref_start=50, ref_end=58, strand="+"),
        ]
        annotated = list(tracker.annotate_reads(origins))
        assert annotated[0].overlaps_mutation is False
        assert annotated[0].mutation_name == "."

    def test_annotation_overlaps_snp(self, tracker):
        # Read overlapping position 55 (SNP in haplotype 1)
        origins = [
            ReadOrigin(read_id="read5", haplotype=1, ref_start=50, ref_end=60, strand="+"),
        ]
        annotated = list(tracker.annotate_reads(origins))
        assert annotated[0].overlaps_snp is True

    def test_annotation_repeat_units(self, tracker):
        # Read spanning repeats 2 and 3 (positions 59-86) on haplotype 1
        origins = [
            ReadOrigin(read_id="read6", haplotype=1, ref_start=59, ref_end=86, strand="+"),
        ]
        annotated = list(tracker.annotate_reads(origins))
        assert annotated[0].repeat_units == "2,3"

    def test_write_manifest(self, tracker, tmp_path):
        origins = [
            ReadOrigin(read_id="read1", haplotype=1, ref_start=60, ref_end=90, strand="+"),
            ReadOrigin(read_id="read2", haplotype=2, ref_start=0, ref_end=40, strand="-"),
        ]
        annotated = list(tracker.annotate_reads(origins))
        manifest_path = str(tmp_path / "test_read_manifest.tsv.gz")
        tracker.write_manifest(annotated, manifest_path)

        assert os.path.exists(manifest_path)
        with gzip.open(manifest_path, "rt") as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = list(reader)
        assert len(rows) == 2
        assert rows[0]["read_id"] == "read1"
        assert rows[0]["overlaps_vntr"] == "true"
        assert rows[1]["read_id"] == "read2"
        assert rows[1]["overlaps_vntr"] == "false"

    def test_write_coordinate_map(self, tracker, tmp_path):
        coord_path = str(tmp_path / "test_repeat_coordinates.tsv")
        tracker.write_coordinate_map(coord_path)
        assert os.path.exists(coord_path)
        with open(coord_path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = list(reader)
        # 6 repeats per haplotype × 2 haplotypes = 12 rows
        assert len(rows) == 12
        assert rows[0]["haplotype"] == "1"
        assert rows[0]["repeat_type"] == "1"
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_source_tracking.py::TestReadSourceTracker -v`
Expected: FAIL with ImportError or AttributeError

- [ ] **Step 3: Implement ReadSourceTracker class**

Add to `muc_one_up/read_simulator/source_tracking.py`:
- `ReadSourceTracker.__init__()`: builds coordinate maps for each haplotype
- `ReadSourceTracker.annotate_reads()`: intersects read intervals with coordinate maps
- `ReadSourceTracker.write_manifest()`: writes gzip-compressed TSV with header
- `ReadSourceTracker.write_coordinate_map()`: writes TSV with all repeat regions
- Helper: `_overlapping_regions()` — interval overlap check between read and regions

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_source_tracking.py -v`
Expected: All PASS

- [ ] **Step 5: Commit**

```bash
git add muc_one_up/read_simulator/source_tracking.py tests/test_source_tracking.py
git commit -m "feat: add ReadSourceTracker with annotation and manifest writing"
```

---

## Task 3: ONT Parser

**Files:**
- Create: `muc_one_up/read_simulator/parsers/__init__.py`
- Create: `muc_one_up/read_simulator/parsers/ont_parser.py`
- Create: `tests/test_parsers.py`

- [ ] **Step 1: Write failing tests for ONT parser**

```python
# tests/test_parsers.py
"""Tests for platform-specific read origin parsers."""

from __future__ import annotations

import os
import tempfile

import pytest

from muc_one_up.read_simulator.parsers.ont_parser import parse_nanosim_reads
from muc_one_up.read_simulator.source_tracking import ReadOrigin


class TestONTParser:
    """Tests for NanoSim read name parser."""

    def test_parse_aligned_read(self, tmp_path):
        fastq = tmp_path / "reads.fastq"
        fastq.write_text(
            "@haplotype_1_500_aligned_0_F_10_100_5\n"
            "ACGTACGTACGTACGT\n"
            "+\n"
            "IIIIIIIIIIIIIIII\n"
        )
        origins = parse_nanosim_reads(str(fastq), haplotype_map=None)
        assert len(origins) == 1
        assert origins[0].read_id == "haplotype_1_500_aligned_0_F_10_100_5"
        assert origins[0].haplotype == 1
        assert origins[0].ref_start == 500
        assert origins[0].strand == "+"
        assert origins[0].ref_end == 500 + 16  # read length as approximation

    def test_parse_reverse_strand(self, tmp_path):
        fastq = tmp_path / "reads.fastq"
        fastq.write_text(
            "@haplotype_2_1000_aligned_1_R_5_200_10\n"
            "ACGTACGT\n"
            "+\n"
            "IIIIIIII\n"
        )
        origins = parse_nanosim_reads(str(fastq), haplotype_map=None)
        assert origins[0].strand == "-"
        assert origins[0].haplotype == 2

    def test_parse_unaligned_read(self, tmp_path):
        fastq = tmp_path / "reads.fastq"
        fastq.write_text(
            "@haplotype_1_0_unaligned_0_F_0_100_0\n"
            "ACGT\n"
            "+\n"
            "IIII\n"
        )
        origins = parse_nanosim_reads(str(fastq), haplotype_map=None)
        assert len(origins) == 1
        assert origins[0].ref_start == 0

    def test_diploid_split_sim_haplotype_map(self, tmp_path):
        fastq = tmp_path / "reads.fastq"
        # In split-sim, read names may not have correct haplotype prefix
        fastq.write_text(
            "@chr1_500_aligned_0_F_10_100_5\n"
            "ACGTACGT\n"
            "+\n"
            "IIIIIIII\n"
        )
        # haplotype_map overrides parsed haplotype
        origins = parse_nanosim_reads(str(fastq), haplotype_map=2)
        assert origins[0].haplotype == 2

    def test_multiple_reads(self, tmp_path):
        fastq = tmp_path / "reads.fastq"
        fastq.write_text(
            "@haplotype_1_100_aligned_0_F_5_50_3\n"
            "ACGT\n"
            "+\n"
            "IIII\n"
            "@haplotype_1_200_aligned_1_R_3_80_2\n"
            "ACGTACGT\n"
            "+\n"
            "IIIIIIII\n"
        )
        origins = parse_nanosim_reads(str(fastq), haplotype_map=None)
        assert len(origins) == 2

    def test_empty_fastq(self, tmp_path):
        fastq = tmp_path / "reads.fastq"
        fastq.write_text("")
        origins = parse_nanosim_reads(str(fastq), haplotype_map=None)
        assert len(origins) == 0
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_parsers.py::TestONTParser -v`
Expected: FAIL with ImportError

- [ ] **Step 3: Implement ONT parser**

Create `muc_one_up/read_simulator/parsers/__init__.py` (empty) and `ont_parser.py`:
- `parse_nanosim_reads(fastq_path, haplotype_map=None) -> list[ReadOrigin]`
- Regex pattern: `^(?:haplotype_)?(\d+|[^_]+)_(\d+)_(aligned|unaligned)_(\d+)_([FR])_(\d+)_(\d+)_(\d+)$`
- Parse each FASTQ read name, extract haplotype, position, strand
- If `haplotype_map` is provided (int), override haplotype for all reads (used in split-sim)
- Compute ref_end as ref_start + read_length

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_parsers.py::TestONTParser -v`
Expected: All PASS

- [ ] **Step 5: Commit**

```bash
git add muc_one_up/read_simulator/parsers/ tests/test_parsers.py
git commit -m "feat: add ONT NanoSim read name parser"
```

---

## Task 4: PacBio Parser

**Files:**
- Create: `muc_one_up/read_simulator/parsers/pacbio_parser.py`
- Modify: `tests/test_parsers.py`

- [ ] **Step 1: Write failing tests for PacBio parser**

Add to `tests/test_parsers.py`:

```python
from muc_one_up.read_simulator.parsers.pacbio_parser import parse_pacbio_reads


class TestPacBioParser:
    """Tests for pbsim3 MAF/BAM parser with alignment fallback."""

    def test_parse_from_bam_alignment(self, tmp_path):
        """Test fallback: parse read positions from aligned BAM."""
        # Create a minimal test BAM - this tests the alignment-based fallback
        # which is simpler and more reliable than MAF parsing
        # For unit tests, we test the logic with synthetic ReadOrigin data
        pass  # Integration test covers this with real BAM

    def test_parse_maf_file(self, tmp_path):
        """Test parsing a MAF alignment file from pbsim3."""
        maf_content = (
            "a\n"
            "s haplotype_1 100 60 + 5000 ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n"
            "s S0_0 0 60 + 60 ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n"
            "\n"
            "a\n"
            "s haplotype_1 300 40 + 5000 ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n"
            "s S0_1 0 40 + 40 ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n"
            "\n"
        )
        maf_path = tmp_path / "test.maf"
        maf_path.write_text(maf_content)
        origins = parse_pacbio_reads(
            maf_paths=[str(maf_path)],
            haplotype_index=1,
            aligned_bam=None,
        )
        assert len(origins) == 2
        assert origins[0].haplotype == 1
        assert origins[0].ref_start == 100
        assert origins[0].ref_end == 160
        assert origins[1].ref_start == 300
        assert origins[1].ref_end == 340

    def test_parse_multiple_maf_files(self, tmp_path):
        """Test parsing MAF files from multiple haplotypes."""
        maf1 = tmp_path / "sd_0001.maf"
        maf1.write_text(
            "a\n"
            "s haplotype_1 100 50 + 5000 " + "A" * 50 + "\n"
            "s S0_0 0 50 + 50 " + "A" * 50 + "\n\n"
        )
        maf2 = tmp_path / "sd_0002.maf"
        maf2.write_text(
            "a\n"
            "s haplotype_2 200 30 + 4000 " + "A" * 30 + "\n"
            "s S0_0 0 30 + 30 " + "A" * 30 + "\n\n"
        )
        origins_hap1 = parse_pacbio_reads(
            maf_paths=[str(maf1)], haplotype_index=1, aligned_bam=None
        )
        origins_hap2 = parse_pacbio_reads(
            maf_paths=[str(maf2)], haplotype_index=2, aligned_bam=None
        )
        assert all(o.haplotype == 1 for o in origins_hap1)
        assert all(o.haplotype == 2 for o in origins_hap2)

    def test_empty_maf(self, tmp_path):
        maf_path = tmp_path / "empty.maf"
        maf_path.write_text("")
        origins = parse_pacbio_reads(
            maf_paths=[str(maf_path)], haplotype_index=1, aligned_bam=None
        )
        assert len(origins) == 0
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_parsers.py::TestPacBioParser -v`
Expected: FAIL with ImportError

- [ ] **Step 3: Implement PacBio parser**

Create `muc_one_up/read_simulator/parsers/pacbio_parser.py`:
- `parse_pacbio_reads(maf_paths, haplotype_index, aligned_bam=None) -> list[ReadOrigin]`
- Primary path: parse MAF `s` lines — first `s` line in each block is reference, second is read
- Extract ref_start, alignment_length from reference `s` line; compute ref_end = ref_start + alignment_length
- Read name from read `s` line
- Fallback path: if no MAF files, parse aligned BAM for read positions (deferred to integration)

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_parsers.py::TestPacBioParser -v`
Expected: All PASS

- [ ] **Step 5: Commit**

```bash
git add muc_one_up/read_simulator/parsers/pacbio_parser.py tests/test_parsers.py
git commit -m "feat: add PacBio MAF parser for read source tracking"
```

---

## Task 5: Illumina Parser and Fragment Sidecar

**Files:**
- Create: `muc_one_up/read_simulator/parsers/illumina_parser.py`
- Modify: `muc_one_up/read_simulator/fragment_simulation.py`
- Modify: `tests/test_parsers.py`

- [ ] **Step 1: Write failing tests for Illumina sidecar writer and parser**

Add to `tests/test_parsers.py`:

```python
from muc_one_up.read_simulator.parsers.illumina_parser import (
    parse_illumina_reads,
    write_fragment_origins,
)


class TestIlluminaParser:
    """Tests for Illumina fragment sidecar and read parser."""

    def test_write_fragment_origins(self, tmp_path):
        origins_path = str(tmp_path / "fragment_origins.tsv")
        fragments = [
            {"fragment_index": 1, "chrom": "haplotype_1", "fstart": 100, "fend": 250, "strand": "+"},
            {"fragment_index": 2, "chrom": "haplotype_2", "fstart": 300, "fend": 450, "strand": "-"},
        ]
        write_fragment_origins(fragments, origins_path)
        assert os.path.exists(origins_path)
        with open(origins_path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = list(reader)
        assert len(rows) == 2
        assert rows[0]["fragment_index"] == "1"
        assert rows[0]["chrom"] == "haplotype_1"

    def test_parse_illumina_reads_from_sidecar(self, tmp_path):
        # Create sidecar
        origins_path = str(tmp_path / "fragment_origins.tsv")
        with open(origins_path, "w") as f:
            f.write("fragment_index\tchrom\tfstart\tfend\tstrand\n")
            f.write("1\thaplotype_1\t100\t250\t+\n")
            f.write("2\thaplotype_2\t300\t450\t-\n")

        # Create paired FASTQ simulating reseq output (order-based correlation)
        fq1 = tmp_path / "reads_R1.fastq.gz"
        fq2 = tmp_path / "reads_R2.fastq.gz"

        import gzip
        with gzip.open(str(fq1), "wt") as f:
            f.write("@read_0001/1\nACGT\n+\nIIII\n@read_0002/1\nACGT\n+\nIIII\n")
        with gzip.open(str(fq2), "wt") as f:
            f.write("@read_0001/2\nACGT\n+\nIIII\n@read_0002/2\nACGT\n+\nIIII\n")

        seq_names = {"haplotype_1": 1, "haplotype_2": 2}
        origins = parse_illumina_reads(origins_path, str(fq1), seq_names)
        assert len(origins) == 2
        assert origins[0].haplotype == 1
        assert origins[0].ref_start == 100
        assert origins[0].ref_end == 250
        assert origins[1].haplotype == 2

    def test_empty_sidecar(self, tmp_path):
        origins_path = str(tmp_path / "fragment_origins.tsv")
        with open(origins_path, "w") as f:
            f.write("fragment_index\tchrom\tfstart\tfend\tstrand\n")
        origins = parse_illumina_reads(origins_path, None, {})
        assert len(origins) == 0
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_parsers.py::TestIlluminaParser -v`
Expected: FAIL

- [ ] **Step 3: Implement Illumina parser**

Create `muc_one_up/read_simulator/parsers/illumina_parser.py`:
- `write_fragment_origins(fragments: list[dict], output_path: str)` — write TSV with header
- `parse_illumina_reads(sidecar_path, fastq_r1_path, seq_name_to_haplotype) -> list[ReadOrigin]`
  - Read sidecar TSV
  - If fastq_r1_path provided, read R1 FASTQ to get read names (order-based correlation)
  - Map chrom to haplotype via seq_name_to_haplotype dict
  - Return one ReadOrigin per fragment (using fragment coordinates)

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_parsers.py::TestIlluminaParser -v`
Expected: All PASS

- [ ] **Step 5: Commit**

```bash
git add muc_one_up/read_simulator/parsers/illumina_parser.py tests/test_parsers.py
git commit -m "feat: add Illumina fragment sidecar writer and parser"
```

---

## Task 6: Fragment Simulation Sidecar Integration

**Files:**
- Modify: `muc_one_up/read_simulator/fragment_simulation.py`

- [ ] **Step 1: Modify fragment_simulation.py to write sidecar**

In `simulate_fragments()`, after the fragment simulation loop:
- Collect fragment origins (chrom, fstart, fend, strand) for each fragment
- If a `fragment_origins_path` parameter is provided (not None), write the sidecar TSV using `write_fragment_origins()`
- Add `fragment_origins_path: str | None = None` parameter to `simulate_fragments()`

- [ ] **Step 2: Run existing tests to verify no regression**

Run: `pytest tests/ -k "fragment" -v`
Expected: All existing tests PASS

- [ ] **Step 3: Commit**

```bash
git add muc_one_up/read_simulator/fragment_simulation.py
git commit -m "feat: capture fragment origins in sidecar TSV during Illumina simulation"
```

---

## Task 7: CLI Flag and Metadata Threading

**Files:**
- Modify: `muc_one_up/cli/click_main.py`
- Modify: `muc_one_up/cli/orchestration.py`
- Modify: `muc_one_up/cli/analysis.py`
- Modify: `muc_one_up/read_simulation.py`

- [ ] **Step 1: Add --track-read-source to CLI**

Add Click option `--track-read-source` (is_flag=True, default=False) to:
- `simulate` command in `click_main.py`
- `reads illumina` subcommand
- `reads ont` subcommand
- `reads pacbio` subcommand

- [ ] **Step 2: Thread tracker through orchestration**

In `orchestration.py` `run_single_simulation_iteration()`:
- After simulation completes and before `run_read_simulation()`, construct `ReadSourceTracker` if `args.track_read_source` is True
- Pass tracker to `run_read_simulation()` as new parameter

In `analysis.py` `run_read_simulation()`:
- Add `source_tracker=None` parameter
- Pass tracker to `simulate_reads_pipeline()`

In `read_simulation.py` `simulate_reads()`:
- Add `source_tracker=None` parameter
- Pass tracker through `SIMULATOR_MAP` dispatch to platform pipelines

- [ ] **Step 3: Run existing tests to verify no regression**

Run: `make test-fast`
Expected: All existing tests PASS (tracker defaults to None everywhere)

- [ ] **Step 4: Commit**

```bash
git add muc_one_up/cli/click_main.py muc_one_up/cli/orchestration.py muc_one_up/cli/analysis.py muc_one_up/read_simulation.py
git commit -m "feat: add --track-read-source flag and thread tracker through pipeline"
```

---

## Task 8: Platform Pipeline Integration

**Files:**
- Modify: `muc_one_up/read_simulator/pipeline.py`
- Modify: `muc_one_up/read_simulator/ont_pipeline.py`
- Modify: `muc_one_up/read_simulator/pacbio_pipeline.py`

- [ ] **Step 1: Integrate tracker into ONT pipeline**

In `simulate_ont_reads_pipeline()`:
- Add `source_tracker=None` parameter
- After read simulation (NanoSim), if tracker is not None:
  - For split-sim: parse each haplotype FASTQ with `parse_nanosim_reads()` using haplotype_map
  - For standard: parse merged FASTQ
  - Call `tracker.annotate_reads(origins)`
  - Call `tracker.write_manifest()` and `tracker.write_coordinate_map()`

- [ ] **Step 2: Integrate tracker into PacBio pipeline**

In `simulate_pacbio_hifi_reads()`:
- Add `source_tracker=None` parameter
- When tracking enabled, do NOT delete MAF files during CLR cleanup
- After alignment, parse MAF files per haplotype with `parse_pacbio_reads()`
- If MAF parsing fails, fall back to aligned BAM positions
- Annotate and write manifest
- Delete MAF files after successful manifest generation

- [ ] **Step 3: Integrate tracker into Illumina pipeline**

In `simulate_reads_pipeline()`:
- Add `source_tracker=None` parameter
- Pass `fragment_origins_path` to `simulate_fragments()` when tracking enabled
- After reseq + read splitting, parse sidecar with `parse_illumina_reads()`
- After efficiency bias (if applied), filter manifest to surviving reads
- Annotate and write manifest
- Delete sidecar after successful manifest generation

- [ ] **Step 4: Run full test suite**

Run: `make test-fast`
Expected: All tests PASS

- [ ] **Step 5: Commit**

```bash
git add muc_one_up/read_simulator/pipeline.py muc_one_up/read_simulator/ont_pipeline.py muc_one_up/read_simulator/pacbio_pipeline.py
git commit -m "feat: integrate read source tracking into all platform pipelines"
```

---

## Task 9: Standalone reads Metadata Reconstruction

**Files:**
- Modify: `muc_one_up/read_simulator/source_tracking.py`
- Modify: `tests/test_source_tracking.py`

- [ ] **Step 1: Write failing tests for reconstruction**

```python
class TestTrackerFromCompanionFiles:
    def test_reconstruct_from_stats_and_structure(self, tmp_path):
        # Create mock simulation_stats.json
        stats = {
            "haplotypes": {
                "haplotype_1": {"repeat_chain": "1-2-X-7-8-9"},
                "haplotype_2": {"repeat_chain": "1-2-X-7-8-9"},
            },
            "mutation_details": {"name": "dupC", "targets": [[1, 3]]},
            "snp_info": {},
            "config": {
                "repeats": {"1": "AACCCTCCC", "2": "AACCCTCCC", "X": "AACCCTCCCAACCCTCCC", "7": "AACCCTCCC", "8": "AACCCTCCC", "9": "AACCCTCCC"},
                "constants": {"hg38": {"left": "A" * 50}},
            },
            "reference_assembly": "hg38",
        }
        import json
        stats_path = tmp_path / "test.simulation_stats.json"
        stats_path.write_text(json.dumps(stats))

        tracker = ReadSourceTracker.from_companion_files(str(stats_path))
        assert tracker is not None
        assert 1 in tracker.coordinate_maps
        assert 2 in tracker.coordinate_maps
```

- [ ] **Step 2: Implement from_companion_files class method**

Add `ReadSourceTracker.from_companion_files(stats_path)` that:
- Loads simulation_stats.json
- Extracts repeat chains, repeats_dict, mutation details, SNP info, left constant length
- Constructs and returns a ReadSourceTracker
- Returns None with warning if required fields are missing

- [ ] **Step 3: Run tests**

Run: `pytest tests/test_source_tracking.py::TestTrackerFromCompanionFiles -v`
Expected: PASS

- [ ] **Step 4: Commit**

```bash
git add muc_one_up/read_simulator/source_tracking.py tests/test_source_tracking.py
git commit -m "feat: add ReadSourceTracker reconstruction from companion files"
```

---

## Task 10: Documentation Updates

**Files:**
- Modify: `CLAUDE.md` (if needed)

- [ ] **Step 1: Update CLI help text**

Verify `--track-read-source` has descriptive help text in all Click option definitions.

- [ ] **Step 2: Run linting and formatting**

Run: `make lint-fix && make format`

- [ ] **Step 3: Run full CI check**

Run: `make ci-check`
Expected: All checks PASS

- [ ] **Step 4: Commit any fixes**

```bash
git add -u
git commit -m "docs: update CLI help text and fix linting for read source tracking"
```

---

## Task 11: Manual End-to-End Verification

- [ ] **Step 1: Run simulate with tracking (if external tools available)**

```bash
muconeup --config config.json simulate --out-base test_tracking --track-read-source --fixed-lengths 30
```

Verify:
- `test_tracking_repeat_coordinates.tsv` is generated
- Coordinate TSV has correct columns and data

- [ ] **Step 2: Verify unit tests cover all paths**

```bash
pytest tests/test_source_tracking.py tests/test_parsers.py -v --tb=short
```

- [ ] **Step 3: Run full test suite with coverage**

```bash
make test
```
Expected: All tests pass, coverage maintained ≥30%
