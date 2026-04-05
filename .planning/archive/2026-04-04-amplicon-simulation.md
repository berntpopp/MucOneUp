# Amplicon Simulation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add PacBio amplicon read simulation with PCR length bias modeling, PBSIM3 template mode, and primer-based amplicon extraction.

**Architecture:** New `amplicon_pipeline.py` orchestrates an 8-stage pipeline: haplotype extraction, primer-based amplicon extraction, PCR bias coverage split, template FASTA generation, PBSIM3 template mode simulation, CCS consensus, BAM merge, and alignment. A shared `sequence_utils.py` is refactored out of `snapshot_validator.py` for primer matching.

**Tech Stack:** Python 3.10+, Click (CLI), BioPython (FASTA I/O), NumPy (stochastic PCR model), PBSIM3 (external), CCS (external), minimap2/samtools (external).

**Spec:** `.planning/specs/2026-04-04-amplicon-simulation-design.md`

---

## File Structure

### New Files
| File | Responsibility |
|------|---------------|
| `muc_one_up/read_simulator/utils/sequence_utils.py` | Shared primer binding site detection (extracted from snapshot_validator) |
| `muc_one_up/read_simulator/utils/amplicon_extractor.py` | Primer-based amplicon extraction from haplotype FASTA |
| `muc_one_up/read_simulator/utils/template_generator.py` | Replicate amplicon N times into template FASTA for PBSIM3 |
| `muc_one_up/read_simulator/pcr_bias.py` | PCR length bias model with presets and stochastic mode |
| `muc_one_up/read_simulator/amplicon_pipeline.py` | 8-stage amplicon simulation pipeline orchestrator |
| `tests/read_simulator/test_sequence_utils.py` | Tests for shared primer matching |
| `tests/read_simulator/test_amplicon_extractor.py` | Tests for amplicon extraction |
| `tests/read_simulator/test_template_generator.py` | Tests for template FASTA generation |
| `tests/read_simulator/test_pcr_bias.py` | Tests for PCR bias model |
| `tests/read_simulator/test_amplicon_pipeline.py` | Integration tests for amplicon pipeline |

### Modified Files
| File | Change |
|------|--------|
| `muc_one_up/analysis/snapshot_validator.py:156-180` | Refactor `_find_binding_sites()` to use `sequence_utils` |
| `muc_one_up/read_simulator/wrappers/pbsim3_wrapper.py` | Add `run_pbsim3_template_simulation()` |
| `muc_one_up/read_simulator/constants.py` | Add amplicon constants |
| `muc_one_up/config.py:285` | Add `"amplicon"` to simulator enum, add `amplicon_params` schema |
| `muc_one_up/read_simulation.py:64-86` | Add `"amplicon"` to `_get_simulator_map()` |
| `muc_one_up/cli/commands/reads.py` | Add `amplicon` subcommand |
| `muc_one_up/exceptions.py` | Add `AmpliconExtractionError` |
| `config.json` | Add `amplicon_params` section |

---

### Task 1: Add `AmpliconExtractionError` to exceptions

**Files:**
- Modify: `muc_one_up/exceptions.py:147`
- Test: `tests/test_exceptions.py` (verify import only)

- [ ] **Step 1: Add the exception class**

Add to `muc_one_up/exceptions.py` after `ReadSimulationError`:

```python
class AmpliconExtractionError(MucOneUpError):
    """Raised when amplicon extraction from a reference fails.

    Examples:
        - Primer binding site not found in haplotype
        - Multiple ambiguous primer binding sites
        - Extracted amplicon outside expected size range
    """
```

- [ ] **Step 2: Verify import works**

Run: `python -c "from muc_one_up.exceptions import AmpliconExtractionError; print('OK')"`
Expected: `OK`

- [ ] **Step 3: Commit**

```bash
git add muc_one_up/exceptions.py
git commit -m "feat: add AmpliconExtractionError exception class"
```

---

### Task 2: Extract shared primer matching into `sequence_utils.py`

**Files:**
- Create: `muc_one_up/read_simulator/utils/sequence_utils.py`
- Modify: `muc_one_up/analysis/snapshot_validator.py:156-180`
- Test: `tests/read_simulator/test_sequence_utils.py`

- [ ] **Step 1: Write failing tests for `find_primer_binding_sites`**

Create `tests/read_simulator/test_sequence_utils.py`:

```python
"""Tests for shared primer binding site detection."""

import pytest

from muc_one_up.read_simulator.utils.sequence_utils import (
    find_primer_binding_sites,
    reverse_complement,
)


class TestFindPrimerBindingSites:
    """Tests for find_primer_binding_sites()."""

    def test_single_match(self):
        template = "AAACCCGGGTTTAAACCC"
        primer = "CCCGGG"
        sites = find_primer_binding_sites(template, primer)
        assert sites == [3]

    def test_no_match(self):
        template = "AAATTTAAATTT"
        primer = "CCCGGG"
        sites = find_primer_binding_sites(template, primer)
        assert sites == []

    def test_multiple_matches(self):
        template = "ACGTACGTACGT"
        primer = "ACGT"
        sites = find_primer_binding_sites(template, primer)
        assert sites == [0, 4, 8]

    def test_overlapping_matches(self):
        template = "AAAA"
        primer = "AAA"
        sites = find_primer_binding_sites(template, primer)
        assert sites == [0, 1]

    def test_case_normalization(self):
        template = "aaaCCCgggTTT"
        primer = "cccGGG"
        sites = find_primer_binding_sites(template, primer)
        assert sites == [3]

    def test_reverse_complement_search(self):
        # Primer is ACGT, RC is ACGT (palindrome)
        template = "TTTACGTAAA"
        primer = "ACGT"
        sites = find_primer_binding_sites(template, primer, reverse_complement=True)
        assert sites == [3]

    def test_reverse_complement_non_palindrome(self):
        # Primer is AAACCC, RC is GGGTTT
        template = "TTTGGGTTTAAA"
        primer = "AAACCC"
        sites = find_primer_binding_sites(template, primer, reverse_complement=True)
        assert sites == [3]

    def test_empty_template(self):
        sites = find_primer_binding_sites("", "ACGT")
        assert sites == []

    def test_primer_longer_than_template(self):
        sites = find_primer_binding_sites("ACG", "ACGTACGT")
        assert sites == []


class TestReverseComplement:
    """Tests for reverse_complement()."""

    def test_simple(self):
        assert reverse_complement("ACGT") == "ACGT"  # palindrome

    def test_non_palindrome(self):
        assert reverse_complement("AAACCC") == "GGGTTT"

    def test_case_preserved(self):
        assert reverse_complement("AcGt") == "aCgT"

    def test_empty(self):
        assert reverse_complement("") == ""
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/read_simulator/test_sequence_utils.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'muc_one_up.read_simulator.utils.sequence_utils'`

- [ ] **Step 3: Implement `sequence_utils.py`**

Create `muc_one_up/read_simulator/utils/sequence_utils.py`:

```python
"""Shared sequence utilities for primer binding site detection.

Provides primer binding site detection used by both the SNaPshot
validation workflow and the amplicon extraction pipeline.
"""

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence.

    Args:
        seq: DNA sequence (ACGT characters, case preserved).

    Returns:
        Reverse complement string.
    """
    return seq.translate(_COMPLEMENT)[::-1]


def find_primer_binding_sites(
    template: str,
    primer: str,
    reverse_complement: bool = False,
) -> list[int]:
    """Find all binding sites for a primer in a template sequence.

    Both template and primer are normalized to uppercase before matching.
    Searches use exact string matching with overlapping allowed.

    Args:
        template: Template DNA sequence.
        primer: Primer sequence (5'->3').
        reverse_complement: If True, search for the reverse complement
            of the primer instead of the primer itself.

    Returns:
        List of 0-based start positions where the primer binds.
    """
    template_upper = template.upper()
    search_seq = primer.upper()

    if reverse_complement:
        search_seq = search_seq.translate(_COMPLEMENT)[::-1]

    sites: list[int] = []
    start = 0
    while True:
        pos = template_upper.find(search_seq, start)
        if pos == -1:
            break
        sites.append(pos)
        start = pos + 1

    return sites
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/read_simulator/test_sequence_utils.py -v`
Expected: All 11 tests PASS

- [ ] **Step 5: Refactor `snapshot_validator.py` to use shared utility**

In `muc_one_up/analysis/snapshot_validator.py`, replace the `_find_binding_sites` method body (lines 156-180). Change the method to delegate to the shared utility:

```python
    def _find_binding_sites(self, template: str, primer: str, reverse: bool = False) -> list[int]:
        """
        Find all binding sites for a primer in template.

        Args:
            template: Template sequence
            primer: Primer sequence (already RC'd for reverse primer in __init__)
            reverse: If True, indicates this is the reverse primer (already RC'd)

        Returns:
            List of binding site positions
        """
        from muc_one_up.read_simulator.utils.sequence_utils import (
            find_primer_binding_sites,
        )

        # Primer is already in correct orientation (RC'd in __init__ if needed)
        # so we search directly without reverse_complement flag
        return find_primer_binding_sites(template, primer, reverse_complement=False)
```

- [ ] **Step 6: Run existing snapshot validator tests to verify no regression**

Run: `pytest tests/test_snapshot_validator.py tests/test_snapshot_validator_integration.py tests/test_snapshot_complete_workflow.py -v`
Expected: All existing tests PASS (behavior unchanged)

- [ ] **Step 7: Run linting**

Run: `ruff check muc_one_up/read_simulator/utils/sequence_utils.py muc_one_up/analysis/snapshot_validator.py tests/read_simulator/test_sequence_utils.py`
Expected: No errors

- [ ] **Step 8: Commit**

```bash
git add muc_one_up/read_simulator/utils/sequence_utils.py tests/read_simulator/test_sequence_utils.py muc_one_up/analysis/snapshot_validator.py
git commit -m "refactor: extract shared primer matching into sequence_utils

Move primer binding site detection logic from PCRSimulator into a
shared utility module. snapshot_validator now delegates to the shared
function. No behavior change — existing tests pass unchanged."
```

---

### Task 3: Implement amplicon extractor

**Files:**
- Create: `muc_one_up/read_simulator/utils/amplicon_extractor.py`
- Test: `tests/read_simulator/test_amplicon_extractor.py`

- [ ] **Step 1: Write failing tests**

Create `tests/read_simulator/test_amplicon_extractor.py`:

```python
"""Tests for primer-based amplicon extraction."""

import pytest

from muc_one_up.exceptions import AmpliconExtractionError
from muc_one_up.read_simulator.utils.amplicon_extractor import (
    AmpliconExtractor,
    AmpliconResult,
)


@pytest.fixture
def muc1_primers():
    """MUC1 VNTR amplification primers (Wenzel et al. 2018)."""
    return {
        "forward": "GGAGAAAAGGAGACTTCGGCTACCCAG",
        "reverse": "GCCGTTGTGCACCAGAGTAGAAGCTGA",
    }


@pytest.fixture
def haplotype_fasta_with_primers(tmp_path, muc1_primers):
    """Create a haplotype FASTA containing primer binding sites."""
    from Bio.Seq import Seq

    fwd = muc1_primers["forward"]
    rev_rc = str(Seq(muc1_primers["reverse"]).reverse_complement())
    # Build: flanking + fwd_primer + VNTR_region + rev_primer_RC + flanking
    vntr_region = "ACGT" * 500  # 2000bp fake VNTR
    sequence = "N" * 100 + fwd + vntr_region + rev_rc + "N" * 100
    fasta = tmp_path / "haplotype.fa"
    fasta.write_text(f">hap1\n{sequence}\n")
    return fasta, len(fwd) + len(vntr_region) + len(rev_rc)


@pytest.fixture
def haplotype_fasta_no_primers(tmp_path):
    """Create a haplotype FASTA without primer binding sites."""
    fasta = tmp_path / "no_primers.fa"
    fasta.write_text(">hap1\n" + "ACGT" * 500 + "\n")
    return fasta


class TestAmpliconExtractor:
    """Tests for AmpliconExtractor."""

    def test_successful_extraction(self, haplotype_fasta_with_primers, muc1_primers, tmp_path):
        extractor = AmpliconExtractor(
            forward_primer=muc1_primers["forward"],
            reverse_primer=muc1_primers["reverse"],
        )
        fasta, expected_len = haplotype_fasta_with_primers
        output = tmp_path / "amplicon.fa"

        result = extractor.extract(str(fasta), str(output))

        assert isinstance(result, AmpliconResult)
        assert result.length == expected_len
        assert result.fasta_path == str(output)
        assert output.exists()

    def test_forward_primer_not_found(self, haplotype_fasta_no_primers, muc1_primers, tmp_path):
        extractor = AmpliconExtractor(
            forward_primer=muc1_primers["forward"],
            reverse_primer=muc1_primers["reverse"],
        )
        output = tmp_path / "amplicon.fa"

        with pytest.raises(AmpliconExtractionError, match="Forward primer.*not found"):
            extractor.extract(str(haplotype_fasta_no_primers), str(output))

    def test_reverse_primer_not_found(self, tmp_path, muc1_primers):
        """Forward primer present but reverse primer missing."""
        fwd = muc1_primers["forward"]
        fasta = tmp_path / "partial.fa"
        fasta.write_text(f">hap1\nNNNN{fwd}ACGTACGT\n")
        output = tmp_path / "amplicon.fa"

        extractor = AmpliconExtractor(
            forward_primer=muc1_primers["forward"],
            reverse_primer=muc1_primers["reverse"],
        )

        with pytest.raises(AmpliconExtractionError, match="Reverse primer.*not found"):
            extractor.extract(str(fasta), str(output))

    def test_multiple_forward_sites_raises(self, tmp_path, muc1_primers):
        """Multiple forward primer binding sites should raise."""
        from Bio.Seq import Seq

        fwd = muc1_primers["forward"]
        rev_rc = str(Seq(muc1_primers["reverse"]).reverse_complement())
        sequence = fwd + "ACGT" * 100 + fwd + "ACGT" * 100 + rev_rc
        fasta = tmp_path / "multi.fa"
        fasta.write_text(f">hap1\n{sequence}\n")
        output = tmp_path / "amplicon.fa"

        extractor = AmpliconExtractor(
            forward_primer=muc1_primers["forward"],
            reverse_primer=muc1_primers["reverse"],
        )

        with pytest.raises(AmpliconExtractionError, match="Multiple forward primer"):
            extractor.extract(str(fasta), str(output))

    def test_product_range_validation_pass(self, haplotype_fasta_with_primers, muc1_primers, tmp_path):
        extractor = AmpliconExtractor(
            forward_primer=muc1_primers["forward"],
            reverse_primer=muc1_primers["reverse"],
            expected_product_range=(100, 10000),
        )
        fasta, _ = haplotype_fasta_with_primers
        output = tmp_path / "amplicon.fa"

        result = extractor.extract(str(fasta), str(output))
        assert result.length > 0  # should succeed

    def test_product_range_validation_fail(self, haplotype_fasta_with_primers, muc1_primers, tmp_path):
        extractor = AmpliconExtractor(
            forward_primer=muc1_primers["forward"],
            reverse_primer=muc1_primers["reverse"],
            expected_product_range=(10, 50),  # too small for the test amplicon
        )
        fasta, _ = haplotype_fasta_with_primers
        output = tmp_path / "amplicon.fa"

        with pytest.raises(AmpliconExtractionError, match="outside expected"):
            extractor.extract(str(fasta), str(output))

    def test_case_insensitive_matching(self, tmp_path, muc1_primers):
        """Primers should match regardless of case in the template."""
        from Bio.Seq import Seq

        fwd = muc1_primers["forward"].lower()  # lowercase in template
        rev_rc = str(Seq(muc1_primers["reverse"]).reverse_complement()).lower()
        sequence = "nnnn" + fwd + "acgt" * 100 + rev_rc + "nnnn"
        fasta = tmp_path / "lowercase.fa"
        fasta.write_text(f">hap1\n{sequence}\n")
        output = tmp_path / "amplicon.fa"

        extractor = AmpliconExtractor(
            forward_primer=muc1_primers["forward"],
            reverse_primer=muc1_primers["reverse"],
        )

        result = extractor.extract(str(fasta), str(output))
        assert result.length > 0
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/read_simulator/test_amplicon_extractor.py -v`
Expected: FAIL — `ModuleNotFoundError`

- [ ] **Step 3: Implement `amplicon_extractor.py`**

Create `muc_one_up/read_simulator/utils/amplicon_extractor.py`:

```python
"""Primer-based amplicon extraction from haplotype FASTA files.

Extracts the PCR amplicon region from a haplotype sequence by locating
forward and reverse primer binding sites. Used by the amplicon simulation
pipeline to define the template region for PBSIM3.
"""

import logging
from pathlib import Path
from typing import NamedTuple

from ...exceptions import AmpliconExtractionError
from .sequence_utils import find_primer_binding_sites, reverse_complement


class AmpliconResult(NamedTuple):
    """Result of amplicon extraction.

    Attributes:
        sequence: Extracted amplicon DNA sequence (includes primers).
        length: Length of extracted amplicon in bp.
        forward_pos: 0-based start position of forward primer in original sequence.
        reverse_pos: 0-based start position of reverse primer RC in original sequence.
        fasta_path: Path to output FASTA file containing the amplicon.
    """

    sequence: str
    length: int
    forward_pos: int
    reverse_pos: int
    fasta_path: str


class AmpliconExtractor:
    """Extract amplicon region from a haplotype FASTA using primer sequences.

    Finds forward and reverse primer binding sites, extracts the region
    between them (inclusive of primer sequences), and writes to a new FASTA.

    Args:
        forward_primer: Forward primer sequence (5'->3').
        reverse_primer: Reverse primer sequence (5'->3', will be RC'd for search).
        expected_product_range: Optional (min, max) inclusive size range in bp.
            If provided, extracted amplicon length is validated against this range.
    """

    def __init__(
        self,
        forward_primer: str,
        reverse_primer: str,
        expected_product_range: tuple[int, int] | None = None,
    ) -> None:
        self.forward_primer = forward_primer
        self.reverse_primer = reverse_primer
        self.expected_product_range = expected_product_range
        self._reverse_primer_rc = reverse_complement(reverse_primer)

    def extract(self, haplotype_fasta: str, output_path: str) -> AmpliconResult:
        """Extract amplicon region from a haplotype FASTA file.

        Reads the first sequence from the FASTA, locates primer binding sites,
        extracts the amplicon region, and writes it to a new FASTA file.

        Args:
            haplotype_fasta: Path to input haplotype FASTA (single sequence).
            output_path: Path for output FASTA containing the extracted amplicon.

        Returns:
            AmpliconResult with extracted sequence details.

        Raises:
            AmpliconExtractionError: If primers not found, ambiguous, or
                amplicon outside expected size range.
        """
        from Bio import SeqIO

        records = list(SeqIO.parse(haplotype_fasta, "fasta"))
        if not records:
            raise AmpliconExtractionError(
                f"No sequences found in {haplotype_fasta}"
            )

        record = records[0]
        template = str(record.seq)
        seq_id = record.id

        # Find forward primer binding sites
        fwd_sites = find_primer_binding_sites(template, self.forward_primer)
        if not fwd_sites:
            raise AmpliconExtractionError(
                f"Forward primer '{self.forward_primer[:20]}...' not found in "
                f"{seq_id}. Ensure the simulated reference includes sufficient "
                f"flanking sequence to contain the primer binding site."
            )
        if len(fwd_sites) > 1:
            raise AmpliconExtractionError(
                f"Multiple forward primer binding sites found at positions "
                f"{fwd_sites} in {seq_id}. Amplicon extraction requires "
                f"unambiguous primer binding."
            )

        # Find reverse primer binding sites (search for RC on forward strand)
        rev_sites = find_primer_binding_sites(
            template, self.reverse_primer, reverse_complement=True
        )
        if not rev_sites:
            raise AmpliconExtractionError(
                f"Reverse primer '{self.reverse_primer[:20]}...' not found in "
                f"{seq_id}. Ensure the simulated reference includes sufficient "
                f"flanking sequence to contain the primer binding site."
            )
        if len(rev_sites) > 1:
            raise AmpliconExtractionError(
                f"Multiple reverse primer binding sites found at positions "
                f"{rev_sites} in {seq_id}. Amplicon extraction requires "
                f"unambiguous primer binding."
            )

        fwd_pos = fwd_sites[0]
        rev_pos = rev_sites[0]

        # Extract amplicon: from forward primer start to reverse primer end (inclusive)
        amplicon_end = rev_pos + len(self._reverse_primer_rc)
        amplicon_seq = template[fwd_pos:amplicon_end]
        amplicon_len = len(amplicon_seq)

        # Validate size range if configured
        if self.expected_product_range is not None:
            min_size, max_size = self.expected_product_range
            if amplicon_len < min_size or amplicon_len > max_size:
                raise AmpliconExtractionError(
                    f"Extracted amplicon length ({amplicon_len} bp) is outside "
                    f"expected product range [{min_size}, {max_size}] bp."
                )

        # Write amplicon to output FASTA
        output = Path(output_path)
        output.parent.mkdir(parents=True, exist_ok=True)

        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord

        amplicon_record = SeqRecord(
            Seq(amplicon_seq),
            id=f"{seq_id}_amplicon",
            description=f"amplicon fwd={fwd_pos} rev={rev_pos} len={amplicon_len}",
        )
        SeqIO.write([amplicon_record], str(output), "fasta")

        logging.info(
            "Extracted amplicon from %s: %d bp (fwd=%d, rev=%d)",
            seq_id, amplicon_len, fwd_pos, rev_pos,
        )

        return AmpliconResult(
            sequence=amplicon_seq,
            length=amplicon_len,
            forward_pos=fwd_pos,
            reverse_pos=rev_pos,
            fasta_path=str(output),
        )
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/read_simulator/test_amplicon_extractor.py -v`
Expected: All 8 tests PASS

- [ ] **Step 5: Run linting**

Run: `ruff check muc_one_up/read_simulator/utils/amplicon_extractor.py tests/read_simulator/test_amplicon_extractor.py`
Expected: No errors

- [ ] **Step 6: Commit**

```bash
git add muc_one_up/read_simulator/utils/amplicon_extractor.py tests/read_simulator/test_amplicon_extractor.py
git commit -m "feat: add primer-based amplicon extractor

Extracts PCR amplicon regions from haplotype FASTAs by locating
forward and reverse primer binding sites. Validates unambiguous
binding and optional product size range."
```

---

### Task 4: Implement template FASTA generator

**Files:**
- Create: `muc_one_up/read_simulator/utils/template_generator.py`
- Test: `tests/read_simulator/test_template_generator.py`

- [ ] **Step 1: Write failing tests**

Create `tests/read_simulator/test_template_generator.py`:

```python
"""Tests for template FASTA generator."""

import pytest
from Bio import SeqIO

from muc_one_up.read_simulator.utils.template_generator import generate_template_fasta


@pytest.fixture
def amplicon_fasta(tmp_path):
    """Create a single-sequence amplicon FASTA."""
    fasta = tmp_path / "amplicon.fa"
    fasta.write_text(">hap1_amplicon\nACGTACGTACGTACGT\n")
    return fasta


class TestGenerateTemplateFasta:
    """Tests for generate_template_fasta()."""

    def test_correct_record_count(self, amplicon_fasta, tmp_path):
        output = tmp_path / "template.fa"
        generate_template_fasta(str(amplicon_fasta), 100, str(output))
        records = list(SeqIO.parse(str(output), "fasta"))
        assert len(records) == 100

    def test_all_records_have_same_sequence(self, amplicon_fasta, tmp_path):
        output = tmp_path / "template.fa"
        generate_template_fasta(str(amplicon_fasta), 10, str(output))
        records = list(SeqIO.parse(str(output), "fasta"))
        seqs = [str(r.seq) for r in records]
        assert all(s == seqs[0] for s in seqs)

    def test_record_naming(self, amplicon_fasta, tmp_path):
        output = tmp_path / "template.fa"
        generate_template_fasta(str(amplicon_fasta), 3, str(output))
        records = list(SeqIO.parse(str(output), "fasta"))
        assert records[0].id == "amplicon_copy_001"
        assert records[1].id == "amplicon_copy_002"
        assert records[2].id == "amplicon_copy_003"

    def test_single_copy(self, amplicon_fasta, tmp_path):
        output = tmp_path / "template.fa"
        generate_template_fasta(str(amplicon_fasta), 1, str(output))
        records = list(SeqIO.parse(str(output), "fasta"))
        assert len(records) == 1

    def test_returns_output_path(self, amplicon_fasta, tmp_path):
        output = tmp_path / "template.fa"
        result = generate_template_fasta(str(amplicon_fasta), 5, str(output))
        assert result == str(output)

    def test_creates_parent_directory(self, amplicon_fasta, tmp_path):
        output = tmp_path / "subdir" / "template.fa"
        generate_template_fasta(str(amplicon_fasta), 5, str(output))
        assert output.exists()

    def test_zero_copies_raises(self, amplicon_fasta, tmp_path):
        output = tmp_path / "template.fa"
        with pytest.raises(ValueError, match="num_copies must be"):
            generate_template_fasta(str(amplicon_fasta), 0, str(output))
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/read_simulator/test_template_generator.py -v`
Expected: FAIL — `ModuleNotFoundError`

- [ ] **Step 3: Implement `template_generator.py`**

Create `muc_one_up/read_simulator/utils/template_generator.py`:

```python
"""Template FASTA generator for PBSIM3 template mode.

Replicates an amplicon sequence N times into a multi-record FASTA file.
PBSIM3 --strategy templ produces one read per FASTA record, so the
number of copies directly controls the number of simulated reads.
"""

import logging
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def generate_template_fasta(
    amplicon_fasta: str,
    num_copies: int,
    output_path: str,
) -> str:
    """Replicate an amplicon sequence into a multi-record template FASTA.

    Args:
        amplicon_fasta: Path to input FASTA with a single amplicon sequence.
        num_copies: Number of copies to generate (one read per copy).
        output_path: Path for the output template FASTA.

    Returns:
        Path to the output template FASTA file.

    Raises:
        ValueError: If num_copies < 1.
    """
    if num_copies < 1:
        raise ValueError(f"num_copies must be >= 1, got {num_copies}")

    # Read the amplicon sequence
    records = list(SeqIO.parse(amplicon_fasta, "fasta"))
    amplicon_seq = str(records[0].seq)

    # Create output directory if needed
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    # Write N copies as separate FASTA records
    template_records = [
        SeqRecord(
            Seq(amplicon_seq),
            id=f"amplicon_copy_{i:03d}",
            description="",
        )
        for i in range(1, num_copies + 1)
    ]

    SeqIO.write(template_records, str(output), "fasta")

    logging.info(
        "Generated template FASTA with %d copies of %d bp amplicon: %s",
        num_copies,
        len(amplicon_seq),
        output_path,
    )

    return str(output)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/read_simulator/test_template_generator.py -v`
Expected: All 7 tests PASS

- [ ] **Step 5: Run linting**

Run: `ruff check muc_one_up/read_simulator/utils/template_generator.py tests/read_simulator/test_template_generator.py`
Expected: No errors

- [ ] **Step 6: Commit**

```bash
git add muc_one_up/read_simulator/utils/template_generator.py tests/read_simulator/test_template_generator.py
git commit -m "feat: add template FASTA generator for PBSIM3 templ mode

Replicates an amplicon sequence N times into a multi-record FASTA.
Each record becomes one read in PBSIM3 --strategy templ mode."
```

---

### Task 5: Implement PCR bias model

**Files:**
- Create: `muc_one_up/read_simulator/pcr_bias.py`
- Modify: `muc_one_up/read_simulator/constants.py`
- Test: `tests/read_simulator/test_pcr_bias.py`

- [ ] **Step 1: Add amplicon constants**

Add to `muc_one_up/read_simulator/constants.py` at the end:

```python
# =============================================================================
# Amplicon Simulation Constants
# =============================================================================

#: Default timeout for amplicon simulation pipeline (1 hour).
DEFAULT_AMPLICON_TIMEOUT: Final[int] = 3600

#: Valid PCR bias preset names.
VALID_PCR_PRESETS: Final[set[str]] = {"default", "no_bias"}

#: Default PCR bias preset.
DEFAULT_PCR_PRESET: Final[str] = "default"
```

- [ ] **Step 2: Write failing tests for PCR bias model**

Create `tests/read_simulator/test_pcr_bias.py`:

```python
"""Tests for PCR length bias model."""

import pytest

from muc_one_up.read_simulator.pcr_bias import PCRBiasModel


class TestPCRBiasModelDeterministic:
    """Tests for deterministic PCR bias computation."""

    def test_equal_length_alleles_get_equal_coverage(self):
        model = PCRBiasModel.from_preset("default")
        n1, n2 = model.compute_coverage_split(1000, 3000, 3000)
        assert n1 == 500
        assert n2 == 500

    def test_shorter_allele_gets_more_reads(self):
        model = PCRBiasModel.from_preset("default")
        n1, n2 = model.compute_coverage_split(1000, 2000, 5000)
        assert n1 > n2
        assert n1 + n2 == 1000

    def test_total_coverage_preserved(self):
        model = PCRBiasModel.from_preset("default")
        n1, n2 = model.compute_coverage_split(1000, 1500, 6000)
        assert n1 + n2 == 1000

    def test_no_bias_preset_always_equal(self):
        model = PCRBiasModel.from_preset("no_bias")
        n1, n2 = model.compute_coverage_split(1000, 1500, 6000)
        assert n1 == 500
        assert n2 == 500

    def test_symmetric_with_swapped_alleles(self):
        model = PCRBiasModel.from_preset("default")
        n1_a, n2_a = model.compute_coverage_split(1000, 2000, 4000)
        n1_b, n2_b = model.compute_coverage_split(1000, 4000, 2000)
        assert n1_a == n2_b
        assert n2_a == n1_b

    def test_odd_total_coverage_rounds_correctly(self):
        model = PCRBiasModel.from_preset("default")
        n1, n2 = model.compute_coverage_split(101, 2000, 4000)
        assert n1 + n2 == 101


class TestPCRBiasModelPresets:
    """Tests for preset loading."""

    def test_default_preset_loads(self):
        model = PCRBiasModel.from_preset("default")
        assert model.e_max == pytest.approx(0.95)
        assert model.cycles == 25

    def test_no_bias_preset_loads(self):
        model = PCRBiasModel.from_preset("no_bias")
        assert model.e_max == pytest.approx(1.0)
        assert model.alpha == pytest.approx(0.0)

    def test_invalid_preset_raises(self):
        with pytest.raises(ValueError, match="Unknown PCR bias preset"):
            PCRBiasModel.from_preset("nonexistent")

    def test_preset_with_overrides(self):
        model = PCRBiasModel.from_preset("default", cycles=30)
        assert model.cycles == 30
        assert model.e_max == pytest.approx(0.95)  # unchanged


class TestPCRBiasModelFromParams:
    """Tests for explicit parameter construction."""

    def test_from_params(self):
        model = PCRBiasModel.from_params(
            e_max=0.90, alpha=0.0002, cycles=20,
        )
        assert model.e_max == pytest.approx(0.90)
        assert model.alpha == pytest.approx(0.0002)
        assert model.cycles == 20

    def test_stochastic_flag(self):
        model = PCRBiasModel.from_params(
            e_max=0.95, alpha=0.0001, cycles=25, stochastic=True,
        )
        assert model.stochastic is True


class TestPCRBiasModelStochastic:
    """Tests for stochastic PCR bias mode."""

    def test_stochastic_reproducible_with_seed(self):
        model = PCRBiasModel.from_params(
            e_max=0.95, alpha=0.0001, cycles=25, stochastic=True,
        )
        r1 = model.compute_coverage_split(1000, 2000, 5000, seed=42)
        r2 = model.compute_coverage_split(1000, 2000, 5000, seed=42)
        assert r1 == r2

    def test_stochastic_different_seeds_differ(self):
        model = PCRBiasModel.from_params(
            e_max=0.95, alpha=0.0001, cycles=25, stochastic=True,
        )
        r1 = model.compute_coverage_split(1000, 2000, 5000, seed=42)
        r2 = model.compute_coverage_split(1000, 2000, 5000, seed=99)
        # Very unlikely to be exactly equal with different seeds
        assert r1 != r2

    def test_stochastic_preserves_total(self):
        model = PCRBiasModel.from_params(
            e_max=0.95, alpha=0.0001, cycles=25, stochastic=True,
        )
        n1, n2 = model.compute_coverage_split(1000, 2000, 5000, seed=42)
        assert n1 + n2 == 1000

    def test_stochastic_shorter_allele_tends_higher(self):
        """Over many runs, shorter allele should average more reads."""
        model = PCRBiasModel.from_params(
            e_max=0.95, alpha=0.0003, cycles=25, stochastic=True,
        )
        ratios = []
        for seed in range(50):
            n1, n2 = model.compute_coverage_split(1000, 2000, 5000, seed=seed)
            ratios.append(n1 / (n1 + n2))
        avg_ratio = sum(ratios) / len(ratios)
        # Shorter allele (2000bp) should get > 50% on average
        assert avg_ratio > 0.5


class TestPCRBiasModelFromConfig:
    """Tests for constructing from config dict."""

    def test_from_config_with_preset(self):
        config = {"preset": "default"}
        model = PCRBiasModel.from_config(config)
        assert model.e_max == pytest.approx(0.95)

    def test_from_config_with_preset_and_overrides(self):
        config = {"preset": "default", "cycles": 30}
        model = PCRBiasModel.from_config(config)
        assert model.cycles == 30

    def test_from_config_without_preset(self):
        config = {"e_max": 0.90, "alpha": 0.0002, "cycles": 20}
        model = PCRBiasModel.from_config(config)
        assert model.e_max == pytest.approx(0.90)

    def test_from_config_empty_uses_default_preset(self):
        config = {}
        model = PCRBiasModel.from_config(config)
        assert model.e_max == pytest.approx(0.95)
```

- [ ] **Step 3: Run tests to verify they fail**

Run: `pytest tests/read_simulator/test_pcr_bias.py -v`
Expected: FAIL — `ModuleNotFoundError`

- [ ] **Step 4: Implement `pcr_bias.py`**

Create `muc_one_up/read_simulator/pcr_bias.py`:

```python
"""PCR length bias model for amplicon simulation.

Models the preferential amplification of shorter alleles in long-range PCR.
Per-cycle efficiency decays exponentially with amplicon length, and the
competitive amplification ratio compounds over PCR cycles.

Mathematical basis:
    E(L) = E_max * exp(-alpha * L)
    ratio = ((1 + E1) / (1 + E2)) ^ cycles

References:
    - Suzuki & Giovannoni 1996 — Competitive PCR model
    - Madritsch et al. 2026 (Sci Rep 16:762) — MUC1 empirical dropout data
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

from .constants import DEFAULT_PCR_PRESET, VALID_PCR_PRESETS

# Preset parameter sets
_PRESETS: dict[str, dict[str, Any]] = {
    "default": {
        "e_max": 0.95,
        "alpha": 0.00035,
        "cycles": 25,
        "denaturation_time": 10.0,
        "stochastic": False,
    },
    "no_bias": {
        "e_max": 1.0,
        "alpha": 0.0,
        "cycles": 25,
        "denaturation_time": 10.0,
        "stochastic": False,
    },
}


@dataclass(frozen=True, slots=True)
class PCRBiasModel:
    """PCR length bias model for computing per-allele coverage splits.

    Attributes:
        e_max: Maximum per-cycle efficiency for short amplicons (0-1).
        alpha: Length decay rate (bp^-1).
        cycles: Number of PCR cycles.
        denaturation_time: Denaturation step duration in seconds.
        stochastic: If True, use Galton-Watson branching process.
    """

    e_max: float
    alpha: float
    cycles: int
    denaturation_time: float = 10.0
    stochastic: bool = False

    @classmethod
    def from_preset(cls, name: str, **overrides: Any) -> PCRBiasModel:
        """Load a named preset with optional parameter overrides.

        Args:
            name: Preset name ("default" or "no_bias").
            **overrides: Individual parameters to override.

        Returns:
            Configured PCRBiasModel.

        Raises:
            ValueError: If preset name is unknown.
        """
        if name not in VALID_PCR_PRESETS:
            raise ValueError(
                f"Unknown PCR bias preset: '{name}'. "
                f"Valid presets: {', '.join(sorted(VALID_PCR_PRESETS))}"
            )
        params = {**_PRESETS[name], **overrides}
        return cls(**params)

    @classmethod
    def from_params(
        cls,
        e_max: float,
        alpha: float,
        cycles: int,
        denaturation_time: float = 10.0,
        stochastic: bool = False,
    ) -> PCRBiasModel:
        """Construct from explicit parameters."""
        return cls(
            e_max=e_max,
            alpha=alpha,
            cycles=cycles,
            denaturation_time=denaturation_time,
            stochastic=stochastic,
        )

    @classmethod
    def from_config(cls, config: dict[str, Any]) -> PCRBiasModel:
        """Construct from a config dictionary.

        If "preset" key is present, loads that preset and applies any
        other keys as overrides. Otherwise, uses explicit parameters
        or falls back to the default preset.

        Args:
            config: PCR bias config dict from amplicon_params.pcr_bias.

        Returns:
            Configured PCRBiasModel.
        """
        if not config:
            return cls.from_preset(DEFAULT_PCR_PRESET)

        preset = config.get("preset")
        if preset is not None:
            overrides = {k: v for k, v in config.items() if k != "preset"}
            return cls.from_preset(preset, **overrides)

        # No preset — need explicit params, fall back to defaults where missing
        defaults = _PRESETS[DEFAULT_PCR_PRESET]
        return cls(
            e_max=config.get("e_max", defaults["e_max"]),
            alpha=config.get("alpha", defaults["alpha"]),
            cycles=config.get("cycles", defaults["cycles"]),
            denaturation_time=config.get(
                "denaturation_time", defaults["denaturation_time"]
            ),
            stochastic=config.get("stochastic", defaults["stochastic"]),
        )

    def _efficiency(self, length: int) -> float:
        """Compute per-cycle efficiency for a given amplicon length."""
        return self.e_max * math.exp(-self.alpha * length)

    def compute_coverage_split(
        self,
        total_coverage: int,
        allele1_length: int,
        allele2_length: int,
        seed: int | None = None,
    ) -> tuple[int, int]:
        """Compute per-allele read counts from total coverage.

        Args:
            total_coverage: Total desired read count (template molecules).
            allele1_length: Length of allele 1 amplicon in bp.
            allele2_length: Length of allele 2 amplicon in bp.
            seed: Random seed (only used in stochastic mode).

        Returns:
            Tuple of (reads_allele1, reads_allele2) summing to total_coverage.
        """
        e1 = self._efficiency(allele1_length)
        e2 = self._efficiency(allele2_length)

        if self.stochastic:
            return self._stochastic_split(total_coverage, e1, e2, seed)
        return self._deterministic_split(total_coverage, e1, e2)

    def _deterministic_split(
        self, total: int, e1: float, e2: float
    ) -> tuple[int, int]:
        """Deterministic coverage split based on yield ratio."""
        yield1 = (1 + e1) ** self.cycles
        yield2 = (1 + e2) ** self.cycles
        total_yield = yield1 + yield2

        # Compute fraction for allele 1
        frac1 = yield1 / total_yield
        n1 = round(total * frac1)
        n2 = total - n1  # Ensures exact total

        return n1, n2

    def _stochastic_split(
        self, total: int, e1: float, e2: float, seed: int | None
    ) -> tuple[int, int]:
        """Stochastic coverage split using Galton-Watson branching process."""
        import numpy as np

        rng = np.random.default_rng(seed)

        # Start with equal copy numbers
        n1 = 1000
        n2 = 1000

        for _ in range(self.cycles):
            n1 += rng.binomial(n1, e1)
            n2 += rng.binomial(n2, e2)

        # Convert to read counts proportional to final molecule counts
        frac1 = n1 / (n1 + n2)
        reads1 = round(total * frac1)
        reads2 = total - reads1

        return reads1, reads2
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `pytest tests/read_simulator/test_pcr_bias.py -v`
Expected: All 18 tests PASS

- [ ] **Step 6: Run linting**

Run: `ruff check muc_one_up/read_simulator/pcr_bias.py muc_one_up/read_simulator/constants.py tests/read_simulator/test_pcr_bias.py`
Expected: No errors

- [ ] **Step 7: Commit**

```bash
git add muc_one_up/read_simulator/pcr_bias.py muc_one_up/read_simulator/constants.py tests/read_simulator/test_pcr_bias.py
git commit -m "feat: add PCR length bias model with presets

Exponential decay model for per-cycle efficiency with amplicon length.
Supports deterministic and stochastic (Galton-Watson) modes, preset
profiles (default, no_bias), and config-based construction."
```

---

### Task 6: Add PBSIM3 template mode wrapper

**Files:**
- Modify: `muc_one_up/read_simulator/wrappers/pbsim3_wrapper.py`
- Test: `tests/read_simulator/test_pbsim3_wrapper.py`

- [ ] **Step 1: Write failing tests for template mode**

Append to `tests/read_simulator/test_pbsim3_wrapper.py`:

```python
"""Tests for PBSIM3 template mode simulation."""

import pytest
from pathlib import Path
from unittest.mock import patch, MagicMock

from muc_one_up.read_simulator.wrappers.pbsim3_wrapper import (
    run_pbsim3_template_simulation,
)
from muc_one_up.exceptions import FileOperationError


class TestRunPbsim3TemplateSimulation:
    """Tests for run_pbsim3_template_simulation()."""

    def test_builds_templ_strategy_command(self, tmp_path):
        """Verify command uses --strategy templ and --template."""
        # Create required files
        template_fa = tmp_path / "template.fa"
        template_fa.write_text(">copy_001\nACGT\n>copy_002\nACGT\n")
        model_file = tmp_path / "test.model"
        model_file.write_text("model data")

        # Create fake output BAM
        output_bam = tmp_path / "out.bam"
        output_bam.write_bytes(b"BAM\x01FAKE")

        with patch(
            "muc_one_up.read_simulator.wrappers.pbsim3_wrapper.run_command"
        ) as mock_run:
            mock_run.return_value = MagicMock(returncode=0)

            run_pbsim3_template_simulation(
                pbsim3_cmd="pbsim",
                samtools_cmd="samtools",
                template_fasta=str(template_fa),
                model_type="qshmm",
                model_file=str(model_file),
                output_prefix=str(tmp_path / "out"),
            )

            cmd = mock_run.call_args[0][0]
            assert "--strategy" in cmd
            idx = cmd.index("--strategy")
            assert cmd[idx + 1] == "templ"
            assert "--template" in cmd
            assert "--genome" not in cmd
            assert "--depth" not in cmd

    def test_rejects_invalid_model_type(self, tmp_path):
        template_fa = tmp_path / "template.fa"
        template_fa.write_text(">copy\nACGT\n")
        model_file = tmp_path / "test.model"
        model_file.write_text("model")

        with pytest.raises(FileOperationError, match="Invalid pbsim3 model type"):
            run_pbsim3_template_simulation(
                pbsim3_cmd="pbsim",
                samtools_cmd="samtools",
                template_fasta=str(template_fa),
                model_type="invalid",
                model_file=str(model_file),
                output_prefix=str(tmp_path / "out"),
            )

    def test_rejects_missing_template(self, tmp_path):
        model_file = tmp_path / "test.model"
        model_file.write_text("model")

        with pytest.raises(FileOperationError, match="Template FASTA.*not found"):
            run_pbsim3_template_simulation(
                pbsim3_cmd="pbsim",
                samtools_cmd="samtools",
                template_fasta=str(tmp_path / "nonexistent.fa"),
                model_type="qshmm",
                model_file=str(model_file),
                output_prefix=str(tmp_path / "out"),
            )

    def test_includes_seed_when_provided(self, tmp_path):
        template_fa = tmp_path / "template.fa"
        template_fa.write_text(">copy\nACGT\n")
        model_file = tmp_path / "test.model"
        model_file.write_text("model")
        output_bam = tmp_path / "out.bam"
        output_bam.write_bytes(b"BAM\x01FAKE")

        with patch(
            "muc_one_up.read_simulator.wrappers.pbsim3_wrapper.run_command"
        ) as mock_run:
            mock_run.return_value = MagicMock(returncode=0)

            run_pbsim3_template_simulation(
                pbsim3_cmd="pbsim",
                samtools_cmd="samtools",
                template_fasta=str(template_fa),
                model_type="qshmm",
                model_file=str(model_file),
                output_prefix=str(tmp_path / "out"),
                seed=42,
            )

            cmd = mock_run.call_args[0][0]
            assert "--seed" in cmd
            idx = cmd.index("--seed")
            assert str(cmd[idx + 1]) == "42"
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/read_simulator/test_pbsim3_wrapper.py::TestRunPbsim3TemplateSimulation -v`
Expected: FAIL — `ImportError: cannot import name 'run_pbsim3_template_simulation'`

- [ ] **Step 3: Implement `run_pbsim3_template_simulation`**

Add to `muc_one_up/read_simulator/wrappers/pbsim3_wrapper.py` after the existing `run_pbsim3_simulation` function (before `validate_pbsim3_parameters`):

```python
def run_pbsim3_template_simulation(
    pbsim3_cmd: str,
    samtools_cmd: str,
    template_fasta: str,
    model_type: str,
    model_file: str,
    output_prefix: str,
    pass_num: int = DEFAULT_PBSIM3_PASS_NUM,
    accuracy_mean: float = DEFAULT_PBSIM3_ACCURACY_MEAN,
    seed: int | None = None,
    timeout: int = DEFAULT_PBSIM3_TIMEOUT,
) -> list[str]:
    """Simulate reads using PBSIM3 template mode.

    In template mode, PBSIM3 produces one read per FASTA record in the
    template file. Reads span the full length of each template sequence
    with realistic error profiles applied.

    Args:
        pbsim3_cmd: Path to the pbsim3 executable.
        samtools_cmd: Path to samtools executable (for SAM->BAM conversion).
        template_fasta: Path to template FASTA (N records = N reads).
        model_type: Simulation model type ("qshmm" or "errhmm").
        model_file: Path to pbsim3 model file.
        output_prefix: Prefix for output files.
        pass_num: Number of passes per molecule (>=2 for HiFi).
        accuracy_mean: Mean per-base accuracy before consensus.
        seed: Random seed for reproducibility.
        timeout: Timeout in seconds.

    Returns:
        List of paths to output BAM files.

    Raises:
        FileOperationError: If model file or template doesn't exist,
            or if model type is invalid.
    """
    # Validate model type
    if model_type not in VALID_PBSIM3_MODEL_TYPES:
        valid_types = ", ".join(sorted(VALID_PBSIM3_MODEL_TYPES))
        raise FileOperationError(
            f"Invalid pbsim3 model type: '{model_type}'. Valid options: {valid_types}"
        )

    # Validate files exist
    model_file_path = Path(model_file)
    if not model_file_path.exists():
        raise FileOperationError(f"pbsim3 model file not found: {model_file}")

    template_path = Path(template_fasta)
    if not template_path.exists():
        raise FileOperationError(f"Template FASTA file not found: {template_fasta}")

    # Build command — template mode uses different flags than WGS
    cmd_args = [
        "--strategy",
        "templ",
        "--method",
        model_type,
        "--qshmm" if model_type == "qshmm" else "--errhmm",
        model_file,
        "--template",
        template_fasta,
        "--pass-num",
        pass_num,
        "--accuracy-mean",
        accuracy_mean,
        "--prefix",
        output_prefix,
    ]

    if seed is not None:
        cmd_args.extend(["--seed", seed])

    cmd = build_tool_command(pbsim3_cmd, *cmd_args)

    logging.info("Running pbsim3 template mode simulation:")
    logging.info("  Template: %s", template_fasta)
    logging.info("  Model: %s (%s)", model_type, model_file)
    logging.info("  Pass number: %d", pass_num)
    logging.info("  Accuracy mean: %.2f", accuracy_mean)
    if seed is not None:
        logging.info("  Seed: %d", seed)

    run_command(cmd, timeout=timeout, stderr_prefix="[pbsim3-templ] ", stderr_log_level=logging.INFO)

    # Handle output files (same logic as WGS mode)
    output_prefix_path = Path(output_prefix)
    multi_bam_pattern = f"{output_prefix_path.name}_*.bam"
    multi_bam_files = sorted(
        str(p) for p in output_prefix_path.parent.glob(multi_bam_pattern)
    )

    output_bams: list[str] = []

    if multi_bam_files:
        for bam_file in multi_bam_files:
            bam_path = Path(bam_file)
            if not bam_path.exists():
                raise FileOperationError(f"pbsim3 output BAM {bam_file} not found")
            if bam_path.stat().st_size == 0:
                raise FileOperationError(f"pbsim3 produced empty BAM: {bam_file}")
            output_bams.append(bam_file)
    else:
        output_bam = f"{output_prefix}.bam"
        output_sam = f"{output_prefix}.sam"

        if Path(output_sam).exists() and not Path(output_bam).exists():
            convert_sam_to_bam(
                samtools_cmd=samtools_cmd,
                input_sam=output_sam,
                output_bam=output_bam,
                threads=4,
                timeout=timeout,
            )
            try:
                Path(output_sam).unlink()
            except Exception as e:
                logging.warning("Could not remove SAM file %s: %s", output_sam, e)

        if not Path(output_bam).exists():
            raise FileOperationError(
                f"pbsim3 template simulation failed: Expected output {output_bam} not found."
            )
        if Path(output_bam).stat().st_size == 0:
            raise FileOperationError(
                f"pbsim3 produced empty BAM: {output_bam}"
            )
        output_bams.append(output_bam)

    logging.info(
        "pbsim3 template simulation complete: %d BAM file(s)", len(output_bams)
    )
    return output_bams
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/read_simulator/test_pbsim3_wrapper.py -v`
Expected: All tests PASS (both existing WGS and new template tests)

- [ ] **Step 5: Run linting**

Run: `ruff check muc_one_up/read_simulator/wrappers/pbsim3_wrapper.py tests/read_simulator/test_pbsim3_wrapper.py`
Expected: No errors

- [ ] **Step 6: Commit**

```bash
git add muc_one_up/read_simulator/wrappers/pbsim3_wrapper.py tests/read_simulator/test_pbsim3_wrapper.py
git commit -m "feat: add PBSIM3 template mode wrapper

New run_pbsim3_template_simulation() uses --strategy templ for
amplicon-like reads. One read per FASTA record, full-length
templates with realistic error profiles."
```

---

### Task 7: Update config schema and config.json

**Files:**
- Modify: `muc_one_up/config.py:282-349`
- Modify: `config.json`

- [ ] **Step 1: Add `"amplicon"` to simulator enum and `amplicon_params` schema**

In `muc_one_up/config.py`, find line 285:

```python
"simulator": {"type": "string", "enum": ["illumina", "ont"]},
```

Change to:

```python
"simulator": {"type": "string", "enum": ["illumina", "ont", "pacbio", "amplicon"]},
```

Then find the `"additionalProperties": False` line that closes the top-level schema (line 472). Before that line, add the `amplicon_params` schema as a new property. Insert it after the `"toxic_protein_detection"` property block (after line 461):

```python
        "amplicon_params": {
            "type": "object",
            "properties": {
                "forward_primer": {"type": "string"},
                "reverse_primer": {"type": "string"},
                "primer_source": {"type": "string"},
                "expected_product_range": {
                    "type": "array",
                    "items": {"type": "number"},
                    "minItems": 2,
                    "maxItems": 2,
                },
                "pcr_bias": {
                    "type": "object",
                    "properties": {
                        "preset": {"type": "string"},
                        "e_max": {"type": "number", "minimum": 0.0, "maximum": 1.0},
                        "alpha": {"type": "number", "minimum": 0.0},
                        "cycles": {"type": "integer", "minimum": 1},
                        "denaturation_time": {"type": "number", "minimum": 0.0},
                        "stochastic": {"type": "boolean"},
                    },
                },
            },
            "required": ["forward_primer", "reverse_primer"],
            "additionalProperties": False,
        },
```

- [ ] **Step 2: Add `amplicon_params` to `config.json`**

Add the following section to `config.json` (after `pacbio_params` or at an appropriate location):

```json
"amplicon_params": {
    "forward_primer": "GGAGAAAAGGAGACTTCGGCTACCCAG",
    "reverse_primer": "GCCGTTGTGCACCAGAGTAGAAGCTGA",
    "primer_source": "Wenzel et al. 2018 (PMID: 29520014)",
    "expected_product_range": [1500, 6000],
    "pcr_bias": {
        "preset": "default"
    }
}
```

- [ ] **Step 3: Verify config loads without validation errors**

Run: `python -c "from muc_one_up.config import load_config; c = load_config('config.json'); print('amplicon_params' in c, c.get('amplicon_params', {}).get('forward_primer', '')[:10])"`
Expected: `True GGAGAAAAG`

- [ ] **Step 4: Run existing config tests to verify no regression**

Run: `pytest tests/test_config.py -v`
Expected: All existing tests PASS

- [ ] **Step 5: Commit**

```bash
git add muc_one_up/config.py config.json
git commit -m "feat: add amplicon_params to config schema and config.json

Extend simulator enum to include 'amplicon' and 'pacbio'. Add
amplicon_params schema with primer sequences, expected product
range, and PCR bias configuration."
```

---

### Task 8: Implement amplicon pipeline orchestrator

**Files:**
- Create: `muc_one_up/read_simulator/amplicon_pipeline.py`
- Test: `tests/read_simulator/test_amplicon_pipeline.py`

- [ ] **Step 1: Write failing integration tests**

Create `tests/read_simulator/test_amplicon_pipeline.py`:

```python
"""Integration tests for amplicon simulation pipeline."""

import pytest
from pathlib import Path
from unittest.mock import patch, MagicMock

from muc_one_up.read_simulator.amplicon_pipeline import simulate_amplicon_reads_pipeline


@pytest.fixture
def muc1_primers():
    return {
        "forward": "GGAGAAAAGGAGACTTCGGCTACCCAG",
        "reverse": "GCCGTTGTGCACCAGAGTAGAAGCTGA",
    }


@pytest.fixture
def diploid_fasta_with_primers(tmp_path, muc1_primers):
    """Create a diploid FASTA with primer sites in both haplotypes."""
    from Bio.Seq import Seq

    fwd = muc1_primers["forward"]
    rev_rc = str(Seq(muc1_primers["reverse"]).reverse_complement())

    # Haplotype 1: shorter VNTR
    vntr1 = "ACGT" * 400  # 1600bp
    seq1 = "N" * 50 + fwd + vntr1 + rev_rc + "N" * 50

    # Haplotype 2: longer VNTR
    vntr2 = "ACGT" * 800  # 3200bp
    seq2 = "N" * 50 + fwd + vntr2 + rev_rc + "N" * 50

    fasta = tmp_path / "diploid.fa"
    fasta.write_text(f">hap1\n{seq1}\n>hap2\n{seq2}\n")
    return fasta


@pytest.fixture
def haploid_fasta_with_primers(tmp_path, muc1_primers):
    """Create a haploid FASTA with primer sites."""
    from Bio.Seq import Seq

    fwd = muc1_primers["forward"]
    rev_rc = str(Seq(muc1_primers["reverse"]).reverse_complement())
    vntr = "ACGT" * 500
    seq = "N" * 50 + fwd + vntr + rev_rc + "N" * 50

    fasta = tmp_path / "haploid.fa"
    fasta.write_text(f">hap1\n{seq}\n")
    return fasta


@pytest.fixture
def amplicon_config(tmp_path, muc1_primers):
    """Minimal config for amplicon pipeline tests."""
    model_file = tmp_path / "test.model"
    model_file.write_text("model data")
    return {
        "tools": {
            "pbsim3": "pbsim",
            "ccs": "ccs",
            "samtools": "samtools",
            "minimap2": "minimap2",
        },
        "amplicon_params": {
            "forward_primer": muc1_primers["forward"],
            "reverse_primer": muc1_primers["reverse"],
            "pcr_bias": {"preset": "default"},
        },
        "pacbio_params": {
            "model_type": "qshmm",
            "model_file": str(model_file),
            "coverage": 100,
            "pass_num": 3,
            "min_passes": 3,
            "min_rq": 0.99,
            "threads": 4,
        },
        "read_simulation": {
            "coverage": 100,
        },
    }


class TestAmpliconPipelineDiploid:
    """Tests for diploid amplicon pipeline."""

    @patch("muc_one_up.read_simulator.amplicon_pipeline.run_pbsim3_template_simulation")
    @patch("muc_one_up.read_simulator.amplicon_pipeline.run_ccs_consensus")
    @patch("muc_one_up.read_simulator.amplicon_pipeline.merge_bam_files")
    @patch("muc_one_up.read_simulator.amplicon_pipeline.convert_bam_to_fastq")
    @patch("muc_one_up.read_simulator.amplicon_pipeline.align_reads_with_minimap2")
    def test_diploid_calls_pbsim3_twice(
        self, mock_align, mock_convert, mock_merge, mock_ccs, mock_pbsim3,
        diploid_fasta_with_primers, amplicon_config, tmp_path,
    ):
        """Diploid input should invoke PBSIM3 template mode once per allele."""
        # Setup mocks to create expected output files
        mock_pbsim3.side_effect = lambda **kw: [kw["output_prefix"] + ".bam"]
        mock_ccs.return_value = str(tmp_path / "hifi.bam")
        mock_merge.return_value = str(tmp_path / "merged.bam")
        mock_convert.return_value = str(tmp_path / "hifi.fastq")
        mock_align.return_value = str(tmp_path / "aligned.bam")

        # Create fake output files that mocks reference
        for name in ["hifi.bam", "merged.bam", "hifi.fastq", "aligned.bam"]:
            (tmp_path / name).write_bytes(b"FAKE")

        simulate_amplicon_reads_pipeline(
            config=amplicon_config,
            input_fa=str(diploid_fasta_with_primers),
            human_reference=str(tmp_path / "ref.fa"),
        )

        assert mock_pbsim3.call_count == 2

    @patch("muc_one_up.read_simulator.amplicon_pipeline.run_pbsim3_template_simulation")
    @patch("muc_one_up.read_simulator.amplicon_pipeline.run_ccs_consensus")
    @patch("muc_one_up.read_simulator.amplicon_pipeline.merge_bam_files")
    @patch("muc_one_up.read_simulator.amplicon_pipeline.convert_bam_to_fastq")
    @patch("muc_one_up.read_simulator.amplicon_pipeline.align_reads_with_minimap2")
    def test_pcr_bias_affects_template_counts(
        self, mock_align, mock_convert, mock_merge, mock_ccs, mock_pbsim3,
        diploid_fasta_with_primers, amplicon_config, tmp_path,
    ):
        """PCR bias should give shorter allele more template copies."""
        template_fastas = []

        def capture_template(**kwargs):
            template_fastas.append(kwargs["template_fasta"])
            return [kwargs["output_prefix"] + ".bam"]

        mock_pbsim3.side_effect = capture_template
        mock_ccs.return_value = str(tmp_path / "hifi.bam")
        mock_merge.return_value = str(tmp_path / "merged.bam")
        mock_convert.return_value = str(tmp_path / "hifi.fastq")
        mock_align.return_value = str(tmp_path / "aligned.bam")

        for name in ["hifi.bam", "merged.bam", "hifi.fastq", "aligned.bam"]:
            (tmp_path / name).write_bytes(b"FAKE")

        simulate_amplicon_reads_pipeline(
            config=amplicon_config,
            input_fa=str(diploid_fasta_with_primers),
            human_reference=str(tmp_path / "ref.fa"),
        )

        # Count records in each template FASTA
        from Bio import SeqIO
        counts = [
            len(list(SeqIO.parse(f, "fasta"))) for f in template_fastas
        ]
        # Shorter allele (hap1) should have more copies than longer (hap2)
        assert counts[0] > counts[1]
        assert sum(counts) == 100  # total coverage


class TestAmpliconPipelineHaploid:
    """Tests for haploid amplicon pipeline."""

    @patch("muc_one_up.read_simulator.amplicon_pipeline.run_pbsim3_template_simulation")
    @patch("muc_one_up.read_simulator.amplicon_pipeline.run_ccs_consensus")
    @patch("muc_one_up.read_simulator.amplicon_pipeline.convert_bam_to_fastq")
    @patch("muc_one_up.read_simulator.amplicon_pipeline.align_reads_with_minimap2")
    def test_haploid_calls_pbsim3_once(
        self, mock_align, mock_convert, mock_ccs, mock_pbsim3,
        haploid_fasta_with_primers, amplicon_config, tmp_path,
    ):
        """Haploid input should invoke PBSIM3 once with full coverage."""
        mock_pbsim3.side_effect = lambda **kw: [kw["output_prefix"] + ".bam"]
        mock_ccs.return_value = str(tmp_path / "hifi.bam")
        mock_convert.return_value = str(tmp_path / "hifi.fastq")
        mock_align.return_value = str(tmp_path / "aligned.bam")

        for name in ["hifi.bam", "hifi.fastq", "aligned.bam"]:
            (tmp_path / name).write_bytes(b"FAKE")

        simulate_amplicon_reads_pipeline(
            config=amplicon_config,
            input_fa=str(haploid_fasta_with_primers),
            human_reference=str(tmp_path / "ref.fa"),
        )

        assert mock_pbsim3.call_count == 1


class TestAmpliconPipelineTrackingRejection:
    """Tests for read-source tracking rejection."""

    def test_track_read_source_raises_error(
        self, haploid_fasta_with_primers, amplicon_config, tmp_path,
    ):
        """--track-read-source should raise for amplicon mode."""
        from muc_one_up.read_simulator.source_tracking import ReadSourceTracker

        tracker = MagicMock(spec=ReadSourceTracker)

        with pytest.raises(RuntimeError, match="not yet supported"):
            simulate_amplicon_reads_pipeline(
                config=amplicon_config,
                input_fa=str(haploid_fasta_with_primers),
                source_tracker=tracker,
            )
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/read_simulator/test_amplicon_pipeline.py -v`
Expected: FAIL — `ModuleNotFoundError`

- [ ] **Step 3: Implement `amplicon_pipeline.py`**

Create `muc_one_up/read_simulator/amplicon_pipeline.py`:

```python
"""Amplicon read simulation pipeline for PacBio.

Orchestrates the complete amplicon simulation workflow:
1. Haplotype extraction (diploid) or pass-through (haploid)
2. Primer-based amplicon extraction per haplotype
3. PCR bias coverage split computation
4. Template FASTA generation (N copies per allele)
5. PBSIM3 template mode simulation per allele
6. CCS consensus generation
7. BAM merge (diploid)
8. Alignment to human reference (optional)
"""

from __future__ import annotations

import logging
import tempfile
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from .output_config import OutputConfig

from ..exceptions import ExternalToolError, FileOperationError
from .constants import MINIMAP2_PRESET_PACBIO_HIFI
from .pcr_bias import PCRBiasModel
from .utils import cleanup_files, write_metadata_file
from .utils.amplicon_extractor import AmpliconExtractor
from .utils.reference_utils import extract_haplotypes, is_diploid_reference
from .utils.template_generator import generate_template_fasta
from .wrappers.ccs_wrapper import run_ccs_consensus
from .wrappers.minimap2_wrapper import align_reads_with_minimap2
from .wrappers.pbsim3_wrapper import run_pbsim3_template_simulation
from .wrappers.samtools_wrapper import convert_bam_to_fastq, merge_bam_files


def simulate_amplicon_reads_pipeline(
    config: dict[str, Any],
    input_fa: str,
    human_reference: str | None = None,
    source_tracker: Any | None = None,
    output_config: OutputConfig | None = None,
) -> str:
    """Run the complete amplicon read simulation pipeline.

    Args:
        config: Configuration dict with amplicon_params, pacbio_params, and tools.
        input_fa: Input FASTA (haploid or diploid).
        human_reference: Human reference for alignment (optional).
        source_tracker: Must be None — read source tracking is not supported
            for amplicon mode in v1.
        output_config: Optional output path configuration.

    Returns:
        Path to final output file (BAM or FASTQ).

    Raises:
        RuntimeError: If source_tracker is provided or pipeline fails.
    """
    # Reject read-source tracking for v1
    if source_tracker is not None:
        raise RuntimeError(
            "Read source tracking is not yet supported for amplicon simulation. "
            "Remove --track-read-source to proceed."
        )

    # Extract configuration
    amplicon_params = config.get("amplicon_params", {})
    pacbio_params = config.get("pacbio_params", {})
    tools = config.get("tools", {})
    rs_config = config.get("read_simulation", {})

    # Tool paths
    pbsim3_cmd = tools.get("pbsim3", "pbsim")
    ccs_cmd = tools.get("ccs", "ccs")
    samtools_cmd = tools.get("samtools", "samtools")
    minimap2_cmd = tools.get("minimap2", "minimap2")

    # PacBio params
    model_type = pacbio_params["model_type"]
    model_file = pacbio_params["model_file"]
    pass_num = pacbio_params.get("pass_num", 3)
    min_passes = pacbio_params.get("min_passes", 3)
    min_rq = pacbio_params.get("min_rq", 0.99)
    threads = pacbio_params.get("threads", 4)
    seed = pacbio_params.get("seed")
    accuracy_mean = pacbio_params.get("accuracy_mean", 0.85)

    # Amplicon params
    forward_primer = amplicon_params["forward_primer"]
    reverse_primer = amplicon_params["reverse_primer"]
    expected_range = amplicon_params.get("expected_product_range")
    expected_range_tuple = tuple(expected_range) if expected_range else None

    # Coverage
    total_coverage = rs_config.get("coverage", pacbio_params.get("coverage", 30))

    # PCR bias model
    pcr_bias_config = amplicon_params.get("pcr_bias", {})
    pcr_model = PCRBiasModel.from_config(pcr_bias_config)

    # Output paths
    input_path = Path(input_fa)
    if output_config is not None:
        output_dir = output_config.out_dir
        output_base = output_config.out_base
    else:
        output_dir = input_path.parent
        output_base = input_path.stem

    output_dir.mkdir(parents=True, exist_ok=True)

    start_time = datetime.now()
    intermediate_files: list[str] = []

    try:
        with tempfile.TemporaryDirectory(prefix="amplicon_sim_") as temp_dir:
            temp_path = Path(temp_dir)

            logging.info("=" * 80)
            logging.info("AMPLICON SIMULATION PIPELINE")
            logging.info("=" * 80)

            # ============================================================
            # STAGE 1: Haplotype detection and extraction
            # ============================================================
            diploid = is_diploid_reference(input_fa)

            if diploid:
                logging.info("STAGE 1: Extracting haplotypes from diploid reference")
                hap1_fa, hap2_fa = extract_haplotypes(
                    input_fa, temp_path, base_name="amplicon_sim"
                )
                haplotype_fastas = [hap1_fa, hap2_fa]
            else:
                logging.info("STAGE 1: Haploid reference — skipping extraction")
                haplotype_fastas = [input_fa]

            # ============================================================
            # STAGE 2: Amplicon extraction per haplotype
            # ============================================================
            logging.info("STAGE 2: Extracting amplicons via primer binding sites")

            extractor = AmpliconExtractor(
                forward_primer=forward_primer,
                reverse_primer=reverse_primer,
                expected_product_range=expected_range_tuple,
            )

            amplicon_results = []
            for i, hap_fa in enumerate(haplotype_fastas, 1):
                amp_out = str(temp_path / f"amplicon_hap{i}.fa")
                result = extractor.extract(hap_fa, amp_out)
                amplicon_results.append(result)
                logging.info(
                    "  Haplotype %d amplicon: %d bp", i, result.length
                )

            # ============================================================
            # STAGE 3: PCR bias coverage split
            # ============================================================
            if diploid:
                logging.info("STAGE 3: Computing PCR bias coverage split")
                n1, n2 = pcr_model.compute_coverage_split(
                    total_coverage,
                    amplicon_results[0].length,
                    amplicon_results[1].length,
                    seed=seed,
                )
                allele_counts = [n1, n2]
                logging.info(
                    "  Coverage split: allele1=%d, allele2=%d (total=%d)",
                    n1, n2, total_coverage,
                )
            else:
                logging.info("STAGE 3: Haploid — full coverage to single allele")
                allele_counts = [total_coverage]

            # ============================================================
            # STAGE 4: Template FASTA generation
            # ============================================================
            logging.info("STAGE 4: Generating template FASTAs")

            template_fastas = []
            for i, (amp_result, count) in enumerate(
                zip(amplicon_results, allele_counts), 1
            ):
                template_out = str(temp_path / f"template_hap{i}.fa")
                generate_template_fasta(amp_result.fasta_path, count, template_out)
                template_fastas.append(template_out)
                logging.info(
                    "  Haplotype %d: %d copies of %d bp amplicon",
                    i, count, amp_result.length,
                )

            # ============================================================
            # STAGE 5: PBSIM3 template mode simulation
            # ============================================================
            logging.info("STAGE 5: Running PBSIM3 template mode simulation")

            clr_bam_groups: list[list[str]] = []
            for i, template_fa in enumerate(template_fastas, 1):
                prefix = str(temp_path / f"clr_hap{i}")
                hap_seed = (seed + i) if seed is not None else None

                bams = run_pbsim3_template_simulation(
                    pbsim3_cmd=pbsim3_cmd,
                    samtools_cmd=samtools_cmd,
                    template_fasta=template_fa,
                    model_type=model_type,
                    model_file=model_file,
                    output_prefix=prefix,
                    pass_num=pass_num,
                    accuracy_mean=accuracy_mean,
                    seed=hap_seed,
                )
                clr_bam_groups.append(bams)
                intermediate_files.extend(bams)
                logging.info("  Haplotype %d: %d CLR BAMs", i, len(bams))

            # ============================================================
            # STAGE 6: CCS consensus generation
            # ============================================================
            logging.info("STAGE 6: Generating HiFi consensus with CCS")

            hifi_bams = []
            for i, clr_bams in enumerate(clr_bam_groups, 1):
                for j, clr_bam in enumerate(clr_bams, 1):
                    hifi_out = str(temp_path / f"hifi_hap{i}_{j:04d}.bam")
                    hap_seed = (seed + i * 100 + j) if seed is not None else None
                    hifi_bam = run_ccs_consensus(
                        ccs_cmd=ccs_cmd,
                        input_bam=clr_bam,
                        output_bam=hifi_out,
                        min_passes=min_passes,
                        min_rq=min_rq,
                        threads=threads,
                        seed=hap_seed,
                    )
                    hifi_bams.append(hifi_bam)
                    intermediate_files.append(hifi_bam)

            # ============================================================
            # STAGE 7: Merge HiFi BAMs
            # ============================================================
            hifi_merged = str(output_dir / f"{output_base}_amplicon_hifi.bam")

            if len(hifi_bams) > 1:
                logging.info("STAGE 7: Merging %d HiFi BAMs", len(hifi_bams))
                hifi_merged = merge_bam_files(
                    samtools_cmd=samtools_cmd,
                    input_bams=hifi_bams,
                    output_bam=hifi_merged,
                    threads=threads,
                )
            else:
                import shutil
                shutil.copy(hifi_bams[0], hifi_merged)

            intermediate_files.append(hifi_merged)

            # Convert to FASTQ
            hifi_fastq = str(output_dir / f"{output_base}_amplicon_hifi.fastq")
            hifi_fastq = convert_bam_to_fastq(
                samtools_cmd=samtools_cmd,
                input_bam=hifi_merged,
                output_fastq=hifi_fastq,
                threads=threads,
            )

            # ============================================================
            # STAGE 8: Alignment (optional)
            # ============================================================
            if human_reference is None:
                logging.info("Amplicon simulation complete (no alignment)")
                logging.info("Final output: %s", hifi_fastq)
                cleanup_files(intermediate_files)
                return hifi_fastq

            intermediate_files.append(hifi_fastq)

            logging.info("STAGE 8: Aligning with minimap2 (map-hifi)")
            aligned_bam = str(output_dir / f"{output_base}_amplicon_aligned.bam")
            aligned_bam = align_reads_with_minimap2(
                minimap2_cmd=minimap2_cmd,
                samtools_cmd=samtools_cmd,
                reference=human_reference,
                reads_fastq=hifi_fastq,
                output_bam=aligned_bam,
                preset=MINIMAP2_PRESET_PACBIO_HIFI,
                threads=threads,
            )

            # Cleanup and metadata
            cleanup_files(intermediate_files)

            end_time = datetime.now()
            duration = end_time - start_time
            logging.info("=" * 80)
            logging.info(
                "Amplicon simulation complete (duration: %s)",
                str(duration).split(".")[0],
            )
            logging.info("Final output: %s", aligned_bam)
            logging.info("=" * 80)

            write_metadata_file(
                output_dir=str(output_dir),
                output_base=output_base,
                config=config,
                start_time=start_time,
                end_time=end_time,
                platform="PacBio-Amplicon",
                tools_used=["pbsim3", "ccs", "minimap2", "samtools"],
            )

            return aligned_bam

    except (ExternalToolError, FileOperationError) as e:
        logging.error("Amplicon pipeline failed: %s", e)
        raise RuntimeError(
            f"Amplicon simulation pipeline failed.\nError: {e}\n\n"
            f"Troubleshooting:\n"
            f"  1. Verify tools installed: pbsim, ccs, samtools, minimap2\n"
            f"  2. Check model file: {model_file}\n"
            f"  3. Ensure primers bind in the reference"
        ) from e
    except Exception as e:
        logging.error("Unexpected error in amplicon pipeline: %s", e)
        raise RuntimeError(f"Amplicon pipeline failed: {e}") from e
    finally:
        try:
            if intermediate_files:
                cleanup_files(intermediate_files)
        except Exception as cleanup_err:
            logging.warning("Cleanup failed (non-fatal): %s", cleanup_err)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/read_simulator/test_amplicon_pipeline.py -v`
Expected: All 4 tests PASS

- [ ] **Step 5: Run linting**

Run: `ruff check muc_one_up/read_simulator/amplicon_pipeline.py tests/read_simulator/test_amplicon_pipeline.py`
Expected: No errors

- [ ] **Step 6: Commit**

```bash
git add muc_one_up/read_simulator/amplicon_pipeline.py tests/read_simulator/test_amplicon_pipeline.py
git commit -m "feat: add amplicon simulation pipeline orchestrator

8-stage pipeline: haplotype extraction, amplicon extraction, PCR bias
split, template generation, PBSIM3 templ mode, CCS, BAM merge, alignment.
Supports diploid and haploid inputs."
```

---

### Task 9: Wire up CLI command and dispatcher

**Files:**
- Modify: `muc_one_up/read_simulation.py:64-86`
- Modify: `muc_one_up/cli/commands/reads.py`

- [ ] **Step 1: Add `"amplicon"` to `_get_simulator_map()` in `read_simulation.py`**

In `muc_one_up/read_simulation.py`, add the import and map entry:

```python
def _get_simulator_map() -> dict[str, Callable[..., str]]:
    """Build simulator map with lazy imports to avoid import-time coupling."""
    from muc_one_up.read_simulator.amplicon_pipeline import (
        simulate_amplicon_reads_pipeline,
    )
    from muc_one_up.read_simulator.ont_pipeline import simulate_ont_reads_pipeline
    from muc_one_up.read_simulator.pacbio_pipeline import simulate_pacbio_hifi_reads
    from muc_one_up.read_simulator.pipeline import (
        simulate_reads_pipeline as simulate_illumina_reads,
    )

    return {
        "illumina": lambda config, input_fa, _, **kw: simulate_illumina_reads(
            config, input_fa, **kw
        ),
        "ont": lambda config, input_fa, human_reference, **kw: simulate_ont_reads_pipeline(
            config, input_fa, human_reference=human_reference, **kw
        ),
        "pacbio": lambda config, input_fa, human_reference, **kw: simulate_pacbio_hifi_reads(
            config, input_fa, human_reference=human_reference, **kw
        ),
        "amplicon": lambda config, input_fa, human_reference, **kw: simulate_amplicon_reads_pipeline(
            config, input_fa, human_reference=human_reference, **kw
        ),
    }
```

Also update the `simulator_names` dict in `simulate_reads()`:

```python
    simulator_names = {
        "illumina": "Illumina read simulation pipeline with reseq/WeSSim",
        "ont": "Oxford Nanopore (ONT) read simulation pipeline with NanoSim",
        "pacbio": "PacBio HiFi read simulation pipeline with pbsim3/CCS",
        "amplicon": "PacBio amplicon read simulation pipeline with pbsim3/CCS (template mode)",
    }
```

And update the alignment warning to include amplicon:

```python
    if simulator in ["ont", "pacbio", "amplicon"] and not human_reference:
```

- [ ] **Step 2: Add `amplicon` CLI command to `reads.py`**

Append to `muc_one_up/cli/commands/reads.py`:

```python
@reads.command()
@click.option(
    "--model-file",
    type=click.Path(exists=True, dir_okay=False),
    default=None,
    help="Path to pbsim3 model file (overrides config if provided).",
)
@click.option(
    "--model-type",
    type=click.Choice(["qshmm", "errhmm"]),
    default=None,
    help="pbsim3 model type (overrides config if provided).",
)
@click.option(
    "--pcr-preset",
    type=click.Choice(["default", "no_bias"]),
    default=None,
    help="PCR bias preset profile (default: from config or 'default').",
)
@click.option(
    "--stochastic-pcr",
    is_flag=True,
    default=False,
    help="Enable stochastic PCR bias (Galton-Watson branching process).",
)
@shared_read_options
@click.pass_context
@cli_error_handler
def amplicon(
    ctx,
    input_fastas,
    out_dir,
    out_base,
    coverage,
    model_type,
    model_file,
    pcr_preset,
    stochastic_pcr,
    seed,
    track_read_source,
):
    """Simulate PacBio amplicon reads from one or more FASTA files.

    Produces full-length amplicon reads using PBSIM3 template mode with
    realistic PCR length bias modeling. Extracts amplicon regions using
    primer binding sites from config.

    \b
    Pipeline:
      1. Extract amplicon region per haplotype (primer-based)
      2. Apply PCR length bias to determine per-allele coverage
      3. Simulate full-length reads (pbsim3 --strategy templ)
      4. Generate HiFi consensus (CCS)
      5. Align to reference (minimap2 map-hifi preset)

    \b
    Examples:
      # Basic amplicon simulation
      muconeup --config X reads amplicon sample.simulated.fa \\
        --model-file /models/QSHMM-SEQUEL.model

      # High coverage with stochastic PCR bias
      muconeup --config X reads amplicon sample.fa \\
        --model-file /models/QSHMM-SEQUEL.model \\
        --coverage 1000 --stochastic-pcr --seed 42

      # No PCR bias (equal coverage per allele)
      muconeup --config X reads amplicon sample.fa \\
        --model-file /models/QSHMM-SEQUEL.model \\
        --pcr-preset no_bias
    """
    require_config(ctx)

    # Reject --track-read-source early
    if track_read_source:
        raise click.ClickException(
            "Read source tracking is not yet supported for amplicon simulation. "
            "Remove --track-read-source to proceed."
        )

    from ...config import load_config_raw

    config = load_config_raw(str(ctx.obj["config_path"]))
    _setup_read_config(config, "amplicon", coverage, seed)

    # Ensure amplicon_params exists
    if "amplicon_params" not in config:
        raise click.ClickException(
            "Missing amplicon_params section in config. "
            "Add forward_primer and reverse_primer to config.json."
        )

    # Ensure pacbio_params exists
    if "pacbio_params" not in config:
        config["pacbio_params"] = {}

    # Override PacBio params from CLI
    if model_type is not None:
        config["pacbio_params"]["model_type"] = model_type
    if model_file is not None:
        config["pacbio_params"]["model_file"] = model_file
    if seed is not None:
        config["pacbio_params"]["seed"] = seed

    # Validate required PacBio params
    required_params = ["model_type", "model_file"]
    missing = [p for p in required_params if p not in config["pacbio_params"]]
    if missing:
        raise click.ClickException(
            f"Missing required PacBio parameters: {', '.join(missing)}. "
            f"Provide via CLI options or config.json pacbio_params section."
        )

    # Apply PCR bias overrides
    if "pcr_bias" not in config["amplicon_params"]:
        config["amplicon_params"]["pcr_bias"] = {}
    if pcr_preset is not None:
        config["amplicon_params"]["pcr_bias"]["preset"] = pcr_preset
    if stochastic_pcr:
        config["amplicon_params"]["pcr_bias"]["stochastic"] = True

    _run_batch_simulation(
        config, input_fastas, out_dir, out_base, "_amplicon", "Amplicon",
        track_read_source=False,
    )
```

- [ ] **Step 3: Verify CLI help renders**

Run: `python -m muc_one_up.cli reads amplicon --help`
Expected: Help text with all options shown

- [ ] **Step 4: Run all existing tests to verify no regression**

Run: `pytest tests/ -x -q --tb=short`
Expected: All existing tests PASS

- [ ] **Step 5: Run linting**

Run: `ruff check muc_one_up/read_simulation.py muc_one_up/cli/commands/reads.py`
Expected: No errors

- [ ] **Step 6: Commit**

```bash
git add muc_one_up/read_simulation.py muc_one_up/cli/commands/reads.py
git commit -m "feat: wire up amplicon CLI command and dispatcher

Add 'reads amplicon' subcommand with --coverage, --model-file,
--pcr-preset, --stochastic-pcr options. Register 'amplicon' in
the simulator map for pipeline dispatch."
```

---

### Task 10: Run full test suite and lint check

**Files:** None (verification only)

- [ ] **Step 1: Run all tests**

Run: `pytest tests/ -v --tb=short`
Expected: All tests PASS including new amplicon tests

- [ ] **Step 2: Run ruff lint on all changed files**

Run: `ruff check muc_one_up/ tests/`
Expected: No errors

- [ ] **Step 3: Run ruff format**

Run: `ruff format muc_one_up/ tests/`
Expected: Files formatted (or already formatted)

- [ ] **Step 4: Run mypy**

Run: `mypy muc_one_up/`
Expected: No new errors from amplicon code

- [ ] **Step 5: Commit any formatting fixes**

```bash
git add -u
git commit -m "style: fix formatting from ruff"
```

(Skip if no changes.)

- [ ] **Step 6: Final verification commit**

If all passes, no additional commit needed. The feature is complete.
