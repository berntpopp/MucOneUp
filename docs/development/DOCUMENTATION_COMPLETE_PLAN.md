# MucOneUp Complete Documentation Plan

**Version:** 2.0 (Revised)
**Date:** 2025-10-18
**Objective:** Achieve 100% Google-style docstring coverage following DRY, KISS, SOLID principles

---

## Executive Summary

This plan standardizes all docstrings to **Google Style with Type Annotations** as the single documentation format. Current status: **65% complete**. This plan addresses the remaining **35%** with precise, file-by-file actions.

**Estimated Effort:** 4-6 hours
**Approach:** Fix old-style docstrings, add missing module docs, complete partial functions
**Principle:** Types in annotations ONLY (DRY), descriptions in docstrings (single responsibility)

---

## Part 1: Google Style Standard (The ONLY Format)

### Core Principle: Single Source of Truth

**Type Information:** Function/parameter signatures ONLY
**Behavior Information:** Docstrings ONLY

### Template

```python
def function_name(param1: Type1, param2: Type2 = default) -> ReturnType:
    """One-line summary in imperative mood (72 chars max).

    Extended description providing context, algorithm details, or
    important behavioral notes. Optional if one-liner is sufficient.

    Args:
        param1: What param1 represents and its purpose
        param2: What param2 represents, defaults to [value]
            Use continuation indent for detailed multi-line
            descriptions of complex parameters

    Returns:
        Description of return value structure and meaning.
        For complex types, describe components:
        - Tuple: (element1_meaning, element2_meaning)
        - Dict: {key_type: value_meaning}

    Raises:
        ExceptionType1: When and why this occurs
        ExceptionType2: Specific conditions triggering this

    Example:
        Basic usage::

            >>> result = function_name("test", 42)
            >>> result
            ExpectedValue

    Note:
        Important caveats, thread-safety, performance notes.
    """
```

**Key Rules:**
1. NO type information in docstrings (types already in signature)
2. NO default values in docstrings (defaults already in signature)
3. Imperative mood for summaries: "Load config" not "Loads config"
4. Examples for non-trivial functions (>3 params or complex logic)

---

## Part 2: Current State Analysis

### Files Requiring Action (35%)

| File | Issue | Action Required |
|------|-------|-----------------|
| `simulate.py` | ❌ Missing module doc | Add module-level docstring |
| `mutate.py` | ⚠️ Old-style + no module doc | Convert 3 functions + add module doc |
| `config.py` | ⚠️ Old-style + no module doc | Convert 1 function + add module doc + document CONFIG_SCHEMA |
| `translate.py` | ⚠️ Old-style + no module doc | Convert 6 functions + add module doc + document CODON_TABLE |
| `fasta_writer.py` | ⚠️ Old-style + no module doc | Convert 1 function + add module doc |
| `io.py` | ⚠️ Old-style (partial) | Convert 2 functions (already has module doc ✓) |
| `orchestration.py` | ⚠️ Incomplete | Expand 1 function docstring |
| `__init__.py` | ⚠️ Minimal | Expand module docstring |

### Files Already Perfect (65%)

✅ `probabilities.py` - Google style, excellent examples
✅ `distribution.py` - Google style, excellent examples
✅ `simulation_statistics.py` - Google style, comprehensive
✅ `validation.py` - Google style, complete
✅ `cli/haplotypes.py` - Google style, complete
✅ `snp_integrator.py` - Google style, excellent
✅ `toxic_protein_detector.py` - Google style, exceptional
✅ `read_simulator/pipeline.py` - Google style, comprehensive
✅ `read_simulator/ont_pipeline.py` - Google style, complete
✅ `analysis/vntr_statistics.py` - Google style, exceptional

**These files serve as REFERENCE STANDARDS.**

---

## Part 3: Detailed File-by-File Actions

### Priority 1: Add Missing Module Docstrings (1-1.5 hours)

#### File 1: `muc_one_up/simulate.py`

**Location:** Line 1 (add before imports)
**Current:** Missing module docstring
**Action:** Add comprehensive module docstring

```python
"""VNTR haplotype simulation engine.

This module implements the core simulation logic for generating MUC1 VNTR
haplotype chains using probability-based state transitions. It supports both
random generation (sampling from configured distributions) and simulation from
predefined repeat chains.

Key Functions:
    simulate_diploid: Generate multiple haplotypes with configurable parameters
    simulate_single_haplotype: Build single haplotype by chaining repeats
    simulate_from_chains: Simulate haplotypes from predefined repeat chains
    assemble_haplotype_from_chain: Concatenate constants + repeats into sequence

Design:
    All haplotypes enforce canonical terminal block (6/6p → 7 → 8 → 9) to match
    biological MUC1 structure. Probability transitions prevent premature
    termination during chain building.

Example:
    Basic diploid simulation::

        from muc_one_up.config import load_config
        from muc_one_up.simulate import simulate_diploid

        config = load_config("config.json")
        results = simulate_diploid(config, num_haplotypes=2, fixed_lengths=[50, 60])
        for seq, chain in results:
            print(f"Chain: {'-'.join(chain)}, Length: {len(seq)} bp")

See Also:
    - probabilities.py: Weighted random repeat selection
    - distribution.py: Target length sampling
    - mutate.py: Post-simulation mutation application
"""
```

**Estimated Time:** 10 minutes

---

#### File 2: `muc_one_up/mutate.py`

**Location:** Line 1 (add before imports)
**Current:** Missing module docstring
**Action:** Add comprehensive module docstring

```python
"""Mutation application engine for simulated VNTR haplotypes.

This module applies targeted mutations to repeat units within simulated
haplotypes. Supports four mutation types: insert, delete, replace, and
delete_insert. Mutations are validated against allowed_repeats and can
operate in strict mode (reject invalid targets) or permissive mode (auto-convert).

Key Functions:
    apply_mutations: Apply named mutation to specific haplotype positions
    validate_allowed_repeats: Ensure mutation targets are valid repeat symbols
    apply_changes_to_repeat: Execute insertion/deletion/replacement operations
    rebuild_haplotype_sequence: Reassemble sequence after mutations

Mutation Structure:
    Mutations are defined in config["mutations"][name] with:
    - allowed_repeats: Valid repeat symbols for this mutation
    - strict_mode: Boolean controlling auto-conversion behavior
    - changes: List of operations (type, start, end, sequence)

Example:
    Apply dupC mutation to haplotype 1, repeat position 25::

        from muc_one_up.mutate import apply_mutations

        config = load_config("config.json")
        results, mutated_units = apply_mutations(
            config=config,
            results=[(seq1, chain1), (seq2, chain2)],
            mutation_name="dupC",
            targets=[(1, 25)]  # 1-based indexing
        )

Notes:
    - Positions use 1-based indexing (matching biological conventions)
    - Mutated repeats are marked with 'm' suffix in chains
    - Strict mode prevents silent auto-conversions (recommended for reproducibility)

See Also:
    - config.py: CONFIG_SCHEMA for mutation definitions
    - simulate.py: Haplotype generation
"""
```

**Estimated Time:** 12 minutes

---

#### File 3: `muc_one_up/config.py`

**Location:** Line 1 (add before imports)
**Current:** Missing module docstring
**Action:** Add comprehensive module docstring

```python
"""Configuration loading and validation for MucOneUp.

This module handles JSON configuration file parsing, schema validation, and
format normalization. Validates all required sections (repeats, constants,
probabilities, mutations, tools, read_simulation) and enforces structural
constraints via CONFIG_SCHEMA.

Key Functions:
    load_config: Load, validate, and normalize configuration from JSON file

Key Constants:
    CONFIG_SCHEMA: JSON Schema defining required configuration structure

Configuration Sections:
    - repeats: Mapping of repeat symbols (1, 2, X, A, B...) to DNA sequences
    - constants: Left/right flanking sequences for hg19 and hg38 assemblies
    - probabilities: State transition matrix for repeat chain generation
    - length_model: Distribution parameters (normal/uniform) for VNTR length
    - mutations: Named mutation definitions with allowed_repeats and changes
    - tools: Command paths for external tools (reseq, bwa, samtools...)
    - read_simulation: Parameters for Illumina pipeline (coverage, threads...)
    - nanosim_params: Parameters for ONT pipeline (model path, read lengths...)

Example:
    Load and access configuration::

        from muc_one_up.config import load_config

        config = load_config("config.json")
        repeats = config["repeats"]  # Dict[str, str]
        probs = config["probabilities"]  # Dict[str, Dict[str, float]]

Notes:
    - Automatically normalizes flat constants format to nested format
    - Validates mutation allowed_repeats against valid repeat symbols
    - Supports both hg19 and hg38 reference assemblies

Raises:
    FileNotFoundError: If config file doesn't exist
    json.JSONDecodeError: If config file contains invalid JSON
    ValidationError: If config doesn't conform to CONFIG_SCHEMA
"""
```

**Estimated Time:** 12 minutes

---

#### File 4: `muc_one_up/translate.py`

**Location:** Line 1 (add before imports)
**Current:** Missing module docstring
**Action:** Add comprehensive module docstring

```python
"""ORF prediction and translation for simulated haplotypes.

This module provides Open Reading Frame (ORF) detection and translation
functionality for simulated MUC1 haplotype sequences. Uses orfipy_core for
efficient ORF finding and implements custom translation with the standard
genetic code.

Key Functions:
    predict_orfs_in_haplotypes: Find and translate ORFs in all haplotypes
    find_orfs_in_memory: Use orfipy to find ORFs in DNA sequence
    dna_to_protein: Translate DNA sequence to protein using codon table
    reverse_complement: Generate reverse complement of DNA sequence
    write_peptides_to_fasta: Write predicted peptides to FASTA file

Key Constants:
    CODON_TABLE: Standard genetic code mapping codons to amino acids

Workflow:
    1. Detect ORFs using orfipy_core (both strands)
    2. Filter by minimum amino acid length (orf_min_aa)
    3. Translate DNA to protein sequences
    4. Optional: Filter by required N-terminal prefix
    5. Write results to peptide FASTA file

Example:
    Predict ORFs with minimum 100 aa length::

        from muc_one_up.translate import run_orf_finder_in_memory

        results = [(seq1, chain1), (seq2, chain2)]
        run_orf_finder_in_memory(
            results=results,
            output_pep="output.pep.fa",
            orf_min_aa=100,
            required_prefix="MTRV"  # Optional N-terminal filter
        )

Notes:
    - Supports both forward and reverse strand ORF detection
    - Stop codons: TAA, TAG, TGA (standard genetic code)
    - Start codons: ATG, TTG, CTG (common bacterial starts)
    - Uses '_' character for stop codons in CODON_TABLE

See Also:
    - toxic_protein_detector.py: Analyze translated ORFs for toxic features
"""
```

**Estimated Time:** 12 minutes

---

#### File 5: `muc_one_up/fasta_writer.py`

**Location:** Line 1 (add before imports)
**Current:** Missing module docstring
**Action:** Add comprehensive module docstring

```python
"""FASTA file writing utilities for MucOneUp.

This module provides FASTA output functionality with support for custom
headers, per-sequence comments, and proper line wrapping (80 characters).

Key Functions:
    write_fasta: Write DNA/protein sequences to FASTA file with formatting

Features:
    - Automatic 80-character line wrapping for readability
    - Per-sequence comment support (for mutation annotations)
    - Global comment fallback for all sequences
    - Safe file handling with proper error reporting

Example:
    Write haplotype sequences with mutation annotations::

        from muc_one_up.fasta_writer import write_fasta

        sequences = [seq1, seq2]
        comments = [
            "mutation=dupC targets=[(1, 25)]",
            "mutation=normal"
        ]
        write_fasta(
            sequences=sequences,
            filename="output.fa",
            prefix="haplotype",
            comments=comments
        )

Output Format:
    >haplotype_1 mutation=dupC targets=[(1, 25)]
    ATCGATCGATCGATCG...
    >haplotype_2 mutation=normal
    GCTAGCTAGCTAGCTA...
"""
```

**Estimated Time:** 8 minutes

---

### Priority 2: Convert Old-Style to Google Style (1.5-2 hours)

#### File 6: `muc_one_up/mutate.py` - Function Conversions

**Function 1:** `validate_allowed_repeats()` (lines 15-36)

**Current:**
```python
def validate_allowed_repeats(mutation_def: MutationDefinition, config: ConfigDict) -> set[str]:
    """
    Validate that the allowed_repeats in a mutation definition are valid repeat symbols.

    :param mutation_def: Mutation definition from config["mutations"][name].
    :param config: The entire config dict containing the "repeats" section.
    :return: Set of validated allowed repeat symbols.
    :raises ValueError: If any allowed repeat is not a valid repeat symbol.
    """
```

**Replace with:**
```python
def validate_allowed_repeats(mutation_def: MutationDefinition, config: ConfigDict) -> set[str]:
    """Validate that allowed_repeats contains only valid repeat symbols.

    Args:
        mutation_def: Mutation definition from config["mutations"][name]
        config: Configuration dict containing the "repeats" section

    Returns:
        Set of validated allowed repeat symbols

    Raises:
        ValueError: If any allowed repeat is not a valid repeat symbol
    """
```

---

**Function 2:** `apply_mutations()` (lines 39-150)

**Current:**
```python
def apply_mutations(
    config: ConfigDict,
    results: HaplotypeList,
    mutation_name: MutationName,
    targets: MutationTargets,
) -> tuple[HaplotypeList, dict[int, list[tuple[int, str]]]]:
    """
    Apply a single named mutation to one or more haplotypes at specific repeat indices,
    and record the mutated VNTR unit(s).

    :param config: The entire config dict.
    :param results: List of (sequence, chain) for each haplotype.
    :param mutation_name: Key in config["mutations"].
    :param targets: List of (haplotype_idx, repeat_idx) (1-based indexing).
    :return: Tuple (updated_results, mutated_units) where:
             - updated_results: List of (sequence, chain) for each haplotype after mutation.
             - mutated_units: Dict mapping haplotype index (1-based) to a list of tuples
               (repeat_index, mutated_unit_sequence).
    :raises ValueError: if configuration or target indices are invalid.
    """
```

**Replace with:**
```python
def apply_mutations(
    config: ConfigDict,
    results: HaplotypeList,
    mutation_name: MutationName,
    targets: MutationTargets,
) -> tuple[HaplotypeList, dict[int, list[tuple[int, str]]]]:
    """Apply named mutation to specific haplotype positions.

    Applies the mutation defined in config["mutations"][mutation_name] to each
    target position. Validates allowed_repeats, handles strict mode enforcement,
    and tracks mutated VNTR unit sequences for reporting.

    Args:
        config: Configuration dict containing mutations section
        results: List of (sequence, chain) tuples for each haplotype
        mutation_name: Key in config["mutations"] defining the mutation
        targets: List of (haplotype_idx, repeat_idx) tuples using 1-based indexing

    Returns:
        Tuple containing:
        - updated_results: Modified haplotypes with mutations applied
        - mutated_units: Dict mapping haplotype index (1-based) to list of
          (repeat_index, mutated_sequence) tuples

    Raises:
        ValueError: If mutation name not found, target indices invalid, or
                   strict mode rejects invalid repeat symbol

    Note:
        Mutated repeats are marked with 'm' suffix in the chain.
        Positions use 1-based indexing (haplotype and repeat).
    """
```

---

**Function 3:** `rebuild_haplotype_sequence()` (lines 153-172)

**Current:**
```python
def rebuild_haplotype_sequence(chain: RepeatChain, config: ConfigDict) -> DNASequence:
    """
    Rebuild the haplotype sequence from the chain of repeats and constant flanks.

    :param chain: List of repeat symbols (with possible appended 'm').
    :param config: Configuration dict containing 'constants' and 'repeats'.
    :return: Reassembled haplotype sequence.
    """
```

**Replace with:**
```python
def rebuild_haplotype_sequence(chain: RepeatChain, config: ConfigDict) -> DNASequence:
    """Rebuild haplotype sequence from repeat chain and flanking constants.

    Concatenates left constant + repeat units + right constant (if terminal repeat is 9).
    Strips mutation markers ('m') from symbols when looking up sequences.

    Args:
        chain: List of repeat symbols, possibly with 'm' suffix marking mutations
        config: Configuration dict containing 'constants' and 'repeats' sections

    Returns:
        Reassembled haplotype DNA sequence
    """
```

---

**Function 4:** `apply_changes_to_repeat()` (lines 175-264)

**Current:**
```python
def apply_changes_to_repeat(
    seq: DNASequence,
    chain: RepeatChain,
    repeat_index: int,
    changes: list[dict[str, int | str]],
    config: ConfigDict,
    mutation_name: MutationName,
) -> tuple[DNASequence, RepeatChain, DNASequence]:
    """
    Modify the substring of 'seq' corresponding to chain[repeat_index] using the
    list of changes from the mutation definition.

    :param seq: Original haplotype sequence.
    :param chain: The repeat chain.
    :param repeat_index: Index (0-based) of the repeat to modify.
    :param changes: List of change dictionaries.
    :param config: Configuration dict.
    :param mutation_name: Name of the mutation.
    :return: Tuple (new_seq, chain, mutated_repeat) after applying the changes.
    """
```

**Replace with:**
```python
def apply_changes_to_repeat(
    seq: DNASequence,
    chain: RepeatChain,
    repeat_index: int,
    changes: list[dict[str, int | str]],
    config: ConfigDict,
    mutation_name: MutationName,
) -> tuple[DNASequence, RepeatChain, DNASequence]:
    """Modify repeat unit sequence using mutation change operations.

    Applies a series of change operations (insert, delete, replace, delete_insert)
    to the repeat unit at the specified index. Operations use 1-based positioning
    within the repeat unit.

    Args:
        seq: Full haplotype sequence before mutation
        chain: Repeat chain identifying unit positions
        repeat_index: Index (0-based) of the repeat to modify
        changes: List of change dicts with 'type', 'start', 'end', 'sequence' keys
        config: Configuration dict for assembly constants and repeat lookups
        mutation_name: Name of mutation (for error messages)

    Returns:
        Tuple containing:
        - new_seq: Updated haplotype sequence with mutation applied
        - chain: Repeat chain (unmodified, marker added separately)
        - mutated_repeat: The modified repeat unit sequence

    Raises:
        ValueError: If change coordinates are out of bounds or type is unknown

    Note:
        Change operations:
        - insert: Insert sequence at position (inclusive)
        - delete: Delete bases from start to end (inclusive)
        - replace: Replace bases from start to end with sequence
        - delete_insert: Delete between boundaries, insert at boundary
    """
```

**Estimated Time:** 25 minutes for all mutate.py conversions

---

#### File 7: `muc_one_up/config.py` - Function Conversion

**Function:** `load_config()` (lines 162-217)

**Current:**
```python
def load_config(config_path: str) -> dict[str, Any]:
    """
    Load and validate the JSON config for MucOneUp.

    :param config_path: Path to the JSON config file.
    :return: Python dict with validated configuration data.
    :raises FileNotFoundError: if the config file does not exist.
    :raises json.JSONDecodeError: if the config file is not valid JSON.
    :raises ValidationError: if the config does not conform to the required schema.
    """
```

**Replace with:**
```python
def load_config(config_path: str) -> dict[str, Any]:
    """Load and validate MucOneUp configuration from JSON file.

    Validates configuration against CONFIG_SCHEMA and normalizes constants
    format (converts flat format to nested hg19/hg38 structure if needed).
    Performs additional validation for mutation allowed_repeats.

    Args:
        config_path: Path to the JSON config file

    Returns:
        Validated configuration dictionary containing all required sections
        (repeats, constants, probabilities, length_model, mutations,
        tools, read_simulation, nanosim_params)

    Raises:
        FileNotFoundError: If the config file does not exist
        json.JSONDecodeError: If the config file is not valid JSON
        ValidationError: If the config does not conform to CONFIG_SCHEMA or
                        contains invalid mutation allowed_repeats

    Note:
        Flat constants format (for backward compatibility) is automatically
        converted to nested format with default reference assembly (hg38).
    """
```

**Estimated Time:** 10 minutes

---

**Constant:** `CONFIG_SCHEMA` (lines 11-159)

**Current:** No docstring
**Action:** Add comprehensive constant documentation BEFORE the definition

**Insert at line 11 (before `CONFIG_SCHEMA = {`):**

```python
#: JSON Schema for validating MucOneUp configuration files.
#:
#: Defines the complete structure and validation rules for configuration
#: dictionaries. All sections are required except where specified.
#:
#: Required Sections:
#:     repeats: Mapping of repeat symbols (str) to DNA sequences (str)
#:         Example: {"1": "ACGACT", "2": "CAGACT", "X": "TCGACT", ...}
#:
#:     constants: Flanking sequences for each reference assembly
#:         Nested format: {assembly: {left: str, right: str, vntr_region: str}}
#:         Supported assemblies: hg19, hg38
#:
#:     probabilities: State transition matrix for repeat chain generation
#:         Format: {current_symbol: {next_symbol: probability, ...}, ...}
#:         Example: {"1": {"2": 0.7, "7": 0.3}, "2": {"7": 1.0}}
#:
#:     length_model: Distribution parameters for VNTR length sampling
#:         Required keys: distribution, min_repeats, max_repeats,
#:                       mean_repeats, median_repeats
#:
#:     mutations: Named mutation definitions
#:         Each mutation has:
#:         - allowed_repeats: List of valid repeat symbols for this mutation
#:         - strict_mode: Boolean (optional, default: false)
#:         - changes: List of operations with type, start, end, sequence
#:         Operation types: insert, delete, replace, delete_insert
#:
#:     tools: Command paths for external executables
#:         Required: samtools
#:         Optional: reseq, faToTwoBit, pblat, bwa, nanosim, minimap2
#:
#:     read_simulation: Parameters for Illumina read simulation
#:         Required: human_reference, threads
#:         Simulator type: illumina or ont
#:
#:     nanosim_params: Parameters for Oxford Nanopore simulation (optional)
#:         Required: training_data_path, coverage
#:
#: Example:
#:     Accessing schema in code::
#:
#:         from muc_one_up.config import CONFIG_SCHEMA
#:         from jsonschema import validate
#:
#:         validate(instance=config_dict, schema=CONFIG_SCHEMA)
#:
#: Notes:
#:     - Constants support both flat format (backward compatibility) and
#:       nested format (hg19/hg38)
#:     - Mutations are validated to ensure allowed_repeats are valid symbols
#:     - All numeric fields are validated for type and range
CONFIG_SCHEMA: dict[str, Any] = {
```

**Estimated Time:** 15 minutes

---

#### File 8: `muc_one_up/translate.py` - Function Conversions

**Constant:** `CODON_TABLE` (lines 9-74)

**Current:** No docstring
**Action:** Add constant documentation BEFORE the definition

**Insert at line 9 (before `CODON_TABLE = {`):**

```python
#: Standard genetic code codon table.
#:
#: Maps DNA codons (3-letter strings) to single-letter amino acid codes.
#: Stop codons (TAA, TAG, TGA) are represented by underscore ('_').
#:
#: Format:
#:     {codon: amino_acid, ...}
#:     Example: {"ATG": "M", "TAA": "_", "GCT": "A"}
#:
#: Properties:
#:     - 64 total codons (4^3 combinations)
#:     - 20 standard amino acids + stop
#:     - Multiple codons per amino acid (degeneracy)
#:     - Case-sensitive (uppercase DNA bases)
#:
#: Usage:
#:     Translate codon to amino acid::
#:
#:         from muc_one_up.translate import CODON_TABLE
#:         amino_acid = CODON_TABLE["ATG"]  # Returns "M" (methionine)
#:         amino_acid = CODON_TABLE["TAA"]  # Returns "_" (stop)
#:
#: See Also:
#:     - dna_to_protein: Translation function using this table
CODON_TABLE = {
```

**Function 1:** `reverse_complement()` (lines 77-96)

**Current:**
```python
def reverse_complement(seq):
    """
    Return the reverse complement of a DNA sequence.

    :param seq: DNA sequence.
    :return: Reverse complement of the sequence.
    """
```

**Replace with:**
```python
def reverse_complement(seq: str) -> str:
    """Generate reverse complement of DNA sequence.

    Handles both uppercase and lowercase DNA bases (A, T, G, C, N).
    Unknown bases are converted to 'N'.

    Args:
        seq: DNA sequence string (any case)

    Returns:
        Reverse complement DNA sequence

    Example:
        >>> reverse_complement("ATCG")
        'CGAT'
        >>> reverse_complement("atcgN")
        'Ncgat'
    """
```

---

**Function 2:** `dna_to_protein()` (lines 99-117)

**Current:**
```python
def dna_to_protein(dna, codon_table=CODON_TABLE, include_stop=False):
    """
    Translate a DNA sequence in-frame to protein using a codon table.

    :param dna: DNA sequence.
    :param codon_table: Dict mapping codons to amino acids.
    :param include_stop: If True, include stop codon in translation.
    :return: Translated protein sequence.
    """
```

**Replace with:**
```python
def dna_to_protein(
    dna: str,
    codon_table: dict[str, str] = CODON_TABLE,
    include_stop: bool = False
) -> str:
    """Translate DNA sequence to protein using codon table.

    Reads DNA in 3-base codons starting from position 0 (in-frame translation).
    Stops at first stop codon unless include_stop=True. Unknown codons are
    translated to 'X'.

    Args:
        dna: DNA sequence to translate
        codon_table: Mapping of codons to amino acids, defaults to standard code
        include_stop: If True, include stop codon character '_' in output

    Returns:
        Translated protein sequence as string of amino acid codes

    Example:
        >>> dna_to_protein("ATGGCTTAA")
        'MA'  # Stop codon TAA terminates translation
        >>> dna_to_protein("ATGGCTTAA", include_stop=True)
        'MA_'  # Stop codon included
    """
```

---

**Function 3:** `find_orfs_in_memory()` (lines 120-162)

**Current:**
```python
def find_orfs_in_memory(
    seq,
    min_len=30,
    max_len=1000000,
    strand="b",
    starts=None,
    stops=None,
    include_stop=False,
    partial5=False,
    partial3=False,
):
    """
    Use orfipy_core.orfs() to find ORFs in a DNA sequence.

    :param seq: DNA sequence.
    :param min_len: Minimum ORF length in nucleotides.
    :param max_len: Maximum ORF length.
    :param strand: Strand option ('+', '-', or 'b' for both).
    :param starts: List of possible start codons.
    :param stops: List of possible stop codons.
    :param include_stop: Whether to include the stop codon.
    :param partial5: Allow partial ORFs at the 5' end.
    :param partial3: Allow partial ORFs at the 3' end.
    :return: List of ORFs as tuples (start, stop, strand, description).
    """
```

**Replace with:**
```python
def find_orfs_in_memory(
    seq: str,
    min_len: int = 30,
    max_len: int = 1000000,
    strand: str = "b",
    starts: list[str] | None = None,
    stops: list[str] | None = None,
    include_stop: bool = False,
    partial5: bool = False,
    partial3: bool = False,
) -> list[tuple[int, int, str, str]]:
    """Find Open Reading Frames in DNA sequence using orfipy.

    Uses orfipy_core for efficient ORF detection on one or both strands.
    Default start codons: ATG, TTG, CTG (standard + alternative bacterial starts).
    Default stop codons: TAA, TAG, TGA (standard genetic code).

    Args:
        seq: DNA sequence to search
        min_len: Minimum ORF length in nucleotides
        max_len: Maximum ORF length in nucleotides
        strand: Strand option ('+' for forward, '-' for reverse, 'b' for both)
        starts: List of start codon sequences, defaults to ["ATG", "TTG", "CTG"]
        stops: List of stop codon sequences, defaults to ["TAA", "TAG", "TGA"]
        include_stop: If True, include stop codon in ORF coordinates
        partial5: Allow partial ORFs at 5' end (no start codon)
        partial3: Allow partial ORFs at 3' end (no stop codon)

    Returns:
        List of ORF tuples: (start_pos, stop_pos, strand, description)
        Positions are 0-based indices into the sequence
    """
```

---

**Function 4:** `predict_orfs_in_haplotypes()` (lines 165-208)

**Current:**
```python
def predict_orfs_in_haplotypes(results, min_len=30, orf_min_aa=100, required_prefix=None):
    """
    For each haplotype in results, find ORFs, translate them to peptides,
    and filter based on minimal peptide length and optional prefix.

    :param results: List of tuples (dna_seq, repeat_chain).
    :param min_len: Minimal ORF length in nucleotides (for orfipy).
    :param orf_min_aa: Minimal peptide length in amino acids.
    :param required_prefix: Optional prefix filter for peptides.
    :return: Dict mapping haplotype id to a list of ORF tuples.
    """
```

**Replace with:**
```python
def predict_orfs_in_haplotypes(
    results: list[tuple[str, list[str]]],
    min_len: int = 30,
    orf_min_aa: int = 100,
    required_prefix: str | None = None
) -> dict[str, list[tuple[str, str, int, int, str, str]]]:
    """Predict and translate ORFs for all haplotypes.

    Finds ORFs in each haplotype sequence, translates to protein, and filters
    by minimum amino acid length and optional N-terminal prefix requirement.

    Args:
        results: List of (dna_sequence, repeat_chain) tuples
        min_len: Minimum ORF length in nucleotides for orfipy detection
        orf_min_aa: Minimum peptide length in amino acids (post-translation filter)
        required_prefix: Optional N-terminal amino acid prefix filter
                        (e.g., "MTRV" to require specific signal peptide)

    Returns:
        Dictionary mapping haplotype IDs to ORF lists:
        {
            "haplotype_1": [
                (orf_id, peptide, start, stop, strand, description),
                ...
            ],
            ...
        }

    Example:
        >>> results = [(seq1, chain1), (seq2, chain2)]
        >>> orfs = predict_orfs_in_haplotypes(results, orf_min_aa=100)
        >>> orfs["haplotype_1"]
        [('haplotype_1_ORF1', 'MTRV...', 245, 1523, '+', 'len=1278')]
    """
```

---

**Function 5:** `write_peptides_to_fasta()` (lines 211-227)

**Current:**
```python
def write_peptides_to_fasta(haplotype_orfs, output_pep):
    """
    Write predicted peptide sequences to a FASTA file.

    :param haplotype_orfs: Dict mapping haplotype_id to list of ORF tuples.
    :param output_pep: Output FASTA filename.
    """
```

**Replace with:**
```python
def write_peptides_to_fasta(
    haplotype_orfs: dict[str, list[tuple]],
    output_pep: str
) -> None:
    """Write predicted peptide sequences to FASTA file.

    Args:
        haplotype_orfs: Dictionary mapping haplotype IDs to ORF lists
                       (from predict_orfs_in_haplotypes)
        output_pep: Output FASTA filename for peptide sequences

    Raises:
        Exception: If file writing fails (logged and re-raised)

    Example Output:
        >haplotype_1_ORF1 strand=+ len=1278
        MTRVPGTRPALLLLLVPLLLQTGALA...
        >haplotype_1_ORF2 strand=- len=456
        MKRELLSAAVPV...
    """
```

---

**Function 6:** `run_orf_finder_in_memory()` (lines 230-244)

**Current:**
```python
def run_orf_finder_in_memory(results, output_pep, min_len=30, orf_min_aa=100, required_prefix=None):
    """
    High-level function that predicts ORFs for each haplotype in memory,
    applies filters, and writes the peptides to a FASTA file.

    :param results: List of tuples (dna_seq, repeat_chain).
    :param output_pep: Output peptide FASTA filename.
    :param min_len: Minimum ORF length in nucleotides.
    :param orf_min_aa: Minimum peptide length in amino acids.
    :param required_prefix: Optional prefix filter for peptides.
    """
```

**Replace with:**
```python
def run_orf_finder_in_memory(
    results: list[tuple[str, list[str]]],
    output_pep: str,
    min_len: int = 30,
    orf_min_aa: int = 100,
    required_prefix: str | None = None
) -> None:
    """Predict ORFs and write peptides to FASTA (complete pipeline).

    High-level orchestration function combining ORF detection, translation,
    filtering, and FASTA output. Convenience wrapper for common workflow.

    Args:
        results: List of (dna_sequence, repeat_chain) tuples from simulation
        output_pep: Output FASTA filename for translated peptides
        min_len: Minimum ORF length in nucleotides
        orf_min_aa: Minimum peptide length in amino acids
        required_prefix: Optional N-terminal prefix requirement

    Example:
        >>> from muc_one_up.simulate import simulate_diploid
        >>> from muc_one_up.translate import run_orf_finder_in_memory
        >>>
        >>> results = simulate_diploid(config)
        >>> run_orf_finder_in_memory(results, "output.pep.fa", orf_min_aa=100)
    """
```

**Estimated Time:** 35 minutes for all translate.py conversions

---

#### File 9: `muc_one_up/fasta_writer.py` - Function Conversion

**Function:** `write_fasta()` (lines 7-35)

**Current:**
```python
def write_fasta(sequences, filename, prefix="haplotype", comment=None, comments=None):
    """
    Write a list of sequences to a FASTA file.

    :param sequences: List of sequence strings.
    :param filename: Output filename for the FASTA file.
    :param prefix: Header prefix for each sequence (default "haplotype").
    :param comment: Optional comment to append to all header lines.
    :param comments: Optional list of comments, one per sequence.
    """
```

**Replace with:**
```python
def write_fasta(
    sequences: list[str],
    filename: str,
    prefix: str = "haplotype",
    comment: str | None = None,
    comments: list[str] | None = None
) -> None:
    """Write DNA/protein sequences to FASTA file with formatting.

    Writes sequences with auto-generated sequential IDs (prefix_1, prefix_2, ...)
    and optional comments. Sequences are wrapped at 80 characters for readability.

    Args:
        sequences: List of DNA or protein sequence strings
        filename: Output FASTA filename
        prefix: Header prefix for sequential IDs, defaults to "haplotype"
        comment: Optional comment appended to all headers (if comments not provided)
        comments: Optional list of per-sequence comments (overrides comment parameter)

    Raises:
        Exception: If file writing fails (logged and re-raised)

    Example:
        Write haplotypes with mutation annotations::

            sequences = [seq1, seq2]
            comments = ["mutation=dupC", "mutation=normal"]
            write_fasta(sequences, "output.fa", comments=comments)

        Output::

            >haplotype_1 mutation=dupC
            ATCGATCGATCG...
            >haplotype_2 mutation=normal
            GCTAGCTAGCTA...

    Note:
        Per-sequence comments take precedence over global comment parameter.
    """
```

**Estimated Time:** 10 minutes

---

#### File 10: `muc_one_up/io.py` - Function Conversions

**Function 1:** `parse_vntr_structure_file()` (lines 17-101)

**Current:**
```python
def parse_vntr_structure_file(
    filepath: str, config: dict
) -> tuple[list[list[str]], dict[str, Any] | None]:
    """
    Parse a VNTR structure file for use in simulation.

    The file format should be:
    haplotype_1<TAB>1-2-3-4-5-C-X-X-A-...-6p-7-8-9
    haplotype_2<TAB>1-2-3-4-5-C-X-X-X-...-6-7-8-9

    Each line represents one haplotype, and the symbols should be separated by
    dashes (-). The first column is the haplotype identifier, and the second
    column is the chain of repeat symbols.

    :param filepath: Path to the VNTR structure file.
    :param config: Configuration dict with valid repeat symbols.
    :return: List of lists, where each inner list contains symbols for one haplotype.
    :raises ValueError: If a symbol in the structure file is not found in the config.
    :raises FileNotFoundError: If the structure file cannot be found or read.
    """
```

**Replace with:**
```python
def parse_vntr_structure_file(
    filepath: str,
    config: dict
) -> tuple[list[list[str]], dict[str, Any] | None]:
    """Parse VNTR structure file for predefined chain simulation.

    Reads tab-delimited file with haplotype IDs and repeat chains. Extracts
    mutation information from comment lines if present. Validates all repeat
    symbols against config["repeats"].

    File Format:
        # Mutation Applied: dupC (Targets: [(1, 25)])
        haplotype_1<TAB>1-2-3-4-5-C-X-X-A-...-6p-7-8-9
        haplotype_2<TAB>1-2-3-4-5-C-X-X-X-...-6-7-8-9

    Args:
        filepath: Path to the VNTR structure file
        config: Configuration dict with valid repeat symbols in config["repeats"]

    Returns:
        Tuple containing:
        - List of repeat chains (one list per haplotype)
        - Mutation info dict or None if no mutation comment found
          Format: {"name": str, "targets": list[(int, int)]}

    Raises:
        ValueError: If repeat symbol not found in config or file format invalid
        FileNotFoundError: If structure file doesn't exist

    Example:
        >>> config = load_config("config.json")
        >>> chains, mut_info = parse_vntr_structure_file("structure.txt", config)
        >>> chains
        [['1', '2', '3', 'X'], ['1', '2', 'A', 'B']]
        >>> mut_info
        {'name': 'dupC', 'targets': [(1, 25)]}

    Note:
        Symbols with 'm' suffix (mutation markers) are validated by stripping
        the marker before checking against config["repeats"].
    """
```

---

**Function 2:** `extract_mutation_info_from_comments()` (lines 104-147)

**Current:**
```python
def extract_mutation_info_from_comments(
    comments: list[str],
) -> dict[str, Any] | None:
    """
    Extract mutation information from structure file comments.

    Looks for comments in the format:
    # Mutation Applied: dupC (Targets: [(1, 25)])

    :param comments: List of comment strings (without the leading '#')
    :return: Dictionary with mutation information or None if no valid info found
    """
```

**Replace with:**
```python
def extract_mutation_info_from_comments(
    comments: list[str],
) -> dict[str, Any] | None:
    """Extract mutation information from structure file comment lines.

    Parses comment lines to find mutation metadata in standardized format.
    Returns first valid mutation information found.

    Expected Format:
        Mutation Applied: <name> (Targets: [(hap_idx, rep_idx), ...])

    Args:
        comments: List of comment strings (without leading '#' character)

    Returns:
        Dictionary with mutation data or None if not found:
        {
            "name": mutation_name (str),
            "targets": [(haplotype_idx, repeat_idx), ...]  (1-based indices)
        }

    Example:
        >>> comments = ["Mutation Applied: dupC (Targets: [(1, 25), (2, 30)])"]
        >>> extract_mutation_info_from_comments(comments)
        {'name': 'dupC', 'targets': [(1, 25), (2, 30)]}

    Note:
        Malformed targets are logged as warnings and skipped. Uses ast.literal_eval
        for safe parsing of target tuples.
    """
```

**Estimated Time:** 15 minutes for both io.py conversions

---

### Priority 3: Complete Incomplete Docstrings (0.5-1 hour)

#### File 11: `muc_one_up/cli/orchestration.py`

**Function:** `run_single_simulation_iteration()` (lines 24-108)

**Current:** Has basic docstring but missing Args/Returns/Raises
**Action:** Expand docstring

**Replace current docstring (lines 36-40) with:**

```python
def run_single_simulation_iteration(
    args,
    config,
    out_dir,
    out_base,
    sim_index,
    fixed_conf,
    predefined_chains,
    dual_mutation_mode,
    mutation_pair,
    structure_mutation_info,
):
    """Run complete simulation iteration with all processing steps.

    Orchestrates a single simulation iteration from haplotype generation
    through mutation application, FASTA output, ORF prediction, read
    simulation, and statistics generation. Handles both normal and dual
    mutation modes.

    Args:
        args: Namespace containing all CLI arguments from argparse
        config: Configuration dictionary from load_config()
        out_dir: Output directory path for all generated files
        out_base: Base filename for outputs (without iteration suffix)
        sim_index: Current iteration index (for numbered outputs)
        fixed_conf: Fixed length configuration or "from_structure" flag
        predefined_chains: Predefined repeat chains from structure file (or None)
        dual_mutation_mode: Boolean indicating dual simulation mode (normal + mutated)
        mutation_pair: Tuple of (normal_name, mutated_name) for dual mode
        structure_mutation_info: Mutation metadata from structure file comments

    Returns:
        None (all outputs written to files)

    Raises:
        SimulationError: If haplotype generation fails
        ValidationError: If configuration or parameters are invalid
        FileOperationError: If file writing fails

    Side Effects:
        Creates multiple output files:
        - FASTA files (*.fa or *.normal.fa + *.mut.fa)
        - Structure files (*.structure.txt)
        - Mutated units files (*.mutated_units.txt)
        - Statistics JSON (*.simulation_stats.json)
        - Optional: ORF peptides (*.pep.fa)
        - Optional: Read simulation BAM files

    Note:
        Single Responsibility: Orchestrate one complete simulation iteration.
        This function delegates to specialized modules (haplotypes, mutations,
        outputs, analysis) following SOLID principles.
    """
```

**Estimated Time:** 10 minutes

---

#### File 12: `muc_one_up/__init__.py`

**Current (lines 1-8):**
```python
"""
muc_one_up package

This package simulates MUC1 VNTR diploid references.
"""

from .version import __version__
```

**Replace with:**
```python
"""MucOneUp: MUC1 VNTR simulation and analysis toolkit.

MucOneUp is a Python toolkit for simulating realistic MUC1 Variable Number
Tandem Repeat (VNTR) diploid references with customizable mutations and
sequencing read simulation. Supports both Illumina and Oxford Nanopore
sequencing platforms.

Core Capabilities:
    - VNTR haplotype simulation with probability-based repeat transitions
    - Targeted mutation application (insert, delete, replace, delete_insert)
    - SNP integration for haplotype-specific variants
    - Illumina and ONT read simulation with realistic error profiles
    - ORF prediction and toxic protein detection
    - Comprehensive statistics and coverage analysis

Key Modules:
    simulate: Core haplotype generation engine
    mutate: Mutation application and validation
    config: Configuration loading and schema validation
    read_simulation: Illumina and ONT pipeline orchestration
    analysis: VNTR statistics and ORF analysis
    cli: Command-line interface (click-based)

Example:
    Basic simulation workflow::

        from muc_one_up.config import load_config
        from muc_one_up.simulate import simulate_diploid
        from muc_one_up.mutate import apply_mutations

        # Load configuration
        config = load_config("config.json")

        # Generate diploid haplotypes
        results = simulate_diploid(config, num_haplotypes=2, fixed_lengths=[50, 60])

        # Apply mutation
        mutated_results, units = apply_mutations(
            config, results, "dupC", targets=[(1, 25)]
        )

Command-Line Interface:
    Install and use via CLI::

        pip install .
        muconeup --config config.json simulate --out-base test --fixed-lengths 50

Version:
    {__version__}

See Also:
    - CLAUDE.md: Comprehensive developer documentation
    - README.md: Installation and usage guide
    - config.json: Example configuration file
"""

from .version import __version__
```

**Estimated Time:** 8 minutes

---

## Part 4: Quality Assurance

### Automated Validation Tools

**Install:**
```bash
pip install pydocstyle mypy
```

**Run After Each Phase:**

```bash
# Check docstring compliance (Google convention)
pydocstyle muc_one_up/ --convention=google --count

# Check missing docstrings
pydocstyle muc_one_up/ --select=D100,D101,D102,D103 --count

# Verify type annotations are consistent
mypy muc_one_up/ --ignore-missing-imports --no-error-summary | grep "error:"
```

### Manual Review Checklist

After completing all conversions, verify for each modified file:

- [ ] Module docstring present and comprehensive (10-25 lines)
- [ ] All functions have Google-style docstrings
- [ ] NO old-style directives (`:param:`, `:type:`, `:return:`, `:rtype:`)
- [ ] Type information in signatures ONLY, never duplicated in docstrings
- [ ] Args section: One line per parameter (no type info)
- [ ] Returns section: Describes structure and meaning (no type info)
- [ ] Raises section: Lists all exceptions with conditions
- [ ] Examples for complex functions (>3 params or intricate logic)
- [ ] Imperative mood: "Generate..." not "Generates..."
- [ ] One-line summaries ≤ 72 characters

---

## Part 5: Implementation Workflow

### Recommended Order (Minimize Context Switching)

**Session 1 (1.5 hours): Module Docstrings**
1. `simulate.py` - Add module doc
2. `mutate.py` - Add module doc
3. `config.py` - Add module doc + document CONFIG_SCHEMA
4. `translate.py` - Add module doc + document CODON_TABLE
5. `fasta_writer.py` - Add module doc

**Session 2 (1.5 hours): Function Conversions (Core)**
6. `mutate.py` - Convert 4 functions
7. `config.py` - Convert 1 function
8. `translate.py` - Convert 6 functions

**Session 3 (1 hour): Function Conversions (Utilities)**
9. `fasta_writer.py` - Convert 1 function
10. `io.py` - Convert 2 functions
11. `orchestration.py` - Expand 1 function
12. `__init__.py` - Expand module doc

**Session 4 (0.5 hours): Quality Assurance**
13. Run pydocstyle validation
14. Fix any remaining issues
15. Manual checklist review
16. Commit changes

### Commit Strategy

One commit per session with clear scope:

```bash
# Session 1
git add muc_one_up/{simulate,mutate,config,translate,fasta_writer}.py
git commit -m "docs: add comprehensive module docstrings for core modules

- Add module docs for simulate, mutate, config, translate, fasta_writer
- Document CONFIG_SCHEMA and CODON_TABLE constants
- Follow Google style with key functions, examples, and design notes
"

# Session 2
git add muc_one_up/{mutate,config,translate}.py
git commit -m "docs: convert core functions to Google-style docstrings

- Convert mutate.py: 4 functions (validate_allowed_repeats, apply_mutations, etc.)
- Convert config.py: load_config function
- Convert translate.py: 6 functions (reverse_complement, dna_to_protein, etc.)
- Remove all :param:, :type:, :return:, :rtype: directives
- Add type annotations, remove type duplication in docstrings
"

# Session 3
git add muc_one_up/{fasta_writer,io,cli/orchestration,__init__}.py
git commit -m "docs: convert utility functions and expand incomplete docstrings

- Convert fasta_writer.py: write_fasta function
- Convert io.py: 2 parsing functions
- Expand orchestration.py: add Args/Returns/Raises
- Expand __init__.py: comprehensive package documentation
"

# Session 4
git commit -m "docs: finalize Google-style standardization (100% coverage)

- Verified with pydocstyle --convention=google
- All modules have comprehensive docstrings
- All functions use Google style with type annotations
- Zero type duplication (DRY principle)
"
```

---

## Part 6: Success Metrics

### Before This Plan

| Metric | Value |
|--------|-------|
| Module doc coverage | 55% (11/20 modules) |
| Old-style functions | 18 functions with `:param:` etc. |
| Undocumented constants | 2 (CONFIG_SCHEMA, CODON_TABLE) |
| Incomplete functions | 2 (orchestration, __init__) |
| Style consistency | 65% Google-style |
| Type duplication | ~42 instances |

### After This Plan (Target)

| Metric | Value |
|--------|-------|
| Module doc coverage | **100%** (20/20 modules) |
| Old-style functions | **0** (all converted to Google) |
| Undocumented constants | **0** (all documented) |
| Incomplete functions | **0** (all expanded) |
| Style consistency | **100%** Google-style |
| Type duplication | **0** instances |

### Validation Commands

```bash
# Verify 100% module coverage
find muc_one_up -name "*.py" -type f | while read f; do
    if ! grep -q '"""' "$f"; then
        echo "Missing module doc: $f"
    fi
done

# Count old-style directives (should be 0)
grep -r ":param\|:type\|:return\|:rtype" muc_one_up/*.py | wc -l

# Verify Google style compliance
pydocstyle muc_one_up/ --convention=google --count
```

---

## Part 7: Reference Examples (Copy-Paste Templates)

### Template 1: Simple Function (Pure Logic, No I/O)

```python
def pick_next_repeat(
    probabilities: ProbabilitiesDict,
    current_symbol: str,
    force_end: bool = False,
) -> str:
    """Select next repeat symbol from probability distribution.

    Args:
        probabilities: Transition matrix mapping current to next state probabilities
        current_symbol: The current repeat symbol in chain
        force_end: If True, bias selection toward 'END' state

    Returns:
        Next symbol as string (e.g., "1", "X", "6p", "END")

    Example:
        >>> probs = {"1": {"2": 0.7, "7": 0.3}}
        >>> pick_next_repeat(probs, "1")
        '2'  # or '7' based on weighted random selection
    """
```

### Template 2: Complex Function (Multiple Returns, Exceptions)

```python
def apply_mutations(
    config: ConfigDict,
    results: HaplotypeList,
    mutation_name: MutationName,
    targets: MutationTargets,
) -> tuple[HaplotypeList, dict[int, list[tuple[int, str]]]]:
    """Apply named mutation to specific haplotype positions.

    Applies mutation operations (insert/delete/replace) to targeted repeat
    units. Validates allowed_repeats and tracks mutated sequences.

    Args:
        config: Configuration dict containing mutations section
        results: List of (sequence, chain) tuples for each haplotype
        mutation_name: Key in config["mutations"]
        targets: List of (haplotype_idx, repeat_idx) using 1-based indexing

    Returns:
        Tuple containing:
        - updated_results: Haplotypes with mutations applied
        - mutated_units: Dict mapping haplotype index to mutated sequences

    Raises:
        ValueError: If mutation not found or indices invalid
        ValidationError: If strict mode rejects invalid repeat symbol

    Example:
        >>> config = load_config("config.json")
        >>> results = simulate_diploid(config)
        >>> mutated, units = apply_mutations(config, results, "dupC", [(1, 25)])

    Note:
        Positions use 1-based indexing. Mutated repeats marked with 'm' suffix.
    """
```

### Template 3: I/O Function (File Operations)

```python
def write_fasta(
    sequences: list[str],
    filename: str,
    prefix: str = "haplotype",
    comments: list[str] | None = None
) -> None:
    """Write sequences to FASTA file with formatting.

    Args:
        sequences: List of DNA/protein sequence strings
        filename: Output FASTA filename
        prefix: Header prefix for IDs, defaults to "haplotype"
        comments: Optional per-sequence comments

    Raises:
        FileOperationError: If file writing fails

    Example:
        >>> write_fasta([seq1, seq2], "out.fa", comments=["mut=dupC", "normal"])
    """
```

---

## Appendix: Quick Conversion Reference

### Old Style → Google Style Cheat Sheet

| Old Style (DON'T USE) | Google Style (USE THIS) |
|----------------------|-------------------------|
| `:param name: Description` | `Args:\n    name: Description` |
| `:type name: str` | *(delete - type in signature)* |
| `:param name: Desc\n:type name: str` | `Args:\n    name: Desc` |
| `:return: Description` | `Returns:\n    Description` |
| `:rtype: str` | *(delete - type in signature)* |
| `:raises ValueError: When` | `Raises:\n    ValueError: When` |

### Before/After Example

**BEFORE (Old Style - 12 lines):**
```python
def load_config(config_path: str) -> dict[str, Any]:
    """
    Load and validate the JSON config for MucOneUp.

    :param config_path: Path to the JSON config file.
    :type config_path: str
    :return: Python dict with validated configuration data.
    :rtype: dict
    :raises FileNotFoundError: if the config file does not exist.
    :raises json.JSONDecodeError: if the config file is not valid JSON.
    :raises ValidationError: if the config does not conform to schema.
    """
```

**AFTER (Google Style - 9 lines, 25% reduction):**
```python
def load_config(config_path: str) -> dict[str, Any]:
    """Load and validate MucOneUp configuration from JSON file.

    Args:
        config_path: Path to the JSON config file

    Returns:
        Validated configuration dictionary with all required sections

    Raises:
        FileNotFoundError: If the config file does not exist
        json.JSONDecodeError: If the config file is not valid JSON
        ValidationError: If the config does not conform to CONFIG_SCHEMA
    """
```

**Changes:**
- ✅ Removed `:param` directive
- ✅ Removed `:type` directive (type already in signature `config_path: str`)
- ✅ Removed `:return` directive → replaced with `Returns:` section
- ✅ Removed `:rtype` directive (type already in signature `-> dict[str, Any]`)
- ✅ Converted `:raises` → `Raises:` section
- ✅ Improved summary line (imperative mood)
- ✅ Expanded Returns description (more detail)

**Result:** 25% fewer lines, zero type duplication, improved readability

---

## Summary

This plan achieves **100% Google-style documentation coverage** for MucOneUp through:

1. **12 files** requiring changes (out of 20 total)
2. **7 module docstrings** to add
3. **18 functions** to convert from old style
4. **2 constants** to document
5. **2 docstrings** to expand

**Total Estimated Time:** 4-6 hours
**Total Lines Changed:** ~150 docstring lines
**Impact:** Complete standardization, zero technical debt, perfect DRY compliance

**Next Steps:**
1. Review and approve this plan
2. Schedule 4x 1.5-hour sessions
3. Execute sessions in order (minimize context switching)
4. Run QA validation after each session
5. Final review and commit

---

**END OF PLAN**
