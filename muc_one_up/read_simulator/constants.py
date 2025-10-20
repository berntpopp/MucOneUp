# muc_one_up/read_simulator/constants.py
"""Constants for read simulation pipelines.

This module centralizes magic numbers, validation constants, and default
values used across Illumina, ONT, and PacBio read simulation pipelines.

Key Constants:
    VALID_PBSIM3_MODEL_TYPES: Valid pbsim3 model types (qshmm, errhmm)
    MINIMAP2_PRESET_*: Minimap2 alignment presets for different technologies
    DEFAULT_*_TIMEOUT: Default timeout values for external tools
    MIN/MAX_COVERAGE: Valid coverage range for read simulation

Example:
    Importing constants:

        from muc_one_up.read_simulator.constants import (
            VALID_PBSIM3_MODEL_TYPES,
            MINIMAP2_PRESET_PACBIO_HIFI,
            DEFAULT_PBSIM3_TIMEOUT,
        )

Notes:
    - Centralization prevents magic numbers and improves maintainability
    - Used by both wrappers (runtime validation) and JSON Schema (config validation)
    - Timeout values are conservative for large genomes and high coverage

References:
    - pbsim3 model types: https://github.com/yukiteruono/pbsim3#model-based-simulation
    - minimap2 presets: https://github.com/lh3/minimap2#usage
"""

from typing import Final

# =============================================================================
# pbsim3 Constants
# =============================================================================

#: Valid pbsim3 model types for PacBio simulation.
#:
#: - qshmm: Quality Score Hidden Markov Model (recommended for most use cases)
#: - errhmm: Error Hidden Markov Model (alternative model type)
#:
#: Reference: https://github.com/yukiteruono/pbsim3#model-based-simulation
VALID_PBSIM3_MODEL_TYPES: Final[set[str]] = {"qshmm", "errhmm"}

#: Default number of passes for multi-pass CLR simulation.
#: HiFi reads require ≥2 passes for high accuracy consensus.
DEFAULT_PBSIM3_PASS_NUM: Final[int] = 3

#: Default accuracy mean for multi-pass CLR simulation (0.0-1.0).
#: Represents mean base accuracy before consensus generation.
DEFAULT_PBSIM3_ACCURACY_MEAN: Final[float] = 0.85

#: Default accuracy standard deviation for CLR simulation (0.0-1.0).
DEFAULT_PBSIM3_ACCURACY_SD: Final[float] = 0.05

#: Default minimum accuracy for CLR simulation (0.0-1.0).
DEFAULT_PBSIM3_ACCURACY_MIN: Final[float] = 0.75

#: Default length mean for CLR simulation (base pairs).
#: PacBio Sequel II produces ~15-20kb reads on average.
DEFAULT_PBSIM3_LENGTH_MEAN: Final[int] = 15000

#: Default length standard deviation for CLR simulation (base pairs).
DEFAULT_PBSIM3_LENGTH_SD: Final[int] = 5000

#: Default minimum length for CLR simulation (base pairs).
DEFAULT_PBSIM3_LENGTH_MIN: Final[int] = 5000

#: Default maximum length for CLR simulation (base pairs).
DEFAULT_PBSIM3_LENGTH_MAX: Final[int] = 30000

# =============================================================================
# CCS (Circular Consensus Sequencing) Constants
# =============================================================================

#: Default minimum number of passes required for CCS HiFi generation.
#: Higher values increase accuracy but reduce yield.
DEFAULT_CCS_MIN_PASSES: Final[int] = 3

#: Default minimum predicted accuracy for CCS HiFi reads (0.0-1.0).
#: Standard HiFi threshold is 0.99 (Q20).
DEFAULT_CCS_MIN_RQ: Final[float] = 0.99

#: Default number of threads for CCS processing.
DEFAULT_CCS_THREADS: Final[int] = 4

#: Minimum file size (bytes) for valid CCS output BAM.
#: Empty or corrupted BAM files are typically <1KB.
MIN_VALID_CCS_OUTPUT_SIZE: Final[int] = 1024

# =============================================================================
# Minimap2 Alignment Presets
# =============================================================================

#: Minimap2 preset for Oxford Nanopore reads.
#: Optimized for ONT error profile with higher indel rate.
MINIMAP2_PRESET_ONT: Final[str] = "map-ont"

#: Minimap2 preset for PacBio CLR reads.
#: Optimized for older PacBio Continuous Long Reads.
MINIMAP2_PRESET_PACBIO_CLR: Final[str] = "map-pb"

#: Minimap2 preset for PacBio HiFi reads.
#: Optimized for high-accuracy HiFi reads (Q20+).
MINIMAP2_PRESET_PACBIO_HIFI: Final[str] = "map-hifi"

#: Default number of threads for minimap2 alignment.
DEFAULT_MINIMAP2_THREADS: Final[int] = 4

# =============================================================================
# Timeout Constants (seconds)
# =============================================================================

#: Default timeout for pbsim3 simulation (2 hours).
#: Large genomes with high coverage may require longer timeouts.
DEFAULT_PBSIM3_TIMEOUT: Final[int] = 7200

#: Default timeout for CCS consensus generation (2 hours).
#: Processing time scales with number of reads and passes.
DEFAULT_CCS_TIMEOUT: Final[int] = 7200

#: Default timeout for minimap2 alignment (1 hour).
#: Alignment is typically faster than simulation/consensus.
DEFAULT_ALIGNMENT_TIMEOUT: Final[int] = 3600

#: Default timeout for samtools operations (30 minutes).
#: SAM/BAM conversion and indexing are usually quick.
DEFAULT_SAMTOOLS_TIMEOUT: Final[int] = 1800

# =============================================================================
# Coverage and Quality Limits
# =============================================================================

#: Minimum valid coverage for read simulation.
#: Values below 0.1x may produce insufficient reads for analysis.
MIN_COVERAGE: Final[float] = 0.1

#: Maximum valid coverage for read simulation.
#: Very high coverage (>1000x) can cause memory issues and slow processing.
MAX_COVERAGE: Final[float] = 10000.0

#: Minimum valid quality score (Phred scale).
MIN_QUALITY_SCORE: Final[int] = 0

#: Maximum valid quality score (Phred scale).
#: Q93 represents 99.9999999% accuracy (theoretical maximum).
MAX_QUALITY_SCORE: Final[int] = 93

# =============================================================================
# File Size Limits
# =============================================================================

#: Minimum reference sequence length (base pairs).
#: Very short sequences may not simulate realistically.
MIN_REFERENCE_LENGTH: Final[int] = 100

#: Maximum reference sequence length (base pairs).
#: Protects against accidentally simulating entire genomes.
MAX_REFERENCE_LENGTH: Final[int] = 1_000_000_000  # 1 Gbp

# =============================================================================
# Pass Number Limits (Multi-pass Simulation)
# =============================================================================

#: Minimum number of passes for multi-pass CLR simulation.
#: HiFi requires ≥2 passes; 1 pass would be equivalent to CLR.
MIN_PASS_NUM: Final[int] = 2

#: Maximum number of passes for multi-pass CLR simulation.
#: Very high pass counts are unrealistic and slow down simulation.
MAX_PASS_NUM: Final[int] = 50

# =============================================================================
# Accuracy Limits (0.0-1.0 scale)
# =============================================================================

#: Minimum valid accuracy value.
MIN_ACCURACY: Final[float] = 0.0

#: Maximum valid accuracy value.
MAX_ACCURACY: Final[float] = 1.0
