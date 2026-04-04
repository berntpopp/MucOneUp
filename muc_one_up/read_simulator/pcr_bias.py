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
        """Load a named preset with optional parameter overrides."""
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
        """
        if not config:
            return cls.from_preset(DEFAULT_PCR_PRESET)

        preset = config.get("preset")
        if preset is not None:
            overrides = {k: v for k, v in config.items() if k != "preset"}
            return cls.from_preset(preset, **overrides)

        defaults = _PRESETS[DEFAULT_PCR_PRESET]
        return cls(
            e_max=config.get("e_max", defaults["e_max"]),
            alpha=config.get("alpha", defaults["alpha"]),
            cycles=config.get("cycles", defaults["cycles"]),
            denaturation_time=config.get("denaturation_time", defaults["denaturation_time"]),
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

    def _deterministic_split(self, total: int, e1: float, e2: float) -> tuple[int, int]:
        """Deterministic coverage split based on yield ratio."""
        yield1 = (1 + e1) ** self.cycles
        yield2 = (1 + e2) ** self.cycles
        total_yield = yield1 + yield2

        frac1 = yield1 / total_yield
        n1 = round(total * frac1)
        n2 = total - n1

        return n1, n2

    def _stochastic_split(
        self, total: int, e1: float, e2: float, seed: int | None
    ) -> tuple[int, int]:
        """Stochastic coverage split using Galton-Watson branching process."""
        import numpy as np

        rng = np.random.default_rng(seed)

        n1 = 1000
        n2 = 1000

        for _ in range(self.cycles):
            n1 += rng.binomial(n1, e1)
            n2 += rng.binomial(n2, e2)

        frac1 = n1 / (n1 + n2)
        reads1 = round(total * frac1)
        reads2 = total - reads1

        return reads1, reads2
