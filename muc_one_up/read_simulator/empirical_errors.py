"""Empirical long-read error channel calibrated to real-data aggregate statistics.

pbsim3's ONT models cannot reach R10.4.1 sup accuracy in template mode (error
floor ~3.5-4% with ~16-22% C7 length errors, against ~2.1% and 11-48% measured
per strand). This channel reproduces measured rates directly:

* mismatches, insertions and deletions at per-base rates outside homopolymer
  runs (runs >= ``hp_min_len`` are left to the molecule model's stutter table,
  so the table equals the measured homopolymer spectrum without deconvolution);
* indel lengths from measured distributions;
* per-read error variation as a log-normal multiplier (mean 1);
* Phred qualities from the read's error rate, scaled by ``q_scale`` because
  basecaller qualities over-promise (observed error ~1.4x Q-predicted).

Pure Python with an injected RNG; deterministic per seed.
"""

from __future__ import annotations

import math
import random
from collections import Counter
from collections.abc import Mapping
from dataclasses import dataclass, fields

from .molecules import homopolymer_runs
from .stutter import MIN_RUN_LEN

_BASES = "ACGT"
_MAX_RATE = 0.2


def _pmf(
    data: Mapping[int, float] | Mapping[str, float], name: str
) -> tuple[tuple[int, float], ...]:
    items = tuple(sorted((int(k), float(v)) for k, v in data.items()))
    if (
        not items
        or any(k < 1 or p < 0 for k, p in items)
        or not math.isclose(sum(p for _, p in items), 1.0, abs_tol=1e-6)
    ):
        raise ValueError(f"{name} must map lengths >= 1 to probabilities that sum to 1")
    return items


def _draw(pmf: tuple[tuple[int, float], ...], rng: random.Random) -> int:
    draw, cumulative = rng.random(), 0.0
    for value, prob in pmf:
        cumulative += prob
        if draw < cumulative:
            return value
    return pmf[-1][0]


@dataclass(frozen=True)
class EmpiricalErrorModel:
    """Per-base error rates for non-homopolymer positions plus quality model."""

    mismatch_rate: float
    insertion_rate: float
    deletion_rate: float
    insertion_len_pmf: Mapping[int, float]
    deletion_len_pmf: Mapping[int, float]
    read_error_sigma: float = 0.45
    q_scale: float = 0.7
    hp_min_len: int = MIN_RUN_LEN

    def __post_init__(self) -> None:
        for name in ("mismatch_rate", "insertion_rate", "deletion_rate"):
            value = getattr(self, name)
            if not 0.0 <= value <= _MAX_RATE:
                raise ValueError(f"{name} must be in [0, {_MAX_RATE}], got {value}")
        if self.read_error_sigma < 0 or not 0 < self.q_scale <= 1 or self.hp_min_len < 2:
            raise ValueError(
                "read_error_sigma >= 0, q_scale in (0, 1] and hp_min_len >= 2 required"
            )
        object.__setattr__(
            self, "insertion_len_pmf", dict(_pmf(self.insertion_len_pmf, "insertion_len_pmf"))
        )
        object.__setattr__(
            self, "deletion_len_pmf", dict(_pmf(self.deletion_len_pmf, "deletion_len_pmf"))
        )

    @classmethod
    def from_dict(cls, data: Mapping[str, object]) -> EmpiricalErrorModel:
        unknown = set(data) - {f.name for f in fields(cls)}
        if unknown:
            raise ValueError(f"unknown error model keys: {sorted(unknown)}")
        return cls(**data)  # type: ignore[arg-type]

    @property
    def total_rate(self) -> float:
        return self.mismatch_rate + self.insertion_rate + self.deletion_rate


def _read_multiplier(model: EmpiricalErrorModel, rng: random.Random) -> float:
    if model.read_error_sigma == 0:
        return 1.0
    sigma = model.read_error_sigma
    return math.exp(rng.gauss(-sigma * sigma / 2.0, sigma))  # log-normal with mean 1


def _quality(rate: float, length: int, model: EmpiricalErrorModel, rng: random.Random) -> str:
    phred = -10.0 * math.log10(max(rate * model.q_scale, 1e-4))
    return "".join(chr(33 + min(40, max(2, round(rng.gauss(phred, 2.0))))) for _ in range(length))


def apply_errors(
    seq: str, model: EmpiricalErrorModel, rng: random.Random, events: Counter[str] | None = None
) -> tuple[str, str]:
    """Return (read, qualities) for one molecule sequence.

    If ``events`` is given, it is incremented with ``mismatch``, ``insertion`` and
    ``deletion`` event counts and ``exposed`` (non-homopolymer bases).
    """
    events = events if events is not None else Counter()
    protected = bytearray(len(seq))
    for start, end, _ in homopolymer_runs(seq, model.hp_min_len):
        protected[start:end] = b"\x01" * (end - start)
    scale = _read_multiplier(model, rng)
    p_del = min(_MAX_RATE, model.deletion_rate * scale)
    p_ins = p_del + min(_MAX_RATE, model.insertion_rate * scale)
    p_mm = p_ins + min(_MAX_RATE, model.mismatch_rate * scale)
    ins_pmf = tuple(model.insertion_len_pmf.items())
    del_pmf = tuple(model.deletion_len_pmf.items())
    out: list[str] = []
    i = 0
    while i < len(seq):
        base = seq[i]
        if protected[i]:
            out.append(base)
            i += 1
            continue
        events["exposed"] += 1
        draw = rng.random()
        if draw < p_del:
            events["deletion"] += 1
            # Stop at the next protected base: homopolymer lengths are set by
            # the stutter table only.
            stop = min(len(seq), i + _draw(del_pmf, rng))
            i += 1
            while i < stop and not protected[i]:
                i += 1
            continue
        if draw < p_ins:
            events["insertion"] += 1
            out.append(base)
            out.extend(rng.choice(_BASES) for _ in range(_draw(ins_pmf, rng)))
        elif draw < p_mm:
            events["mismatch"] += 1
            out.append(rng.choice(_BASES.replace(base, "") or _BASES))
        else:
            out.append(base)
        i += 1
    read = "".join(out)
    return read, _quality(model.total_rate * scale, len(read), model, rng)
