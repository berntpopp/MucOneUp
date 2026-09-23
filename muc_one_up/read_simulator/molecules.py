"""Template molecule model for truth-tracked long-read simulation.

Builds the molecules that pbsim3 template mode turns into reads. Each molecule
records its haplotype, product kind, strand, source interval and any injected
homopolymer edits, so every simulated read has exact truth.

The module is pure: no I/O, no external tools, and all randomness comes from an
injected :class:`random.Random`. Default :class:`MoleculeModel` values are a
no-op (full-length, forward-strand molecules without edits).

Strand convention: ``"+"`` means the molecule is read in the source orientation
(for MUC1 haplotypes, motif 1 -> motif 9 with C-runs read as C). Stutter keys
use the base and run length in source orientation plus the read strand, e.g.
``"C7|+"``.
"""

from __future__ import annotations

import math
import random
from collections.abc import Iterator, Mapping, Sequence
from dataclasses import dataclass, field, fields
from typing import NamedTuple

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
_STRANDS = ("+", "-")


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return seq.translate(_COMPLEMENT)[::-1]


@dataclass(frozen=True)
class HpEdit:
    """One injected homopolymer length change (positions in source orientation)."""

    pos: int
    base: str
    true_len: int
    new_len: int


class SourceInterval(NamedTuple):
    """A contiguous stretch of haplotype ``hap`` starting at 0-based ``start``."""

    hap: int
    start: int
    seq: str


@dataclass(frozen=True)
class Molecule:
    """A simulated template molecule; ``seq`` is in sequencing orientation."""

    id: int
    hap: int  # 1-based haplotype of origin
    kind: str  # full | smear | chimera | concatemer | offtarget | fragment
    strand: str
    seq: str
    src_start: int
    src_end: int
    hp_edits: tuple[HpEdit, ...] = ()
    detail: str = ""


@dataclass(frozen=True)
class StutterTable:
    """Per-(base, length, strand) probability mass over homopolymer length deltas."""

    pmfs: Mapping[str, tuple[tuple[int, float], ...]]
    min_len: int = 3

    @classmethod
    def from_dict(cls, data: Mapping[str, Mapping[str, float]], min_len: int = 3) -> StutterTable:
        pmfs: dict[str, tuple[tuple[int, float], ...]] = {}
        for key, pmf in data.items():
            _validate_stutter_key(key)
            items = tuple(sorted((int(delta), float(p)) for delta, p in pmf.items()))
            if any(p < 0 for _, p in items) or not math.isclose(
                sum(p for _, p in items), 1.0, abs_tol=1e-6
            ):
                raise ValueError(f"stutter pmf for '{key}' must be non-negative and sum to 1")
            pmfs[key] = items
        return cls(pmfs, min_len)

    def sample(self, base: str, length: int, strand: str, rng: random.Random) -> int:
        """Draw a length delta; 0 when the table has no entry for this run."""
        pmf = self.pmfs.get(f"{base}{length}|{strand}") or self.pmfs.get(f"{base}{length}|both")
        if not pmf:
            return 0
        draw = rng.random()
        cumulative = 0.0
        for delta, prob in pmf:
            cumulative += prob
            if draw < cumulative:
                return delta
        return pmf[-1][0]


def _validate_stutter_key(key: str) -> None:
    run, sep, strand = key.partition("|")
    if (
        not sep
        or strand not in (*_STRANDS, "both")
        or len(run) < 2
        or run[0] not in "ACGT"
        or not run[1:].isdigit()
    ):
        raise ValueError(f"invalid stutter key '{key}': expected e.g. 'C7|+', 'G4|-' or 'A5|both'")


@dataclass(frozen=True)
class MoleculeModel:
    """Library/molecule artefact rates. Defaults reproduce legacy full-length reads."""

    forward_frac: float = 1.0
    smear_rate: float = 0.0
    smear_junction_beta: tuple[float, float] = (2.0, 5.0)
    smear_min_keep: float = 0.15
    chimera_rate: float = 0.0
    concatemer_rate: float = 0.0
    offtarget_frac: float = 0.0
    offtarget_median_bp: int = 370
    offtarget_sigma: float = 0.5
    stutter: StutterTable | None = field(default=None)

    def __post_init__(self) -> None:
        for name in (
            "forward_frac",
            "smear_rate",
            "smear_min_keep",
            "chimera_rate",
            "concatemer_rate",
        ):
            value = getattr(self, name)
            if not 0.0 <= value <= 1.0:
                raise ValueError(f"{name} must be in [0, 1], got {value}")
        if not 0.0 <= self.offtarget_frac < 1.0:
            raise ValueError(f"offtarget_frac must be in [0, 1), got {self.offtarget_frac}")
        if self.smear_rate + self.chimera_rate + self.concatemer_rate > 1.0:
            raise ValueError("smear_rate + chimera_rate + concatemer_rate must not exceed 1")
        if (
            min(self.smear_junction_beta) <= 0
            or self.offtarget_median_bp < 1
            or self.offtarget_sigma <= 0
        ):
            raise ValueError(
                "smear_junction_beta, offtarget_median_bp and offtarget_sigma must be positive"
            )

    @classmethod
    def from_dict(cls, data: Mapping[str, object]) -> MoleculeModel:
        """Build from a JSON-like mapping; unknown keys are rejected."""
        known = {f.name for f in fields(cls)}
        unknown = set(data) - known
        if unknown:
            raise ValueError(f"unknown molecule model keys: {sorted(unknown)}")
        kwargs = dict(data)
        if "stutter" in kwargs and kwargs["stutter"] is not None:
            kwargs["stutter"] = StutterTable.from_dict(kwargs["stutter"])  # type: ignore[arg-type]
        if "smear_junction_beta" in kwargs:
            kwargs["smear_junction_beta"] = tuple(kwargs["smear_junction_beta"])  # type: ignore[arg-type]
        return cls(**kwargs)  # type: ignore[arg-type]


def homopolymer_runs(seq: str, min_len: int) -> Iterator[tuple[int, int, str]]:
    """Yield (start, end, base) of runs of one base with length >= min_len."""
    i = 0
    while i < len(seq):
        j = i
        while j < len(seq) and seq[j] == seq[i]:
            j += 1
        if j - i >= min_len:
            yield i, j, seq[i]
        i = j


def apply_stutter(
    seq: str, table: StutterTable, strand: str, rng: random.Random
) -> tuple[str, tuple[HpEdit, ...]]:
    """Apply stutter to every homopolymer run of ``seq`` (source orientation)."""
    parts: list[str] = []
    edits: list[HpEdit] = []
    last = 0
    for start, end, base in homopolymer_runs(seq, table.min_len):
        delta = table.sample(base, end - start, strand, rng)
        if delta == 0:
            continue
        new_len = max(1, end - start + delta)
        parts.append(seq[last:start])
        parts.append(base * new_len)
        edits.append(HpEdit(start, base, end - start, new_len))
        last = end
    parts.append(seq[last:])
    return "".join(parts), tuple(edits)


def _finish(
    mol_id: int,
    hap: int,
    kind: str,
    source_seq: str,
    span: tuple[int, int],
    model: MoleculeModel,
    rng: random.Random,
    detail: str = "",
) -> Molecule:
    """Assign strand, apply stutter in source orientation, then orient for sequencing."""
    strand = "+" if rng.random() < model.forward_frac else "-"
    seq: str = source_seq
    edits: tuple[HpEdit, ...] = ()
    if model.stutter is not None:
        seq, edits = apply_stutter(source_seq, model.stutter, strand, rng)
    if strand == "-":
        seq = reverse_complement(seq)
    return Molecule(mol_id, hap, kind, strand, seq, span[0], span[1], edits, detail)


def _smear(seq: str, model: MoleculeModel, rng: random.Random) -> tuple[str, str]:
    a, b = model.smear_junction_beta
    start = int(rng.betavariate(a, b) * len(seq))
    keep = rng.uniform(model.smear_min_keep, 0.95)
    resume = min(len(seq), start + round((1.0 - keep) * len(seq)))
    return seq[:start] + seq[resume:], f"deletion:{start}-{resume}"


def _chimera(first: str, second: str, rng: random.Random) -> tuple[str, str]:
    frac = rng.uniform(0.1, 0.9)
    cut_a, cut_b = int(frac * len(first)), int(frac * len(second))
    return first[:cut_a] + second[cut_b:], f"junction:{cut_a}|{cut_b}"


def _choose_product(
    amplicon: str, partner: str | None, model: MoleculeModel, rng: random.Random
) -> tuple[str, str, str]:
    """Pick the PCR product kind for one molecule slot; returns (kind, seq, detail)."""
    draw = rng.random()
    if draw < model.smear_rate:
        seq, detail = _smear(amplicon, model, rng)
        return "smear", seq, detail
    draw -= model.smear_rate
    if draw < model.chimera_rate:
        if partner is None:  # chimeras need two haplotypes; the slot stays full-length
            return "full", amplicon, ""
        seq, detail = _chimera(amplicon, partner, rng)
        return "chimera", seq, detail
    draw -= model.chimera_rate
    if draw < model.concatemer_rate:
        return "concatemer", amplicon + amplicon, ""
    return "full", amplicon, ""


def _lognormal_length(median: float, sigma: float, rng: random.Random) -> int:
    return max(1, round(rng.lognormvariate(math.log(median), sigma)))


def build_amplicon_molecules(
    amplicons: Sequence[str],
    counts: Sequence[int],
    model: MoleculeModel,
    rng: random.Random,
    *,
    offtarget_sources: Sequence[SourceInterval] = (),
) -> list[Molecule]:
    """Build PCR amplicon molecules: ``counts[i]`` per haplotype plus artefacts.

    Each full-length slot becomes a smear, chimera (only with two haplotypes) or
    concatemer product with the configured probabilities. Off-target products
    are added on top (``offtarget_frac`` of the final total). Each is cut from
    one of ``offtarget_sources`` (chosen weighted by length), so it never
    crosses an interval end, and records that haplotype and its coordinates.
    """
    if len(amplicons) != len(counts):
        raise ValueError("amplicons and counts must have equal length")
    sources = [s for s in offtarget_sources if s.seq]
    if model.offtarget_frac and not sources:
        raise ValueError("offtarget_sources are required when offtarget_frac > 0")
    molecules: list[Molecule] = []
    for hap_index, (amplicon, count) in enumerate(zip(amplicons, counts, strict=True), start=1):
        partner = amplicons[1 - (hap_index - 1)] if len(amplicons) == 2 else None
        for _ in range(count):
            mol_id = len(molecules) + 1
            kind, seq, detail = _choose_product(amplicon, partner, model, rng)
            molecules.append(
                _finish(mol_id, hap_index, kind, seq, (0, len(amplicon)), model, rng, detail)
            )
    if model.offtarget_frac:
        n_off = round(len(molecules) * model.offtarget_frac / (1.0 - model.offtarget_frac))
        weights = [len(s.seq) for s in sources]
        for _ in range(n_off):
            source = rng.choices(sources, weights=weights)[0]
            length = min(
                len(source.seq),
                _lognormal_length(model.offtarget_median_bp, model.offtarget_sigma, rng),
            )
            start = rng.randrange(0, len(source.seq) - length + 1)
            piece = source.seq[start : start + length]
            span = (source.start + start, source.start + start + length)
            molecules.append(
                _finish(len(molecules) + 1, source.hap, "offtarget", piece, span, model, rng)
            )
    return molecules


def build_fragment_molecules(
    sources: Sequence[str],
    n_molecules: int,
    length_median: float,
    length_sigma: float,
    model: MoleculeModel,
    rng: random.Random,
) -> list[Molecule]:
    """Sample genomic fragments: haplotype uniform, log-normal length, uniform start."""
    if not sources or n_molecules < 0 or length_median <= 0 or length_sigma <= 0:
        raise ValueError("sources must be non-empty and length parameters positive")
    molecules: list[Molecule] = []
    for mol_id in range(1, n_molecules + 1):
        hap = rng.randrange(len(sources)) + 1
        source = sources[hap - 1]
        start = rng.randrange(len(source))
        end = min(len(source), start + _lognormal_length(length_median, length_sigma, rng))
        molecules.append(
            _finish(mol_id, hap, "fragment", source[start:end], (start, end), model, rng)
        )
    return molecules
