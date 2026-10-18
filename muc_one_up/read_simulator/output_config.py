"""Output configuration for read simulation pipelines."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class OutputConfig:
    """Configuration for output file placement."""

    out_dir: Path
    out_base: str

    def __post_init__(self) -> None:
        """Validate output configuration."""
        if not self.out_base or not self.out_base.strip():
            raise ValueError("out_base cannot be empty")
        if (
            "/" in self.out_base
            or "\\" in self.out_base
            or self.out_base in (".", "..")
            or Path(self.out_base).name != self.out_base
            or Path(self.out_base).is_absolute()
        ):
            raise ValueError(
                f"Invalid out_base: '{self.out_base}'. Must be a single filename component without path separators or directory traversal."
            )

    def derive_path(self, suffix: str) -> Path:
        """Return out_dir / (out_base + suffix)."""
        return self.out_dir / f"{self.out_base}{suffix}"

    def derive_path_str(self, suffix: str) -> str:
        """Return string version of derive_path."""
        return str(self.derive_path(suffix))

    def ensure_dir(self) -> None:
        """Create out_dir (and parents) if it does not exist."""
        self.out_dir.mkdir(parents=True, exist_ok=True)

    @classmethod
    def from_input_fasta(
        cls,
        input_fa: str,
        out_dir: str | None,
        out_base: str | None,
        suffix: str,
    ) -> OutputConfig:
        """Derive OutputConfig from an input FASTA path with optional overrides.

        Parameters
        ----------
        input_fa : str
            Path to the input FASTA file.
        out_dir : str | None
            Output directory override. If None, uses input_fa's parent.
        out_base : str | None
            Output base name override. If None, uses input_fa stem + suffix.
        suffix : str
            Suffix appended to stem when out_base is not provided.

        Returns
        -------
        OutputConfig
            Resolved output configuration.
        """
        input_path = Path(input_fa)
        resolved_dir = Path(out_dir) if out_dir else input_path.parent
        resolved_base = out_base if out_base else input_path.stem + suffix
        return cls(out_dir=resolved_dir, out_base=resolved_base)
