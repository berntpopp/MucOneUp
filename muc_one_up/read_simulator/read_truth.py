"""Per-read truth for pbsim3 template simulations.

pbsim3 template mode writes one read (or one ZMW of subreads) per template
record and a MAF alignment pairing each read with its template name. Templates
are named ``m{molecule_id:07d}``; this module maps simulated reads back to
their :class:`~muc_one_up.read_simulator.molecules.Molecule`, renames them to
unique ``{base}_h{hap}_m{molecule_id}`` identifiers and writes a gzipped TSV
truth manifest.
"""

from __future__ import annotations

import gzip
from collections.abc import Iterable, Iterator, Mapping
from pathlib import Path
from typing import IO

from ..exceptions import ReadSimulationError
from .molecules import Molecule

TRUTH_COLUMNS = (
    "read_id",
    "hap",
    "molecule",
    "kind",
    "strand",
    "src_start",
    "src_end",
    "n_hp_edits",
    "hp_edits",
    "detail",
)


def template_id(molecule_id: int) -> str:
    """Template record name for a molecule id."""
    return f"m{molecule_id:07d}"


def molecule_index(template: str) -> int:
    """Molecule id encoded in a template record name."""
    return int(template.lstrip("m"))


def _open_text(path: Path) -> IO[str]:
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path)


def parse_maf_read_templates(maf_path: Path) -> dict[str, str]:
    """Map read name -> template name from a pbsim3 MAF (template line, then read line)."""
    mapping: dict[str, str] = {}
    block: list[str] = []
    with _open_text(Path(maf_path)) as handle:
        for line in [*handle, "\n"]:
            if line.startswith("s "):
                block.append(line.split()[1])
            elif block:
                if len(block) != 2:
                    raise ReadSimulationError(
                        f"Unexpected MAF block in {maf_path}: expected template and read, got {block}"
                    )
                mapping[block[1]] = block[0]
                block = []
    return mapping


def _read_fastq(path: Path) -> Iterator[tuple[str, str, str]]:
    with _open_text(path) as handle:
        while header := handle.readline():
            seq, _, qual = handle.readline(), handle.readline(), handle.readline()
            yield header[1:].split()[0], seq.rstrip("\n"), qual.rstrip("\n")


def _resolve_template(read_name: str, mapping: Mapping[str, str], zmw: Mapping[str, str]) -> str:
    if read_name in mapping:
        return mapping[read_name]
    movie_zmw = read_name.rsplit("/", 1)[0]  # CCS reads: '{prefix}/{zmw}/ccs'
    if movie_zmw in zmw:
        return zmw[movie_zmw]
    raise ReadSimulationError(f"Simulated read '{read_name}' has no template in the pbsim3 MAF")


def _truth_row(read_id: str, molecule: Molecule) -> list[str]:
    edits = ",".join(f"{e.pos}:{e.base}:{e.true_len}>{e.new_len}" for e in molecule.hp_edits)
    return [
        read_id,
        str(molecule.hap),
        str(molecule.id),
        molecule.kind,
        molecule.strand,
        str(molecule.src_start),
        str(molecule.src_end),
        str(len(molecule.hp_edits)),
        edits,
        molecule.detail,
    ]


def relabel_reads(
    fastqs: Iterable[Path],
    read_templates: Mapping[str, str],
    molecules: Mapping[int, Molecule],
    base: str,
    out_fastq: Path,
    truth_tsv: Path,
) -> int:
    """Rename simulated reads to unique truth-bearing ids and write the truth manifest.

    Returns the number of reads written. Raises if any read cannot be traced to a
    molecule, so truth is never silently incomplete.
    """
    zmw = {name.rsplit("/", 1)[0]: tmpl for name, tmpl in read_templates.items() if "/" in name}
    count = 0
    with open(out_fastq, "w") as fq_out, gzip.open(truth_tsv, "wt") as truth_out:
        truth_out.write("\t".join(TRUTH_COLUMNS) + "\n")
        for fastq in fastqs:
            for name, seq, qual in _read_fastq(Path(fastq)):
                molecule = molecules[molecule_index(_resolve_template(name, read_templates, zmw))]
                read_id = f"{base}_h{molecule.hap}_{template_id(molecule.id)}"
                fq_out.write(f"@{read_id}\n{seq}\n+\n{qual}\n")
                truth_out.write("\t".join(_truth_row(read_id, molecule)) + "\n")
                count += 1
    return count
