from __future__ import annotations

from pathlib import Path
from typing import List, Tuple
import re

from .utils import verify_bed_file, verify_fasta_file


Coord = Tuple[str, int, int]


def _parse_coords(query: str) -> List[Coord]:
    """Return a list of ``(chrom, start, end)`` tuples from ``query``.

    ``query`` may be a path to a BED-like file or a coordinate string in the
    form ``chr:start-end``. Multiple coordinates can be separated by commas.
    """
    path = Path(query)
    coords: List[Coord] = []
    if path.exists():
        verify_bed_file(str(path))
        with open(path) as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                chrom, start, end = line.strip().split()[:3]
                coords.append((chrom, int(start), int(end)))
        return coords

    for token in query.split(','):
        m = re.match(r"^([^:]+):(\d+)-(\d+)$", token.strip())
        if not m:
            raise ValueError(f"Invalid coordinate: {token}")
        chrom, start, end = m.groups()
        coords.append((chrom, int(start), int(end)))
    return coords


def run_module3(query: str, output: str, extra_fasta: str | None = None) -> None:
    """Extract sequences for L1 insertions at ``query``.

    ``query`` may be a BED file or ``chr:start-end`` string. ``output`` is the
    path for the resulting FASTA file. If ``extra_fasta`` is provided, those
    sequences will be appended to the extracted FASTA.
    """
    coords = _parse_coords(query)
    print(
        "Module 3 running with:\n"
        f"  Coordinates: {coords}\n"
        f"  Output: {output}"
    )
    if extra_fasta:
        verify_fasta_file(extra_fasta)
        print(f"  Extra FASTA: {extra_fasta}")
    # TODO: implement extraction logic

