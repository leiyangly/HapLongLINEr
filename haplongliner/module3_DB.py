from __future__ import annotations

from pathlib import Path
from typing import DefaultDict, Dict, Iterable, List, Tuple
import re
import sys
import io
import zipfile

from pyfaidx import Fasta

from .utils import verify_bed_file, verify_fasta_file


def _verify_zip_file(path: Path) -> None:
    """Exit if ``path`` is missing or appears to be a Git LFS placeholder."""
    if not path.exists():
        sys.exit(f"Error: data file '{path}' not found.")
    if path.stat().st_size < 200:
        with open(path) as fh:
            first_line = fh.readline()
        if first_line.startswith("version https://git-lfs.github.com"):
            sys.exit(
                "Error: data files appear to be placeholders. Run 'git lfs pull' to fetch them."
            )


def _revcomp(seq: str) -> str:
    complement = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(complement)[::-1]


def _best_alignment(query: str, reference: str) -> Tuple[str, str]:
    """Return orientation (+/-) and CIGAR string for best alignment."""
    import edlib

    plus = edlib.align(query, reference, mode="NW", task="path")
    minus = edlib.align(_revcomp(query), reference, mode="NW", task="path")
    if minus["editDistance"] < plus["editDistance"]:
        return "-", minus["cigar"]
    return "+", plus["cigar"]


def _load_sequences(zip_path: Path) -> Dict[str, List[Tuple[str, str]]]:
    """Return mapping ``L1_name -> [(id, sequence), ...]`` from ``zip_path``."""
    _verify_zip_file(zip_path)
    seqs: Dict[str, List[Tuple[str, str]]] = {}
    with zipfile.ZipFile(zip_path) as zf:
        for name in zf.namelist():
            if not name.endswith(".fa") or name.startswith("__MACOSX/"):
                continue
            l1_name = Path(name).stem
            with zf.open(name) as fh:
                record_id = None
                lines: List[str] = []
                for line in io.TextIOWrapper(fh):
                    line = line.strip()
                    if not line:
                        continue
                    if line.startswith(">"):
                        if record_id is not None:
                            seqs.setdefault(l1_name, []).append(
                                (record_id, "".join(lines))
                            )
                        record_id = line[1:].split()[0]
                        lines = []
                    else:
                        lines.append(line)
                if record_id is not None:
                    seqs.setdefault(l1_name, []).append((record_id, "".join(lines)))
    return seqs


def _load_database() -> Dict[str, List[Tuple[str, str]]]:
    base = Path("data")
    seqs = {}
    for fname in ["HPRC_L1_seq_by_site_v2.zip", "HPRC_L1_seq_by_site_v2fl.zip"]:
        seqs_zip = _load_sequences(base / fname)
        for k, vals in seqs_zip.items():
            seqs.setdefault(k, []).extend(vals)
    return seqs


def _load_reference(path: Path) -> str:
    fa = Fasta(str(path))
    first = next(iter(fa.records.values()))
    return str(first[:])


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

    base = Path("data")
    bed_path = base / "HPRC_L1_hs1_v2_v2fl.bed"
    if not bed_path.exists():
        sys.exit("Error: hs1 reference BED not found in data directory")

    # Load reference coordinates
    by_chrom: DefaultDict[str, List[Tuple[int, int, str, str]]] = DefaultDict(list)
    with open(bed_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, name, _len, strand = line.strip().split()[:6]
            by_chrom[chrom].append((int(start), int(end), name, strand))

    # Load sequence database and reference sequence
    seq_db = _load_database()
    ref_seq = _load_reference(base / "L1rp.fa")

    with open(output, "w") as out_f:
        for chrom, qstart, qend in coords:
            hits = [e for e in by_chrom.get(chrom, []) if not (qend <= e[0] or qstart >= e[1])]
            if not hits:
                print(
                    f"No L1 in hs1 at {chrom}:{qstart}-{qend}. Use -e/--extra to add sequences."
                )
                continue
            for start, end, name, strand in hits:
                records = seq_db.get(name)
                if not records:
                    print(f"Warning: sequences for {name} not found in data archives")
                    continue
                for rid, seq in records:
                    orient, cigar = _best_alignment(seq, ref_seq)
                    header = f">{rid};{chrom}:{start}-{end};{orient};{cigar}"
                    out_f.write(header + "\n" + seq + "\n")

        if extra_fasta:
            with open(extra_fasta) as fh:
                for line in fh:
                    out_f.write(line)

