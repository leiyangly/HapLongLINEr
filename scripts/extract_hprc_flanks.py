#!/usr/bin/env python3
"""Extract flanking +/-2kb sequences for L1s listed in HPRC_L1_hs_v2_v2fl.bed.

This script mimics the flank extraction behaviour of Module 1 but operates on
``data/HPRC_L1_hs_v2_v2fl.bed``. It requires a reference FASTA from which the
flanking sequences will be retrieved.
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Tuple


def parse_bed(path: Path) -> List[Tuple[str, int, int, str, str]]:
    """Return list of entries from ``path`` as
    ``(chrom, start, end, name, strand)``.
    """
    entries = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split()[:5]
            if len(fields) < 5:
                continue
            chrom, start, end, name, strand = fields
            entries.append((chrom, int(start), int(end), name, strand))
    return entries


def load_fasta(path: Path) -> Dict[str, str]:
    """Return mapping ``sequence_name -> sequence`` from FASTA ``path``."""
    seqs: Dict[str, List[str]] = {}
    name = None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def write_fasta(entries: List[Tuple[str, str]], path: Path) -> None:
    """Write ``entries`` as FASTA to ``path``."""
    with open(path, "w") as out:
        for header, seq in entries:
            out.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i+60] + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract +/-2kb flanks for HPRC L1 insertions")
    parser.add_argument("reference", help="Reference FASTA containing chromosomes referenced by the BED file")
    parser.add_argument("output_prefix", help="Prefix for output FASTA files")
    args = parser.parse_args()

    bed_path = Path("data") / "HPRC_L1_hs_v2_v2fl.bed"
    entries = parse_bed(bed_path)
    sequences = load_fasta(Path(args.reference))

    minus_entries: List[Tuple[str, str]] = []
    plus_entries: List[Tuple[str, str]] = []

    for chrom, start, end, name, strand in entries:
        seq = sequences.get(chrom)
        if seq is None:
            continue
        up_start = max(0, start - 2000)
        up_end = start
        down_start = end
        down_end = end + 2000
        minus_seq = seq[up_start:up_end]
        plus_seq = seq[down_start:down_end]
        minus_entries.append((f"{name}_-2kb", minus_seq))
        plus_entries.append((f"{name}_+2kb", plus_seq))

    write_fasta(minus_entries, Path(f"{args.output_prefix}-2kb.fa"))
    write_fasta(plus_entries, Path(f"{args.output_prefix}+2kb.fa"))


if __name__ == "__main__":
    main()
