# -*- coding: utf-8 -*-
"""Create HPRC L1 database from BED and sequence archives.

This script parses the `HPRC_L1_hs1_master_v2*.bed` files along with the
corresponding sequence archives (`HPRC_L1_seq_by_site_v2*.zip`) and
produces a tab separated text file summarising all insertions.
Sequences are stored relative to the L1rp reference using orientation and
CIGAR strings to save space.
"""

from __future__ import annotations

import csv
import io
from pathlib import Path
from typing import Dict, Iterable, Tuple
import zipfile

from Bio import SeqIO

# Reuse the alignment helper from the package
from haplongliner.store_l1_diffs import _best_alignment


def _load_reference(path: Path) -> str:
    """Return sequence string from FASTA ``path``."""
    record = next(SeqIO.parse(str(path), "fasta"))
    return str(record.seq)


def _load_sequences(zip_path: Path) -> Dict[Tuple[str, str, str], str]:
    """Return mapping ``(l1_name, sample, hap) -> sequence`` from ``zip_path``."""
    seqs: Dict[Tuple[str, str, str], str] = {}
    with zipfile.ZipFile(zip_path) as zf:
        for name in zf.namelist():
            if not name.endswith(".fa") or name.startswith("__MACOSX/"):
                continue
            l1_name = Path(name).stem
            with zf.open(name) as fh:
                for record in SeqIO.parse(io.TextIOWrapper(fh), "fasta"):
                    parts = record.id.split("#")
                    if len(parts) < 2:
                        continue
                    sample = parts[0]
                    hap = parts[1]
                    seqs[(l1_name, sample, hap)] = str(record.seq)
    return seqs


def _parse_hs1_coord(field: str) -> Tuple[str, int, int, str]:
    """Return ``(chrom, start, end, strand)`` from ``chr1_123_456_+`` format."""
    chrom, start, end, strand = field.split("_")
    return chrom, int(start), int(end), strand


def _parse_bed(
    bed: Path,
    sequences: Dict[Tuple[str, str, str], str],
    ref_seq: str,
) -> Iterable[Dict[str, str]]:
    """Yield dictionaries for each BED row with alignment info added."""
    hs1_ranges: Dict[str, Tuple[str, int, int, str]] = {}
    rows = []

    with open(bed) as fh:
        for line in fh:
            if not line.strip():
                continue
            fields = line.rstrip().split()
            if len(fields) < 9:
                continue
            l1_name = fields[0]
            sample = fields[1]
            hap_status = fields[2]
            status = fields[3]
            lineage = fields[4]
            site_lineage = fields[5]
            present_freq = fields[6]
            intact_freq = fields[7]
            assembly_info = fields[8]
            hs1_coord = fields[9] if len(fields) > 9 else ""

            hap = "1" if "paternal" in hap_status else "2"
            seq = sequences.get((l1_name, sample, hap))
            orient = cigar = ""
            if seq:
                orient, cigar = _best_alignment(seq, ref_seq)

            if hs1_coord:
                chrom, start, end, strand = _parse_hs1_coord(hs1_coord)
                rng = hs1_ranges.get(l1_name)
                if rng is None:
                    hs1_ranges[l1_name] = (chrom, start, end, strand)
                else:
                    c, s, e, st = rng
                    hs1_ranges[l1_name] = (
                        c,
                        min(s, start),
                        max(e, end),
                        st,
                    )

            rows.append(
                {
                    "l1_name": l1_name,
                    "sample": sample,
                    "haplotype": hap_status,
                    "status": status,
                    "lineage": lineage,
                    "site_lineage": site_lineage,
                    "present_freq": present_freq,
                    "intact_freq": intact_freq,
                    "assembly_info": assembly_info,
                    "hs1_coord": hs1_coord,
                    "orientation": orient,
                    "cigar": cigar,
                }
            )

    # Update hs1 coordinate ranges to use largest span
    for row in rows:
        rng = hs1_ranges.get(row["l1_name"])
        if rng:
            chrom, start, end, strand = rng
            row["hs1_coord"] = f"{chrom}_{start}_{end}_{strand}"
        yield row


def create_database(output: Path) -> None:
    """Create the combined database and write it to ``output``."""
    ref_seq = _load_reference(Path("data") / "L1rp.fa")

    # Load all sequences
    seq_v2 = _load_sequences(Path("data") / "HPRC_L1_seq_by_site_v2.zip")
    seq_v2fl = _load_sequences(Path("data") / "HPRC_L1_seq_by_site_v2fl.zip")

    rows: Iterable[Dict[str, str]] = []
    rows = list(_parse_bed(Path("HPRC_L1_hs1_master_v2.bed"), seq_v2, ref_seq)) + \
        list(_parse_bed(Path("HPRC_L1_hs1_master_v2fl.bed"), seq_v2fl, ref_seq))

    fieldnames = [
        "l1_name",
        "sample",
        "haplotype",
        "status",
        "lineage",
        "site_lineage",
        "present_freq",
        "intact_freq",
        "assembly_info",
        "hs1_coord",
        "orientation",
        "cigar",
    ]

    with open(output, "w", newline="") as out_f:
        writer = csv.DictWriter(out_f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Create HPRC L1 database")
    parser.add_argument(
        "-o",
        "--output",
        default="hprc_l1_db.txt",
        help="Output TSV file",
    )
    args = parser.parse_args()
    create_database(Path(args.output))

