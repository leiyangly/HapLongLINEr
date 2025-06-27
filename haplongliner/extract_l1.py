"""Utility to extract full-length L1 records from a BED file."""

import sys


def extract_l1_from_bed(infile: str, outfile: str | None = None, min_length: int = 5000) -> None:
    """Write BED records for L1s at least ``min_length`` bp long."""
    out = open(outfile, "w") if outfile else sys.stdout
    with open(infile) as f:
        for line in f:
            if line.strip() == "" or line.startswith("#"):
                continue
            fields = line.strip().split()
            if len(fields) < 4:
                continue
            name = fields[3]
            start = int(fields[1])
            end = int(fields[2])
            if name.startswith("L1") and (end - start) >= min_length:
                print(line, end="", file=out)
    if outfile:
        out.close()


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Extract full-length L1s from BED file."
    )
    parser.add_argument("infile", help="Input BED file")
    parser.add_argument("-o", "--outfile", help="Output file (default: stdout)")
    parser.add_argument(
        "-l",
        "--length",
        type=int,
        default=5000,
        help="Minimum L1 length to keep (default: 5000)",
    )
    args = parser.parse_args()
    extract_l1_from_bed(args.infile, args.outfile, args.length)
