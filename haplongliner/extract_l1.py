"""Utility to extract full-length transposable element records from a BED file."""

import sys
from collections.abc import Iterable


def _expand_te_names(te_list: Iterable[str]) -> set[str]:
    """Return the set of RepeatMasker names for the requested TEs.

    ``te_list`` contains user supplied identifiers such as ``"L1"`` or
    ``"SVA"``.  The mapping is as follows::

        L1    -> {"L1HS", "L1PA2"}
        SVA   -> {"SVA_E", "SVA_F"}
        ALU   -> {"AluY"}
        HERVK -> {"HERVK-int", "LTR5_Hs"}

    Any other value is used as-is, allowing custom TE names.
    """

    mapping = {
        "L1": {"L1HS", "L1PA2"},
        "SVA": {"SVA_E", "SVA_F"},
        "ALU": {"AluY"},
        "HERVK": {"HERVK-int", "LTR5_Hs"},
    }
    expanded: set[str] = set()
    for te in te_list:
        te_up = te.upper()
        expanded.update(mapping.get(te_up, {te}))
    return expanded


def extract_te_from_bed(
    infile: str,
    outfile: str | None = None,
    min_length: int = 5000,
    te: str | Iterable[str] = ("L1", "L1PA3"),
) -> None:
    """Write BED records for selected TEs at least ``min_length`` bp long."""

    if isinstance(te, str):
        te_list = [t for t in te.split(",") if t]
    else:
        te_list = list(te)
    allowed = _expand_te_names(te_list)

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
            if name in allowed and (end - start) >= min_length:
                print(line, end="", file=out)
    if outfile:
        out.close()


# Backward compatible alias
def extract_l1_from_bed(
    infile: str,
    outfile: str | None = None,
    min_length: int = 5000,
    te: str | Iterable[str] = ("L1", "L1PA3"),
) -> None:
    extract_te_from_bed(infile, outfile, min_length, te)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Extract full-length TEs from BED file.",
    )
    parser.add_argument("infile", help="Input BED file")
    parser.add_argument("-o", "--outfile", help="Output file (default: stdout)")
    parser.add_argument(
        "-l",
        "--length",
        type=int,
        default=5000,
        help="Minimum length to keep (default: 5000)",
    )
    parser.add_argument(
        "-t",
        "--te",
        default="L1,L1PA3",
        help=(
            "Comma-separated list of TE families. Shortcuts: "
            "L1=L1HS,L1PA2; SVA=SVA_E,SVA_F; ALU=AluY; "
            "HERVK=HERVK-int,LTR5_Hs"
        ),
    )
    args = parser.parse_args()
    extract_te_from_bed(args.infile, args.outfile, args.length, args.te)

