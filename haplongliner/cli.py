import argparse
import sys
from .module1_RM import run_module1
from .module2_SV import run_module2
from .module3_DB import run_module3
from .utils import check_dependencies

__version__ = "0.1.0"

HS1_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz"
HG38_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"

def main():
    parser = argparse.ArgumentParser(
        prog="haplongliner",
        description=f"HapLongLINEr v{__version__}: Full-length L1 discovery pipeline",
        add_help=False
    )
    parser.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS,
                        help="Show this help message and exit.")
    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {__version__}",
                        help="Show program's version number and exit.")

    subparsers = parser.add_subparsers(
        dest="command",
        title="Pipeline module to run",
    )

    # Module 1: RepeatMasker-based
    parser_rm = subparsers.add_parser("rm", help="Module 1: RepeatMasker-based L1 discovery", add_help=False)
    parser_rm.add_argument("-i", "--in", dest="input", required=True, help="Input haploid assembly FASTA")
    parser_rm.add_argument("-m", "--mask", required=True, help="RepeatMasker BED or .out file")

    parser_rm.add_argument(
        "--log-skipped",
        dest="log_skipped",
        help="File to log malformed RepeatMasker lines",
    )

    ref_group = parser_rm.add_mutually_exclusive_group(required=True)
    ref_group.add_argument("-r", "--reference", choices=["hs1", "hg38"], help="Reference genome: 'hs1' or 'hg38' (remote)")
    ref_group.add_argument("-c", "--custom", help="Custom reference FASTA or gzipped FASTA (local path)")

    parser_rm.add_argument("-o", "--out", dest="output", required=True, help="Output directory for intermediate files")
    parser_rm.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS,
                           help="Show this help message and exit.")

    # Module 2: SV-based
    parser_sv = subparsers.add_parser("sv", help="Module 2: SV-based L1 discovery", add_help=False)
    parser_sv.add_argument("-i", "--in", dest="input", required=True, help="Input haploid assembly FASTA")
    parser_sv.add_argument("-s", "--sv", required=True, help="Structural variant callset")
    parser_sv.add_argument("-l", "--l1ref", required=True, help="Pangenome-level L1 reference FASTA")
    parser_sv.add_argument("-o", "--out", dest="output", required=True, help="Output BED file")
    parser_sv.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS,
                           help="Show this help message and exit.")

    # Module 3: Database
    parser_db = subparsers.add_parser("db", help="Module 3: L1 sequence repository", add_help=False)
    parser_db.add_argument("-o", "--out", dest="output", required=True, help="Output directory or file")
    parser_db.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS,
                           help="Show this help message and exit.")

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    if args.command == "rm":
        # Determine reference path/URL
        if args.reference == "hs1":
            reference = HS1_URL
        elif args.reference == "hg38":
            reference = HG38_URL
        else:
            reference = args.custom
        run_module1(
            args.input,
            args.mask,
            reference,
            args.output,
            log_skipped=args.log_skipped,
        )
    elif args.command == "sv":
        run_module2(args.input, args.sv, args.l1ref, args.output)
    elif args.command == "db":
        run_module3(args.output)

if __name__ == "__main__":
    check_dependencies()
    main()