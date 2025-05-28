import argparse
import shutil
import sys
from .module1_RM import run_module1
from .module2_SV import run_module2
from .module3_DB import run_module3
from .utils import check_dependencies

__version__ = "0.1.0"

def main():
    parser = argparse.ArgumentParser(
        prog="haplongliner",
        description=f"HapLongLINEr v{__version__}: Full-length L1 discovery pipeline",
        add_help=False
    )
    # Custom help and version options
    parser.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS,
                        help="Show this help message and exit.")
    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {__version__}",
                        help="Show program's version number and exit.")

    # Custom subparser title and description
    subparsers = parser.add_subparsers(
        dest="command",
        title="Pipeline module to run",
    )

    # Module 1: RepeatMasker-based
    parser_rm = subparsers.add_parser("rm", help="Module 1: RepeatMasker-based L1 discovery", add_help=False)
    parser_rm.add_argument("-i", "-in", dest="input", required=True, help="Input haploid assembly FASTA")
    parser_rm.add_argument("-r", "--repeatmasker", required=True, help="RepeatMasker BED or .out file")
    parser_rm.add_argument("-o", "-out", dest="output", required=True, help="Output BED file")
    parser_rm.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS,
                           help="Show this help message and exit.")

    # Module 2: SV-based
    parser_sv = subparsers.add_parser("sv", help="Module 2: SV-based L1 discovery", add_help=False)
    parser_sv.add_argument("-i", "-in", dest="input", required=True, help="Input haploid assembly FASTA")
    parser_sv.add_argument("-s", "--sv", required=True, help="Structural variant callset")
    parser_sv.add_argument("-l", "--l1ref", required=True, help="Pangenome-level L1 reference FASTA")
    parser_sv.add_argument("-o", "-out", dest="output", required=True, help="Output BED file")
    parser_sv.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS,
                           help="Show this help message and exit.")

    # Module 3: Database
    parser_db = subparsers.add_parser("db", help="Module 3: L1 sequence repository", add_help=False)
    parser_db.add_argument("-o", "-out", dest="output", required=True, help="Output directory or file")
    parser_db.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS,
                           help="Show this help message and exit.")

    args = parser.parse_args()

    # If no arguments are provided, print help and exit
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    if args.command == "rm":
        run_module1(args.input, args.repeatmasker, args.output)
    elif args.command == "sv":
        run_module2(args.input, args.sv, args.l1ref, args.output)
    elif args.command == "db":
        run_module3(args.output)

if __name__ == "__main__":
    check_dependencies()
    main()