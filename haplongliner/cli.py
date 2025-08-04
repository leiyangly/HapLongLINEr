import argparse
import sys
from textwrap import dedent


class _Tee:
    """Simple stdout splitter."""

    def __init__(self, *streams):
        self.streams = streams

    def write(self, data):
        for s in self.streams:
            s.write(data)

    def flush(self):
        for s in self.streams:
            s.flush()

from .module1_RM import run_module1
from .module2_SV import run_module2
from .module3_DB import run_module3
from .utils import check_dependencies

__version__ = "0.1.0"

HS1_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz"
HG38_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"

def main():
    # Top-level parser mimicking the seqtk style help output
    description = dedent(
        f"""haplongliner

        Usage:   haplongliner <command> <arguments>
        Version: {__version__}

        Command: rm        RepeatMasker-based TE discovery
                 sv        SV-based L1 discovery
                 db        L1 sequence repository
        """
    )

    parser = argparse.ArgumentParser(
        prog="haplongliner",
        usage=argparse.SUPPRESS,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=description,
        add_help=False,
    )

    parser.add_argument(
        "-h",
        "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit.",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
        help="Show program's version number and exit.",
    )
    parser.add_argument(
        "-g",
        "--log",
        dest="log_file",
        metavar="FILE",
        help="Write console output to FILE",
    )

    subparsers = parser.add_subparsers(dest="command")

    # Module 1: RepeatMasker-based
    parser_rm = subparsers.add_parser(
        "rm",
        add_help=False,
        usage=argparse.SUPPRESS,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Usage:   haplongliner rm [options]",
        help="Module 1: RepeatMasker-based TE discovery",
    )
    parser_rm._positionals.title = ""
    parser_rm._optionals.title = "Options"
    parser_rm.add_argument("-i", "--in", dest="input", required=True, help="Input haploid assembly FASTA")
    parser_rm.add_argument("-m", "--mask", required=True, help="RepeatMasker BED or .out file")

    parser_rm.add_argument(
        "-l",
        "--length",
        dest="length",
        type=int,
        default=0,
        help="Minimum TE length (default: 0)",
    )

    parser_rm.add_argument(
        "-t",
        "--te",
        dest="te",
        default="L1,L1PA3",
        help=(
            "Comma-separated TE families to search. Shortcuts: "
            "L1=L1HS,L1PA2; SVA=SVA_E,SVA_F; ALU=AluY; "
            "HERVK=HERVK-int,LTR5_Hs (default: L1,L1PA3)"
        ),
    )

    parser_rm.add_argument(
        "-g",
        "--log",
        dest="log",
        help="File to log malformed RepeatMasker lines",
    )
    parser_rm.add_argument(
        "-a",
        "--asm",
        dest="asm",
        type=int,
        choices=[5, 10, 20],
        default=10,
        help="minimap2 asm preset (5, 10, or 20; default: 10)",
    )

    parser_rm.add_argument(
        "-v",
        "--liftover",
        dest="liftover",
        choices=["fl", "2kb"],
        default="fl",
        help="Liftover mode: 'fl' (whole-genome, default) or '2kb' (2kb flanks)",
    )

    ref_group = parser_rm.add_mutually_exclusive_group(required=True)
    ref_group.add_argument(
        "-r",
        "--ref",
        dest="ref",
        choices=["hs1", "hg38"],
        help="Reference genome: 'hs1' or 'hg38' (downloaded to data/ if missing)"
    )
    ref_group.add_argument("-c", "--custom", help="Custom reference FASTA or gzipped FASTA (local path)")

    parser_rm.add_argument("-o", "--out", dest="output", required=True, help="Output directory for intermediate files")
    parser_rm.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS,
                           help="Show this help message and exit.")

    # Module 2: SV-based
    parser_sv = subparsers.add_parser(
        "sv",
        add_help=False,
        usage=argparse.SUPPRESS,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Usage:   haplongliner sv [options]",
        help="Module 2: SV-based L1 discovery",
    )
    parser_sv._positionals.title = ""
    parser_sv._optionals.title = "Options"
    parser_sv.add_argument("-i", "--in", dest="input", required=True, help="Input haploid assembly FASTA")
    parser_sv.add_argument("-s", "--sv", required=True, help="Structural variant callset")
    l1_group = parser_sv.add_mutually_exclusive_group(required=True)
    l1_group.add_argument(
        "-1",
        "--l1ref",
        dest="l1ref",
        choices=["hprc"],
        help="Use built-in HPRC L1 reference",
    )
    l1_group.add_argument(
        "-2",
        "--l1cus",
        dest="l1cus",
        help="Custom L1 BED file to generate flanks",
    )
    parser_sv.add_argument(
        "-r",
        "--ref",
        dest="ref",
        choices=["hs1", "hg38"],
        required=True,
        help="Reference genome for liftover flanks",
    )
    parser_sv.add_argument(
        "-l",
        "--length",
        dest="length",
        type=int,
        default=5000,
        help="Minimum L1 length (default: 5000)",
    )
    parser_sv.add_argument(
        "-g",
        "--log",
        dest="log",
        help="File to log malformed SV lines",
    )
    parser_sv.add_argument(
        "-a",
        "--asm",
        dest="asm",
        type=int,
        choices=[5, 10, 20],
        default=10,
        help="minimap2 asm preset (5, 10, or 20; default: 10)",
    )
    parser_sv.add_argument(
        "-o",
        "--out",
        dest="output",
        required=True,
        help="Output directory for intermediate files",
    )
    parser_sv.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS,
                           help="Show this help message and exit.")

    # Module 3: Database
    parser_db = subparsers.add_parser(
        "db",
        add_help=False,
        usage=argparse.SUPPRESS,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Usage:   haplongliner db [options]",
        help="Module 3: Retrieve L1 sequences by coordinate (hs1 only)",
    )
    parser_db._positionals.title = ""
    parser_db._optionals.title = "Options"
    parser_db.add_argument(
        "-q",
        "--query",
        dest="query",
        required=True,
        help="Coordinate on hs1 as chr:start-end or BED file",
    )
    parser_db.add_argument(
        "-e",
        "--extra",
        dest="extra",
        help="Additional FASTA with non-HPRC L1s",
    )
    parser_db.add_argument(
        "-o", "--out", dest="output", required=True, help="Output FASTA file"
    )
    parser_db.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS,
                           help="Show this help message and exit.")

    # Hide the autogenerated positional arguments section in the main help
    parser._action_groups = [g for g in parser._action_groups if g.title != "positional arguments"]

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    if len(sys.argv) == 2:
        if sys.argv[1] == "rm":
            parser_rm.print_help(sys.stderr)
            sys.exit(1)
        elif sys.argv[1] == "sv":
            parser_sv.print_help(sys.stderr)
            sys.exit(1)
        elif sys.argv[1] == "db":
            parser_db.print_help(sys.stderr)
            sys.exit(1)

    args = parser.parse_args()

    log_handle = None
    if getattr(args, "log_file", None):
        log_handle = open(args.log_file, "w")
        sys.stdout = _Tee(sys.stdout, log_handle)

    if args.command == "rm":
        check_dependencies()
        # Determine reference path/URL
        if args.ref == "hs1":
            reference = HS1_URL
        elif args.ref == "hg38":
            reference = HG38_URL
        else:
            reference = args.custom
        run_module1(
            args.input,
            args.mask,
            reference,
            args.output,
            log=args.log,
            min_length=args.length,
            asm=args.asm,
            liftover=args.liftover,
            te=args.te,
        )
    elif args.command == "sv":
        # Determine reference path/URL for generating flanks when needed
        if args.ref == "hs1":
            reference = HS1_URL
        else:
            reference = HG38_URL
        check_dependencies()
        run_module2(
            args.input,
            args.sv,
            reference,
            args.output,
            l1ref=args.l1ref,
            l1cus=args.l1cus,
            log=args.log,
            min_length=args.length,
            asm=args.asm,
        )
    elif args.command == "db":
        run_module3(args.query, args.output, extra_fasta=args.extra)

    if log_handle:
        log_handle.close()

if __name__ == "__main__":
    main()
