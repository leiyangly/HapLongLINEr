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

from .rm_mode import run_rm_mode
from .sv_mode import run_sv_mode
from .sequence_retrieval_function import run_sequence_retrieval_function
from .liftover_paf import liftover_paf as run_liftover_paf
from .utils import check_dependencies
from .extract_l1 import _expand_te_names

__version__ = "0.1.0"

HS1_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz"
HG38_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"

def main():
    # Top-level parser mimicking the seqtk style help output
    description = dedent(
        f"""haplongliner

        Usage:   haplongliner <command> <arguments>
        Version: {__version__}

        Command: rm        RM mode (RepeatMasker-based TE discovery)
                 sv        SV mode (SV-based L1 discovery)
                 seq       Sequence retrieval function
                 lf        Liftover PAF alignments
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

    # RM mode: RepeatMasker-based
    parser_rm = subparsers.add_parser(
        "rm",
        add_help=False,
        usage=argparse.SUPPRESS,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Usage:   haplongliner rm [options]",
        help="RM mode: RepeatMasker-based TE discovery",
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
        default=None,
        help="Minimum TE length (default: 5000 for L1-only, otherwise 0)",
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
        "-y",
        "--legacy",
        dest="legacy",
        choices=["yes", "no"],
        default="no",
        help="Use legacy 2kb flank liftover ('yes') or full-genome liftover ('no', default)",
    )

    parser_rm.add_argument(
        "-r",
        "--ref",
        dest="ref",
        default="hs1",
        help="Reference genome: hs1 (default), hg38, or path to FASTA",
    )

    parser_rm.add_argument("-o", "--out", dest="output", required=True, help="Output directory for intermediate files")
    parser_rm.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS,
                           help="Show this help message and exit.")

    # SV mode
    parser_sv = subparsers.add_parser(
        "sv",
        add_help=False,
        usage=argparse.SUPPRESS,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Usage:   haplongliner sv [options]",
        help="SV mode: SV-based L1 discovery",
    )
    parser_sv._positionals.title = ""
    parser_sv._optionals.title = "Options"
    parser_sv.add_argument("-i", "--in", dest="input", required=True, help="Input haploid assembly FASTA")
    parser_sv.add_argument("-s", "--sv", required=True, help="Structural variant callset")
    parser_sv.add_argument(
        "-f",
        "--teref",
        dest="teref",
        default="hprc",
        help="TE reference: hprc, hs1, hg38, or path to BED (default: hprc)",
    )
    parser_sv.add_argument(
        "-r",
        "--ref",
        dest="ref",
        default="hs1",
        help="Reference genome for liftover flanks: hs1 (default), hg38, or path to FASTA",
    )
    parser_sv.add_argument(
        "-l",
        "--length",
        dest="length",
        type=int,
        default=None,
        help="Minimum TE length (default: 5000 for L1-only, otherwise 0)",
    )
    parser_sv.add_argument(
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

    # Sequence retrieval function
    parser_seq = subparsers.add_parser(
        "seq",
        add_help=False,
        usage=argparse.SUPPRESS,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Usage:   haplongliner seq [options]",
        help="Sequence retrieval function: Retrieve L1 sequences by coordinate (hs1 only)",
    )
    parser_seq._positionals.title = ""
    parser_seq._optionals.title = "Options"
    parser_seq.add_argument(
        "-q",
        "--query",
        dest="query",
        required=True,
        help="Coordinate on hs1 as chr:start-end or BED file",
    )
    parser_seq.add_argument(
        "-e",
        "--extra",
        dest="extra",
        help="Additional FASTA with non-HPRC L1s",
    )
    parser_seq.add_argument(
        "-o", "--out", dest="output", required=True, help="Output FASTA file"
    )
    parser_seq.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS,
                           help="Show this help message and exit.")

    # Liftover function
    parser_lf = subparsers.add_parser(
        "lf",
        add_help=False,
        usage=argparse.SUPPRESS,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Usage:   haplongliner lf <aln.paf> <query.bed> [options]",
        help="Liftover PAF alignments",
    )
    parser_lf._positionals.title = ""
    parser_lf._optionals.title = "Options"
    parser_lf.add_argument("aln_paf", help="PAF alignment file or '-' for stdin")
    parser_lf.add_argument("query_bed", help="BED file on query sequences")
    parser_lf.add_argument("-m", action="store_true", help="merge BED intervals")
    parser_lf.add_argument("-q", type=int, default=0, help="min mapping quality [0]")
    parser_lf.add_argument("-l", type=int, default=0, help="min alignment length [0]")
    parser_lf.add_argument(
        "-d",
        type=float,
        default=2.0,
        help="max sequence divergence (>=1 to disable) [2.0]",
    )
    parser_lf.add_argument(
        "-h",
        "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit.",
    )

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
        elif sys.argv[1] == "seq":
            parser_seq.print_help(sys.stderr)
            sys.exit(1)
        elif sys.argv[1] == "lf":
            parser_lf.print_help(sys.stderr)
            sys.exit(1)

    args = parser.parse_args()

    log_handle = None
    if getattr(args, "log_file", None):
        log_handle = open(args.log_file, "w")
        sys.stdout = _Tee(sys.stdout, log_handle)

    if args.command == "rm":
        if args.length is None:
            te_list = [t for t in args.te.split(",") if t]
            expanded_te = _expand_te_names(te_list)
            if all(t.upper().startswith("L1") for t in expanded_te):
                args.length = 5000
            else:
                args.length = 0
        check_dependencies()
        # Determine reference path/URL
        if args.ref == "hs1":
            reference = HS1_URL
        elif args.ref == "hg38":
            reference = HG38_URL
        else:
            reference = args.ref
        liftover_mode = "flank2kb" if args.legacy == "yes" else "full"
        run_rm_mode(
            args.input,
            args.mask,
            reference,
            args.output,
            log=args.log,
            min_length=args.length,
            asm=args.asm,
            liftover=liftover_mode,
            te=args.te,
        )
    elif args.command == "sv":
        if args.length is None:
            te_list = [t for t in args.te.split(",") if t]
            expanded_te = _expand_te_names(te_list)
            if all(t.upper().startswith("L1") for t in expanded_te):
                args.length = 5000
            else:
                args.length = 0
        # Determine reference path/URL for generating flanks when needed
        if args.ref == "hs1":
            reference = HS1_URL
        elif args.ref == "hg38":
            reference = HG38_URL
        else:
            reference = args.ref
        check_dependencies()
        run_sv_mode(
            args.input,
            args.sv,
            reference,
            args.output,
            teref=args.teref,
            te=args.te,
            log=args.log,
            min_length=args.length,
            asm=args.asm,
        )
    elif args.command == "seq":
        run_sequence_retrieval_function(args.query, args.output, extra_fasta=args.extra)
    elif args.command == "lf":
        try:
            run_liftover_paf(args.aln_paf, args.query_bed, args.q, args.l, args.d, args.m)
        except BrokenPipeError:
            pass

    if log_handle:
        log_handle.close()

if __name__ == "__main__":
    main()
