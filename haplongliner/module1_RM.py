import subprocess
from pathlib import Path
import gzip
import re
import urllib.request
import shutil
import os

from pyfaidx import Fasta

from .process_orf import process_orf_fasta
from .find_longest_orf import find_longest_orf
from .find_intact_orf import find_intact_orf
from .combine_table import combine_table
from .utils import (
    verify_blast_db,
    run_quiet,
    verify_fasta_file,
    verify_repeatmasker_file,
)


def _revcomp(seq: str) -> str:
    complement = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(complement)[::-1]


def parse_repeatmasker(input_path, output_path, log_path=None):
    """
    Parse RepeatMasker BED, BED.gz, .out, or .out.gz file and write a unified
    BED-like file using 0-based half-open coordinates::

        chrom  start  end  name  length  strand

    ``log_path`` optionally records skipped malformed lines.
    """
    # Open plain or gzipped file
    opener = gzip.open if str(input_path).endswith(".gz") else open
    skipped = []
    with opener(input_path, "rt") as fin, open(output_path, "w") as fout:
        lines = fin.readlines()
        # Detect .out header (skip first 4 lines if header detected)
        if any("SW" in l and "perc" in l for l in lines[:4]):
            lines = lines[4:]

        for line in lines:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = re.split(r"\s+", line.strip())

            is_out = (
                len(fields) >= 14
                and fields[0].replace(".", "", 1).isdigit()
                and fields[5].isdigit()
                and fields[6].isdigit()
            )

            try:
                if is_out:
                    chrom = fields[4]
                    # RepeatMasker .out uses 1-based inclusive coordinates
                    # Convert to 0-based half-open
                    start = int(fields[5]) - 1
                    end = int(fields[6])
                    name = fields[9]
                    strand = fields[8]
                    strand = "-" if strand == "C" else "+"
                elif len(fields) >= 5:
                    chrom = fields[0]
                    start = int(fields[1])
                    end = int(fields[2])
                    name = fields[3]
                    strand = fields[5] if len(fields) >= 6 else fields[4]
                else:
                    raise ValueError
            except Exception:
                skipped.append(line.rstrip())
                continue

            length = end - start
            fout.write(f"{chrom}\t{start}\t{end}\t{name}\t{length}\t{strand}\n")

    if log_path and skipped:
        with open(log_path, "w") as logf:
            logf.write("\n".join(skipped) + "\n")
    print(f"Skipped {len(skipped)} malformed lines")


def download_if_needed(url, local_path):
    """
    Download the file from url to local_path if it does not exist.
    """
    local_path = Path(local_path)
    if local_path.exists():
        print(f"[INFO] Reference genome already exists at {local_path}.")
        return str(local_path)
    print(f"[INFO] Downloading reference genome from {url} ...")
    local_path.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url) as response, open(local_path, "wb") as out_file:
        shutil.copyfileobj(response, out_file)
    print(f"[INFO] Download complete: {local_path}")
    return str(local_path)


def _fix_getorf_headers(fa_path: Path) -> None:
    """Rewrite ``getorf`` FASTA headers for easier downstream parsing.

    ``getorf`` produces headers in the form ``>name_<idx> [<s> - <e>]``.  This
    function converts them to ``>name,<idx>,<s>,<e>`` so the different fields can
    be reliably retrieved by simply splitting on commas.  ``name`` itself may
    contain underscores (e.g. scaffold identifiers), so only the separators
    introduced by ``getorf`` are replaced.
    """

    tmp = fa_path.with_suffix(".tmp")
    pattern = re.compile(r"^>([^\s]+)_([0-9]+)\s+\[([0-9]+)\s*-\s*([0-9]+)\]")
    with open(fa_path) as fin, open(tmp, "w") as out:
        for line in fin:
            if line.startswith(">"):
                m = pattern.match(line)
                if m:
                    name, idx, start, end = m.groups()
                    line = f">{name},{idx},{start},{end}\n"
            out.write(line)
    os.replace(tmp, fa_path)


def _extract_fasta(fa: Fasta, bed: Path, out: Path) -> None:
    """Write sequences for regions in ``bed`` to ``out`` using ``fa``."""

    with open(bed) as b, open(out, "w") as out_f:
        for line in b:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.strip().split()
            if len(fields) < 6:
                continue
            chrom, start, end, *_rest, strand = fields[:6]
            start_i = int(start)
            end_i = int(end)
            seq = fa[chrom][start_i:end_i].seq.upper()
            if strand == "-":
                seq = _revcomp(seq)
            out_f.write(f">{chrom},{start_i},{end_i},{strand}\n{seq}\n")


def _write_rm_sequences(fasta: Path, bed: Path, out_fa: Path) -> None:
    """Write sequences from ``fasta`` for each BED entry in ``bed``."""

    fa = Fasta(str(fasta))
    with open(bed) as b, open(out_fa, "w") as out:
        for line in b:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.strip().split()
            if len(fields) < 6:
                continue
            chrom, start, end, *_rest, strand = fields[:6]
            start_i = int(start)
            end_i = int(end)
            seq = fa[chrom][start_i:end_i].seq.upper()
            if strand == "-":
                seq = _revcomp(seq)
            header = f"{chrom},{start_i},{end_i},{end_i - start_i},{strand},RPM"
            out.write(f">{header}\n{seq}\n")


def run_module1(
    input_fasta,
    repeatmasker_file,
    reference_fasta,
    output_dir="module1_output",
    log=None,
    min_length: int = 5000,
    asm: int = 10,
):
    """
    RepeatMasker-based L1 discovery pipeline.
    Downloads remote reference if needed.
    Handles RepeatMasker BED, BED.gz, .out, or .out.gz input.
    ``log`` specifies a file to log malformed RepeatMasker lines. If not
    provided, the ``HAPLOGLINER_LOG`` environment variable is checked.
    ``min_length`` controls the minimum L1 length to retain (default 5000 bp).
    ``asm`` sets the minimap2 assembly preset (5, 10, or 20; default 10).
    """
    if log is None:
        log = os.getenv("HAPLOGLINER_LOG")
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Validate input files
    verify_fasta_file(input_fasta)
    verify_repeatmasker_file(repeatmasker_file)
    if not reference_fasta.startswith("http://") and not reference_fasta.startswith(
        "https://"
    ):
        verify_fasta_file(reference_fasta)

    # If reference_fasta is a URL, download it to the data folder
    if reference_fasta.startswith("http://") or reference_fasta.startswith("https://"):
        data_dir = Path("data")
        data_dir.mkdir(exist_ok=True)
        ref_local = data_dir / Path(reference_fasta).name
        reference_fasta = download_if_needed(reference_fasta, ref_local)

    print(
        "Module 1 running with:\n"
        f"  Input: {input_fasta}\n"
        f"  RepeatMasker: {repeatmasker_file}\n"
        f"  Reference: {reference_fasta}\n"
        f"  Output Dir: {outdir}\n"
        f"  ASM preset: asm{asm}"
    )

    print()

    print("[STEP 1] Parsing input and extracting sequences")
    # 1-3. Parse RepeatMasker file to unified BED6 and obtain candidate sequences
    parsed_bed = outdir / "parsed_repeatmasker.bed"
    parse_repeatmasker(repeatmasker_file, parsed_bed, log)

    # Extract full-length L1s from parsed BED
    candidate_bed = outdir / "candidate.bed"
    run_quiet(
        [
            "python3",
            "-m",
            "haplongliner.extract_l1",
            parsed_bed,
            "-o",
            str(candidate_bed),
            "-l",
            str(min_length),
        ],
        check=True,
    )

    # Extract the sequence of the full-length L1s using pyfaidx
    candidate_fa = outdir / "candidate.fa"
    fa = Fasta(str(input_fasta))
    _extract_fasta(fa, candidate_bed, candidate_fa)

    print("\n[STEP 2] Performing liftover based on 2kb flanking sequences")
    # 4-6. Extract flanking 2kb regions, obtain their sequences and map them to the reference genome
    # Extract flanking 2kb regions (upstream and downstream)
    candidate_minus2kb_bed = outdir / "candidate_minus2kb.bed"
    candidate_plus2kb_bed = outdir / "candidate_plus2kb.bed"
    # Upstream
    run_quiet(
        f"""awk 'BEGIN{{OFS=\"\t\"}} {{$3=$2; $2=$2-2000; print $0}}' {candidate_bed} > {candidate_minus2kb_bed}""",
        shell=True,
        check=True,
    )
    # Downstream
    run_quiet(
        f"""awk 'BEGIN{{OFS=\"\t\"}} {{$2=$3; $3=$3+2000; print $0}}' {candidate_bed} > {candidate_plus2kb_bed}""",
        shell=True,
        check=True,
    )

    # Extract sequences for flanking regions
    candidate_minus2kb_fa = outdir / "candidate_minus2kb.fa"
    candidate_plus2kb_fa = outdir / "candidate_plus2kb.fa"
    _extract_fasta(fa, candidate_minus2kb_bed, candidate_minus2kb_fa)
    _extract_fasta(fa, candidate_plus2kb_bed, candidate_plus2kb_fa)

    # Map flanking regions to reference genome with minimap2 (using local FASTA)
    candidate_minus2kb_paf = outdir / "candidate_minus2kb.paf"
    candidate_plus2kb_paf = outdir / "candidate_plus2kb.paf"
    run_quiet(
        f"minimap2 -x asm{asm} {reference_fasta} {candidate_minus2kb_fa} > {candidate_minus2kb_paf}",
        shell=True,
        check=True,
    )
    run_quiet(
        f"minimap2 -x asm{asm} {reference_fasta} {candidate_plus2kb_fa} > {candidate_plus2kb_paf}",
        shell=True,
        check=True,
    )

    print("\n[STEP 3] Detecting intact ORFs")
    # 7-8. Detect ORFs and identify intact ones
    # Detect ORFs and choose the longest ORF1/ORF2 per locus
    orf_fa = outdir / "candidate_orf.fa"
    run_quiet(
        [
            "getorf",
            "-sequence",
            str(candidate_fa),
            "-find",
            "1",
            "-outseq",
            str(orf_fa),
        ],
        check=True,
    )
    _fix_getorf_headers(orf_fa)
    orf_bed = outdir / "candidate_orf.bed"
    process_orf_fasta(orf_fa, orf_bed)
    blastp_out = outdir / "candidate_orf.blastp"
    db_prefix = Path("data") / "L1rpORF12p.fa"
    verify_blast_db(db_prefix)
    run_quiet(
        [
            "blastp",
            "-db",
            str(db_prefix),
            "-query",
            str(orf_fa),
            "-outfmt",
            "6 std qlen slen sacc",
            "-out",
            str(blastp_out),
        ],
        check=True,
    )
    longest_orf_out = outdir / "candidate_orf.combine.blastp"
    find_longest_orf(blastp_out, longest_orf_out)

    # Identify intact ORFs
    intact_out = outdir / "candidate_orf.intact.blastp"
    find_intact_orf(longest_orf_out, intact_out)

    print("\n[STEP 4] Preparing output files")
    # 9-10. Integrate ORF status, liftover information and write candidate sequences
    # Integrate ORF status and liftover information
    combined_out = outdir / "haplongliner_rm.bed"
    combine_table(
        candidate_plus2kb_paf,
        candidate_minus2kb_paf,
        intact_out,
        candidate_bed,
        combined_out,
        min_length=min_length,
    )

    rm_fa = outdir / "haplongliner_rm.fa"
    _write_rm_sequences(Path(input_fasta), candidate_bed, rm_fa)

    # Final output table
    print(f"Module 1 completed. Results in {combined_out} and {rm_fa}")

    # Remove large intermediate files to save space
    # for tmp in [
    #     blastp_out,
    #     orf_fa,
    #     fl_fa,
    #     parsed_bed,
    #     fl_minus2kb_fa,
    #     fl_plus2kb_fa,
    # ]:
    #     try:
    #         os.remove(tmp)
    #     except FileNotFoundError:
    #         pass
