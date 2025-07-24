import subprocess
from pathlib import Path
import gzip
import re
import urllib.request
import shutil
import os

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
            fout.write(
                f"{chrom}\t{start}\t{end}\t{name}\t{length}\t{strand}\n"
            )

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
    with urllib.request.urlopen(url) as response, open(local_path, 'wb') as out_file:
        shutil.copyfileobj(response, out_file)
    print(f"[INFO] Download complete: {local_path}")
    return str(local_path)


def _rename_fa_with_bed(fa_path: Path, bed_path: Path) -> None:
    """Rename FASTA headers based on BED coordinates."""
    records = []
    with open(bed_path) as bed:
        for line in bed:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.strip().split()
            if len(fields) < 5:
                continue
            if len(fields) == 5:
                chrom, start, end, _name, strand = fields
            else:
                chrom, start, end, *_rest, strand = fields[:6]
            records.append((chrom, int(start) + 1, int(end), strand))

    tmp = fa_path.with_suffix(".tmp")
    with open(fa_path) as fin, open(tmp, "w") as out:
        idx = 0
        for line in fin:
            if line.startswith(">"):
                chrom, s, e, strand = records[idx]
                out.write(f">{chrom};{s};{e};{strand}\n")
                idx += 1
            else:
                out.write(line)
    os.replace(tmp, fa_path)

def run_module1(
    input_fasta,
    repeatmasker_file,
    reference_fasta,
    output_dir="module1_output",
    log=None,
    min_length: int = 5000,
):
    """
    RepeatMasker-based L1 discovery pipeline.
    Downloads remote reference if needed.
    Handles RepeatMasker BED, BED.gz, .out, or .out.gz input.
    ``log`` specifies a file to log malformed RepeatMasker lines. If not
    provided, the ``HAPLOGLINER_LOG`` environment variable is checked.
    ``min_length`` controls the minimum L1 length to retain (default 5000 bp).
    """
    if log is None:
        log = os.getenv("HAPLOGLINER_LOG")
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Validate input files
    verify_fasta_file(input_fasta)
    verify_repeatmasker_file(repeatmasker_file)
    if not reference_fasta.startswith("http://") and not reference_fasta.startswith("https://"):
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
    )

    print()

    print("[STEP 1] Parsing RepeatMasker output")
    # 1. Parse RepeatMasker file to unified BED6
    parsed_bed = outdir / "parsed_repeatmasker.bed"
    parse_repeatmasker(repeatmasker_file, parsed_bed, log)

    print("\n[STEP 2] Extracting full-length L1s")
    # 2. Extract full-length L1s from parsed BED
    fl_bed = outdir / "FL.bed"
    run_quiet([
        "python3",
        "-m",
        "haplongliner.extract_l1",
        parsed_bed,
        "-o",
        str(fl_bed),
        "-l",
        str(min_length),
    ], check=True)

    print("\n[STEP 3] Extracting full-length L1 sequences")
    # 3. Extract the sequence of the full-length L1s (plus and minus strand)
    fl_fa = outdir / "FL.fa"
    with open(fl_fa, "w") as out_fa:
        # Plus strand
        plus_cmd = (
            f"awk '$6==\"+\"' {fl_bed} | "
            f"seqtk subseq {input_fasta} - | "
            f"seqtk seq -U -l 0 - | "
            "sed -e '/^>/ s/:/;/' -e '/^>/ s/-/;/' -e '/^>/ s/$/;+/'"
        )
        run_quiet(plus_cmd, shell=True, stdout=out_fa, check=True)
        # Minus strand
        minus_cmd = (
            f"awk '$6==\"-\"' {fl_bed} | "
            f"seqtk subseq {input_fasta} - | "
            f"seqtk seq -U -r -l 0 - | "
            "sed -e '/^>/ s/:/;/' -e '/^>/ s/-/;/' -e '/^>/ s/$/;-/'"
        )
        run_quiet(minus_cmd, shell=True, stdout=out_fa, check=True)

        
    print("\n[STEP 4] Extracting 2kb flanking regions")
    # 4. Extract flanking 2kb regions (upstream and downstream)
    fl_minus2kb_bed = outdir / "FL-2kb.bed"
    fl_plus2kb_bed = outdir / "FL+2kb.bed"
    # Upstream
    run_quiet(
        f"""awk 'BEGIN{{OFS=\"\t\"}} {{$3=$2; $2=$2-2000; print $0}}' {fl_bed} > {fl_minus2kb_bed}""",
        shell=True, check=True
    )
    # Downstream
    run_quiet(
        f"""awk 'BEGIN{{OFS=\"\t\"}} {{$2=$3; $3=$3+2000; print $0}}' {fl_bed} > {fl_plus2kb_bed}""",
        shell=True, check=True
    )

    print("\n[STEP 5] Getting sequences for flanking regions")
    # 5. Extract sequences for flanking regions
    fl_minus2kb_fa = outdir / "FL-2kb.fa"
    fl_plus2kb_fa = outdir / "FL+2kb.fa"
    run_quiet(
        f"seqtk subseq {input_fasta} {fl_minus2kb_bed} | seqtk seq -U -l 0 - > {fl_minus2kb_fa}",
        shell=True,
        check=True,
    )
    _rename_fa_with_bed(fl_minus2kb_fa, fl_minus2kb_bed)
    run_quiet(
        f"seqtk subseq {input_fasta} {fl_plus2kb_bed} | seqtk seq -U -l 0 - > {fl_plus2kb_fa}",
        shell=True,
        check=True,
    )
    _rename_fa_with_bed(fl_plus2kb_fa, fl_plus2kb_bed)

    print("\n[STEP 6] Mapping flanks to reference genome")
    # 6. Map flanking regions to reference genome with minimap2 (using local FASTA)
    fl_minus2kb_paf = outdir / "FL-2kb.paf"
    fl_plus2kb_paf = outdir / "FL+2kb.paf"
    run_quiet(
        f"minimap2 -x asm5 {reference_fasta} {fl_minus2kb_fa} > {fl_minus2kb_paf}",
        shell=True,
        check=True,
    )
    run_quiet(
        f"minimap2 -x asm5 {reference_fasta} {fl_plus2kb_fa} > {fl_plus2kb_paf}",
        shell=True,
        check=True,
    )

    print("\n[STEP 7] Detecting ORFs")
    # 7. Detect ORFs and choose the longest ORF1/ORF2 per locus
    orf_fa = outdir / "FLAllORF.fa"
    run_quiet([
        "getorf",
        "-sequence",
        str(fl_fa),
        "-find",
        "1",
        "-outseq",
        str(orf_fa),
    ], check=True)
    orf_bed = outdir / "FLAllORF.bed"
    process_orf_fasta(orf_fa, orf_bed)
    blastp_out = outdir / "FLAllORF.blastp"
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
    longest_orf_out = outdir / "FLAllORF.combine.blastp"
    find_longest_orf(blastp_out, longest_orf_out)

    print("\n[STEP 8] Identifying intact ORFs")
    # 8. Identify intact ORFs
    intact_out = outdir / "FLAllORF.intact.blastp"
    find_intact_orf(longest_orf_out, intact_out)

    print("\n[STEP 9] Integrating ORF status and liftover info")
    # 9. Integrate ORF status and liftover information
    combined_out = outdir / "HapLongLINErRM.txt"
    combine_table(
        fl_plus2kb_paf,
        fl_minus2kb_paf,
        intact_out,
        fl_bed,
        combined_out,
        min_length=min_length,
    )

    # Final output table
    print(f"Module 1 completed. Results in {combined_out}")

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

