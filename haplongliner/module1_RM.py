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

def _ensure_minimap_index(reference_fasta: str) -> str:
    """Return a minimap2 index for the reference, creating it if needed."""
    ref_path = Path(reference_fasta)
    # Determine index path (hs1.fa.gz -> hs1.fa.gz.mmi)
    if ref_path.suffix == ".gz":
        index_path = ref_path.with_suffix(ref_path.suffix + ".mmi")
    else:
        index_path = ref_path.with_suffix(".mmi")
    if not index_path.exists():
        print(f"[INFO] Building minimap2 index for {ref_path} ...")
        subprocess.run([
            "minimap2",
            "-x",
            "asm5",
            "-d",
            str(index_path),
            str(reference_fasta),
        ], check=True)
    else:
        print(f"[INFO] Using existing minimap2 index at {index_path}.")
    return str(index_path)


def run_module1(
    input_fasta,
    repeatmasker_file,
    reference_fasta,
    output_dir="module1_output",
    log_skipped=None,
):
    """
    RepeatMasker-based L1 discovery pipeline.
    Downloads remote reference if needed.
    Handles RepeatMasker BED, BED.gz, .out, or .out.gz input.
    ``log_skipped`` specifies a file to log malformed RepeatMasker lines. If not
    provided, the ``HAPLOGLINER_LOG_SKIPPED`` environment variable is checked.
    """
    if log_skipped is None:
        log_skipped = os.getenv("HAPLOGLINER_LOG_SKIPPED")
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

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

    print("[STEP 1] Parsing RepeatMasker output")
    # 1. Parse RepeatMasker file to unified BED6
    parsed_bed = outdir / "parsed_repeatmasker.bed"
    parse_repeatmasker(repeatmasker_file, parsed_bed, log_skipped)

    print("[STEP 2] Extracting full-length L1s")
    # 2. Extract full-length L1s (>=5000bp) from parsed BED
    fl_bed = outdir / "FL.bed"
    subprocess.run([
        "python3", "-m", "haplongliner.extract_l1", parsed_bed, "-o", str(fl_bed)
    ], check=True)

    print("[STEP 3] Extracting full-length L1 sequences")
    # 3. Extract the sequence of the full-length L1s (plus and minus strand)
    fl_fa = outdir / "FL.fa"
    with open(fl_fa, "w") as out_fa:
        # Plus strand
        plus_cmd = (
            f"awk '$6==\"+\"' {fl_bed} | "
            f"seqtk subseq {input_fasta} - | "
            f"seqtk seq -U -l 0 - | "
            "sed '/^>/ s/$/(+)/'"
        )
        subprocess.run(plus_cmd, shell=True, stdout=out_fa, check=True)
        # Minus strand
        minus_cmd = (
            f"awk '$6==\"-\"' {fl_bed} | "
            f"seqtk subseq {input_fasta} - | "
            f"seqtk seq -U -r -l 0 - | "
            "sed '/^>/ s/$/(-)/'"
        )
        subprocess.run(minus_cmd, shell=True, stdout=out_fa, check=True)

    # Sanitize FASTA headers for getorf compatibility
    fl_rename_fa = outdir / "FL.rename.fa"
    with open(fl_fa) as fin, open(fl_rename_fa, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                header = line.strip()
                # Temporarily replace strand annotation to protect the minus sign
                header = (
                    header.replace("(+)", "_plus")
                    .replace("(-)", "_minus")
                    .replace(":", "_")
                    .replace("-", "_")
                    .replace("_plus", "_+")
                    .replace("_minus", "_-")
                )
                fout.write(header + "\n")
            else:
                fout.write(line)

    print("[STEP 4] Extracting 2kb flanking regions")
    # 4. Extract flanking 2kb regions (upstream and downstream)
    fl_minus2kb_bed = outdir / "FL-2kb.bed"
    fl_plus2kb_bed = outdir / "FL+2kb.bed"
    # Upstream
    subprocess.run(
        f"""awk 'BEGIN{{OFS=\"\t\"}} {{$3=$2; $2=$2-2000; print $0}}' {fl_bed} > {fl_minus2kb_bed}""",
        shell=True, check=True
    )
    # Downstream
    subprocess.run(
        f"""awk 'BEGIN{{OFS=\"\t\"}} {{$2=$3; $3=$3+2000; print $0}}' {fl_bed} > {fl_plus2kb_bed}""",
        shell=True, check=True
    )

    print("[STEP 5] Getting sequences for flanking regions")
    # 5. Extract sequences for flanking regions
    fl_minus2kb_fa = outdir / "FL-2kb.fa"
    fl_plus2kb_fa = outdir / "FL+2kb.fa"
    subprocess.run(f"seqtk subseq {input_fasta} {fl_minus2kb_bed} | seqtk seq -U -l 0 - > {fl_minus2kb_fa}", shell=True, check=True)
    subprocess.run(f"seqtk subseq {input_fasta} {fl_plus2kb_bed} | seqtk seq -U -l 0 - > {fl_plus2kb_fa}", shell=True, check=True)

    print("[STEP 6] Mapping flanks to reference genome")
    # 6. Map flanking regions to reference genome with minimap2 (using local FASTA)
    fl_minus2kb_minimap = outdir / "FL-2kb.minimap.txt"
    fl_plus2kb_minimap = outdir / "FL+2kb.minimap.txt"
    index_fa = _ensure_minimap_index(reference_fasta)
    subprocess.run(
        f"minimap2 -x asm5 {index_fa} {fl_minus2kb_fa} > {fl_minus2kb_minimap}",
        shell=True,
        check=True,
    )
    subprocess.run(
        f"minimap2 -x asm5 {index_fa} {fl_plus2kb_fa} > {fl_plus2kb_minimap}",
        shell=True,
        check=True,
    )

    print("[STEP 7] Detecting ORFs")
    # 7. Detect ORFs and choose the longest ORF1/ORF2 per locus
    orf_fa = outdir / "FLAllORF.fa"
    subprocess.run([
        "getorf",
        "-sequence",
        fl_rename_fa,
        "-find",
        "1",
        "-outseq",
        str(orf_fa),
    ], check=True)
    orf_bed = outdir / "FLAllORF.bed"
    process_orf_fasta(orf_fa, orf_bed)
    blastp_out = outdir / "FLAllORF.blastp"
    subprocess.run(
        [
            "blastp",
            "-db",
            str(Path("data") / "L1rpORF12p.fa"),
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

    print("[STEP 8] Identifying intact ORFs")
    # 8. Identify intact ORFs
    intact_out = outdir / "FLAllORF.intact.blastp"
    find_intact_orf(longest_orf_out, intact_out)

    print("[STEP 9] Integrating ORF status and liftover info")
    # 9. Integrate ORF status and liftover information
    combined_out = outdir / "HapLongLINErRM.txt"
    combine_table(
        fl_plus2kb_minimap,
        fl_minus2kb_minimap,
        intact_out,
        fl_bed,
        combined_out,
    )

    # Final output table
    print(f"Module 1 completed. Results in {combined_out}")

    # Remove large intermediate files to save space
    for tmp in [
        blastp_out,
        orf_fa,
        fl_rename_fa,
        fl_fa,
        parsed_bed,
        fl_minus2kb_fa,
        fl_plus2kb_fa,
    ]:
        try:
            os.remove(tmp)
        except FileNotFoundError:
            pass
