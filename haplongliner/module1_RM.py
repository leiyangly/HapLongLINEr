import subprocess
from pathlib import Path
import gzip
import re
import urllib.request
import shutil

from .process_orf import process_orf_fasta
from .find_longest_orf import find_longest_orf

def parse_repeatmasker(input_path, output_path):
    """
    Parse RepeatMasker BED, BED.gz, .out, or .out.gz file and write a unified BED-like file:
    chrom  start  end  name  .  strand
    """
    # Open plain or gzipped file
    opener = gzip.open if str(input_path).endswith(".gz") else open
    with opener(input_path, "rt") as fin, open(output_path, "w") as fout:
        lines = fin.readlines()
        # Detect .out header (skip first 4 lines if header detected)
        if any("SW" in l and "perc" in l for l in lines[:4]):
            lines = lines[4:]

        for line in lines:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = re.split(r'\s+', line.strip())
            # Try .out format
            if len(fields) >= 10 and fields[4] and fields[5].isdigit() and fields[6].isdigit():
                chrom = fields[4]
                start = int(fields[5]) - 1  # .out is 1-based, BED is 0-based
                end = int(fields[6])
                name = fields[9]
                strand = fields[8]
                strand = "-" if strand == "C" else "+"
            # Otherwise, treat as BED
            elif len(fields) >= 5:
                chrom = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                name = fields[3]
                strand = fields[4]
            else:
                continue
            fout.write(f"{chrom}\t{start}\t{end}\t{name}\t.\t{strand}\n")

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

def run_module1(input_fasta, repeatmasker_file, reference_fasta, output_bed="module1_output.bed"):
    """
    RepeatMasker-based L1 discovery pipeline.
    Downloads remote reference if needed.
    Handles RepeatMasker BED, BED.gz, .out, or .out.gz input.
    """
    output_bed = Path(output_bed)
    outdir = output_bed.parent if output_bed.parent != Path("") else Path(".")
    outdir.mkdir(exist_ok=True)

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
        f"  Output BED: {output_bed}\n"
    )

    # 1. Parse RepeatMasker file to unified BED6
    parsed_bed = outdir / "parsed_repeatmasker.bed"
    parse_repeatmasker(repeatmasker_file, parsed_bed)

    # 2. Extract full-length L1s (>=5000bp) from parsed BED
    fl_bed = outdir / "FL.bed"
    subprocess.run([
        "python3", "-m", "haplongliner.extract_l1", parsed_bed, "-o", str(fl_bed)
    ], check=True)

    # 3. Extract the sequence of the full-length L1s (plus and minus strand)
    fl_fa = outdir / "FL.fa"
    with open(fl_fa, "w") as out_fa:
        # Plus strand
        plus_cmd = (
            f"awk '$6==\"+\"' {fl_bed} | "
            f"seqtk subseq {input_fasta} - | "
            f"seqtk seq -U -l 0 -"
        )
        subprocess.run(plus_cmd, shell=True, stdout=out_fa)
        # Minus strand
        minus_cmd = (
            f"awk '$6==\"-\"' {fl_bed} | "
            f"seqtk subseq {input_fasta} - | "
            f"seqtk seq -U -r -l 0 -"
        )
        subprocess.run(minus_cmd, shell=True, stdout=out_fa)

    # 4. Extract flanking 2kb regions (upstream and downstream)
    fl_minus2kb_bed = outdir / "FL-2kb.bed"
    fl_plus2kb_bed = outdir / "FL+2kb.bed"
    # Upstream
    subprocess.run(
        f"""awk 'BEGIN{{OFS="\\t"}} {{$2=$2-2000; $3=$2+2000; print $0}}' {fl_bed} > {fl_minus2kb_bed}""",
        shell=True, check=True
    )
    # Downstream
    subprocess.run(
        f"""awk 'BEGIN{{OFS="\\t"}} {{$3=$3+2000; print $0}}' {fl_bed} > {fl_plus2kb_bed}""",
        shell=True, check=True
    )

    # 5. Extract sequences for flanking regions
    fl_minus2kb_fa = outdir / "FL-2kb.fa"
    fl_plus2kb_fa = outdir / "FL+2kb.fa"
    subprocess.run(f"seqtk subseq {input_fasta} {fl_minus2kb_bed} | seqtk seq -U -l 0 - > {fl_minus2kb_fa}", shell=True, check=True)
    subprocess.run(f"seqtk subseq {input_fasta} {fl_plus2kb_bed} | seqtk seq -U -l 0 - > {fl_plus2kb_fa}", shell=True, check=True)

    # 6. Map flanking regions to reference genome with minimap2 (using local FASTA)
    fl_minus2kb_minimap = outdir / "FL-2kb.minimap.txt"
    fl_plus2kb_minimap = outdir / "FL+2kb.minimap.txt"
    subprocess.run(
        f"minimap2 -x asm5 {reference_fasta} {fl_minus2kb_fa} > {fl_minus2kb_minimap}",
        shell=True,
        check=True,
    )
    subprocess.run(
        f"minimap2 -x asm5 {reference_fasta} {fl_plus2kb_fa} > {fl_plus2kb_minimap}",
        shell=True,
        check=True,
    )

    # 7. Detect ORFs and choose the longest ORF1/ORF2 per locus
    orf_fa = outdir / "FLAllORF.fa"
    subprocess.run(["getorf", "-sequence", fl_fa, "-find", "1", "-outseq", str(orf_fa)], check=True)
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

    # Final output BED
    shutil.copy(fl_bed, output_bed)
    print(f"Module 1 completed. Results in {output_bed}")
