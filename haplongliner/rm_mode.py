import subprocess
from pathlib import Path
import gzip
import re
import urllib.request
import shutil
import os

from pyfaidx import Fasta
from collections import Counter

from .find_longest_orf import find_longest_orf
from .find_intact_orf import find_intact_orf
from .combine_table import combine_table, _read_intact
from .extract_l1 import _expand_te_names
from .utils import (
    verify_blast_db,
    run_quiet,
    verify_fasta_file,
    verify_repeatmasker_file,
    _fix_blast_query_names,
    read_nonempty_fasta,
    sort_bed,
    append_fasta,
)
from .liftover_paf import liftover_paf


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
    print(f"[INFO] Skipped {len(skipped)} malformed lines")


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


def _resolve_repeatmasker_file(
    input_fasta: str,
    repeatmasker_file: str | None,
    outdir: Path,
) -> Path:
    """Return a validated RepeatMasker output, generating one when needed."""

    if repeatmasker_file:
        verify_repeatmasker_file(repeatmasker_file)
        return Path(repeatmasker_file)

    fasta_path = Path(input_fasta)
    if str(fasta_path).endswith(".gz"):
        raise RuntimeError(
            "Automatic RepeatMasker execution requires an uncompressed FASTA input. "
            "Provide --mask or decompress the FASTA first."
        )

    rm_out = outdir / f"{fasta_path.name}.out"
    if rm_out.exists():
        verify_repeatmasker_file(str(rm_out))
        print(f"[INFO] Reusing existing RepeatMasker output: {rm_out}")
        return rm_out

    print(f"[INFO] No RepeatMasker file provided; masking {fasta_path} first")
    try:
        run_quiet(
            [
                "RepeatMasker",
                "-e",
                "rmblast",
                "-species",
                "human",
                "-dir",
                str(outdir),
                str(fasta_path),
            ],
            check=True,
        )
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            "RepeatMasker failed while masking the input assembly. "
            "Ensure RepeatMasker and the rmblast engine are installed and "
            "configured correctly."
        ) from exc

    verify_repeatmasker_file(str(rm_out))
    print(f"[INFO] RepeatMasker output written to {rm_out}")
    return rm_out


def _fix_getorf_headers(fa_path: Path, orig_fa: Path | None = None) -> None:
    """Rewrite ``getorf`` FASTA headers for easier downstream parsing.

    ``getorf`` produces headers in the form ``>name_<idx> [<s> - <e>]`` and also
    replaces commas in ``name`` with underscores.  This function converts the
    headers to ``>name,<idx>,<s>,<e>`` so the different fields can be reliably
    retrieved by simply splitting on commas.  If ``orig_fa`` is provided, the
    original headers are used to restore commas that ``getorf`` replaced with
    underscores.
    """

    mapping: Dict[str, str] = {}
    if orig_fa:
        with open(orig_fa) as fh:
            for line in fh:
                if line.startswith(">"):
                    name = line[1:].strip()
                    mapping[name.replace(",", "_")] = name

    tmp = fa_path.with_suffix(".tmp")
    pattern = re.compile(r"^>([^\s]+)_([0-9]+)\s+\[([0-9]+)\s*-\s*([0-9]+)\]")
    with open(fa_path) as fin, open(tmp, "w") as out:
        for line in fin:
            if line.startswith(">"):
                m = pattern.match(line)
                if m:
                    name, idx, start, end = m.groups()
                    orig = mapping.get(name, name)
                    line = f">{orig},{idx},{start},{end}\n"
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


def _combine_full_liftover(
    lifted_bed: Path,
    intact_file: Path,
    cand_bed: Path,
    out_file: Path,
    *,
    min_length: int = 0,
) -> None:
    """Integrate ORF status and liftover information for full-mode liftover."""

    lifted: dict[tuple[str, int, int, str], tuple[str, str, str, str]] = {}
    with open(lifted_bed) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.strip().split()
            if len(fields) < 6:
                continue
            chr_ref, start_ref, end_ref, name, _score, out_strand = fields[:6]
            parts = name.split(",")
            if len(parts) < 5:
                continue
            qchrom, qstart, qend, qstrand, _ = parts[:5]
            try:
                lifted[(qchrom, int(qstart), int(qend), qstrand)] = (
                    chr_ref,
                    start_ref,
                    end_ref,
                    out_strand,
                )
            except ValueError:
                continue

    intact = _read_intact(intact_file)

    with open(cand_bed) as bed, open(out_file, "w") as out:
        for line in bed:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.strip().split()
            if len(f) < 6:
                continue
            chrom, start, end, name, length, strand = f[:6]
            start_i = int(start)
            end_i = int(end)
            key = (chrom, start_i, end_i, strand)
            ikey = f"{chrom}_{start_i}_{end_i}"
            status = "intact" if ikey in intact else "present"

            chr_ref, start_ref, end_ref, out_strand = lifted.get(
                key, ("NA", "NA", "NA", "NA")
            )

            ref_len = "NA"
            if start_ref != "NA" and end_ref != "NA":
                try:
                    rl = int(end_ref) - int(start_ref)
                    if min_length == 0 or rl <= 2 * min_length:
                        ref_len = str(rl)
                    else:
                        chr_ref, start_ref, end_ref, ref_len = "NA", "NA", "NA", "NA"
                except ValueError:
                    chr_ref, start_ref, end_ref, ref_len = "NA", "NA", "NA", "NA"

            scaffold_info = f"{chrom},{start},{end},{length},{strand},RPM"
            out.write(
                f"{chr_ref}\t{start_ref}\t{end_ref}\t{name}\t{ref_len}\t{out_strand}\t{status}\t{scaffold_info}\n"
            )

def run_rm_mode(
    input_fasta,
    reference_fasta,
    repeatmasker_file=None,
    output_dir="rm_mode_output",
    log=None,
    min_length: int = 0,
    asm: int = 10,
    liftover: str = "full",
    te: str = "L1,L1PA3",
):
    """
    RepeatMasker-based TE discovery pipeline.
    Downloads remote reference if needed.
    Handles RepeatMasker BED, BED.gz, .out, or .out.gz input. When
    ``repeatmasker_file`` is omitted, RepeatMasker is first run on
    ``input_fasta`` and the generated ``.out`` file is used automatically.
    ``log`` specifies a file to log malformed RepeatMasker lines. If not
    provided, the ``HAPLOGLINER_LOG`` environment variable is checked.
    ``min_length`` controls the minimum TE length to retain (default 0 bp).
    ``asm`` sets the minimap2 assembly preset (5, 10, or 20; default 10).
    ``liftover`` chooses between 'full' (whole-genome liftover) and
    'flank2kb' (2kb flanking sequence liftover). Default is 'full'.
    ``te`` is a comma-separated list of TE families to search for.  Shortcuts:
    ``L1``=L1HS/L1PA2, ``SVA``=SVA_E/SVA_F, ``ALU``=AluY, ``HERVK``=HERVK-int/LTR5_Hs.
    """
    if log is None:
        log = os.getenv("HAPLOGLINER_LOG")
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Validate input files
    verify_fasta_file(input_fasta)
    repeatmasker_path = _resolve_repeatmasker_file(input_fasta, repeatmasker_file, outdir)
    ref_is_url = reference_fasta.startswith("http://") or reference_fasta.startswith(
        "https://"
    )
    if not ref_is_url:
        verify_fasta_file(reference_fasta)

    te_list = [t for t in te.split(",") if t]
    expanded_te = _expand_te_names(te_list)
    perform_orf = (
        min_length >= 5000
        and bool(expanded_te & {"L1HS", "L1PA2", "L1PA3"})
    )

    info_lines = [
        "[INFO] RM mode running with:",
        f"[INFO]   Input: {input_fasta}",
        f"[INFO]   RepeatMasker: {repeatmasker_path}",
        f"[INFO]   Reference: {reference_fasta}",
        f"[INFO]   Output Dir: {outdir}",
        f"[INFO]   TE types: {te} (min length {min_length})",
        f"[INFO]   ASM preset: asm{asm}",
    ]
    print("\n".join(info_lines))

    if ref_is_url:
        data_dir = Path("data")
        data_dir.mkdir(exist_ok=True)
        ref_local = data_dir / Path(reference_fasta).name
        reference_fasta = download_if_needed(reference_fasta, ref_local)
        verify_fasta_file(reference_fasta)
    if not perform_orf:
        print(
            "[INFO] Skipping intact ORF detection and associated BLASTP/sequence extraction steps "
            "(requires --length ≥5000 and --te including L1HS/L1PA2/L1PA3)"
        )

    print()

    print("[STEP 1] Parsing input and extracting sequences")
    # 1-3. Parse RepeatMasker file to unified BED6 and obtain candidate sequences
    parsed_bed = outdir / "parsed_repeatmasker.bed"
    parse_repeatmasker(repeatmasker_path, parsed_bed, log)

    # Extract full-length TEs from parsed BED
    candidate_bed = outdir / "cand.bed"
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
            "-t",
            te,
        ],
        check=True,
    )

    # Extract the sequence of the full-length TEs using pyfaidx when ORF detection is enabled
    fa = Fasta(str(input_fasta))
    candidate_lift_fa = outdir / "cand_lift.fa"
    candidate_ins_fa = outdir / "cand_ins.fa"
    candidate_fa = outdir / "cand.fa"
    if perform_orf:
        _extract_fasta(fa, candidate_bed, candidate_lift_fa)
        candidate_ins_fa.write_text("")
        candidate_fa.write_text("")
        append_fasta(candidate_fa, candidate_lift_fa)
        append_fasta(candidate_fa, candidate_ins_fa)

    # Summary for Step 1
    total_entries = 0
    te_counts: Counter[str] = Counter()
    with open(parsed_bed) as fh:
        for line in fh:
            if not line.strip():
                continue
            total_entries += 1
    with open(candidate_bed) as fh:
        for line in fh:
            if not line.strip():
                continue
            te = line.split()[3]
            te_counts[te] += 1
    num_candidates = sum(te_counts.values())
    te_type_count = len(te_counts)
    print(
        f"[SUM] Parsed {total_entries} entries; {num_candidates} candidates after filtering across {te_type_count} TE types"
    )

    if liftover == "flank2kb":
        print("\n[STEP 2] Performing liftover based on 2kb flanking sequences")
        # 4-6. Extract flanking 2kb regions, obtain their sequences and map them to the reference genome
        # Extract flanking 2kb regions (upstream and downstream)
        candidate_minus2kb_bed = outdir / "cand_minus2kb.bed"
        candidate_plus2kb_bed = outdir / "cand_plus2kb.bed"
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
        candidate_minus2kb_fa = outdir / "cand_minus2kb.fa"
        candidate_plus2kb_fa = outdir / "cand_plus2kb.fa"
        _extract_fasta(fa, candidate_minus2kb_bed, candidate_minus2kb_fa)
        _extract_fasta(fa, candidate_plus2kb_bed, candidate_plus2kb_fa)

        # Map flanking regions to reference genome with minimap2 (using local FASTA)
        candidate_minus2kb_paf = outdir / "cand_minus2kb.paf"
        candidate_plus2kb_paf = outdir / "cand_plus2kb.paf"
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
        print(f"[SUM] Generated flanking mappings for {num_candidates} candidates")
    else:
        print("\n[STEP 2] Performing liftover using whole-genome alignment")
        aln_paf = outdir / "genome_alignment.paf"
        with open(aln_paf, "w") as out:
            run_quiet(
                [
                    "minimap2",
                    "-x",
                    f"asm{asm}",
                    "-c",
                    "--cs",
                    str(reference_fasta),
                    str(input_fasta),
                ],
                check=True,
                stdout=out,
            )
        lifted_bed = outdir / "lifted.bed"
        import contextlib
        with open(lifted_bed, "w") as out, contextlib.redirect_stdout(out):
            liftover_paf(str(aln_paf), str(candidate_bed), 0, 0, 2.0, False)
        lifted_count = sum(1 for line in open(lifted_bed) if line.strip())
        print(f"[SUM] Lifted coordinates for {lifted_count} candidates")

    if perform_orf:
        print("\n[STEP 3] Detecting intact ORFs")
        print("[INFO] Using getorf -minsize 810 to retain ORFs ≥270 aa (≥80% of L1 ORF1)")
        # 7-8. Detect ORFs and identify intact ones
        # Detect ORFs and choose the longest ORF1/ORF2 per locus
        orf_fa = outdir / "cand_orf.fa"
        run_quiet(
            [
                "getorf",
                "-sequence",
                str(candidate_fa),
                "-find",
                "1",
                "-outseq",
                str(orf_fa),
                "-minsize",
                "810",
            ],
            check=True,
        )
        _fix_getorf_headers(orf_fa, candidate_fa)
        blastp_out = outdir / "cand_orf.blastp"
        db_prefix = Path("data") / "L1rpORF12p.fa"
        verify_blast_db(db_prefix)
        filtered_fasta = read_nonempty_fasta(orf_fa)
        if filtered_fasta:
            run_quiet(
                [
                    "blastp",
                    "-db",
                    str(db_prefix),
                    "-query",
                    "-",
                    "-outfmt",
                    "6 std qlen slen sacc",
                    "-out",
                    str(blastp_out),
                ],
                input=filtered_fasta,
                text=True,
                check=True,
            )
            _fix_blast_query_names(blastp_out)
        else:
            blastp_out.touch()
        longest_orf_out = outdir / "cand_orf_combine.blastp"
        find_longest_orf(blastp_out, longest_orf_out)
        _fix_blast_query_names(longest_orf_out)

        # Identify intact ORFs
        intact_out = outdir / "cand_orf_intact.blastp"
        find_intact_orf(longest_orf_out, intact_out)
        _fix_blast_query_names(intact_out)
        append_fasta(orf_fa, candidate_fa)
    else:
        print(
            "\n[STEP 3] Skipping intact ORF detection and BLASTP steps "
            "(requires --length ≥5000 and --te including L1HS/L1PA2/L1PA3)"
        )
        intact_out = outdir / "cand_orf_intact.blastp"

    # Summary for Step 3
    intact_records = _read_intact(intact_out)
    print(f"[SUM] {len(intact_records)}/{num_candidates} candidates with intact ORFs")

    print("\n[STEP 4] Preparing output files")
    # 9-10. Integrate ORF status, liftover information and write candidate sequences
    # Integrate ORF status and liftover information
    combined_out = outdir / "haplongliner_rm.bed"
    if liftover == "flank2kb":
        combine_table(
            candidate_plus2kb_paf,
            candidate_minus2kb_paf,
            intact_out,
            candidate_bed,
            combined_out,
            min_length=min_length,
        )
    else:
        _combine_full_liftover(
            lifted_bed,
            intact_out,
            candidate_bed,
            combined_out,
            min_length=min_length,
        )

    sort_bed(combined_out)

    rm_fa = outdir / "haplongliner_rm.fa"
    _write_rm_sequences(Path(input_fasta), candidate_bed, rm_fa)

    # Summary for Step 4
    total_final = 0
    intact_final = 0
    final_te_counts: Counter[str] = Counter()
    with open(combined_out) as fh:
        for line in fh:
            if not line.strip():
                continue
            total_final += 1
            fields = line.strip().split()
            final_te_counts[fields[3]] += 1
            if len(fields) > 6 and fields[6] == "intact":
                intact_final += 1
    te_type_count = len(final_te_counts)
    print(
        f"[SUM] Final table contains {total_final} entries; {intact_final} intact across {te_type_count} TE types"
    )

    # Final output table
    print(
        f"[INFO] RM mode completed. Detected {te} elements ≥{min_length} bp. "
        f"Results in {combined_out} and {rm_fa}"
    )

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
