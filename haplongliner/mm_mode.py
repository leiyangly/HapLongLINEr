import contextlib
import os
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

from pyfaidx import Fasta

from .combine_table import _read_intact, combine_table
from .extract_l1 import _expand_te_names
from .find_intact_orf import find_intact_orf
from .find_longest_orf import find_longest_orf
from .liftover_paf import liftover_paf
from .repeatmasker import read_collapsed_repeatmasker_out
from .rm_mode import (
    _combine_full_liftover,
    _extract_fasta,
    _fix_getorf_headers,
    _write_rm_sequences,
    download_if_needed,
)
from .utils import (
    _fix_blast_query_names,
    append_fasta,
    read_nonempty_fasta,
    run_quiet,
    sort_bed,
    verify_blast_db,
    verify_fasta_file,
)


L1RP_FASTA = Path("data") / "L1rp.fa"
MM_MERGE_GAP = 200


def _validate_mm_te(te: str) -> set[str]:
    """Return expanded TE names and ensure they are L1-only."""

    te_list = [t for t in te.split(",") if t]
    expanded = _expand_te_names(te_list)
    if not expanded or any(not t.upper().startswith("L1") for t in expanded):
        sys.exit("Error: mm mode only supports L1-family targets.")
    return expanded


def _mm_seed_min_span(min_length: int) -> int:
    """Return the minimum L1rp-aligned query span used for seeding."""

    if min_length <= 0:
        return 1000
    return max(1000, min_length - 1000)


def _collect_mm_candidates(
    paf_file: Path,
    out_bed: Path,
    *,
    min_span: int,
    merge_gap: int = MM_MERGE_GAP,
) -> int:
    """Collect and merge query intervals that align to L1rp."""

    by_scaffold: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
    with open(paf_file) as fh:
        for line in fh:
            if not line.strip():
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 12:
                continue
            qname = fields[0]
            try:
                qstart = int(fields[2])
                qend = int(fields[3])
            except ValueError:
                continue
            strand = fields[4]
            qspan = qend - qstart
            if qspan < min_span:
                continue
            by_scaffold[(qname, strand)].append((qstart, qend))

    count = 0
    with open(out_bed, "w") as out:
        for (qname, strand), intervals in sorted(by_scaffold.items()):
            intervals.sort()
            cur_start, cur_end = intervals[0]
            for start, end in intervals[1:]:
                if start <= cur_end + merge_gap:
                    cur_end = max(cur_end, end)
                    continue
                out.write(
                    f"{qname}\t{cur_start}\t{cur_end}\tL1rp_seed\t{cur_end - cur_start}\t{strand}\n"
                )
                count += 1
                cur_start, cur_end = start, end
            out.write(
                f"{qname}\t{cur_start}\t{cur_end}\tL1rp_seed\t{cur_end - cur_start}\t{strand}\n"
            )
            count += 1
    return count


def _shorten_fasta_headers(fasta_path: Path) -> tuple[Path, Path]:
    """Write a copy of ``fasta_path`` with short headers and record mapping."""

    short_fa = fasta_path.with_name(fasta_path.stem + "_short.fa")
    list_file = fasta_path.with_name(fasta_path.stem + ".list")
    with open(fasta_path) as src, open(short_fa, "w") as dst, open(list_file, "w") as lst:
        idx = 1
        for line in src:
            if line.startswith(">"):
                orig = line[1:].strip()
                name = f"seq{idx}"
                dst.write(f">{name}\n")
                lst.write(f"{name}\t{orig}\n")
                idx += 1
            else:
                dst.write(line)
    return short_fa, list_file


def _restore_rm_out_names(short_out: Path, list_file: Path, out_file: Path) -> None:
    """Restore original sequence names in RepeatMasker ``.out`` file."""

    mapping: dict[str, str] = {}
    with open(list_file) as fh:
        for line in fh:
            if not line.strip():
                continue
            short, full = line.rstrip().split("\t", 1)
            mapping[short] = full

    with open(short_out) as src, open(out_file, "w") as dst:
        for line in src:
            parts = line.split()
            if (
                len(parts) > 4
                and parts[0].replace(".", "", 1).isdigit()
                and parts[4] in mapping
            ):
                parts[4] = mapping[parts[4]]
                line = " ".join(parts) + "\n"
            dst.write(line)


def _run_candidate_repeatmasker(candidate_fa: Path) -> Path:
    """Run RepeatMasker on candidate FASTA and return restored ``.out`` path."""

    rm_out = candidate_fa.with_suffix(candidate_fa.suffix + ".out")
    if not candidate_fa.exists() or candidate_fa.stat().st_size == 0:
        rm_out.touch()
        return rm_out

    short_fa, list_file = _shorten_fasta_headers(candidate_fa)
    try:
        run_quiet(
            [
                "RepeatMasker",
                "-e",
                "rmblast",
                "-species",
                "human",
                short_fa.name,
            ],
            cwd=candidate_fa.parent,
            check=True,
        )
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            "RepeatMasker failed while annotating minimap2-seeded candidate sequences. "
            "Ensure RepeatMasker and the rmblast engine are installed and "
            "configured correctly."
        ) from exc

    short_out = candidate_fa.parent / (short_fa.name + ".out")
    _restore_rm_out_names(short_out, list_file, rm_out)
    return rm_out


def _project_mm_repeatmasker(
    out_file: Path,
    out_bed: Path,
    *,
    min_length: int,
    te: str,
) -> int:
    """Project candidate-level RepeatMasker hits back to assembly coordinates."""

    allowed = _validate_mm_te(te)
    collapsed = read_collapsed_repeatmasker_out(out_file)

    written = 0
    with open(out_bed, "w") as out:
        for hit in collapsed:
            if not hit.is_l1 or not any(hit.family.startswith(name) for name in allowed):
                continue
            header = hit.query_name.split(",")
            if len(header) < 4:
                continue
            chrom = header[0]
            try:
                seed_start = int(header[1])
                seed_end = int(header[2])
            except ValueError:
                continue
            seed_strand = header[3]
            if seed_strand not in {"+", "-"}:
                continue
            local_start = hit.query_start
            local_end = hit.query_end

            if seed_strand == "+":
                start = seed_start + local_start
                end = seed_start + local_end
            else:
                start = seed_end - local_end
                end = seed_end - local_start
            length = end - start
            if length < min_length:
                continue
            strand = seed_strand if hit.strand == "+" else ("-" if seed_strand == "+" else "+")
            out.write(f"{chrom}\t{start}\t{end}\t{hit.family}\t{length}\t{strand}\n")
            written += 1
    return written


def run_mm_mode(
    input_fasta,
    reference_fasta,
    output_dir="mm_mode_output",
    log=None,
    min_length: int = 0,
    asm: int = 10,
    liftover: str = "full",
    te: str = "L1,L1PA3",
):
    """Minimap2-seeded L1 discovery followed by candidate-only RepeatMasker."""

    if log is None:
        log = os.getenv("HAPLOGLINER_LOG")
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    verify_fasta_file(input_fasta)
    ref_is_url = reference_fasta.startswith("http://") or reference_fasta.startswith("https://")
    if not ref_is_url:
        verify_fasta_file(reference_fasta)
    verify_fasta_file(str(L1RP_FASTA))

    expanded_te = _validate_mm_te(te)
    perform_orf = min_length >= 5000 and bool(expanded_te & {"L1HS", "L1PA2", "L1PA3"})

    info_lines = [
        "[INFO] MM mode running with:",
        f"[INFO]   Input: {input_fasta}",
        f"[INFO]   L1 seed reference: {L1RP_FASTA}",
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
    print("[STEP 1] Discovering candidate L1 loci with minimap2 and candidate-level RepeatMasker")

    seed_paf = outdir / "genome_vs_l1rp.paf"
    with open(seed_paf, "w") as out:
        run_quiet(
            [
                "minimap2",
                "-x",
                f"asm{asm}",
                "--secondary=yes",
                "-N",
                "1000",
                str(L1RP_FASTA),
                str(input_fasta),
            ],
            check=True,
            stdout=out,
        )

    seed_bed = outdir / "seed.bed"
    seed_min_span = _mm_seed_min_span(min_length)
    seed_count = _collect_mm_candidates(seed_paf, seed_bed, min_span=seed_min_span)

    seed_fa = outdir / "seed.fa"
    _extract_fasta(Fasta(str(input_fasta)), seed_bed, seed_fa)
    seed_rm_out = _run_candidate_repeatmasker(seed_fa)

    candidate_bed = outdir / "cand.bed"
    candidate_count = _project_mm_repeatmasker(
        seed_rm_out,
        candidate_bed,
        min_length=min_length,
        te=te,
    )

    te_counts: Counter[str] = Counter()
    with open(candidate_bed) as fh:
        for line in fh:
            if not line.strip():
                continue
            te_counts[line.split()[3]] += 1
    num_candidates = sum(te_counts.values())
    print(
        f"[SUM] Seeded {seed_count} candidate loci from minimap2; {candidate_count} candidates remained after RepeatMasker lineage assignment across {len(te_counts)} TE types"
    )

    fa = Fasta(str(input_fasta))
    candidate_lift_fa = outdir / "cand_lift.fa"
    candidate_ins_fa = outdir / "cand_ins.fa"
    candidate_fa = outdir / "cand.fa"
    candidate_ins_fa.write_text("")
    if perform_orf:
        _extract_fasta(fa, candidate_bed, candidate_lift_fa)
        candidate_fa.write_text("")
        append_fasta(candidate_fa, candidate_lift_fa)
        append_fasta(candidate_fa, candidate_ins_fa)

    if liftover == "flank2kb":
        print("\n[STEP 2] Performing liftover based on 2kb flanking sequences")
        candidate_minus2kb_bed = outdir / "cand_minus2kb.bed"
        candidate_plus2kb_bed = outdir / "cand_plus2kb.bed"
        run_quiet(
            f"""awk 'BEGIN{{OFS=\"\t\"}} {{$3=$2; $2=$2-2000; print $0}}' {candidate_bed} > {candidate_minus2kb_bed}""",
            shell=True,
            check=True,
        )
        run_quiet(
            f"""awk 'BEGIN{{OFS=\"\t\"}} {{$2=$3; $3=$3+2000; print $0}}' {candidate_bed} > {candidate_plus2kb_bed}""",
            shell=True,
            check=True,
        )

        candidate_minus2kb_fa = outdir / "cand_minus2kb.fa"
        candidate_plus2kb_fa = outdir / "cand_plus2kb.fa"
        _extract_fasta(fa, candidate_minus2kb_bed, candidate_minus2kb_fa)
        _extract_fasta(fa, candidate_plus2kb_bed, candidate_plus2kb_fa)

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
        with open(lifted_bed, "w") as out, contextlib.redirect_stdout(out):
            liftover_paf(str(aln_paf), str(candidate_bed), 0, 0, 2.0, False)
        lifted_count = sum(1 for line in open(lifted_bed) if line.strip())
        print(f"[SUM] Lifted coordinates for {lifted_count} candidates")

    if perform_orf and candidate_fa.exists() and candidate_fa.stat().st_size > 0:
        print("\n[STEP 3] Detecting intact ORFs")
        print("[INFO] Using getorf -minsize 810 to retain ORFs ≥270 aa (≥80% of L1 ORF1)")
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
        intact_out.touch()

    intact_records = _read_intact(intact_out)
    print(f"[SUM] {len(intact_records)}/{num_candidates} candidates with intact ORFs")

    print("\n[STEP 4] Preparing output files")
    combined_out = outdir / "haplongliner_mm.bed"
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

    mm_fa = outdir / "haplongliner_mm.fa"
    _write_rm_sequences(Path(input_fasta), candidate_bed, mm_fa)

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
    print(
        f"[SUM] Final table contains {total_final} entries; {intact_final} intact across {len(final_te_counts)} TE types"
    )
    print(
        f"[INFO] MM mode completed. Detected {te} elements ≥{min_length} bp. "
        f"Results in {combined_out} and {mm_fa}"
    )
