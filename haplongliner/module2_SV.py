import re
from pathlib import Path
from typing import Dict, List, Tuple, Iterable, Set

from .utils import (
    run_quiet,
    verify_fasta_file,
    verify_sv_file,
    verify_bed_file,
)
from .module1_RM import download_if_needed


def _read_paf(path: Path) -> Dict[str, List[str]]:
    """Return mapping ``query_name -> fields`` from a minimap2 PAF."""
    hits: Dict[str, List[str]] = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            fields = line.rstrip().split('\t')
            qname = fields[0]
            if qname not in hits:
                hits[qname] = fields
    return hits


def _load_master_coords(paths: Iterable[Path]) -> Dict[str, Tuple[str, int, int, str]]:
    """Return mapping ``name -> (chrom, start, end, strand)`` from master BED files."""
    coords: Dict[str, Tuple[str, int, int, str]] = {}
    for path in paths:
        with open(path) as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.strip().split()
                if len(fields) < 5:
                    continue
                chrom, start, end, name, strand = fields[:5]
                try:
                    coords[name] = (chrom, int(start), int(end), strand)
                except ValueError:
                    continue
    return coords


def _liftover_l1s(
    minus_paf: Path,
    plus_paf: Path,
    ref_bed: Path,
    min_length: int,
    master_files: Iterable[Path] | None = None,
) -> List[Tuple[str, int, int, str, int, str, str, str]]:
    """Infer reference coordinates for each L1 by lifting from scaffolds.

    Returns tuples of ``(chrom, start, end, name, length, strand, plus, minus)``
    where ``plus`` and ``minus`` are ``scaffold;start;end;strand`` strings
    describing the originating PAF alignments.
    """
    if master_files is None:
        master_files = [Path("data") / "HPRC_L1_hs1_v2_v2fl.bed"]
    ref_coords = _load_master_coords(master_files)

    minus = _read_paf(minus_paf)
    plus = _read_paf(plus_paf)

    lifted: List[Tuple[str, int, int, str, int, str, str, str]] = []
    for name, (rchrom, rstart, rend, rstrand) in ref_coords.items():
        m = minus.get(f"{name}_-2kb")
        p = plus.get(f"{name}_+2kb")
        if not m or not p:
            continue
        if m[5] != p[5]:
            continue

        start_t: int | None = None
        end_t: int | None = None
        if int(m[7]) <= int(p[8]):
            start_t = int(m[8])
            end_t = int(p[7])
        if int(p[7]) <= int(m[8]):
            start_t = int(p[8])
            end_t = int(m[7])
        if start_t is None or end_t is None:
            continue
        if end_t < start_t:
            start_t, end_t = end_t, start_t

        length = rend - rstart
        if length >= min_length:
            plus_info = f"{p[5]};{p[7]};{p[8]};{p[4]}"
            minus_info = f"{m[5]};{m[7]};{m[8]};{m[4]}"
            lifted.append((rchrom, rstart, rend, name, length, rstrand, plus_info, minus_info))
    return lifted


def _write_bed(entries: List[Tuple[str, int, int, str, int, str, str, str]], path: Path) -> None:
    with open(path, 'w') as out:
        for chrom, start, end, name, length, strand, plus_info, minus_info in entries:
            out.write(
                f"{chrom}\t{start}\t{end}\t{name}\t{length}\t{strand}\t{plus_info}\t{minus_info}\n"
            )


def _parse_sv(
    sv_path: Path, log_path: Path | None = None
) -> Tuple[List[Tuple[str, int, int]], List[Tuple[str, int, int, str]]]:
    """Parse a simple VCF or BED SV file and return deletion and insertion regions.

    For insertions the inserted sequence is also returned if present. Unrecognized
    lines are written to ``log_path`` if provided."""
    deletions: List[Tuple[str, int, int]] = []
    insertions: List[Tuple[str, int, int, str]] = []
    is_vcf = sv_path.suffix.lower().endswith('vcf') or sv_path.suffix.lower() == '.gz'
    import gzip

    skipped: List[str] = []
    opener = gzip.open if str(sv_path).endswith(".gz") else open
    with opener(sv_path, "rt") as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip().split('\t')
            if is_vcf or len(fields) > 5:
                try:
                    chrom = fields[0]
                    pos = int(fields[1]) - 1
                    info = fields[7] if len(fields) > 7 else ''
                    svtype = None
                    end = None
                    for item in info.split(';'):
                        if item.startswith('SVTYPE='):
                            svtype = item.split('=', 1)[1]
                        elif item.startswith('END='):
                            end = int(item.split('=', 1)[1]) - 1
                    if svtype == 'DEL' and end is not None:
                        deletions.append((chrom, pos, end))
                    elif svtype == 'INS':
                        alt = fields[4]
                        seq = '' if alt.startswith('<') else alt
                        insertions.append((chrom, pos, pos + 1, seq))
                except Exception:
                    skipped.append(line.rstrip())
            else:
                try:
                    chrom, start, end = fields[:3]
                    svtype = fields[3].upper() if len(fields) > 3 else ''
                    if svtype == 'DEL':
                        deletions.append((chrom, int(start), int(end)))
                    elif svtype == 'INS':
                        seq = fields[4] if len(fields) > 4 else ''
                        insertions.append((chrom, int(start), int(end), seq))
                except Exception:
                    skipped.append(line.rstrip())
    if log_path and skipped:
        with open(log_path, 'w') as logf:
            logf.write('\n'.join(skipped) + '\n')
    return deletions, insertions


def _bedtools_intersect(a: Path, b: Path, output: Path) -> None:
    with open(output, 'w') as out:
        run_quiet(['bedtools', 'intersect', '-wa', '-wb', '-a', str(a), '-b', str(b)], check=True, stdout=out)


def _classify_sv(
    lifted: List[Tuple[str, int, int, str, int, str, str, str]],
    deletions: List[Tuple[str, int, int]],
    insertions: List[Tuple[str, int, int, str]],
    outdir: Path,
) -> Tuple[Dict[str, str], Dict[str, str]]:
    """Write SV BED files and intersect with lifted coordinates.

    Returns a status dictionary and a mapping of insertion sequences by L1 name."""

    lift_bed = outdir / "lifted.bed"
    del_bed = outdir / "sv_del.bed"
    ins_bed = outdir / "sv_ins.bed"

    _write_bed(lifted, lift_bed)

    with open(del_bed, "w") as out:
        for chrom, start, end in deletions:
            out.write(f"{chrom}\t{start}\t{end}\n")

    with open(ins_bed, "w") as out:
        for chrom, start, end, seq in insertions:
            out.write(f"{chrom}\t{start}\t{end}\t{seq}\n")

    inter_del = outdir / "intersect_del.bed"
    inter_ins = outdir / "intersect_ins.bed"
    _bedtools_intersect(lift_bed, del_bed, inter_del)
    _bedtools_intersect(lift_bed, ins_bed, inter_ins)

    status: Dict[str, str] = {name: "present" for _, _, _, name, _, _, _, _ in lifted}
    ins_seqs: Dict[str, str] = {}
    with open(inter_del) as fh:
        for line in fh:
            f = line.strip().split("\t")
            if len(f) < 9:
                continue
            a_start = int(f[1])
            a_end = int(f[2])
            name = f[3]
            d_start = int(f[6])
            d_end = int(f[7])
            overlap = max(0, min(a_end, d_end) - max(a_start, d_start))
            del_len = d_end - d_start
            cov = overlap / del_len if del_len else 0
            if cov >= 0.95:
                status[name] = "missing"
            else:
                status[name] = "absent"

    with open(inter_ins) as fh:
        for line in fh:
            f = line.strip().split("\t")
            if len(f) >= 10:
                name = f[3]
                seq = f[9]
                ins_seqs[name] = seq

    return status, ins_seqs


def _extract_sequences(
    fasta: Path,
    lifted: List[Tuple[str, int, int, str, int, str, str, str]],
    status: Dict[str, str],
    ins_seqs: Dict[str, str],
    out_fa: Path,
    min_length: int,
) -> None:
    bed_path = out_fa.with_suffix('.bed')
    with open(bed_path, 'w') as bed:
        for chrom, start, end, name, length, _, _, _ in lifted:
            if status.get(name) == 'missing' and abs((end - start) - length) / length < 0.1 and length >= min_length:
                bed.write(f"{chrom}\t{start}\t{end}\t{name}\n")
    tmp_del = out_fa.with_suffix('.del.fa')
    if bed_path.stat().st_size > 0:
        cmd = f"seqtk subseq {fasta} {bed_path} | seqtk seq -U -l 0 - > {tmp_del}"
        run_quiet(cmd, shell=True, check=True)
    else:
        tmp_del.touch()

    with open(out_fa, 'w') as out:
        if tmp_del.stat().st_size > 0:
            out.write(Path(tmp_del).read_text())
        for chrom, start, end, name, length, _, _, _ in lifted:
            seq = ins_seqs.get(name)
            if seq and len(seq) >= min_length:
                out.write(f">{name}\n{seq}\n")

    tmp_del.unlink(missing_ok=True)
    bed_path.unlink(missing_ok=True)


def _parse_repeatmasker(out_file: Path, l1_len: int, cov_thresh: float) -> Set[str]:
    """Return names of sequences covering ``cov_thresh`` of the L1 reference."""
    coverage: Dict[str, int] = {}
    with open(out_file) as fh:
        for line in fh:
            if line.startswith(' '):
                parts = line.split()
                if len(parts) >= 14 and re.search(r'L1', parts[9]):
                    name = parts[4]
                    rep_start = int(parts[11].strip('()'))
                    rep_end = int(parts[12].strip('()'))
                    cov = abs(rep_end - rep_start) + 1
                    coverage[name] = coverage.get(name, 0) + cov
    return {n for n, c in coverage.items() if c / l1_len >= cov_thresh}


def run_module2(
    input_fasta: str,
    sv_file: str,
    reference_fasta: str,
    output_dir: str,
    *,
    l1ref: str | None = None,
    l1cus: str | None = None,
    log: str | None = None,
    min_length: int = 5000,
) -> None:
    """RepeatMasker-free L1 discovery using structural variants.

    ``log`` specifies a file to record malformed SV lines if provided.
    ``reference_fasta`` can be a local path or URL (downloaded if needed).
    """

    print(
        "Module 2 running with:\n"
        f"  Input: {input_fasta}\n"
        f"  SV: {sv_file}\n"
        f"  Reference: {reference_fasta}\n"
        f"  Output Dir: {output_dir}\n"
        f"  Min Length: {min_length}"
    )

    print()

    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    verify_fasta_file(input_fasta)
    verify_sv_file(sv_file)

    print("[STEP 1] Preparing reference L1 flanking sequences")

    if l1cus:
        verify_bed_file(l1cus)
        ref_path = reference_fasta
        if ref_path.startswith("http://") or ref_path.startswith("https://"):
            data_dir = Path("data")
            data_dir.mkdir(exist_ok=True)
            ref_local = data_dir / Path(ref_path).name
            ref_path = download_if_needed(ref_path, ref_local)
        else:
            verify_fasta_file(ref_path)

        minus_bed = outdir / "custom_minus2kb.bed"
        plus_bed = outdir / "custom_plus2kb.bed"
        run_quiet(
            f"awk 'BEGIN{{OFS=\"\t\"}} {{$3=$2; $2=$2-2000; print $0}}' {l1cus} > {minus_bed}",
            shell=True,
            check=True,
        )
        run_quiet(
            f"awk 'BEGIN{{OFS=\"\t\"}} {{$2=$3; $3=$3+2000; print $0}}' {l1cus} > {plus_bed}",
            shell=True,
            check=True,
        )
        minus_fa = outdir / "-2kb.fa"
        plus_fa = outdir / "+2kb.fa"
        run_quiet(
            f"seqtk subseq {ref_path} {minus_bed} | seqtk seq -U -l 0 - > {minus_fa}",
            shell=True,
            check=True,
        )
        run_quiet(
            f"seqtk subseq {ref_path} {plus_bed} | seqtk seq -U -l 0 - > {plus_fa}",
            shell=True,
            check=True,
        )
        ref_bed = Path(l1cus)
    else:
        if l1ref != "hprc":
            raise ValueError("Unsupported L1 reference")
        minus_fa = Path("data") / "-2kb.fa"
        plus_fa = Path("data") / "+2kb.fa"
        ref_bed = Path("data") / "HPRC_L1_hs1_v2_v2fl.bed"

    print("\n[STEP 2] Mapping L1 flanks to input assembly")

    minus_paf = outdir / "minus2kb.paf"
    plus_paf = outdir / "plus2kb.paf"

    run_quiet(
        f"minimap2 -x asm5 {input_fasta} {minus_fa} > {minus_paf}",
        shell=True,
        check=True,
    )
    run_quiet(
        f"minimap2 -x asm5 {input_fasta} {plus_fa} > {plus_paf}",
        shell=True,
        check=True,
    )

    print("\n[STEP 3] Lifting over candidate L1 coordinates")

    lifted = _liftover_l1s(minus_paf, plus_paf, ref_bed, min_length)

    print("\n[STEP 4] Parsing structural variants")

    deletions, insertions = _parse_sv(Path(sv_file), Path(log) if log else None)

    print("\n[STEP 5] Validating candidate L1 sequences")

    status, ins_seqs = _classify_sv(lifted, deletions, insertions, outdir)

    candidate_fa = outdir / "candidates.fa"
    _extract_sequences(Path(input_fasta), lifted, status, ins_seqs, candidate_fa, min_length)

    if candidate_fa.stat().st_size > 0:
        run_quiet(["RepeatMasker", "-lib", str(Path("data") / "L1rp.fa"), str(candidate_fa)], check=True)
        rm_out = candidate_fa.with_suffix(".fa.out")
        l1_names = _parse_repeatmasker(rm_out, 6019, 0.9)
    else:
        l1_names = set()

    print("\n[STEP 6] Writing output table")

    out_table = outdir / "HapLongLINErSV.txt"
    with open(out_table, "w") as out:
        for chrom, start, end, name, length, strand, _, _ in lifted:
            stat = status.get(name, "present")
            l1flag = "L1" if name in l1_names else "NA"
            out.write(
                f"{chrom}\t{start}\t{end}\t{name}\t{length}\t{strand}\t{stat}\t{l1flag}\n"
            )

    print(f"Module 2 completed. Results in {out_table}")
