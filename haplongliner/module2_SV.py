import re
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple


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


def _liftover_l1s(minus_paf: Path, plus_paf: Path, ref_bed: Path) -> List[Tuple[str, int, int, str, int, str]]:
    """Infer target assembly coordinates for each L1 listed in ``ref_bed``."""
    minus = _read_paf(minus_paf)
    plus = _read_paf(plus_paf)

    lifted: List[Tuple[str, int, int, str, int, str]] = []
    with open(ref_bed) as fh:
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            chrom, start, end, name, strand = line.strip().split()[:5]
            m = minus.get(f"{name}_-2kb")
            p = plus.get(f"{name}_+2kb")
            if not m or not p:
                continue
            if m[5] != p[5] or m[4] != p[4]:
                continue
            tname = m[5]
            orient = m[4]
            if orient == '+':
                start_t = int(m[8])
                end_t = int(p[7])
            else:
                start_t = int(p[8])
                end_t = int(m[7])
            if end_t < start_t:
                start_t, end_t = end_t, start_t
            length = int(end) - int(start)
            lifted.append((tname, start_t, end_t, name, length, orient))
    return lifted


def _write_bed(entries: List[Tuple[str, int, int, str, int, str]], path: Path) -> None:
    with open(path, 'w') as out:
        for chrom, start, end, name, length, strand in entries:
            out.write(f"{chrom}\t{start}\t{end}\t{name}\t{length}\t{strand}\n")


def _parse_sv(sv_path: Path) -> Tuple[List[Tuple[str, int, int]], List[Tuple[str, int, int]]]:
    """Parse a simple VCF or BED SV file and return deletion and insertion regions."""
    deletions: List[Tuple[str, int, int]] = []
    insertions: List[Tuple[str, int, int]] = []
    is_vcf = sv_path.suffix.lower().endswith('vcf') or sv_path.suffix.lower() == '.gz'
    with open(sv_path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip().split('\t')
            if is_vcf or len(fields) > 5:
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
                    insertions.append((chrom, pos, pos + 1))
            else:
                chrom, start, end = fields[:3]
                svtype = fields[3].upper() if len(fields) > 3 else ''
                if svtype == 'DEL':
                    deletions.append((chrom, int(start), int(end)))
                elif svtype == 'INS':
                    insertions.append((chrom, int(start), int(end)))
    return deletions, insertions


def _bedtools_intersect(a: Path, b: Path, output: Path) -> None:
    with open(output, 'w') as out:
        subprocess.run(['bedtools', 'intersect', '-wa', '-wb', '-a', str(a), '-b', str(b)], check=True, stdout=out)


def _classify_deletions(lifted: List[Tuple[str, int, int, str, int, str]], deletions: List[Tuple[str, int, int]], outdir: Path) -> Dict[str, str]:
    lift_bed = outdir / 'lifted.bed'
    del_bed = outdir / 'sv_del.bed'
    _write_bed(lifted, lift_bed)
    with open(del_bed, 'w') as out:
        for chrom, start, end in deletions:
            out.write(f"{chrom}\t{start}\t{end}\n")
    inter_file = outdir / 'intersect.bed'
    _bedtools_intersect(lift_bed, del_bed, inter_file)

    status: Dict[str, str] = {name: 'present' for _, _, _, name, _, _ in lifted}
    with open(inter_file) as fh:
        for line in fh:
            f = line.strip().split('\t')
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
                status[name] = 'missing'
            else:
                status[name] = 'absent'
    return status


def _extract_sequences(fasta: Path, lifted: List[Tuple[str, int, int, str, int, str]], status: Dict[str, str], out_fa: Path) -> None:
    bed_path = out_fa.with_suffix('.bed')
    with open(bed_path, 'w') as bed:
        for chrom, start, end, name, length, _ in lifted:
            if status.get(name) in {'missing', 'absent'} and abs((end - start) - length) / length < 0.1:
                bed.write(f"{chrom}\t{start}\t{end}\t{name}\n")
    if bed_path.stat().st_size > 0:
        cmd = f"seqtk subseq {fasta} {bed_path} | seqtk seq -U -l 0 - > {out_fa}"
        subprocess.run(cmd, shell=True, check=True)
    else:
        out_fa.touch()


def _parse_repeatmasker(out_file: Path) -> List[str]:
    hits: List[str] = []
    with open(out_file) as fh:
        for line in fh:
            if line.startswith(' '):
                parts = line.split()
                if len(parts) >= 10 and re.search(r'L1', parts[9]):
                    hits.append(parts[4])
    return hits


def run_module2(input_fasta: str, sv_file: str, l1ref_fasta: str, output_bed: str) -> None:
    print(
        f"Module 2 running with:\n  Input: {input_fasta}\n  SV: {sv_file}\n  L1 Reference: {l1ref_fasta}\n  Output: {output_bed}"
    )

    out_path = Path(output_bed)
    outdir = out_path.parent
    outdir.mkdir(parents=True, exist_ok=True)

    minus_fa = Path('data') / '-2kb.fa'
    plus_fa = Path('data') / '+2kb.fa'

    minus_paf = outdir / 'minus2kb.paf'
    plus_paf = outdir / 'plus2kb.paf'

    subprocess.run(f"minimap2 -x asm5 {input_fasta} {minus_fa} > {minus_paf}", shell=True, check=True)
    subprocess.run(f"minimap2 -x asm5 {input_fasta} {plus_fa} > {plus_paf}", shell=True, check=True)

    ref_bed = Path('data') / 'HPRC_L1_hs_v2_v2fl.bed'
    lifted = _liftover_l1s(minus_paf, plus_paf, ref_bed)

    deletions, insertions = _parse_sv(Path(sv_file))
    status = _classify_deletions(lifted, deletions, outdir)

    candidate_fa = outdir / 'candidates.fa'
    _extract_sequences(Path(input_fasta), lifted, status, candidate_fa)

    if candidate_fa.stat().st_size > 0:
        subprocess.run(['RepeatMasker', str(candidate_fa)], check=True)
        rm_out = candidate_fa.with_suffix('.fa.out')
        l1_names = set(_parse_repeatmasker(rm_out))
    else:
        l1_names = set()

    with open(output_bed, 'w') as out:
        for chrom, start, end, name, length, strand in lifted:
            stat = status.get(name, 'present')
            l1flag = 'L1' if name in l1_names else 'NA'
            out.write(f"{chrom}\t{start}\t{end}\t{name}\t{length}\t{strand}\t{stat}\t{l1flag}\n")

    print(f"Module 2 completed. Results in {output_bed}")
