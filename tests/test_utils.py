import gzip
import gzip
from pathlib import Path

import pytest

from haplongliner.utils import (
    verify_fasta_file,
    verify_repeatmasker_file,
    verify_sv_file,
    _fix_blast_query_names,
    _candidate_key_from_query,
    sort_bed,
    append_fasta,
    ensure_pyfaidx_fasta,
)


def test_verify_fasta_file_accepts_gz(tmp_path):
    src = Path('tests/test.genome.HG00410.1.fa')
    gz = tmp_path / 'test.fa.gz'
    with open(src, 'rb') as f_in, gzip.open(gz, 'wb') as f_out:
        f_out.write(f_in.read())
    verify_fasta_file(str(src))
    verify_fasta_file(str(gz))


def test_ensure_pyfaidx_fasta_decompresses_standard_gz(tmp_path):
    src = tmp_path / "test.fa.gz"
    with gzip.open(src, "wt") as out:
        out.write(">chr1\nACGT\n")
    out_path = ensure_pyfaidx_fasta(src)
    assert out_path == tmp_path / "test.fa"
    assert out_path.read_text() == ">chr1\nACGT\n"


def test_ensure_pyfaidx_fasta_keeps_plain_fasta(tmp_path):
    src = tmp_path / "plain.fa"
    src.write_text(">chr1\nACGT\n")
    out_path = ensure_pyfaidx_fasta(src)
    assert out_path == src


def test_verify_repeatmasker_file_accepts_out_and_bed(tmp_path):
    out = Path('tests/test.genome.HG00410.1.out')
    gz = tmp_path / 'test.out.gz'
    with open(out, 'rb') as f_in, gzip.open(gz, 'wb') as f_out:
        f_out.write(f_in.read())
    verify_repeatmasker_file(str(out))
    verify_repeatmasker_file(str(gz))
    bed = tmp_path / 'test.bed'
    bed.write_text('chr1\t1\t2\tname\t1\t+\n')
    verify_repeatmasker_file(str(bed))


def test_verify_repeatmasker_file_bad(tmp_path):
    bad = tmp_path / 'bad.txt'
    bad.write_text('malformed line')
    with pytest.raises(SystemExit):
        verify_repeatmasker_file(str(bad))


def test_verify_sv_file_accepts_gz_and_plain(tmp_path):
    src_gz = Path('tests/test.genome.HG00410.vcf.gz')
    plain = tmp_path / 'test.vcf'
    with gzip.open(src_gz, 'rb') as f_in, open(plain, 'wb') as f_out:
        f_out.write(f_in.read())
    verify_sv_file(str(src_gz))
    verify_sv_file(str(plain))


def test_verify_sv_file_bad(tmp_path):
    bad = tmp_path / 'bad.txt'
    bad.write_text('bad line')
    with pytest.raises(SystemExit):
        verify_sv_file(str(bad))


def test_fix_blast_query_names(tmp_path):
    path = tmp_path / "test.blast"
    path.write_text(
        "011138A1_186_phaseblock_2_1087928_1093929_+_6_876_1889\tsub\n"
    )
    _fix_blast_query_names(path)
    fixed = path.read_text().strip()
    assert fixed.startswith(
        "011138A1,186_phaseblock_2,1087928,1093929,+,6,876,1889\tsub"
    )


@pytest.mark.parametrize(
    "query,expected",
    [
        (
            "chr1:1-2,scaf,10,20,+,1,5,100",
            "chr1:1-2,scaf,10",
        ),
        (
            "chr1:1-2_scaf_10_20_+_1_5_100",
            "chr1:1-2,scaf,10",
        ),
        (
            "INS_chr1_100_chr1_100_200_._1_20_200",
            "INS,chr1_100_chr1,100",
        ),
        (
            "simple_name",
            "simple_name",
        ),
    ],
)
def test_candidate_key_from_query(query, expected):
    assert _candidate_key_from_query(query) == expected


def test_sort_bed(tmp_path):
    path = tmp_path / "unsorted.bed"
    path.write_text(
        "chr2\t50\t60\nchr1\t30\t40\nchr1\t10\t20\n"
    )
    sort_bed(path)
    assert path.read_text().splitlines() == [
        "chr1\t10\t20",
        "chr1\t30\t40",
        "chr2\t50\t60",
    ]


def test_sort_bed_with_na_start(tmp_path):
    path = tmp_path / "unsorted_na.bed"
    path.write_text("chr1\tNA\t20\nchr1\t10\t20\n")
    sort_bed(path)
    assert path.read_text().splitlines() == [
        "chr1\t10\t20",
        "chr1\tNA\t20",
    ]


def test_append_fasta(tmp_path):
    dst = tmp_path / "orf.fa"
    dst.write_text(">a,1,2,+,1,1,2\nAA\n")
    src = tmp_path / "cand.fa"
    src.write_text(">a,1,2,+\nAAAA\n")
    append_fasta(dst, src)
    assert dst.read_text().splitlines() == [
        ">a,1,2,+,1,1,2",
        "AA",
        ">a,1,2,+",
        "AAAA",
    ]
