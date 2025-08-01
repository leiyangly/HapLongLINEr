from pathlib import Path
from haplongliner.module1_RM import _write_rm_sequences, _fix_getorf_headers


def test_write_rm_sequences(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">chr1\n" + "A" * 100 + "\n")
    bed = tmp_path / "sites.bed"
    bed.write_text("chr1\t10\t20\tL1a\t10\t+\nchr1\t30\t40\tL1b\t10\t-\n")
    out = tmp_path / "rm.fa"
    _write_rm_sequences(fa, bed, out)
    headers = [l.strip() for l in open(out) if l.startswith(">")]
    assert headers == [">chr1,10,20,10,+,RPM", ">chr1,30,40,10,-,RPM"]


def test_fix_getorf_headers(tmp_path):
    fa = tmp_path / "orf.fa"
    fa.write_text(
        ">scaf_A_B_1 [5 - 10]\nAAAA\n>name_2 [1 - 2]\nTT\n"
    )
    _fix_getorf_headers(fa)
    headers = [l.strip() for l in open(fa) if l.startswith(">")]
    assert headers == [">scaf_A_B,1,5,10", ">name,2,1,2"]


def test_fix_getorf_headers_restore_commas(tmp_path):
    orig = tmp_path / "cand.fa"
    orig.write_text(">a,b,c,d,-\nAA\n")
    orf = tmp_path / "orf.fa"
    orf.write_text(">a_b_c_d_-_1 [1 - 2]\nAA\n")
    _fix_getorf_headers(orf, orig)
    headers = [l.strip() for l in open(orf) if l.startswith(">")]
    assert headers == [">a,b,c,d,-,1,1,2"]
