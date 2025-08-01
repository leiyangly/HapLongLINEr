from pathlib import Path
from haplongliner.process_orf import process_orf_fasta


def test_process_orf_fasta(tmp_path):
    fa = tmp_path / "orf.fa"
    fa.write_text(
        ">chr1,10,20,10,+,RPM,1,5,10\nAA\n"
        ">chr2,30,40,10,-,RPM,2,3,8\nTT\n"
        ">chr3,50,60,10,+,RPM,3,9,4\nCC\n"
    )
    bed = tmp_path / "out.bed"
    process_orf_fasta(fa, bed)
    lines = [l.strip() for l in open(bed)]
    assert lines == [
        "chr1_10_20\t4\t10\t+\t6\t+",
        "chr2_30_40\t2\t8\t+\t6\t-",
        "chr3_50_60\t3\t9\t-\t6\t+",
    ]


def test_process_orf_fasta_invalid(tmp_path):
    """Headers missing required fields should be ignored."""
    fa = tmp_path / "bad_orf.fa"
    fa.write_text(">bad_header\nAA\n")
    bed = tmp_path / "out.bed"
    process_orf_fasta(fa, bed)
    assert bed.read_text() == ""
