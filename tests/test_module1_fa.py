from pathlib import Path
from haplongliner.module1_RM import _write_rm_sequences


def test_write_rm_sequences(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">chr1\n" + "A" * 100 + "\n")
    bed = tmp_path / "sites.bed"
    bed.write_text("chr1\t10\t20\tL1a\t10\t+\nchr1\t30\t40\tL1b\t10\t-\n")
    out = tmp_path / "rm.fa"
    _write_rm_sequences(fa, bed, out)
    headers = [l.strip() for l in open(out) if l.startswith(">")]
    assert headers == [">chr1;10;20;10;+;ALN", ">chr1;30;40;10;-;ALN"]
