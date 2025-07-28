from pathlib import Path
from haplongliner.module2_SV import _extract_sequences

def test_extract_sequences_names(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">chr1\n" + "A" * 100 + "\n")
    lifted = [
        ("chr1", 10, 20, "L1a", 10, "+", "chr1;10;20;+", "chr1;10;20;+", 10, 20),
        ("chr1", 30, 40, "L1b", 10, "-", "chr1;30;40;-", "chr1;30;40;-", 30, 40),
    ]
    status = {"L1a": "present", "L1b": "present"}
    out = tmp_path / "out.fa"
    _extract_sequences(fa, lifted, status, {}, out, 5)
    headers = [l.strip() for l in open(out) if l.startswith(">")]
    assert headers == [">L1a;chr1;10;20;+", ">L1b;chr1;30;40;-"]
