from pathlib import Path
from haplongliner.sv_mode import _extract_sequences


def test_extract_sequences_names(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">chr1\n" + "A" * 100 + "\n")
    lifted = [
        ("chr1", 10, 20, "L1a", 10, "+", "chr1,10,20,+", "chr1,10,20,+", 10, 20),
        ("chr1", 30, 40, "L1b", 10, "-", "chr1,30,40,-", "chr1,30,40,-", 30, 40),
    ]
    out_lift = tmp_path / "lift.fa"
    out_ins = tmp_path / "ins.fa"
    _extract_sequences(fa, lifted, out_lift, out_ins, [], 5, 1000)
    headers = [l.strip() for l in open(out_lift) if l.startswith(">")]
    assert headers == [">L1a,chr1,10,20,+", ">L1b,chr1,30,40,-"]
    assert out_ins.read_text() == ""


def test_extract_sequences_with_insertions(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">chr1\n" + "A" * 100 + "\n")
    lifted: list = []
    insertions = [("chr1", 50, 70, "T" * 20)]
    out_lift = tmp_path / "lift.fa"
    out_ins = tmp_path / "extra.fa"
    _extract_sequences(fa, lifted, out_lift, out_ins, insertions, 5, 1000)
    headers = [l.strip() for l in open(out_ins) if l.startswith(">")]
    assert headers == [">INS_chr1_50,chr1,50,70,."]
    assert out_lift.read_text() == ""


def test_extract_sequences_filters_short_insertions(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">chr1\n" + "A" * 100 + "\n")
    lifted = [
        ("chr1", 10, 20, "L1a", 10, "+", "chr1,10,20,+", "chr1,10,20,+", 10, 20)
    ]
    insertions = [("chr1", 80, 82, "G" * 2)]
    out_lift = tmp_path / "lift.fa"
    out_ins = tmp_path / "ins.fa"
    _extract_sequences(fa, lifted, out_lift, out_ins, insertions, 50, 1000)
    lines = [l.strip() for l in open(out_lift)]
    assert lines == [">L1a,chr1,10,20,+", "A" * 10]
    assert out_ins.read_text() == ""


def test_extract_sequences_filters_long_insertions(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">chr1\n" + "A" * 100 + "\n")
    insertions = [("chr1", 50, 90, "T" * 40)]
    out_lift = tmp_path / "lift.fa"
    out_ins = tmp_path / "long.fa"
    _extract_sequences(fa, [], out_lift, out_ins, insertions, 5, 30)
    assert out_lift.read_text() == ""
    assert out_ins.read_text() == ""
