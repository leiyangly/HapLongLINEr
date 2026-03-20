from pathlib import Path

from haplongliner.sv_mode import _extract_sequences


def test_extract_sequences_names(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">chr1\n" + "A" * 100 + "\n")
    lifted = [
        ("chr1", 10, 20, "chr1:10-20", 10, "+", "chr1,10,20,+", "chr1,10,20,+", 10, 20),
        ("chr1", 30, 40, "chr1:30-40", 10, "-", "chr1,30,40,-", "chr1,30,40,-", 30, 40),
    ]
    out_lift = tmp_path / "lift.fa"
    out_ins = tmp_path / "ins.fa"
    _extract_sequences(fa, lifted, out_lift, out_ins, [], 5, 1000)
    headers = [l.strip() for l in open(out_lift) if l.startswith(">")]
    assert headers == [
        ">chr1,10,20,+,chr1:10-20,chr1,10,20,+",
        ">chr1,30,40,-,chr1:30-40,chr1,30,40,-",
    ]
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
    assert headers == [">chr1,50,70,+,INS_chr1_50"]
    assert out_lift.read_text() == ""


def test_extract_sequences_filters_short_insertions(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">chr1\n" + "A" * 100 + "\n")
    lifted = [
        ("chr1", 10, 20, "chr1:10-20", 10, "+", "chr1,10,20,+", "chr1,10,20,+", 10, 20)
    ]
    insertions = [("chr1", 80, 82, "G" * 2)]
    out_lift = tmp_path / "lift.fa"
    out_ins = tmp_path / "ins.fa"
    _extract_sequences(fa, lifted, out_lift, out_ins, insertions, 50, 1000)
    lines = [l.strip() for l in open(out_lift)]
    assert lines == [">chr1,10,20,+,chr1:10-20,chr1,10,20,+", "A" * 100]
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


def test_extract_sequences_expands_short_lifted_loci_for_orf_rescue(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">chr1\n" + "A" * 200 + "\n")
    lifted = [
        (
            "chr1",
            100,
            102,
            "chr1:100-102",
            2,
            "+",
            "chr1,100,102,+",
            "chr1,100,102,+",
            100,
            102,
        )
    ]
    out_lift = tmp_path / "lift.fa"
    out_ins = tmp_path / "ins.fa"

    eval_intervals = _extract_sequences(fa, lifted, out_lift, out_ins, [], 50, 1000)

    lines = [l.strip() for l in open(out_lift)]
    assert lines == [
        ">chr1,100,102,+,chr1:100-102,chr1,100,102,+",
        "A" * 150,
    ]
    assert eval_intervals == {"chr1:100-102": (26, 176)}
