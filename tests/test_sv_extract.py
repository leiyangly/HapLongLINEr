from pathlib import Path
from haplongliner.sv_mode import _extract_sequences


def test_extract_sequences_names(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">chr1\n" + "A" * 100 + "\n")
    lifted = [
        ("chr1", 10, 20, "L1a", 10, "+", "chr1,10,20,+", "chr1,10,20,+", 10, 20),
        ("chr1", 30, 40, "L1b", 10, "-", "chr1,30,40,-", "chr1,30,40,-", 30, 40),
    ]
    status = {"L1a": "present", "L1b": "present"}
    out = tmp_path / "out.fa"
    _extract_sequences(fa, lifted, status, {}, out, 5)
    headers = [l.strip() for l in open(out) if l.startswith(">")]
    assert headers == [">L1a,chr1,10,20,+", ">L1b,chr1,30,40,-"]


def test_extract_sequences_with_extras(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">chr1\n" + "A" * 100 + "\n")
    lifted: list = []
    status = {}
    extras = [("chr1", 50, 70, "T" * 20, "INS_chr1_50")]
    out = tmp_path / "extra.fa"
    _extract_sequences(fa, lifted, status, {}, out, 5, extras)
    headers = [l.strip() for l in open(out) if l.startswith(">")]
    assert headers == [">INS_chr1_50,chr1,50,70,."]


def test_extract_sequences_insertions_unfiltered(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">chr1\n" + "A" * 100 + "\n")
    lifted: list = []
    status = {}
    ins = {"L1a": "T" * 10}
    extras = [("chr1", 80, 82, "G" * 2, "INS_chr1_80")]
    out = tmp_path / "ins.fa"
    _extract_sequences(fa, lifted, status, ins, out, 50, extras)
    lines = [l.strip() for l in open(out)]
    assert lines == [">L1a_ins", "T" * 10, ">INS_chr1_80,chr1,80,82,.", "G" * 2]
