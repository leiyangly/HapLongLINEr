from haplongliner.sv_mode import _write_sv_sequences


def _make_lifted_entry() -> list:
    return [
        (
            "chr1",
            0,
            10,
            "chr1:0-10",
            10,
            "+",
            "chr1,0,10,+",
            "chr1,0,10,+",
            0,
            10,
        )
    ]


def test_write_sv_sequences_filters_short_insertions(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">chr1\n" + "A" * 100 + "\n")
    out = tmp_path / "sv.fa"
    lifted = _make_lifted_entry()
    status = {"chr1:0-10": "present"}
    ins_seqs = {"chr1:0-10": "A" * 4}  # shorter than min_length
    _write_sv_sequences(fa, lifted, status, ins_seqs, out, 5, 1000)
    assert out.read_text() == ""


def test_write_sv_sequences_writes_long_insertions(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">chr1\n" + "A" * 100 + "\n")
    out = tmp_path / "sv_long.fa"
    lifted = _make_lifted_entry()
    status = {"chr1:0-10": "present"}
    ins_seqs = {"chr1:0-10": "T" * 6}
    _write_sv_sequences(fa, lifted, status, ins_seqs, out, 5, 1000)
    lines = [l.strip() for l in out.read_text().splitlines() if l.strip()]
    assert lines == [">chr1,0,10,10,+,INS", "T" * 6]


def test_write_sv_sequences_filters_long_insertions(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">chr1\n" + "A" * 100 + "\n")
    out = tmp_path / "sv_long_filter.fa"
    lifted = _make_lifted_entry()
    status = {"chr1:0-10": "present"}
    ins_seqs = {"chr1:0-10": "T" * 40}
    _write_sv_sequences(fa, lifted, status, ins_seqs, out, 5, 30)
    assert out.read_text() == ""

