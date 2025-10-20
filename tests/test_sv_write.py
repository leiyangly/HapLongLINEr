import sys
import types


if "Bio" not in sys.modules:
    bio = types.ModuleType("Bio")

    class _SeqIOStub:
        @staticmethod
        def parse(*_args, **_kwargs):  # pragma: no cover - simple stub
            return []

    bio.SeqIO = _SeqIOStub()
    sys.modules["Bio"] = bio


from haplongliner.sv_mode import _write_sv_sequences, _infer_te_type


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


def test_infer_te_type_handles_various_names():
    assert _infer_te_type("5398_L1HS_intact") == "L1HS"
    assert _infer_te_type("INS_chr1_200") == "INS"
    assert _infer_te_type("L1PA3") == "L1PA3"
    assert _infer_te_type("L1HS_3end") == "L1HS_3end"
    assert _infer_te_type("chr1:0-10") == "chr1:0-10"

