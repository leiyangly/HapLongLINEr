from pathlib import Path
from haplongliner.sv_mode import _write_sv_sequences


def test_write_sv_sequences(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">chr1\n" + "A" * 100 + "\n")
    lifted = [
        ("chr1", 10, 20, "L1a", 10, "+", "chr1,10,20,+", "chr1,10,20,+", 10, 20),
        ("chr1", 30, 40, "L1b", 10, "+", "chr1,30,40,+", "chr1,30,40,+", 30, 40),
    ]
    status = {"L1a": "present", "L1b": "present"}
    ins_seqs = {"L1b": "T" * 10}
    out = tmp_path / "sv.fa"
    _write_sv_sequences(fa, lifted, status, ins_seqs, out, 5)
    headers = [l.strip() for l in open(out) if l.startswith(">")]
    assert headers == [">chr1,10,20,10,+,ALN", ">chr1,30,40,10,+,INS"]
