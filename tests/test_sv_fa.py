from pathlib import Path
from haplongliner.sv_mode import _write_sv_sequences


def test_write_sv_sequences(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">chr1\n" + "A" * 100 + "\n")
    lifted = [
        ("chr1", 10, 20, "chr1:10-20", 10, "+", "chr1,10,20,+", "chr1,10,20,+", 10, 20),
        ("chr1", 30, 40, "chr1:30-40", 10, "+", "chr1,30,40,+", "chr1,30,40,+", 30, 40),
    ]
    status = {"chr1:10-20": "disrupted", "chr1:30-40": "disrupted"}
    ins_seqs = {"chr1:30-40": "T" * 10}
    out = tmp_path / "sv.fa"
    _write_sv_sequences(
        fa, lifted, status, ins_seqs, out, 5, 1000, present={"chr1_30_40"}
    )
    headers = [l.strip() for l in open(out) if l.startswith(">")]
    assert headers == [">chr1,10,20,10,+,ALN", ">chr1,30,40,10,+,INS"]
