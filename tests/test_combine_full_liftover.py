from haplongliner.rm_module import _combine_full_liftover


def test_combine_full_liftover(tmp_path):
    lifted = tmp_path / "lifted.bed"
    lifted.write_text("chr1\t1000\t1100\tscaf1,100,200,+,L1HS\t0\t+\n")
    intact = tmp_path / "intact.blastp"
    intact.write_text("scaf1,100,200,+\n")
    cand = tmp_path / "cand.bed"
    cand.write_text("scaf1\t100\t200\tL1HS\t100\t+\n")
    out = tmp_path / "out.bed"
    _combine_full_liftover(lifted, intact, cand, out, min_length=50)
    fields = out.read_text().strip().split("\t")
    assert fields[0:7] == ["chr1", "1000", "1100", "L1HS", "100", "+", "intact"]
    assert fields[7] == "scaf1,100,200,100,+,RPM"


def test_combine_full_liftover_minlen_zero(tmp_path):
    lifted = tmp_path / "lifted.bed"
    lifted.write_text("chr1\t1000\t1100\tscaf1,100,200,+,L1HS\t0\t+\n")
    intact = tmp_path / "intact.blastp"
    intact.write_text("scaf1,100,200,+\n")
    cand = tmp_path / "cand.bed"
    cand.write_text("scaf1\t100\t200\tL1HS\t100\t+\n")
    out = tmp_path / "out.bed"
    _combine_full_liftover(lifted, intact, cand, out, min_length=0)
    fields = out.read_text().strip().split("\t")
    assert fields[0:7] == ["chr1", "1000", "1100", "L1HS", "100", "+", "intact"]
    assert fields[7] == "scaf1,100,200,100,+,RPM"
