from haplongliner.sv_mode import _parse_repeatmasker


def test_parse_repeatmasker_normalizes_names(tmp_path):
    out = tmp_path / "cand.out"
    # Minimal RepeatMasker-like line with enough columns and an L1 family
    out.write_text(
        "   500   5.0  0.0  0.0  seq1,chr1,100,200,+  1  100  (0)  +  L1HS  LINE/L1  1  100  (0)  1\n"
    )
    seq_lengths = {("seq1", "chr1", "100"): 100}
    present, annotations = _parse_repeatmasker(out, seq_lengths, cov_thresh=0.0)
    assert present == {"seq1_chr1_100"}
    hit = annotations["seq1_chr1_100"]
    assert hit.family == "L1HS"
    assert hit.identity == 95.0
    assert hit.coverage == 100.0


def test_parse_repeatmasker_accepts_unindented_lines(tmp_path):
    out = tmp_path / "cand.out"
    out.write_text(
        "500 5.0 0.0 0.0 seq1,chr1,100,200,+ 1 100 (0) + L1HS LINE/L1 1 100 (0) 1\n"
    )
    seq_lengths = {("seq1", "chr1", "100"): 100}
    present, annotations = _parse_repeatmasker(out, seq_lengths, cov_thresh=0.0)
    assert present == {"seq1_chr1_100"}
    hit = annotations["seq1_chr1_100"]
    assert hit.family == "L1HS"
    assert hit.identity == 95.0
    assert hit.coverage == 100.0
