from haplongliner.sv_mode import _parse_repeatmasker


def test_parse_repeatmasker_normalizes_names(tmp_path):
    out = tmp_path / "cand.out"
    out.write_text(
        "19083 1.2 0.0 0.0 seq1,chr1,100,6096,+ 1 2127 (3869) + L1HS LINE/L1 10 2136 (3896) 1\n"
        "28563 1.6 0.0 0.0 seq1,chr1,100,6096,+ 1978 5996 (0) + L1PA2 LINE/L1 2110 6128 (27) 2\n"
    )
    seq_lengths = {("seq1", "chr1", "100"): 5996}
    present, annotations = _parse_repeatmasker(out, seq_lengths, cov_thresh=0.95)
    assert present == {"seq1_chr1_100"}
    hit = annotations["seq1_chr1_100"]
    assert hit.family == "L1PA2"
    assert hit.identity > 98.0
    assert hit.coverage == 100.0
