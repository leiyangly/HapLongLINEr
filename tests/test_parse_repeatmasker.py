from haplongliner.sv_mode import _parse_repeatmasker


def test_parse_repeatmasker_normalizes_names(tmp_path):
    out = tmp_path / "cand.out"
    # Minimal RepeatMasker-like line with enough columns and an L1 family
    out.write_text(
        "    1 0 0 0 seq1,chr1,100,200,+ 0 0 0 L1HS L1 1 100 0 0\n"
    )
    names = _parse_repeatmasker(out, l1_len=100, cov_thresh=0.0)
    assert names == {"seq1,chr1,100"}

