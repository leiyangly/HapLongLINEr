import pytest

from haplongliner.mm_mode import (
    _collect_mm_candidates,
    _mm_seed_min_span,
    _project_mm_repeatmasker,
    _validate_mm_te,
)


def test_mm_seed_min_span_defaults():
    assert _mm_seed_min_span(0) == 1000
    assert _mm_seed_min_span(2000) == 1000
    assert _mm_seed_min_span(5000) == 4000


def test_validate_mm_te_rejects_non_l1():
    with pytest.raises(SystemExit, match="L1-family"):
        _validate_mm_te("L1,SVA")


def test_collect_mm_candidates_merges_close_hits(tmp_path):
    paf = tmp_path / "seed.paf"
    paf.write_text(
        "scaf1\t10000\t100\t2500\t+\tL1rp\t6019\t0\t2400\t2300\t2400\t60\n"
        "scaf1\t10000\t2550\t5200\t+\tL1rp\t6019\t0\t2600\t2550\t2650\t60\n"
        "scaf1\t10000\t7000\t8900\t-\tL1rp\t6019\t0\t1900\t1800\t1900\t60\n"
        "scaf1\t10000\t9000\t9500\t-\tL1rp\t6019\t0\t500\t450\t500\t60\n"
    )
    out = tmp_path / "seed.bed"

    count = _collect_mm_candidates(paf, out, min_span=1500, merge_gap=100)

    assert count == 2
    lines = out.read_text().strip().splitlines()
    assert lines == [
        "scaf1\t100\t5200\tL1rp_seed\t5100\t+",
        "scaf1\t7000\t8900\tL1rp_seed\t1900\t-",
    ]


def test_project_mm_repeatmasker_projects_plus_and_reverse_candidates(tmp_path):
    out = tmp_path / "seed.fa.out"
    out.write_text(
        "500 5.0 0.0 0.0 scaf1,100,220,+ 1 120 (0) + L1HS LINE/L1 1 120 (0) 1\n"
        "500 2.0 0.0 0.0 scaf2,500,620,- 21 120 (0) + L1PA3 LINE/L1 1 100 (0) 1\n"
    )
    bed = tmp_path / "cand.bed"

    count = _project_mm_repeatmasker(out, bed, min_length=100, te="L1,L1PA3")

    assert count == 2
    assert bed.read_text().strip().splitlines() == [
        "scaf1\t100\t220\tL1HS\t120\t+",
        "scaf2\t500\t600\tL1PA3\t100\t-",
    ]


def test_project_mm_repeatmasker_filters_non_requested_family(tmp_path):
    out = tmp_path / "seed.fa.out"
    out.write_text(
        "500 5.0 0.0 0.0 scaf1,100,220,+ 1 120 (0) + L1PA4 LINE/L1 1 120 (0) 1\n"
    )
    bed = tmp_path / "cand.bed"

    count = _project_mm_repeatmasker(out, bed, min_length=100, te="L1,L1PA3")

    assert count == 0
    assert bed.read_text() == ""


def test_project_mm_repeatmasker_merges_fragmented_young_l1(tmp_path):
    out = tmp_path / "seed.fa.out"
    out.write_text(
        "19083 1.2 0.0 0.0 scaf1,100,6096,+ 1 2127 (3869) + L1HS LINE/L1 10 2136 (3896) 1\n"
        "28563 1.6 0.0 0.0 scaf1,100,6096,+ 1978 5996 (0) + L1PA2 LINE/L1 2110 6128 (27) 2\n"
    )
    bed = tmp_path / "cand.bed"

    count = _project_mm_repeatmasker(out, bed, min_length=5000, te="L1,L1PA3")

    assert count == 1
    assert bed.read_text().strip() == "scaf1\t100\t6096\tL1PA2\t5996\t+"
