from pathlib import Path

import pytest

from haplongliner.mm_mode import (
    _collect_mm_candidates,
    _mm_seed_min_span,
    _project_mm_repeatmasker,
    _run_candidate_repeatmasker,
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
def test_run_candidate_repeatmasker_uses_default_pa(monkeypatch, tmp_path):
    fa = tmp_path / "seed.fa"
    fa.write_text(">scaf1,100,220,+\n" + "A" * 120 + "\n")

    def fake_run_quiet(cmd, check=True, cwd=None, **kwargs):
        assert cmd[0:7] == [
            "RepeatMasker",
            "-e",
            "rmblast",
            "-pa",
            "4",
            "-species",
            "human",
        ]
        short_out = Path(cwd) / "seed_short.fa.out"
        short_out.write_text("")

    monkeypatch.setattr("haplongliner.mm_mode.run_quiet", fake_run_quiet)
    out = _run_candidate_repeatmasker(fa, 4)
    assert out == tmp_path / "seed.fa.out"


def test_run_candidate_repeatmasker_respects_pa_override(monkeypatch, tmp_path):
    fa = tmp_path / "seed.fa"
    fa.write_text(">scaf1,100,220,+\n" + "A" * 120 + "\n")

    def fake_run_quiet(cmd, check=True, cwd=None, **kwargs):
        assert cmd[0:7] == [
            "RepeatMasker",
            "-e",
            "rmblast",
            "-pa",
            "8",
            "-species",
            "human",
        ]
        short_out = Path(cwd) / "seed_short.fa.out"
        short_out.write_text("")

    monkeypatch.setattr("haplongliner.mm_mode.run_quiet", fake_run_quiet)
    _run_candidate_repeatmasker(fa, 8)
