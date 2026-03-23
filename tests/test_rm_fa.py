from pathlib import Path
import pytest

from haplongliner.rm_mode import (
    _write_rm_sequences,
    _fix_getorf_headers,
    _resolve_repeatmasker_file,
)


def test_write_rm_sequences(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">chr1\n" + "A" * 100 + "\n")
    bed = tmp_path / "sites.bed"
    bed.write_text("chr1\t10\t20\tL1a\t10\t+\nchr1\t30\t40\tL1b\t10\t-\n")
    out = tmp_path / "rm.fa"
    _write_rm_sequences(fa, bed, out)
    headers = [l.strip() for l in open(out) if l.startswith(">")]
    assert headers == [">chr1,10,20,10,+,RPM", ">chr1,30,40,10,-,RPM"]


def test_fix_getorf_headers(tmp_path):
    fa = tmp_path / "orf.fa"
    fa.write_text(
        ">scaf_A_B_1 [5 - 10]\nAAAA\n>name_2 [1 - 2]\nTT\n"
    )
    _fix_getorf_headers(fa)
    headers = [l.strip() for l in open(fa) if l.startswith(">")]
    assert headers == [">scaf_A_B,1,5,10", ">name,2,1,2"]


def test_fix_getorf_headers_restore_commas(tmp_path):
    orig = tmp_path / "cand.fa"
    orig.write_text(">a,b,c,d,-\nAA\n")
    orf = tmp_path / "orf.fa"
    orf.write_text(">a_b_c_d_-_1 [1 - 2]\nAA\n")
    _fix_getorf_headers(orf, orig)
    headers = [l.strip() for l in open(orf) if l.startswith(">")]
    assert headers == [">a,b,c,d,-,1,1,2"]


def test_resolve_repeatmasker_file_uses_explicit_mask(tmp_path, monkeypatch):
    fa = tmp_path / "test.fa"
    fa.write_text(">chr1\nAAAA\n")
    mask = tmp_path / "mask.out"
    mask.write_text(
        "500 5.0 0.0 0.0 chr1 1 4 (0) + L1HS LINE/L1 1 4 (0) 1\n"
    )

    def fail(*args, **kwargs):
        raise AssertionError("RepeatMasker should not run when --mask is provided")

    monkeypatch.setattr("haplongliner.rm_mode.run_quiet", fail)
    got = _resolve_repeatmasker_file(str(fa), str(mask), tmp_path)
    assert got == mask


def test_resolve_repeatmasker_file_runs_repeatmasker_when_mask_missing(tmp_path, monkeypatch):
    fa = tmp_path / "test.fa"
    fa.write_text(">chr1\nAAAA\n")
    calls = []

    def fake_run_quiet(cmd, **kwargs):
        calls.append(cmd)
        (tmp_path / "test.fa.out").write_text(
            "500 5.0 0.0 0.0 chr1 1 4 (0) + L1HS LINE/L1 1 4 (0) 1\n"
        )

    monkeypatch.setattr("haplongliner.rm_mode.run_quiet", fake_run_quiet)
    got = _resolve_repeatmasker_file(str(fa), None, tmp_path)
    assert got == tmp_path / "test.fa.out"
    assert calls == [[
        "RepeatMasker",
        "-e",
        "rmblast",
        "-species",
        "human",
        "-dir",
        str(tmp_path),
        str(fa),
    ]]


def test_resolve_repeatmasker_file_reuses_existing_generated_mask(tmp_path, monkeypatch):
    fa = tmp_path / "test.fa"
    fa.write_text(">chr1\nAAAA\n")
    auto_mask = tmp_path / "test.fa.out"
    auto_mask.write_text(
        "500 5.0 0.0 0.0 chr1 1 4 (0) + L1HS LINE/L1 1 4 (0) 1\n"
    )

    def fail(*args, **kwargs):
        raise AssertionError("RepeatMasker should not rerun when cached output exists")

    monkeypatch.setattr("haplongliner.rm_mode.run_quiet", fail)
    got = _resolve_repeatmasker_file(str(fa), None, tmp_path)
    assert got == auto_mask


def test_resolve_repeatmasker_file_rejects_gz_input_without_mask(tmp_path):
    fa = tmp_path / "test.fa.gz"
    fa.write_text(">chr1\nAAAA\n")
    with pytest.raises(RuntimeError, match="uncompressed FASTA"):
        _resolve_repeatmasker_file(str(fa), None, tmp_path)
