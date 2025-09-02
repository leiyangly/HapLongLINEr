from pathlib import Path
import subprocess

import pytest

from haplongliner.sv_mode import _validate_presence
import haplongliner.sv_mode as sv_mode


def test_empty_sequences_are_filtered(monkeypatch, tmp_path):
    cand = tmp_path / "cand.fa"
    cand.write_text(">a\nACGT\n>b\n>c\nTTTT\n")

    def fake_run_quiet(cmd, check=True, cwd=None, **kwargs):  # pragma: no cover - behavior verified via assertions
        query = Path(cmd[-1])
        mapping = query.with_name("cand.list").read_text()
        assert "a" in mapping
        assert "c" in mapping
        assert "b" not in mapping
        out_path = query.with_suffix(query.suffix + ".out")
        out_path.write_text("")

    monkeypatch.setattr("haplongliner.sv_mode.run_quiet", fake_run_quiet)
    assert _validate_presence(cand, min_length=0) == set()


def test_all_empty_sequences_skip_repeatmasker(monkeypatch, tmp_path):
    cand = tmp_path / "cand.fa"
    cand.write_text(">a\n>b\n")

    called = {"flag": False}

    def fake_run_quiet(cmd, check=True, cwd=None, **kwargs):  # pragma: no cover - should not be called
        called["flag"] = True

    monkeypatch.setattr("haplongliner.sv_mode.run_quiet", fake_run_quiet)
    assert _validate_presence(cand, min_length=0) == set()
    assert not called["flag"]


def test_repeatmasker_runs_in_output_dir(monkeypatch, tmp_path):
    cand = tmp_path / "cand.fa"
    cand.write_text(">a\nACGT\n")

    def fake_run_quiet(cmd, check=True, cwd=None, **kwargs):
        assert cwd == tmp_path
        query = Path(cmd[-1])
        out_path = query.with_suffix(query.suffix + ".out")
        out_path.write_text("")

    monkeypatch.setattr("haplongliner.sv_mode.run_quiet", fake_run_quiet)
    assert _validate_presence(cand, min_length=0) == set()

def test_validate_presence_reports_repeatmasker_failure(tmp_path, monkeypatch):
    fa = tmp_path / "cand.fa"
    fa.write_text(">seq\nAAAA\n")

    def fake_run_quiet(cmd, **kwargs):
        raise subprocess.CalledProcessError(255, cmd)

    monkeypatch.setattr(sv_mode, "run_quiet", fake_run_quiet)

    with pytest.raises(RuntimeError) as excinfo:
        sv_mode._validate_presence(fa)

    assert "RepeatMasker" in str(excinfo.value)
