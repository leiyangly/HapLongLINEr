from pathlib import Path
import subprocess

import pytest

from haplongliner.sv_mode import _validate_presence
import haplongliner.sv_mode as sv_mode


def test_empty_sequences_are_filtered(monkeypatch, tmp_path):
    cand = tmp_path / "cand.fa"
    cand.write_text(">a\nACGT\n>b\n>c\nTTTT\n")

    def fake_run_quiet(cmd, check=True, cwd=None, **kwargs):  # pragma: no cover - behavior verified via assertions
        query = Path(cwd) / cmd[-1]
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
        assert Path(cwd).resolve() == tmp_path
        query = Path(cwd) / cmd[-1]
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


def test_validate_presence_relative_path_subdir(monkeypatch, tmp_path):
    sub = tmp_path / "sub"
    sub.mkdir()
    cand = sub / "cand.fa"
    cand.write_text(">a\nACGT\n")

    data_dir = tmp_path / "data"
    data_dir.mkdir()
    (data_dir / "L1rp.fa").write_text(">ref\nAAAA\n")

    def fake_run_quiet(cmd, check=True, cwd=None, **kwargs):  # pragma: no cover - behavior validated via assertions
        assert Path(cwd).resolve() == sub
        assert cmd[-1] == "cand_short.fa"
        query = Path(cwd) / cmd[-1]
        out_path = query.with_suffix(query.suffix + ".out")
        out_path.write_text("")

    monkeypatch.setattr("haplongliner.sv_mode.run_quiet", fake_run_quiet)
    monkeypatch.chdir(tmp_path)
    assert _validate_presence(Path("sub") / "cand.fa", min_length=0) == set()
    assert (sub / "cand.fa.out").exists()
