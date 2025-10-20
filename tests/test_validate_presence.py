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
    present, annotations = _validate_presence(cand, min_length=0)
    assert present == set()
    assert annotations == {}


def test_all_empty_sequences_skip_repeatmasker(monkeypatch, tmp_path):
    cand = tmp_path / "cand.fa"
    cand.write_text(">a\n>b\n")

    called = {"flag": False}

    def fake_run_quiet(cmd, check=True, cwd=None, **kwargs):  # pragma: no cover - should not be called
        called["flag"] = True

    monkeypatch.setattr("haplongliner.sv_mode.run_quiet", fake_run_quiet)
    present, annotations = _validate_presence(cand, min_length=0)
    assert present == set()
    assert annotations == {}
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
    present, annotations = _validate_presence(cand, min_length=0)
    assert present == set()
    assert annotations == {}

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
    present, annotations = _validate_presence(Path("sub") / "cand.fa", min_length=0)
    assert present == set()
    assert annotations == {}
    assert (sub / "cand.fa.out").exists()


def test_validate_presence_returns_annotations(monkeypatch, tmp_path):
    cand = tmp_path / "cand.fa"
    cand.write_text(
        ">chr1,10,110,+,L1HS\n"
        + "A" * 100
        + "\n"
        + ">chr1,200,320,-,L1PA3\n"
        + "C" * 120
        + "\n"
    )

    data_dir = tmp_path / "data"
    data_dir.mkdir()
    (data_dir / "L1rp.fa").write_text(">ref\n" + "A" * 6000 + "\n")

    def fake_run_quiet(cmd, check=True, cwd=None, **kwargs):
        query = Path(cwd) / cmd[-1]
        out_path = query.with_suffix(query.suffix + ".out")
        out_path.write_text(
            "   500   5.0  0.0  0.0  chr1,10,110,+,L1HS      1    90   (10)  + L1HS    LINE/L1    1   90  (0)   1\n"
            "   400  10.0  0.0  0.0  chr1,200,320,-,L1PA3   10   120  (0)  C L1PA3   LINE/L1  (100)  310 210  2\n"
        )

    monkeypatch.setattr("haplongliner.sv_mode.run_quiet", fake_run_quiet)
    monkeypatch.chdir(tmp_path)

    present, annotations = _validate_presence(Path("cand.fa"), min_length=0)
    assert present == {"chr1_10_110", "chr1_200_320"}
    assert set(annotations) == {"chr1_10_110", "chr1_200_320"}

    hit = annotations["chr1_10_110"]
    assert hit.family == "L1HS"
    assert pytest.approx(hit.identity, rel=1e-3) == 95.0
    assert pytest.approx(hit.coverage, rel=1e-3) == 90.0

    hit = annotations["chr1_200_320"]
    assert hit.family == "L1PA3"
    assert pytest.approx(hit.identity, rel=1e-3) == 90.0
    # second sequence is 120 bp long with 111 bp covered (positions 10-120)
    assert pytest.approx(hit.coverage, rel=1e-3) == pytest.approx(111 / 120 * 100, rel=1e-3)
