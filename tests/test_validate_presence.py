from pathlib import Path

from haplongliner.sv_mode import _validate_presence


def test_empty_sequences_are_filtered(monkeypatch, tmp_path):
    cand = tmp_path / "cand.fa"
    cand.write_text(">a\nACGT\n>b\n>c\nTTTT\n")

    def fake_run_quiet(cmd, check=True, cwd=None):  # pragma: no cover - behavior verified via assertions
        query = Path(cmd[-1])
        content = query.read_text()
        assert ">a" in content
        assert ">c" in content
        assert ">b" not in content
        out_path = query.with_suffix(query.suffix + ".out")
        out_path.write_text("")

    monkeypatch.setattr("haplongliner.sv_mode.run_quiet", fake_run_quiet)
    assert _validate_presence(cand, min_length=0) == set()


def test_all_empty_sequences_skip_repeatmasker(monkeypatch, tmp_path):
    cand = tmp_path / "cand.fa"
    cand.write_text(">a\n>b\n")

    called = {"flag": False}

    def fake_run_quiet(cmd, check=True, cwd=None):  # pragma: no cover - should not be called
        called["flag"] = True

    monkeypatch.setattr("haplongliner.sv_mode.run_quiet", fake_run_quiet)
    assert _validate_presence(cand, min_length=0) == set()
    assert not called["flag"]


def test_repeatmasker_runs_in_output_dir(monkeypatch, tmp_path):
    cand = tmp_path / "cand.fa"
    cand.write_text(">a\nACGT\n")

    def fake_run_quiet(cmd, check=True, cwd=None):
        assert cwd == tmp_path
        query = Path(cmd[-1])
        out_path = query.with_suffix(query.suffix + ".out")
        out_path.write_text("")

    monkeypatch.setattr("haplongliner.sv_mode.run_quiet", fake_run_quiet)
    assert _validate_presence(cand, min_length=0) == set()

