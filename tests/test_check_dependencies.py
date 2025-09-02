import builtins
import sys
import pytest
import haplongliner.utils as utils


def test_reports_missing_repeatmasker(monkeypatch):
    def fake_which(tool):
        return None if tool == "RepeatMasker" else "/usr/bin/" + tool

    monkeypatch.setattr(utils.shutil, "which", fake_which)

    with pytest.raises(SystemExit) as excinfo:
        utils.check_dependencies()

    assert "RepeatMasker" in str(excinfo.value)


def test_reports_missing_h5py(monkeypatch):
    monkeypatch.setattr(utils.shutil, "which", lambda tool: "/usr/bin/" + tool)

    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "h5py":
            raise ImportError
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)

    with pytest.raises(SystemExit) as excinfo:
        utils.check_dependencies()

    assert "h5py" in str(excinfo.value)


def test_reports_missing_dfam_partition(monkeypatch, tmp_path):
    rm_dir = tmp_path / "RepeatMasker"
    (rm_dir / "Libraries" / "famdb").mkdir(parents=True)
    rm_bin = rm_dir / "RepeatMasker"
    rm_bin.write_text("")

    def fake_which(tool):
        return str(rm_bin) if tool == "RepeatMasker" else "/usr/bin/" + tool

    monkeypatch.setattr(utils.shutil, "which", fake_which)
    monkeypatch.setitem(sys.modules, "h5py", object())

    with pytest.raises(SystemExit) as excinfo:
        utils.check_dependencies()

    assert "partition 7" in str(excinfo.value)
