import pytest
import haplongliner.utils as utils


def test_reports_missing_blastn(monkeypatch):
    def fake_which(tool):
        return None if tool == "blastn" else "/usr/bin/" + tool

    monkeypatch.setattr(utils.shutil, "which", fake_which)

    with pytest.raises(SystemExit) as excinfo:
        utils.check_dependencies()

    assert "blastn" in str(excinfo.value)
