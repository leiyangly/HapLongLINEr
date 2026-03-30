from pathlib import Path

from haplongliner.resources import get_cache_dir, get_data_path


def test_get_data_path_works_outside_repo_cwd(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    path = get_data_path("L1rp.fa")
    assert path.exists()
    assert path.name == "L1rp.fa"


def test_get_data_path_prefers_env_override(monkeypatch, tmp_path):
    override = tmp_path / "custom-data"
    override.mkdir()
    (override / "L1rp.fa").write_text(">ref\nACGT\n")
    monkeypatch.setenv("HAPLONGLINER_DATA_DIR", str(override))

    assert get_data_path("L1rp.fa") == override / "L1rp.fa"


def test_get_cache_dir_uses_override(monkeypatch, tmp_path):
    cache_dir = tmp_path / "cache"
    monkeypatch.setenv("HAPLONGLINER_CACHE_DIR", str(cache_dir))

    assert get_cache_dir() == cache_dir
    assert cache_dir.exists()


def test_get_data_path_raises_for_missing_file():
    missing = Path("definitely_missing_resource.txt")
    try:
        get_data_path(missing)
    except FileNotFoundError as exc:
        assert str(missing) in str(exc)
    else:
        raise AssertionError("Expected FileNotFoundError for missing resource")
