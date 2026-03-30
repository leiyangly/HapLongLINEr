import os
from pathlib import Path


PACKAGE_ROOT = Path(__file__).resolve().parent
PACKAGE_DATA_DIR = PACKAGE_ROOT / "data"
REPO_ROOT_DATA_DIR = PACKAGE_ROOT.parent / "data"


def data_search_paths() -> list[Path]:
    """Return candidate directories searched for bundled data files."""

    paths: list[Path] = []
    override = os.environ.get("HAPLONGLINER_DATA_DIR")
    if override:
        paths.append(Path(override).expanduser())
    paths.append(PACKAGE_DATA_DIR)
    if REPO_ROOT_DATA_DIR != PACKAGE_DATA_DIR:
        paths.append(REPO_ROOT_DATA_DIR)

    seen: set[Path] = set()
    ordered: list[Path] = []
    for path in paths:
        if path in seen:
            continue
        seen.add(path)
        ordered.append(path)
    return ordered


def get_data_path(name: str | Path, *, required: bool = True) -> Path:
    """Return the resolved path for a packaged data asset."""

    candidate = Path(name).expanduser()
    if candidate.is_absolute():
        if candidate.exists() or not required:
            return candidate
        raise FileNotFoundError(f"Data file '{candidate}' not found")

    for base in data_search_paths():
        path = base / candidate
        if path.exists():
            return path

    if not required:
        return data_search_paths()[0] / candidate

    searched = ", ".join(str(p) for p in data_search_paths())
    raise FileNotFoundError(f"Data file '{candidate}' not found in: {searched}")


def get_cache_dir() -> Path:
    """Return a writable cache directory for downloaded references."""

    override = os.environ.get("HAPLONGLINER_CACHE_DIR")
    if override:
        cache_dir = Path(override).expanduser()
    else:
        xdg_cache = os.environ.get("XDG_CACHE_HOME")
        base = Path(xdg_cache).expanduser() if xdg_cache else Path.home() / ".cache"
        cache_dir = base / "haplongliner"
    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir
