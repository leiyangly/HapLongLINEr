import shutil
import sys
from pathlib import Path

def check_dependencies():
    """Ensure required external tools are available."""
    tools = ["seqtk", "minimap2", "getorf", "blastp"]
    missing = [tool for tool in tools if shutil.which(tool) is None]
    if missing:
        sys.exit(
            f"Error: The following required tools are missing from your PATH: {', '.join(missing)}"
        )


def verify_blast_db(db_prefix):
    """Ensure the BLAST database exists and is not a Git LFS placeholder."""
    db_path = Path(db_prefix)
    if not db_path.exists():
        sys.exit(f"Error: BLAST database '{db_prefix}' not found.")

    if db_path.stat().st_size < 200:
        with open(db_path) as fh:
            first_line = fh.readline()
        if first_line.startswith("version https://git-lfs.github.com"):
            sys.exit(
                "Error: BLAST database appears to be a Git LFS placeholder. "
                "Run 'git lfs pull' to download the required data files."
            )


