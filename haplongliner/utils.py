import shutil
import sys


def check_dependencies():
    """Ensure required external tools are available."""
    tools = ["seqtk", "minimap2", "getorf", "blastp"]
    missing = [tool for tool in tools if shutil.which(tool) is None]
    if missing:
        sys.exit(
            f"Error: The following required tools are missing from your PATH: {', '.join(missing)}"
        )
