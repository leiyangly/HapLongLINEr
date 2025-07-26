import shutil
import sys
import os
from pathlib import Path
from typing import Dict, List

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


def run_quiet(cmd, **kwargs):
    """Run a subprocess suppressing output from successful commands."""
    import subprocess

    kwargs.setdefault("stdout", subprocess.PIPE)
    kwargs.setdefault("stderr", subprocess.PIPE)
    result = subprocess.run(cmd, **kwargs)

    def _to_str(stream):
        if stream is None:
            return ""
        if isinstance(stream, bytes):
            return stream.decode()
        return str(stream)

    if result.returncode != 0:
        if result.stderr:
            sys.stderr.write(_to_str(result.stderr))
        if kwargs.get("stdout") == subprocess.PIPE and result.stdout:
            sys.stderr.write(_to_str(result.stdout))
        result.check_returncode()
    return result


def shift_fasta_start(fa_path: Path, delta: int = -1) -> None:
    """Adjust start coordinate in FASTA headers by ``delta``.

    Headers must follow ``chrom;start;end;...``. Lines that do not conform are
    left unchanged.
    """
    tmp = fa_path.with_suffix(".tmp")
    with open(fa_path) as fin, open(tmp, "w") as out:
        for line in fin:
            if line.startswith(">"):
                parts = line[1:].strip().split(";")
                if len(parts) >= 2:
                    try:
                        parts[1] = str(int(parts[1]) + delta)
                        line = ">" + ";".join(parts) + "\n"
                    except ValueError:
                        pass
            out.write(line)
    os.replace(tmp, fa_path)


def verify_fasta_file(path: str) -> None:
    """Exit if ``path`` is not a readable FASTA file."""
    import gzip

    fpath = Path(path)
    if not fpath.exists():
        sys.exit(f"Error: FASTA file '{path}' not found.")

    opener = gzip.open if str(fpath).endswith(".gz") else open
    try:
        with opener(fpath, "rt") as fh:
            first = fh.readline()
        if not first.startswith(">"):
            raise ValueError("Missing FASTA header")
    except Exception as exc:
        sys.exit(f"Error: FASTA file '{path}' appears malformed: {exc}")


def _check_bed_line(line: str) -> None:
    fields = line.rstrip().split()
    if len(fields) < 3:
        raise ValueError("fewer than 3 columns")
    int(fields[1])
    int(fields[2])


def verify_bed_file(path: str) -> None:
    """Exit if ``path`` is not a readable BED-like file."""
    import gzip

    bpath = Path(path)
    if not bpath.exists():
        sys.exit(f"Error: file '{path}' not found.")

    opener = gzip.open if str(bpath).endswith(".gz") else open
    try:
        with opener(bpath, "rt") as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                _check_bed_line(line)
                break
            else:
                raise ValueError("no data lines found")
    except Exception as exc:
        sys.exit(f"Error: BED file '{path}' appears malformed: {exc}")


def verify_sv_file(path: str) -> None:
    """Exit if ``path`` is not a readable SV/VCF/BED file."""
    import gzip

    spath = Path(path)
    if not spath.exists():
        sys.exit(f"Error: file '{path}' not found.")

    opener = gzip.open if str(spath).endswith(".gz") else open
    with opener(spath, "rt") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            if len(line.rstrip().split("\t")) >= 3:
                return
        sys.exit(f"Error: SV file '{path}' appears malformed or empty.")


def read_paf(path: str | Path) -> Dict[str, List[str]]:
    """Return mapping ``query_name -> fields`` from a minimap2 PAF.

    Alignments shorter than 200 bp are ignored and, if multiple alignments are
    present for a query, the longest one is retained.
    """

    hits: Dict[str, List[str]] = {}
    lengths: Dict[str, int] = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            fields = line.rstrip().split('\t')
            if len(fields) < 4:
                continue
            qname = fields[0]
            try:
                aln_len = int(fields[3]) - int(fields[2])
            except ValueError:
                continue
            if aln_len < 200:
                continue
            prev_len = lengths.get(qname, -1)
            if aln_len > prev_len:
                hits[qname] = fields
                lengths[qname] = aln_len

    return hits


def verify_repeatmasker_file(path: str) -> None:
    """Exit if ``path`` is a readable RepeatMasker BED or .out file."""
    import gzip
    import re

    rpath = Path(path)
    if not rpath.exists():
        sys.exit(f"Error: file '{path}' not found.")

    opener = gzip.open if str(rpath).endswith(".gz") else open
    try:
        with opener(rpath, "rt") as fh:
            lines = []
            for _ in range(10):
                line = fh.readline()
                if not line:
                    break
                if line.startswith(("#", "track", "browser")):
                    continue
                lines.append(line)
            if not lines:
                raise ValueError("no data lines found")
            if any("SW" in l and "perc" in l for l in lines[:4]):
                lines = lines[4:]
                if not lines:
                    raise ValueError("header only")
            line0 = lines[0].strip()
            fields = re.split(r"\s+", line0)
            is_out = (
                len(fields) >= 14
                and fields[0].replace(".", "", 1).isdigit()
                and fields[5].isdigit()
                and fields[6].isdigit()
            )
            if is_out:
                int(fields[5])
                int(fields[6])
            elif len(fields) >= 3:
                int(fields[1])
                int(fields[2])
            else:
                raise ValueError("unrecognized format")
    except Exception as exc:
        sys.exit(f"Error: RepeatMasker file '{path}' appears malformed: {exc}")
