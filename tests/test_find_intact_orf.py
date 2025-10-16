import importlib.util
from pathlib import Path


def _load_find_intact_orf():
    module_path = Path(__file__).resolve().parents[1] / "haplongliner" / "find_intact_orf.py"
    spec = importlib.util.spec_from_file_location("_find_intact_orf", module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.find_intact_orf


find_intact_orf = _load_find_intact_orf()


def run_find_intact(tmp_path: Path, contents: str) -> Path:
    infile = tmp_path / "cand_orf_combine.blastp"
    outfile = tmp_path / "cand_orf_intact.blastp"
    infile.write_text(contents)
    find_intact_orf(str(infile), str(outfile))
    return outfile


def test_marks_reverse_subject_as_intact(tmp_path: Path) -> None:
    row = (
        "cand1,extra,fields\tL1rpORF1p\t99.0\t338\t0\t0\t1\t338\t338\t1\t0.0\t600\t338\t338\tL1rpORF1p\t"
        "cand1,extra,fields\tL1rpORF2p\t99.0\t1275\t0\t0\t1\t1275\t1275\t1\t0.0\t2000\t1275\t1275\tL1rpORF2p\n"
    )
    outfile = run_find_intact(tmp_path, row)
    assert outfile.read_text() == row


def test_skips_truncated_orf(tmp_path: Path) -> None:
    truncated = (
        "cand2\tL1rpORF1p\t99.0\t300\t0\t0\t1\t300\t1\t300\t0.0\t500\t338\t338\tL1rpORF1p\t"
        "cand2\tL1rpORF2p\t99.0\t1000\t0\t0\t1\t1000\t1\t1000\t0.0\t1500\t1275\t1275\tL1rpORF2p\n"
    )
    outfile = run_find_intact(tmp_path, truncated)
    assert outfile.read_text() == ""
