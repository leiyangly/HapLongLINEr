from pathlib import Path
import importlib.util
import sys
import pathlib

spec = importlib.util.spec_from_file_location(
    "liftover_paf",
    pathlib.Path(__file__).resolve().parents[1] / "scripts" / "liftover_paf.py",
)
liftover_paf = importlib.util.module_from_spec(spec)
spec.loader.exec_module(liftover_paf)

def test_liftover_respects_strand(tmp_path, capsys):
    paf = tmp_path / "test.paf"
    bed = tmp_path / "test.bed"
    paf.write_text("q1\t1000\t100\t200\t+\tchr1\t2000\t1000\t1100\t100\t100\t60\ttp:A:P\tcg:Z:100M\n")
    bed.write_text("q1\t110\t130\tname1\t0\t+\nq1\t140\t150\tname2\t0\t-\n")
    liftover_paf.liftover_paf(str(paf), str(bed), 0, 0, 2.0, False)
    out = capsys.readouterr().out.strip().splitlines()
    assert out[0].endswith("\t+")
    assert out[1].endswith("\t-")
