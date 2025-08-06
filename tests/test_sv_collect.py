from pathlib import Path
from haplongliner.sv_module import _collect_long_insertions


def test_collect_long_insertions(tmp_path):
    inter = tmp_path / "inter.bed"
    inter.write_text("chr1\t0\t10\tL1\t6000\t+\tinfo\tchr1\t100\t110\tAAAA\n")
    insertions = [
        ("chr1", 100, 110, "A" * 6000),
        ("chr1", 200, 210, "T" * 6000),
    ]
    extras = _collect_long_insertions(insertions, inter, 5000)
    assert extras == [("chr1", 200, 210, "T" * 6000, "INS_chr1_200")]
