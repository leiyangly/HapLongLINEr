from haplongliner.extract_l1 import extract_te_from_bed


def test_extract_te_pattern_matching(tmp_path):
    bed = tmp_path / "input.bed"
    bed.write_text(
        "chr1\t0\t6000\tL1HS\t6000\t+\n"
        "chr1\t0\t6000\tL1HS_3end\t6000\t+\n"
        "chr1\t0\t6000\tL1PA3\t6000\t+\n"
    )
    out = tmp_path / "out.bed"
    extract_te_from_bed(str(bed), str(out), 5000, "L1HS")
    lines = out.read_text().strip().splitlines()
    names = [line.split("\t")[3] for line in lines]
    assert names == ["L1HS", "L1HS_3end"]
