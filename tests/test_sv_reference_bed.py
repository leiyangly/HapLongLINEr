from haplongliner.sv_module import _filter_reference_bed


def test_filter_reference_hprc(tmp_path):
    raw = tmp_path / "hprc.bed"
    raw.write_text(
        "chr1\t0\t6000\tname1\t6000\t+\n" "chr1\t0\t4000\tname2\t4000\t+\n"
    )
    out = _filter_reference_bed(raw, 5000, "L1", "hprc")
    lines = out.read_text().strip().splitlines()
    assert len(lines) == 1
    assert lines[0].split("\t")[3] == "name1"


def test_filter_reference_nonhprc(tmp_path):
    raw = tmp_path / "ref.bed"
    raw.write_text(
        "chr1\t0\t6000\tL1HS\t6000\t+\n" "chr1\t0\t6000\tSVA_E\t6000\t+\n"
    )
    out = _filter_reference_bed(raw, 5000, "L1", "hg38")
    lines = out.read_text().strip().splitlines()
    assert len(lines) == 1
    assert lines[0].split("\t")[3] == "L1HS"
