from haplongliner.sv_mode import _create_lifted_reorg, _filter_reference_bed

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

def test_filter_reference_pattern_match(tmp_path):
    raw = tmp_path / "ref.bed"
    raw.write_text(
        "chr1\t0\t6000\tL1HS_3end\t6000\t+\n" "chr1\t0\t6000\tL1PA2\t6000\t+\n",
    )
    out = _filter_reference_bed(raw, 5000, "L1HS", "hg38")
    lines = out.read_text().strip().splitlines()
    assert len(lines) == 1
    assert lines[0].split("\t")[3] == "L1HS_3end"


def test_create_lifted_reorg_preserves_reference_status(tmp_path):
    lifted_bed = tmp_path / "lifted.bed"
    lifted_bed.write_text("scaf1\t30\t40\tchr1,0,10,+,site1\t0\t+\n")
    ref_bed = tmp_path / "ref.bed"
    ref_bed.write_text("chr1\t0\t10\tsite1\t10\t+\tabsent\n")
    out = tmp_path / "lifted_reorg.bed"

    _create_lifted_reorg(lifted_bed, ref_bed, out)

    fields = out.read_text().strip().split("\t")
    assert fields[:7] == ["chr1", "0", "10", "chr1:0-10", "10", "+", "scaf1,30,40,10,+"]
    assert fields[7] == "absent"
