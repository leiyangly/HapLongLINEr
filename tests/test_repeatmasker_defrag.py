from haplongliner.repeatmasker import read_collapsed_repeatmasker_out
from haplongliner.rm_mode import parse_repeatmasker


def test_repeatmasker_defrag_merges_fragmented_young_l1(tmp_path):
    out = tmp_path / "cand.out"
    out.write_text(
        "19083 1.2 0.0 0.0 scaf1,100,6096,+ 1 2127 (3869) + L1HS LINE/L1 10 2136 (3896) 1\n"
        "28563 1.6 0.0 0.0 scaf1,100,6096,+ 1978 5996 (0) + L1PA2 LINE/L1 2110 6128 (27) 2\n"
    )

    loci = read_collapsed_repeatmasker_out(out)

    assert len(loci) == 1
    hit = loci[0]
    assert hit.query_name == "scaf1,100,6096,+"
    assert (hit.query_start, hit.query_end) == (0, 5996)
    assert hit.family == "L1PA2"
    assert hit.strand == "+"
    assert hit.covered_bp == 5996


def test_repeatmasker_defrag_does_not_merge_tandem_full_length_l1s(tmp_path):
    out = tmp_path / "cand.out"
    out.write_text(
        "27925 1.4 0.0 0.0 scaf1 1 6000 (6000) + L1HS LINE/L1 124 6155 (0) 1\n"
        "27925 1.4 0.0 0.0 scaf1 6001 12000 (0) + L1HS LINE/L1 124 6155 (0) 2\n"
    )

    loci = read_collapsed_repeatmasker_out(out)

    assert len(loci) == 2


def test_rm_parse_repeatmasker_defragments_out_input(tmp_path):
    rm_out = tmp_path / "genome.fa.out"
    rm_out.write_text(
        "19083 1.2 0.0 0.0 chr1 101 2227 (3869) + L1HS LINE/L1 10 2136 (3896) 1\n"
        "28563 1.6 0.0 0.0 chr1 2078 6096 (0) + L1PA2 LINE/L1 2110 6128 (27) 2\n"
        "500 5.0 0.0 0.0 chr1 9000 9100 (0) + AluY SINE/Alu 1 100 (0) 3\n"
    )
    bed = tmp_path / "parsed.bed"

    parse_repeatmasker(rm_out, bed)

    assert bed.read_text().strip().splitlines() == [
        "chr1\t100\t6096\tL1PA2\t5996\t+",
        "chr1\t8999\t9100\tAluY\t101\t+",
    ]
