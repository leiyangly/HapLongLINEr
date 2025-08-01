from haplongliner.find_longest_orf import find_longest_orf
import tempfile, os


def test_find_longest_orf_basic(tmp_path):
    input_file = tmp_path / "input.blastp"
    output_file = tmp_path / "out.blastp"
    input_file.write_text(
        "\n".join([
            "011138A1,186_phaseblock_2,1087928,1093929,+,6,876,1889\tL1rpORF1p\t98.817\t338\t4\t0\t1\t338\t1\t338\t0.0\t680\t338\t338\tL1rpORF1p",
            "011138A1,186_phaseblock_2,1087928,1093929,+,26,1956,5780\tL1rpORF2p\t99.137\t1275\t11\t0\t1\t1275\t1\t1275\t0.0\t2602\t1275\t1275\tL1rpORF2p",
        ])
    )
    find_longest_orf(input_file, output_file)
    out_lines = [l.strip() for l in open(output_file) if l.strip()]
    assert len(out_lines) == 1
    fields = out_lines[0].split()
    assert fields[0].split(',')[5] == '6'
    assert fields[15].split(',')[5] == '26'
