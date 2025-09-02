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
    # The query identifiers include ORF-specific fields at the end.  The index
    # of the ORF is stored three positions from the end, regardless of whether
    # a scaffold name is present.
    assert fields[0].split(',')[-3] == '6'
    assert fields[15].split(',')[-3] == '26'


def test_find_longest_orf_groups_by_candidate(tmp_path):
    """Multiple ORFs for the same candidate should be collapsed."""
    inp = tmp_path / "input.blastp"
    outp = tmp_path / "out.blastp"
    inp.write_text(
        "\n".join(
            [
                "scaf,1,100,+,1,5,50\tL1rpORF1p\t0\t50",
                "scaf,1,100,+,2,10,60\tL1rpORF1p\t0\t100",
                "scaf,1,100,+,3,20,80\tL1rpORF2p\t0\t150",
            ]
        )
    )
    find_longest_orf(inp, outp)
    lines = [l.strip() for l in open(outp) if l.strip()]
    assert len(lines) == 1
    f = lines[0].split()
    assert f[0].split(',')[-3] == '2'
    assert f[1] == 'L1rpORF1p'
    assert f[4].split(',')[-3] == '3'
    assert f[5] == 'L1rpORF2p'


def test_find_longest_orf_sv_headers(tmp_path):
    """Headers with fewer comma-separated fields should still be grouped."""
    inp = tmp_path / "input.blastp"
    outp = tmp_path / "out.blastp"
    base = "22822069-22828100_122_phaseblock_26_660301_666332_-"
    q1 = f"{base},1,906,1919"
    q2 = f"{base},2,1986,5810"
    inp.write_text("\n".join([
        f"{q1}\tL1rpORF1p\t0\t338",
        f"{q2}\tL1rpORF2p\t0\t1275",
    ]))
    find_longest_orf(inp, outp)
    lines = [l.strip() for l in open(outp) if l.strip()]
    assert len(lines) == 1
    fields = lines[0].split('\t')
    assert fields[0] == q1
    assert fields[4] == q2

