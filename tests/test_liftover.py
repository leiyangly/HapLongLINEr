from pathlib import Path

from haplongliner.sv_module import _liftover_l1s, _read_paf


def write_paf(path: Path, lines: list[str]) -> None:
    path.write_text("\n".join(lines) + "\n")


def test_read_paf_filters_and_selects_longest(tmp_path):
    paf = tmp_path / "test.paf"
    lines = [
        "\t".join(
            [
                "q1",
                "1000",
                "0",
                "180",
                "+",
                "chr1",
                "10000",
                "0",
                "180",
                "180",
                "180",
                "60",
            ]
        ),
        "\t".join(
            [
                "q1",
                "1000",
                "0",
                "400",
                "+",
                "chr1",
                "10000",
                "200",
                "600",
                "400",
                "400",
                "60",
            ]
        ),
        "\t".join(
            [
                "q2",
                "1000",
                "0",
                "150",
                "+",
                "chr1",
                "10000",
                "0",
                "150",
                "150",
                "150",
                "60",
            ]
        ),
    ]
    write_paf(paf, lines)
    hits = _read_paf(paf)
    assert "q1" in hits
    assert hits["q1"][3] == "400"
    assert "q2" not in hits


def test_liftover_l1s_plus_and_minus(tmp_path):
    paf = tmp_path / "aln.paf"
    ref_bed = tmp_path / "ref.bed"
    master = tmp_path / "master.bed"

    paf_lines = [
        "\t".join(
            [
                "chrA",
                "1000",
                "0",
                "1000",
                "+",
                "chrX",
                "10000",
                "1200",
                "2200",
                "1000",
                "1000",
                "60",
                "tp:A:P",
                "cg:Z:1000M",
            ]
        ),
        "\t".join(
            [
                "chrB",
                "1000",
                "0",
                "100",
                "+",
                "chrY",
                "10000",
                "300",
                "400",
                "100",
                "100",
                "60",
                "tp:A:P",
                "cg:Z:100M",
            ]
        ),
    ]

    write_paf(paf, paf_lines)

    ref_bed.write_text("chrA\t100\t900\tL1a\t+\n" + "chrB\t0\t100\tL1b\t-\n")

    master.write_text("chrA\t100\t1000\tL1a\t+\n" + "chrB\t0\t100\tL1b\t-\n")

    lifted = _liftover_l1s(paf, ref_bed, min_length=50, master_files=[master])
    lifted.sort(key=lambda x: x[3])

    assert lifted[0] == (
        "chrA",
        100,
        1000,
        "L1a",
        900,
        "+",
        "chrX,1300,2100,+",
        "chrX,1300,2100,+",
        1300,
        2100,
    )
    assert lifted[1] == (
        "chrB",
        0,
        100,
        "L1b",
        100,
        "-",
        "chrY,300,400,+",
        "chrY,300,400,+",
        300,
        400,
    )


def test_liftover_respects_min_length(tmp_path):
    paf = tmp_path / "aln.paf"
    ref_bed = tmp_path / "ref.bed"
    master = tmp_path / "master.bed"

    paf.write_text(
        "\t".join(
            [
                "chrC",
                "1000",
                "0",
                "80",
                "+",
                "scaf1",
                "10000",
                "2000",
                "2080",
                "80",
                "80",
                "60",
                "tp:A:P",
                "cg:Z:80M",
            ]
        )
        + "\n"
    )

    ref_bed.write_text("chrC\t100\t180\tL1c\t+\n")
    master.write_text("chrC\t100\t180\tL1c\t+\n")

    lifted = _liftover_l1s(paf, ref_bed, min_length=100, master_files=[master])
    assert lifted == []
