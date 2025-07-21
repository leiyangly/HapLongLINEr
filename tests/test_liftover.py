from pathlib import Path

from haplongliner.module2_SV import _liftover_l1s


def write_paf(path: Path, lines: list[str]) -> None:
    path.write_text("\n".join(lines) + "\n")


def test_liftover_l1s_plus_and_minus(tmp_path):
    minus_paf = tmp_path / "minus.paf"
    plus_paf = tmp_path / "plus.paf"
    ref_bed = tmp_path / "ref.bed"

    minus_lines = [
        "\t".join([
            "L1a_-2kb",
            "1000",
            "0",
            "1000",
            "+",
            "chrX",
            "10000",
            "1000",
            "1200",
            "100",
            "200",
            "60",
        ]),
        "\t".join([
            "L1b_-2kb",
            "1000",
            "0",
            "1000",
            "-",
            "chrY",
            "10000",
            "500",
            "700",
            "100",
            "200",
            "60",
        ]),
    ]

    plus_lines = [
        "\t".join([
            "L1a_+2kb",
            "1000",
            "0",
            "1000",
            "+",
            "chrX",
            "10000",
            "1300",
            "1500",
            "100",
            "200",
            "60",
        ]),
        "\t".join([
            "L1b_+2kb",
            "1000",
            "0",
            "1000",
            "-",
            "chrY",
            "10000",
            "300",
            "480",
            "100",
            "200",
            "60",
        ]),
    ]

    write_paf(minus_paf, minus_lines)
    write_paf(plus_paf, plus_lines)

    ref_bed.write_text(
        "chrX\t0\t900\tL1a\t+\n" +
        "chrY\t0\t100\tL1b\t-\n"
    )

    lifted = _liftover_l1s(minus_paf, plus_paf, ref_bed, min_length=50)
    lifted.sort(key=lambda x: x[3])

    assert lifted[0] == ("chrX", 1200, 1300, "L1a", 900, "+")
    assert lifted[1] == ("chrY", 480, 500, "L1b", 100, "-")
