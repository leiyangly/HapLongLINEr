from pathlib import Path

import pytest

from haplongliner.sv_mode import _classify_sv


@pytest.fixture()
def lifted_entry(tmp_path: Path) -> tuple[list[tuple], Path]:
    lifted = [
        (
            "chr1",
            0,
            100,
            "chr1:0-100",
            100,
            "+",
            "scaf1,0,100,100,+",
            "scaf1,0,100,100,+",
            0,
            100,
        )
    ]
    lifted_bed = tmp_path / "lifted_reorg.bed"
    lifted_bed.write_text(
        "chr1\t0\t100\tchr1:0-100\t100\t+\tscaf1,0,100,100,+\n"
    )
    return lifted, lifted_bed


@pytest.fixture(autouse=True)
def fake_bedtools(monkeypatch):
    def _fake_intersect(a: Path, b: Path, output: Path) -> None:
        with open(a) as fa, open(b) as fb, open(output, "w") as out:
            a_lines = [line.strip() for line in fa if line.strip()]
            b_lines = [line.strip() for line in fb if line.strip()]
            if not a_lines or not b_lines:
                return
            for a_line in a_lines:
                for b_line in b_lines:
                    out.write(f"{a_line}\t{b_line}\n")

    monkeypatch.setattr("haplongliner.sv_mode._bedtools_intersect", _fake_intersect)


def test_classify_sv_marks_absent_only_with_reciprocal_coverage(
    tmp_path: Path, lifted_entry: tuple[list[tuple], Path]
) -> None:
    lifted, lifted_bed = lifted_entry

    status, _ = _classify_sv(
        lifted,
        [("chr1", 0, 100)],
        [],
        tmp_path,
        lifted_bed,
    )
    assert status["chr1:0-100"] == "absent"

    status, _ = _classify_sv(
        lifted,
        [("chr1", 0, 120)],
        [],
        tmp_path,
        lifted_bed,
    )
    assert status["chr1:0-100"] == "present"

    status, _ = _classify_sv(
        lifted,
        [("chr1", 5, 95)],
        [],
        tmp_path,
        lifted_bed,
    )
    assert status["chr1:0-100"] == "present"
