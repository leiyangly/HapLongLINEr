from pathlib import Path

import pytest

from haplongliner.sv_mode import _classify_sv


def _make_lifted(
    tmp_path: Path,
    *,
    ref_len: int = 100,
    lifted_len: int = 100,
    info_suffix: str = "+",
    ref_status: str = "disrupted",
) -> tuple[list[tuple], Path]:
    lifted = [
        (
            "chr1",
            0,
            ref_len,
            "chr1:0-100",
            ref_len,
            "+",
            f"scaf1,0,{lifted_len},{lifted_len},{info_suffix}",
            f"scaf1,0,{lifted_len},{lifted_len},{info_suffix}",
            0,
            lifted_len,
        )
    ]
    lifted_bed = tmp_path / "lifted_reorg.bed"
    lifted_bed.write_text(
        f"chr1\t0\t{ref_len}\tchr1:0-100\t{ref_len}\t+\t"
        f"scaf1,0,{lifted_len},{lifted_len},{info_suffix}\t{ref_status}\n"
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


def test_classify_sv_ref_present_plus_del_marks_absent(tmp_path: Path) -> None:
    lifted, lifted_bed = _make_lifted(tmp_path, ref_status="disrupted")

    status, _, _ = _classify_sv(
        lifted,
        [("chr1", 0, 100)],
        [],
        tmp_path,
        lifted_bed,
    )
    assert status["chr1:0-100"] == "absent"


def test_classify_sv_ref_absent_without_ins_stays_absent(tmp_path: Path) -> None:
    lifted, lifted_bed = _make_lifted(tmp_path, ref_status="absent")

    status, _, _ = _classify_sv(
        lifted,
        [],
        [],
        tmp_path,
        lifted_bed,
    )
    assert status["chr1:0-100"] == "absent"


def test_classify_sv_ref_absent_plus_ins_marks_disrupted(tmp_path: Path) -> None:
    lifted, lifted_bed = _make_lifted(tmp_path, ref_status="absent")

    status, _, _ = _classify_sv(
        lifted,
        [],
        [("chr1", 0, 1, "A" * 50)],
        tmp_path,
        lifted_bed,
    )
    assert status["chr1:0-100"] == "disrupted"
