from haplongliner.module3_DB import _parse_coords
import pytest


def test_parse_coords_string():
    coords = _parse_coords("chr1:10-20")
    assert coords == [("chr1", 10, 20)]


def test_parse_coords_file(tmp_path):
    bed = tmp_path / "sites.bed"
    bed.write_text("chr2\t30\t40\nchr3\t50\t60\n")
    coords = _parse_coords(str(bed))
    assert coords == [("chr2", 30, 40), ("chr3", 50, 60)]


def test_parse_coords_bad():
    with pytest.raises(ValueError):
        _parse_coords("badcoord")
