import gzip
from pathlib import Path

from haplongliner.sv_module import _parse_sv


def test_parse_sv_handles_gz(tmp_path):
    sv_gz = Path('tests/test.genome.HG00410.vcf.gz')
    deletions, insertions = _parse_sv(sv_gz)
    assert deletions[0] == ('chr1', 10828, 10858)
    assert insertions[0][:3] == ('chr1', 66365, 66366)
    assert insertions[0][3].startswith('TATTATTATATA')
    assert len(deletions) == 20629
    assert len(insertions) == 23382

    # also ensure plain VCF works
    plain = tmp_path / 'test.vcf'
    with gzip.open(sv_gz, 'rb') as src, open(plain, 'wb') as out:
        out.write(src.read())
    del2, ins2 = _parse_sv(plain)
    assert del2 == deletions
    assert ins2 == insertions

