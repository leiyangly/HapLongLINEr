import gzip
from pathlib import Path
import pytest

from haplongliner.utils import verify_fasta_file, verify_repeatmasker_file


def test_verify_fasta_file_accepts_gz(tmp_path):
    src = Path('tests/test.genome.HG00410.1.fa')
    gz = tmp_path / 'test.fa.gz'
    with open(src, 'rb') as f_in, gzip.open(gz, 'wb') as f_out:
        f_out.write(f_in.read())
    verify_fasta_file(str(src))
    verify_fasta_file(str(gz))


def test_verify_repeatmasker_file_accepts_out_and_bed(tmp_path):
    out = Path('tests/test.genome.HG00410.1.out')
    gz = tmp_path / 'test.out.gz'
    with open(out, 'rb') as f_in, gzip.open(gz, 'wb') as f_out:
        f_out.write(f_in.read())
    verify_repeatmasker_file(str(out))
    verify_repeatmasker_file(str(gz))
    bed = tmp_path / 'test.bed'
    bed.write_text('chr1\t1\t2\tname\t1\t+\n')
    verify_repeatmasker_file(str(bed))


def test_verify_repeatmasker_file_bad(tmp_path):
    bad = tmp_path / 'bad.txt'
    bad.write_text('malformed line')
    with pytest.raises(SystemExit):
        verify_repeatmasker_file(str(bad))
