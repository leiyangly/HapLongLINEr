import gzip
from pathlib import Path

import pytest

from haplongliner.sv_mode import _parse_sv, run_sv_mode


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


def test_reuse_existing_alignment(monkeypatch, tmp_path):
    inp = tmp_path / "in.fa"
    inp.write_text(">chr1\nACGT\n")
    sv = tmp_path / "sv.bed"
    sv.write_text("chr1\t0\t1\n")
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n")

    outdir = tmp_path / "out"
    outdir.mkdir()
    paf = outdir / "genome_alignment.paf"
    paf.write_text("paf\n")

    called = {"flag": False}

    def fake_run_quiet(cmd, check=True, cwd=None, **kwargs):  # pragma: no cover - should not be called
        called["flag"] = True

    monkeypatch.setattr("haplongliner.sv_mode.run_quiet", fake_run_quiet)
    def fake_liftover_paf(*args, **kwargs):  # pragma: no cover - stops execution
        raise RuntimeError("stop")

    monkeypatch.setattr("haplongliner.sv_mode.liftover_paf", fake_liftover_paf)

    with pytest.raises(RuntimeError, match="stop"):
        run_sv_mode(str(inp), str(sv), str(ref), str(outdir), exist="yes")

    assert not called["flag"]

