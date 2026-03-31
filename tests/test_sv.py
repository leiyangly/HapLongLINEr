import gzip
from pathlib import Path

import pytest

from haplongliner.sv_mode import (
    _parse_sv,
    _resolve_sv_haplotype_selection,
    run_sv_mode,
)


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


def test_parse_sv_filters_phased_haplotype_vcf(tmp_path):
    vcf = tmp_path / "calls.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample\n"
        "chr1\t11\tdel1\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=21\tGT\t1|0\n"
        "chr1\t31\tdel2\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=41\tGT\t0|1\n"
        "chr1\t51\tdel3\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=61\tGT\t1|1\n"
        "chr1\t71\tdel4\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=81\tGT\t0|0\n"
        "chr1\t91\tins1\tA\tACGT\t.\tPASS\tSVTYPE=INS;END=91\tGT\t0/1\n"
    )

    del1, ins1 = _parse_sv(vcf, hap="1")
    del2, ins2 = _parse_sv(vcf, hap="2")
    delb, insb = _parse_sv(vcf, hap="both")

    assert del1 == [("chr1", 10, 20), ("chr1", 50, 60)]
    assert del2 == [("chr1", 30, 40), ("chr1", 50, 60)]
    assert delb == [("chr1", 10, 20), ("chr1", 30, 40), ("chr1", 50, 60)]
    assert ins1 == [("chr1", 90, 91, "ACGT")]
    assert ins2 == [("chr1", 90, 91, "ACGT")]
    assert insb == [("chr1", 90, 91, "ACGT")]


def test_resolve_sv_haplotype_selection_auto_from_input_name():
    resolved, note = _resolve_sv_haplotype_selection(
        "auto",
        "/tmp/HG00410.LSK110.R9_hapdup_phased_1.fasta",
    )

    assert resolved == "1"
    assert "phased_1" in note


def test_resolve_sv_haplotype_selection_auto_falls_back_to_both():
    resolved, note = _resolve_sv_haplotype_selection("auto", "/tmp/sample.fa")

    assert resolved == "both"
    assert "fallback" in note
