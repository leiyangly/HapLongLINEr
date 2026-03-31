import sys

import haplongliner.cli as cli


def test_cli_mm_repeatmasker_pa_default(monkeypatch):
    captured = {}

    monkeypatch.setattr(cli, "check_dependencies", lambda: None)

    def fake_run_mm_mode(*args, **kwargs):
        captured["kwargs"] = kwargs

    monkeypatch.setattr(cli, "run_mm_mode", fake_run_mm_mode)
    monkeypatch.setattr(
        sys,
        "argv",
        ["haplongliner", "mm", "-i", "in.fa", "-o", "outdir"],
    )

    cli.main()

    assert captured["kwargs"]["repeatmasker_pa"] == 4


def test_cli_mm_repeatmasker_pa_override(monkeypatch):
    captured = {}

    monkeypatch.setattr(cli, "check_dependencies", lambda: None)

    def fake_run_mm_mode(*args, **kwargs):
        captured["kwargs"] = kwargs

    monkeypatch.setattr(cli, "run_mm_mode", fake_run_mm_mode)
    monkeypatch.setattr(
        sys,
        "argv",
        ["haplongliner", "mm", "-i", "in.fa", "-o", "outdir", "-pa", "8"],
    )

    cli.main()

    assert captured["kwargs"]["repeatmasker_pa"] == 8


def test_cli_sv_repeatmasker_pa_default(monkeypatch):
    captured = {}

    monkeypatch.setattr(cli, "check_dependencies", lambda: None)

    def fake_run_sv_mode(*args, **kwargs):
        captured["kwargs"] = kwargs

    monkeypatch.setattr(cli, "run_sv_mode", fake_run_sv_mode)
    monkeypatch.setattr(
        sys,
        "argv",
        ["haplongliner", "sv", "-i", "in.fa", "-s", "calls.vcf", "-o", "outdir"],
    )

    cli.main()

    assert captured["kwargs"]["repeatmasker_pa"] == 4


def test_cli_sv_repeatmasker_pa_override(monkeypatch):
    captured = {}

    monkeypatch.setattr(cli, "check_dependencies", lambda: None)

    def fake_run_sv_mode(*args, **kwargs):
        captured["kwargs"] = kwargs

    monkeypatch.setattr(cli, "run_sv_mode", fake_run_sv_mode)
    monkeypatch.setattr(
        sys,
        "argv",
        ["haplongliner", "sv", "-i", "in.fa", "-s", "calls.vcf", "-o", "outdir", "-pa", "12"],
    )

    cli.main()

    assert captured["kwargs"]["repeatmasker_pa"] == 12


def test_cli_sv_haplotype_option(monkeypatch):
    captured = {}

    monkeypatch.setattr(cli, "check_dependencies", lambda: None)

    def fake_run_sv_mode(*args, **kwargs):
        captured["kwargs"] = kwargs

    monkeypatch.setattr(cli, "run_sv_mode", fake_run_sv_mode)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "haplongliner",
            "sv",
            "-i",
            "HG00410.LSK110.R9_hapdup_phased_2.fasta",
            "-s",
            "calls.vcf",
            "-o",
            "outdir",
            "--hap",
            "2",
        ],
    )

    cli.main()

    assert captured["kwargs"]["hap"] == "2"
