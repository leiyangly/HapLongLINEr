from pathlib import Path

from haplongliner.sv_mode import _validate_orfs, _append_liftover_orfs


def test_validate_orfs_outputs_only_orfs(monkeypatch, tmp_path):
    cand = tmp_path / "cand.fa"
    cand.write_text(">L1,chr1,0,24,+\nATGGCCATTGTAATGGGCCGCTGAA\n")

    def fake_run_quiet(cmd, check=True, cwd=None, **kwargs):
        if cmd[0] == "getorf":
            out = Path(cmd[cmd.index("-outseq") + 1])
            out.write_text("")
        elif cmd[0] == "blastp":
            out = Path(cmd[cmd.index("-out") + 1])
            out.write_text("")

    monkeypatch.setattr("haplongliner.sv_mode.run_quiet", fake_run_quiet)
    monkeypatch.setattr("haplongliner.sv_mode.verify_blast_db", lambda x: None)
    monkeypatch.setattr(
        "haplongliner.sv_mode.find_longest_orf", lambda inp, out: Path(out).write_text("")
    )
    monkeypatch.setattr(
        "haplongliner.sv_mode.find_intact_orf", lambda inp, out: Path(out).write_text("")
    )

    _validate_orfs(cand)
    orf_fa = cand.with_name("cand_orf.fa")
    assert orf_fa.read_text() == ""


def test_append_liftover_orfs(monkeypatch, tmp_path):
    cand = tmp_path / "cand.fa"
    cand.write_text(">INS_chr1_50,chr1,50,70,.\nATGGCC\n")
    sv_fa = tmp_path / "haplongliner_sv.fa"
    sv_fa.write_text(
        f">chr1,0,810,810,+,ALN\n{'ATG' * 270}\n"
        f">chr1,50,830,780,+,INS\n{'ATG' * 270}\n"
    )

    class Dummy:
        def __init__(self, stdout=""):
            self.stdout = stdout

    def fake_run_quiet(cmd, check=True, cwd=None, text=False, **kwargs):
        if cmd[0] == "getorf":
            outseq = cmd[cmd.index("-outseq") + 1]
            if outseq == "stdout":
                prot = "M" * 270
                stdout = (
                    f">chr1_0_810_810_+_ALN_1 [1 - 810]\n{prot}\n"
                    f">chr1_50_830_780_+_INS_1 [1 - 780]\n{prot}\n"
                )
                return Dummy(stdout=stdout)
            else:
                Path(outseq).write_text("")
                return Dummy()
        elif cmd[0] == "blastp":
            out = Path(cmd[cmd.index("-out") + 1])
            out.write_text("")
            return Dummy()

    monkeypatch.setattr("haplongliner.sv_mode.run_quiet", fake_run_quiet)
    monkeypatch.setattr("haplongliner.sv_mode.verify_blast_db", lambda x: None)
    monkeypatch.setattr(
        "haplongliner.sv_mode.find_longest_orf", lambda inp, out: Path(out).write_text("")
    )
    monkeypatch.setattr(
        "haplongliner.sv_mode.find_intact_orf", lambda inp, out: Path(out).write_text("")
    )

    _validate_orfs(cand)
    orf_fa = cand.with_name("cand_orf.fa")
    _append_liftover_orfs(orf_fa, sv_fa)

    content = orf_fa.read_text()
    assert ">INS_chr1_50,chr1,50,70,.,1,1,6" not in content
    assert ">chr1,0,810,810,+,ALN,1,1,810" in content
    assert ">chr1,50,830,780,+,INS,1,1,780" not in content


def test_cand_orf_intact_pipeline(monkeypatch, tmp_path):
    cand = tmp_path / "cand.fa"
    cand.write_text(">L1,chr1,0,24,+\nATGGCCATTGTAATGGGCCGCTGAA\n")

    def fake_run_quiet(cmd, check=True, cwd=None, input=None, text=False, **kwargs):
        if cmd[0] == "getorf":
            out = Path(cmd[cmd.index("-outseq") + 1])
            out.write_text(">L1,chr1,0,24,+,1,1,24\nATG\n")
        elif cmd[0] == "blastp":
            out = Path(cmd[cmd.index("-out") + 1])
            out.write_text("blast\n")

    monkeypatch.setattr("haplongliner.sv_mode.run_quiet", fake_run_quiet)
    monkeypatch.setattr("haplongliner.sv_mode.verify_blast_db", lambda x: None)

    calls = []

    def fake_find_longest(inp, out):
        calls.append(("longest", inp, out))
        Path(out).write_text("longest\n")

    def fake_find_intact(inp, out):
        calls.append(("intact", inp, out))
        Path(out).write_text("intact\n")

    monkeypatch.setattr("haplongliner.sv_mode.find_longest_orf", fake_find_longest)
    monkeypatch.setattr("haplongliner.sv_mode.find_intact_orf", fake_find_intact)

    _validate_orfs(cand)

    assert calls[0][0] == "longest"
    assert calls[1][0] == "intact"
    assert calls[1][1] == calls[0][2]
    intact_file = cand.with_name("cand_orf_intact.blastp")
    assert intact_file.read_text() == "intact\n"


def test_validate_orfs_sv_headers(monkeypatch, tmp_path):
    cand = tmp_path / "cand.fa"
    cand.write_text(">dummy\nATG\n")
    base = "22822069-22828100_122_phaseblock_26_660301_666332_-"
    q1 = f"{base}_1_906_1919"
    q2 = f"{base}_2_1986_5810"

    def fake_run_quiet(cmd, check=True, cwd=None, input=None, text=False, **kwargs):
        if cmd[0] == "getorf":
            out = Path(cmd[cmd.index("-outseq") + 1])
            out.write_text(f">{q1}\nM\n>{q2}\nM\n")
        elif cmd[0] == "blastp":
            out = Path(cmd[cmd.index("-out") + 1])
            out.write_text(
                f"{q1}\tL1rpORF1p\t0\t338\n{q2}\tL1rpORF2p\t0\t1275\n"
            )

    monkeypatch.setattr("haplongliner.sv_mode.run_quiet", fake_run_quiet)
    monkeypatch.setattr("haplongliner.sv_mode.verify_blast_db", lambda x: None)
    from haplongliner.find_longest_orf import find_longest_orf
    monkeypatch.setattr("haplongliner.sv_mode.find_longest_orf", find_longest_orf)
    monkeypatch.setattr(
        "haplongliner.sv_mode.find_intact_orf",
        lambda inp, out: Path(out).write_text(Path(inp).read_text()),
    )

    intact = _validate_orfs(cand)
    assert intact == {"22822069-22828100,122_phaseblock_26,660301"}
