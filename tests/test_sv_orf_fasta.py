from pathlib import Path

from haplongliner.sv_mode import _validate_orfs, _append_liftover_orfs


def test_validate_orfs_outputs_only_orfs(monkeypatch, tmp_path):
    cand = tmp_path / "cand.fa"
    cand.write_text(">L1,chr1,0,24,+\n" "ATGGCCATTGTAATGGGCCGCTGAA\n")

    def fake_run_quiet(cmd, check=True, cwd=None, **kwargs):
        # Simulate getorf and blastp calls by creating expected output files
        if cmd[0] == "getorf":
            out = Path(cmd[-1])
            # Name is converted by getorf (commas -> underscores)
            out.write_text(
                ">L1_chr1_0_24_+_1 [1 - 24]\nMAIVMGR*\n"
            )
        elif cmd[0] == "blastp":
            out = Path(cmd[-1])
            out.write_text("")

    monkeypatch.setattr("haplongliner.sv_mode.run_quiet", fake_run_quiet)
    monkeypatch.setattr("haplongliner.sv_mode.verify_blast_db", lambda x: None)
    monkeypatch.setattr("haplongliner.sv_mode.find_longest_orf", lambda inp, out: Path(out).write_text(""))
    monkeypatch.setattr("haplongliner.sv_mode.find_intact_orf", lambda inp, out: Path(out).write_text(""))

    _validate_orfs(cand)
    orf_fa = cand.with_name("cand_orf.fa")
    content = orf_fa.read_text()
    assert ">L1,chr1,0,24,+\n" not in content
    assert ">L1,chr1,0,24,+,1,1,24" in content


def test_append_liftover_orfs(monkeypatch, tmp_path):
    cand = tmp_path / "cand.fa"
    cand.write_text(">INS_chr1_50,chr1,50,70,.\n" "ATGGCC\n")
    sv_fa = tmp_path / "haplongliner_sv.fa"
    sv_fa.write_text(
        ">chr1,0,24,24,+,ALN\nATGGCCATTGTAATGGGCCGCTGAA\n"
        ">chr1,50,70,20,+,INS\nATGGCC\n"
    )

    def fake_run_quiet(cmd, check=True, cwd=None, **kwargs):
        out = Path(cmd[-1])
        if cmd[0] == "getorf":
            if out.name == "cand_orf.fa":
                out.write_text(
                    ">INS_chr1_50_chr1_50_70_._1 [1 - 6]\nMAIVMG\n"
                )
            else:
                out.write_text(
                    ">chr1_0_24_24_+_ALN_1 [1 - 24]\nMAIVMGR*\n"
                    ">chr1_50_70_20_+_INS_1 [1 - 20]\nMAIVMGR*\n"
                )
        elif cmd[0] == "blastp":
            out.write_text("")

    monkeypatch.setattr("haplongliner.sv_mode.run_quiet", fake_run_quiet)
    monkeypatch.setattr("haplongliner.sv_mode.verify_blast_db", lambda x: None)
    monkeypatch.setattr("haplongliner.sv_mode.find_longest_orf", lambda inp, out: Path(out).write_text(""))
    monkeypatch.setattr("haplongliner.sv_mode.find_intact_orf", lambda inp, out: Path(out).write_text(""))

    _validate_orfs(cand)
    orf_fa = cand.with_name("cand_orf.fa")
    _append_liftover_orfs(orf_fa, sv_fa)

    content = orf_fa.read_text()
    assert ">INS_chr1_50,chr1,50,70,.,1,1,6" in content
    assert ">chr1,0,24,24,+,ALN,1,1,24" in content
    assert ">chr1,50,70,20,+,INS,1,1,20" not in content
