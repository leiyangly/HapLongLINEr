from pathlib import Path

from haplongliner.sv_mode import _validate_orfs


def test_validate_orfs_outputs_only_orfs(monkeypatch, tmp_path):
    cand = tmp_path / "cand.fa"
    cand.write_text(">L1,chr1,0,24,+\n" "ATGGCCATTGTAATGGGCCGCTGAA\n")

    def fake_run_quiet(cmd, check=True, cwd=None):
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
