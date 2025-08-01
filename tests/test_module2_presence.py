from pathlib import Path
from haplongliner.module2_SV import _validate_presence


def test_validate_presence_ignores_empty(tmp_path, monkeypatch):
    fa = tmp_path / "cand.fa"
    # one entry with no sequence, one valid entry
    fa.write_text(">empty\n>ok\nATGC\n")

    calls = {}

    def fake_run(cmd, **kwargs):
        calls['cmd'] = cmd
        out_idx = cmd.index('-out') + 1
        Path(cmd[out_idx]).write_text("")
        class Dummy:
            returncode = 0
            stdout = b""
            stderr = b""
            def check_returncode(self):
                pass
        return Dummy()

    monkeypatch.setattr("haplongliner.module2_SV.run_quiet", fake_run)
    names = _validate_presence(fa, min_length=1)
    assert names == set()
    assert 'cmd' in calls
    q_idx = calls['cmd'].index('-query') + 1
    q_path = Path(calls['cmd'][q_idx])
    q_text = q_path.read_text()
    assert '>empty' not in q_text
    assert '>ok' in q_text
