import sys


def find_intact_orf(in_file, out_file):
    """Filter BLASTP combined results for intact ORFs.

    A line is kept only if the ORF1 alignment spans positions 1-338 and
    the ORF2 alignment spans positions 1-1275. This mirrors the behaviour
    of the legacy ``FindIntactORF.pl`` script.
    """
    with open(in_file) as fin, open(out_file, "w") as fout:
        for raw_line in fin:
            if not raw_line.strip():
                continue
            fields = raw_line.strip().split()
            if len(fields) < 30:
                # ``cand_orf_combine.blastp`` stores an ORF1 and ORF2 hit on the
                # same line.  Fewer than 30 fields means one of the alignments
                # is missing, so the candidate cannot be intact.
                continue

            def _is_intact(chunk_start: int, expected_len: int) -> bool:
                """Return ``True`` when the BLAST chunk spans the expected ORF."""

                chunk = fields[chunk_start:chunk_start + 15]
                if len(chunk) < 15:
                    return False
                try:
                    qstart = int(chunk[6])
                    qend = int(chunk[7])
                    sstart = int(chunk[8])
                    send = int(chunk[9])
                    qlen = int(chunk[12])
                    slen = int(chunk[13])
                except ValueError:
                    return False

                q_min, q_max = sorted((qstart, qend))
                s_min, s_max = sorted((sstart, send))

                return (
                    q_min == 1
                    and q_max == expected_len
                    and qlen == expected_len
                    and s_min == 1
                    and s_max == expected_len
                    and slen == expected_len
                )

            if _is_intact(0, 338) and _is_intact(15, 1275):
                fout.write(raw_line)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("Usage: python -m haplongliner.find_intact_orf <in> <out>")
    find_intact_orf(sys.argv[1], sys.argv[2])
