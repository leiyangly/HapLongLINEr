import sys


def find_intact_orf(in_file, out_file):
    """Filter BLASTP combined results for intact ORFs.

    A line is kept only if the ORF1 alignment spans positions 1-338 and
    the ORF2 alignment spans positions 1-1275. This mirrors the behaviour
    of the legacy ``FindIntactORF.pl`` script.
    """
    with open(in_file) as fin, open(out_file, "w") as fout:
        for line in fin:
            if not line.strip():
                continue
            fields = line.strip().split()
            # Ensure there are enough columns
            if len(fields) < 25:
                continue
            if (fields[8] == "1" and fields[9] == "338" and
                    fields[23] == "1" and fields[24] == "1275"):
                fout.write(line)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("Usage: python -m haplongliner.find_intact_orf <in> <out>")
    find_intact_orf(sys.argv[1], sys.argv[2])
