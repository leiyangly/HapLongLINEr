import re


def find_longest_orf(blastp_file, out_file):
    """Select the longest ORF1 and ORF2 alignment for each L1.

    Mirrors the behaviour of ``FindLongestORF.pl``. ``blastp_file`` should be
    generated with ``-outfmt '6 std qlen slen sacc'`` against the ORF1/2
    reference database.
    """

    len_dict = {}
    info = {}
    index = set()
    with open(blastp_file) as fh:
        for line in fh:
            if not line.strip():
                continue
            F = line.strip().split()
            name = re.sub(r"_[0-9]+$", "", F[0])
            subject = F[1]
            aln_len = int(F[3])
            if aln_len >= len_dict.get(name, {}).get(subject, 0):
                len_dict.setdefault(name, {})[subject] = aln_len
                info.setdefault(name, {})[subject] = line.strip()
            index.add(name)

    with open(out_file, "w") as out:
        for key in index:
            if "L1rpORF1p" in info.get(key, {}) and "L1rpORF2p" in info.get(key, {}):
                out.write(
                    f"{info[key]['L1rpORF1p']}\t{info[key]['L1rpORF2p']}\n"
                )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Select longest ORF1/ORF2 alignments from BLASTP results"
    )
    parser.add_argument("blastp_file", help="BLASTP tabular output")
    parser.add_argument("out_file", help="Output file")
    args = parser.parse_args()
    find_longest_orf(args.blastp_file, args.out_file)
