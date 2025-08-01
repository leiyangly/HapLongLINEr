import re


def process_orf_fasta(in_fasta, out_bed):
    """Parse getorf FASTA output and convert it to a BED-like table.

    Mirrors the behaviour of ``ProcessORF.pl`` from the legacy pipeline.
    Only header lines from ``getorf`` are examined. The output columns are::

        l1_id  start  end  strand  length  l1_strand
    """

    pattern_new = re.compile(r"^(.+),(\d+),(\d+),(\d+)$")
    pattern_old = re.compile(r"^(.+?)[,_](\d+)[,_](\d+)[,_]([+-])(?:[,_].*)?$")

    with open(in_fasta) as fin, open(out_bed, "w") as fout:
        for line in fin:
            if not line.startswith(">"):
                continue
            line = line.strip()

            header = line[1:]
            m_new = pattern_new.match(header)
            if m_new:
                base, idx, pos_start, pos_end = m_new.groups()
                pos_start = int(pos_start)
                pos_end = int(pos_end)
                parts = base.split(",")
                if len(parts) >= 5:
                    chrom, lstart, lend, _length, l1_strand = parts[:5]
                else:
                    # Fallback: derive values from the tail of the list
                    chrom, lstart, lend = parts[0], parts[1], parts[2]
                    l1_strand = parts[-1]
            else:
                fields = line.split()
                if len(fields) < 4:
                    continue
                pos_start = int(fields[1].lstrip("["))
                pos_end = int(fields[3].rstrip("]"))
                header = fields[0][1:]
                m_old = pattern_old.match(header)
                if not m_old:
                    continue
                chrom, lstart, lend, l1_strand = m_old.groups()

            # getorf reports coordinates as 1-based inclusive. Convert to
            # 0-based half-open to match the rest of the pipeline.
            start = min(pos_start, pos_end) - 1
            end = max(pos_start, pos_end)
            strand = "+" if pos_end >= pos_start else "-"

            l1_id = f"{chrom}_{lstart}_{lend}"
            length = end - start
            fout.write(f"{l1_id}\t{start}\t{end}\t{strand}\t{length}\t{l1_strand}\n")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Process getorf FASTA output to BED format"
    )
    parser.add_argument("in_fasta", help="FASTA file from getorf")
    parser.add_argument("out_bed", help="Output BED file")
    args = parser.parse_args()
    process_orf_fasta(args.in_fasta, args.out_bed)
