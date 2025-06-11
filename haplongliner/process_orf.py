import re


def process_orf_fasta(in_fasta, out_bed):
    """Parse getorf FASTA output and convert it to a BED-like table.

    Mirrors the behaviour of ``ProcessORF.pl`` from the legacy pipeline.
    Only header lines from ``getorf`` are examined. The output columns are::

        l1_id  start  end  strand  length  l1_strand
    """

    with open(in_fasta) as fin, open(out_bed, "w") as fout:
        for line in fin:
            if not line.startswith(">"):
                continue
            fields = line.strip().split()
            if len(fields) < 4:
                continue
            start = int(fields[1].lstrip("["))
            end = int(fields[3].rstrip("]"))

            header = fields[0][1:]
            m = re.match(r"^(.+?)_(\d+)_(\d+)_([+-])(?:_.*)?$", header)
            if not m:
                continue
            chrom, lstart, lend, l1_strand = m.groups()
            l1_id = f"{chrom}_{lstart}_{lend}"
            strand = "+" if end >= start else "-"
            length = abs(end - start)
            fout.write(
                f"{l1_id}\t{start}\t{end}\t{strand}\t{length}\t{l1_strand}\n"
            )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Process getorf FASTA output to BED format"
    )
    parser.add_argument("in_fasta", help="FASTA file from getorf")
    parser.add_argument("out_bed", help="Output BED file")
    args = parser.parse_args()
    process_orf_fasta(args.in_fasta, args.out_bed)
