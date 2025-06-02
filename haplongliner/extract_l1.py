def extract_l1_from_bed(infile, outfile=None):
    import sys
    out = open(outfile, "w") if outfile else sys.stdout
    with open(infile) as f:
        for line in f:
            if line.strip() == "" or line.startswith("#"):
                continue
            fields = line.strip().split()
            if len(fields) < 4:
                continue
            name = fields[3]
            start = int(fields[1])
            end = int(fields[2])
            if name.startswith("L1") and (end - start) >= 5000:
                print(line, end="", file=out)
    if outfile:
        out.close()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Extract full-length L1s (>=5000bp) from BED file.")
    parser.add_argument("infile", help="Input BED file")
    parser.add_argument("-o", "--outfile", help="Output file (default: stdout)")
    args = parser.parse_args()
    extract_l1_from_bed(args.infile, args.outfile)