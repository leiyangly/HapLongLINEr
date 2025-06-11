import sys
import re
from typing import Dict


def _read_minimap(file_path: str) -> Dict[str, str]:
    result = {}
    with open(file_path) as fh:
        for line in fh:
            if not line.strip():
                continue
            fields = line.strip().split()
            if len(fields) < 4:
                continue
            key = fields[0]
            aln_len = int(fields[3]) - int(fields[2])
            if aln_len >= 200 and (key not in result or "_" in result[key]):
                result[key] = line.strip()
    return result


def _read_intact(file_path: str) -> Dict[str, str]:
    result = {}
    with open(file_path) as fh:
        for line in fh:
            if not line.strip():
                continue
            fields = line.strip().split()
            m = re.match(r"^(.+?)_(\d+)_(\d+)_([+-])(?:_.*)?$", fields[0])
            if not m:
                continue
            chrom, start, end, _ = m.groups()
            key = f"{chrom}_{start}_{end}"
            result[key] = line.strip()
    return result


def combine_table(plus_file: str, minus_file: str, intact_file: str, fl_bed: str, out_file: str) -> None:
    plus = _read_minimap(plus_file)
    minus = _read_minimap(minus_file)
    intact = _read_intact(intact_file)

    with open(fl_bed) as bed, open(out_file, "w") as out:
        for line in bed:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.strip().split()
            if len(f) < 6:
                continue
            chrom, start, end, name, dot, strand = f[:6]
            start_i = int(start)
            end_i = int(end)
            # Convert to the coordinate system used by the legacy pipeline.
            m = start_i - 1999
            p = end_i + 2000
            px = end_i + 1
            mx = start_i + 1
            mkey = f"{chrom}:{m}-{start_i}"
            pkey = f"{chrom}:{px}-{p}"
            # ORF headers store 1-based coordinates
            ikey = f"{chrom}_{mx}_{end_i}"

            m_val = minus.get(mkey, "").split("\t") if minus.get(mkey) else []
            p_val = plus.get(pkey, "").split("\t") if plus.get(pkey) else []

            status = "present"
            if ikey in intact:
                status = "intact"

            chr_ref = "NA"
            if len(m_val) > 5 and len(p_val) > 5 and m_val[5] == p_val[5]:
                chr_ref = m_val[5]

            out_strand = "NA"
            if len(m_val) > 4 and len(p_val) > 4:
                if m_val[4] == p_val[4] and m_val[4] == strand:
                    out_strand = "+"
                elif m_val[4] == p_val[4] and m_val[4] != strand:
                    out_strand = "-"

            start_ref = "NA"
            end_ref = "NA"
            if len(m_val) > 8 and len(p_val) > 8 and len(m_val) > 8 and m_val[5] == p_val[5]:
                if int(m_val[7]) <= int(p_val[8]):
                    start_ref = m_val[8]
                    end_ref = p_val[7]
                if int(p_val[7]) <= int(m_val[8]):
                    start_ref = p_val[8]
                    end_ref = m_val[7]
                if start_ref != "NA" and end_ref != "NA" and int(end_ref) < int(start_ref):
                    start_ref, end_ref = end_ref, start_ref

            out.write(
                f"{chrom}_{start}_{end}_{strand}_{dot}_{name}_{status}\t{chr_ref}_{start_ref}_{end_ref}_{out_strand}\n"
            )


if __name__ == "__main__":
    if len(sys.argv) != 6:
        sys.exit(
            "Usage: python -m haplongliner.combine_table <plus> <minus> <intact> <bed> <out>"
        )
    combine_table(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
