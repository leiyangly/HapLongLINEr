#!/usr/bin/env python3
"""Simplistic liftover for PAF alignments.

This script replicates the ``liftover`` command implemented in
``scripts/paftools.js`` from the minimap2 package.  It reads a PAF
alignment file and a BED file containing query intervals and outputs the
lifted coordinates on the target assembly as BED6.
"""
from __future__ import annotations

import argparse
import sys
from typing import Dict, Iterable, List, Tuple


# ---------------------------------------------------------------------------
# interval utilities (port of Interval.* helpers in paftools.js)
# ---------------------------------------------------------------------------

Interval = List[List]


def interval_sort(a: Interval) -> None:
    a.sort(key=lambda x: (x[0], x[1]))


def interval_merge(a: Interval) -> None:
    if not a:
        return
    interval_sort(a)
    k = 0
    for i in range(1, len(a)):
        if a[k][1] >= a[i][0]:
            if a[i][1] > a[k][1]:
                a[k][1] = a[i][1]
        else:
            k += 1
            a[k] = a[i][:]
    del a[k + 1 :]


def interval_index_end(a: Interval) -> None:
    if not a:
        return
    interval_sort(a)
    a[0].append(0)
    k = 0
    k_en = a[0][1]
    for i in range(1, len(a)):
        if k_en <= a[i][0]:
            k += 1
            while k < i and a[k][1] <= a[i][0]:
                k += 1
            k_en = a[k][1]
        a[i].append(k)


def interval_find_intv(a: Interval, x: int) -> int:
    left = -1
    right = len(a)
    while right - left > 1:
        mid = (left + right) // 2
        if a[mid][0] > x:
            right = mid
        elif a[mid][0] < x:
            left = mid
        else:
            return mid
    return left


def interval_find_ovlp(a: Interval, st: int, en: int) -> Interval:
    if not a or st >= en:
        return []
    l = interval_find_intv(a, st)
    k = 0 if l < 0 else a[l][-1]
    b: Interval = []
    for i in range(k, len(a)):
        if a[i][0] >= en:
            break
        if st < a[i][1]:
            b.append(a[i])
    return b


# ---------------------------------------------------------------------------
# BED parsing helper
# ---------------------------------------------------------------------------


def read_bed(path: str | None, to_merge: bool) -> Dict[str, Interval]:
    if path is None:
        return {}
    bed: Dict[str, Interval] = {}
    fh = sys.stdin if path == "-" else open(path)
    with fh:
        for line in fh:
            if not line.strip():
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 3:
                continue
            chrom, start, end = fields[:3]
            name = fields[3] if len(fields) > 3 else ""
            score = fields[4] if len(fields) > 4 else "0"
            strand = fields[5] if len(fields) > 5 else "+"
            try:
                s = int(start)
                e = int(end)
            except ValueError:
                continue
            bed.setdefault(chrom, []).append([s, e, strand, name, score])
    for chrom in list(bed.keys()):
        interval_sort(bed[chrom])
        if to_merge:
            interval_merge(bed[chrom])
        interval_index_end(bed[chrom])
    return bed


# ---------------------------------------------------------------------------
# liftover implementation
# ---------------------------------------------------------------------------


def liftover_paf(
    paf_path: str,
    bed_path: str,
    min_mapq: int,
    min_len: int,
    max_div: float,
    merge: bool,
) -> None:
    bed = read_bed(bed_path, merge)
    fh = sys.stdin if paf_path == "-" else open(paf_path)
    re_cigar = re_tag = None
    import re as _re

    re_cigar = _re.compile(r"(\d+)([MID])")
    re_tag = _re.compile(r"^(\S\S):([AZif]):(\S+)$")

    with fh:
        for line in fh:
            if not line.strip():
                continue
            t = line.rstrip().split("\t")
            if t[0] not in bed:
                continue
            # parse tp and cg tags
            tp = None
            cg = None
            for field in t[12:]:
                m = re_tag.match(field)
                if not m:
                    continue
                if m.group(1) == "tp":
                    tp = m.group(3)
                elif m.group(1) == "cg":
                    cg = m.group(3)
            if tp not in {"P", "I"}:
                continue
            if cg is None:
                raise RuntimeError("unable to find the 'cg' tag")
            try:
                for i in range(1, 4):
                    t[i] = int(t[i])
                for i in range(6, 12):
                    t[i] = int(t[i])
            except ValueError:
                continue
            if t[11] < min_mapq or t[10] < min_len:
                continue
            regs = interval_find_ovlp(bed[t[0]], t[2], t[3])
            if not regs:
                continue
            if 0.0 <= max_div < 1.0:
                n_gaps = 0
                n_opens = 0
                for m in re_cigar.finditer(cg):
                    if m.group(2) in {"I", "D"}:
                        n_gaps += int(m.group(1))
                        n_opens += 1
                n_mm = t[10] - t[9] - n_gaps
                n_diff2 = n_mm + n_opens
                if n_diff2 / (n_diff2 + t[9]) > max_div:
                    continue
            # extract start and end positions
            a: List[List[int]] = []
            r: List[List[int]] = []
            strand = t[4]
            for i, (s, e, *_) in enumerate(regs):
                if strand == "+":
                    a.append([s, 0, i, -2])
                    a.append([e - 1, 1, i, -2])
                else:
                    a.append([t[1] - e, 0, i, -2])
                    a.append([t[1] - s - 1, 1, i, -2])
                r.append([-2, -2])
            a.sort(key=lambda x: x[0])

            k = 0
            x = t[7]
            y = t[2] if strand == "+" else t[1] - t[3]
            for m in re_cigar.finditer(cg):
                length = int(m.group(1))
                op = m.group(2)
                if op == "D":
                    x += length
                    continue
                while k < len(a) and a[k][0] < y:
                    k += 1
                for i in range(k, len(a)):
                    if y <= a[i][0] < y + length:
                        a[i][3] = x + (a[i][0] - y) if op == "M" else x
                    else:
                        break
                y += length
                if op == "M":
                    x += length
            if (
                x != t[8]
                or (strand == "+" and y != t[3])
                or (strand == "-" and y != t[1] - t[2])
            ):
                raise RuntimeError("CIGAR is inconsistent with mapping coordinates")

            for qpos, flag, idx, coord in a:
                if flag == 0:
                    r[idx][0] = coord
                else:
                    r[idx][1] = coord + 1
            for i, (s, e, b_strand, b_name, b_score, *_) in enumerate(regs):
                name = f"{t[0]},{s},{e},{b_strand},{b_name}"
                start = r[i][0]
                end = r[i][1]
                if start < 0:
                    name += ",t5"
                    start = t[7]
                if end < 0:
                    name += ",t3"
                    end = t[8]
                out_strand = t[4] if b_strand == "+" else ("+" if t[4] == "-" else "-")
                print(
                    "\t".join(
                        [t[5], str(start), str(end), name, str(b_score), out_strand]
                    )
                )


# ---------------------------------------------------------------------------
# command line interface
# ---------------------------------------------------------------------------


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Simplistic liftOver for PAF alignments. "
            "By default this script uses permissive filtering (-q 0 -l 0 -d 2.0). "
            "To emulate the more stringent paftools.js behaviour, use "
            "'-q 5 -l 50000 -d 2.0'."
        )
    )
    parser.add_argument("aln_paf", help="PAF alignment file or '-' for stdin")
    parser.add_argument("query_bed", help="BED file on query sequences")
    parser.add_argument("-m", action="store_true", help="merge BED intervals")
    parser.add_argument(
        "-q",
        type=int,
        default=0,
        help="min mapping quality [0]",
    )
    parser.add_argument(
        "-l",
        type=int,
        default=0,
        help="min alignment length [0]",
    )
    parser.add_argument(
        "-d",
        type=float,
        default=2.0,
        help="max sequence divergence (>=1 to disable) [2.0]",
    )
    args = parser.parse_args(argv)

    try:
        liftover_paf(args.aln_paf, args.query_bed, args.q, args.l, args.d, args.m)
    except BrokenPipeError:
        pass
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
