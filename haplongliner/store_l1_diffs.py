import sqlite3
from typing import Tuple

import edlib
from Bio import SeqIO


def _revcomp(seq: str) -> str:
    complement = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(complement)[::-1]


def _best_alignment(query: str, reference: str) -> Tuple[str, str]:
    """Return orientation (+/-) and cigar string for best alignment."""
    plus = edlib.align(query, reference, mode="NW", task="path")
    minus = edlib.align(_revcomp(query), reference, mode="NW", task="path")
    if minus["editDistance"] < plus["editDistance"]:
        return "-", minus["cigar"]
    return "+", plus["cigar"]


def store_diffs(fasta: str, reference: str = "data/L1rp.fa", db: str = "l1rp_diff.db") -> None:
    """Align sequences in *fasta* to *reference* and store differences in *db*."""
    ref_record = next(SeqIO.parse(reference, "fasta"))
    ref_seq = str(ref_record.seq)

    conn = sqlite3.connect(db)
    cur = conn.cursor()
    cur.execute(
        "CREATE TABLE IF NOT EXISTS diffs (name TEXT PRIMARY KEY, orientation TEXT, cigar TEXT)"
    )

    for record in SeqIO.parse(fasta, "fasta"):
        orient, cigar = _best_alignment(str(record.seq), ref_seq)
        cur.execute(
            "INSERT OR REPLACE INTO diffs (name, orientation, cigar) VALUES (?, ?, ?)",
            (record.id, orient, cigar),
        )

    conn.commit()
    conn.close()


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Align sequences to L1rp reference and store differences"
    )
    parser.add_argument("fasta", help="Input FASTA file")
    parser.add_argument(
        "-r", "--reference", default="data/L1rp.fa", help="Reference FASTA"
    )
    parser.add_argument(
        "-d", "--db", default="l1rp_diff.db", help="SQLite database path"
    )
    args = parser.parse_args()
    store_diffs(args.fasta, args.reference, args.db)
