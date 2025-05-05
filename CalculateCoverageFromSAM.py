#!/usr/bin/env python3

import argparse
import csv
import re
import sys


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Compute per-base coverage for a single reference from a SAM file"
    )
    parser.add_argument("sam_path", help="Path to the input SAM file")
    parser.add_argument("ref_name", help="Reference sequence name to tally coverage for")
    return parser.parse_args()


def get_ref_length_and_header(sam_path, ref_name):
    """
    Read SAM headers to find the length of the specified reference.
    Exits with an error if the reference is not found.
    """
    with open(sam_path) as f:
        for line in f:
            if not line.startswith("@SQ"):
                continue
            # Header fields: SN:<name>  LN:<length>
            fields = dict(kv.split(":", 1) for kv in line.strip().split("\t")[1:])
            if fields.get("SN") == ref_name:
                return int(fields["LN"])
    sys.exit(f"ERROR: reference '{ref_name}' not found in SAM headers")


def increment_coverage(coverage, pos, cigar):
    """
    Update coverage list based on one alignment.
    
    Args:
        coverage: list of ints, per-base counts (0-based index)
        pos:       1-based leftmost alignment position on reference
        cigar:     CIGAR string describing alignment
    """
    ref_pos = pos  # 1-based coordinate
    # Parse CIGAR operations and lengths
    for length_str, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar):
        length = int(length_str)
        if op in ("M", "=", "X"):
            # Count matches/mismatches toward coverage
            start = ref_pos - 1
            end = start + length
            for i in range(start, end):
                coverage[i] += 1
            ref_pos += length
        elif op in ("D", "N"):
            # Deletion or skipped region advances reference pos
            ref_pos += length
        # Other operations (I, S, H, P) do not consume reference bases


def main():
    """Main workflow: determine reference length, tally coverage, write CSV."""
    args = parse_args()
    sam_path = args.sam_path
    ref_name = args.ref_name

    # Step 1: Determine reference length from SAM header
    ref_length = get_ref_length_and_header(sam_path, ref_name)
    coverage = [0] * ref_length

    # Step 2: Scan alignments and increment coverage
    with open(sam_path) as f:
        for line in f:
            if line.startswith("@"):
                continue
            cols = line.rstrip().split("\t")
            if cols[2] != ref_name:
                continue
            flag = int(cols[1])
            if flag & 0x4:
                continue  # skip unmapped reads
            pos = int(cols[3])
            cigar = cols[5]
            if cigar != "*":
                increment_coverage(coverage, pos, cigar)

    # Step 3: Write per-site coverage to CSV
    out_csv = f"{ref_name}_coverage.csv"
    with open(out_csv, "w", newline="") as outf:
        writer = csv.writer(outf)
        writer.writerow(["Site", "Coverage"])
        for site, count in enumerate(coverage, start=1):
            writer.writerow([site, count])

    print(f"Wrote coverage for '{ref_name}' ({ref_length} bp) to '{out_csv}'")


if __name__ == "__main__":
    main()
