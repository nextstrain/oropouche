"""
Parse a Nextclade TSV file provided on STDIN and write out a simplified TSV to
STDOUT of the strain names which which the Nextclade TSV indicates were
successfully mapped. The resulting TSV (on STDOUT) has three columns:
"accession" (corresponding to 'seqName' in the Nextclade TSV), "segmment_{SEGMENT}"
(where {SEGMENT} is provided via --segment) and "qc_{SEGMENT}" which corresponds to
'qc.overallStatus' in the Nextclade TSV.
"""

from sys import stdin,stdout,stderr
import argparse
import csv

if __name__=="__main__":
    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument('segment', type=str, help="segment name")
    args = parser.parse_args()

    n_lines = 0
    n_hits = 0
    output_columns = ["accession", f"segment_{args.segment}", f"qc_{args.segment}"]

    reader = csv.DictReader(stdin, delimiter='\t')
    writer = csv.DictWriter(stdout, output_columns, extrasaction='ignore', delimiter='\t', lineterminator='\n')
    writer.writeheader()

    for row in reader:
        n_lines += 1
        # empty QC means no alignment, which we interpret as "not this segment"
        if not row['qc.overallStatus']:
            continue
        writer.writerow(dict(zip(output_columns, [row['seqName'], "1", row['qc.overallStatus']])))
        n_hits += 1

    print(f"Nextclade aligned {n_hits}/{n_lines} sequences to segment {args.segment}", file=stderr)
