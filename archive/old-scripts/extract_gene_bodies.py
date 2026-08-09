# scripts/extract_gene_bodies.py
"""
Extract gene body coordinates from a GENCODE GTF into sorted BED format.

Selects gene-level records of a given gene_type on a given set of chromosomes
and writes them as BED (0-based half-open) for intersection with per-site
methylation calls.

Usage:
    python extract_gene_bodies.py \
        --gtf /path/to/gencode.vM25.annotation.gtf.gz \
        --output /path/to/gene_bodies.protein_coding.bed
"""

import argparse
import gzip
import re
import sys
from pathlib import Path

DEFAULT_CHROMS = [f"chr{i}" for i in range(1, 20)] + ["chrX", "chrY"]

GTF_FEATURE_COL = 2
GTF_CHROM_COL = 0
GTF_START_COL = 3
GTF_END_COL = 4
GTF_STRAND_COL = 6
GTF_ATTR_COL = 8

ATTR_PATTERN = re.compile(r'(\S+)\s+"([^"]*)"')

OUTPUT_COLUMNS = [
    "chr", "start", "end", "gene_name", "gene_id",
    "strand", "gene_type", "gene_length",
]


def open_maybe_gzip(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "rt")


def parse_attributes(attr_field: str) -> dict:
    """Parse a GTF attribute column into a dict of key -> value."""
    return dict(ATTR_PATTERN.findall(attr_field))


def chrom_sort_key(chrom: str) -> tuple:
    """Sort chr1..chr19 numerically, then chrX, chrY."""
    suffix = chrom[3:]
    if suffix.isdigit():
        return (0, int(suffix))
    return (1, suffix)


def extract_genes(gtf_path: Path, gene_type: str, chroms: set) -> list:
    """Read a GENCODE GTF and return gene-body records matching the filters."""
    genes = []
    seen_gene_ids = set()

    with open_maybe_gzip(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if fields[GTF_FEATURE_COL] != "gene":
                continue

            chrom = fields[GTF_CHROM_COL]
            if chrom not in chroms:
                continue

            attrs = parse_attributes(fields[GTF_ATTR_COL])
            if attrs.get("gene_type") != gene_type:
                continue

            gene_id = attrs["gene_id"]
            if gene_id in seen_gene_ids:
                raise ValueError(f"Duplicate gene_id in GTF: {gene_id}")
            seen_gene_ids.add(gene_id)

            # GTF is 1-based closed; BED is 0-based half-open.
            start = int(fields[GTF_START_COL]) - 1
            end = int(fields[GTF_END_COL])

            genes.append({
                "chr": chrom,
                "start": start,
                "end": end,
                "gene_name": attrs["gene_name"],
                "gene_id": gene_id,
                "strand": fields[GTF_STRAND_COL],
                "gene_type": attrs["gene_type"],
                "gene_length": end - start,
            })

    return genes


def write_bed(genes: list, output_path: Path):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    genes_sorted = sorted(genes, key=lambda g: (chrom_sort_key(g["chr"]), g["start"], g["end"]))

    with open(output_path, "w") as fh:
        for gene in genes_sorted:
            fh.write("\t".join(str(gene[col]) for col in OUTPUT_COLUMNS) + "\n")


def report_summary(genes: list, gene_type: str):
    """Print gene counts and length distribution to stderr."""
    per_chrom = {}
    for gene in genes:
        per_chrom[gene["chr"]] = per_chrom.get(gene["chr"], 0) + 1

    lengths = sorted(g["gene_length"] for g in genes)
    n = len(lengths)

    def quantile(frac):
        return lengths[int(frac * (n - 1))]

    print(f"Gene type:  {gene_type}", file=sys.stderr)
    print(f"Total genes: {n:,}", file=sys.stderr)
    print("", file=sys.stderr)
    print("Per-chromosome counts:", file=sys.stderr)
    for chrom in sorted(per_chrom, key=chrom_sort_key):
        print(f"  {chrom:<8} {per_chrom[chrom]:>6,}", file=sys.stderr)
    print("", file=sys.stderr)
    print("Gene length distribution (bp):", file=sys.stderr)
    print(f"  min    {lengths[0]:>12,}", file=sys.stderr)
    print(f"  q25    {quantile(0.25):>12,}", file=sys.stderr)
    print(f"  median {quantile(0.50):>12,}", file=sys.stderr)
    print(f"  q75    {quantile(0.75):>12,}", file=sys.stderr)
    print(f"  max    {lengths[-1]:>12,}", file=sys.stderr)
    print("", file=sys.stderr)
    for threshold in (10_000, 50_000, 100_000):
        count = sum(1 for length in lengths if length >= threshold)
        print(f"  genes >= {threshold:>7,} bp: {count:>6,}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Extract gene body coordinates from a GENCODE GTF into sorted BED")
    parser.add_argument("--gtf", type=Path, required=True,
                        help="Input GENCODE GTF (plain or gzipped)")
    parser.add_argument("--output", type=Path, required=True,
                        help="Output BED path")
    parser.add_argument("--gene-type", default="protein_coding",
                        help="gene_type attribute to keep (default: protein_coding)")
    parser.add_argument("--chroms", default=",".join(DEFAULT_CHROMS),
                        help="Comma-separated chromosomes to keep "
                             "(default: chr1-chr19,chrX,chrY)")
    args = parser.parse_args()

    if not args.gtf.exists():
        sys.exit(f"ERROR: GTF not found: {args.gtf}")

    chroms = set(args.chroms.split(","))

    genes = extract_genes(args.gtf, args.gene_type, chroms)
    if not genes:
        sys.exit(f"ERROR: no genes matched gene_type={args.gene_type} on chroms={sorted(chroms)}")

    write_bed(genes, args.output)
    report_summary(genes, args.gene_type)
    print(f"Wrote {len(genes):,} genes to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
