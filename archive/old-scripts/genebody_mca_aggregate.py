# scripts/genebody_mca_aggregate.py
"""
Aggregate CA-context methylation over gene bodies for one sample.

Combines the CA-filtered CHH and CHG calls, intersects them with gene body
coordinates, and sums coverage and methylated reads per gene. A CA site inside
overlapping genes is counted for every gene containing it, since mCA is read by
MeCP2 regardless of which annotation the site falls under.

Usage:
    python genebody_mca_aggregate.py \
        --ca-chh /path/to/sample_CA_CHH.methylKit.gz \
        --ca-chg /path/to/sample_CA_CHG.methylKit.gz \
        --gene-bed /path/to/gene_bodies.protein_coding.bed \
        --sample-id ctrl_M1 \
        --output /path/to/ctrl_M1.genebody_mca.tsv
"""

import argparse
import gzip
import sys
from collections import defaultdict
from pathlib import Path

import pybedtools

DEFAULT_CHROMS = [f"chr{i}" for i in range(1, 20)] + ["chrX", "chrY"]

METHYLKIT_HEADER = ["chrBase", "chr", "base", "strand", "coverage", "freqC", "freqT"]
MK_CHROM_COL = 1
MK_POS_COL = 2
MK_COVERAGE_COL = 4
MK_FREQC_COL = 5

# Layout of `bedtools intersect -wa -wb` output: 5 CA columns then 8 gene columns.
ISECT_COVERAGE_COL = 3
ISECT_FREQC_COL = 4
ISECT_GENE_NAME_COL = 8

GENE_BED_COLUMNS = [
    "chr", "start", "end", "gene_name", "gene_id",
    "strand", "gene_type", "gene_length",
]

OUTPUT_COLUMNS = [
    "gene_name", "gene_id", "chr", "start", "end", "strand",
    "gene_length", "gene_type",
    "n_ca_sites", "total_coverage", "meth_reads", "mca_rate", "sample_id",
]


def load_gene_bed(gene_bed_path: Path) -> list:
    """Read the gene body BED into a list of dicts, preserving file order."""
    genes = []
    with open(gene_bed_path) as fh:
        for line_num, line in enumerate(fh, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != len(GENE_BED_COLUMNS):
                raise ValueError(
                    f"{gene_bed_path}:{line_num}: expected {len(GENE_BED_COLUMNS)} "
                    f"columns, found {len(fields)}")
            genes.append(dict(zip(GENE_BED_COLUMNS, fields)))

    if not genes:
        raise ValueError(f"Gene BED is empty: {gene_bed_path}")
    return genes


def validate_methylkit_header(header_line: str, source: Path):
    columns = header_line.rstrip("\n").split("\t")
    if columns != METHYLKIT_HEADER:
        raise ValueError(
            f"{source}: unexpected methylKit header.\n"
            f"  expected: {METHYLKIT_HEADER}\n"
            f"  found:    {columns}")


def write_ca_bed(ca_paths: list, bed_path: Path, chroms: set) -> int:
    """Write CA sites from the given methylKit files as BED, returning the count.

    methylKit positions are 1-based; BED is 0-based half-open.
    """
    n_sites = 0

    with open(bed_path, "w") as out:
        for ca_path in ca_paths:
            with gzip.open(ca_path, "rt") as fh:
                validate_methylkit_header(fh.readline(), ca_path)

                for line in fh:
                    fields = line.rstrip("\n").split("\t")
                    chrom = fields[MK_CHROM_COL]
                    if chrom not in chroms:
                        continue

                    pos = int(fields[MK_POS_COL])
                    out.write(f"{chrom}\t{pos - 1}\t{pos}\t"
                              f"{fields[MK_COVERAGE_COL]}\t{fields[MK_FREQC_COL]}\n")
                    n_sites += 1

    if n_sites == 0:
        raise ValueError(
            f"No CA sites on {sorted(chroms)} found in {[str(p) for p in ca_paths]}")
    return n_sites


def aggregate_by_gene(intersect_path: str) -> dict:
    """Sum CA site count, coverage, and methylated reads per gene."""
    totals = defaultdict(lambda: {"n_ca_sites": 0, "total_coverage": 0, "meth_reads": 0.0})

    with open(intersect_path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            coverage = int(fields[ISECT_COVERAGE_COL])
            freq_c = float(fields[ISECT_FREQC_COL])

            gene_totals = totals[fields[ISECT_GENE_NAME_COL]]
            gene_totals["n_ca_sites"] += 1
            gene_totals["total_coverage"] += coverage
            gene_totals["meth_reads"] += coverage * freq_c / 100.0

    return totals


def write_output(genes: list, totals: dict, sample_id: str, output_path: Path):
    """Write one row per gene, including genes with no covered CA sites."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as fh:
        fh.write("\t".join(OUTPUT_COLUMNS) + "\n")

        for gene in genes:
            gene_totals = totals.get(
                gene["gene_name"],
                {"n_ca_sites": 0, "total_coverage": 0, "meth_reads": 0.0})

            total_coverage = gene_totals["total_coverage"]
            meth_reads = gene_totals["meth_reads"]
            mca_rate = meth_reads / total_coverage if total_coverage > 0 else 0.0

            fh.write("\t".join([
                gene["gene_name"], gene["gene_id"], gene["chr"],
                gene["start"], gene["end"], gene["strand"],
                gene["gene_length"], gene["gene_type"],
                str(gene_totals["n_ca_sites"]), str(total_coverage),
                f"{meth_reads:.4f}", f"{mca_rate:.8f}", sample_id,
            ]) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Aggregate CA-context methylation over gene bodies for one sample")
    parser.add_argument("--ca-chh", type=Path, required=True,
                        help="CA-filtered CHH methylKit file")
    parser.add_argument("--ca-chg", type=Path, required=True,
                        help="CA-filtered CHG methylKit file")
    parser.add_argument("--gene-bed", type=Path, required=True,
                        help="Gene body BED from extract_gene_bodies.py")
    parser.add_argument("--sample-id", required=True,
                        help="Sample identifier written into the output")
    parser.add_argument("--output", type=Path, required=True,
                        help="Output TSV path")
    parser.add_argument("--work-dir", type=Path, default=None,
                        help="Directory for intermediate files "
                             "(default: <output parent>/work)")
    parser.add_argument("--chroms", default=",".join(DEFAULT_CHROMS),
                        help="Comma-separated chromosomes to keep "
                             "(default: chr1-chr19,chrX,chrY)")
    parser.add_argument("--keep-intermediate", action="store_true",
                        help="Keep the intermediate CA BED file for inspection")
    args = parser.parse_args()

    for path in (args.ca_chh, args.ca_chg, args.gene_bed):
        if not path.exists():
            sys.exit(f"ERROR: input not found: {path}")

    chroms = set(args.chroms.split(","))
    work_dir = args.work_dir if args.work_dir is not None else args.output.parent / "work"
    work_dir.mkdir(parents=True, exist_ok=True)

    # Keep bedtools scratch inside the project tree rather than the system temp.
    pybedtools.set_tempdir(str(work_dir))

    genes = load_gene_bed(args.gene_bed)
    print(f"Loaded {len(genes):,} genes from {args.gene_bed}", flush=True)

    ca_bed_path = work_dir / f"{args.sample_id}.ca_sites.bed"
    print(f"Writing CA sites to {ca_bed_path} ...", flush=True)
    n_sites = write_ca_bed([args.ca_chh, args.ca_chg], ca_bed_path, chroms)
    print(f"  {n_sites:,} CA sites on {len(chroms)} chromosomes", flush=True)

    print("Intersecting CA sites with gene bodies ...", flush=True)
    ca_bed = pybedtools.BedTool(str(ca_bed_path))
    intersection = ca_bed.intersect(str(args.gene_bed), wa=True, wb=True)

    print("Aggregating per gene ...", flush=True)
    totals = aggregate_by_gene(intersection.fn)

    covered_genes = sum(1 for gene in genes if gene["gene_name"] in totals)
    sites_in_genes = sum(entry["n_ca_sites"] for entry in totals.values())
    print(f"  {covered_genes:,} / {len(genes):,} genes have CA sites", flush=True)
    print(f"  {sites_in_genes:,} site-gene assignments", flush=True)

    write_output(genes, totals, args.sample_id, args.output)
    print(f"Wrote {args.output}", flush=True)

    pybedtools.cleanup()
    if not args.keep_intermediate:
        ca_bed_path.unlink()
        print(f"Removed intermediate {ca_bed_path}", flush=True)


if __name__ == "__main__":
    main()
