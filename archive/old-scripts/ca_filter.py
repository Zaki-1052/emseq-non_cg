# scripts/ca_filter.py
"""
CA-context filtering for EM-seq methylKit files.

Reads CHH and CHG methylKit files, looks up each position in the reference
genome to determine trinucleotide context, writes CA-only sites to new files,
and produces per-file summary statistics as JSON.

Requires pysam and a FASTA reference with .fai index.

Usage:
    python ca_filter.py \
        --ref /path/to/mm10+controls.fa \
        --input-dir /path/to/methylKit/files \
        --output-dir ./ca_filtered/ \
        [--file-index N]
"""

import argparse
import gzip
import json
import sys
from collections import defaultdict
from pathlib import Path

import pysam

AUTOSOMES = {f"chr{i}" for i in range(1, 20)}
SEX_CHROMS = {"chrX", "chrY"}
SPIKE_INS = {"phage_lambda", "plasmid_puc19c", "phage_T4", "phage_Xp12"}

# Marks files this script has already produced, so a rerun pointed at an
# output directory does not re-filter its own output.
CA_FILTERED_MARKER = "_CA_"


def discover_input_files(input_dir: Path) -> list:
    """Return sorted CHH/CHG methylKit inputs, excluding CA-filtered output."""
    candidates = sorted(input_dir.glob("*_CH[HG].methylKit.gz"))
    files = [path for path in candidates if CA_FILTERED_MARKER not in path.name]
    if not files:
        raise FileNotFoundError(
            f"No *_CHH.methylKit.gz or *_CHG.methylKit.gz files in {input_dir}")
    return files


def categorize_chrom(chrom: str) -> str:
    if chrom in AUTOSOMES:
        return "autosomes"
    if chrom in SEX_CHROMS:
        return "sex"
    if chrom in SPIKE_INS:
        return chrom
    return "other"


def is_ca_context(fasta: pysam.FastaFile, chrom: str, pos_1based: int, strand: str) -> bool:
    """Check if a cytosine site is in CA context.

    methylKit positions are 1-based. pysam.fetch uses 0-based half-open intervals.

    For forward strand (F): C at pos, next base at pos+1
        -> fetch(chrom, pos, pos+1) to get the base immediately after C
        -> CA if next_base == 'A'

    For reverse strand (R): G on forward strand at pos (complement of C)
        -> The C on reverse strand 'sees' the base to its 3' side on the reverse strand,
           which is the base to the 5' side on the forward strand (pos - 1 on forward, 1-based)
        -> fetch(chrom, pos-2, pos-1) to get the forward-strand base at position pos-1
        -> CA on reverse strand if that forward base == 'T' (complement of A)
    """
    try:
        if strand == "F":
            # Next base after C on forward strand (0-based: pos to pos+1)
            next_base = fasta.fetch(chrom, pos_1based, pos_1based + 1).upper()
            return next_base == "A"
        elif strand == "R":
            # Base before G on forward strand = base at pos-1 (1-based)
            # 0-based fetch: pos_1based - 2 to pos_1based - 1
            if pos_1based < 2:
                return False
            prev_base = fasta.fetch(chrom, pos_1based - 2, pos_1based - 1).upper()
            return prev_base == "T"
    except (ValueError, KeyError):
        return False
    return False


def process_file(fasta: pysam.FastaFile, input_path: Path, output_path: Path) -> dict:
    """Filter a single methylKit file for CA-context sites.

    Returns summary stats dict.
    """
    stats = {
        "total_sites": 0,
        "ca_sites": 0,
        "categories": defaultdict(lambda: {
            "total_sites": 0, "ca_sites": 0,
            "total_cov": 0, "ca_cov": 0,
            "total_meth_reads": 0.0, "ca_meth_reads": 0.0,
            "per_chr": defaultdict(lambda: {
                "total_sites": 0, "ca_sites": 0,
                "total_cov": 0, "ca_cov": 0,
                "total_meth_reads": 0.0, "ca_meth_reads": 0.0,
            }),
        }),
    }

    with gzip.open(input_path, "rt") as fin, gzip.open(output_path, "wt") as fout:
        header = fin.readline()
        fout.write(header)

        for line_num, line in enumerate(fin, start=1):
            parts = line.rstrip("\n").split("\t")
            chrom = parts[1]
            pos = int(parts[2])
            strand = parts[3]
            coverage = int(parts[4])
            freq_c = float(parts[5])
            meth_reads = coverage * freq_c / 100.0

            cat = categorize_chrom(chrom)
            stats["total_sites"] += 1
            stats["categories"][cat]["total_sites"] += 1
            stats["categories"][cat]["total_cov"] += coverage
            stats["categories"][cat]["total_meth_reads"] += meth_reads
            stats["categories"][cat]["per_chr"][chrom]["total_sites"] += 1
            stats["categories"][cat]["per_chr"][chrom]["total_cov"] += coverage
            stats["categories"][cat]["per_chr"][chrom]["total_meth_reads"] += meth_reads

            ca = is_ca_context(fasta, chrom, pos, strand)
            if ca:
                fout.write(line)
                stats["ca_sites"] += 1
                stats["categories"][cat]["ca_sites"] += 1
                stats["categories"][cat]["ca_cov"] += coverage
                stats["categories"][cat]["ca_meth_reads"] += meth_reads
                stats["categories"][cat]["per_chr"][chrom]["ca_sites"] += 1
                stats["categories"][cat]["per_chr"][chrom]["ca_cov"] += coverage
                stats["categories"][cat]["per_chr"][chrom]["ca_meth_reads"] += meth_reads

            if line_num % 10_000_000 == 0:
                pct_ca = stats["ca_sites"] / stats["total_sites"] * 100 if stats["total_sites"] > 0 else 0
                print(f"    {line_num:>12,} lines processed, {stats['ca_sites']:>10,} CA sites ({pct_ca:.1f}%)",
                      flush=True)

    return stats


def write_stats_json(stats: dict, input_name: str, output_name: str, stats_path: Path):
    """Write one file's CA statistics as JSON.

    Written per input file so SLURM array tasks never contend for a shared
    output. merge_ca_stats.py combines these into the run-level summary.
    """
    payload = {
        "input_file": input_name,
        "output_file": output_name,
        "total_sites": stats["total_sites"],
        "ca_sites": stats["ca_sites"],
        "categories": {
            cat_name: {
                "total_sites": cat_data["total_sites"],
                "ca_sites": cat_data["ca_sites"],
                "total_cov": cat_data["total_cov"],
                "ca_cov": cat_data["ca_cov"],
                "total_meth_reads": cat_data["total_meth_reads"],
                "ca_meth_reads": cat_data["ca_meth_reads"],
                "per_chr": {
                    chrom: dict(chr_data)
                    for chrom, chr_data in sorted(cat_data["per_chr"].items())
                },
            }
            for cat_name, cat_data in sorted(stats["categories"].items())
        },
    }

    with open(stats_path, "w") as fh:
        json.dump(payload, fh, indent=2)


def main():
    parser = argparse.ArgumentParser(description="Filter methylKit files for CA-context sites")
    parser.add_argument("--ref", type=Path, required=True,
                        help="Reference FASTA (must have .fai index)")
    parser.add_argument("--input-dir", type=Path, required=True,
                        help="Directory containing *_CHH/*_CHG methylKit.gz files")
    parser.add_argument("--output-dir", type=Path, required=True,
                        help="Output directory for CA-filtered files and per-file stats")
    parser.add_argument("--file-index", type=int, default=None,
                        help="Process only this file from the discovered list (0-based), "
                             "for SLURM array tasks. Omit to process every file.")
    parser.add_argument("--list-files", action="store_true",
                        help="Print the discovered files with their indices and exit")
    args = parser.parse_args()

    input_files = discover_input_files(args.input_dir)

    if args.list_files:
        for index, path in enumerate(input_files):
            print(f"{index}\t{path.name}")
        return

    if not args.ref.exists():
        sys.exit(f"ERROR: Reference not found: {args.ref}")
    fai_path = Path(str(args.ref) + ".fai")
    if not fai_path.exists():
        sys.exit(f"ERROR: Reference index not found: {fai_path} (run samtools faidx)")

    if args.file_index is not None:
        if not 0 <= args.file_index < len(input_files):
            sys.exit(f"ERROR: --file-index {args.file_index} out of range "
                     f"[0, {len(input_files)}) for {args.input_dir}")
        input_files = [input_files[args.file_index]]

    args.output_dir.mkdir(parents=True, exist_ok=True)
    fasta = pysam.FastaFile(str(args.ref))

    for input_path in input_files:
        out_basename = (input_path.name
                        .replace("_CHH.", "_CA_CHH.")
                        .replace("_CHG.", "_CA_CHG."))
        output_path = args.output_dir / out_basename
        stats_path = args.output_dir / f"{input_path.name}.ca_stats.json"

        print(f"\nProcessing {input_path.name} ...", flush=True)
        stats = process_file(fasta, input_path, output_path)

        ca_pct = stats["ca_sites"] / stats["total_sites"] * 100 if stats["total_sites"] > 0 else 0
        print(f"  Total sites: {stats['total_sites']:,}, "
              f"CA sites: {stats['ca_sites']:,} ({ca_pct:.1f}%)")

        write_stats_json(stats, input_path.name, out_basename, stats_path)
        print(f"  Stats written to {stats_path}")

    fasta.close()
    print("\nDone.")


if __name__ == "__main__":
    main()

