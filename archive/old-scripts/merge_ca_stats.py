# scripts/merge_ca_stats.py
"""
Combine per-file CA statistics from ca_filter.py into run-level summaries.

Produces the per-file/per-category/per-chromosome summary table, a per-sample
spike-in table giving the conversion noise floor measured on unmethylated
phage lambda, and runs QC gates that stop the pipeline before gene-body
aggregation if any input looks wrong.

Usage:
    python merge_ca_stats.py \
        --stats-dir /path/to/ca_filtered \
        --summary-out /path/to/ca_summary_stats.tsv \
        --spike-in-out /path/to/spike_in_rates.tsv
"""

import argparse
import json
import re
import sys
from pathlib import Path

LAMBDA_CATEGORY = "phage_lambda"

EXPECTED_STATS_FILES = 16
EXPECTED_SAMPLES = 8

SAMPLE_PATTERN = re.compile(
    r"_index_(?P<well>[A-H]\d)_"
    r".*?_(?P<genotype>ctrl|mut)_(?P<sex>[MF])(?P<replicate>\d+)_"
    r"S(?P<snum>\d+)_(?P<context>CHH|CHG)\.methylKit\.gz$"
)

SUMMARY_COLUMNS = [
    "file", "category", "chromosome",
    "total_sites", "ca_sites", "ca_fraction",
    "total_cov", "ca_cov",
    "total_rate_pct", "ca_rate_pct",
]

SPIKE_IN_COLUMNS = [
    "sample_id", "genotype", "sex",
    "lambda_ca_sites", "lambda_ca_coverage", "lambda_ca_meth_reads",
    "lambda_ca_rate",
]


def parse_sample_metadata(filename: str) -> dict:
    """Derive sample identity and covariates from a methylKit filename.

    Batch is the trailing digit of the index well: E1/F1/G1/H1 are plate 1,
    A2/B2/C2/D2 are plate 2.
    """
    match = SAMPLE_PATTERN.search(filename)
    if match is None:
        raise ValueError(f"Cannot parse sample metadata from filename: {filename}")

    genotype = match.group("genotype")
    sex = match.group("sex")
    replicate = match.group("replicate")
    well = match.group("well")

    return {
        "sample_id": f"{genotype}_{sex}{replicate}",
        "genotype": genotype,
        "sex": sex,
        "context": match.group("context"),
    }


def load_stats(stats_dir: Path) -> list:
    """Load every *.ca_stats.json in stats_dir."""
    paths = sorted(stats_dir.glob("*.ca_stats.json"))
    if len(paths) != EXPECTED_STATS_FILES:
        raise ValueError(
            f"Expected {EXPECTED_STATS_FILES} .ca_stats.json files in {stats_dir}, "
            f"found {len(paths)}: {[p.name for p in paths]}")

    records = []
    for path in paths:
        with open(path) as fh:
            stats = json.load(fh)
        stats["metadata"] = parse_sample_metadata(stats["input_file"])
        records.append(stats)
    return records


def rate_pct(meth_reads: float, coverage: int) -> float:
    return meth_reads / coverage * 100 if coverage > 0 else 0.0


def fraction(numerator: int, denominator: int) -> float:
    return numerator / denominator if denominator > 0 else 0.0


def write_summary(records: list, summary_path: Path):
    """Write the per-file, per-category, per-chromosome CA summary table."""
    summary_path.parent.mkdir(parents=True, exist_ok=True)

    with open(summary_path, "w") as fh:
        fh.write("\t".join(SUMMARY_COLUMNS) + "\n")

        for stats in sorted(records, key=lambda r: r["input_file"]):
            filename = stats["input_file"]

            for cat_name, cat_data in sorted(stats["categories"].items()):
                fh.write("\t".join([
                    filename, cat_name, "ALL",
                    str(cat_data["total_sites"]), str(cat_data["ca_sites"]),
                    f"{fraction(cat_data['ca_sites'], cat_data['total_sites']):.6f}",
                    str(cat_data["total_cov"]), str(cat_data["ca_cov"]),
                    f"{rate_pct(cat_data['total_meth_reads'], cat_data['total_cov']):.6f}",
                    f"{rate_pct(cat_data['ca_meth_reads'], cat_data['ca_cov']):.6f}",
                ]) + "\n")

                for chrom, chr_data in sorted(cat_data["per_chr"].items()):
                    fh.write("\t".join([
                        filename, cat_name, chrom,
                        str(chr_data["total_sites"]), str(chr_data["ca_sites"]),
                        f"{fraction(chr_data['ca_sites'], chr_data['total_sites']):.6f}",
                        str(chr_data["total_cov"]), str(chr_data["ca_cov"]),
                        f"{rate_pct(chr_data['total_meth_reads'], chr_data['total_cov']):.6f}",
                        f"{rate_pct(chr_data['ca_meth_reads'], chr_data['ca_cov']):.6f}",
                    ]) + "\n")


def aggregate_spike_in(records: list) -> list:
    """Sum lambda CA counts across CHH and CHG for each sample."""
    per_sample = {}

    for stats in records:
        meta = stats["metadata"]
        sample_id = meta["sample_id"]

        if LAMBDA_CATEGORY not in stats["categories"]:
            raise ValueError(
                f"No {LAMBDA_CATEGORY} data in {stats['input_file']}; the conversion "
                f"noise floor cannot be estimated for sample {sample_id}")

        lambda_data = stats["categories"][LAMBDA_CATEGORY]
        entry = per_sample.setdefault(sample_id, {
            "sample_id": sample_id,
            "genotype": meta["genotype"],
            "sex": meta["sex"],
            "lambda_ca_sites": 0,
            "lambda_ca_coverage": 0,
            "lambda_ca_meth_reads": 0.0,
        })

        entry["lambda_ca_sites"] += lambda_data["ca_sites"]
        entry["lambda_ca_coverage"] += lambda_data["ca_cov"]
        entry["lambda_ca_meth_reads"] += lambda_data["ca_meth_reads"]

    for entry in per_sample.values():
        entry["lambda_ca_rate"] = fraction(
            entry["lambda_ca_meth_reads"], entry["lambda_ca_coverage"])

    return sorted(per_sample.values(), key=lambda e: e["sample_id"])


def write_spike_in(spike_in_rows: list, spike_in_path: Path):
    spike_in_path.parent.mkdir(parents=True, exist_ok=True)

    with open(spike_in_path, "w") as fh:
        fh.write("\t".join(SPIKE_IN_COLUMNS) + "\n")
        for row in spike_in_rows:
            fh.write("\t".join([
                row["sample_id"], row["genotype"], row["sex"],
                str(row["lambda_ca_sites"]), str(row["lambda_ca_coverage"]),
                f"{row['lambda_ca_meth_reads']:.2f}", f"{row['lambda_ca_rate']:.8f}",
            ]) + "\n")


def check_ca_fractions(records: list, min_fraction: float, max_fraction: float) -> list:
    """CA sites should be a stable share of all non-CG sites in every file."""
    failures = []
    print("CA fraction per file:")
    for stats in sorted(records, key=lambda r: r["input_file"]):
        observed = fraction(stats["ca_sites"], stats["total_sites"])
        status = "ok" if min_fraction <= observed <= max_fraction else "FAIL"
        print(f"  {stats['input_file']:<70} {observed:>7.4f}  {status}")
        if status == "FAIL":
            failures.append(
                f"{stats['input_file']}: CA fraction {observed:.4f} outside "
                f"[{min_fraction}, {max_fraction}]")
    print("")
    return failures


def check_lambda_rates(spike_in_rows: list, max_rate: float) -> list:
    """Apparent methylation on unmethylated lambda is the conversion noise floor."""
    failures = []
    print("Lambda CA rate per sample (conversion noise floor):")
    for row in spike_in_rows:
        observed = row["lambda_ca_rate"]
        status = "ok" if observed <= max_rate else "FAIL"
        print(f"  {row['sample_id']:<10} {observed:>9.6f}  "
              f"({observed * 100:.4f}%)  cov={row['lambda_ca_coverage']:>10,}  {status}")
        if status == "FAIL":
            failures.append(
                f"{row['sample_id']}: lambda CA rate {observed:.6f} exceeds {max_rate}")
    print("")
    return failures


def check_sample_completeness(spike_in_rows: list) -> list:
    """Every sample must contribute both a CHH and a CHG file."""
    failures = []
    if len(spike_in_rows) != EXPECTED_SAMPLES:
        failures.append(
            f"Expected {EXPECTED_SAMPLES} samples, found {len(spike_in_rows)}: "
            f"{[r['sample_id'] for r in spike_in_rows]}")

    genotype_counts = {}
    for row in spike_in_rows:
        genotype_counts[row["genotype"]] = genotype_counts.get(row["genotype"], 0) + 1
    print(f"Samples: {len(spike_in_rows)}  {dict(sorted(genotype_counts.items()))}")
    print("")
    return failures


def main():
    parser = argparse.ArgumentParser(
        description="Combine CA filter statistics and extract spike-in noise floors")
    parser.add_argument("--stats-dir", type=Path, required=True,
                        help="Directory containing *.ca_stats.json from ca_filter.py")
    parser.add_argument("--summary-out", type=Path, required=True,
                        help="Output path for the combined summary TSV")
    parser.add_argument("--spike-in-out", type=Path, required=True,
                        help="Output path for the per-sample spike-in rate TSV")
    parser.add_argument("--min-ca-fraction", type=float, default=0.25,
                        help="Lower bound on the CA share of non-CG sites (default: 0.25)")
    parser.add_argument("--max-ca-fraction", type=float, default=0.50,
                        help="Upper bound on the CA share of non-CG sites (default: 0.50)")
    parser.add_argument("--max-lambda-rate", type=float, default=0.02,
                        help="Maximum tolerated lambda CA rate (default: 0.02)")
    args = parser.parse_args()

    if not args.stats_dir.is_dir():
        sys.exit(f"ERROR: stats directory not found: {args.stats_dir}")

    records = load_stats(args.stats_dir)
    spike_in_rows = aggregate_spike_in(records)

    write_summary(records, args.summary_out)
    write_spike_in(spike_in_rows, args.spike_in_out)

    print(f"Summary written to  {args.summary_out}")
    print(f"Spike-in written to {args.spike_in_out}")
    print("")
    print("=== QC gates ===")
    print("")

    failures = []
    failures += check_sample_completeness(spike_in_rows)
    failures += check_ca_fractions(records, args.min_ca_fraction, args.max_ca_fraction)
    failures += check_lambda_rates(spike_in_rows, args.max_lambda_rate)

    if failures:
        print("QC FAILED:", file=sys.stderr)
        for failure in failures:
            print(f"  {failure}", file=sys.stderr)
        sys.exit(1)

    print("All QC gates passed.")


if __name__ == "__main__":
    main()
