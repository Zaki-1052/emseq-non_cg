# scripts/utils/extract_direction_flipped.sh
#
# Extract the 760 genes where edgeR calls hypo (logFC < 0, FDR < 0.05)
# but the raw mch_diff is positive (mut rate > ctrl rate).

INPUT="results/03_differential/mch_differential_results.tsv"
OUTPUT="results/03_differential/direction_flipped_genes.tsv"

awk -F'\t' '
  NR == 1 {
    for (i = 1; i <= NF; i++) h[$i] = i
    print $0
    next
  }
  $h["sig_fdr005"] == "TRUE" &&
  $h["edger_logFC"] + 0 < 0 &&
  $h["mch_diff"] + 0 > 0 {
    print $0
  }
' "$INPUT" > "$OUTPUT"

n=$(tail -n +2 "$OUTPUT" | wc -l | tr -d ' ')
echo "$n direction-flipped genes written to $OUTPUT"
