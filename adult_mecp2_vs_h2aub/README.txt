LOG2LOG2 QUADRANT PLOT
1. MeCP2 and H2AK119Ub diffbind output filtered for only Conc >=4. All peaks that do not map to a clear genomic region were removed.
2. Each peak in both files was designated to a gene. If not in gene body, it was not used.
3. There may be multiple peaks per gene in each file. So, I used this calculation on both files to come up with one FC value for MeCP2 and one FC value for H2AUb per gene:
    * Fold=∑(peak Fold×peak weight) / ∑(peak weight)
        * peak weight= | Fold Change | × −log10 (FDR)
4. Significance (pink label) was given for genes that at least one sig peak in BOTH H2aUb and MeCP2

TRIED, BUT MAY NOT BE GOOD
Significance was weighted:
    * gene FDR=∑(FDR ×peak weight) / ∑(peak weight)
        * peak weight= | Fold Change | × −log10 (FDR)
    1. Genes were labeled pink if they were sig in BOTH MeCP2 and H2AUb diffbind results.

Effect size
* abs(mecp2_Fold)
* How much MeCP2 binding changed.
* Bigger fold change = larger weight.
Statistical confidence
* -log10(FDR)
* How significant the peak is.
* Smaller FDR = larger weight.

GO ANALYSIS
1. ONLY selected genes in quadrant 1 where they were statistically significant.
2. enrichGO function and ENTREZIDs
3. BP, MF, and CC were all run