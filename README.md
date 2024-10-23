# Pea Embryo
## Source Data
**PeaEmbryo.xlsx and PeaEndosperm.xlsx** are metabolomics data collected from pea:

* 2 genotypes: Frisson (green pea, all sample names start with “sF” ) and CDC Meadow (yellow pea, all sample names start with “sM”)
* 4 embryo development stages: D1, D2, D3 and D4
* 2 tissue types: embryo (“E” in sample name) and endosperm (“En” in sample name)
* 3 reps  

Total: 24 samples for embryos; 24 samples for endosperm

**psat.csv** is the pathway information for Pisum sativum (pea) is retrieved from [KEGG Pisum sativum](https://www.kegg.jp/kegg-bin/show_organism?menu_type=pathway_maps&org=psat)

## Script
* **name mapping**: find KEGG IDs for metabolites lacking external IDs
* **dam**: differentially accumulated metabolites and volcano plots
* **enrichment**: enrichment analysis based on differentially accumulated metabolites annotation information
* **mummichog**: [Mummichog](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003123) is a method that predicts metabolic network activity directly from untargeted metabolomics data without the need for prior metabolite identification.  
* **compare_enrichment_mummichog**: compare result from enrichment analysis and mummchiog 






## Methods
Enrichment analysis was conducted by using two different methods: a commonly used approach and the Mummichog (detailed parameters are provided in the next section). Afterward, the hypergeometric test was performed to calculate the p-value for the overlap in results between these two methods.

The results show that the median p-value across all comparisons is approximately 0.0103, while the p-value calculated using the average number of pathways and their overlap is 0.0074 (did not directly average all the p-values, as a single p-value of 1 could severely distort the mean, even if the remaining seven values are close to 0, pushing the average up to 0.125).

These findings indicate a statistically significant overlap between the two methods, suggesting that the shared elements between them are unlikely to be due to random chance, we can get a similar conclusion from the two methods.