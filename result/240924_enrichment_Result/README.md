

## New KEGG Database
|              |Num. of Comp.|Previous KEGG Coverage| Latest KEGG Coverage |
|---|---|---|---|
| Pea Embryo      |873|150 (17%)|425 (48%)|   
| Pea Endosperm     |1159|225 (19%)|581 (50%)|


## Cross-validation of Two Enrichment Methods
Enrichment analysis was conducted by using two different methods: a commonly used approach and the Mummichog (detailed parameters are provided in the next section). Afterward, the hypergeometric test was performed to calculate the p-value for the overlap in results between these two methods.

The results show that the median p-value across all comparisons is approximately **0.0103**, while the p-value calculated using the average number of pathways and their overlap is **0.0074** (did not directly average all the p-values, as a single p-value of 1 could severely distort the mean, even if the remaining seven values are close to 0, pushing the average up to 0.125).

These findings indicate a statistically significant overlap between the two methods, suggesting that the shared elements between them are unlikely to be due to random chance, we can get a similar conclusion from the two methods.

The complete results are shown below:

| Group            | Number of Mummichog | Number of Enrichment | Number of Common Elements | Percentage of Common Elements in Enrichment (p < 0.05) | p-value     |
|------------------|---------------------|----------------------|---------------------------|-------------------------------------------------------|-------------|
| sMEnD1\_vs\_sMED1  | 19                  | 15                   | 8                         | 53.33                                                 | 0.001012    |
| sMEnD2\_vs\_sMED2  | 16                  | 18                   | 9                         | 50.00                                                 | 0.053551    |
| sMEnD3\_vs\_sMED3  | 18                  | 10                   | 4                         | 40.00                                                 | 0.153645    |
| sMEnD4\_vs\_sMED4  | 19                  | 10                   | 4                         | 40.00                                                 | 0.153645    |
| sFEnD2\_vs\_sFEnD1 | 23                  | 14                   | 6                         | 42.86                                                 | 0.051054    |
| sFEnD4\_vs\_sFEnD3 | 14                  | 14                   | 6                         | 42.86                                                 | 0.051054    |
| sFEnD3\_vs\_sFEnD2 | 17                  | 16                   | 6                         | 37.50                                                 | 0.109860    |
| sMED3\_vs\_sMED2   | 19                  | 7                    | 3                         | 42.86                                                 | 0.307746    |
| sMED4\_vs\_sMED3   | 11                  | 10                   | 2                         | 20.00                                                 | 0.555607    |
| sMED2\_vs\_sMED1   | 10                  | 14                   | 5                         | 35.71                                                 | 0.231563    |
| sFED4\_vs\_sMED4   | 14                  | 16                   | 7                         | 43.75                                                 | 0.082763    |
| sFED1\_vs\_sMED1   | 14                  | 8                    | 5                         | 62.50                                                 | 0.015671    |
| sFED2\_vs\_sMED2   | 24                  | 9                    | 2                         | 22.22                                                 | 0.789981    |
| sFED3\_vs\_sMED3   | 30                  | 7                    | 3                         | 42.86                                                 | 0.249143    |
| sFEnD3\_vs\_sMEnD3 | 26                  | 9                    | 4                         | 44.44                                                 | 0.053551    |
| sFEnD1\_vs\_sMEnD1 | 12                  | 10                   | 5                         | 50.00                                                 | 0.012391    |
| sFEnD2\_vs\_sMEnD2 | 12                  | 7                    | 0                         | 0.00                                                  | 1.000000    |
| sFEnD4\_vs\_sMEnD4 | 21                  | 12                   | 5                         | 41.67                                                 | 0.087911    |
| sMEnD2\_vs\_sMEnD1 | 23                  | 19                   | 11                        | 57.89                                                 | 0.008362    |
| sMEnD3\_vs\_sMEnD2 | 19                  | 17                   | 7                         | 41.18                                                 | 0.077271    |
| sMEnD4\_vs\_sMEnD3 | 18                  | 19                   | 7                         | 36.84                                                 | 0.077271    |
| sFEnD2\_vs\_sFED2  | 13                  | 8                    | 2                         | 25.00                                                 | 0.151012    |
| sFEnD3\_vs\_sFED3  | 26                  | 9                    | 4                         | 44.44                                                 | 0.053551    |
| sFEnD1\_vs\_sFED1  | 21                  | 12                   | 8                         | 66.67                                                 | 0.000017    |
| sFEnD4\_vs\_sFED4  | 14                  | 11                   | 1                         | 9.09                                                  | 0.683662    |


## psat.csv
The pathway information for Pisum sativum (pea) is retrieved from [KEGG Pisum sativum](https://www.kegg.jp/kegg-bin/show_organism?menu_type=pathway_maps&org=psat)
## mummichog\_pathway\_enrichment.csv
[Mummichog](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003123) is a method that predicts metabolic network activity directly from untargeted metabolomics data without the need for prior metabolite identification.  

The exact mass tolerance is 10 ppm, and the pathway library used is thale cress (pea is not provided in this software), other parameters are default.  

Since Mummichog does not perform metabolite annotation in advance, it cannot be used to generate an enrichment analysis network. Therefore, it is only utilized here for cross-validation with the commonly used method.
## enrichment.csv
**Num**: the total number of metabolites in that pathway  
**Num\_N**: the number of detected metabolites in that pathway  
**Num_n**: the number of differentially accumulated metabolites in that pathway   
**P value**: are calculated using hypergeometric test, and only keep the top 15 for plotting  

## enrichment.pdf
Each node represents a pathway; the node size reflects the number of differentially accumulated metabolites,
and the color represents -log10(p-value).
  
The edges’ width and color indicate the number of shared differentially accumulated metabolites between
pathways.  

The layout used is the Kamada-Kawai layout.
degree.csv shows the sharing of metabolites between pathways.

## volcano_plot.pptx
All the volcano plot
