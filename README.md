# Genetic similarity among 178 disease phenotypes predicts therapeutic and side effects for 1,711 drugs

## Panagiotis N. Lalagkas & Rachel D. Melamed
Link to the [preprint](https://www.medrxiv.org/content/10.1101/2025.05.13.25327511v1).

This repository contains the data and source code used in the analyses for the above manuscript. The repository is organized as follows:
- [data](https://github.com/lalagkaspn/genetic_similarity_drug_discovery/tree/main/data) includes data required to reproduce the analyses. When raw data cannot be publicly shared, we provide either a summary version or processed outputs.
- [src](https://github.com/lalagkaspn/genetic_similarity_drug_discovery/tree/main/src) contains all source code for data preparation, model training and evaluation, and figure/table generation.
- [figures](https://github.com/lalagkaspn/genetic_similarity_drug_discovery/tree/main/figures) stores all generated figures, both main and supplementary.
- [results](https://github.com/lalagkaspn/genetic_similarity_drug_discovery/tree/main/results) includes the predicted drug indications and side effects, as well as results from the direct genetic evidence analysis for our predictions (used for Supplementary Figure 10)
- [tables](https://github.com/lalagkaspn/genetic_similarity_drug_discovery/tree/main/tables) contains two supplementary tables:
    - `Supplementary table 1`: clean list of predicted drug indications
    - `Supplementary table 2`: clean list of predicted drug side effects
    - Both tables include the predicted probabiliies and the corrresponding known drug indications and side effects whose genetic similarity to the test disease was used for making the prediction.

If you have any questions, please reach out to [panagiotis.lalagkas@gmail.com](mailto:panagiotis.lalagkas@gmail.com)

All analyses were performed in `R 4.2.3`. R packages used (CRAN/Bioconductor):
- binom: v1.1.1.1
- biomaRt: v2.54.1
- data.table: v.1.16.0
- DescTools: v0.99.48
- dplyr: v1.1.4
- GenomicRanges: v1.50.2
- ggplot2: v3.4.2
- ggpubr: v0.6.0
- ggrepel: v0.9.3
- ggsignif: v0.6.4
- glue: v1.6.2
- grid: v4.2.3
- gridExtra: v2.3
- lsa: v0.73.3
- openxlsx: v4.2.5.2
- parallel: v4.2.3
- pROC: v1.18.0
- reshape2: v1.4.4
- rstan: v2.32.7
- stringr: v1.5.1
- tidyr: v1.3.0
- viridisLite: v0.4.2
