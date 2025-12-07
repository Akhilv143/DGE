# Differential Gene Expression Analysis of the *airway* RNA-seq Dataset

This repository contains an end-to-end differential gene expression (DGE) workflow using the Bioconductor **DESeq2** package on the human airway smooth muscle RNA-seq dataset from the **airway** experiment package.

## Overview

- Perform quality filtering and normalization of gene-level count data.
- Fit a negative binomial model with DESeq2 to compare dexamethasone-treated vs untreated samples.
- Visualize results with MA plot, PCA, heatmap, and volcano plot.
- Extract significant genes and run Gene Ontology enrichment with **clusterProfiler**.

The goal is to provide a minimal, reproducible example of an RNA-seq DGE workflow that others can clone and adapt to their own datasets.

## Data and experimental design

The analysis uses the Bioconductor **airway** package, which summarizes an RNA-seq experiment on four human airway smooth muscle cell lines, each profiled in untreated and dexamethasone-treated conditions.

- 8 samples total: 4 cell lines × 2 conditions (treated vs untreated).
- Input count matrix: `counts_data.csv` (genes × samples) derived from the `airway` RangedSummarizedExperiment object.
- Sample metadata: `sample_info.csv`, containing cell line IDs and dexamethasone treatment status.

The design formula in DESeq2 is `~ dexamethasone`, testing for differences between treated and untreated samples.

## Environment and dependencies

This workflow was run in R using Bioconductor packages for count-based RNA-seq analysis.

Core R packages:

- `DESeq2` - differential expression modelling and statistics.
- `airway` - example RNA-seq dataset.
- `tidyverse` - data wrangling and plotting helpers.
- `EnhancedVolcano` - volcano plot visualization.
- `pheatmap` - clustered heatmaps of transformed counts.
- `org.Hs.eg.db` and `clusterProfiler` - gene annotation and GO enrichment.

Install required Bioconductor packages (run once):

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c(
  "DESeq2",
  "airway",
  "org.Hs.eg.db",
  "clusterProfiler"
))
install.packages(c("tidyverse", "pheatmap", "EnhancedVolcano"))
```

## Reproducing the analysis

1. **Clone the repository**

```bash
git clone https://github.com/Akhilv143/DGE.git
cd DGE
```

2. **Open the R project / folder in RStudio**

- Make sure working directory points to the folder containing `DGE_script.R`, `counts_data.csv`, and `sample_info.csv`.

3. **Run the analysis script**

- Open `DGE_script.R`.
- Source the script line-by-line or click "Source" to run the full workflow.
- The script performs:
  - Construction of a `DESeqDataSet` from counts and sample metadata.
  - Pre-filtering of low-count genes and releveling of the treatment factor.
  - Normalization, dispersion estimation, and model fitting via `DESeq()`.
  - Extraction of differential expression results for treated vs untreated (`res`).
  - Filtering for significant genes (padj and log2 fold-change thresholds).
  - Generation of MA plot, PCA plot, clustered heatmap, and volcano plot.
  - Annotation of significant genes and GO enrichment of Biological Process terms.

4. **Check outputs**

- CSV tables:
  - `DESeq2_results_all.csv` - all genes with statistics.
  - `DESeq2_results_significant.csv` and `DESeq2_significant_genes.csv` - filtered sets of differentially expressed genes.
- Figures:
  - `MA_plot.jpeg` - global log2 fold change vs mean expression.
  - `PCA_plot.jpeg` - sample clustering by dexamethasone treatment.
  - `heatmap.jpeg` - top variable or DE genes across samples.
  - `volcano_plot.jpeg` - log2 fold change vs -log10 adjusted p-value.

## Interpretation and next steps

- Inspect volcano plot and ranked tables to identify strongly up- and down-regulated genes; notable candidates in the original airway study included glucocorticoid-responsive genes linked to asthma biology.
- Use GO enrichment results to highlight pathways associated with steroid response, inflammation, and airway smooth muscle function.
- Adapt the script for more complex experimental designs (e.g. multiple factors or interactions) by extending the DESeq2 design formula and contrast specification.

## References

- Love et al., "RNA-Seq workflow: gene-level exploratory analysis and differential expression," F1000Research, describing the *airway* dataset and DESeq2-based analysis.
- Bioconductor DESeq2 vignette for detailed explanations of model assumptions, design formulas, and advanced use cases.
- Bioconductor **airway** package documentation for details on the experimental setup and raw data sources.

### Useful Resources

- [RNA-seq Analysis Course](https://statomics.github.io/SGA/airway_salmon_edgeR.html)
- [DESeq2 DGE Analysis on Kaggle](https://www.kaggle.com/datasets/mannekuntanagendra/deseq2-dge-analysis-airway-dataset)
- [Bioconductor RNA-seq Course Materials](https://www.bioconductor.org/help/course-materials/2015/Uruguay2015/V6-RNASeq.html)
- [RNA-seq Analysis Guide](https://bdsr.stephenturner.us/rnaseq.html)
- [F1000Research Article](https://f1000research.com/articles/4-1070)
- [Advanced RNA-seq Lab](https://www.bioconductor.org/help/course-materials/2015/BioC2015/V6-Lab-Advanced-RnaSeq.html)
- [DESeq2 Full Vignette](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
- [RNA-seq Simulation Guide](https://cran.r-project.org/web/packages/seqgendiff/vignettes/simulate_rnaseq.html)
