# TCGA-PAAD Molecular Subtyping Analysis

[![R](https://img.shields.io/badge/R-%3E%3D4.0-blue)](https://www.r-project.org/)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-%3E%3D3.14-green)](https://bioconductor.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Overview

This project performs **unsupervised molecular subtyping** of Pancreatic Ductal Adenocarcinoma (PAAD) using publicly available TCGA (The Cancer Genome Atlas) gene expression data. The pipeline identifies distinct molecular subtypes and characterizes them through survival analysis, differential gene expression, and pathway enrichment.

Pancreatic cancer has one of the lowest 5-year survival rates (~12%) among all cancers. Identifying molecular subtypes with different prognostic outcomes is critical for understanding disease heterogeneity and guiding treatment strategies (Bailey et al., 2016; Moffitt et al., 2015).

---

## Analysis Pipeline

```
TCGA-PAAD Gene Expression Data
           │
           ▼
    Data Preprocessing
  (ID matching, deduplication,
   zero-variance gene removal)
           │
           ▼
     PCA + K-means (k=3)
    Molecular Subtyping
           │
           ▼
   Kaplan-Meier Survival
     Analysis (log-rank)
           │
           ▼
  Differential Expression
      Analysis (limma)
           │
           ▼
  GO & KEGG Pathway
     Enrichment
```

---

## Methods

| Step | Method | Tool / Package | Reference |
|------|--------|---------------|-----------|
| Dimensionality reduction | PCA | Base R `prcomp` | — |
| Clustering | K-means (k=3) | Base R `kmeans` | Moffitt et al., 2015 |
| Survival analysis | Kaplan-Meier + log-rank | `survival`, `survminer` | Kassambara et al., 2021 |
| Differential expression | Empirical Bayes (limma) | `limma` | Ritchie et al., 2015 |
| Visualization | Volcano plot, Heatmap | `ggplot2`, `ComplexHeatmap` | — |
| Pathway enrichment | GO-BP + KEGG | `clusterProfiler` | Yu et al., 2012 |

**Cluster number justification:** k=3 was selected based on the elbow method (WSS plot) and is consistent with published PAAD molecular subtype literature (Moffitt et al. identified classical/basal subtypes; Bailey et al. identified 4 subtypes).

---

## Dataset

- **Source:** [TCGA-PAAD](https://portal.gdc.cancer.gov/projects/TCGA-PAAD) via The Cancer Genome Atlas
- **Data type:** RNA-seq gene expression (normalized counts)
- **Clinical data:** Overall survival, vital status, patient barcodes
- **Access:** Open-access TCGA data (no patient consent required for download)

> **Note:** Raw data files are not included in this repository due to size.  
> Download from [GDC Data Portal](https://portal.gdc.cancer.gov/) or [UCSC Xena](https://xenabrowser.net/).

---

## Repository Structure

```
TCGA-PAAD-Subtyping/
├── analysis_pipeline.R          # Main analysis script (fully documented)
├── README.md                    # This file
├── LICENSE                      # MIT License
├── data/
│   ├── Resolved_Clinical_Data.csv       # TCGA clinical metadata
│   └── TCGA_PAAD_Gene_Expression.csv    # Gene expression matrix
├── results/
│   ├── Top50_DEGs_Cluster2_vs_Cluster1.csv
│   ├── Top50_DEGs_Cluster3_vs_Cluster1.csv
│   ├── Top50_DEGs_Cluster3_vs_Cluster2.csv
│   ├── GO_Enrichment_Results.csv
│   ├── KEGG_Enrichment_Results.csv
│   └── session_info.txt
└── figures/
    ├── elbow_plot.png            # Optimal k selection
    ├── pca_clusters.png          # PCA coloured by cluster
    ├── survival_curves.png       # Kaplan-Meier curves
    ├── volcano_C2vsC1.png        # Volcano plot
    ├── heatmap_top50.png         # Expression heatmap
    ├── GO_dotplot.png            # GO enrichment
    └── KEGG_dotplot.png          # KEGG enrichment
```

## Key Results

| Output | Description |
|--------|-------------|
| `figures/elbow_plot.png` | WSS elbow curve justifying k=3 |
| `figures/pca_clusters.png` | PCA scatter plot with 3 molecular subtypes |
| `figures/survival_curves.png` | KM curves with log-rank p-value |
| `figures/volcano_C2vsC1.png` | DEGs between Cluster 2 and Cluster 1 |
| `figures/heatmap_top50.png` | Top 50 DEG expression heatmap |
| `results/GO_Enrichment_Results.csv` | Significant GO-BP terms |
| `results/KEGG_Enrichment_Results.csv` | Significant KEGG pathways |

---

## References

1. **TCGA PAAD Consortium** (2017). Integrated Genomic Characterization of Pancreatic Ductal Adenocarcinoma. *Cancer Cell*, 32(2), 185–203. https://doi.org/10.1016/j.ccell.2017.07.007

2. **Moffitt, R.A. et al.** (2015). Virtual microdissection identifies distinct tumor- and stroma-specific subtypes of pancreatic ductal adenocarcinoma. *Nature Genetics*, 47, 1168–1178. https://doi.org/10.1038/ng.3398

3. **Bailey, P. et al.** (2016). Genomic analyses identify molecular subtypes of pancreatic cancer. *Nature*, 531, 47–52. https://doi.org/10.1038/nature16965

4. **Ritchie, M.E. et al.** (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Research*, 43(7), e47. https://doi.org/10.1093/nar/gkv007

5. **Yu, G. et al.** (2012). clusterProfiler: an R Package for Comparing Biological Themes Among Gene Clusters. *OMICS: A Journal of Integrative Biology*, 16(5), 284–287. https://doi.org/10.1089/omi.2011.0118

6. **Kassambara, A. et al.** (2021). survminer: Drawing Survival Curves using 'ggplot2'. R package version 0.4.9. https://CRAN.R-project.org/package=survminer

7. **Gu, Z. et al.** (2016). Complex heatmaps reveal patterns and correlations in multidimensional genomic data. *Bioinformatics*, 32(18), 2847–2849. https://doi.org/10.1093/bioinformatics/btw313

---

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.

---

## Author

**Malak Hatem**  
Nile University  
M.hatem2261@nu.edu.eg  
