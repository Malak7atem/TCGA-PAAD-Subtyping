# =============================================================================
# TCGA-PAAD Molecular Subtyping Analysis Pipeline
# =============================================================================
# Author      : [Your Name]
# Date        : 2025
# Description : Unsupervised molecular subtyping of pancreatic ductal
#               adenocarcinoma (PAAD) using TCGA gene expression data.
#               Pipeline includes: PCA, K-means clustering, Kaplan-Meier
#               survival analysis, differential expression (limma), and
#               GO/KEGG pathway enrichment (clusterProfiler).
#
# References  :
#   - TCGA PAAD: Cancer Cell (2017). doi:10.1016/j.ccell.2017.07.007
#   - limma    : Ritchie et al., Nucleic Acids Res (2015). doi:10.1093/nar/gkv007
#   - clusterProfiler: Yu et al., OMICS (2012). doi:10.1089/omi.2011.0118
#   - Moffitt subtypes: Nat Genet (2015). doi:10.1038/ng.3398
#   - survminer: Kassambara et al., CRAN (2021)
# =============================================================================

# ── 0. Package Installation & Loading ────────────────────────────────────────

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

cran_packages  <- c("dplyr", "tibble", "ggplot2", "survival", "survminer")
bioc_packages  <- c("limma", "ComplexHeatmap", "clusterProfiler", "org.Hs.eg.db")

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, ask = FALSE)
  library(pkg, character.only = TRUE)
}

# Output directory
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

cat("All packages loaded successfully.\n")


# ── 1. Data Loading & Preprocessing ──────────────────────────────────────────

# Load clinical and gene expression data
clinical_data   <- read.csv("data/Resolved_Clinical_Data.csv")
gene_expression <- read.csv("data/TCGA_PAAD_Gene_Expression.csv", row.names = 1)

# Standardize gene expression sample IDs to 12-character TCGA barcodes
colnames(gene_expression) <- gsub("\\.", "-", substr(colnames(gene_expression), 1, 12))

# Remove duplicated samples from both datasets
clinical_data   <- clinical_data[!duplicated(clinical_data$bcr_patient_barcode), ]
gene_expression <- gene_expression[, !duplicated(colnames(gene_expression))]

# Keep only samples present in both datasets
common_samples  <- intersect(clinical_data$bcr_patient_barcode, colnames(gene_expression))
clinical_data   <- clinical_data[clinical_data$bcr_patient_barcode %in% common_samples, ]
gene_expression <- gene_expression[, common_samples]

# Align row order of clinical data to column order of gene expression matrix
clinical_data <- clinical_data[match(common_samples, clinical_data$bcr_patient_barcode), ]

# Remove zero-variance genes (uninformative for PCA / clustering)
gene_expression <- gene_expression[apply(gene_expression, 1, var) > 0, ]

cat("Clinical data dimensions  :", dim(clinical_data), "\n")
cat("Gene expression dimensions:", dim(gene_expression), "\n")
cat("Common samples            :", length(common_samples), "\n")


# ── 2. PCA & Optimal Cluster Number ──────────────────────────────────────────

# Principal Component Analysis (samples as rows)
pca_result <- prcomp(t(gene_expression), scale. = TRUE)
cat("PCA completed. Variance explained by PC1 & PC2:",
    round(summary(pca_result)$importance[2, 1:2] * 100, 1), "%\n")

# Elbow method: determine optimal number of clusters (k)
# Justification: K-means WSS plot following standard practice
set.seed(123)
wss <- sapply(1:8, function(k) {
  kmeans(pca_result$x[, 1:2], centers = k, nstart = 25)$tot.withinss
})

png("figures/elbow_plot.png", width = 700, height = 500, res = 120)
plot(1:8, wss, type = "b", pch = 19,
     xlab = "Number of Clusters (k)",
     ylab = "Total Within-Cluster Sum of Squares",
     main = "Elbow Method for Optimal k",
     col = "#2C7BB6")
abline(v = 3, lty = 2, col = "red")
dev.off()
cat("Elbow plot saved to figures/elbow_plot.png\n")


# ── 3. K-means Clustering ─────────────────────────────────────────────────────

# K-means with k=3 (supported by elbow plot and PAAD literature: Moffitt et al.)
set.seed(123)
kmeans_result <- kmeans(pca_result$x[, 1:2], centers = 3, nstart = 25)
clinical_data$subtype <- as.factor(kmeans_result$cluster)

cat("\nCluster distribution:\n")
print(table(clinical_data$subtype))

# PCA scatter plot coloured by cluster
pca_df <- data.frame(
  PC1     = pca_result$x[, 1],
  PC2     = pca_result$x[, 2],
  Cluster = as.factor(kmeans_result$cluster)
)

pca_plot <- ggplot(pca_df, aes(PC1, PC2, colour = Cluster)) +
  geom_point(size = 2.5, alpha = 0.8) +
  scale_colour_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A"),
                      labels = paste("Cluster", 1:3)) +
  labs(title    = "PCA of TCGA-PAAD Samples",
       subtitle = "K-means clustering (k = 3)",
       x        = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "%)"),
       y        = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "%)"),
       colour   = "Molecular Subtype") +
  theme_classic(base_size = 13)

ggsave("figures/pca_clusters.png", pca_plot, width = 7, height = 5, dpi = 150)
cat("PCA plot saved to figures/pca_clusters.png\n")


# ── 4. Survival Analysis (Kaplan-Meier) ──────────────────────────────────────

# Build overall survival variables from TCGA clinical columns
clinical_data$OS.time <- as.numeric(as.character(clinical_data$last_contact_days_to))
clinical_data$OS      <- ifelse(clinical_data$vital_status == "Alive", 0, 1)

# Remove patients with missing survival information
surv_data <- clinical_data[!is.na(clinical_data$OS.time) & !is.na(clinical_data$OS), ]
cat("\nPatients with complete survival data:", nrow(surv_data), "\n")

# Kaplan-Meier survival object & fit
surv_object <- Surv(time = surv_data$OS.time, event = surv_data$OS)
surv_fit    <- survfit(surv_object ~ subtype, data = surv_data)

# Plot & save
png("figures/survival_curves.png", width = 900, height = 700, res = 120)
km_plot <- ggsurvplot(
  surv_fit,
  data         = surv_data,
  pval         = TRUE,
  conf.int     = TRUE,
  risk.table   = TRUE,
  ggtheme      = theme_classic(base_size = 13),
  palette      = c("#E41A1C", "#377EB8", "#4DAF4A"),
  legend.title = "Molecular Subtype",
  legend.labs  = paste("Cluster", 1:3),
  title        = "Overall Survival by Molecular Subtype (TCGA-PAAD)"
)
print(km_plot)
dev.off()
cat("Survival curves saved to figures/survival_curves.png\n")


# ── 5. Differential Expression Analysis (limma) ──────────────────────────────

# Design matrix: one column per cluster, no intercept
design <- model.matrix(~ 0 + subtype, data = clinical_data)
colnames(design) <- c("Cluster1", "Cluster2", "Cluster3")

# Fit linear model
lm_fit <- lmFit(gene_expression, design)

# Contrast matrix: all pairwise comparisons
contrast_matrix <- makeContrasts(
  Cluster2_vs_Cluster1 = Cluster2 - Cluster1,
  Cluster3_vs_Cluster1 = Cluster3 - Cluster1,
  Cluster3_vs_Cluster2 = Cluster3 - Cluster2,
  levels = design
)

lm_fit2 <- contrasts.fit(lm_fit, contrast_matrix)
lm_fit2 <- eBayes(lm_fit2)

cat("\nDifferential expression summary:\n")
print(summary(decideTests(lm_fit2)))

# Save top 50 DEGs for each contrast
for (contrast in colnames(contrast_matrix)) {
  top <- topTable(lm_fit2, coef = contrast, number = 50, sort.by = "P")
  write.csv(top, file = paste0("results/Top50_DEGs_", contrast, ".csv"))
}
cat("DEG tables saved to results/\n")

# ── 5a. Volcano Plot ──────────────────────────────────────────────────────────

top_genes_c2v1 <- topTable(lm_fit2, coef = "Cluster2_vs_Cluster1",
                            number = Inf, sort.by = "P")
top_genes_c2v1$gene      <- rownames(top_genes_c2v1)
top_genes_c2v1$sig       <- with(top_genes_c2v1,
                                  ifelse(adj.P.Val < 0.05 & abs(logFC) > 1,
                                         ifelse(logFC > 1, "Up", "Down"), "NS"))

volcano <- ggplot(top_genes_c2v1, aes(logFC, -log10(P.Value), colour = sig)) +
  geom_point(alpha = 0.5, size = 1.2) +
  scale_colour_manual(values = c(Up = "#E41A1C", Down = "#377EB8", NS = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "black") +
  labs(title    = "Volcano Plot: Cluster 2 vs Cluster 1",
       x        = "log2 Fold Change",
       y        = expression(-log[10](P~value)),
       colour   = "Expression") +
  theme_classic(base_size = 13)

ggsave("figures/volcano_C2vsC1.png", volcano, width = 7, height = 5, dpi = 150)
cat("Volcano plot saved to figures/volcano_C2vsC1.png\n")

# ── 5b. Heatmap ───────────────────────────────────────────────────────────────

top50_genes  <- rownames(topTable(lm_fit2, coef = "Cluster2_vs_Cluster1", number = 50))
heatmap_data <- as.matrix(gene_expression[top50_genes, ])

# Column annotation: cluster label
col_annot <- HeatmapAnnotation(
  Cluster = clinical_data$subtype,
  col     = list(Cluster = c("1" = "#E41A1C", "2" = "#377EB8", "3" = "#4DAF4A"))
)

png("figures/heatmap_top50.png", width = 1000, height = 800, res = 120)
draw(Heatmap(
  heatmap_data,
  name                = "Expression",
  top_annotation      = col_annot,
  column_title        = "Samples",
  row_title           = "Top 50 DEGs (Cluster2 vs Cluster1)",
  show_row_names      = TRUE,
  show_column_names   = FALSE,
  cluster_rows        = TRUE,
  cluster_columns     = TRUE,
  row_names_gp        = gpar(fontsize = 7)
))
dev.off()
cat("Heatmap saved to figures/heatmap_top50.png\n")


# ── 6. Pathway Enrichment (GO & KEGG) ────────────────────────────────────────

# Strip Ensembl version numbers (e.g. ENSG00000141510.11 -> ENSG00000141510)
top50_clean <- gsub("\\..*", "", top50_genes)

# Map Ensembl IDs to Entrez IDs
entrez_ids <- bitr(
  top50_clean,
  fromType = "ENSEMBL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)
cat("\nMapped", nrow(entrez_ids), "of", length(top50_clean), "genes to Entrez IDs.\n")

# GO Biological Process enrichment
go_enrich <- enrichGO(
  gene          = entrez_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

if (!is.null(go_enrich) && nrow(go_enrich@result) > 0) {
  go_dot <- dotplot(go_enrich, showCategory = 20) +
    ggtitle("GO Enrichment – Biological Process") +
    theme_classic(base_size = 11)
  ggsave("figures/GO_dotplot.png", go_dot, width = 9, height = 7, dpi = 150)
  write.csv(go_enrich@result, "results/GO_Enrichment_Results.csv", row.names = FALSE)
  cat("GO enrichment results saved.\n")
} else {
  cat("No significant GO terms found.\n")
}

# KEGG pathway enrichment
kegg_enrich <- enrichKEGG(
  gene          = entrez_ids$ENTREZID,
  organism      = "hsa",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)

if (!is.null(kegg_enrich) && nrow(kegg_enrich@result) > 0) {
  kegg_dot <- dotplot(kegg_enrich, showCategory = 20) +
    ggtitle("KEGG Pathway Enrichment") +
    theme_classic(base_size = 11)
  ggsave("figures/KEGG_dotplot.png", kegg_dot, width = 9, height = 7, dpi = 150)
  write.csv(kegg_enrich@result, "results/KEGG_Enrichment_Results.csv", row.names = FALSE)
  cat("KEGG enrichment results saved.\n")
} else {
  cat("No significant KEGG pathways found.\n")
}


# ── 7. Session Info (Reproducibility) ────────────────────────────────────────

cat("\n=== Session Information ===\n")
writeLines(capture.output(sessionInfo()), "results/session_info.txt")
sessionInfo()

cat("\n=== Pipeline Complete ===\n")
cat("All figures saved to:  figures/\n")
cat("All results saved to:  results/\n")
