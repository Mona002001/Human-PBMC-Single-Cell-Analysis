
# ============================================================
# Single-Cell RNA-seq Analysis of Human PBMCs Using Seurat
# Dataset: 10x Genomics Human PBMC 10k (v2)
# ============================================================


# 1. Load required packages

library(Seurat)
library(ggplot2)
library(dplyr)
library(tibble)
library(patchwork)
library(DoubletFinder)


# 2. Define project directories

data_dir <- "C:/Users/aline/OneDrive/Desktop/filtered_feature_bc_matrix"

dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)


# 3. Load 10x Genomics data

pbmc.data <- Read10X(data.dir = data_dir)

# Select gene expression data
pbmc.counts <- pbmc.data[["Gene Expression"]]

# Create Seurat object
pbmc <- CreateSeuratObject(
  counts = pbmc.counts,
  project = "Human_PBMC",
  min.cells = 3,
  min.features = 200
)

print(pbmc)


# 4. Quality control

# Calculate mitochondrial gene percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(
  pbmc,
  pattern = "^MT-"
)

# Save QC summary
qc.summary <- data.frame(
  Metric = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  Minimum = c(
    min(pbmc$nFeature_RNA),
    min(pbmc$nCount_RNA),
    min(pbmc$percent.mt)
  ),
  Median = c(
    median(pbmc$nFeature_RNA),
    median(pbmc$nCount_RNA),
    median(pbmc$percent.mt)
  ),
  Mean = c(
    mean(pbmc$nFeature_RNA),
    mean(pbmc$nCount_RNA),
    mean(pbmc$percent.mt)
  ),
  Maximum = c(
    max(pbmc$nFeature_RNA),
    max(pbmc$nCount_RNA),
    max(pbmc$percent.mt)
  )
)

write.csv(
  qc.summary,
  "results/QC_summary.csv",
  row.names = FALSE
)

# Visualize QC metrics
p1 <- VlnPlot(
  pbmc,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0.1
)

ggsave(
  "figures/QC_violin_plots.png",
  p1,
  width = 12,
  height = 5
)

p2 <- FeatureScatter(
  pbmc,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA"
)

ggsave(
  "figures/QC_feature_scatter.png",
  p2,
  width = 6,
  height = 5
)


# 5. Filter low-quality cells

pbmc.filtered <- subset(
  pbmc,
  subset =
    nFeature_RNA > 200 &
    nFeature_RNA < 5000 &
    nCount_RNA < 25000 &
    percent.mt < 10
)

qc.retention <- data.frame(
  Stage = c("Before_QC", "After_QC"),
  Cells = c(ncol(pbmc), ncol(pbmc.filtered))
)

qc.retention$Percentage <- (
  qc.retention$Cells / ncol(pbmc)
) * 100

write.csv(
  qc.retention,
  "results/QC_cell_retention.csv",
  row.names = FALSE
)

print(pbmc.filtered)


# 6. Normalize data

pbmc.filtered <- NormalizeData(
  pbmc.filtered,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)


# 7. Identify highly variable genes

pbmc.filtered <- FindVariableFeatures(
  pbmc.filtered,
  selection.method = "vst",
  nfeatures = 2000
)

top10 <- head(
  VariableFeatures(pbmc.filtered),
  10
)

p3 <- VariableFeaturePlot(pbmc.filtered)

ggsave(
  "figures/Variable_features.png",
  p3,
  width = 8,
  height = 6
)

write.csv(
  data.frame(Gene = VariableFeatures(pbmc.filtered)),
  "results/Variable_genes.csv",
  row.names = FALSE
)


# 8. Scale data and run PCA

pbmc.filtered <- ScaleData(
  pbmc.filtered,
  features = VariableFeatures(pbmc.filtered)
)

pbmc.filtered <- RunPCA(
  pbmc.filtered,
  features = VariableFeatures(pbmc.filtered)
)

p4 <- ElbowPlot(
  pbmc.filtered,
  ndims = 30
)

ggsave(
  "figures/PCA_elbow_plot.png",
  p4,
  width = 7,
  height = 5
)


# 9. Find neighbors and clusters

pbmc.filtered <- FindNeighbors(
  pbmc.filtered,
  dims = 1:20
)

pbmc.filtered <- FindClusters(
  pbmc.filtered,
  resolution = 0.5
)

# Save cluster cell counts
cluster.counts <- as.data.frame(
  table(Idents(pbmc.filtered))
)

colnames(cluster.counts) <- c("Cluster", "Cell_Count")

write.csv(
  cluster.counts,
  "results/Cluster_cell_counts.csv",
  row.names = FALSE
)


# 10. Run UMAP

pbmc.filtered <- RunUMAP(
  pbmc.filtered,
  dims = 1:20
)

p5 <- DimPlot(
  pbmc.filtered,
  reduction = "umap",
  label = TRUE,
  repel = TRUE
) + ggtitle("Human PBMC Clusters")

ggsave(
  "figures/UMAP_clusters.png",
  p5,
  width = 9,
  height = 7
)


# 11. Find cluster marker genes

pbmc.markers <- FindAllMarkers(
  pbmc.filtered,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write.csv(
  pbmc.markers,
  "results/All_cluster_markers.csv",
  row.names = FALSE
)

top.markers <- pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 10,
    with_ties = FALSE
  )

write.csv(
  top.markers,
  "results/Top10_markers_per_cluster.csv",
  row.names = FALSE
)


# 12. Visualize canonical immune markers

canonical.markers <- c(
  "CD3D", "CD3E", "IL7R", "CCR7",
  "CD8A", "NKG7", "GNLY",
  "MS4A1", "CD79A", "CD14",
  "FCGR3A", "FCN1", "LST1",
  "FCER1A", "CLEC10A", "CLEC4C",
  "LILRA4", "PPBP", "PF4"
)

available.markers <- intersect(
  canonical.markers,
  rownames(pbmc.filtered)
)

p6 <- DotPlot(
  pbmc.filtered,
  features = available.markers,
  scale = FALSE
) + RotatedAxis()

ggsave(
  "figures/Canonical_marker_DotPlot.png",
  p6,
  width = 14,
  height = 7
)


# 13. Preliminary cluster annotation

cluster.annotations <- data.frame(
  Cluster = as.character(0:17),
  Preliminary_Identity = c(
    "Classical monocytes",
    "T cells, possibly CD4-like",
    "Activated CD4-like T cells",
    "CD8 T cells",
    "Naive B cells",
    "Inflammatory monocytes",
    "NK cells",
    "Effector/memory T cells",
    "Gamma delta T cells",
    "Mixed T/monocyte-like population",
    "Plasma-like B cells",
    "cDC2-like / myeloid cells",
    "C1Q-positive monocyte-like cells",
    "Mixed T/B marker signature",
    "pDC-like cells",
    "Platelet-like population",
    "Uncharacterized population",
    "Basophil-like / immature granulocyte-like"
  ),
  Confidence = c(
    "High", "Medium", "Medium", "High",
    "High", "High", "High", "Medium",
    "High", "Ambiguous", "Medium", "Medium",
    "Medium", "Ambiguous", "Uncertain",
    "High", "Low", "Low"
  )
)

write.csv(
  cluster.annotations,
  "results/Preliminary_cluster_annotations.csv",
  row.names = FALSE
)


# 14. DoubletFinder analysis

# Estimate expected doublet number using an assumed 5% rate
# This is an assumption, not a measured loading rate.

set.seed(123)

sweep.res <- paramSweep(
  pbmc.filtered,
  PCs = 1:20,
  sct = FALSE
)

sweep.stats <- summarizeSweep(
  sweep.res,
  GT = FALSE
)

bcmvn <- find.pK(sweep.stats)

write.csv(
  bcmvn,
  "results/DoubletFinder_parameter_sweep.csv",
  row.names = FALSE
)

# Use the previously selected pK value
selected.pK <- 0.25

expected.doublets <- round(
  0.05 * ncol(pbmc.filtered)
)

# Run DoubletFinder
pbmc.filtered <- doubletFinder(
  pbmc.filtered,
  PCs = 1:20,
  pN = 0.25,
  pK = selected.pK,
  nExp = expected.doublets,
  reuse.pANN = NULL,
  sct = FALSE
)

# Identify the metadata columns generated by DoubletFinder
df.columns <- grep(
  "^DF.classifications_",
  colnames(pbmc.filtered[[]]),
  value = TRUE
)

df.column <- tail(df.columns, 1)

# Save DoubletFinder summary
df.summary <- as.data.frame(
  table(pbmc.filtered[[]][[df.column]])
)

colnames(df.summary) <- c("Classification", "Cell_Count")

write.csv(
  df.summary,
  "results/DoubletFinder_summary.csv",
  row.names = FALSE
)

# Save cluster-level DoubletFinder summary
df.cluster.summary <- pbmc.filtered[[]] %>%
  mutate(
    Cluster = as.character(seurat_clusters),
    Classification = .data[[df.column]]
  ) %>%
  count(Cluster, Classification) %>%
  tidyr::pivot_wider(
    names_from = Classification,
    values_from = n,
    values_fill = 0
  )

write.csv(
  df.cluster.summary,
  "results/DoubletFinder_by_cluster.csv",
  row.names = FALSE
)


# 15. Create a separate singlet object for exploratory analysis

# Extract metadata from the filtered Seurat object
metadata <- pbmc.filtered[[]]

# Identify Singlet cell barcodes
singlet.cells <- rownames(metadata)[
  metadata[[df.column]] == "Singlet"
]

# Create a separate Seurat object containing Singlet cells
pbmc.singlets.test <- subset(
  pbmc.filtered,
  cells = singlet.cells
)

# Display the Singlet Seurat object
print(pbmc.singlets.test)

# Save singlet cell count
singlet.summary <- data.frame(
  Population = c("All QC-passed cells", "Predicted singlets"),
  Cell_Count = c(
    ncol(pbmc.filtered),
    ncol(pbmc.singlets.test)
  )
)

# Save the singlet summary
write.csv(
  singlet.summary,
  "results/Singlet_cell_summary.csv",
  row.names = FALSE
)

# Save the Singlet Seurat object locally
saveRDS(
  pbmc.singlets.test,
  "pbmc_singlets.rds"
)

# 16. Save final Seurat objects locally

# Save the QC-filtered Seurat object
saveRDS(
  pbmc.filtered,
  "results/pbmc_filtered.rds"
)

# Save the Singlet Seurat object
saveRDS(
  pbmc.singlets.test,
  "results/pbmc_singlets.rds"
)

# Verify that both files were saved successfully
file.exists("results/pbmc_filtered.rds")
file.exists("results/pbmc_singlets.rds")

# 17. Save session information

# Capture R session information
session.info <- capture.output(sessionInfo())

# Save session information to a text file
writeLines(
  session.info,
  "results/sessionInfo.txt"
)

# Check whether the file was saved
file.exists("results/sessionInfo.txt")


# 18. Final report

cat("Analysis completed successfully.\n")
cat("QC-passed cells:", ncol(pbmc.filtered), "\n")
cat("Predicted singlets:", ncol(pbmc.singlets.test), "\n")
cat("Number of clusters:", length(levels(Idents(pbmc.filtered))), "\n")

# ============================================================
# End of analysis
# ============================================================
