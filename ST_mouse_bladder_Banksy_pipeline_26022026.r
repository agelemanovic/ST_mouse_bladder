
set.seed(1234)

# ------------------------------
# USER INPUT
# ------------------------------
install.packages("BiocManager")
BiocManager::install(c(
  "Banksy",
  "SpatialExperiment",
  "scuttle",
  "scater",
  "scran"
))

#!/usr/bin/env Rscript

############################################################
# BANKSY Nucleus-Level Spatial Clustering Pipeline
############################################################

library(Banksy)
library(SpatialExperiment)
library(scuttle)
library(scater)
library(scran)
library(ggplot2)
counts_file  <- "banksy_counts.csv"
coords_file  <- "banksy_coords.csv"
output_dir   <- "banksy_output"

k_geom       <- 25
alpha        <- 0.15
resolution   <- 0.55
npcs         <- 40

dir.create(output_dir, showWarnings = FALSE)

# ------------------------------
# Load data
# ------------------------------
counts <- as.matrix(read.csv(counts_file, row.names = 1))
coords <- as.matrix(read.csv(coords_file, row.names = 1))

# Remove nuclei with zero transcripts
keep <- colSums(counts) > 0
counts <- counts[, keep]
coords <- coords[keep, , drop = FALSE]

# ------------------------------
# Create SpatialExperiment
# ------------------------------
se <- SpatialExperiment(
  assays = list(counts = counts),
  spatialCoords = coords
)

# ------------------------------
# Normalization
# ------------------------------
se <- computeLibraryFactors(se)
assay(se, "normcounts") <- normalizeCounts(se, log = FALSE)

# ------------------------------
# BANKSY smoothing
# ------------------------------
se <- computeBanksy(
  se,
  assay_name = "normcounts",
  k_geom = k_geom,
  alpha = alpha
)

# ------------------------------
# PCA
# ------------------------------
se <- runBanksyPCA(se, npcs = npcs)

# ------------------------------
# Clustering
# ------------------------------
se <- clusterBanksy(
  se,
  resolution = resolution
)

cluster_col <- grep("clust_", colnames(colData(se)), value = TRUE)[1]
colData(se)$BANKSY_cluster <- colData(se)[[cluster_col]]

# ------------------------------
# UMAP
# ------------------------------
se <- runBanksyUMAP(se)

pdf(file.path(output_dir, "BANKSY_UMAP.pdf"), width = 6, height = 5)
print(
  plotUMAP(se, colour_by = "BANKSY_cluster") +
    theme_classic()
)
dev.off()

# ------------------------------
# Marker detection
# ------------------------------
markers <- findMarkers(
  assay(se, "normcounts"),
  groups = colData(se)$BANKSY_cluster
)

top10 <- do.call(rbind, lapply(names(markers), function(cl) {
  df <- as.data.frame(markers[[cl]])
  df$gene <- rownames(df)
  df$cluster <- cl
  head(df[order(df$Top), ], 10)
}))

write.csv(
  top10,
  file.path(output_dir, "BANKSY_top10_markers.csv"),
  row.names = FALSE
)

saveRDS(se, file.path(output_dir, "BANKSY_object.rds"))

message("BANKSY pipeline finished successfully.")

