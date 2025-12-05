################################################################################

############ Spatial transcriptomics analysis pipeline using Seurat ############


################################################################################
## 0 | Setup
################################################################################

# Install packages
install.packages(c("Seurat", "ggplot2", "dplyr", "patchwork",
                   "arrow", "png", "viridis", "ggrepel", "ggsci",
                   "tidyr", "pheatmap", "hdf5r", "remotes"))
  
# LoupeR must be installed manually
os <- sub("Darwin", "macOS", Sys.info()["sysname"])
url <- paste0("https://github.com/10XGenomics/loupeR/releases/latest/download/loupeR_", os, ".tar.gz")
install.packages(url, repos = NULL, type = "source")

# Load packages
packages <- c("Seurat",
              "ggplot2",
              "dplyr",
              "patchwork",
              "arrow",
              "png",
              "viridis",
              "ggrepel",
              "ggsci",
              "hdf5r",
              "loupeR",
              "tidyr",
              "pheatmap")
invisible(lapply(packages, library, character.only = TRUE))

# Set seed for reproducibility
set.seed(1234)


### Pipeline parameters

# Create output folder
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# User-defined paths
data_dir   <- "path/to/visium/data"       
image_dir  <- "path/to/visium/data/spatial"
output_dir <- "path/to/output"

# Analysis parameters
min_nCount_Spatial <- 40
max_percent_mt     <- 10
min_cells_expr     <- 100
sct_method         <- "glmGamPoi"
sct_assay          <- "Spatial"
clust_resolution   <- 0.6
pca_dims           <- 1:8
nCount_field       <- "nCount_Spatial"
nFeature_field     <- "nFeature_Spatial"
min.pct            <- 0.3
logfc.threshold    <- 0.25


################################################################################
## 1 | Load Visium data
################################################################################

seurat_obj <- Load10X_Spatial(data.dir = data_dir,
                              image = "hires",
                              filter.matrix = TRUE)
  
print(seurat_obj)

# OPTIONAL: Load hires image if available
if (file.exists(file.path(image_dir, "tissue_hires_image.png"))) {
  hires_image <- Read10X_Image(image.dir = image_dir,
                               image.name = "tissue_hires_image.png",
                               assay = "Spatial",
                               filter.matrix = TRUE)
  hires_image <- hires_image[Cells(seurat_obj)]
  seurat_obj[["slice1"]] <- hires_image
}

################################################################################
## 2 | QC & Filtering
################################################################################

# Select cells
seurat_obj_filtered <- subset(seurat_obj, cells = Cells(seurat_obj))

# Mitochondrial percentage
seurat_obj_filtered[["percent.mt"]] <- PercentageFeatureSet(seurat_obj_filtered,
                                                            pattern = "^mt-")
  

# Compute nFeature and nCount manually
seurat_obj_filtered[[nFeature_field]] <- Matrix::colSums(seurat_obj_filtered@assays$Spatial$counts > 0)
seurat_obj_filtered[[nCount_field]] <- Matrix::colSums(seurat_obj_filtered@assays$Spatial$counts)

# Filter bins by defined parameters
seurat_obj_filtered <- subset(seurat_obj_filtered,
                              subset = nCount_Spatial > min_nCount_Spatial & 
                                percent.mt < max_percent_mt)

print(seurat_obj_filtered)

# Remove lowly-expressed genes
gene_counts_nonzero <- rowSums(seurat_obj_filtered@assays$Spatial$counts > 0)
seurat_obj_filtered <- seurat_obj_filtered[gene_counts_nonzero >= min_cells_expr, ]

print(seurat_obj_filtered)


################################################################################
## 3 | Normalization & Dimensionality reduction
################################################################################

seurat_obj_filtered <- SCTransform(seurat_obj_filtered,
                                   assay   = sct_assay,
                                   method  = sct_method,
                                   verbose = TRUE)
  
seurat_obj_filtered <- RunPCA(seurat_obj_filtered, 
                              assay = "SCT")
ElbowPlot(seurat_obj_filtered, ndims = 15)

seurat_obj_filtered <- FindNeighbors(seurat_obj_filtered, 
                                     dims = pca_dims)

seurat_obj_filtered <- FindClusters(seurat_obj_filtered, 
                                    resolution = clust_resolution)

seurat_obj_filtered <- RunUMAP(seurat_obj_filtered, 
                               dims = pca_dims)


################################################################################
## 4 | Cluster renaming / ordering
################################################################################

# Custom order
new_order <- c(
  "1",   # Fibroblasts
  "0",   # Smooth muscle A
  "2",   # Smooth muscle B
  "4",   # Smooth muscle C
  "11",  # Smooth muscle D
  "7",   # Stromal cells (Apoe+)
  "8",   # Stromal cells (Cd74+)
  "3",   # Stromal cells (Dcn+)
  "6",   # Stromal cells (Mgp+)
  "5",   # Urothelium - Basal cells
  "12",  # Urothelium - Intermediate cells
  "10",  # Urothelium - Superficial cells
  "9"    # Urothelium - Umbrella cells
)

seurat_obj_filtered$seurat_clusters <- factor(as.character(seurat_obj_filtered$seurat_clusters),
                                              levels = new_order,
                                              ordered = TRUE)
Idents(seurat_obj_filtered) <- seurat_obj_filtered$seurat_clusters

# Custom naming
cluster_names <- c("1"  = "Fibroblasts (Col1a1+)",
                   "0"  = "Smooth muscle A",
                   "2"  = "Smooth muscle B",
                   "4"  = "Smooth muscle C",
                   "11" = "Smooth muscle D",
                   "7"  = "Stromal cells (Apoe+)",
                   "8"  = "Stromal cells (Cd74+)",
                   "3"  = "Stromal cells (Dcn+)",
                   "6"  = "Stromal cells (Mgp+)",
                   "5"  = "Urothelium - Basal cells",
                   "12" = "Urothelium - Intermediate cells",
                   "10" = "Urothelium - Superficial cells",
                   "9"  = "Urothelium - Umbrella cells")

seurat_obj_filtered <- RenameIdents(seurat_obj_filtered, 
                                    cluster_names)

# Preserve order after renaming
Idents(seurat_obj_filtered) <- factor(Idents(seurat_obj_filtered), 
                                      levels = cluster_names, 
                                      ordered = TRUE)

print("Final cluster order:")
print(levels(seurat_obj_filtered))


################################################################################
## 5 | Marker genes
################################################################################

# Find all marker genes per cluster
markers <- FindAllMarkers(seurat_obj_filtered,
                          assay = "SCT",
                          only.pos = TRUE,
                          min.pct = min.pct,
                          logfc.threshold = logfc.threshold)
write.csv(markers, file.path(output_dir, "all_markers.csv"))

# Top 10 per cluster
top10 <- markers %>% 
  group_by(cluster) %>% 
  filter(pct.1 > 0.3, p_val_adj < 0.05) %>%
  slice_max(avg_log2FC, n = 10)
write.csv(top10, file.path(output_dir, "top10_markers.csv"))

# Top 1 per cluster
top_per_cluster <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = pct.1, n = 1)
write.csv(top_per_cluster, file.path(output_dir, "top1_markers.csv"))


################################################################################
## 6 | Visualization
################################################################################

# Consistent color palette
cluster_colors <- ggsci::pal_d3("category20")(length(levels(seurat_obj_filtered)))
names(cluster_colors) <- levels(seurat_obj_filtered)


# DimPlot
p_umap <- DimPlot(seurat_obj_filtered, 
                  reduction = "umap", 
                  label = TRUE)


# Spatial DimPlot
p_spatial <- SpatialDimPlot(seurat_obj_filtered,
                            image.scale = "hires",
                            image.alpha = 0.6,
                            alpha = c(0.7, 0.7),
                            shape = 21,
                            pt.size.factor = 2.5,
                            label = FALSE,
                            cols = cluster_colors) +
  labs(fill = "") +
  theme(legend.text = element_text(size = 14),
        legend.key.size = unit(1.2, "cm")) +
  guides(fill = guide_legend(override.aes = list(size = 6)))


# DotPlot
p_dot <- DotPlot(seurat_obj_filtered, 
                 features = unique(top10$gene), 
                 assay = "SCT",
                 cols = c("lightgrey", "blue"), 
                 col.min = -2.5, 
                 col.max = 2.5, 
                 dot.min = 0, 
                 dot.scale = 6, 
                 cluster.idents = FALSE) + 
  RotatedAxis()


# Heatmap
p_heatmap <- DoHeatmap(seurat_obj_filtered, 
                       features = top10$gene) + 
  scale_fill_gradientn(colors = viridis::viridis(10)) +
  NoLegend()


# FeaturePlot
p_feature_plot <- FeaturePlot(seurat_obj_filtered, 
                              features = top_per_cluster$gene, 
                              reduction = "umap", 
                              label = TRUE)


# SpatialFeaturePlot
p_spatial_feature_plot <- SpatialFeaturePlot(seurat_obj_filtered,
                                             features = top_per_cluster$gene,
                                             image.alpha = 0.6,
                                             image.scale = "hires") +
  theme(legend.text = element_text(size = 14),
        legend.key.size = unit(1.2, "cm"))
    
  
# Save plots
ggsave(file.path(output_dir, "UMAP_clusters.png"), p_umap, dpi = 300)
ggsave(file.path(output_dir, "Spatial_clusters.png"), p_spatial, dpi = 300)
ggsave(file.path(output_dir, "DotPlot_top10.png"), p_dot, dpi = 300)
ggsave(file.path(output_dir, "Heatmap_top10.png"), p_heatmap, dpi = 300)
ggsave(file.path(output_dir, "FeaturePlot_top1.png"), p_feature_plot, dpi = 300)
ggsave(file.path(output_dir, "SpatialFeaturePlot_top1.png"), p_spatial_feature_plot, dpi = 300)



################################################################################
## 7 | ADDITIONAL ANALYSES
################################################################################

### Module scoring

fibroblast_markers      <- c("Col1a1", "Dcn", "Lum", "Pdgfra", "Col3a1", "Car3")
umbrella_cells_markers  <- c("Upk1a", "Upk1b", "Upk2", "Upk3a", "Krt20")
basal_cells_markers     <- c("Krt5", "Shh", "Trp63")
smooth_muscle_markers   <- c("Acta2", "Myh11", "Tagln", "Myl9", "Cnn1", "Mylk", "Actg2")
endothelial_markers     <- c("Pecam1", "Cdh5", "Vwf", "Kdr", "Tek")
pericytes               <- c("Kcnj8", "Pdgfrb", "Rgs5")
vSMC                    <- c("Tesc", "Pln", "Wtip")
immune_markers2         <- c("Ptprc","Epor","Pecam","Itgam","Ly","Adgre","Cd","Cd","Nos","Cd",
                             "H","Mrc","Arg","Chil","Cd","Itgax","Cd","Thbd","Cd","Mcpt",
                             "Mcpt","Mcpt","Cpa","Mcpt","Tpsb","Kit","Ly","Cd","Itgav","Pdgfra",
                             "Mpo","Csf","Cd","Klrb","Ncr","Ncam","Fcgr","Cd","Cd","Cd","Foxp",
                             "Trdc","Klrk","Cd","Ms","Sdc")

seurat_obj_filtered <- AddModuleScore(seurat_obj_filtered,
  features = list(
    fibroblast_markers,
    umbrella_cells_markers,
    basal_cells_markers,
    smooth_muscle_markers,
    endothelial_markers,
    pericytes,
    vSMC,
    immune_markers2
  ),
  name = c("Fibroblast", "Umbrella_cells", "Basal_cells", "SmoothMuscle",
           "Endothelial", "Pericytes", "vSMC", "Immune"),
  nbin = 15
)
  
cluster_scores <- seurat_obj_filtered@meta.data %>%
  group_by(seurat_clusters) %>%
  summarize(
    Fibroblast = mean(Fibroblast1),
    Umbrella   = mean(Umbrella_cells2),
    Basal      = mean(Basal_cells3),
    Smooth     = mean(SmoothMuscle4),
    Endo       = mean(Endothelial5),
    Pericytes  = mean(Pericytes6),
    vSMC       = mean(vSMC7),
    Immune     = mean(Immune8)
  )

write.csv(cluster_scores, file.path(output_dir, "cluster_scores_module_scores.csv"), row.names = FALSE)

############################################################

### Make pairwise comparison between each group Smooth muscle / Urothelium

# List of cluster groups you want to compare
group_list <- list(
  Smooth_muscle = c(0, 2, 4, 11),
  Urothelium = c(5, 9, 10, 12),
)

# Function to run pairwise FindMarkers within each group

library(dplyr)
library(purrr)

pairwise_markers_within_group <- function(seurat_obj, clusters,
                                          assay = "SCT",
                                          only_pos = TRUE,
                                          min_pct = 0.2,
                                          logfc_threshold = 0.25) {
  # Generate all unique pairs of clusters
  cluster_pairs <- combn(clusters, 2, simplify = FALSE)
  
  # List to store results
  results_list <- list()
  
  for (pair in cluster_pairs) {
    cluster1 <- pair[1]
    cluster2 <- pair[2]
    
    # Subset Seurat object to only those clusters
    cells_use <- WhichCells(seurat_obj, idents = c(cluster1, cluster2))
    seurat_sub <- subset(seurat_obj, cells = cells_use)
    
    # Set Idents to clusters in subset
    Idents(seurat_sub) <- seurat_sub$seurat_clusters
    
    # Run FindMarkers cluster1 vs cluster2
    markers_1_vs_2 <- FindMarkers(
      seurat_sub,
      ident.1 = cluster1,
      ident.2 = cluster2,
      assay = assay,
      only.pos = only_pos,
      min.pct = min_pct,
      logfc.threshold = logfc_threshold
    ) %>%
      rownames_to_column("gene") %>%
      mutate(
        comparison = paste0(cluster1, "_vs_", cluster2),
        cluster1 = cluster1,
        cluster2 = cluster2,
        direction = "cluster1_up"
      )
    
    # Run FindMarkers cluster2 vs cluster1
    markers_2_vs_1 <- FindMarkers(
      seurat_sub,
      ident.1 = cluster2,
      ident.2 = cluster1,
      assay = assay,
      only.pos = only_pos,
      min.pct = min_pct,
      logfc.threshold = logfc_threshold
    ) %>%
      rownames_to_column("gene") %>%
      mutate(
        comparison = paste0(cluster2, "_vs_", cluster1),
        cluster1 = cluster2,
        cluster2 = cluster1,
        direction = "cluster2_up"
      )
    
    # Combine both directions
    results_list[[paste0(cluster1, "_vs_", cluster2)]] <- bind_rows(markers_1_vs_2, markers_2_vs_1)
  }
  
  # Bind all pairs into one data frame
  all_results <- bind_rows(results_list, .id = "pair_id")
  
  return(all_results)
}


# Run this function on your groups
all_pairwise_results <- lapply(names(group_list), function(group_name) {
  clusters <- group_list[[group_name]]
  df <- pairwise_markers_within_group(seurat_obj_filtered, clusters)
  df$group <- group_name
  return(df)
})

# Combine all groups results
all_pairwise_results_df <- bind_rows(all_pairwise_results)

# View
head(all_pairwise_results_df)



library(dplyr)

# Extract top 10 per group with filtering
top10_markers_pairwise <- all_pairwise_results_df %>%
  dplyr::filter(pct.1 > 0.4, p_val_adj < 0.05) %>%
  group_by(group, comparison) %>%
  dplyr::top_n(n = 10, wt = avg_log2FC) %>%
  ungroup() %>%
  as.data.frame()

head(top10_markers_pairwise)

# Extract top 1 per group with filtering
top1_markers_pairwise <- all_pairwise_results_df %>%
  group_by(group, comparison) %>%
  top_n(n = 1, wt = pct.1) %>%
  ungroup() %>%
  as.data.frame()

head(top1_markers_pairwise)


# Save all created dataframes
write.xlsx(all_pairwise_results_df, file.path(output_dir, "DE_markers_all_pairwiseWithinGroup.xlsx"), asTable = F, overwrite = T)
write.xlsx(top10_markers_pairwise, file.path(output_dir, "DE_markers_top10_pairwiseWithinGroup.xlsx"), asTable = F, overwrite = T)
write.xlsx(top1_markers_pairwise, file.path(output_dir,"DE_markers_top1_pairwiseWithinGroup.xlsx"), asTable = F, overwrite = T)

############################################################

### Visualize pairwise results as heatmap

# Load input table (EDIT PATH)

df <- read.csv(
  # You will get this by exporting as.csv from our Supplementary Tables sheet S3
  "path/to/your/Supplementary_TableS3.csv",    
  sep = ";",
  check.names = FALSE)
  
tbl <- as.table(as.matrix(df))

# Make sure comparison is a factor
df$comparison <- as.factor(df$comparison)

df$comparison <- as.factor(df$comparison)
df$neg_log10_p <- -log10(df$p_val_adj)

# Clean dataframe (remove NA / Inf)
df_clean <- df %>%
  filter(!is.na(avg_log2FC),
         !is.na(p_val_adj),
         !is.na(gene),
         gene != "",
         !is.infinite(avg_log2FC),
         !is.infinite(neg_log10_p))

############################################################

## HEATMAP 1 — Smooth Muscle Clusters

smooth_muscle_clusters <- c(0, 2, 4, 11)

# Filter
df_sm <- df_clean %>%
  filter(cluster1 %in% smooth_muscle_clusters & cluster2 %in% smooth_muscle_clusters)

# Define cluster names
smooth_muscle_names <- c(
  "0"  = "Smooth muscle A",
  "2"  = "Smooth muscle B",
  "4"  = "Smooth muscle C",
  "11" = "Smooth muscle D"
)

# Create readable comparison names
df_sm <- df_sm %>%
  mutate(
    comparison_name = case_when(
      comparison == "0_vs_2"  ~ paste(smooth_muscle_names["0"], "vs", smooth_muscle_names["2"]),
      comparison == "0_vs_4"  ~ paste(smooth_muscle_names["0"], "vs", smooth_muscle_names["4"]),
      comparison == "0_vs_11" ~ paste(smooth_muscle_names["0"], "vs", smooth_muscle_names["11"]),
      comparison == "2_vs_4"  ~ paste(smooth_muscle_names["2"], "vs", smooth_muscle_names["4"]),
      comparison == "2_vs_11" ~ paste(smooth_muscle_names["2"], "vs", smooth_muscle_names["11"]),
      comparison == "4_vs_11" ~ paste(smooth_muscle_names["4"], "vs", smooth_muscle_names["11"]),
      comparison == "11_vs_0" ~ paste(smooth_muscle_names["11"], "vs", smooth_muscle_names["0"]),
      comparison == "11_vs_2" ~ paste(smooth_muscle_names["11"], "vs", smooth_muscle_names["2"]),
      comparison == "11_vs_4" ~ paste(smooth_muscle_names["11"], "vs", smooth_muscle_names["4"]),
      comparison == "4_vs_0"  ~ paste(smooth_muscle_names["4"], "vs", smooth_muscle_names["0"]),
      comparison == "4_vs_2"  ~ paste(smooth_muscle_names["4"], "vs", smooth_muscle_names["2"]),
      TRUE ~ comparison
    )
  )

# Significant genes
signif_genes_sm <- df_sm %>%
  filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25) %>%
  pull(gene) %>%
  unique()

# Pivot to heatmap matrix
heatmap_data_sm <- df_sm %>%
  select(gene, comparison_name, avg_log2FC) %>%
  pivot_wider(
    names_from = comparison_name,
    values_from = avg_log2FC,
    values_fill = 0
  )

heatmap_matrix_sm <- as.matrix(heatmap_data_sm[,-1])
rownames(heatmap_matrix_sm) <- heatmap_data_sm$gene

# Keep only significant genes
heatmap_matrix_sm <- heatmap_matrix_sm[rownames(heatmap_matrix_sm) %in% signif_genes_sm, ]

# Draw heatmap
pheatmap(
  heatmap_matrix_sm,
  scale = "row",
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  color = colorRampPalette(c("blue", "white", "firebrick3"))(50),
  main = "Smooth Muscle – Differential Expression Heatmap"
)


############################################################

## HEATMAP 2 — Urothelial Clusters

uro_clusters <- c(5, 9, 10, 12)

# Filter
df_uro <- df_clean %>%
  filter(cluster1 %in% uro_clusters & cluster2 %in% uro_clusters)

# Define cluster names
uro_names <- c(
  "5"  = "Urothelium - Basal",
  "9"  = "Urothelium - Umbrella",
  "10" = "Urothelium - Superficial",
  "12" = "Urothelium - Intermediate"
)

# Create readable comparison names
df_uro <- df_uro %>%
  mutate(
    comparison_name = case_when(
      comparison == "5_vs_9"  ~ paste(uro_names["5"],  "vs", uro_names["9"]),
      comparison == "9_vs_5"  ~ paste(uro_names["9"],  "vs", uro_names["5"]),
      comparison == "5_vs_10" ~ paste(uro_names["5"],  "vs", uro_names["10"]),
      comparison == "10_vs_5" ~ paste(uro_names["10"], "vs", uro_names["5"]),
      comparison == "5_vs_12" ~ paste(uro_names["5"],  "vs", uro_names["12"]),
      comparison == "12_vs_5" ~ paste(uro_names["12"], "vs", uro_names["5"]),
      comparison == "9_vs_10" ~ paste(uro_names["9"],  "vs", uro_names["10"]),
      comparison == "10_vs_9" ~ paste(uro_names["10"], "vs", uro_names["9"]),
      comparison == "9_vs_12" ~ paste(uro_names["9"],  "vs", uro_names["12"]),
      comparison == "12_vs_9" ~ paste(uro_names["12"], "vs", uro_names["9"]),
      comparison == "10_vs_12" ~ paste(uro_names["10"], "vs", uro_names["12"]),
      comparison == "12_vs_10" ~ paste(uro_names["12"], "vs", uro_names["10"]),
      TRUE ~ comparison
    )
  )

# Significant genes
signif_genes_uro <- df_uro %>%
  filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25) %>%
  pull(gene) %>%
  unique()

# Pivot to matrix
heatmap_data_uro <- df_uro %>%
  select(gene, comparison_name, avg_log2FC) %>%
  pivot_wider(
    names_from = comparison_name,
    values_from = avg_log2FC,
    values_fill = 0
  )

heatmap_matrix_uro <- as.matrix(heatmap_data_uro[,-1])
rownames(heatmap_matrix_uro) <- heatmap_data_uro$gene

heatmap_matrix_uro <- heatmap_matrix_uro[rownames(heatmap_matrix_uro) %in% signif_genes_uro, ]

# Draw heatmap
pheatmap(
  heatmap_matrix_uro,
  scale = "row",
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Urothelium – Differential Expression Heatmap"
)

################################################################################
## 8| OPTIONAL: Loupe export
################################################################################
assay <- seurat_obj[["Spatial"]]

# You get this by exporting bins from your raw Loupe file
feature_ids <- read_feature_ids_from_tsv("path/to/your/tissue_positions.csv") 
counts <- counts_matrix_from_assay(assay)

projections <- select_projections(seurat_obj_filtered)

count_mat <- seurat_obj_filtered@assays$SCT@counts

clusters <- select_clusters(seurat_obj_filtered)

create_loupe(count_mat = count_mat,                                          # you will get separate loupe file with the same UMAP as in Seurat  
             clusters = clusters,                                            # when you generate Loupe file - 1. export this UMAP projection as .csv
             projections = projections,                                      # 2. export clusters from Loupe as .csv
             output_name = "bladder_loupe_<dimensions>_<resolution>_<date>") # 3.import both .csv files in appropriate field in your raw Loupe browser file
# Your Seurat clusters will now be presented spatially the same both in R and in Loupe browser

################################################################################
## 9 | END
################################################################################

message("Pipeline finished. Results saved to: ", normalizePath(output_dir))

############################################################