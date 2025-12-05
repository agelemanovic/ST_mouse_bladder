# Spatial Transcriptomics Analysis — Mouse Bladder (10x Visium, Seurat v5)

This repository contains the complete R pipeline used to analyze **10x Genomics Visium spatial transcriptomics data** from healthy mouse bladder tissue.  
The workflow includes QC, filtering, SCTransform normalization, clustering, marker detection, visualizations, additional analyses, and optional Loupe Browser export.

Run the analysis with the provided R script "ST_spatial_bladder_pipeline.R".

## Requirements (R package versions used)
R ≥ 4.2  
Seurat 5.3.1  
dplyr 1 1 4  
tidyr 1 3 1  
ggplot2 4.0.0  
ggrepel 0 9 6  
ggsci 4 0 0  
pheatmap 1.0.13  
patchwork 1 3 2  
remotes 2.5.0  
hdf5r  1.3.12  
loupeR 1.1.4  
arrow 21.0.0.1
