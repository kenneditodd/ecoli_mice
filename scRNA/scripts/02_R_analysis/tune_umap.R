# Load libraries
library(ggplot2)
library(Seurat)
library(gridExtra)

# Set working directory
setwd(".")

# Source thresholds and output paths
source("../../refs/thresholds_and_outs.R")
out <- paste0(out, "pass2/")

# Load integrated Seurat object
mouse.merged <- readRDS(paste0("../../rObjects/", filtering_method, "_pass2_harmony_seurat_obj.rds"))

# Parse SLURM-passed arguments
args <- commandArgs(trailingOnly = TRUE)
dm <- as.numeric(args[1])       # Number of dimensions
md <- as.numeric(args[2])       # min.dist
nn <- as.numeric(args[3])       # n.neighbors

# FeaturePlot markers
marker_genes <- c("Tmem119", "Mrc1", "Slc1a3", 
                  "Flt1", "Pdgfrb", "Snap25", 
                  "Col1a2", "Tnr", "Mag", 
                  "Ttr", "Dcn")


# Run UMAP with specified parameters
mouse.merged <- RunUMAP(
  object = mouse.merged,
  dims = 1:dm,
  reduction = "harmony",
  min.dist = md,
  n.neighbors = nn,
  n.components = 3,
  seed.use = 42,
  verbose = FALSE
)

# Create DimPlot
u <- DimPlot(mouse.merged, group.by = "group", shuffle = TRUE) +
  NoLegend() +
  ggtitle(paste0("Dims: 1-", dm, ", Min.dist: ", md, ", nn: ", nn))

# Create 11 FeaturePlots
fplots <- lapply(marker_genes, function(gene) {
  FeaturePlot(mouse.merged,
              reduction = "umap",
              features = gene) +
    scale_colour_gradientn(colours = c("blue", "lightblue", "yellow", "orange", "red")) +
    ggtitle(gene) +
    theme(legend.position = "none")
})

# Combine all plots (1 DimPlot + 11 FeaturePlots = 12 total)
all_plots <- c(list(u), fplots)
layout <- matrix(1:12, nrow = 3, ncol = 4, byrow = TRUE)

# Save as PDF
path <- paste0(out, "clustering_QC/tune_umap/UMAP_dims", dm, "_md", md, "_nn", nn, ".pdf")
pdf(path, width = 12, height = 9)
grid.arrange(grobs = all_plots, layout_matrix = layout)
dev.off()
