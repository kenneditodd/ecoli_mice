# load libraries
library(dotenv)
library(dplyr)
library(Seurat)
setwd(".")

# thresholds
nCount.min <- 400
nCount.max <- 20000
nFeature.min <- 250
complexity.cutoff <- 0.8
mt.cutoff <- 10
hb.cutoff <- 1

# load the environment variables
load_dot_env(file = "../../refs/.env")
ann_dir <- Sys.getenv("ANNOTATION_REFS")

# read in annotation file
if (file.exists("../../rObjects/annotation.rds")) {
  genes <- readRDS("../../rObjects/annotation.rds")
} else {
  file <- paste0(ann_dir, "/refdata-gex-GRCm39-2024-A/genes/genes.gtf.gz")
  genes <- rtracklayer::import(file)
  genes <- as.data.frame(genes)
  saveRDS(genes, "../../rObjects/annotation.rds")
}

# load data
mouse <- readRDS("../../rObjects/seurat_obj_before_filtering.rds")

# reorder columns
new_order <- mouse@meta.data %>%
  arrange(treatment, timepoint_days, age, sex, animal_id) %>%
  rownames()
mouse <- mouse[, new_order]

# filter
mouse.filtered <- subset(mouse, subset = (nCount_RNA > nCount.min) &
                           (nCount_RNA < nCount.max) &
                           (nFeature_RNA > nFeature.min) &
                           (cell_complexity > complexity.cutoff) &
                           (percent_mt < mt.cutoff) &
                           (percent_hb < hb.cutoff))

# print cells removed
print(paste0(dim(mouse)[2] - dim(mouse.filtered)[2]," cells removed"))

# filter genes
mouse.filtered <- JoinLayers(mouse.filtered)
counts <- GetAssayData(object = mouse.filtered, layer = "counts")
nonzero <- counts > 0  # produces logical
keep <- Matrix::rowSums(nonzero) >= 10  # sum the true/false
counts.filtered <- counts[keep,]  # keep certain genes

# overwrite mouse.filtered
mouse.filtered <- CreateSeuratObject(counts.filtered, 
                                     meta.data = mouse.filtered@meta.data)

# print features removed
print(paste0(dim(counts)[1] - dim(counts.filtered)[1], " features removed"))

# get mt.genes
gene.names <- genes[genes$type == "gene",]
gene.names <- gene.names[gene.names$gene_type == "protein_coding",]
mt.genes <- gene.names[gene.names$seqnames == "chrM", "gene_name"]

# remove mt.genes
counts <- GetAssayData(object = mouse.filtered, layer = "counts")
keep <- !rownames(counts) %in% mt.genes # false when mt.gene
counts.filtered <- counts[keep,]

# overwrite mouse.filtered
mouse.filtered <- CreateSeuratObject(counts.filtered,
                                     meta.data = mouse.filtered@meta.data)

# print features removed
print(paste0(dim(counts)[1] - dim(counts.filtered)[1], " features removed"))

# cleanup data
remove(mouse,counts,counts.filtered,nonzero)

# save
saveRDS(mouse.filtered, "../../rObjects/seurat_obj_filtered.rds")
