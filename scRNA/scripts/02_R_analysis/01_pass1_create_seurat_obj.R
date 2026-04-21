# load libraries
library(dotenv)      # load_dot_env()
library(dplyr)       # left_join()
library(rtracklayer) # import()
library(scCustomize) # Merge_Seurat_List()
library(Seurat)      # CreateSeuratObject()
library(stringr)     # str_match()

# load the environment variables
load_dot_env(file = "../../refs/.env")
ann_dir <- Sys.getenv("ANNOTATION_REFS")

# get sample names
samples <- readLines("../../refs/sample_list.tsv")
samples <- gtools::mixedsort(samples)

# read in annotation file
if (file.exists("../../rObjects/annotation.rds")) {
  genes <- readRDS("../../rObjects/annotation.rds")
} else {
  file <- paste0(ann_dir, "/refdata-gex-GRCm39-2024-A/genes/genes.gtf")
  genes <- rtracklayer::import(file)
  genes <- as.data.frame(genes)
  saveRDS(genes, "../../rObjects/annotation.rds")
}

# read in meta
meta <- readRDS("../../rObjects/meta.rds")

# Initialize an empty list
seurat_obj_list <- list()
  
# create list of individual seurat objects
for (i in 1:length(samples)) {
  print(i)
  sample <- samples[i]
    
  # Create Seurat object with Cell Ranger output
  obj <- CreateSeuratObject(
    Read10X(data.dir = paste0("../../counts/", sample, "/outs/filtered_feature_bc_matrix"))
  )
  
  # Add sample as prefix to cell names
  obj <- RenameCells(obj, add.cell.id = sample)
    
  # Add Seurat object to the list with the sample name as the key
  seurat_obj_list[[sample]] <- obj
    
  # cleanup - helps with memory
  remove(obj)
  gc()
}
  
# Merge all Seurat objects
mouse <- merge(seurat_obj_list[[1]], y = seurat_obj_list[-1])
  
# Set project name
mouse@project.name <- "E.coli Mice scRNAseq"
  
# Join layers
mouse$orig.ident <- colnames(mouse)
mouse <- JoinLayers(mouse)
  
# Extract animal_id
mouse$animal_id <- str_match(colnames(mouse), "[oldyung]+_[fesm]+_([0-9]+)_.+")[,2]
  
# Check
table(mouse$animal_id)
  
# Add meta
mouse@meta.data <- left_join(x = mouse@meta.data, y = meta, by = "animal_id")
rownames(mouse@meta.data) <- mouse$orig.ident
  
# cleanup
remove(meta, seurat_obj_list)
gc()

# preview
mouse

# set idents
Idents(mouse) <- "sample"

# add meta columns
# cell.complexity
mouse$cell_complexity <- log10(mouse$nFeature_RNA) / log10(mouse$nCount_RNA)
# percent.mt
mt.genes <- genes[genes$type == "gene",]
mt.genes <- mt.genes[mt.genes$gene_type == "protein_coding",]
mt.genes <- mt.genes[mt.genes$seqnames == "chrM",]
mt.genes <- mt.genes$gene_name
mouse$percent_mt <- PercentageFeatureSet(mouse, features = mt.genes)
mt.genes
# percent.ribo
# ribosomal proteins begin with 'Rps' or 'Rpl' in this annotation file
# mitochondrial ribosomes start with 'Mrps' or 'Mrpl'
ribo <- genes[genes$type == "gene",]
ribo <- ribo[ribo$gene_type == "protein_coding",]
mt.ribo <- ribo[grep("^Mrp[sl]", ribo$gene_name),"gene_name"]
ribo <- ribo[grep("^Rp[sl]", ribo$gene_name), "gene_name"]
ribo.combined <- c(mt.ribo,ribo)
ribo.combined <- ribo.combined[ribo.combined %in% rownames(mouse)]
mouse$percent_ribo <- PercentageFeatureSet(mouse, features = ribo.combined)
ribo.combined
# percent.hb
hb.genes <- c("Hba-x","Hba-a1","Hba-a2","Hbb-bt","Hbb-bs","Hbb-bh2","Hbb-bh1","Hbb-y")
mouse$percent_hb <- PercentageFeatureSet(mouse, features = hb.genes)
hb.genes

# set sample levels
new_order <- mouse@meta.data %>%
  arrange(timepoint_days, age, treatment, sex, animal_id) %>%
  select(sample)
new_order <- new_order$sample
new_order <- unique(new_order)
mouse$sample <- factor(mouse$sample, levels = new_order)

# reorder cell columns
new_order <- mouse@meta.data %>%
  arrange(timepoint_days, age,treatment, sex, animal_id) %>%
  rownames()
mouse <- mouse[, new_order]

# save with compression
saveRDS(mouse, "../../rObjects/cellranger_merged_before_filtering.rds", compress = FALSE)
