# Kennedi Todd
# Jan 3, 2024

# Set .libPaths() to the paths in R_LIBS_USER
lib_paths <- strsplit(Sys.getenv("R_LIBS_USER"), ":")[[1]]
.libPaths(lib_paths)

# Verify the library paths
print(.libPaths())

# set wd
setwd(".")

# load packages
library(dplyr)        # ungroup()
library(gtools)       # smartbind()
library(parallel)     # detectCores()
library(Seurat)       # DimPlot()
library(stringr)      # str_match()
library(tidyr)        # %>%

# variables
samples <- readLines("../../refs/sample_names.txt")
samples <- gtools::mixedsort(samples)
out <- "../../results/all_samples/"

# single cell functions
files <- list.files("../../functions", full.names = TRUE)
invisible(lapply(files, source))

# load object
human.annotated <- readRDS("../../rObjects/seurat_obj_annotated.rds")
human.annotated$group <- str_match(human.annotated$donor, "([HDRRP]+).+")[,2]
human.annotated$group2 <- paste0(human.annotated$group, "_", human.annotated$sex)

############################### DE: RRP vs HD ##################################

print("RRP vs HD")

Idents(human.annotated) <- "Azimuth_L1_predictions"
DE_within_each_cluster(obj = human.annotated,
                       outDir = paste0(out, "DEGs/DEG_tables"),
                       clusterCol = "Azimuth_L1_predictions",
                       groupCol = "group",
                       group1 = "RRP",
                       group2 = "HD")
gc()

############################### DE: RRP_F vs HD_F ##############################

print("RRP_F vs HD_F")

DE_within_each_cluster(obj = human.annotated,
                       outDir = paste0(out, "DEGs/DEG_tables"),
                       clusterCol = "Azimuth_L1_predictions",
                       groupCol = "group2",
                       group1 = "RRP_F",
                       group2 = "HD_F")
gc()

############################### DE: RRP_M vs HD_M ##############################

print("RRP_M vs HD_M")

DE_within_each_cluster(obj = human.annotated,
                       outDir = paste0(out, "DEGs/DEG_tables"),
                       clusterCol = "Azimuth_L1_predictions",
                       groupCol = "group2",
                       group1 = "RRP_M",
                       group2 = "HD_M")
gc()
