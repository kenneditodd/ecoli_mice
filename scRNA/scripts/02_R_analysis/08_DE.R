# Set working directory
setwd(".")

# Load packages
library(dplyr)
library(gtools)
library(MAST)
library(parallel)
library(Seurat)
library(stringr)
library(tidyr)

# single cell functions
files <- list.files("../../functions", full.names = TRUE)
invisible(lapply(files, source))

# Set output directory
out <- "../../results/without_integration/pass2/"

# Load Seurat object
mouse.annotated <- readRDS("../../rObjects/seurat_obj_annotated_pass2.rds")
DefaultAssay(mouse.annotated) <- "RNA"
Idents(mouse.annotated) <- "annotated_clusters"

################################################################################

# Define comparisons with the Male/Female variable
comparisons_with_sex <- list(
  c("E2DYF", "E2DYM"), c("E2DYF", "S2DYF"), c("E2DYM", "S2DYM"),
  c("E2DOF", "E2DOM"), c("E2DOF", "S2DOF"), c("E2DOM", "S2DOM"),
  c("E18DYF", "E18DYM"), c("E18DYF", "S18DYF"), c("E18DYM", "S18DYM"),
  c("E18DOF", "E18DOM"), c("E18DOF", "S18DOF"), c("E18DOM", "S18DOM"),
  c("E2DOF", "E2DYF"), c("E2DOM", "E2DYM"), c("S2DOF", "S2DYF"), c("S2DOM", "S2DYM"),
  c("E18DOF", "E18DYF"), c("E18DOM", "E18DYM"), c("S18DOF", "S18DYF"), c("S18DOM", "S18DYM"),
  c("E18DYF", "E2DYF"), c("E18DYM", "E2DYM"), c("E18DOF", "E2DOF"), c("E18DOM", "E2DOM"),
  c("S18DYF", "S2DYF"), c("S18DYM", "S2DYM"), c("S18DOF", "S2DOF"), c("S18DOM", "S2DOM"),
  c("S2DYF", "S2DYM"), c("S2DOF", "S2DOM"), c("S18DYF", "S18DYM"), c("S18DOF", "S18DOM")
)

# Loop through comparisons and execute DE
for (comparison in comparisons_with_sex) {
  group1 <- comparison[1]
  group2 <- comparison[2]
  
  message(paste("Comparing", group1, "vs.", group2))
  
  DE_within_each_cluster(
    obj = mouse.annotated,
    outDir = paste0(out, "DEGs/DEG_tables"),
    clusterCol = "annotated_clusters",
    groupCol = "group2",
    group1 = group1,
    group2 = group2
  )
  
  gc()  # Garbage collection to free memory
}

################################################################################

# Define comparisons without the Male/Female variable
comparisons_no_sex <- list(
  c("E2DY", "S2DY"), c("E18DY","S18DY"), c("E2DO","S2DO"), c("E18DO","S18DO"),
  c("E18DY","E2DY"), c("S18DY","S2DY"), c("E18DO","E2DO"), c("S18DO","S2DO"),
  c("E2DO","E2DY"), c("S2DO","S2DY"), c("E18DO","E18DY"), c("S18DO","S18DY")
)

# Loop through comparisons and execute DE with groupCol = "group"
for (comparison in comparisons_no_sex) {
  group1 <- comparison[1]
  group2 <- comparison[2]
  
  message(paste("Comparing", group1, "vs.", group2, "with groupCol = 'group'"))
  
  DE_within_each_cluster(
    obj = mouse.annotated,
    outDir = paste0(out, "DEGs/DEG_tables"),
    clusterCol = "annotated_clusters",
    groupCol = "group",  # Different group column for these comparisons
    group1 = group1,
    group2 = group2
  )
  
  gc()  # Garbage collection to free memory
}

################################################################################

comparisons <- list(c("old","young"))

for (comparison in comparisons) {
  group1 <- comparison[1]
  group2 <- comparison[2]
  
  message(paste("Comparing", group1, "vs.", group2, "with groupCol = 'age'"))
  
  DE_within_each_cluster(
    obj = mouse.annotated,
    outDir = paste0(out, "DEGs/DEG_tables"),
    clusterCol = "annotated_clusters",
    groupCol = "age",  # Different group column for these comparisons
    group1 = group1,
    group2 = group2
  )
  
  gc()  # Garbage collection to free memory
}

################################################################################

comparisons <- list(c("Female","Male"))

for (comparison in comparisons) {
  group1 <- comparison[1]
  group2 <- comparison[2]
  
  message(paste("Comparing", group1, "vs.", group2, "with groupCol = 'sex'"))
  
  DE_within_each_cluster(
    obj = mouse.annotated,
    outDir = paste0(out, "DEGs/DEG_tables"),
    clusterCol = "annotated_clusters",
    groupCol = "sex",  # Different group column for these comparisons
    group1 = group1,
    group2 = group2
  )
  
  gc()  # Garbage collection to free memory
}

################################################################################

comparisons <- list(c("18","2"))

for (comparison in comparisons) {
  group1 <- comparison[1]
  group2 <- comparison[2]
  
  message(paste("Comparing", group1, "vs.", group2, "with groupCol = 'timepoint_days'"))
  
  DE_within_each_cluster(
    obj = mouse.annotated,
    outDir = paste0(out, "DEGs/DEG_tables"),
    clusterCol = "annotated_clusters",
    groupCol = "timepoint_days",
    group1 = group1,
    group2 = group2
  )
  
  gc()  # Garbage collection to free memory
}

################################################################################

comparisons <- list(c("ecoli","saline"))

for (comparison in comparisons) {
  group1 <- comparison[1]
  group2 <- comparison[2]
  
  message(paste("Comparing", group1, "vs.", group2, "with groupCol = 'treatment'"))
  
  DE_within_each_cluster(
    obj = mouse.annotated,
    outDir = paste0(out, "DEGs/DEG_tables"),
    clusterCol = "annotated_clusters",
    groupCol = "treatment",
    group1 = group1,
    group2 = group2
  )
  
  gc()  # Garbage collection to free memory
}
