# load libraries
library(dotenv)
library(DropletUtils)
library(Seurat)
library(scCustomize)

# source environment variables
load_dot_env(file = "../../refs/.env")
project_dir <- Sys.getenv("PROJECT_DIR")

# read list of samples
samples <- readLines(paste0(project_dir, "/refs/sample_names.txt"))

###################### Copy 10x Genomics MEX files over ########################

# loop through samples and copy matrix files to new folder
for (i in samples) {
  # print the sample you're on
  print(i)

  # set paths and dir names
  input_dir <- paste0(project_dir, "/counts/", i, "/outs/filtered_feature_bc_matrix")
  output_base_dir <- paste0(project_dir, "/matrices")
  output_new_dir_name <- i
  
  # copy barcodes
  file.copy(from = paste0(input_dir, "/barcodes.tsv.gz"), to = output_base_dir)
  file.rename(file.path(output_base_dir, "barcodes.tsv.gz"), 
              paste0(output_base_dir, "/", i, "_barcodes.tsv.gz"))
  
  # copy genes
  file.copy(from = paste0(input_dir, "/features.tsv.gz"), to = output_base_dir)
  file.rename(file.path(output_base_dir, "features.tsv.gz"), 
              paste0(output_base_dir, "/", i, "_features.tsv.gz"))
  
  # copy matrix
  file.copy(from = paste0(input_dir, "/matrix.mtx.gz"), to = output_base_dir)
  file.rename(file.path(output_base_dir, "matrix.mtx.gz"), 
              paste0(output_base_dir, "/", i, "_matrix.mtx.gz"))
  
}
