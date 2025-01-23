# load libraries
library(stringr)
library(dotenv)

# load the environment variables
load_dot_env(file = "../../refs/.env")

# access the environment variables
project_dir <- Sys.getenv("PROJECT_DIR")

# file locations
samples <- c("E3CF1","E3CF2","E3CM1","E3CM2","E3PF1","E3PF2","E3PM1","E3PM2",
             "E4CF1","E4CF2","E4CM1","E4CM2","E4PF1","E4PF2","E4PM1","E4PM2")
locations <- paste0(project_dir, "/counts/", samples, "/outs/metrics_summary.csv")

# initialize df and loop through files
df <- data.frame()
for (i in 1:length(locations)) {
  if (i == 1) {
    df <- read.csv(locations[i])
  } else {
    row <- read.csv(locations[i])[1,]
    df <- rbind(df,row)
  }
}

rownames(df) <- samples
c.names <- c("estimated_cells", "mean_reads", "median_genes", "number_reads",
             "valid_barcodes", "sequencing_saturation", "Q30_bases_barcode",
             "Q30_bases_read", "Q30_bases_UMI", "reads_mapped_genome", "confident_reads_mapped_genome",
             "confident_intergenic_reads_mapped", "confident_intronic_reads_mapped",
             "confident_exonic_reads_mapped", "confident_reads_mapped_transcriptome",
             "reads_mapped_antisense", "fraction_reads", "total_genes", "median_UMI")
colnames(df) <- tolower(gsub("\\.", "_", colnames(df)))

write.table(df, 
            paste0(project_dir, "/counts/web_summaries/overall_metrics.tsv"),
            sep = "\t",
            quote = FALSE)
