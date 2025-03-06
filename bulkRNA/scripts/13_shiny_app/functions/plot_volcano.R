source("functions/assign_colors.R")

# For testing
#deg_table <- "E.2D.O.F_vs_E.2D.O.M_FDRq_1.00_LFC_0.00.tsv"
#fdrq <- 0.05
#lfc <- 0
#labels <- NULL
#dynamic <- TRUE

plot_volcano <- function(fdrq, lfc, deg_table, labels = NULL, dynamic = FALSE) {
  path <- paste0("data/DEGs/", deg_table)
  data <- readr::read_tsv(file = path, show_col_types = FALSE)
  
  # Ensure gene_name column exists and is character
  if (!"gene_name" %in% colnames(data)) {
    stop("The 'gene_name' column does not exist in the dataset.")
  }
  data$gene_name <- as.character(data$gene_name)
  
  # Assign colors based on FDRq and logFC
  data <- assign_colors(data, fdrq, lfc)
  
  # Select top 50 genes
  passed_genes <- data %>% filter(adj.P.Val < fdrq & (logFC > lfc | logFC < -lfc))
  top_genes <- rbind(passed_genes %>% arrange(adj.P.Val) %>% head(20),
                     passed_genes %>% arrange(logFC) %>% head(10),
                     passed_genes %>% arrange(-logFC) %>% head(10))
  top_genes <- unique(top_genes)
  
  hadjpval <- (-log10(max(data$P.Value[data$adj.P.Val < fdrq], na.rm=TRUE)))
  
  # Define plot
  p <- ggplot(data, aes(x = logFC, y = -log10(adj.P.Val), color = color_adjpval)) +
    geom_point() +
    scale_color_manual(values = c("up-regulated" = "red", "down-regulated" = "blue", "not significant" = "gray")) +
    labs(x = "log2(Fold Change)", y = "-log10(FDRq)", title = gsub("_FDRq_1.00_LFC_0.00.tsv","",deg_table)) +
    geom_hline(yintercept = hadjpval,  #  horizontal line
               colour = "#000000",
               linetype = "dashed") +
    geom_vline(xintercept = lfc,  #  vertical line
               colour = "#000000",
               linetype = "dashed") +
    geom_vline(xintercept = -lfc,  #  vertical line
               colour = "#000000",
               linetype = "dashed") +
    theme_minimal()+
    theme(
      axis.text = element_text(size = 15, color = "black"),
      axis.title = element_text(size = 15, color = "black"),
      plot.title = element_text(size = 20, color = "black", face = "bold"),
      legend.title = element_text(size = 15, color = "black"),
      legend.text = element_text(size = 15, color = "black"))
  
  # Always label top genes in static mode
  if (!dynamic && nrow(top_genes) > 0) {
    p <- p + geom_text_repel(data = top_genes, aes(label = gene_name), size = 4, color = "black")
  } 
  
  # Label user-selected gene only if provided
  if (!dynamic && !is.null(labels) && labels %in% data$gene_name) {
    p <- p + geom_text_repel(data = data %>% filter(gene_name == labels), aes(label = gene_name), size = 4, color = "black")
  }
  
  if (dynamic) {
    return(ggplotly(p))  # Convert ggplot to interactive plotly object
  } else {
    return(p)  # Return static ggplot object
  }
}
