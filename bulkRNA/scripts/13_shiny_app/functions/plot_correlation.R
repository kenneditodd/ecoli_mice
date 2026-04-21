# For testing
#deg_list1 <- "E.2D.O.F_vs_S.2D.O.F_FDRq_1.00_LFC_0.00.tsv"
#deg_list2 <- "E.2D.Y.F_vs_S.2D.Y.F_FDRq_1.00_LFC_0.00.tsv"
#fdrq1 <- 0.05
#lfc1 <- 0
#fdrq2 <- NULL
#lfc2 <- NULL

plot_correlation <- function(deg_list1, fdrq1, lfc1, deg_list2, fdrq2 = NULL, lfc2 = NULL) {
  
  # read tables
  data1 <- readr::read_tsv(paste0("data/DEGs/", deg_list1), show_col_types = FALSE)
  data2 <- readr::read_tsv(paste0("data/DEGs/", deg_list2), show_col_types = FALSE)

  # Filter based on criteria - AND/OR
  if (is.null(fdrq2) && is.null(lfc2) && !is.null(fdrq1) && !is.null(lfc1)) {
    # filter by OR
    genes_to_keep <- unique(c(
      data1 %>% filter(adj.P.Val < fdrq1 & (logFC > lfc1 | logFC < -lfc1)) %>% pull(gene_name),
      data2 %>% filter(adj.P.Val < fdrq1 & (logFC > lfc1 | logFC < -lfc1)) %>% pull(gene_name)
    ))
    data1 <- data1 %>% filter(gene_name %in% genes_to_keep)
    data2 <- data2 %>% filter(gene_name %in% genes_to_keep)
  } else if (!is.null(fdrq2) && !is.null(lfc2)) {
    # filter by AND
    data1 <- data1 %>% filter(adj.P.Val < fdrq1 & (logFC > lfc1 | logFC < -lfc1))
    data2 <- data2 %>% filter(adj.P.Val < fdrq2 & (logFC > lfc2 | logFC < -lfc2))
  }
  
  # Find shared genes
  shared_genes <- intersect(data1$gene_name_unique, data2$gene_name_unique)
  
  # Filter shared genes
  data1 <- data1 %>% filter(gene_name_unique %in% shared_genes)
  data2 <- data2 %>% filter(gene_name_unique %in% shared_genes)
  
  shared <- left_join(data1, data2, by = "gene_name_unique", suffix = c(".x", ".y"))
  
  # Check for sufficient data
  if (nrow(shared) < 2) {
    print(length(shared))
    warning("Not enough data points to plot correlation")
    return(NULL)
  }
  
  # Spearman's correlation
  res <- cor.test(shared$logFC.x, shared$logFC.y, method = "spearman")
  p_value <- round(res$p.value, 5)
  rho_value <- round(res$estimate, 3)
  
  # Linear model
  lm_fit <- lm(logFC.y ~ logFC.x, data = shared)
  slope <- coef(lm_fit)[2]
  intercept <- coef(lm_fit)[1]
  
  # Calculate the limits for the plot
  x_limits <- range(shared$logFC.x, na.rm = TRUE)
  y_limits <- range(shared$logFC.y, na.rm = TRUE)
  
  # Create the plot
  p <- ggplot(data = shared, 
              aes(x = logFC.x, y = logFC.y, text = gene_name_unique)) +
    annotate("rect", 
             xmin = 0, xmax = max(shared$logFC.x, na.rm = TRUE), 
             ymin = 0, ymax = max(shared$logFC.y, na.rm = TRUE), 
             fill = "lightpink3", alpha = 0.5) +  
    annotate("rect", 
             xmin = min(shared$logFC.x, na.rm = TRUE), xmax = 0, 
             ymin = min(shared$logFC.y, na.rm = TRUE), ymax = 0, 
             fill = "cadetblue3", alpha = 0.5) +
    geom_point(size = 2) +
    geom_abline(slope = slope, intercept = intercept, color = "gray", 
                size = 0.8, linetype = "dashed") +  
    annotate("text", 
             x = 0, y = 0, 
             label = paste0("rho = ", rho_value), 
             color = "red", size = 5) +
    scale_x_continuous(limits = x_limits) +  
    scale_y_continuous(limits = y_limits) +  
    labs(title = paste("Spearman's rho:", rho_value),
         x = paste0(basename(deg_list1), " LogFC"),
         y = paste0(basename(deg_list2), " LogFC")) +
    theme_minimal()
  
  return(ggplotly(p))
}