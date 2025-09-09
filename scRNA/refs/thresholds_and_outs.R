# variables
filtering_method <- "thresh1_cellbender"
out <- paste0("../../results/", filtering_method, "/")

# thresholds
if (filtering_method == "thresh1_cellbender") {
  nCount.min <- 2500
  nCount.max <- 28000
  nFeature.min <- 1200
  complexity.cutoff <- 0.8
  mt.cutoff <- 10
  hb.cutoff <- 1
}
