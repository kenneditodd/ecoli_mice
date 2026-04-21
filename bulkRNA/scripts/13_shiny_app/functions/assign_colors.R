assign_colors <- function(data, fdrq, lfc) {
  data %>%
    mutate(color_adjpval = case_when(
      adj.P.Val < fdrq & logFC > lfc ~ "up-regulated",
      adj.P.Val < fdrq & logFC < -lfc ~ "down-regulated",
      TRUE ~ "not significant"
    ))
}
