plot_boxplot <- function(counts, gene) {
  
  # Create df
  names <- colnames(counts)
  group <- stringr::str_match(names, "([ES\\.]+\\.[21D]+\\.[OY]+\\.[MF])")[,2]
  value <- as.numeric(subset(counts, rownames(counts) == gene))
  df <- data.frame(sample = names,
                   group = factor(group, levels = c("S.2D.Y.M","S.2D.Y.F","E.2D.Y.M","E.2D.Y.F",
                                                    "S.2D.O.M","S.2D.O.F","E.2D.O.M","E.2D.O.F",
                                                    "S.21D.Y.M","S.21D.Y.F","E.21D.Y.M","E.21D.Y.F",
                                                    "S.21D.O.M","S.21D.O.F","E.21D.O.M","E.21D.O.F")),
                   value = as.vector(value))
  df <- df[order(df$group),]
  rownames(df) <- 1:nrow(df)
  
  # Define colors for groups
  group_colors <- c("#f55d3f", "#c44532",  # Red shades
                    "#f5ae3f", "#c48632",  # Orange shades
                    "#EFF964", "#C4D14D",  # Yellow shades
                    "#59C751", "#3A8A39",  # Green shades
                    "#60ABF6", "#3A78B8",  # Blue shades
                    "#ffaee0", "#cc8ab2",  # Pink shades
                    "#9D54DE", "#6D3BA0",  # Purple shades
                    "#bbbbbb", "#888888")  # Gray shades
  
  # Generate boxplot
  ggplot(df, aes(x = group, y = value, fill = group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust = 0)) +
    ggtitle(gene) +
    scale_color_manual(values = group_colors) +
    scale_fill_manual(values = group_colors) +
    ylab("CPM") + xlab("Group") +
    theme(
      axis.text = element_text(size = 15, color = "black"),
      axis.title = element_text(size = 15, color = "black"),
      plot.title = element_text(size = 20, color = "black", face = "bold"),
      legend.title = element_text(size = 15, color = "black"),
      legend.text = element_text(size = 15, color = "black")
    )
}
