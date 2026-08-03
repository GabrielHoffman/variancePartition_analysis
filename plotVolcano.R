
library(data.table)
library(ggrepel)

.plotVolcano = function(tab, nGenes = 5, size = 12, minp = 1.0e-310, cutoff = 0.05, ncol = 3, ...){

  if ("ID" %in% colnames(tab)) {
    df_combine <- data.table(tab)
  } else {
    df_combine <- data.table(ID = rownames(tab), tab)
  }

  xmax <- max(abs(df_combine$logFC))
  ymax <- -log10(min(df_combine$P.Value))

  # Pass R CMD check
  .SD <- logFC <- P.Value <- isSignif <- ID <- NULL

  # check for 0 p-values
  if (!is.finite(ymax)) {
    nzero <- length(df_combine$P.Value == 0)
    txt <- paste0("There are ", nzero, " features with p-value of 0. Plotting will be affected")
    warning(txt)
  }

  df_combine$isSignif <- c("no", "yes")[(df_combine$FDR < cutoff) + 1]
  df_combine$P.Value <- pmax(minp, df_combine$P.Value)

  # top significant genes in each cell type
  df2 <- df_combine[, head(.SD, nGenes)]

  # reverse order to plot significant points last
  ggplot(df_combine[seq(nrow(df_combine), 1), , drop = FALSE], aes(logFC, -log10(P.Value), color = isSignif)) +
    geom_point() +
    theme_classic(size) +
    theme(aspect.ratio = 1, legend.position = "none", plot.title = element_text(hjust = 0.5)) +
    xlab(bquote(log[2] ~ fold ~ change)) +
    ylab(bquote(-log[10] ~ P)) +
    scale_color_manual(values = c("grey", "darkred")) +
    scale_y_continuous(limits = c(0, ymax * 1.02), expand = c(0, 0)) +
    geom_text_repel(data = df2, aes(logFC, -log10(P.Value), label = ID), segment.size = .5, segment.color = "black", color = "black", force = 1, nudge_x = .005, nudge_y = .5)
}
