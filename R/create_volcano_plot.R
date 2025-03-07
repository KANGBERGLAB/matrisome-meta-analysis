create_volcano_plot <- function(deseq_result, title = "Volcano Plot", fc_threshold = 1, p_threshold = 0.05, top.n = 10) {
  deseq_result <- within(deseq_result, {
    DEG <- "Not significant"
    DEG[avg_log2FC > fc_threshold & p_val_adj < p_threshold] <- "Upregulated"
    DEG[avg_log2FC < -fc_threshold & p_val_adj < p_threshold] <- "Downregulated"
    gene <- rownames(deseq_result)
  })
  top_up <- deseq_result %>%
    filter(DEG == "Upregulated") %>%
    top_n(top.n, avg_log2FC) %>%
    pull(gene)

  top_down <- deseq_result %>%
    filter(DEG == "Downregulated") %>%
    top_n(top.n, -avg_log2FC) %>%
    pull(gene)
  deseq_result$delabel <- ifelse(deseq_result$gene %in% c(top_up, top_down), deseq_result$gene, NA)
  mycolors <- c("slateblue1", "tomato2", "grey")
  names(mycolors) <- c("Downregulated", "Upregulated", "Not significant")
  p <- ggplot(data=deseq_result, aes(x=avg_log2FC, y=-log10(p_val_adj), col=DEG, label=delabel)) +
    geom_point(size =1.5, alpha = 0.5) + 
    theme_classic() +
    geom_text_repel(max.overlaps = Inf, size = 3) +
    scale_color_manual(values=mycolors) +
    geom_vline(xintercept=c(-fc_threshold, fc_threshold), col = "grey", linetype = 'dashed') +
    geom_hline(yintercept=-log10(p_threshold), col = "grey", linetype = 'dashed') +
    labs(x = expression("log"[2]*"FC"), y = "-log10(Adj. p-value)") +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5))

  return(p)
}
