function(merged_seurat2, project_dir, 
                                    split.by = c("NSC", "IPC", "Imm.GExN", "M.GExN",
                                               "Unknown", "GABA IN", "Astrocyte", "OPC", 
                                               "Microglia", "Endothelial"),
                                    trend.in = "trimester",
                                    trend.order = c("first", "early", "mid"),
                                    pseudobulk.by = "donor",
                                    suppress.warnings = TRUE) {
  if(suppress.warnings) {
    options(warn = -1)
  }
  
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggrepel)
  
  # Create project directory
  dir.create(project_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Initialize dataframes
  all_plot_data <- data.frame(
    gene = character(),
    trimester = factor(levels = trend.order),
    expression = numeric(),
    cell_type = character(),
    type = character(),
    donor_value = numeric(),
    stringsAsFactors = FALSE
  )
  all_sig_genes <- data.frame()
  
  for (cell_type in split.by) {
    # Create pseudo-bulk object
    print(paste0("Making ",cell_type, " pseudobulk matrix by ", trend.in, " & ", pseudobulk.by))
    cell_seurat <- subset(merged_seurat2, sctype2 == cell_type)
    cell_seurat[["RNA"]] <- JoinLayers(cell_seurat[["RNA"]])
    pseudo_obj <- AggregateExpression(cell_seurat, 
                                    assays = "RNA", 
                                    return.seurat = TRUE, 
                                    group.by = c(trend.in, pseudobulk.by))
    pseudo_obj <- NormalizeData(pseudo_obj)
    
    # Calculate correlations
    expr_matrix <- GetAssayData(pseudo_obj, layer = "data")
    time_points <- 1:length(trend.order)
    
    # Calculate correlation and p-values
    print(paste0("Calculating Pearson's correlation coefficient of genes expression in ", trend.in," and p-values for ",cell_type))
    cor_results <- data.frame(
      gene = rownames(expr_matrix),
      cellType = cell_type,
      cor = sapply(1:nrow(expr_matrix), function(i) {
        expr_values <- sapply(trend.order, function(trim) {
          mean(expr_matrix[i, which(pseudo_obj[[trend.in]] == trim)])
          })
        tryCatch(cor.test(time_points, expr_values)$estimate,
            error = function(e) NA)
        }),
      pval = sapply(1:nrow(expr_matrix), function(i) {
        expr_values <- sapply(trend.order, function(trim) {
          mean(expr_matrix[i, which(pseudo_obj[[trend.in]] == trim)])
          })
        tryCatch(cor.test(time_points, expr_values)$p.value,
                 error = function(e) NA)
        }),
      slope = sapply(1:nrow(expr_matrix), function(i) {
        expr_values <- sapply(trend.order, function(trim) {
          mean(expr_matrix[i, which(pseudo_obj[[trend.in]] == trim)])
          })
        tryCatch(coef(lm(expr_values ~ time_points))[2],
                 error = function(e) NA)
        })
      )

    
    # Filter matrisome genes and get significant ones
    print(paste0("Filtering for significantly increasing & decreasing genes for ", cell_type))
    sig_genes <- cor_results[cor_results$pval < 0.05 & abs(cor_results$cor) > 0.5,]
    sig_genes <- sig_genes [complete.cases(sig_genes),]
    all_sig_genes <- rbind(all_sig_genes, sig_genes)
    
    # Get top genes
    top_increase <- head(sig_genes[sig_genes$slope > 0,][order(-sig_genes[sig_genes$slope > 0,]$slope),], 10)
    top_decrease <- head(sig_genes[sig_genes$slope < 0,][order(sig_genes[sig_genes$slope < 0,]$slope),], 10)
    
    # Prepare plot data
    print(paste0("Preparing plot data for ", cell_type))
    for (gene in c(top_increase$gene, top_decrease$gene)) {
      for (trim in trend.order) {
        expr_values <- expr_matrix[gene, which(pseudo_obj[[trend.in]] == trim)]
        all_plot_data <- rbind(all_plot_data, data.frame(
          gene = gene,
          trimester = factor(trim, levels = trend.order),
          expression = mean(expr_values),
          cell_type = cell_type,
          type = ifelse(gene %in% top_increase$gene, "Increasing", "Decreasing"),
          donor_value = expr_values
        ))
      }
    }
  }
  
  # Create plot with unique labels
  print(paste0("Generating a trend plot"))
  p <- ggplot(all_plot_data, aes(x = trimester, y = expression, group = gene)) +
    facet_wrap(~factor(cell_type, levels = split.by), 
              scales = "free_y", ncol = 5) +
    geom_point(aes(y = donor_value, color = type), alpha = 0.5, size = 0.8) +
    geom_line(aes(color = type), alpha = 0.8) +
    geom_text_repel(
      data = subset(all_plot_data, 
                   trimester == trend.order[1] & type == "Decreasing") %>%
        group_by(cell_type) %>%
        distinct(gene, .keep_all = TRUE) %>%
        slice_max(order_by = expression, n = 3),
      aes(label = gene),
      direction = "y",
      hjust = 1,
      nudge_x = -0.2,
      size = 2.5,
      box.padding = 0.5,
      force = 2
    ) +
    geom_text_repel(
      data = subset(all_plot_data, 
                   trimester == trend.order[length(trend.order)] & type == "Increasing") %>%
        group_by(cell_type) %>%
        distinct(gene, .keep_all = TRUE) %>%
        slice_max(order_by = expression, n = 3),
      aes(label = gene),
      direction = "y",
      hjust = 0,
      nudge_x = 0.2,
      size = 2.5,
      box.padding = 0.5,
      force = 2) +
    scale_color_manual(values = c("Increasing" = "tomato", "Decreasing" = "slateblue")) +
    theme_minimal() +
    theme(
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(size = 10, face = "bold"),
      axis.text = element_text(size = 8),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0.5, "lines")
    ) +
    labs(x = "Trimester", 
         y = "Expression Level (log)", 
         color = "Trend")
  
  # Save plot
  print(paste0("Saving results to ", project_dir))
  ggsave(file.path(project_dir, "temporal_trends_all_celltypes.pdf"), 
         p, width = 20, height = 9, dpi = 300)
  
  # Save significant genes table
  write.csv(all_sig_genes %>%
              select(gene, cellType, cor, pval) %>%
              arrange(cellType, desc(abs(cor))),
            file.path(project_dir, "all_significant_matrisome_genes.csv"),
            row.names = FALSE)
  
  # Calculate and save cell type rankings
  print(paste0("Ranking ", split.by, " by number of genes with high Pearson's correlation coefficent"))
  cell_rankings <- all_sig_genes %>%
    group_by(cellType) %>%
    summarize(
      total_sig_genes = n(),
      pos_cor = sum(cor > 0, na.rm = TRUE),
      neg_cor = sum(cor < 0, na.rm = TRUE),
      pos_neg_cor_ratio = pos_cor/neg_cor,
      strong_cor_genes = sum(abs(cor) > 0.5, na.rm = TRUE),
      percent_strong = round(strong_cor_genes/total_sig_genes * 100, 2),
      pos_coef = sum(slope > 0, na.rm = TRUE),
      neg_coef = sum(slope < 0, na.rm = TRUE),
      pos_neg_coef_ratio = pos_coef/neg_coef,
      strong_coef_genes = sum(abs(slope) > 0.5, na.rm = TRUE),
      percent_strong_coef = round(strong_coef_genes/total_sig_genes * 100, 2)
    ) %>%
    arrange(desc(strong_cor_genes))
  
  write.csv(cell_rankings,
            file.path(project_dir, "cell_type_rankings.csv"),
            row.names = FALSE)
  
  return(list(
    plot = p,
    data = all_plot_data,
    significant_genes = all_sig_genes,
    rankings = cell_rankings
  ))
}
