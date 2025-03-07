function(seurat_object, batch = "batch", cluster.anno = "cell_type", reduction = "umap", sample.size = 20000) {
  # Load required libraries
  library(Seurat)
  library(lisi)
  
  # Check if the object has more than sample.size cells
  if (ncol(seurat_object) > sample.size) {
    set.seed(42)  # For reproducibility
    cells_to_use <- sample(colnames(seurat_object), sample.size)
    seurat_object <- subset(seurat_object, cells = cells_to_use)
  }
  
  # Extract the integrated embedding
  integrated_embedding <- Embeddings(seurat_object, reduction = reduction)
  
  # Prepare metadata
  metadata <- data.frame(
    batch = seurat_object[[batch]][,1],
    cell_type = seurat_object[[cluster.anno]][,1]
  )
  
  # Calculate LISI scores
  lisi_results <- compute_lisi(integrated_embedding, metadata, c("batch", "cell_type"))
  
  # Rename columns for clarity
  colnames(lisi_results) <- c("iLISI", "cLISI")
  
  # Calculate mean scores
  mean_iLISI <- mean(lisi_results$iLISI)
  mean_cLISI <- mean(lisi_results$cLISI)
  
  # Return results
  return(list(
    lisi_scores = lisi_results,
    mean_iLISI = mean_iLISI,
    mean_cLISI = mean_cLISI
  ))
}
