calculate_integration_separation_scores <- function(object, batch, cluster.anno, seurat_clusters, sample_size = 10000, sampling.size = 3) {
  # Original data extraction with seurat clusters
  umap_data <- data.frame(
    UMAP1 = Embeddings(object, "umap")[,1],
    UMAP2 = Embeddings(object, "umap")[,2],
    batch = object@meta.data[[batch]],
    cluster = object@meta.data[[cluster.anno]],
    seurat_cluster = object@meta.data[[seurat_clusters]]
  )
  
  # Sampling step for large datasets
  if (nrow(umap_data) > sample_size) {
    set.seed(42)
    sampled_indices <- sample(1:nrow(umap_data), sample_size)
    umap_data <- umap_data[sampled_indices, ]
  }
  
  # Create random sampling for baseline
  set.seed(42)
  umap_data$random_batch <- sample(1:sampling.size, nrow(umap_data), replace = TRUE)
  
  # Calculate pairwise distances
  dist_matrix <- as.matrix(dist(umap_data[, c("UMAP1", "UMAP2")]))
  
  # Modified calculate_scores function with local structure penalty
  calculate_scores <- function(group_var) {
    unique_groups <- unique(umap_data[[group_var]])
    scores <- list()
    
    # Get total number of unique Seurat clusters
    total_clusters <- length(unique(umap_data$seurat_cluster))
    
    for (i in seq_along(unique_groups)) {
      current_group <- unique_groups[i]
      other_groups <- unique_groups[-i]
      
      current_indices <- which(umap_data[[group_var]] == current_group)
      other_indices <- which(umap_data[[group_var]] %in% other_groups)
      
      # Calculate distinct clusters between current batch and others
      current_clusters <- unique(umap_data$seurat_cluster[current_indices])
      other_clusters <- unique(umap_data$seurat_cluster[other_indices])
      distinct_clusters <- length(setdiff(current_clusters, other_clusters))
      
      # Calculate local structure penalty
      cluster_penalty <- distinct_clusters / total_clusters
      
      internal_dist <- mean(dist_matrix[current_indices, current_indices])
      external_dist <- mean(dist_matrix[current_indices, other_indices])
      
      # Calculate base score and apply penalty
      base_score <- abs(internal_dist - external_dist) / (internal_dist + external_dist)
      final_score <- base_score * (1 + cluster_penalty)  # Increase score proportionally to distinct clusters
      
      scores[[i]] <- data.frame(
        group = current_group,
        internal_dist = internal_dist,
        external_dist = external_dist,
        distinct_clusters = distinct_clusters,
        cluster_penalty = cluster_penalty,
        base_score = base_score,
        score = final_score
      )
    }
    
    do.call(rbind, scores)
  }
  
  # Calculate scores for actual data and random sampling
  integration_scores <- calculate_scores("batch")
  separation_scores <- calculate_scores("cluster")
  random_integration_scores <- calculate_scores("random_batch")
  
  # Calculate overall scores
  overall_integration_score <- mean(integration_scores$score)
  overall_separation_score <- mean(separation_scores$score)
  best_integration_score <- mean(random_integration_scores$score)
  worst_separation_score <- mean(random_integration_scores$score)
  
  # Normalize scores based on random baseline
  normalized_integration_score <- (overall_integration_score - best_integration_score) / (1 - best_integration_score)
  normalized_separation_score <- (overall_separation_score - worst_separation_score) / (1 - worst_separation_score)
  
  return(list(
    integration_scores = integration_scores,
    separation_scores = separation_scores,
    overall_integration_score = overall_integration_score,
    overall_separation_score = overall_separation_score,
    best_integration_score = best_integration_score,
    worst_separation_score = worst_separation_score,
    normalized_integration_score = normalized_integration_score,
    normalized_separation_score = normalized_separation_score
  ))
}
