suppressMessages({
  library(dplyr)
  library(ggplot2)
})

# After visual inspection of statistics from consensusClusterPlus and 
# after having decided the number of clusters (in our case 6), extract the
# consensus attribution

cluster_path <- "path/to/folder/with/clustering_runs"
output_path <- "path/to/output/folder"

## Read all iterations and extract the grouping 
cl_runs <- NULL

for(i in 1:20) {
  cl_obj <- readRDS(paste0(cluster_path, "hc_eu.", i, ".rds"))
  cl_df <- data.frame(
    sample = names(cl_obj[[6]]$consensusClass),
    cluster = as.vector(cl_obj[[6]]$consensusClass),
    row.names = NULL
  )
  colnames(cl_df)[2] <- paste0("cluster", i)  # Rename cluster column
  
  # If cl_runs is NULL (first iteration), assign cl_df to cl_runs
  if (is.null(cl_runs)) {
    cl_runs <- cl_df
  } else {
    cl_runs <- inner_join(cl_runs, cl_df, by = "sample")
  }
  
  print(paste0("Joined ", i))
}
write.table(cl_runs, file.path(output_path, "clustering_runs_comparison.tsv"), sep="\t", col.names = T)

## Detect changing samples anmd the proportion of runs per each cluster in these samples
changing_samples <- data.frame()

for (i in 1:nrow(cl_runs)) {
  unique_clusters <- unique(as.numeric(cl_runs[i, 2:ncol(cl_runs)]))
  
  if (length(unique_clusters) > 1) {
    # Extract the first and second cluster
    first_cluster <- unique_clusters[1]
    second_cluster <- unique_clusters[2]
    
    # Compute proportions of each cluster
    cluster_counts <- table(as.numeric(cl_runs[i, 2:ncol(cl_runs)]))

    first_cluster_prop <- cluster_counts[as.character(first_cluster)] / (ncol(cl_runs)-1)
    second_cluster_prop <- cluster_counts[as.character(second_cluster)] / (ncol(cl_runs)-1)
    
    # Append to the dataframe
    changing_samples <- rbind(
      changing_samples,
      data.frame(
        sample = cl_runs[i, 1],
        first_cluster = first_cluster,
        second_cluster = second_cluster,
        first_cluster_prop = first_cluster_prop,
        second_cluster_prop = second_cluster_prop
      )
    )
  }
}
write.table(changing_samples, file.path(output_path, "changing_samples.tsv"), sep="\t", col.names = T)

## Create final clustering taking the highest proportion cluster for the changing samples
max_cluster_df <- changing_samples %>%
  mutate(
    cluster = ifelse(first_cluster_prop >= second_cluster_prop, first_cluster, second_cluster) # Choose the max proportion cluster
  ) %>%
  select(sample, cluster)

# Samples with unique clustering attribution
clusters <- cl_runs %>%
  filter(!sample %in% changing_samples$sample) %>%  # Exclude changing samples
  mutate(cluster = as.numeric(cluster1)) %>%  # Take the only clustering option
  select(sample, cluster) %>%
  rbind(max_cluster_df) %>%
  mutate(order = match(sample, cl_runs$sample)) %>%  # Create an ordering column
  arrange(order) %>%  # Sort by this column
  select(-order)
write.table(clusters, file.path(output_path, "clusters.tsv"), sep="\t", col.names = T)


