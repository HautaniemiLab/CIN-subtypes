suppressMessages({
  library(tidyverse)
  library(ggplot2)
})

# After visual inspection of statistics from consensusClusterPlus and 
# after having decided the number of clusters (in our case 6), extract the
# consensus attribution

cluster_path <- "/results/multiple_clustering" # path to clustering runs
output_path <- "/results"

## Read all iterations and extract the grouping 
n <- 6 #set the number of clusters after inspection of consensusClusterPlus output
cl_runs <- NULL

for(i in 1:50) {
  cl_obj <- readRDS(paste0(cluster_path, "/hc_eu.", i, ".rds"))
  cl_df <- data.frame(
    sample = names(cl_obj[[n]]$consensusClass),
    cluster = as.vector(cl_obj[[n]]$consensusClass),
    row.names = NULL
  )
  colnames(cl_df)[2] <- paste0("run", i)  # Rename clustering column
  
  # If cl_runs is NULL (first iteration), assign cl_df to cl_runs
  if (is.null(cl_runs)) {
    cl_runs <- cl_df
  } else {
    cl_runs <- inner_join(cl_runs, cl_df, by = "sample")
  }
  
  print(paste0("Joined ", i))
}
write.table(cl_runs, file.path(output_path, "clustering_runs_comparison.tsv"), sep="\t", col.names = T)

## Statistics about the iterative clustering
cl_stats <- cl_runs |>
  pivot_longer(cols = starts_with("run"), names_to = "cluster_id", values_to = "cluster_value") |>
  group_by(sample, cluster_value) |>
  summarise(count = n(), .groups = "drop") |>
  group_by(sample) |>
  mutate(max_count = max(count)) |>
  filter(count == max_count) |>
  arrange(sample, cluster_value) |>
  reframe(
    majority_cluster = first(cluster_value),
    tied_with = if (n() > 1) paste(cluster_value[-1], collapse = ",") else NA_character_,
    majority_pct = round(first(count) / 20, 2),
    disagreement = majority_pct < 1,
    tie = n() > 1
  )
write.table(cl_stats, file.path(output_path, "clustering_statistics.tsv"), sep="\t", col.names = T)

## Selection of the run with the best clustering
# NB: ties are handled by taking the first cluster assignment and not counted when selecting the run with the best match. 
# It is possible to change or exclude these samples manually.
best_match <- cl_runs |>
  left_join(
    cl_stats |> filter(majority_pct != 0.5) |> select(sample, majority_cluster),
    by = "sample"
  ) |>
  pivot_longer(cols = starts_with("run"), names_to = "run", values_to = "cluster_value") |>
  mutate(match = cluster_value == majority_cluster) |>
  group_by(run) |>
  summarise(
    agreement_rate = mean(match, na.rm = TRUE),
    n_matches = sum(match, na.rm = TRUE)
  ) |>
  arrange(desc(agreement_rate)) |>
  slice(1)

## Extraction of the final clustering
cl_final <- cl_runs |>
  select("sample", best_match$run) |>
  setNames(c("sample", "cluster"))

write.table(cl_final, file.path(output_path, "clusters.tsv"), sep="\t", col.names = T, row.names = F)