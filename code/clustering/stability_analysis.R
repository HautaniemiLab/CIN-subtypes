suppressMessages({
  library(tidyverse)
  library(ggplot2)
})

setwd("path/CIN-subtypes")
output_path <- "path/to/output_folder"

## Analysis of patient stability: are all samples from the same patient?

# Load data: cluster data frame from cluster_analysis_and_selection.R and sample_data.tsv
sample_data <- read.table(file.path("path/to/sample_data.tsv"), sep="\t", header=T)
clusters <- read.table(file.path("path/to/clusters.tsv"), sep="\t", header=T)

# Prepare the data-frame for the analysis:
# - join folders
# - add a column with the sample timing information (p = pre-treatment, r = relapse, ...)

clusters <- clusters |>
  inner_join(sample_data) |>
  select(sample, patient, cluster) |>
  mutate(phase = sub( '^(\\w+\\d+)_([piro])(.*)', '\\2', sample), # add phase information
         cluster_pretreatment = NA, 
         cluster_raw = cluster) %>%
  rename()

# Determine for each patient the list of phases
stability_list <- list()

# Iterate over each unique patient
for(i in unique(clusters$patient)) {
  p <- clusters |> filter(patient == i) |> select(cluster, phase)
  unique_clusters <- unique(p$cluster)
  cl_unique <- paste(unique_clusters, collapse = ", ")
  
  # Case 1: General Stability (All samples fall in one cluster)
  if(length(unique_clusters) == 1) {
    stability_list[[i]] <- list(patient = i, stability_general = TRUE, stability_pretreatment = TRUE, all_clusters = cl_unique)
    clusters[clusters$patient == i, c("cluster", "cluster_pretreatment")] <- unique_clusters
    next
  } 
  
  # Case 2: General Stability (One differing sample)
  if(nrow(p) >= 6) {
    cl_counts <- p %>%
      group_by(cluster) %>%
      summarise(count = n(), .groups = "drop")
    
    max_cl <- cl_counts$cluster[which.max(cl_counts$count)]
    
    if(nrow(cl_counts) == 2 && any(cl_counts$count == 1)) {
      stability_list[[i]] <- list(patient = i, stability_general = TRUE, stability_pretreatment = TRUE, all_clusters = cl_unique)
      clusters[clusters$patient == i, c("cluster", "cluster_pretreatment")] <- max_cl
      next
    }
  }
  
  # Case 3: Pretreatment Stability (Only phase "p" samples fall in one cluster)
  pp <- p |> filter(phase == "p")
  unique_clusters_p <- unique(pp$cluster)
  
  if(length(unique_clusters_p) == 1) {
    stability_list[[i]] <- list(patient = i, stability_general = FALSE, stability_pretreatment = TRUE, all_clusters = cl_unique)
    clusters[clusters$patient == i & clusters$phase == "p", "cluster_pretreatment"] <- unique(pp$cluster)
    next
  } 
  
  # Case 4: Unstable
  stability_list[[i]] <- list(patient = i, stability_general = FALSE, stability_pretreatment = FALSE, all_clusters = cl_unique)
  clusters[clusters$patient == i & clusters$phase == "p", "cluster_pretreatment"] <- clusters[clusters$patient == i & clusters$phase == "p", "cluster_raw"]
}

# Convert the list to a data frame efficiently
stability <- bind_rows(stability_list)

# For convenience, join stability info in clusters table
clusters <- clusters |>
  inner_join(stability) |>
  select(sample, patient, phase, cluster, stability_general, cluster_pretreatment, stability_pretreatment, all_clusters, cluster_raw)

write.table(clusters, file.path(output_path, "clusters_with_stability.tsv"), sep="\t", col.names = T, row.names = F)
