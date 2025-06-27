suppressMessages({
  library(dplyr)
  library(ConsensusClusterPlus)
  library(tidyr)
})

###### Clustering of signatures activities ###### 
##' Input: 
##'   - activities_file.txt : signature activities from sigProfilerAssignment or sigProfilerExtractor
##'   - samples_info.txt : minimal information required is sample and patient
##'   - output_path
##'   
##' Run: Rscript iterative_clustering.R activities_file.txt patient_filter = T

# Retrieve command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript iterative_clustering.R <activities_file.txt> <samples_info.txt> <output_path>")
}

activities_file <- args[1]
samples_info_file <- args[2]
output_path <- args[3]

# Read input files
act <- read.table(activities_file, sep="\t", header=TRUE) 
colnames(act) <- gsub("\\.", "-", colnames(act))
sample_info <- read.table(samples_info_file, sep="\t", header=TRUE)

# Transform activities in proportions 
props <- act |>
  mutate(across(-Samples, ~ .x / rowSums(across(-Samples)))) %>% # Normalize row-wise
  inner_join(sample_info, by = c("Samples" = "sample")) %>%
  rename_with(~ gsub("CH44", "SCN-", .x)) %>%  # Rename columns
  filter(`SCN-G` <= 0.4) %>% # Filter out high signature G samples
  select(Samples, patient, everything()) %>%
  rename("sample" = "Samples")

# Create matrix
data_matrix <- as.matrix(props[,3:ncol(props)])
rownames(data_matrix) <- props$sample
input_mat <- t(data_matrix)

# Define different reproducible seeds for each run
set.seed(17171)
rand_seeds <- sample(10000:99999, 50)

# Run ConsensusClusterPlus
for(i in 1:50)
{
  r <- ConsensusClusterPlus(input_mat, maxK = 10, 
                            reps=1000, 
                            pItem=0.9, 
                            pFeature=1, 
                            title = paste0(output_path, "/euW2_", i), 
                            clusterAlg="hc", 
                            distance="euclidean", 
                            innerLinkage="ward.D2", 
                            finalLinkage="ward.D2", 
                            seed=rand_seeds[i], 
                            plot="png")
  saveRDS(r, paste0(output_path, "/hc_eu.", i, ".rds"))
  print(paste0("Done iteration ", i))
}
print("Finished")