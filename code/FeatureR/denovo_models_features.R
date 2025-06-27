suppressMessages({
  library(tidyverse)
  library(StructuralVariantAnnotation)
  library(ShatterSeek)
  library(mclust)
  library(BAMMtools)
  library(flexmix)
}) 

# Load functions
source(file.path("FeatureR", "CNVfeatures_functions.R"))
source(file.path("FeatureR", "SVfeatures_functions.R"))
source(file.path("FeatureR", "utils.R"))

##### De novo extraction of features and models #####
## Input data 
#' Data frame constituted by the fields "sample", "patient" and "path"
#'  - sample: sample name
#'  - patient: patient name
#'  - path: path to the directory containing PURPLE and LINX results relative to the specified sample
sample_data <- read.table("/data/example/sample_data.tsv", sep="\t", header=T) |>
  mutate(path = paste0("/", path))

## Define output path 
output_path <- "/results"

#### Create segmentation file ####
# Optional: modify the parameter purityCut in extract_segments, default = 0.2
segments <- extract_segments(sample_data, verbose =T)
write.table(segments, file.path(output_path, "segmentation.tsv"), sep="\t", col.names = T, row.names = F)
print("Segmentation written")

# If samples have been filtered out by purityCut, remove them from sample_data
sample_data <- sample_data %>% filter(sample %in% unique(segments$sample))

#### De novo SV features #### 
print("############## Starting SV features de novo extraction ##############")

seed = 17171

## Complex and single SVs
cosim_feat <- cosim(sample_data)

## Intra-chromosomal and inter-chromosomal SVs
trater_feat <- trater(sample_data)

## Complex and sparse SVs
breakpoints <- bk_distance(sample_data)
compact_sparse_feat <- extract_compact_sparse(breakpoints, seed) 

#Deletion and duplication lenghts
deldups_long <- deldup(sample_data)
dels_feat <- extract_classes(deldups_long$deletions, "deletion", 6, seed, plot = F)
dups_feat <- extract_classes(deldups_long$duplications, "duplication", 5, seed, plot = F)

# SVs types: insertions, unbalanced translocation, reciprocal translocation, reciprocal inversions, LINE insertions
sv_types <-SV_categories(sample_data)

# Magnitudes: quantiles 0.25 and 0.75 of the logR
magnitudes_feat <- logR_quantiles(sample_data, segments)

# Chromotripsis using Shatterseek
chrmtrps_feat <- is.chromotripsis(sample_data, segments)

#### De novo CNV features ####
print("############## Starting CNV features de novo extraction ##############")

ss <- segments %>%
  mutate(length = end - start) %>%
  filter(length>0) %>%
  filter(sample %in% sample_data$sample) %>%
  dplyr::select(sample, length)
ss_feat <- extract_classes(ss, "segsize", 4, seed, plot = F)

bp5MB_counts <- bp5MB_extract(sample_data, segments)
bp5MB_feat <- bp5MB_discretize(bp5MB_counts, seed, plot = F)

oscil_counts<- oscil_extract(sample_data, segments)
oscil_feat <- extract_classes(oscil_counts, "oscillation", 5, seed, plot = F)

bpArm_counts <- bpArm_extract(sample_data, segments)
bpArm_feat <- extract_classes(bpArm_counts, "arm_breakpoints", 5, seed, plot = F)

changep_counts <- changp_extract(sample_data, segments)
changep_feat <- extract_classes(changep_counts, "changepoint", 6, seed, plot = F)

#### Aggregation in one feature dataframe and one model list ####
features <- Reduce(function(df1, df2) merge(df1, df2, by = "sample", all = TRUE),
                   list(cosim_feat, trater_feat, compact_sparse_feat[[2]], dels_feat[[2]], dups_feat[[2]], sv_types,
                        magnitudes_feat, chrmtrps_feat, ss_feat[[2]], bp5MB_feat[[2]], oscil_feat[[2]],
                        bpArm_feat[[2]], changep_feat[[2]]))
write.table(features, file.path(output_path, "denovo_features.tsv"), sep='\t', col.names = TRUE)

models <- list(compact_sparse = compact_sparse_feat[[1]], deletions = dels_feat[[1]], duplications = dups_feat[[1]],
               segsize = ss_feat[[1]], bp5MB = bp5MB_feat[[1]], oscill = oscil_feat[[1]], bpArm = bpArm_feat[[1]],
               changepoint = changep_feat[[1]])
saveRDS(models, file.path(output_path, "denovo_models.rds"))

#### Once created the feature data frame, it can be used for signature extraction (de novo signatures) or for signature 
#### attribution (calculating activities from existing signatures). In any case, the extracted features the feature data
#### frame must be reshaped. After this step proceed with:
#### - de novo signature extraction: run_sigProfilerExtractor.py
#### - signature attribution: run_sigProfilerAssignment.py
ex_mat <- t(as.matrix(features[, -1]))
colnames(ex_mat) <- features$sample
ex_mat[is.na(ex_mat)] <- 0

write.table(ex_mat, file.path(output_path, "denovo_extraction_df.tsv"), sep="\t", col.names = T, row.names =T)

print("denovo_features.tsv, denovo_models.rds and denovo_extraction_df.tsv have been correctly saved in the results folder")
