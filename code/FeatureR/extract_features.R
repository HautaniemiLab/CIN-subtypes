library(tidyverse)
library(StructuralVariantAnnotation)
library(ShatterSeek)

# Set working directory to this repo
setwd("path/to/CIN-subtypes")

# Load functions
source(file.path("code", "FeatureR", "CNVfeatures_functions.R"))
source(file.path("code", "FeatureR", "SVfeatures_functions.R"))
source(file.path("code", "FeatureR", "utils.R"))

##### Extraction of features based on pre-existing models #####
## Input data 
#' Data frame constituted by the fields "sample", "patient" and "path"
#'  - sample: sample name
#'  - patient: patient name
#'  - path: path to the directory containing PURPLE and LINX results relative to the specified sample
sample_data <- read.table("path/to/sample_data.tsv", sep="\t", header=T)

## Load models
models <- readRDS(file.path("data", "discretization_models.rds"))

## Define output path 
output_path <- "output/path"

#### Create segmentation file ####
# Optional: modify the parameter purityCut in extract_segments, default = 0.2
segments <- extract_segments(sample_data, verbose =T)
write.table(segments, file.path(output_path, "segmentation.tsv"), sep="\t", col.names = T, row.names = F)

# If samples have been filtered out by purityCut, remove them from sample_data
sample_data <- sample_data %>% filter(sample %in% unique(segments$sample))

#### Extract SV features #### 
## Complex and single SVs
cosim_feat <- cosim(sample_data)

## Intra-chromosomal and inter-chromosomal SVs
trater_feat <- trater(sample_data)

## Complex and sparse SVs
breakpoints <- bk_distance(sample_data)
compact_sparse_feat <- MM_fromModel(breakpoints, "compact_sparse",models[["compact_sparse"]])

#Deletion and duplication lenghts
deldups_long <- deldup(sample_data)
dels_feat <- classes_fromModel(deldups_long$deletions, "deletion", models[["deletions"]])
dups_feat <- classes_fromModel(deldups_long$duplications, "duplication", models[["duplications"]])

# SVs types: insertions, unbalanced translocation, reciprocal translocation, reciprocal inversions, LINE insertions
sv_types <-SV_categories(sample_data)

# Magnitudes: quantiles 0.25 and 0.75 of the logR
magnitudes_feat <- logR_quantiles(sample_data, segments)

# Chromotripsis using Shatterseek
chrmtrps_feat <- is.chromotripsis(sample_data, segments)

#### Extract CNV features ####
ss <- segments %>%
  mutate(length = end - start) %>%
  filter(length>0) %>%
  filter(sample %in% sample_data$sample) %>%
  dplyr::select(sample, length)
ss_feat <- classes_fromModel(ss, "segsize", models[["segsize"]])

bp5MB_counts <- bp5MB_extract(sample_data, segments)
bp5MB_feat <- MM_fromModel(bp5MB_counts, "bp5MB", models[["bp5MB"]])

oscil_counts<- oscil_extract(sample_data, segments)
oscil_feat <- classes_fromModel(oscil_counts, "oscillation", models[["oscill"]])

bpArm_counts <- bpArm_extract(sample_data, segments)
bpArm_feat <- classes_fromModel(bpArm_counts, "arm_breakpoints", models[["bpArm"]])

changep_counts <- changp_extract(sample_data, segments)
changep_feat <- classes_fromModel(changep_counts, "changepoint", models[["changepoint"]])

#### Aggregation in one feature dataframe ####
features <- Reduce(function(df1, df2) merge(df1, df2, by = "sample", all = TRUE),
                   list(cosim_feat, trater_feat, compact_sparse_feat, dels_feat, dups_feat, sv_types,
                        magnitudes_feat, chrmtrps_feat, ss_feat, bp5MB_feat, oscil_feat,
                        bpArm_feat, changep_feat))
write.table(features, file.path(output_path, "features.tsv"), sep='\t', col.names = TRUE)

#### Once created the feature data frame, it can be used for signature extraction (de novo signatures) or for signature 
#### attribution (calculating activities from existing signatures). In any case, the extracted features the feature data
#### frame must be reshaped. After this step proceed with:
#### - de novo signature extraction: run_sigProfilerExtractor.py
#### - signature attribution: run_sigProfilerAssignment.py
ex_mat <- t(as.matrix(features[, -1]))
colnames(ex_mat) <- features$sample
ex_mat[is.na(ex_mat)] <- 0

write.table(ex_mat, file.path(output_path, "extraction_df.tsv"), sep="\t", col.names = T, row.names =F)
