# Code to reproduce Fig3 from Micoli et al. 2025
# Assumption: same features and number of clusters as in the original work

suppressMessages({
  library(ComplexHeatmap)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(circlize)
  library(dendextend)
})

# Dara required:
# - activities from sigProfilerExtractor/Assignment
# - metadata containing:
#     - HR: homologous recombination status (HRD/HRP)
#     - CDK12mut: mutational status of CDK12 (mut/none)
#     - CCNE1amp: copy number status of CCNE1 (amp/none)
#     - MYCamp: copy number status of MYC (amp/none)
#     - AKT2amp: copy number status of AKT2 (amp/none)
#     - BRCA1mut: mutational status of BRCA1 (mut/none)
#     - BRCA2mut: mutational status of BRCA1 (mut/none)
#     - complex structural variations counts (range = 0-5): chromoplexy, rigma, pyrgo, tyfonas, bfb, chromothripsis
#     - breaks: number of breakpoints in the genomes
#     - TMB: tumor mutational burden value
#     - LINE: LINE events (range = 0-100)
# - clustering object of the selected run (rds object)

# Load data
act <- read.table(file.path("~/your/path", "signature_activities.tsv"), sep="\t", header=T)
metadata <- read.table(file.path("~/your/path", "figure_metadata.tsv"), sep="\t", header=T)
cl_obj <- readRDS(file.path("~/your/path", "hc_eu.n.rds"))
output_path <- "~/your/output/path"

# Transform activities into proportions
props <- act %>%
  rowwise() %>%
  mutate(across(`SCN-A`:`SCN-K`, ~ .x / sum(c_across(`SCN-A`:`SCN-K`)))) %>%
  ungroup() %>%
  filter(sample %in% metadata$sample) # check to have a consistent sample set 

# Into matrix
mat <- t(as.matrix(props[,3:ncol(props)]))
colnames(mat) <- props$sample

# Retrieve cluster assignment
cl_assignment <- data.frame(sample = names(cl_obj[[6]]$consensusClass), 
                            cluster = cl_obj[[6]]$consensusClass)

# Join cluster assignment and proportion and order the data according to the consensus tree order
ordered_props <- props %>%
  inner_join(cl_assignment, by = c("sample" = "sample")) %>%
  mutate(
    cluster = factor(cluster,
                     levels = 1:6,
                     labels = c("Core HRP", "Proliferative HRP", "BRCA1-like HRD",
                                "EMT HRP", " ", "BRCA2-like HRD"))
  ) %>%
  slice(match(colnames(mat)[cl_obj[[6]]$consensusTree$order], sample)) %>%
  mutate(sample = factor(sample, levels = colnames(mat)[cl_obj[[6]]$consensusTree$order]))

# Create the matrix with the complex structural variation data
svs_mat <- as.matrix(t(metadata[, 13:18]))
colnames(svs_mat) <- metadata$sample
svs_mat <- ifelse(svs_mat > 5, 5, svs_mat) # Rescale extreme values

##### Plot #####
# Define color vectors
genes_colors <- c("TRUE"="red", "FALSE"="white")
amp_colors <- c("TRUE"="blue", "FALSE"="white")
signature_colors <- c("SCN-A" = "#ff477e", "SCN-B" = "#ff99ac", "SCN-E" = "#f7cad0", "SCN-I" = "#68d8d6", "SCN-J" = "#c4fff9",
                      "SCN-C" = "#9b72cf", "SCN-D" = "#ffc60a", "SCN-F" = "#60bf60", "SCN-H" = "#007fff", "SCN-G" = "#ab967d", "SCN-K" = "#b9b3af")
cluster_colors <- c("Core HRP"="#540d6e", "Proliferative HRP"="#ee4266", "EMT HRP" = "#a2d2ff", "BRCA1-like HRD"="#ffd23f", " " = "#b9b3af", "BRCA2-like HRD"="#0ead69")
HR_colors <- c("HRD"="#fb5607", "HRP"= "#1b4965")

## Dendrogram 
dhc <- as.dendrogram(cl_obj[[6]]$consensusTree)
dhc = color_branches(dhc, k=6, col = cluster_colors[c(5,3,1,2,6,4)])

## Bottom annotation
# Reshape metadata 
metadata$LINE <- ifelse(metadata$LINE >=100, 100, metadata$LINE)
metadata$logTMB <- log(metadata$tumorBurden)

# Annotation object
b <- HeatmapAnnotation(
  breaks = anno_barplot(metadata$breaks, 
                        bar_width = 1,
                        gp = gpar(fill = "#8632e6", col =NA),
                        axis_param = list(direction = "reverse",  at = c(0, 2000, 4000), labels = c("0", "2e3", "4e3"), side = "left")), 
  TMB = anno_lines(metadata$tumorBurden,
                   gp = gpar(col = 2), 
                   add_point = T, 
                   pt_gp = gpar(col = 5), 
                   pch = 1, 
                   size= unit(1, "mm"),
                   axis_param = list(at = c(0, 5e-06, 1e-05, 1.5e-05), labels = c("0", "5e-6", "1e-5", "1.5e-5"))),
  LINE = metadata$LINE,
  col = list(LINE = colorRamp2(c(0, max(metadata$LINE)), c("#fff1fa", "#ff0a54"))),
  annotation_height = unit(c(1,2,1), "cm"),
  gap = unit(c(1.5,1.5, 1.5), "mm"),
  border = c(breaks = T, TMB = T, LINE = T),
  which = "col",
  show_legend =F
)

## Top annotations
t <- HeatmapAnnotation(
  clusters = anno_block(gp = gpar(fill = cluster_colors)[c(5,3,1,2,6,4)],
                        labels = names(cluster_colors)[c(5,3,1,2,6,4)], 
                        labels_gp = gpar(col = "white", fontsize = 10)),
  sigs = anno_barplot(ordered_props[, 4:14], gp = gpar(fill = signature_colors, col = NA),
                      bar_width = 1, height = unit(10, "cm"), axis_param=list(at = NA, labels =NA, side = "left")),
  HR = metadata$HR,
  CDK12 = metadata$CDK12mut,
  CCNE1 = metadata$CCNE1amp,
  MYC = metadata$MYCamp,
  AKT2 = metadata$AKT2amp,
  BRCA1 = metadata$BRCA1mut,
  BRCA2 = metadata$BRCA2mut,
  
  col = list(HR = HR_colors,
             #wgdStatus = c("TRUE"="#fa9f42", "FALSE"="white"),
             CDK12 = genes_colors,
             CCNE1 = amp_colors,
             MYC = amp_colors,
             AKT2 = amp_colors,
             BRCA1 = genes_colors, 
             BRCA2 = genes_colors),
  na_col = "#ccc5b9",
  
  show_annotation_name = c(sigs = FALSE, clusters = FALSE), # only turn off `bar`
  border = c(HR = TRUE, CDK12 = T, CCNE1=T, MYC=T, AKT2=T, BRCA1=T, BRCA2=T, sigs=F),
  gap = unit(c(1.5, 1.5, 1.5, 0, 0, 0,0, 0), "mm"),
  annotation_height = unit(c(0.5, 8, rep(0.5,7)), "cm"),
  
  show_legend =F
)

## Heatmap
ht <- Heatmap(svs_mat, name = "count SVs", 
              col = c("#fbfaff","#C7ECD0", "#54CABE", "#0090BA", "#2F327D" ), 
              cluster_rows = FALSE, cluster_columns = dhc, border = T,
              column_dend_height = unit(2.5, "cm"), show_column_names = FALSE,
              column_title = NULL,
              top_annotation = t,
              bottom_annotation = b,
              column_split = 6,
              height = unit(3, "cm"),
              heatmap_legend_param = list(at = c(0, 2.5, 5), 
                                          labels = c("0", "2.5", ">5"), 
                                          labels_gp = gpar(size = 2),
                                          border = "black", 
                                          title_position = "topcenter")
)

## Legends
anno_legends <- list(
  Legend(at = c(0, 50, 100), labels = c("0", "50", ">100"), labels_gp = gpar(size = 2),col_fun = colorRamp2(c(0, 100), c("#fff1fa", "#ff0a54")), border = "black", 
         title_position = "topcenter", title = "LINE"),
  Legend(labels = names(signature_colors), title = "signatures",
         legend_gp = gpar(fill = signature_colors), border = "black", nrow = 4, by_row = F, title_position = "topleft"),
  Legend(labels = c("HRD", "HRP"), title="HRstatus", legend_gp = gpar(size=2, fill = HR_colors), border="black", title_position = "topcenter"),
  Legend(labels = c("mut", "amp", "none"), title= "genes", legend_gp = gpar(fill=c("red", "blue", "white")),
         border="black", title_position = "topcenter")
)

spacer <- Legend(labels = " ", legend_gp = gpar(fill = "white", col = "white"), title = NULL)#, width = unit(5, "mm"))
anno_legends_with_spacers <- list(
  spacer,
  anno_legends[[1]], spacer,
  anno_legends[[2]], spacer,
  anno_legends[[3]], spacer,
  anno_legends[[4]]
)
packed_legends_with_spacers <- ComplexHeatmap::packLegend(anno_legends_with_spacers, direction = "horizontal")

## Print
png(file.path(output_path, "clustering_heatmap.png"), height = 40, width = 35, units = "cm", res = 1000)
draw(ht, annotation_legend_list = packed_legends_with_spacers, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", 
     merge_legend = T)
dev.off()

