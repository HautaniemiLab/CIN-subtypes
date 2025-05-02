library(tidyverse)
library(pheatmap)

setwd("path/to/CIN-subtypes")

######## Analysis of extracted or quantified signatures ######## 
### 1. Signatures definition
# Load signatures
sigs <- read.table(file.path("data", "SCN_Signatures.txt"), sep="\t", header=T)

# Convert in matrix
sigsM <- as.matrix(t(sigs[, 2:ncol(sigs)]))
colnames(sigsM) <- sigs$MutationType

# Multiply by 1000 to make the weigths readable
sigsM <- round(sigsM * 1000, digits = 0)
mode(sigsM) <- "numeric"

pdf(file.path("output/path", "heatmap_signatures.pdf"), width = 10, height = 4)
pheatmap(sigsM, cluster_rows = F, cluster_cols = F, fontsize = 7, display_numbers = formatC(sigsM, format = "f", digits = 0))
dev.off()

### 2. Intra-patient analysis
# Load activities
... ongoing
activities <- read.table(file.path(act_path, "activities_set1-15.tsv"), sep="\t", header=T)
sample_info <- read.table(file.path(home, "SCNA_Purple/resources/sample_info_extended_231004.csv"), sep="\t", header=T)
activities_new <- activities |>
  inner_join(sample_info[, c("sample","set")]) |>
  mutate(set_type = ifelse(grepl("14|15", set), "new", "old")) |>
  select(-set) |>
  group_by(patient) |>
  filter(any(set_type == "new")) |>
  ungroup() |>
  as.data.frame()

complete_colors <- c('#4ad66d','#9be564', '#ffe119', '#9a6324',  '#f58231', '#d80032','#FFC0CB','#f032e6', '#911eb4', '#4363d8', '#48cae4', '#ade8f4')
for (i in unique(activities_new$patient))
{
  curr_df <- activities_new |>
    filter(patient == i) |>
    gather(signature, value, 3:13) |>
    mutate(sample = factor(sample, levels = unique(sample)))
  
  sig_plot <- ggplot(curr_df, aes(x=sample, y=value, fill=signature, color = set_type)) +
    labs(title = paste0("Signatures in patient ", i)) +
    geom_bar(position="stack", stat="identity") +
    scale_fill_manual(values = complete_colors) +
    scale_color_manual(values = c("new" = "black", old = "white")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=6),
          plot.title = element_text(hjust=0.5, size = 10))
  ggsave(paste0(out_path, "/patients_plots/", i, ".png"), sig_plot, width = 20, height = 15, units = "cm", dpi="retina")
}
