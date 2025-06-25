library(tidyverse)
library(pheatmap)

setwd("path/to/CIN-subtypes")
output_path <- "path/to/output_dir"

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

pdf(file.path(output_path, "heatmap_signatures.pdf"), width = 10, height = 4)
pheatmap(sigsM, cluster_rows = F, cluster_cols = F, fontsize = 7, display_numbers = formatC(sigsM, format = "f", digits = 0))
dev.off()

### 2. Intra-patient analysis
# Load activities
activities <- read.table(file.path("data/example", "De_Novo_Solution/Activities/De_Novo_Activities.txt"), sep="\t", header=T) |>
  rename(sample = Samples) |>
  mutate(patient = sub( '^(\\w+\\d+)_([piro])(.*)', '\\1', sample)) |>
  select(sample, patient, 2:12)
colnames(activities) <- gsub("\\.", "-", colnames(activities))

signature_colors <- c("SCN-A" = "#ff477e", "SCN-B" = "#ff99ac", "SCN-E" = "#f7cad0", "SCN-I" = "#68d8d6", "SCN-J" = "#c4fff9",
                      "SCN-C" = "#9b72cf", "SCN-D" = "#ffc60a", "SCN-F" = "#60bf60", "SCN-H" = "#007fff", "SCN-G" = "#ab967d", "SCN-K" = "#b9b3af")

for (i in unique(activities$patient))
{
  curr_df <- activities |>
    filter(patient == i) |>
    gather(signature, value, 3:13) |>
    mutate(sample = factor(sample, levels = unique(sample)),
           signature = factor(signature, levels = c("SCN-A", "SCN-B", "SCN-E", "SCN-I", "SCN-J", "SCN-C", "SCN-D", "SCN-F", "SCN-H", "SCN-G", "SCN-K")))
  
  sig_plot <- ggplot(curr_df, aes(x=sample, y=value, fill=signature)) +
    labs(title = paste0("Signatures in patient ", i)) +
    geom_bar(position="stack", stat="identity") +
    scale_fill_manual(values = signature_colors) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=6),
          plot.title = element_text(hjust=0.5, size = 10))
  ggsave(paste0(out_path, "/patients_plots/", i, ".png"), sig_plot, width = 20, height = 15, units = "cm", dpi="retina")
}

### 3. Variability of each signature in the cohort
activities |>
  gather(signature, value, 3:13) |>
  mutate(signature = factor(signature, levels = c("SCN-A", "SCN-B", "SCN-E", "SCN-I", "SCN-J", "SCN-C", "SCN-D", "SCN-F", "SCN-H", "SCN-G", "SCN-K"))) |>
  ggplot(aes(signature, value, fill = signature)) +
  geom_boxplot() +
  scale_fill_manual(values = signature_colors) +
  labs(title = "Signature variabily in the cohort", y = "SCN activity") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black"))

ggsave(file.path(out_path, "variability_signatures.png"), width = 20, height = 15, units = "cm", dpi="retina")
