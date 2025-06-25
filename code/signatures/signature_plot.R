# Code to reproduce Fig1a from Micoli et al. 2025
# Assumption: same features and number of signatures as in the original work

library(tidyverse)
library(gridExtra)
library(ggpubr)

setwd("path/to/CIN-subtypes")
output_path <- "path/to/output_dir"

# Load signatures definition
sigs_definition <- read.table(file.path("data", "SCN_Signatures.txt"), sep="\t", header=T)
colnames(sigs_definition) <- gsub("\\.", "-", colnames(sigs_definition))

# Color vectors
gpt_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#001c7f")
gpt_colors_faded <- c('#62a0cb', '#ffa556', '#6bbc6b', "#e26868", "#b495d1", "#eba0d4", "#a5a5a5", "#d0d164", "#5dd2dd", "#4d60a5") 

#Set the upper limit for the y scale
ceiling_dec <- function(x, level=1) round(x + 5*10^(-level-1), level)
up_limit <- ceiling_dec(max(apply(sigs_definition[,2:12], 2, range)), level = 2 )

#Add the class variable for grouping features
sigs_definition$Class <- rep(c(
  "Event type",
  "Deletions' length",
  "Duplications' length",
  "SV type",
  "CNV weight",
  "Segment size",
  "Segments per 5Mb",
  "Oscillation length",
  "Segments per ChrArm",
  "Changepoint"
), times = c(6, 6, 3, 6, 2, 3, 5, 4, 4, 5))
sigs_definition$Class <- factor(sigs_definition$Class, levels=unique(sigs_definition$Class))

#Gather the signature dataframe
melt_sigs <- gather(sigs_definition, key="Signature", value = "Value", 2:12)
melt_sigs$Class<-factor(melt_sigs$Class, levels=unique(melt_sigs$Class))

# Divide in single df per each signature
subsets <- lapply(unique(melt_sigs$Signature), function(sig) {
  melt_sigs %>% filter(Signature == sig)
})

# Create plot list
plots <- list()
for (i in 1:length(subsets)) {
  p <- ggplot(subsets[[i]], aes(x = interaction(MutationType, Class), y = Value, fill = Class, color=Class, group=Class)) +
    geom_bar(stat = "identity", position = position_dodge(), show.legend = F) +
    ylim(0, up_limit) +
    scale_fill_manual(values=gpt_colors_faded) +
    scale_color_manual(values = gpt_colors) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          axis.text.x = element_blank(), 
          #plot.title = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.line.x = element_line(color="black", linewidth = 0.3), 
          axis.line.y = element_line(color = "black", linewidth=0.3, ),
          axis.title = element_blank(), 
          strip.background = element_blank(),
          strip.text = element_blank()) +
    theme(axis.text.y = element_text(color = ifelse(i %% 4 == 1, "black", "transparent")),
          axis.ticks.y = element_line(color = ifelse(i %% 4 == 1, "black", "transparent"))) 
  
  positions <- data.frame(xpos = 6, ypos = 0.32)
  positions$text <- unique(subsets[[i]]$Signature)
  
  p <- p + geom_text(data=positions, aes(x = xpos, y = ypos,
                                         label = text), inherit.aes = F, size = 3)
  
  
  # Store each plot in the list
  plots[[i]] <- p
}
# Arrange the plots in a grid
final_p <- grid.arrange(grobs = plots, ncol = 4)

# Get legend through a fake plot
fake_plot <- ggplot(subsets[[1]], aes(x = interaction(MutationType, Class), y = Value, fill = Class, color=Class, group=Class)) +
  geom_bar(stat = "identity", position = position_dodge(), show.legend = T) +
  scale_fill_manual(values=gpt_colors_faded) +
  scale_color_manual(values = gpt_colors) +
  guides(fill = guide_legend(ncol = 2)) + # Set legend to two columns
  theme(
    legend.key.size = unit(0.5, "cm"), # Adjust legend item size
    legend.text = element_text(size =9) # Adjust legend text size if needed
  )

leg <- get_legend(fake_plot)

# Arrange plots and legend
ncol_grid <- 4
nrow_grid <- 3

grid_plots <- c(plots, rep(list(NULL), nrow_grid * ncol_grid - length(plots) - 1), list(leg))
sigs_plot <- do.call(grid.arrange, c(grid_plots, ncol = ncol_grid))
sigs_plot
ggsave(file.path(output_path, "signaturesSCN.pdf"), sigs_plot, width =14, height =7, dpi = "retina")
