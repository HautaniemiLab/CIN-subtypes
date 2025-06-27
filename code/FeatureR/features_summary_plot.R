# Code to reproduce Fig1a from Micoli et al. 2025
suppressMessages({
  library(StructuralVariantAnnotation)
  library(tidyverse)
  library(flexmix)
  library(gridExtra)
  library(ggpubr)
}) 

output_path <- "/results"

############ Features plot ##############
#### Barplot
# Load features and combine in a unique dataframe
features <- read.table("/results/features.tsv", sep="\t", header=T)

# Normalize features for better visualization
normalize <- function(x) {
  x[is.na(x)] <- 0
  return ((x - min(x)) / (max(x) - min(x)))
}

features_normalized <- apply(features[, 2:ncol(features)], 2, normalize ) |>
  as.data.frame() |>
  mutate(sample = features$sample) |>
  select(sample, 1:44)

# Calculate medians and quantiles per features
medians <- as.data.frame(t(apply(features_normalized[,2:ncol(features_normalized)], 2, median)))
quantiles <- as.data.frame(t(apply(features_normalized[,2:ncol(features_normalized)], 2, quantile)))

stats <- data.frame(FeatureType = names(medians),
                    Value = t(medians),
                    q25 = quantiles$`25%`,
                    q75 = quantiles$`75%`)

# Add the class column that groups features in categories
stats$Class <- rep(c(
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
stats$Class <- factor(stats$Class, levels=unique(stats$Class))

# Color vectors
gpt_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#001c7f")
gpt_colors_faded <- c('#62a0cb', '#ffa556', '#6bbc6b', "#e26868", "#b495d1", "#eba0d4", "#a5a5a5", "#d0d164", "#5dd2dd", "#4d60a5") 

# Position of class labels
positions <- data.frame(xmin = c(0.5, 6.5, 11.5, 15.5, 21.5, 23.5, 26.5, 31.5, 35.5, 39.5),
                        xmax = c(6.4, 11.4, 15.4, 21.4, 23.4, 26.4, 31.4, 35.4, 39.4, 44.5),
                        ymin = round(max(stats$q75), 2) + 0.02,
                        ymax= round(max(stats$q75), 2) + 0.12)
labels = c("Event type", "Deletions'\\nlength", "Duplications'\\nlength", "SV type", "CNV\\nweights", "Segment\\nsize", "Segments\\nper 5Mb", "Oscillation\\nlength", "Segments\\nper ChrArm", "Changepoint")
positions$text <- labels

# Build barplot
patient_plot <- ggplot(stats, aes(x=interaction(FeatureType, Class), y=Value, fill=Class, color=Class, group=Class)) +
  geom_col(position="dodge", show.legend = F)+
  geom_errorbar(aes(ymin = q25, ymax = q75), show.legend = F, width = 0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=9), legend.position = "right",
        legend.text = element_text(size=7), legend.title = element_text(size=8),
        legend.key.size = unit(0.5, "cm"), plot.margin = unit(c(0,0,2,0), 'lines'))+
  scale_x_discrete(labels=stats$FeatureType)+
  scale_fill_manual(values= gpt_colors_faded) +
  scale_color_manual(values = gpt_colors) +
  ylab("")+
  xlab("")+
  coord_cartesian(ylim=c(0,unique(positions$ymax)+0.001))

patient_plot <- patient_plot +
  geom_rect(data = positions, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),inherit.aes = F, color=gpt_colors, fill = gpt_colors_faded, alpha=0.7) +
  geom_text(data=positions, aes(x = xmin + (xmax - xmin) / 2, y = ymin + (ymax - ymin) / 2,
                                label = gsub("\\\\n", "\n", text)), inherit.aes = F, size = 2.5)

#patient_plot

#### Distributions
# Ricalculate non-discretized features 
# Load functions
source(file.path("FeatureR", "CNVfeatures_functions.R"))
source(file.path("FeatureR", "SVfeatures_functions.R"))
source(file.path("FeatureR", "utils.R"))

# Use the same table used to extract features
sample_data <- read.table("/data/example/sample_data.tsv", sep="\t", header=T)

# Quantify non-discretized features
segments <- read.table("/results/segmentation.tsv", sep="\t", header=T)

deldups_long <- deldup(sample_data)
ss <- segments %>%
  mutate(length = end - start) %>%
  filter(length>0) %>%
  filter(sample %in% sample_data$sample) %>%
  dplyr::select(sample, length)
bp5MB_counts <- bp5MB_extract(sample_data, segments)
oscil_counts <- oscil_extract(sample_data, segments)
bpArm_counts <- bpArm_extract(sample_data, segments)
changep_counts <- changp_extract(sample_data, segments)

non_discretized_features <- list(deldups_long$deletions, deldups_long$duplications, ss, bp5MB_counts, oscil_counts, bpArm_counts, changep_counts)

# Load models
models <- readRDS("/data/discretization_models.rds")
names(non_discretized_features) <- names(models)[c(6,7,1:5)]

# Get division lines for 5MB breaks
set.seed(17171)
dat <- bp5MB_counts$value
control<-new("FLXcontrol")
control@minprior<-0.001
control@iter.max<-1000

fit<-stepFlexmix(dat ~ 1,model = flexmix::FLXMCmvpois(),k=1:10,nrep=1,control=control)
fit<-getModel(fit,which="BIC")
posteriors <- flexmix::posterior(fit)

fiveMB_df <- data.frame(value = dat, group = as.character(fit@cluster)) %>% 
  group_by(group) %>% 
  summarise(mean= mean(value), min=min(value), max=max(value), .groups = "drop") %>% 
  arrange(min) %>% 
  mutate(line = (min - lag(max, default = 0))/2 + max)
lines <- log1p(fiveMB_df$line)[-5]
lines[1] <- lines[1]+0.4

# Colors and labels
color_vector <- c("#ff7f0e", "#2ca02c", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#001c7f")
filling_vector <- c('#ffa556', '#6bbc6b', "#eba0d4", "#a5a5a5", "#d0d164", "#5dd2dd", "#4d60a5")
x_label <- "log"

plot_list <- list()
for(f in 1:length(non_discretized_features))
{
  current_feature <- non_discretized_features[[f]]
  name_current_feature <- names(non_discretized_features)[f]
  
  current_model <- models[[name_current_feature]]
  current_model <- current_model[1:nrow(current_model)-1, ]
  
  if(grepl("log", colnames(current_model)[2]))
  {
    if(grepl("log1p", colnames(current_model)[2]))
    {
      current_feature$transf <- log1p(current_feature[,2])
    } else if (grepl("log$", colnames(current_model)[2]))
    {
      current_feature$transf <- log(current_feature[,2])
    }
    
    if(name_current_feature == "oscill")
    {
      p <- ggplot(current_feature, aes(x=transf)) +
        labs(#title = paste0(names_features[f]), 
          x = x_label) +
        geom_histogram(color=color_vector[f], fill=filling_vector[f], bins = 10) +
        theme_classic() +
        theme(plot.title = element_text(hjust=0.5, size = 10), axis.title.y = element_blank(), 
              axis.title.x = element_text(color = "#666666", size = 8), axis.text.x = element_text(size=8))
    } else {
      p <- ggplot(current_feature, aes(x=transf)) +
        labs(#title = paste0(names_features[f]), 
          x = x_label) +
        geom_density(color=color_vector[f], fill=filling_vector[f]) +
        theme_classic() +
        theme(plot.title = element_text(hjust=0.5, size = 10), axis.title.y = element_blank(), 
              axis.title.x = element_text(color = "#666666", size = 8), axis.text.x = element_text(size=8))
    }
    
    # Adding vertical lines
    for (x_val in current_model[, 2]) {
      p <- p + geom_vline(xintercept = x_val, linetype = "longdash", color = "red")
    }
    plot_list[[f]] <- p
  } else 
  {
    current_feature$transf <- log1p(current_feature[,2])
    
    p <- ggplot(current_feature, aes(x=transf)) +
      labs(#title = paste0(names_features[f]), 
        x = x_label) +
      geom_density(color=color_vector[f], fill=filling_vector[f]) +
      theme_classic() +
      theme(plot.title = element_text(hjust=0.5, size= 10), axis.title.y = element_blank(), 
            axis.title.x = element_text(color = "#666666", size = 8), axis.text.x = element_text(size=8))
    
    for (x_val in lines) {
      p <- p + geom_vline(xintercept = x_val, linetype = "longdash", color = "red")
    }
    
    plot_list[[f]] <- p
  }
}

density_plots <- grid.arrange(grobs = plot_list, ncol = 7)

features_overview <- ggarrange(patient_plot, density_plots, ncol = 1, heights = c(3,1)) 
ggsave(file.path(output_path, "features_overview.pdf"), features_overview, width =12, height = 8, dpi = "retina")
