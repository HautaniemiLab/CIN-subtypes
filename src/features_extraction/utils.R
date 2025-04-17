#' Annotation of segments from Purple segmentation files in a dataframe
#'
#' @param sample name
#' @param patient name
#' @param path to the Purple folder
#'
#' @return data frame with the segmentation of the specified sample added with
#'    sample and patient fields
#' @export
#'
read_and_annotate <- function(sample, patient, path) {
  file <- paste0(path, "/", sample, ".purple.cnv.somatic.tsv")
  if (file.exists(file)) {
    df <- read.table(file, sep="\t", header=T) %>%
      mutate(sample = sample, patient = patient)
    return(df)
  } else {
    return(NULL)
  }
}

#' Creation and manipulation of a segmentation file from Purple single files
#'
#' @param samples_data : data frame constituted by the fields "sample", "patient" and "path"
#'  sample: sample name
#'  patient: patient name
#'  path: path to the directory containing Purple results relative to the specified sample
#' @param purityCut : threshold for sample selection. All samples with purity < purityCut
#'  will be discarded from the analysis
#' @param verbose : default = F, if T speicifes which steps the function is performing
#'
#' @return segmentation of all the samples merged in a single data frame, filtered according to
#'  purityCut, corrected for purity bias and added with logR and Loh fields
#' @export
extract_segments <- function(samples_data, purityCut = -0.1, verbose = F)
{
  # retrieve purity and ploidy for each sample
  if(verbose)
  {
    print("Starting to add purity and ploidy information")
  }

  for (i in 1:nrow(samples_data)) {
    file <- paste0(samples_data[i, "path"], "/", samples_data[i, "sample"], ".purple.purity.tsv")
    if (file.exists(file)) {
      pur_pl <-read.table(file, sep="\t", header=T)
      samples_data[i, "ploidy"] <- pur_pl[1, "ploidy"]
      samples_data[i, "non_curated_purity"] <- pur_pl[1, "purity"]

      #Add checking of aberrancy using the purity-range method with a threshold of 0.5
      samples_data[i, "aberrant"] <- ifelse(pur_pl$maxPurity - pur_pl$minPurity <= 0.5, T, F)
    } else {
      print(paste0("The Purple purity file for the sample ", samples_data[i, "sample"], " doesn't exist."))
    }
  }

  # Join all segments retrieved from PURPLE results
  if(verbose) {  print("Starting to add segmentation") }

  all_segs <- map2_dfr(samples_data$sample, samples_data$patient, ~read_and_annotate(.x, .y, samples_data$path[which(samples_data$sample == .x & samples_data$patient == .y)]))

  #LogR, Loh, corrected purity and filtering according to the parameter purityCut
  if(verbose) { print("Addition of logR and Loh, purity correction and purity filtering") }
  suppressMessages({
    all_segs <- all_segs %>%
      inner_join(samples_data) %>%
      mutate(logR = ifelse(copyNumber<=0, -10, log(copyNumber/ploidy)),
             Loh = abs(baf - 0.5) * 2,
             purity = ifelse(aberrant == T, non_curated_purity, 0.00))}) #purity correction for non-aberrant samples
  if(length(unique(all_segs[all_segs$purity <= purityCut, "sample"]))>0)
  {
    print(paste0("The samples ", paste(unique(all_segs[all_segs$purity <= purityCut, "sample"]), collapse=", "), " have their purity <= the purity cut-off ", purityCut,
                 ", so they will be discarded from the analysis."))
  }

  all_segs <- all_segs %>%
    filter(purity >= purityCut) %>%
    dplyr::select(sample, patient, 1:16, "ploidy", "purity", 20:23)

  return(all_segs)
}

#

#' Discretization of features with Jenks natural breaks and plot
#'
#' @param df : data frame with the sample names and the corresponding feature values
#' @param feature_name : name of the feature. Possibilities: deletion, duplication, segsize, oscillation, arm_breakpoints, changepoint.
#' @param breaks : number of breaks to extract (the first and last breaks extracted are always
#'  the lowest and highest values of the distribution respectively)
#' @param seed : set a random number as seed to keep reproducibility
#' @param plot : default = F, possibility to plot the distribution of the feature with the
#'  discretization breaks
#'
#'  @export
extract_classes <- function(df, feature_name, breaks, seed = 1234, plot=F) {
  if(!is.na(seed))
  {
    set.seed(seed)
  }

  possible_features <- data.frame(name = c("deletion", "duplication", "segsize", "oscillation", "arm_breakpoints", "changepoint"),
                                  nickname = c("del", "dup", "segsize", "osc", "bpArm", "chg"),
                                  fill = c('#ffa556', '#6bbc6b', "#eba0d4", "#d0d164", "#5dd2dd", "#4d60a5"),
                                  color =  c("#ff7f0e", "#2ca02c", "#e377c2", "#bcbd22", "#17becf", "#001c7f"))
  if(!feature_name %in% possible_features$name) {
    print(paste0("Error: the possible feature names are: ", paste(possible_features, collapse = ", ")))
  }

  colnames(df) <- c("sample", "value")
  if (feature_name %in% c("deletion", "duplication", "segsize"))
  {
    df <- df %>% mutate(log = log(value)) %>% filter(!is.na(value))
  } else {
    df <- df %>% mutate(log = log1p(value)) %>% filter(!is.na(value))
  }

  jb <- getJenksBreaks(df$log, breaks)

  #Check for duplicates
  if(sum(duplicated(jb))>0)
  {
    print(paste0("There are duplicated break values that will be removed. Your duplicated values: ", jb[duplicated(jb)]))
    jb <- jb[!duplicated(jb)]
    breaks <- breaks -1
  }

  df$group <- cut(df$log, jb, labels = F, include.lowest = T)

  #Calculate the model of the fit (JB breaks in log)
  metrics_df <- data.frame(group = paste0(rep(possible_features[possible_features$name == feature_name, "nickname"], length(jb)-1), seq(1, length(jb)-1)),
                           break_log = jb[-1])
  #Discretize
  classes <- df %>% group_by(sample, group) %>% summarise(count = n(), .groups = "drop") %>% as.data.frame() %>% spread(group, count, fill = 0)
  colnames(classes) <- c("sample", metrics_df$group)

  results <- c(list(metrics_df), list(classes))

  #Plot distribution and discretization if asked
  if(plot)
  {
    p <- ggplot(df, aes(x=log)) +
      labs(title = paste0(feature_name, " ", breaks, " breaks"), x = "log") +
      #ggtitle(paste0(names_features[f])) +
      geom_density(color=possible_features[possible_features$name==feature_name, "color"], fill=possible_features[possible_features$name==feature_name, "fill"]) +
      theme_classic() +
      theme(plot.title = element_text(hjust=0.5, size = 10), axis.title.y = element_blank(),
            axis.title.x = element_text(color = "#666666", size = 8), axis.text.x = element_text(size=8))
    for (x_val in metrics_df[, 2][-(breaks-1)]) {
      p <- p + geom_vline(xintercept = x_val, linetype = "longdash", color = "red")
    }
    print(p)
  }

  return(results)
}

#' Discretization of features from eisting models
#'
#' @param df : data frame with the sample names and the corresponding feature values
#' @param feature_name : name of the feature. Possibilities: deletion, duplication, segsize, oscillation, arm_breakpoints, changepoint.
#' @param model : data frame with the parameters for the discretization deriving from the function `extract_classes`
#'
#' @export
classes_fromModel <- function(df, feature_name, model)
{
  possible_features <- data.frame(name = c("deletion", "duplication", "segsize", "oscillation", "arm_breakpoints", "changepoint"),
                                  nickname = c("del", "dup", "segsize", "osc", "bpArm", "chg"),
                                  transformation = c("log", "log", "log", "log1p", "log1p", "log1p"))
  #fill = c('#ffa556', '#6bbc6b', "#eba0d4", "#d0d164", "#5dd2dd", "#4d60a5"),
  #color =  c("#ff7f0e", "#2ca02c", "#e377c2", "#bcbd22", "#17becf", "#001c7f"))
  if(!feature_name %in% possible_features$name) {
    print(paste0("Error: the possible feature names are: ", paste(possible_features$name, collapse = ", ")))
  }

  colnames(df) <- c("sample", "value")

  if (grepl("1p", possible_features[possible_features$name == feature_name, "transformation"]))
  {
    df$transf <- log1p(df$value)
  } else {
    df$transf <- log(df$value)
  }
  df <- df[!is.na(df$value),]

  #Extract breaks from models
  breaks <- model[, 2]
  breaks <- append(0, breaks)

  #Cutting and retrieving the classes
  df$group <- cut(df$transf, breaks, labels = F, include.lowest = T)
  df[df$transf>breaks[length(breaks)], "group"] <- nrow(model) #if higher than the maximum, top group
  classes <- df %>%
    group_by(sample, group) %>%
    summarise(count = n(), .groups = "drop") %>%
    as.data.frame() %>% #count the elements per group
    spread(group, count, fill = 0) #spread making different columns

  colnames(classes)[2:ncol(classes)] <- model[,1]

  return(classes)
}

MM_fromModel <- function(df, feature_name, model)
{
  possible_features <- c("compact_sparse", "bp5MB")

  if(!feature_name %in% possible_features) {
    print(paste0("Error: the possible feature names are: ", paste(possible_features, collapse = ", ")))
  }

  if(feature_name == "compact_sparse")
  {
    df <- df[, c("sample", "dist")]
  }

  colnames(df) <- c("sample", "value")
  dat <- df$value

  if(feature_name == "compact_sparse")
  {
    postDatUnscaled = sapply(1:nrow(model), function(x) dnorm(x=dat, model[[x ,"mean"]], model[[x ,"sd"]]))
  } else {
    postDatUnscaled = sapply(1:nrow(model), function(x) dpois(x=dat, lambda = model[[x ,"mean"]]))
    
    # postDatUnscaled = sapply(1:nrow(model), function(x) {
    #   model$prior[x] * dpois(x = dat, lambda = model$mean[x])
    #   })
  }

  postDatScaled = data.frame( postDatUnscaled / rowSums(postDatUnscaled) )
  postDatScaled$sample = df$sample
  data_discr = aggregate(postDatScaled, . ~ sample, sum)
  data_discr[, 2:ncol(data_discr)] <- round(data_discr[, 2:ncol(data_discr)])
  colnames(data_discr)[2:ncol(data_discr)] <- model$group

  return(data_discr)
}


