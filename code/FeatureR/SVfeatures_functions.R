#' Complex and simple events
#' @description
#' Calculates the number of complex events and simple events per sample
#'
#' @param samples_data : data frame constituted by the fields "sample", "patient" and "path"
#'  sample: sample name
#'  patient: patient name
#'  path: path to the directory containing Purple results relative to the specified sample
#'
#' @export
#'
cosim <- function(samples_data)
{
  cosim_df <- data.frame()
  for (i in 1:nrow(samples_data))
  {
    clusters <- read.table(paste0(samples_data[i, "path"], "/", samples_data[i, "sample"], ".linx.clusters.tsv"), header=T, sep="\t")
    simple_ev <- clusters %>% filter(category=="SIMPLE") %>% nrow()
    complex_ev <- clusters %>% filter(category!="SIMPLE") %>% nrow()

    cosim_df <- rbind(cosim_df, c(samples_data[i, "sample"], simple_ev, complex_ev))
  }
  colnames(cosim_df) <- c("sample", "simple", "complex")

  return(cosim_df)
}

#' Interchromosomal and intrachromosomal events
#' @description
#' Calculates the number of intra/inter-chromosomal events per sample according to the pairs of SVs laying on the same chromosome or not
#'
#'
#' @param samples_data :  data frame constituted by the fields "sample", "patient" and "path"
#'  sample: sample name
#'  patient: patient name
#'  path: path to the directory containing Purple results relative to the specified sample
#'
#' @export
#'
trater <- function(samples_data)
{
  trater_df <- data.frame()
  for (i in 1:nrow(samples_data))
  {
    bks <- read.table(paste0(samples_data[i, "path"], "/", samples_data[i, "sample"], ".linx.svs.tsv"), header=T, sep="\t")
    svs <- readVcf(paste0(samples_data[i, "path"], "/", samples_data[i, "sample"], ".purple.sv.vcf.gz"))
    suppressWarnings({gr <- breakpointRanges(svs, inferMissingBreakends=F)})
    couples_bks <- breakpointgr2bedpe(gr)

    interchr <- couples_bks %>% filter(chrom1 != chrom2) %>% nrow()
    intrachr <- couples_bks %>% filter(chrom1 == chrom2) %>% nrow()
    trater_df <- rbind(trater_df, c(samples_data[i, "sample"], interchr, intrachr))
  }
  colnames(trater_df) <- c("sample", "interchr", "intrachr")

  return(trater_df)
}

#' Breakpoints distance
#' @description
#' In complex events, it calculated the distance between the breakpoints.
#'
#' @param samples_data :  data frame constituted by the fields "sample", "patient" and "path"
#'  sample: sample name
#'  patient: patient name
#'  path: path to the directory containing Purple results relative to the specified sample
#'
#' @export
#'
bk_distance <- function(samples_data)
{
  breakpoints <- data.frame()
  for(i in 1:nrow(samples_data))
  {
    clusters <- read.table(paste0(samples_data[i, "path"], "/", samples_data[i, "sample"], ".linx.clusters.tsv"), header=T, sep="\t")
    clusters <- clusters %>% filter(category=="COMPLEX")
    cl_ids <- clusters$clusterId
    Bkp <- read.table(paste0(samples_data[i, "path"], "/", samples_data[i, "sample"], ".linx.svs.tsv"), header=T, sep="\t")
    Bkp <- Bkp %>% filter(clusterId %in% cl_ids) %>% dplyr::select(vcfId, clusterId)

    SVs <- readVcf(paste0(samples_data[i, "path"], "/", samples_data[i, "sample"], ".purple.sv.vcf.gz"))
    suppressWarnings({gr <- breakpointRanges(SVs)})
    pairs <- breakpointgr2bedpe(gr)
    pairs <- pairs %>% filter(pairs$name %in% Bkp$vcfId)

    if (nrow(pairs)!=0)
    {
      final <- inner_join(pairs, Bkp, by=c("name"="vcfId"))

      #divide in single breakends
      bk_left <- final %>% dplyr::select(chrom1, start1, end1, clusterId) %>% dplyr::rename("chrom"="chrom1", "start"="start1", "end"="end1")
      bk_right <- final %>% dplyr::select(chrom2, start2, end2, clusterId) %>% dplyr::rename("chrom"="chrom2", "start"="start2", "end"="end2")
      bks <- rbind(bk_left, bk_right)
      bks$clusterId <- as.character(bks$clusterId)

      #if there is at least one breakpoint in the cluster spanning 2 chromosomes it is extrachromsomal
      inorex <- bks %>% group_by(clusterId) %>% summarize(extra = ifelse(n_distinct(chrom)>1, T, F))
      inorex <- inorex %>% as.data.frame()

      #select only intrachromosomal clusters
      intra_clusters <- inorex[inorex$extra==F, "clusterId"]
      intra <- bks %>% filter(clusterId %in% intra_clusters)

      #group_by clusterId and sort by chromosome and start position
      intra <- intra[order(intra[, "clusterId"], intra[,"chrom"], intra[,"start"]),]

      #calculate the interbreakpoint distance
      intra <- intra %>%
        dplyr::group_by(clusterId) %>%
        dplyr::mutate(dist = end - lag(end))

      intra$sample <- samples_data[i, "sample"]
      breakpoints <- rbind(breakpoints, intra)
      #print(paste0(samples_data[i, "sample"], " done"))
    } else {
      print(paste0(samples_data[i, "sample"], " skipped"))
    }
  }
  colnames(breakpoints) <- c("chr", "start", "end", "clusterId", "dist", "sample")
  breakpoints <- breakpoints[!is.na(breakpoints$dist), ] %>% as.data.frame() %>% filter (dist>0)

  return(breakpoints)
}

#' Compact/sparse events
#' @description
#' From the breakpoint distances it defines the threshold between a distant pair of breakpoints and a compact one.
#' The threshold is set using Mclust package, so extracting two Gaussian functions.
#'
#'
#' @param bkps : data frame coming from the function `bk_distance`
#' @param seed : set a random number as seed to keep reproducibility
#'
#' @export
extract_compact_sparse <- function(bkps, seed)
{
  if(!is.na(seed))
  {
    set.seed(seed)
  }

  #calculation of the GMM
  dist_gmm <- suppressMessages(Mclust(bkps$dist, G=2, verbose = F))
  bkps$group <- dist_gmm$classification

  #identify which group corresponds to sparse and compact
  if(min(bkps[bkps$group==1, "dist"]) < min(bkps[bkps$group==2, "dist"]))
  {
    bkps[bkps$group==1, "group"] <- "compact"
    bkps[bkps$group==2, "group"] <- "sparse"
  } else {
    bkps[bkps$group==1, "group"] <- "sparse"
    bkps[bkps$group==2, "group"] <- "compact"
  }

  #Calculate the model of the fit (Gaussian)
  metrics_df <- bkps %>% group_by(group) %>% summarize(mean=mean(dist), sd=sd(dist), .groups = 'drop') %>% as.data.frame()

  #count for each sample the number of compact and sparse events
  complex_evs <- bkps %>% dplyr::select(sample, group) %>% group_by(sample, group) %>% summarise(count = n(), .groups = 'drop')
  compact_sparse <- spread(complex_evs, key = "group", value = "count") %>% as.data.frame()
  compact_sparse[is.na(compact_sparse$compact), "compact"] <- 0
  compact_sparse[is.na(compact_sparse$sparse), "sparse"] <- 0

  results <- c(list(metrics_df), list(compact_sparse))
  return(results)
}

#' Deletion and duplication lengths
#' @description
#' Extracts all the deletions and duplications from all samples clusters (output from Linx)
#'
#' @param samples_data :  data frame constituted by the fields "sample", "patient" and "path"
#'  sample: sample name
#'  patient: patient name
#'  path: path to the directory containing Purple results relative to the specified sample
#'
#' @export
#'
deldup <- function(samples_data) {

  deletions <- data.frame()
  duplications <- data.frame()

  for (i in 1:nrow(samples_data)) {
    sample_i <- samples_data[i, "sample"]

    #Define the paths to the files
    path_clusters <- paste0(samples_data[i, "path"], "/", sample_i, ".linx.clusters.tsv")
    path_linx_svs <- paste0(samples_data[i, "path"], "/", sample_i, ".linx.svs.tsv")
    path_purple_svs <- paste0(samples_data[i, "path"], "/", sample_i, ".purple.sv.vcf.gz")

    #Extraction clusters, breakpoints and SVs
    clusters <- read.table(path_clusters, sep="\t", header=T)
    bks <- read.table(path_linx_svs, sep="\t", header=T)
    svs <- readVcf(path_purple_svs)

    #transform svs VCF to pairs df
    suppressWarnings({svs <- breakpointRanges(svs, inferMissingBreakends = TRUE)})
    pairs <- breakpointgr2bedpe(svs)

    #rename column vcfID and remove the direction of the breakend
    pairs <- pairs %>% as.data.frame() %>% dplyr::rename("vcfId" = "name")
    pairs$vcfId <- substr(pairs$vcfId, 1, nchar(pairs$vcfId) - 1)

    #select duplications and deletions from the clusters
    dups_cl <- clusters %>% filter(resolvedType == "DUP") %>% dplyr::select(clusterId)
    dels_cl <- clusters %>% filter(resolvedType == "DEL" & synthetic == "false") %>% dplyr::select(clusterId)

    #extract for duplications the position in the genome from pairs and calculate the length
    if (nrow(dups_cl) > 0) {
      dups_bk <- bks %>% filter(clusterId %in% dups_cl$clusterId) %>% dplyr::select(vcfId, clusterId)
      dups_bk$vcfId <- substr(dups_bk$vcfId, 1, nchar(dups_bk$vcfId) - 1)
      suppressMessages({dups <- right_join(pairs, dups_bk) %>% mutate(sample = sample_i, length = abs(end2 - start1)) %>% dplyr::select(sample, length)})
      duplications <- rbind(duplications, dups)
    }
    #extract for deletions the position in the genome from pairs and calculate the length
    if (nrow(dels_cl) > 0) {
      dels_bk <- bks %>% filter(clusterId %in% dels_cl$clusterId) %>% dplyr::select(vcfId, clusterId)
      dels_bk$vcfId <- substr(dels_bk$vcfId, 1, nchar(dels_bk$vcfId) - 1)
      suppressMessages({dels <- right_join(pairs, dels_bk) %>% mutate(sample = sample_i, length = abs(end2 - start1)) %>% dplyr::select(sample, length)})
      deletions <- rbind(deletions, dels)
    }

    #print(paste0(sample_i, " done!"))
  }

  return(list(deletions = deletions, duplications = duplications))
}

#' SV categories
#' @description
#' Extracts from Linx cluster results the number of insertions, LINE events, reciprocal inversions, reciprocal translocations
#' and unbalanced translocations.
#'
#' @param samples_data : data frame constituted by the fields "sample", "patient" and "path"
#'  sample: sample name
#'  patient: patient name
#'  path: path to the directory containing Purple results relative to the specified sample
#'
#' @export
#'
SV_categories <- function(samples_data)
{
  SV_classes <- data.frame()
  for(i in 1:nrow(samples_data))
  {
    clusters <- read.table(paste0(samples_data[i, "path"], "/", samples_data[i, "sample"], ".linx.clusters.tsv"), sep='\t', header=T)
    clusters_filtered <- clusters %>% filter(category!= "ARTIFACT" & category != "INCOMPLETE")

    INS <- clusters_filtered %>% filter(resolvedType=="INS") %>% nrow()
    RECIP_INV <-clusters_filtered[grepl("RECIP_INV", clusters_filtered$resolvedType),] %>% nrow()
    RECIP_TRANS <-clusters_filtered[grepl("RECIP_TRANS", clusters_filtered$resolvedType),] %>% nrow()
    UNBAL_TRANS <-clusters_filtered %>% filter(resolvedType=="UNBAL_TRANS") %>% nrow()
    LINE <- clusters_filtered %>% filter(resolvedType=="LINE") %>% nrow()

    SV_classes <- rbind(SV_classes, c(samples_data[i, "sample"], INS, RECIP_INV, RECIP_TRANS, UNBAL_TRANS, LINE))
  }
  colnames(SV_classes) <- c("sample", "ins", "recip_inv", "recip_trans", "unbal_trans", "LINE")
  return(SV_classes)
}

#' LogR quantiles
#' @description
#' Extracts the quantiles 20 and 80 of the logR for each sample
#'
#'
#' @param samples_data : data frame constituted by the fields "sample", "patient" and "path"
#'  sample: sample name
#'  patient: patient name
#'  path: path to the directory containing Purple results relative to the specified sample
#' @param segmentation : data frame derived from the function `extract_segments` containing the CN segmentation for all samples specified in the paramenter `samples_data`.
#'  Minimum requirement of fields:
#'    sample: sample name
#'    chromosome
#'    start: starting point of the segment
#'    end: ending point of the segment
#'    copyNumber
#'
#' @export
logR_quantiles <- function(samples_data, segmentation)
{
  quantiles <- segmentation %>%
    group_by(sample) %>%
    summarise(logR_q20 = quantile(logR, probs=0.2),
              logR_q80 = quantile(logR, probs=0.8))

  # logR_q20 is mainly negative so the idea is turning the signs and shifting
  # all the values so that they are all positive (0 -> -1)
  quantiles$logR_q20neg <- quantiles$logR_q20 * -1
  quantiles$logR_q20neg <- quantiles$logR_q20neg + 1

  # shift logR_q80 in the same way so that they are comparable
  quantiles$logR_q80neg <- quantiles$logR_q80 + 1

  quantiles <- quantiles %>% dplyr::rename("dup_magnitude" = "logR_q80neg",
                                           "del_magnitude" = "logR_q20neg") %>%
    dplyr::select(sample, del_magnitude, dup_magnitude)

  logRs <- quantiles %>% as.data.frame() %>% filter(sample %in% samples_data$sample)
  return(logRs)
}

#' Chromothripsis events
#' @description
#' Function for identification of the presence of chromotripsis in each sample. It uses Shatterseek for identifying
#' the regions eligible for chromotripsis and their criteria for the classification of the events as "high-confidence" or "low-confidence"
#' The number of high-confidence regions for each sample is annotated.
#'
#' @param samples_data : data frame constituted by the fields "sample", "patient" and "path"
#'  sample: sample name
#'  patient: patient name
#'  path: path to the directory containing Purple results relative to the specified sample
#' @param segmentation : data frame derived from the function `extract_segments` containing the CN segmentation for all samples specified in the paramenter `samples_data`.
#'  Minimum requirement of fields:
#'    sample: sample name
#'    chromosome
#'    start: starting point of the segment
#'    end: ending point of the segment
#'    copyNumber
#'
#' @export
#'
is.chromotripsis <- function(samples_data, segmentation)
{
  chromo <- data.frame()
  for (j in 1: nrow(samples_data))
  {
    suppressMessages(suppressWarnings({
      SV_calls <- readVcf(paste0(samples_data[j,"path"], "/", samples_data[j,"sample"], ".purple.sv.vcf.gz"), "hg38")
      gr_calls <- breakpointRanges(SV_calls, inferMissingBreakends=F)
    }))
    bedpe_calls <- breakpointgr2bedpe(gr_calls)

    #column with the type of SV depending on the reciprocal position of the calls
    bedpe_calls <- bedpe_calls %>%
      mutate(svtype=ifelse(chrom1!=chrom2, "TRA",
                           ifelse(strand1=="+" & strand2=="-", "DEL",
                                  ifelse(strand1=="-" & strand2=="+", "DUP",
                                         ifelse(strand1=="-" & strand2=="-", "t2tINV", "h2hINV")))))

    bedpe_calls <- bedpe_calls %>% filter(chrom1 != "chrY" & chrom2!= "chrY")
    #create the Sv_data object
    SV_data <- SVs(chrom1=gsub("chr", "", bedpe_calls$chrom1),
                   pos1=as.numeric(bedpe_calls$start1),
                   chrom2=gsub("chr", "", bedpe_calls$chrom2),
                   pos2=as.numeric(bedpe_calls$end2),
                   SVtype=as.character(bedpe_calls$svtype),
                   strand1=as.character(bedpe_calls$strand1),
                   strand2=as.character(bedpe_calls$strand2))
    #create the CN_data object
    CNVs_calls <- segmentation %>%
      filter(chromosome != "chrY") %>%
      filter(sample == samples_data[j, "sample"])%>%
      mutate(total_cn=round(copyNumber)) %>%
      dplyr::select(chromosome, start, end, total_cn)

    CN_data <- CNVsegs(chrom=gsub("chr", "", CNVs_calls$chromosome),
                       start=CNVs_calls$start,
                       end=CNVs_calls$end,
                       total_cn=CNVs_calls$total_cn)

    suppressWarnings({chr <- shatterseek(SV.sample=SV_data,
                                         seg.sample=CN_data,
                                         genome="hg38")})

    table_chr <- chr@chromSummary
    #removing NA values from max_oscillation_CN_2_states
    table_chr <- table_chr[!is.na(table_chr$max_number_oscillating_CN_segments_2_states),]
    #removing NA from p_val columns
    if (sum(is.na(table_chr$pval_fragment_joins)) >0)
    {
      table_chr[is.na(table_chr$pval_fragment_joins), "pval_fragment_joins"] <- 1
    }
    if (sum(is.na(table_chr$chr_breakpoint_enrichment)) >0)
    {
      table_chr[is.na(table_chr$chr_breakpoint_enrichment), "chr_breakpoint_enrichment"] <- 1
    }
    if(sum(is.na(table_chr$pval_exp_cluster)) >0)
    {
      table_chr[is.na(table_chr$pval_exp_cluster), "pval_exp_cluster"] <- 1
    }

    if(nrow(table_chr)!=0)
    {

      table_chr$category <- ""
      for (i in 1:nrow(table_chr))
      {
        sum <- sum(as.numeric(table_chr[i, c("number_DEL", "number_DUP","number_h2hINV","number_t2tINV")]))
        if((table_chr[i, "max_number_oscillating_CN_segments_2_states"] >= 7 && sum >=6) &&
           (table_chr[i, "pval_fragment_joins"] < 0.05 && (table_chr[i, "chr_breakpoint_enrichment"] < 0.05 || table_chr[i, "pval_exp_cluster"] < 0.05)))
        {
          categ <- "high-confidence"
        }
        else if ((table_chr[i, "max_number_oscillating_CN_segments_2_states"] >= 7 & sum >=3) &
                 (table_chr[i, "pval_fragment_joins"] < 0.05 & table_chr[i, "number_TRA"]>=4))
        {
          categ <- "high-confidence"
        } else if ((table_chr[i, "max_number_oscillating_CN_segments_2_states"] >= 6 & sum >=4) &
                   (table_chr[i, "pval_fragment_joins"] < 0.05 & (table_chr[i, "chr_breakpoint_enrichment"] < 0.05 | table_chr[i, "pval_exp_cluster"] < 0.05)))
        {
          categ <- "low-confidence"
        } else
        {
          categ <- ""
        }
        table_chr[i, "category"] <- categ
      }
    }

    esit <- sum(table_chr$category=="high-confidence")
    chromo <- rbind(chromo, c(samples_data[j, "sample"], esit))
  }
  colnames(chromo) <- c("sample", "chrtps")

  return(chromo)
}
