make_output_directories <- function(outputDir) {
  if (!dir.exists(outputDir)) {
    dir.create(outputDir, recursive = T)
  }
  if (!dir.exists(file.path(outputDir, "objs"))) {
    dir.create(file.path(outputDir, "objs"), recursive = T)
  }
  if (!dir.exists(file.path(outputDir, "reports", "figs"))) {
    dir.create(file.path(outputDir, "reports", "figs"), recursive = T)
  }
  if (!dir.exists(file.path(outputDir, "reports", "de_genes"))) {
    dir.create(file.path(outputDir, "reports", "de_genes"), recursive = T)
  }
  
}


make_dge <- function(counts_path, countFilter, region="cds") {
  # function to make dge object for Ribo-seq or RNA-seq and output the counts table and CPM table
  # Args: 
  #   counts_path: path for the output count tables from htseq-count
  #   countFilter: common pattern of the htseq file names in the data type
  #   region: cds, utr3, or utr5
  # Returns:
  #   null. But will save the dge object 
  
  files <- dir(path = counts_path, pattern = "*.count_reduced$")
  
  ## sort the files in asc order
  files <- sort(files)
  counts <- readDGE(files, path = counts_path)$counts
  # filter out metadata fields in counts 
  noint <- rownames(counts) %in% c("__no_feature","__ambiguous","__too_low_aQual",
                                   "__not_aligned","__alignment_not_unique")
  counts <- counts[!noint,]

  # get gene information from ensembl
  ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL", version = "Ensembl Genes 90", dataset = "csabaeus_gene_ensembl", host = "http://aug2017.archive.ensembl.org")
  stopifnot(length(unique(rownames(counts))) == nrow(counts))
  gtf.ens <- getBM(attributes=c('ensembl_gene_id','external_gene_name',"gene_biotype"),
                   filters = 'ensembl_gene_id', values = rownames(counts), mart = ensembl)
  saveRDS(gtf.ens, file.path(file.path(OUTPUT, "objs"), paste0(region, "_genes.rds")))
  
  genes <- readRDS(file.path(file.path(OUTPUT, "objs"), paste0(region, "_genes.rds")))
  genes <- genes[match(rownames(counts), genes$ensembl_gene_id),]
  genes <- na.omit(genes)
  counts <- counts[genes$ensembl_gene_id, ]
  stopifnot(genes$ensembl_gene_id == rownames(counts))

  ## get subset of columns from count table
  counts_sub <- counts[, grep(countFilter, colnames(counts))]

  ## get samples information
  samples <- read_csv("samples.csv")
  # make sure order of the colnames in counts is the same as samples
  counts_sub <- counts_sub[,paste0(samples$sample, "")]
  stopifnot(colnames(counts_sub) == samples$sample)
  write.csv(counts_sub, file.path(file.path(OUTPUT, "reports"), paste0("gene_counts_", region, ".csv")))
  colnames(counts_sub) <- paste(samples$experiment, samples$independent_var, samples$experiment_type, sep = "_")
  dge <- DGEList(counts = counts_sub, 
                 group = paste("time", samples$independent_var, samples$experiment_type, sep = "_"),
                 genes = genes)
  dge$samples$experiment <- samples$experiment
  dge$samples$independent_var <- samples$independent_var
  dge$samples$experiment_type <- samples$experiment_type
  # store only protein coding genes for later analysis
  dge <- dge[dge$genes$gene_biotype=="protein_coding",]
  dge <- na.omit(dge)
  dge <- calcNormFactors(dge)
  saveRDS(dge, file.path(file.path(OUTPUT, "objs"), paste0("dge_", region, "_protein_coding.rds")))
  cpm <- cpm(dge, log = T)
  cpm <- as.data.frame(cpm)
  cpm$ensembl_gene_id <- rownames(cpm)
  write.csv(as.data.frame(cpm), file.path(file.path(OUTPUT, "reports"), paste0("cpm_", region, ".csv")), row.names = F)
  write.csv(colSums(dge$counts), file.path(file.path(OUTPUT, "reports"), paste0("colsums", region, ".csv")))
  write.csv(head(dge$counts,n=15), file.path(file.path(OUTPUT, "reports"), paste0("test_", region, ".csv")))
  
}

generate_reprod_plot <- function(cpmTable, dataType) {
  # function to generate reproducibility plot
  # Args:
  # cpmTable: cpm table created by make_dge
  # dataType: experiment name either rnaseq or riboseq
  
  cor_ <- function(df, col1, col2) {
    stopifnot(c(col1, col2) %in% names(df))
    return(paste0("r = ", format(cor(df[[col1]], df[[col2]]), digits = 3)))
  }
  colors<-list("rnaseq" = "deepskyblue", "riboseq" = "darkorchid2")
  samples <- paste(SAMPLES, dataType, sep = "_")
  comb <- matrix(samples, nrow=2)
  plots <- list()
  for(i in 1:ncol(comb)) { 
    df.cor <- cor_(cpmTable, comb[1, i], comb[2, i])
    plots[[i]] <- ggplot(cpmTable, aes_string(comb[1, i], comb[2, i])) +
      geom_point(alpha=.1, colour = colors[[dataType]]) +
      annotate("text", x=9, y = 0, label = df.cor, size = 4) +
      labs(x=paste(dataType, comb[1, i], sep = " "), y=paste(dataType, comb[2, i], sep = " ")) +
      theme_bw(base_size = 8) 
  }
  if (save_figs) pdf(file.path(OUTPUT, "reports", "figs", paste0(dataType, "_reproducibility.pdf")))
  do.call(grid.arrange, c(plots, ncol=4))
  if (save_figs) dev.off()
}

generate_mds <- function(dge, dataType) {
  # function to generate MDS plot
  # Args:
  #   dge: DGE object create by make_dge
  #   dataType: experiment name either rnaseq or riboseq
  
  colors<-c('deepskyblue', 'darkorchid2', 'seagreen2', 'midnightblue')
  dge_sub <- dge[, grep(dataType, colnames(dge))]
  if (save_figs) pdf(file.path(file.path(OUTPUT, "reports", "figs"), paste0(dataType, "_MDS_dim_1_2.pdf")))
  plotMDS(dge_sub,col=colors, dim.plot = c(1,2))
  #legend('topright', legend=dge_sub$samples$group, col=colors, cex=1, pch=18)
  if (save_figs) dev.off()
  
  if (save_figs) pdf(file.path(file.path(OUTPUT, "reports", "figs"), paste0(dataType, "_MDS_dim_1_3.pdf")))
  plotMDS(dge_sub,col=colors, dim.plot = c(1,3))
  #legend('bottomright', legend=dge_sub$samples$group, col=colors, cex=1, pch=18)
  if (save_figs) dev.off()
  
  if (save_figs) pdf(file.path(file.path(OUTPUT, "reports", "figs"), paste0(dataType, "_MDS_dim_2_3.pdf")))
  plotMDS(dge_sub,col=colors, dim.plot = c(2, 3))
  #legend('bottomright', legend=dge_sub$samples$group, col=colors, cex=1, pch=18)
  if (save_figs) dev.off()
}

plot_CV_vs_CPM <- function(cpmTable, repList, dataType) {
  # function to plot the coefficient variation vs the mean log2 CPM of replications 
  # Args:
  #   cpmTable: CPM table generated by make_dge
  #   repList: specify the list of replications
  #   dataType: either rnaseq or riboseq
  for (reps in repList) {
    reps <- paste(reps, dataType, sep = "_")
    tb <- cpmTable %>% 
      mutate(mean=abs(rowMeans(.[reps])), sd = rowSds(.[reps]), cv = sd/mean,
             reps = paste(reps, collapse = "_vs_"))
    reps_name <- paste(reps, collapse = "_and_")
    if (save_figs) pdf(file.path(file.path(OUTPUT, "reports", "figs"), paste0(reps_name, "_CV_vs_CPM.pdf")))
    pt <- ggplot(tb, aes(x=mean, y=cv)) +
      geom_point(size=0.5) +
      ylim(c(0, 25)) +
      xlim(c(0, 25)) +
    labs(x = "Mean log2 CPM", y = "Coefficient variantion", title = reps_name)
    if (save_figs) dev.off()
    
    print(pt)
  }
  
}

filter_dge <- function(dge, region, threshold, numRep) {
  # function to filter out genes with very low counts across all libraries
  # Args
  #   dge: DGE obejct created by make_dge
  #   region: cds, utr3 or utr5
  #   threshold: cut-off value of log2 CPM for filtering
  #   numRep: minimum number of replicates
  print(dim(dge))
  print(head(cpm(dge, log=T),n=10))
  keep <- rowSums(cpm(dge)>threshold) >= numRep 
  print(table(keep))
  dge <- dge[keep, ,keep.lib.sizes=FALSE]
  print(dim(dge))
  print(head(dge))
  
  # normalizes 
  dge <- calcNormFactors(dge)
  saveRDS(dge, file.path(file.path(OUTPUT, "objs"), paste0("dge_", region, "_filt_norm.rds")))
  cpm <- cpm(dge, log = T)
  cpm <- as.data.frame(cpm)
  cpm$ensembl_gene_id <- rownames(cpm)
  write.csv(cpm, file.path(file.path(OUTPUT, "reports"), paste0("cpm_", region, "_filt.csv")), row.names = F)
}