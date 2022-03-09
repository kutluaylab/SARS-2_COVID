### Manual delta TE Calculation Using a Combined Table of RNA-seq and Ribo-seq Log2CPM ###

library(edgeR)


## Loading Data ##

# Paths to Ribo-seq and RNA-seq count files of all replicates and all time points #
# Count files should be named as such: "experiment_condition_timepoint". For example, 162_infected_96hpi and 162_mock_4h.
ribo_count_path <- "/path/to/riboseq/counts"
rseq_count_path <- "/path/to/rnaseq/counts"

# gtf file for gene name conversion from ensembl id to gene symbols #
# This gtf file needs to be gene-only and have the same number of genes as the count files;
# this can be achieved by only selecting rows which have "gene" as the third field (instead of "exon" or anything else).
# The genes in this gtf also need to be in the same order as the genes in count files.
gtf_path <- "/path/to/gene_only_gtf_file"

outputpath <- "/output/path/of/choice"

# Timepoint in which TE is being calculated #
# HBEC timepoints: 24hpi, 48hpi, 72hpi, 96hpi; Vero timepoints: 2hpi, 6hpi, 12hpi, 24hpi.
time <- "96hpi"

# Number of replicates #
numrep <- 4 # 4 for HBEC, 3 for Vero


ribo_infected <- list.files(ribo_count_path, pattern = paste0("*infected_", time), full.names = T)
ribo_mock <- list.files(ribo_count_path, pattern = "*mock_96h", full.names = T) # <<< change timepoint accordingly

rseq_infected <- list.files(rseq_count_path, pattern = paste0("*infected_", time), full.names = T)
rseq_mock <- list.files(rseq_count_path, pattern = "*mock_96h", full.names = T) # <<< change timepoint accordingly


# Vero
# labels <- c("ribo_176_mock", "ribo_RP3_mock", "ribo_RP9_mock",
#             "ribo_176_infected", "ribo_RP3_infected", "ribo_RP9_infected",
#             "rseq_176_mock", "rseq_RP3_mock", "rseq_RP9_mock",
#             "rseq_176_infected", "rseq_RP3_infected", "rseq_RP9_infected")

# HBEC
labels <- c("ribo_162-1_mock", "ribo_162-2_mock", "ribo_163-1_mock", "ribo_163-2_mock",
            "ribo_162-1_infected", "ribo_162-2_infected", "ribo_163-1_infected", "ribo_163-2_infected",
            "rseq_162-1_mock", "rseq_162-2_mock", "rseq_163-1_mock", "rseq_163-2_mock",
            "rseq_162-1_infected", "rseq_162-2_infected", "rseq_163-1_infected", "rseq_163-2_infected")
dge <- readDGE(c(ribo_mock, ribo_infected, rseq_mock, rseq_infected), header = F, labels = labels)


## Human Gene Name Conversion from Ensembl to Gene Names ##
if (grepl(pattern = "ENSG", rownames(dge$counts)[1])) {
  gtf <- read.table(gtf_path, sep = " ")
  gene_names <- as.character(gtf[, 6])
  gene_names <- substr(gene_names, start = 1, stop = nchar(gene_names) - 1) # getting rid of the ";" at the end
  stopifnot(length(gene_names) == nrow(dge))
  rownames(dge) <- gene_names
}


## Filtering ##
keep <- rowSums(cpm(dge) > 1) >= numrep # filtering by expression level
print(paste0("before filtering: ", nrow(dge)))

dge <- dge[keep, ,keep.lib.sizes=FALSE]
print(paste0("after filtering: ", nrow(dge)))


## Deduplication ##
duplicate_genes <- row.names(dge)[duplicated(row.names(dge))] # getting rid of duplicate genes in cpm table, including the original one
dge <- dge[!(row.names(dge) %in% duplicate_genes), ]
print(paste0("after deduplication: ", nrow(dge)))


## CalcNorm and CPM ##
dge <- calcNormFactors(dge)

cpm <- cpm(dge, log = T, prior.count = 1)

write.csv(cpm, file.path(outputpath, paste0(time, "_combined_cpm.csv")))


## TE Calculation ##
te_table <- NULL

for (i in seq(numrep * 2)) {
  te <- cpm[, i] - cpm[, numrep * 2 + i] # substracting rseq cpm from ribo cpm to get TE
  te_table <- cbind(te_table, te)
}

# Vero
# colnames(te_table) <- paste0(c("176_mock", "RP3_mock", "RP9_mock",
#                                "176_infected", "RP3_infected", "RP9_infected"),
#                              "_te")

# HBEC
colnames(te_table) <- paste0(c("162-1_mock", "162-2_mock", "163-1_mock", "163-2_mock",
                               "162-1_infected", "162-2_infected", "163-1_infected", "163-2_infected"),
                             "_te")


## Mean TE and delta TE Calculation ##
te_table <- data.frame(te_table)
te_table[, "mock_mean"] <- rowMeans(te_table[, seq(numrep)])
te_table[, "infected_mean"] <- rowMeans(te_table[, seq(numrep + 1, numrep * 2)])
te_table[, "delta_te"] <- te_table[, "infected_mean"] - te_table[, "mock_mean"]

write.csv(te_table, file = file.path(outputpath, paste0(time, "_dTE.csv")))
