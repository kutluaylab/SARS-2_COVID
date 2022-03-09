library(tidyverse)
library(edgeR)
library(biomaRt)
library(gridExtra)
library(genefilter)


OUTPUT <- "exp1_results"
## global variable: group of samples 
SAMPLES <- c("RP3_mock", "RP9_mock", "RP3_12hpi", "RP9_12hpi")

## create counts table with all 24 samples
## combined_counts folder should contain all counts_reduced files generated with 3'UTR, 5'UTR and CDS gtfs
counts_path <- file.path(OUTPUT, "combined_counts")
files <- dir(path=counts_path, pattern="*.count_reduced$")

## sort files - asc order
files <- sort(files)
counts <- readDGE(files, path=counts_path)$counts

## filter out metadata fields
# noint <- rownames(counts) %in% c("__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique")
# counts <- counts[!noint,]

## get gene info
ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL", version = "Ensembl Genes 75", dataset = "hsapiens_gene_ensembl", host = "http://feb2014.archive.ensembl.org")
stopifnot(length(unique(rownames(counts))) == nrow(counts))
gtf.ens <- getBM(attributes=c('ensembl_gene_id','external_gene_id',"gene_biotype"),filters = 'ensembl_gene_id', values = rownames(counts), mart = ensembl)
saveRDS(gtf.ens, file.path(file.path(OUTPUT, "objs"), paste0("all", "_genes.rds")))
  
genes <- readRDS(file.path(file.path(OUTPUT, "objs"), paste0("all", "_genes.rds")))
genes <- genes[match(rownames(counts), genes$ensembl_gene_id),]
genes <- na.omit(genes)
counts <- counts[genes$ensembl_gene_id, ]

stopifnot(genes$ensembl_gene_id == rownames(counts))

write.csv(counts, file.path(file.path(OUTPUT, "reports"), paste0("raw_counts_combined", ".csv")))


## saving for each samples
sums <- colSums(counts)

table_12hpi <- sums[grepl("12hpi", names(sums))]
table_12hpi <- matrix((table_12hpi),ncol=4, byrow=TRUE)
colnames(table_12hpi) <- c("Ribo_E32","Ribo_E37","Rseq_E32","Rseq_E37")
rownames(table_12hpi) <- c("cds","utr3","utr5")
write.csv(table_12hpi, file.path("/path/to/output/count_tables/table_12hpi.csv"))

table_mock <- sums[grepl("mock", names(sums))]
table_mock <- matrix((table_mock),ncol=4, byrow=TRUE)
colnames(table_mock) <- c("Ribo_E32","Ribo_E37","Rseq_E32","Rseq_E37")
rownames(table_mock) <- c("cds","utr3","utr5")
write.csv(table_mock, file.path("/path/to/output/count_tables/table_mock.csv"))
