### Scatterplot of edgeR RNA-seq logFC vs. manual delta TE, labeling interesting genes ###

library(ggplot2)
library(ggrepel)


## Loading Data ##
time <- "96hpi"
diffex <- read.delim(file = "C:/kutluaylab/data/HBEC_162_163/rnaseq/edgeR/96hpi_vs_mock_96h/96hpi_vs_mock_96h_diffexgenes_full")

# dTE <- read.csv(file = "C:/kutluaylab/data/Vero/Yating/24hpi/exp1_results/reports/cds_table_all_te.csv",
#                 row.names = 1)
dTE <- read.csv(file = "C:/kutluaylab/data/HBEC_162_163/overall/TE/96hpi_dTE.csv",
                row.names = 1)

output_path <- "C:/kutluaylab/data/HBEC_162_163/overall/rseq_logFC_vs_manual_dTE"
output_name <- "96hpi_logFC_dTE"

# genes <- intersect(rownames(diffex), dTE$external_gene_id)
genes <- intersect(row.names(diffex), row.names(dTE))

# diffex_dTE <- data.frame(diffex[genes, ]$logFC, diffex[genes, ]$FDR, dTE[match(genes, dTE$external_gene_id), ]$delta_te)
diffex_dTE <- data.frame(diffex[genes, ]$logFC,
                         diffex[genes, ]$FDR,
                         dTE[genes, ]$delta_te)

row.names(diffex_dTE) <- genes
colnames(diffex_dTE) <- c("RNAseq_logFC", "RNAseq_FDR", "delta_TE")

diffex_dTE$significant <- if_else(condition = (abs(diffex_dTE[, "RNAseq_logFC"]) > 2 & diffex_dTE[, "RNAseq_FDR"] < 0.1),
                                  true = "true", false = "false")

write.csv(diffex_dTE, file.path(output_path, paste0(output_name, ".csv")))




## Labeling Interesting Genes ##

## Vero ##
# label_2hpi <- c("EGR3", "NFKBIA", "FOSB")
# label_6hpi <- c("CXCL11", "CXCL3", "CXCL1", "IFIT1", "EGR1", "TNFAIP3")
# label_12hpi <- c("CXCL11", "OASL", "NR4A3", "IFNL1", "CXCL8", "NFKBIA",
#                  "CXCL3", "MX1", "CCL5", "CXCL1", "IFIT1", "IFIT2", "IFIT3",
#                  "USP18", "IL1A", "IRF1", "ISG15", "IFIT5")
# label_24hpi <- c("NR4A3", "OASL", "FOS", "CXCL8", "CXCL11", "NFKBIA",
#                  "CXCL3", "IFI44", "ISG15", "EGR1", "EGR2", "IFIT1",
#                  "IFIT2", "IFIT3", "CXCR4", "IL6", "IRF1", "IL11")

## HBEC ##
label_24hpi <- c("OASL", "CXCL11", "CXCL10", "MX2", "ISG15", "IFI27",
                 "CXCL9", "RSAD2", "IFI6", "OAS3", "CMPK2", "HELZ2",
                 "IDO1", "CCL2", "OAS1", "IFIT1", "IFIT2", "IFIT3")

label_48hpi <- c("CXCL10", "CXCL11", "OASL", "CXCL9", "MX2", "IFNB1",
                 "CMPK2", "OAS3", "IFIT1", "IFI44L", "IFIT2", "IFI44",
                 "OAS1", "IFI6", "IFIT3", "STAT2", "MX1", "ISG15", "IDO1")

label_72hpi_old <- c("CXCL10", "OASL", "CXCL9", "IFNB1", "MX2", "IFNL1",
                 "OAS3", "CXCL11", "CXCL6", "RSAD2", "IFIT2", "IFIT1",
                 "OAS1", "IFIT3", "CMPK2", "OAS2", "IL8", "MX1",
                 "IDO1", "IFI44", "IFI6", "IFI44L", "ISG15", "STAT1",
                 "STAT2", "CXCL2", "KCNE4", "LGALS9B", "SPTBN5", "NR4A1",
                 "EGR3", "CD34", "NR1D1", "NR4A3")

label_72hpi <- c("CXCL10", "OASL", "CXCL9", "IFNB1", "MX2", "IFNL1",
                 "OAS3", "CXCL11", "CXCL6", "RSAD2", "IFIT2", "IFIT1",
                 "OAS1", "IFIT3", "CMPK2", "OAS2", "IL8", "MX1",
                 "IDO1", "IFI44", "IFI6", "IFI44L", "ISG15", "STAT1",
                 "STAT2", "CXCL2", "NR4A1", "EGR3", "NR4A3")


label_96hpi <- c("IFI6", "XAF1", "IFITM1", "IFI44L", "AIM2", "OAS1",
                 "IFI44", "CCL20", "IFIT1", "OAS2", "OAS3", "IFIH1",
                 "IFIT3", "CMPK2", "OASL", "IDO1", "ISG15", "IFNL1",
                 "MX2", "IFIT5", "IL6", "CXCL11", "CXCL12", "CXCL3",
                 "CXCL6", "CXCL9", "CXCL2", "SELL", "CARD17", "NOS2", "TNFAIP6")

label_current <- label_96hpi

time <- "96hpi"
diffex_dTE <- read.csv(file = paste0("C:/kutluaylab/data/HBEC_162_163/overall/rseq_logFC_vs_manual_dTE/",
                                     time, "_logFC_dTE.csv"),
                       row.names = 1)

diffex_dTE$label <- ""
diffex_dTE[label_current, ]$label <- row.names(diffex_dTE[label_current, ])


# Removing unwanted rows #
# row.names.remove <- "AF127936.5" # HBEC 72hpi
row.names.remove <- c("ZEB1", "ZNF783", "AF127936.5") # HBEC 96hpi
diffex_dTE <- diffex_dTE[!row.names(diffex_dTE) %in% row.names.remove, ]


## Scatterplot ##

output_name <- paste0(time, "_logFC_dTE")


p <- (ggplot(data = diffex_dTE, mapping = aes(x = RNAseq_logFC, y = delta_TE, label = label))
  + geom_point(mapping = aes(color = significant),
               size = 1,
               alpha = 0.8)
  + geom_text_repel(size = 2.5, max.overlaps = 100, max.iter = 10000, max.time = 2)
  + xlab(label = "RNAseq logFC")
  + ylab(label = "delta TE")
  + scale_color_manual(labels = c("insignificant change in transcription", "significant change in transcription"),
                       values = c("turquoise", "orange"))
  + theme_classic()
  + theme(legend.title = element_blank(),
          legend.position = c(0.75, 0.95))
)


png(filename = file.path(output_path, paste0(output_name, ".png")), height = 1500, width = 1500, res = 300)
print(p)
dev.off()
