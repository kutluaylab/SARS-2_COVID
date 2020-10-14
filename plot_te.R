library(edgeR)
library(tidyverse)
library(wordcloud)
library(org.Hs.eg.db)

setwd("/directory/with/input/files")

# input count files 
files <- c("ribo_12hpi", "rna_12hpi")
hpi <- strsplit(files, "_")[[1]][2]
DG <- readDGE(files, header=FALSE);
keep <- rowSums(cpm(DG)>1) >= 1;
DG <- DG[keep, ,keep.lib.sizes=FALSE];

if(grepl("ENSG", rownames(DG$counts)[1])) {
	symbols <- mapIds(org.Hs.eg.db, keys = rownames(DG$counts), keytype = "ENSEMBL", column="SYMBOL")
	rownames(DG) <- symbols
	symbols_to_keep <- na.omit(symbols)
	DG <- DG[symbols_to_keep,]
}

DG <- calcNormFactors(DG);

cpm <- cpm(DG, log=T)
dat <- as.data.frame(cpm)
colnames(dat) <- files
dat$te <- dat[,1] - dat[,2] 

genes <- dat[order(abs(dat$te), decreasing=TRUE), ]
te_save <- paste(direc, "te_genes", sep="")
write.table(genes, te_save, quote=FALSE)

#ribo_de <- read.table("DE_48hpi_ribo", row.names=1, header=TRUE)
ribo_de <- read.table("ribo_IFIT1_diffexgenes_full", row.names=1, header=TRUE)
ribo_de <- rownames(ribo_de[ribo_de$FDR < .05 & abs(ribo_de$logFC) > 1,])
#rseq_de <- read.table("DE_48hpi_rseq", row.names=1, header=TRUE)
rseq_de <- read.table("rseq_IFIT1_diffexgenes_full", row.names=1, header=TRUE)
rseq_de <- rownames(rseq_de[rseq_de$FDR < .05 & abs(rseq_de$logFC) > 1,])

de <- union(ribo_de, rseq_de)
rna_de <- rseq_de[!(rseq_de %in% ribo_de)]

write.table(dat, "dat", quote=FALSE)

	teplot <- paste(direc, "translational_efficiency.png", sep="")
	png(teplot, width=700, height=700)
	par(mgp=c(2.5,1,0))  
	plot(dat[1:(nrow(dat) - 9),2], dat[1:(nrow(dat) - 9),3], xlab="logCPM (RNASeq) ", main="Vero Translational Efficiencies, 48hpi", ylab="logCPM (Riboseq / RNASeq) ", ylim=c(-10,12), xlim=c(-5,20), pch=20, col="gray")  
	abline(h=colMeans(dat)[3], col="blue", lty=2, lwd=3)

	# manually selected gene list
	te_genes <- c("PDPR", "MGAT5", "JUN", "CBL", "TRAM2", "FOS", "ATF3", "CLK1", "NLRP1", "CXCL1", "TLE3", "NR4A3", "DDX25", "RPS2", "TMA7", "RPL12", "RPL36A", "OASL", "IRGQ", "NSA2", "CXCL3", "CXCL8", "CXCL10", "CXCL11", "PMAIP1", "NFKBIA", "EEF1A1", "REPS2", "IFIH1", "DDX58", "IL1A")

	points(dat[de,2], dat[de,3], col="orange", pch=20)
	
	te_genes <- rownames(genes)[1:150]
	genes.sub7 <- genes[te_genes,]
	write.table(genes.sub7, "top_150", quote=FALSE)
	te_genes <- intersect(te_genes, rownames(dat))

	points(dat[te_genes,2], dat[te_genes,3], col="blue", pch=20)
	textplot(dat[te_genes,2], dat[te_genes,3], rownames(dat[te_genes, ]),ylim=c(-10,12), show.lines=F, xlim=c(-5,20), new=F, pos=1, cex=1)

	viral_points <- dat[(nrow(dat) - 9):nrow(dat), ]
	points(viral_points[,2], viral_points[,3], col="red", pch=19)
	textplot(viral_points[,2], viral_points[,3], rownames(viral_points),ylim=c(-10,12), show.lines=F, xlim=c(-5,20), new=F, pos=2, cex=1)


	dev.off()
