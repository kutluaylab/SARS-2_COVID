library(edgeR)
library(tidyverse)
library(wordcloud)
library(org.Hs.eg.db)

#setwd("/scratch/skutluay/Files/Outputs/covid_te/even_more_updated/48hpi/")
setwd("/scratch/skutluay/Files/Outputs/Experiment32-37_IFIT/te/de_genes_list/")
direc <- '/scratch/skutluay/Files/Outputs/Experiment32-37_IFIT/te/'

#files <- c("ribo_48hpi", "rna_48hpi")
files <- c("ribo_IFIT1", "rseq_IFIT1")
hpi <- strsplit(files, "_")[[1]][2]
DG <- readDGE(files, header=FALSE);
#keep <- rowSums(cpm(DG)>1) >= 1;
#DG <- DG[keep, ,keep.lib.sizes=FALSE];

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
	#plot(dat[1:(nrow(dat) - 9),2], dat[1:(nrow(dat) - 9),3], xlab="logCPM (RNASeq) ", main="Vero Translational Efficiencies, 48hpi", ylab="logCPM (Riboseq / RNASeq) ", ylim=c(-10,12), xlim=c(-5,20), pch=20, col="gray")  
	plot(dat[,2], dat[,3], xlab="logCPM (RNASeq) ", main=paste0("293T Translational Efficiencies, ", hpi), ylab="logCPM (Riboseq / RNASeq) ", ylim=c(-10,12), xlim=c(-5,20), pch=20, col="gray", cex.lab=1.5, cex.axis=1.2, cex.main=2)  
	abline(h=colMeans(dat)[3], col="blue", lty=2, lwd=3)

	#te_genes <- c("PDPR", "MGAT5", "JUN", "CBL", "TRAM2", "FOS", "ATF3", "CLK1", "NLRP1", "CXCL1", "TLE3", "NR4A3", "DDX25", "RPS2", "TMA7", "RPL12", "RPL36A", "OASL", "IRGQ", "NSA2", "CXCL3", "CXCL8", "CXCL10", "CXCL11", "PMAIP1", "NFKBIA", "EEF1A1", "REPS2", "IFIH1", "DDX58", "IL1A")
	#te_genes <- c("PDPR", "JUN", "TRAM2", "FOS", "TLE3", "NR4A3", "DDX25", "RPS2", "TMA7", "RPL12", "RPL36A", "OASL", "IRGQ", "NSA2", "CXCL10", "PMAIP1", "NFKBIA", "EEF1A1", "REPS2", "IFIH1", "DDX58", "IL1A")

	points(dat[de,2], dat[de,3], col="orange", pch=20)
	
	#te_genes <- c("IL23A", "CCL25", "MGAT5", "JUN", "CBL", "TRAM2", "FOS", "ATF3", "LAMP3", "SSTR5", "CXCL1", "IFITM3", "NR4A3", "DDX25", "RPS2", "TMA7", "RPL12", "RPL36A", "OASL", "IRGQ", "NSA2", "CXCL3", "CXCL8", "CXCL10", "CXCL11", "PMAIP1", "NFKBIA", "EEF1A1", "REPS2", "IFIH1", "DDX58", "IL1A", "IL20RA", "CXCL1")
	te_genes <- rownames(genes)[1:150]
	genes.sub7 <- genes[te_genes,]
	write.table(genes.sub7, "top_150", quote=FALSE)
	te_genes <- intersect(te_genes, rownames(dat))

	points(dat[te_genes,2], dat[te_genes,3], col="blue", pch=20)
	#textplot(dat[te_genes,2], dat[te_genes,3], rownames(dat[te_genes, ]),ylim=c(-10,12), show.lines=F, xlim=c(-5,20), new=F, pos=1, cex=1)



	#viral_points <- dat[(nrow(dat) - 9):nrow(dat), ]
	#points(viral_points[,2], viral_points[,3], col="red", pch=19)
	#textplot(viral_points[,2], viral_points[,3], rownames(viral_points),ylim=c(-10,12), show.lines=F, xlim=c(-5,20), new=F, pos=2, cex=1)


	dev.off()



# 1. call the top 150 te genes of LPCX
# 2. find the overlap of 1 and IFIT top 150 te genes
# 3. plot them in blue dots



























# infection_table <- read_csv("cpm_cds_te.csv")
# infection_table <- as.data.frame(infection_table)


# plot(infection_table[,4], infection_table[,7], xlab="logCPM (RNASeq) ", main="COVID and Vero Translational Efficiencies", ylab="logCPM (Riboseq / RNASeq) ", ylim=c(-5,10), xlim=c(-1,20), pch=20, col="gray")  
# points(infection_table[,5], infection_table[,8], pch=20, col="gray") 
# points(infection_table[,6], infection_table[,9], pch=20, col="gray") 


# te_table <- infection_table[,grep("te", colnames(infection_table))]  
# te_table <- na.omit(te_table)

# average_te <- mean(colMeans(infection_table[,7:9])) #te_table wasn't working
# print(average_te)
# abline(h=average_te, col="blue", lty=2, lwd=3)

# covid_counts <- read.table("cpm_covid_reduced_all.csv", header=TRUE, row.names=1, sep=",")
# print(covid_counts)
# #covid_counts <- data.frame(covid_counts[,2:4], row.names=covid_counts[,1])

# points(covid_counts[,2], covid_counts[,3], col="red", pch=19)
# text(covid_counts[,2], covid_counts[,3], labels=row.names(covid_counts), cex=0.8, pos=2)

# dev.off()

# scatter <- ggplot(infection_table, aes(names(infection_table)[2], names(infection_table)[4])) + geom_point() + ylim(-5,10) + xlim(-1,15)

# for (i in seq(5,length(colnames(infection_table)), 3)) {
# 	scatter <- scatter + geom_point(data=data.frame(names(infection_table)[i], names(infection_table)[(i+2)]))
# }



# scatter <- scatter + geom_point(data=covid_counts, aes(data.frame(covid_counts[,1], covid_counts[,3])), color="red")


# scatter <- scatter + geom_hline(aes(yintercept=average_te), color="blue", linetype="dashed", size=1)



