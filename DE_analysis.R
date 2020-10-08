library(edgeR)
library(statmod)
library(RColorBrewer)
library(gplots)
library(org.Hs.eg.db)


directory <- "/counts_files" # directory containing counts files
files <- c("RP3_12hpi","RP3_mock","RP9_12hpi","RP9_mock")
donor_list <- c("RP3","RP3","RP9","RP9")
conditions <- c("12hpi","mock","12hpi","mock")
n <- 200  # number of genes to display on output plots
prefix <- "12hpi_"


setwd(directory)


## Create a DGE list
labels <- paste(donor_list, conditions, sep="_")
DG <- readDGE(files, header=FALSE, labels=labels)
keep <- rowSums(cpm(DG)>1) >= 2
print(cpm(DG)[1:20,1:4])
DG <- DG[keep, ,keep.lib.sizes=FALSE]


## Convert ensembl gene ids to gene names
if(grepl("ENSG", rownames(DG$counts)[1])) {
	symbols <- mapIds(org.Hs.eg.db, keys = rownames(DG$counts), keytype = "ENSEMBL", column="SYMBOL")
	rownames(DG) <- symbols
	symbols_to_keep <- na.omit(symbols)

	DG <- DG[symbols_to_keep,]
}


## Normalizing factor
DG <- calcNormFactors(DG)
DGgroups <- c(conditions)
DG <- DGEList(counts=DG, group=factor(DGgroups))
DG <- calcNormFactors(DG)
Donor <- factor(c(donor_list))
Condition <- factor(c(conditions))


## MDS plot
mdsfilename <- paste(prefix, "_MDSplot.pdf", sep="")

pdf(mdsfilename)

if(length(labels) > 2) {
  colors <- brewer.pal(length(labels),"Set1")
} else {
  colors <- c('deepskyblue', 'red')
}

names(colors) <- labels
plotMDS(DG,col=colors, pch=18, cex=2)
legend('topright', title="Experiment", legend=names(colors), col=colors, cex=.75, pch=18, xpd=T, bty="n", inset=c(0, -.162))

dev.off()


## Dispersion
data.frame(Sample=labels, Donor, Condition)
design <- model.matrix(~Donor+Condition)
DG <- estimateDisp(DG, design, robust=TRUE)


## BCV plot
bcv_file_name <- paste(prefix, "_BCVplot.pdf", sep="")
pdf(bcv_file_name) #name of file to write bcv graph to
plotBCV(DG) #!labels appears to be not working!
dev.off()


## DE anlaysis
de <- exactTest(DG, pair=rev(levels(Condition))) 
de.genes <- rownames(topTags(de, n=n)$table)
diffex_file <- paste(prefix, "_diffexgenes", sep="")
diffex_file_full <- paste(prefix, "_diffexgenes_full", sep="")
write.table(de.genes, diffex_file, quote=FALSE, row.names=FALSE, col.names=TRUE) #file to write de.genes to
write.table(topTags(de, n=n)$table, diffex_file_full, quote=FALSE, row.names=TRUE, col.names=TRUE)
diffextable <- topTags(de, n=n)$table


## Smear plot
smearplotname <- paste(prefix, "_SmearPlot.pdf", sep="")
pdf(smearplotname)
plotSmear(DG, de.tags=de.genes) #Generates a smear plot for differentially expressed genes, upregulated and downregulated
dev.off()


## Volcano Plot
detags <- topTags(de, n=n) 

volcanodata <- cbind(detags$table$logFC, -log10(detags$table$PValue), detags$table$FDR)
rownames(volcanodata) <- rownames(detags$table)
colnames(volcanodata) <- c("logFC", "negLogpVal", "FDR")

volcanodata_insig <- volcanodata[volcanodata[,2] < -log10(.05),]
volcanodata_sig <- volcanodata[volcanodata[,2] >= -log10(.05),]
volcano_logFC <- volcanodata_sig[abs(volcanodata_sig[,1]) >= 2,] 
volcanoplotname <- paste(prefix, "_VolcanoPlot.png", sep="")

windowsize=10
upperlim=10

png(volcanoplotname, height=1000, width=1000)
par(mar = c(5,5,5,5))
plot(volcanodata_insig, col="red", pch=19, cex=1.5, ylim=c(0,upperlim), xlim=c(-windowsize, windowsize), cex.lab=2, cex.axis=1.5, font.lab=2, font=2) #Generates a volcano plot of the DE data
points(volcanodata_sig, col="green", cex=1.5, pch=19)
points(volcano_logFC, col="blue", cex=1.5, pch=19)

print(dim(volcanodata))
sig_FDR <- volcanodata[volcanodata[,3] < .05 & volcanodata[,2] < upperlim, ]
sig_FDR <- sig_FDR[order(abs(sig_FDR[,1]), decreasing=TRUE), ]
sig_FDR <- sig_FDR[sig_FDR[,1] < windowsize, ]

num_gene = min(15, length(rownames(sig_FDR)))
text(sig_FDR[1:num_gene,1], sig_FDR[1:num_gene,2] , labels=row.names(sig_FDR[1:num_gene,]), cex=1.2, pos=2, font=2) #Change cex to 2.5 for covid
abline(h=-log10(.05), col="blue", lty=2, lwd=3)
dev.off()


## Heatmap with top DE genes
y <- cpm(DG, log=TRUE, prior.count = 1)
write.table(y, "cpm", quote=FALSE, row.names=TRUE, col.names=TRUE)
volcdata_sub = volcanodata[rownames(volcanodata) %in% rownames(y),]
max_rows <- rownames(volcdata_sub[order(abs(volcdata_sub[volcdata_sub[,3] < .05, ][,1]), decreasing=TRUE), ])
head(max_rows)
selY <- y[max_rows, ]
selY <- selY[1:min(40, length(rownames(selY))),]

write.table(selY, "heatmap_table", quote=FALSE)

heatmapname <- paste(prefix, "_HeatMap.png", sep="")
png(heatmapname, height=1000, width=1000) 
heatmap.2(selY, col=bluered(100), , trace="none", density.info="none", symkey=F, scale="none", cexRow=1.5,cexCol=2,margins=c(14,16), Colv=FALSE)
dev.off()


## Collapsed Heatmap
cols <- unique(conditions)
collapsed <- matrix(nrow=nrow(selY), ncol=length(cols))
if("mock" %in% cols) {
	cols <- c("mock", cols[cols != "mock"])
}
colnames(collapsed) <- cols
cols <- colnames(collapsed)
rownames(collapsed) <- rownames(selY)

for(i in seq(1, length(cols))) {
	collapsed[,i] <- rowMeans(selY[,grepl(cols[i], colnames(selY))])
}

png('collapsed_heatmap.png', height=1000, width=1000)
heatmap.2(collapsed, col=bluered(100), , trace="none", density.info="none", symkey=F, scale="none", cexRow=1.5,cexCol=2,margins=c(14,16), Colv=FALSE)
dev.off()
