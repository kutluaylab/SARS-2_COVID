library(VennDiagram)


setwd("C:/Users/NLEE/Desktop/mbio_venndiagram/")

outdir <- "C:/Users/NLEE/Desktop/mbio_venndiagram/outputs/"
dir.create(outdir)

samples <- c("ribo_1","ribo_2","ribo_3","ribo_4",
             "rseq_1","rseq_2","rseq_3","rseq_4")

readAndFilter <- function(species, expr, tp) {
  if (species == "vero") time <- c("2hpi","6hpi","12hpi","24hpi")
  else if (species == "hbec") time <- c("24hpi", "48hpi", "72hpi", "96hpi")
  else (print("Wrong species!"))
  
  fn <- paste(expr,species,time[tp],"diffexgenes","full", sep="_")
  # print(fn)
  table <- read.table(fn, header=T)
  table <- table[(table$FDR < 0.05 & abs(table$logFC) > 1), ]
  print(paste0(dim(table)[1]," DE genes in ",fn))
  return(rownames(table))
}

drawDiagram <- function(vero, hbec, title) {
  plotname <- unlist(strsplit(title, split = "_"))
  venn.diagram(
    x = list(vero, hbec), 
    category.names = c("Vero" , "HBEC"),
    filename = paste0(outdir,title,'.png'), 
    # output specification
    imagetype="png", margin = 1.5,
    height = 500, width = 500, resolution = 300, 
    compression = "lzw",
    # main title
    main = paste0(plotname[1]," timepoint #",plotname[2]),
    main.pos =  c(0.47, 0.625), main.just = c(0.5, 1), main.cex = 0.6,
    main.fontface = "bold", main.fontfamily = "sans",
    # circles
    lwd = 0.5, lty = "blank", col = c("tomato1","steelblue3"), 
    fill = c("tomato1","steelblue3"), 
    # numbers
    fontfamily = "sans", cex = .6,
    ext.line.lwd = 0.3,
    ext.length = 0.85,
    ext.percent = c(0,0,0.1),
    # sets
    cat.cex = 0.6, cat.default.pos = "outer", cat.fontfamily = "sans",
    cat.pos = c(155, 225), cat.dist = c(0.055, 0.055)
  )
}

plotVennDiag <- function(exper_tp) {
  print(exper_tp)
  fname <- unlist(strsplit(exper_tp, split="_"))
  
  vero <- readAndFilter("vero", fname[1], as.numeric(fname[2]))
  hbec <- readAndFilter("hbec", fname[1], as.numeric(fname[2]))

  drawDiagram(vero, hbec, exper_tp)
  
  overlapGenes <- intersect(vero, hbec)
  print(paste0(">> total overlapping genes: ", length(overlapGenes)))
  return(overlapGenes)
}

convListToDf <- function(l, ml) {
  length(l) <- ml
  return(l)
} 

overlapTable <- sapply(samples, function(x) plotVennDiag(x))
ml <- max(sapply(overlapTable, function(x) length(x)))
overlapTable <- sapply(overlapTable, function(x) convListToDf(x, ml))
overlapTable <- as.data.frame(overlapTable)

write.csv(overlapTable, file = paste0(outdir,"overlapGenes.csv"), 
            quote = F, row.names = F, sep = "\t", na = "")