# ---
# title: "viral_density_plot"
# author: "NaKyung Lee"
# date: "3/12/2021"
# output: html_document
# ---

library(tidyverse)
library(RColorBrewer)
library(pals)
library(ggthemes)
library(ggrepel)
library(gtable)
library(grid)
library(gridExtra)


setwd('C:/Users/nlee1/Documents/lab_data/COVID/HBEC_ribo/162-1/162-1_repaired_counts/')
dir.create('plots')

### Read in the annotaion file
anno = read_tsv('C:/Users/nlee1/Documents/lab_data/scripts/viral_density_test/COVID-19_sepORF1AB.gff')
anno = data.frame(anno)

cds = anno[anno$Feature != "gene", 1:8]
cds = cds[cds$Strand == "+",]
cds$Frame = as.numeric(cds$Frame)

cds$Feature[which(cds$Feature=='five_prime_UTR')] = "5'UTR"
cds$Feature[which(cds$Feature=='three_prime_UTR')] = "3'UTR"
cds$Name[which(cds$Name == 'RefSeq')] = cds$Feature[which(cds$Name=='RefSeq')]


### Change the range of ORF1A
### Note: there is an overlapping region between ORF1A and ORF1B
### the overlapping region will be labeled ORF1B for the purpose of plotting
cds$End[2] = cds$Start[3]-1


### Read in the counts_repaired file
filelist = list.files(pattern='*_sense_counts_repaired')
#filelist = paste('repaired_counts',filelist,sep='/')

for (i in 1:length(filelist)) {
filename = unlist(strsplit(filelist[i],'_'))
filename = paste(filename[1],filename[2],filename[3],sep='_')
savepath = paste('plots',filename,sep='/')
print(filename)

sense = read.table(filelist[i], header=TRUE, sep='\t')
sense = data.frame(sense$bp, rowSums(sense[2:5]))
colnames(sense) = c('bp','count')

### Build Plotting Data
plotting.data = data.frame(bp = 1:tail(cds$End,1))

plotting.data = merge(plotting.data, sense, all=TRUE)
colnames(plotting.data)[2] = 'sense_count'

plotting.data$sense_count[is.na(plotting.data$sense_count)] <- 0
plotting.data$gene = ''

for (j in 1:length(cds$Name)) {
  plotting.data$gene[cds$Start[j]:cds$End[j]] = cds$Name[j]
}

y_max = 1.2*max(sense$count) 
y_lim = max(log(plotting.data$sense_count+1))


### PLOT TYPE I

#mixPalette = c('#FFFFFF', sample(as.vector(tableau20(14))))
#colPalette = mixPalette
colPalette = c("#FFFFFF", "#AEC7E87F", "#98DF8A7F", "#FF98967F", 
               "#C49C947F", "#C5B0D57F", "#FFBB787F", "#D627287F", 
               "#FF7F0E7F", "#9467BD7F", "#8C564B7F", "#E377C27F", 
               "#F7B6D27F", "#2CA02C7F", "#1F77B47F")


### Plotting Sense-counts
ps <- (ggplot(data = plotting.data, aes(x = bp, y = log10(sense_count+1))) 
      #+ geom_line(aes(color = gene)) 
      + geom_area(aes(fill = gene))
      + scale_fill_manual(breaks = c('',cds$Name),
                          values = colPalette,
                          labels = c('',cds$Name))
      + theme_classic()
      + xlim(0,length(plotting.data$bp))
      + labs(y = 'log(counts)')
      + theme(legend.position = 'none',
              axis.title.x = element_blank(),
              # axis.title.y = element_text(family='sans', size=12),
              axis.title.y = element_text(size=12),
              axis.text.x = element_text(size=12),
              axis.text.y = element_text(size=12))
      + scale_x_continuous(breaks=seq(0,30000,5000), labels=seq(0,30000,5000))
      )

psl <- (ggplot(data = plotting.data, aes(x = bp, y = sense_count)) 
       #+ geom_line(aes(color = gene)) 
       + geom_area(aes(fill = gene))
       + scale_fill_manual(breaks = c('',cds$Name),
                           values = colPalette,
                           labels = c('',cds$Name))
       + theme_classic()
       + xlim(0,length(plotting.data$bp))
       + labs(y = 'counts')
       + theme(legend.position = 'none',
               axis.title.x = element_blank(),
               # axis.title.y = element_text(family='sans', size=12),
               axis.title.y = element_text(size=12),
               axis.text.x = element_text(size=12),
               axis.text.y = element_text(size=12))
       + scale_x_continuous(breaks=seq(0,30000,5000), labels=seq(0,30000,5000))
)

# ps

### Plotting Genome CDS regions
gme <- (ggplot(data = plotting.data, aes(x = bp))
      + theme_tufte()
      + xlim(0,length(plotting.data$bp))
      + ylim(0,1)
      + geom_rect(aes(NULL,NULL,xmin=bp-1,xmax=bp,fill=gene),
                  ymin=0, ymax=0.1)
      # + geom_rect(data=cds, aes(NULL,NULL,xmin=Start,xmax=End,fill=Name),
      #             ymin=0, ymax=0.1, alpha=0.5)
      + scale_fill_manual(breaks = (c('',cds$Name)),
                          values = colPalette,
                          labels = (c('',cds$Name)))
      + labs(title = 'Viral Density')
      + theme(plot.title = element_text(family='sans', size=16, face='bold', hjust=0.5),
      # + theme(plot.title = element_text(size=16, face='bold', hjust=0.5),
              axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              legend.position = c(0.5,0.3),
              legend.direction = 'horizontal',
              legend.title = element_blank(),
              legend.key.size = unit(0.5,'cm'),
              legend.text = element_text(family='sans', size=12))
              # legend.text = element_text(size=12))
      + guides(fill = guide_legend(nrow = 1, order = 2))
      )



### Combining plots
g1 <- ggplotGrob(gme)
g2 <- ggplotGrob(ps)
#g3 <- ggplotGrob(psl)


### log-normed plot
g <- rbind(g1, g2, size = "first")
g$widths <- unit.pmax(g1$widths, g2$widths)

png(paste(savepath,'_vgenome.png',sep=''), 12, 3, units='in', res=300)
grid.newpage()
grid.draw(g)
dev.off()

### linear (raw) plot
# g <- rbind(g1, g3, size = "first")
# g$widths <- unit.pmax(g1$widths, g3$widths)
# 
# png(paste(savepath,'_vgenome_linear.png',sep=''), 12, 3, units='in', res=300)
# grid.newpage()
# grid.draw(g)
# dev.off()
# 
# ### put toghether
# g <- rbind(g1, g3, g2, size = "first")
# g$widths <- unit.pmax(g1$widths, g3$widths, g2$widths)
# 
# png(paste(savepath,'_vgenome_tgh.png',sep=''), 12, 5, units='in', res=300)
# grid.newpage()
# grid.draw(g)
# dev.off()

}
