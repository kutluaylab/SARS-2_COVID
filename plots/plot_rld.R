library(ggplot2)


args <- commandArgs()

# directory to save the output plots
wd <- args[7]
# input distribution information: must be full path if not in wd
file <- args[6]
outputname <- args[8]

setwd(wd)


## Plotting with raw counts

dat <- read.table(file)
colnames(dat) <- c('count','length')

head(dat)


if (nrow(dat) > 31) {
  dat_select <- dat[1:31, ]
} else {
  dat_select <- dat
}  # selecting part of dat to plot

p <- (ggplot(data = dat_select, mapping = aes(x = length, y = count))
  + geom_bar(stat = 'identity', fill = 'black', width = 0.9)
  + scale_x_continuous(breaks = min(dat_select$length):max(dat_select$length))
  + xlab('Read Length')
  + ylab('Count')
  + ggtitle(outputname)
  + theme_classic()
  + theme(plot.title = element_text(hjust = 0.5, size = 12),
          axis.text = element_text(size = 7))
  )

outplotname <- paste0(outputname, '.png')
png(outplotname, width = 1500, height = 900, res = 300)
p
dev.off()


## Calculating the counts percentage per read length (Normalizing)

totsum <- sum(dat[,1])
new <- data.frame(length = dat[,2], count = dat[,1], percentage = dat[,1] / totsum)
head(new)

outtablename <- paste0(outputname, '_percentage')
write.table(new, file = outtablename, sep = "\t", quote = F, row.names = F)
