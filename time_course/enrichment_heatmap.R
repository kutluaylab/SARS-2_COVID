# Gene set over-representation analysis for gene clusters from time course heatmap


library(dplyr)
library(msigdbr)
library(clusterProfiler)
library(tidyr)
library(stringr)
library(circlize)
library(ComplexHeatmap)



### SETUP ###

# The base path should contain an enrichment_heatmap/ folder and a time_course_heatmap/ folder, which was created in the time course heatmap script.
# The enrichment_heatmap/ folder should contain the following subfolders:
# enrichment_results/, heatmap/, and selected_term/.

basepath <- "/path/tp/basepath/"

path <- file.path(basepath, "enrichment_heatmap")
cpm_path <- file.path(basepath, "time_course_heatmap", "sig_zscore")

prefix <- "experiment_name" # ex. HBEC_ribo



### DATA PREPARATION ###

# filtered cpm table from time course cpm heatmap for gene universe
cpm <- read.delim(file.path(cpm_path, paste0(prefix, "_TCHeatmap_cpm")),
                  row.names = NULL)

outputname <- paste0(prefix, "_ER")

# cpm table of DE genes for cluster information
sig_cpm <- read.delim(file.path(cpm_path, paste0(prefix, "_TCHeatmap_sig_cpm_cluster")))



gene_cluster <- na.omit(data.frame(row.names(sig_cpm), sig_cpm$cluster))
colnames(gene_cluster) <- c("gene_name", "cluster")

# total number of clusters
num_cluster <- max(gene_cluster$cluster)

# clusters selected to be plotted in the enrichment heatmap, will modify the variable later
selected_clusters <- seq(num_cluster)


cpm <- rename(cpm, gene_name = row.names)


# deduplication
# getting rid of duplicate genes in cpm table, including the original one
duplicate_genes <- cpm$gene_name[duplicated(cpm$gene_name)]
cpm <- cpm[!(cpm$gene_name %in% duplicate_genes), ]
print(paste0(length(duplicate_genes) + length(unique(duplicate_genes)),
             " duplicate genes removed."))



geneNames.universe <- cpm$gene_name


# getting rid of genes in the duplicate list that may be present in clustered genes
gene_cluster <- gene_cluster[!gene_cluster$gene_name %in% duplicate_genes, ]
stopifnot(gene_cluster$gene_name %in% geneNames.universe)


### just for Vero riboseq ###
selected_clusters <- c(1, 2, 4, 5, 6)
gene_cluster_orig <- gene_cluster
gene_cluster <- gene_cluster[gene_cluster$cluster %in% selected_clusters, ]
num_cluster <- length(selected_clusters)
#############################


## Retrieving info from Molecular Signatures database ##
msig_bp <- msigdbr(species = "Homo sapiens") %>%
  filter(gs_cat == "H" | (gs_cat == "C5" & gs_subcat == "GO:BP")) %>%
  select(gs_name, human_gene_symbol)

save(msig_bp, file = file.path(path, "enrichment_results", "ER_heatmap_msig.RData"))
load(file = file.path(path, "enrichment_results", "ER_heatmap_msig.RData"))




### OVER-REPRESENTATION ANALYSIS ###

enrich <- function(cluster_table, universe, database, selected_clusters) {
  enrichment.results <- NULL

  ## Over-representation analysis for each cluster, filtering q-value <= 0.1 ##
  for (cluster.id in selected_clusters) {

    geneNames.cluster <- cluster_table[cluster_table$cluster == cluster.id, ]$gene_name

    em <- enricher(gene = geneNames.cluster,
                   pAdjustMethod = "fdr",
                   universe = universe,
                   TERM2GENE = database,
                   minGSSize = 10,
                   qvalueCutoff = 0.1)


    if (!is.null(em)) {
      result <- em@result
      result$cluster <- cluster.id
      result <- filter(result, qvalue <= 0.1)
      if (nrow(result) > 0) {
        enrichment.results <- rbind(enrichment.results, result)
      }
    }
  }

  return(enrichment.results)
}


enrichment.results <- enrich(cluster_table = gene_cluster,
                             universe = geneNames.universe,
                             database = msig_bp,
                             selected_clusters = selected_clusters)



enrichment.results$cluster <- as.integer(enrichment.results$cluster)

# clusters that contain terms with q-value <= 0.1
enriched_cluster <- unique(enrichment.results$cluster)

enrichment.results <- rename(enrichment.results, genes = geneID)

enrichment.results$category <- ""
enrichment.results$category[str_detect(enrichment.results$Description, pattern = "^HALLMARK")] <- "Hallmark"
enrichment.results$category[str_detect(enrichment.results$Description, pattern = "^GOBP")] <- "Biological Process"



## Getting number of genes present in each cluster ##
num.in.cluster <- data.frame(table(gene_cluster$cluster))
colnames(num.in.cluster) <- c("cluster", "n")
write.csv(num.in.cluster,
          file.path(path, "enrichment_results", paste0(outputname, "_num_in_cluster.csv")),
          row.names = F)


## Getting the ratio of genes in a cluster that belong to a term against the total number of genes in the cluster ##
enrichment.results$ratio <- enrichment.results$Count / num.in.cluster[match(enrichment.results$cluster, num.in.cluster$cluster), ]$n
write.csv(enrichment.results,
          file.path(path, "enrichment_results", paste0(outputname, "_results.csv")),
          row.names = F)


## "Flattening" the results by showing results for each cluster in parallel, merging terms shared by different clusters ##
enrichment.results_wide <- pivot_wider(enrichment.results,
                                       id_cols = c(Description, category),
                                       names_from = cluster,
                                       values_from = c(pvalue, p.adjust, qvalue, Count,
                                                       GeneRatio, BgRatio, ratio, genes))
enrichment.results_wide <- arrange(enrichment.results_wide, category)
write.csv(enrichment.results_wide,
          file.path(path, "enrichment_results", paste0(outputname, "_results_wide.csv")),
          row.names = F)
save(enrichment.results_wide,
     file = file.path(path, "enrichment_results", paste0(outputname, "_results_wide.RData")))




## Manual selection of interesting terms ##

# The terms shown on the eventual heatmap are manually selected from the terms in enrichment_results_wide.
# The indices of the selected terms should be specified in a "selected_term_index.txt" file in the "selected_term" folder.
# Each index should take up one line in the file.

load(file.path(path, "enrichment_results", paste0(outputname, "_results_wide.RData")))


selected_term_index <- read.table(file = file.path(path, "selected_term", "selected_term_index.txt"),
                                  header = F)[, 1]

enrichment.results.pruned <- enrichment.results_wide[selected_term_index, ]

enrichment.results.pruned <- arrange(enrichment.results.pruned, category)


# clusters that have selected terms
# Note that it is possible that not all clusters with significantly enriched terms have a term that is selected to be in the final heatmap.
# In other words, enriched_cluster and selected_clusters may not be the same.
selected_clusters <- c(2, 5)


## Changing the style of the Description column ##
# convert to lower case; replace "_" with " "; get rid of the first word (for example, "gobp").
selected_terms <- tolower(enrichment.results.pruned$Description) %>%
  gsub(pattern = "_", replacement = " ") %>%
  substr(start = nchar(word(., start = 1)) + 2, stop = nchar(.))

write.table(data.frame(selected_terms),
            file = file.path(path, "selected_term", paste0(outputname, "_selected_terms.txt")),
            col.names = F, row.names = F, quote = F)

# Here, some manual style editing of the selected terms is required.
# For example, change the first letter of each word to upper case, and change all abbreviations to upper case.
# Save the edited selected term file to "[outputname]_selected_terms_edited.txt".

## After manual editing of descriptions ##
selected_terms_edited <- read.table(file = file.path(path, "selected_term", paste0(outputname, "_selected_terms_edited.txt")),
                                    header = F,
                                    sep = "\t")[, 1]

enrichment.results.pruned$Description <- selected_terms_edited


## Labeling hallmark terms as a special case of biological process terms and reordering the terms (only used for HBEC)
# Hallmark term indices may need to be modified to fit your run.

## Only for HBEC rseq ##
enrichment.results.pruned[26, "Description"] <- "Coagulation (HM)"
enrichment.results.pruned[26, "category"] <- "Biological Process"

custom.order <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 26, 18, 19, 20, 21, 22, 23, 24, 25)
enrichment.results.pruned <- enrichment.results.pruned[custom.order, ]
########################

## Only for HBEC ribo ##
# labeling hallmark terms as a special case of biological process terms and reordering the terms
enrichment.results.pruned[22, "Description"] <- "Coagulation (HM)"
enrichment.results.pruned[22, "category"] <- "Biological Process"
enrichment.results.pruned[23, "Description"] <- "KRAS Signaling Up (HM)"
enrichment.results.pruned[23, "category"] <- "Biological Process"

custom.order <- c(1, 2, 3, 4, 5, 6, 7, 8, 22, 23, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21)
enrichment.results.pruned <- enrichment.results.pruned[custom.order, ]
########################


save(enrichment.results.pruned,
     file = file.path(path, "enrichment_results", paste0(outputname, "_pruned.RData")))
load(file = file.path(path, "enrichment_results", paste0(outputname, "_pruned.RData")))

write.csv(enrichment.results.pruned,
          file = file.path(path, "enrichment_results", paste0(outputname, "_pruned.csv")),
          row.names = F)



## Selecting q-values for plotting ##
em.q.values <- select(enrichment.results.pruned, grep(pattern = "qvalue",
                                                      colnames(enrichment.results.pruned),
                                                      value = T))
enrichment.results.pruned$Description[20] <- "Inflammatory Response " # Vero rseq only

# getting rid of q-value columns whose cluster does not contain selected terms
em.q.values <- em.q.values[, selected_clusters %in% enriched_cluster]

em.q.values <- as.data.frame(em.q.values)
row.names(em.q.values) <- as.character(enrichment.results.pruned$Description)
colnames(em.q.values) <- selected_clusters


# taking log10 of the q-values
em.q.values <- -log10(as.matrix(em.q.values))


write.csv(em.q.values,
          file = file.path(path, "enrichment_results", paste0(outputname, "_qvalues.csv")),
          row.names = F)

save(em.q.values,
     file = file.path(path, "enrichment_results", paste0(outputname, "_qvalues.RData")))
load(file = file.path(path, "enrichment_results", paste0(outputname, "_qvalues.RData")))




## Selecting ratio values for plotting ##
ratio.values <- select(enrichment.results.pruned, grep(pattern = "ratio",
                                                       colnames(enrichment.results.pruned),
                                                       value = T))

ratio.values <- ratio.values[, selected_clusters %in% enriched_cluster]

ratio.values <- as.data.frame(ratio.values)
row.names(ratio.values) <- as.character(enrichment.results.pruned$Description)
colnames(ratio.values) <- selected_clusters


write.csv(ratio.values,
          file = file.path(path, "enrichment_results", paste0(outputname, "_ratio.csv")),
          row.names = F)
save(ratio.values,
     file = file.path(path, "enrichment_results", paste0(outputname, "_ratio.RData")))
load(file = file.path(path, "enrichment_results", paste0(outputname, "_ratio.RData")))






### ENRICHMENT HEATMAP ###

## col_fun function for heatmap to plot different shades of colors of dots according to q-values ##
max.q <- max(em.q.values, na.rm = T)
min.q <- min(em.q.values, na.rm = T)

col_fun <- colorRamp2(c(min.q, max.q), c("black", "red"))


## Column annotation blocks with cluster information ##
colPalette <- c("#AEC7E87F", "#98DF8A7F", "#FF98967F",
                "#C49C947F", "#C5B0D57F", "#FFBB787F",
                "turquoise", "pink", "maroon", "wheat3")


col.annot <- HeatmapAnnotation(cluster = anno_simple(selected_clusters,
                                                     height = unit(3, "mm"),
                                                     col = structure(colPalette[selected_clusters],
                                                                     names = selected_clusters)
                                                     ),
                               show_annotation_name = F,
                               show_legend = F)



## Heatmap ##

# Graphic parameters can be customized for the atheistics of the plot.

h <- Heatmap(em.q.values,
             # labels
             column_title = "Enriched Gene Sets",
             column_title_gp = gpar(fontsize = 9),
             column_names_gp = gpar(fontsize = 10),
             column_names_rot = 0,
             row_title_gp = gpar(fontsize = 12),
             row_names_gp = gpar(fontsize = 9),

             # annotation
             bottom_annotation = col.annot,

             # legends (custom legends are introduced later)
             show_heatmap_legend = F,

             show_row_names = T,
             show_column_names = T,

             # clustering
             cluster_rows = F,
             cluster_columns = F,

             ## Vero rseq and Vero ribo only (uncomment) ##
             # split = enrichment.results.pruned$category,
             ##############################################

             # circles
             cell_fun = function(j, i, x, y, width, height, fill) {
               grid.rect(x = x, y = y, width = width, height = height,
                         gp = gpar(fill = "white", col = "#EEEEEE")) # white rectangles in the background
               grid.circle(x = x, y = y,
                           r = sqrt(ratio.values[i, j]) * 8, # <<< change size of circles according to ratio distribution, so that all circles fit in the grids
                           default.units = "mm",

                           # coloring the circles according to the q-values
                           gp = gpar(fill = col_fun(em.q.values[i, j]), col = NA))
             }
)



## Q-value legend ##

q_legend <- Legend(col_fun = col_fun,
                   title = "-log10(q-value)",
                   title_gp = gpar(fontsize = 9),
                   labels_gp = gpar(fontsize = 8),
                   grid_height = unit(4, "mm"),
                   direction = "horizontal")



## Ratio value legend ##

# ratio legend breaks: please only uncomment one of the following for the corresponding experiment.
ratio_breaks <- c(0.1, 0.075, 0.05, 0.025) # HBEC rseq and HBEC ribo
# ratio_breaks <- c(0.25, 0.2, 0.15, 0.1, 0.05) # Vero rseq
# ratio_breaks <- c(0.8, 0.6, 0.4, 0.2) # Vero ribo

ratio_legend <- Legend(labels = ratio_breaks,
                       grid_height = unit(6, units = "mm"),
                       grid_width = unit(6, units = "mm"),
                       title = "Percent\nof Cluster",
                       title_gp = gpar(fontsize = 9),
                       labels_gp = gpar(fontsize = 8),
                       type = "points",
                       pch = 1, # circle

                       # circle size of ratio legend breaks
                       # IMPORTANT: change the multiplier number below so that the legend circle sizes match the sizes of circles in the heatmap grids.
                       size = unit(sqrt(ratio_breaks) * 20, units = "mm"),

                       background = "white"
)



## Graphic Devices Settings ##
# Can customize the numbers for aesthetics
# Only uncomment one of the following blocks at a time for the given experiment.

# Vero rseq only
# png(filename = file.path(path, "heatmap", paste0(outputname, ".png")), height = 2000, width = 1350, res = 300)
# draw(h, padding = unit(c(2, 5, 8, 20), "mm")) # padding: empty space around heatmap, for legends, row names...
# draw(q_legend, x = unit(0.7, "npc"), y = unit(0.95, "npc")) # drawing the legend at (x, y) position
# draw(ratio_legend, x = unit(0.9, "npc"), y = unit(0.86, "npc"))
# dev.off()

# Vero ribo only
# png(filename = file.path(path, "heatmap", paste0(outputname, ".png")), height = 2000, width = 1550, res = 300)
# draw(h, padding = unit(c(2, 5, 8, 15), "mm")) # padding: empty space around heatmap, for legends, row names...
# draw(q_legend, x = unit(0.7, "npc"), y = unit(0.95, "npc")) # drawing the legend at (x, y) position
# draw(ratio_legend, x = unit(0.9, "npc"), y = unit(0.88, "npc"))
# dev.off()

# HBEC rseq only
# png(filename = file.path(path, "heatmap", paste0(outputname, ".png")), height = 2000, width = 1400, res = 300)
# draw(h, padding = unit(c(2, 5, 8, 20), "mm")) # padding: empty space around heatmap, for legends, row names...
# draw(q_legend, x = unit(0.7, "npc"), y = unit(0.95, "npc")) # drawing the legend at (x, y) position
# draw(ratio_legend, x = unit(0.9, "npc"), y = unit(0.89, "npc"))
# dev.off()

# HBEC ribo only
png(filename = file.path(path, "heatmap", paste0(outputname, ".png")), height = 2000, width = 1350, res = 300)
draw(h, padding = unit(c(2, 5, 8, 20), "mm")) # padding: empty space around heatmap, for legends, row names...
draw(q_legend, x = unit(0.7, "npc"), y = unit(0.95, "npc")) # drawing the legend at (x, y) position
draw(ratio_legend, x = unit(0.9, "npc"), y = unit(0.86, "npc"))
dev.off()

