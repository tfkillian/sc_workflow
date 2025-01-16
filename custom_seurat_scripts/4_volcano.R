## Visualizations
## this script reads all the DE gene results files in a given directory,
## performs QC on the results and creates a volcano plot for each file

library("dplyr")
library("ggplot2")

## path to directory containing Seurat DE .txt files
# path <- "/Users/u0112671/Documents/tmp/pancreatic_melanie/results/complete3/DE_results/"
path <- paste0("/Users/u0112671/Documents/tmp//")
# proj_name <- "complete"
setwd(path)
fl <- list.files(path)
### need to add a function that removes non-txt files
# dl <- list()
dl <- lapply(fl, function(i) read.table(file = i, header = TRUE, sep = "\t"))
# names(dl) <- gsub("_complete.txt", "", fl)
names(dl) <- gsub("_rerun.txt", "", fl)

## for loop sets max p_val_adj to 1.0e-300 and creates threshold factor
for (i in 1:length(dl)){
  dl[[i]]$p_val_adj[dl[[i]]$p_val_adj == 0] <- 1.0e-300
  dl[[i]]$threshold <- as.factor(abs(dl[[i]]$avg_log2FC) > 1 & 
                                            dl[[i]]$p_val_adj < 0.05)
  dl[[i]] <- dl[[i]] %>% dplyr::rename(gene = X)
}

## for loop creates a volcano plot for each file in the file list
gglist <- list()
for (i in 1:length(dl)){
#pdf(paste0("Volcano_plots", gsub("_", " ", names(dl)[i]), ".pdf"))
ggplot2::ggplot(
  data = dl[[i]],
  aes(x = avg_log2FC, y = -log10(p_val_adj), color = threshold)) +
  geom_text(label = dl[[i]]$gene, nudge_x = 0.25, nudge_y = 0.25,
            check_overlap = TRUE) +
  theme(legend.position = "none") +
  geom_point(alpha = 0.4, size = 0.5) +
  xlab("log2 fold change") + ylab("-log10 padj-value") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("#000000", "#FF0000")) +
  ggtitle(paste0(gsub("_", " ", names(dl)[i]))) -> gglist[[i]]
#dev.off()
}

## for loop saves a volcano plot pdf for each file in the file list
pdf(paste0("Vol_plot_all.pdf"))
for (i in 1:length(dl)){
print(gglist[[i]])
}
dev.off()

# for (i in 1:length(dl)){
# pdf(paste0("Vol_plot_", names(dl)[i], ".pdf"))
# print(gglist[[i]])
# dev.off()
# }
