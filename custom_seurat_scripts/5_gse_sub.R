## Visualizations
## this script reads all the DE gene results files in a given directory,
## performs GSE the results

library("biomaRt")
library("dplyr")
library("org.Mm.eg.db")
library("GO.db")
library("hypeR")
library("ggplot2")
library("tibble")

## you may want to try other GSE packages that are explicitly used for GSEA like:
## http://www.bioconductor.org/packages/release/bioc/vignettes/escape/inst/doc/vignette.html

## path to directory containing Seurat DE .txt files
# path <- "/Users/u0112671/Documents/tmp/pancreatic_melanie/results/complete3/DE_results/"
# path <- paste0("/Users/u0112671/Documents/tmp/pancreatic_melanie/results/raw1000dr/cancer_treatment_markers/marker_genes_cancer_clusters/")
path <- paste0("/Users/u0112671/Documents/tmp/pancreatic_melanie/gsea/myeloid/myeloid_de_genes/")
# # proj_name <- "complete"
setwd(path)
file_list <- list.files(path, pattern = "*.txt")
### need to add a function that removes non-txt files
rank_list <- list()
data_list <- lapply(file_list, function(i) read.table(file = i, header = TRUE, sep = "\t"))
# names(data_list) <- gsub("_complete.txt", "", file_list) 
# names(data_list) <- gsub("_pancreatic_mice.txt", "", file_list)
names(data_list) <- gsub("_rerun.txt", "", file_list)

## query biomart for gene ids
getBM(attributes = c("ensembl_gene_id", "external_gene_name", "entrezgene_id"),
      mart = useDataset("mmusculus_gene_ensembl", useMart("ensembl"))) %>% 
  dplyr::rename(ensembl_id = ensembl_gene_id,
                gene_name = external_gene_name,
                entrez_id = entrezgene_id) -> mouse_ensembl_2_geneName

## Load Pathway Gene Sets. Gene sets that are mapped to gene ontology (GO) and
## KEGG pathway terms are loaded from the `msigdb` database.
# hypeR::msigdb_info()
GO_BP <- hypeR::msigdb_gsets(species = "Mus musculus", category = "C5",
                             subcategory = "BP")
GO_CC <- hypeR::msigdb_gsets(species = "Mus musculus", category = "C5",
                             subcategory = "CC")
GO_MF <- hypeR::msigdb_gsets(species = "Mus musculus", category = "C5",
                             subcategory = "MF")
KEGG <- hypeR::msigdb_gsets(species = "Mus musculus", category = "C2",
                            subcategory = "CP:KEGG")

## for loop sets max p_val_adj to 1.0e-300 and creates threshold factor
for (i in 1:length(data_list)){
  data_list[[i]]$p_val_adj[data_list[[i]]$p_val_adj == 0] <- 1.0e-300
  # data_list[[i]]$threshold <- as.factor(abs(data_list[[i]]$avg_log2FC) > 1 & 
  #                                           data_list[[i]]$p_val_adj < 0.05)
  data_list[[i]] %>%
    dplyr::rename(gene_name = X) %>%
    dplyr::left_join(mouse_ensembl_2_geneName, by = "gene_name") %>%
    dplyr::filter(!is.na(p_val_adj)) %>%
    dplyr::filter(!is.na(entrez_id)) %>%
    dplyr::filter(!duplicated(entrez_id)) %>%
    dplyr::filter(!is.na(gene_name)) %>%
    dplyr::filter(!duplicated(gene_name)) %>% 
    # dplyr::select(gene_name, avg_logFC) %>%
    # dplyr::arrange(desc(abs(avg_logFC))) %>%
    dplyr::select(gene_name, avg_log2FC) %>%
    dplyr::arrange(desc(abs(avg_log2FC))) %>%
    tibble::deframe() -> rank_list[[i]]
}

## GSEA GO BP 
GO_BP_GSEA <- list()
for (i in 1:length(rank_list)){
fgsea::fgseaMultilevel(pathways = GO_BP$genesets, stats = rank_list[[i]],
                       minSize = 1, maxSize = Inf, eps = 0, nPermSimple = 1000
                       ) -> GO_BP_GSEA[[i]]
}

## GSEA GO CC
GO_CC_GSEA <- list()
for (i in 1:length(rank_list)){
fgsea::fgseaMultilevel(pathways = GO_CC$genesets, stats = rank_list[[i]],
                       minSize = 1, maxSize = Inf, eps = 0, nPermSimple = 1000
                       ) -> GO_CC_GSEA[[i]]
}

## GSEA GO MF
GO_MF_GSEA <- list()
for (i in 1:length(rank_list)){
fgsea::fgseaMultilevel(pathways = GO_MF$genesets, stats = rank_list[[i]],
                       minSize = 1, maxSize = Inf, eps = 0, nPermSimple = 1000
                       ) -> GO_MF_GSEA[[i]]
}

## GSEA KEGG
KEGG_GSEA <- list()
for (i in 1:length(rank_list)){
fgsea::fgseaMultilevel(pathways = KEGG$genesets, stats = rank_list[[i]],
                       minSize = 1, maxSize = Inf, eps = 0, nPermSimple = 1000
                       ) -> KEGG_GSEA[[i]]
}

## save results
for (i in 1:length(rank_list)){
GO_BP_GSEA[[i]] %>% as.data.frame() %>%
                dplyr::mutate(leadingEdge = paste0(leadingEdge)) %>%
write.table(file = paste0(path, names(data_list)[i], "_GO_BP_genes.csv"),
            sep = ";", quote = FALSE, row.names = FALSE)
GO_BP_GSEA[[i]] %>% as.data.frame() %>%
  dplyr::mutate(leadingEdge = paste0(leadingEdge)) %>%
  dplyr::select(-leadingEdge) %>%
write.table(file = paste0(path, names(data_list)[i], "_GO_BP.csv"),
            sep = ";", quote = FALSE, row.names = FALSE)

GO_CC_GSEA[[i]] %>% as.data.frame() %>%
                dplyr::mutate(leadingEdge = paste0(leadingEdge)) %>%
write.table(file = paste0(path, names(data_list)[i], "_GO_CC_genes.csv"),
            sep = ";", quote = FALSE, row.names = FALSE)
GO_CC_GSEA[[i]] %>% as.data.frame() %>%
  dplyr::mutate(leadingEdge = paste0(leadingEdge)) %>%
  dplyr::select(-leadingEdge) %>%
write.table(file = paste0(path, names(data_list)[i], "_GO_CC.csv"),
            sep = ";", quote = FALSE, row.names = FALSE)

GO_MF_GSEA[[i]] %>% as.data.frame() %>%
                dplyr::mutate(leadingEdge = paste0(leadingEdge)) %>%
write.table(file = paste0(path, names(data_list)[i], "_GO_MF_genes.csv"),
            sep = ";", quote = FALSE, row.names = FALSE)
GO_MF_GSEA[[i]] %>% as.data.frame() %>%
  dplyr::mutate(leadingEdge = paste0(leadingEdge)) %>%
  dplyr::select(-leadingEdge) %>%
write.table(file = paste0(path, names(data_list)[i], "_GO_MF.csv"),
            sep = ";", quote = FALSE, row.names = FALSE)

KEGG_GSEA[[i]] %>% as.data.frame() %>%
    dplyr::mutate(leadingEdge = paste0(leadingEdge)) %>%
write.table(file = paste(path, names(data_list)[i], "_KEGG_genes.csv"),
            sep = ";", quote = FALSE, row.names = FALSE)
KEGG_GSEA[[i]] %>% as.data.frame() %>%
  dplyr::mutate(leadingEdge = paste0(leadingEdge)) %>%
  dplyr::select(-leadingEdge) %>%
write.table(file = paste0(path, names(data_list)[i], "_KEGG.csv"),
            sep = ";", quote = FALSE, row.names = FALSE)
}

## bar charts
GO_BP_plot_list <- list()
for (i in 1:length(rank_list)){
PPInfer::GSEA.barplot(GO_BP_GSEA[[i]], category = "pathway", score = "NES",
                      top = 30, pvalue = "padj", sort = "padj", numChar = 100) +
  geom_abline(slope = 0, intercept = -log10(0.05), col = "red", linetype = "dashed") +
  theme(plot.margin = margin(10, 10, 10, 70)) +
  ggtitle(paste0(gsub("_", " ", names(data_list)[i]))) -> GO_BP_plot_list[[i]]
}

GO_CC_plot_list <- list()
for (i in 1:length(rank_list)){
PPInfer::GSEA.barplot(GO_CC_GSEA[[i]], category = "pathway", score = "NES",
                      top = 30, pvalue = "padj", sort = "padj", numChar = 100) +
  geom_abline(slope = 0, intercept = -log10(0.05), col = "red", linetype = "dashed") +
  theme(plot.margin = margin(10, 10, 10, 70)) +
  ggtitle(paste0(gsub("_", " ", names(data_list)[i]))) -> GO_CC_plot_list[[i]]
}

## bar charts
GO_MF_plot_list <- list()
for (i in 1:length(rank_list)){
PPInfer::GSEA.barplot(GO_MF_GSEA[[i]], category = "pathway", score = "NES",
                      top = 30, pvalue = "padj", sort = "padj", numChar = 100) +
  geom_abline(slope = 0, intercept = -log10(0.05), col = "red", linetype = "dashed") +
  theme(plot.margin = margin(10, 10, 10, 70)) +
  ggtitle(paste0(gsub("_", " ", names(data_list)[i]))) -> GO_MF_plot_list[[i]]
}

## bar charts
KEGG_plot_list <- list()
for (i in 1:length(rank_list)){
PPInfer::GSEA.barplot(KEGG_GSEA[[i]], category = "pathway", score = "NES",
                      top = 30, pvalue = "padj", sort = "padj", numChar = 100) +
  geom_abline(slope = 0, intercept = -log10(0.05), col = "red", linetype = "dashed") +
  theme(plot.margin = margin(10, 10, 10, 70)) +
  ggtitle(paste0(gsub("_", " ", names(data_list)[i]))) -> KEGG_plot_list[[i]]
}

## for loop saves a plot pdf for each file in the file list
# for (i in 1:length(GO_BP_plot_list)){
# pdf(paste0("GO_BP_GSEA_", names(data_list)[i], ".pdf"), width = 16, height = 8)
# print(GO_BP_plot_list[[i]])
# dev.off()
# }
# 
# for (i in 1:length(GO_CC_plot_list)){
# pdf(paste0("GO_CC_GSEA_", names(data_list)[i], ".pdf"), width = 16, height = 8)
# print(GO_CC_plot_list[[i]])
# dev.off()
# }
# for (i in 1:length(GO_MF_plot_list)){
# pdf(paste0("GO_MF_GSEA_", names(data_list)[i], ".pdf"), width = 16, height = 8)
# print(GO_MF_plot_list[[i]])
# dev.off()
# }
# for (i in 1:length(KEGG_plot_list)){
# pdf(paste0("KEGG_GSEA_", names(data_list)[i], ".pdf"), width = 16, height = 8)
# print(KEGG_plot_list[[i]])
# dev.off()
# }

##
pdf(paste0("GO_BP_GSEA_all.pdf"), width = 16, height = 8)
for (i in 1:length(GO_BP_plot_list)){
print(GO_BP_plot_list[[i]])
}
dev.off()

pdf(paste0("GO_CC_GSEA_all.pdf"), width = 16, height = 8)
for (i in 1:length(GO_CC_plot_list)){
print(GO_CC_plot_list[[i]])
}
dev.off()

pdf(paste0("GO_MF_GSEA_all.pdf"), width = 16, height = 8)
for (i in 1:length(GO_MF_plot_list)){
print(GO_MF_plot_list[[i]])
}
dev.off()

pdf(paste0("KEGG_GSEA_all.pdf"), width = 16, height = 8)
for (i in 1:length(KEGG_plot_list)){
print(KEGG_plot_list[[i]])
}
dev.off()
