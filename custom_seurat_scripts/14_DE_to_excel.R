## the script takes a list of DE results from the Seurat FindAllMarkers function
## and converts it to an Excel spreadsheet, with each cluster on a new sheet

## load libraries
library("dplyr")
library("ggplot2")
library("writexl")
library("biomaRt")

## read seurat object
# file_path <- "~/Documents/tmp/pancreatic_melanie/"
cell_type <- "tcell"
# file_path <- paste0("/Users/u0112671/Documents/tmp/pancreatic_melanie/gsea/",
#                     cell_type, "/", cell_type, "_de_genes/")
file_path <- paste0("/Users/u0112671/Documents/tmp/pancreatic_melanie/",
                    cell_type, "/")
proj_name <- "rerun"
setwd(file_path)
file_list <- list.files(file_path, pattern = "*.txt")
seu_markers <- list()
data_list <- lapply(file_list, function(i) read.table(file = i, header = TRUE, sep = "\t"))
names(data_list) <- gsub("_rerun.txt", "", file_list)

## query biomart for gene ids
getBM(attributes = c("ensembl_gene_id", "external_gene_name", "entrezgene_id"),
      mart = useDataset("mmusculus_gene_ensembl", useMart("ensembl"))) %>% 
  dplyr::rename(ensembl_id = ensembl_gene_id,
                gene_name = external_gene_name,
                entrez_id = entrezgene_id) -> mouse_ensembl_2_geneName

## for loop sets max p_val_adj to 1.0e-300 and creates threshold factor
for (i in 1:length(data_list)){
  data_list[[i]]$p_val_adj[data_list[[i]]$p_val_adj == 0] <- 1.0e-300
  data_list[[i]]$p_val[data_list[[i]]$p_val == 0] <- 1.0e-300
  data_list[[i]]$threshold <- as.factor(abs(data_list[[i]]$avg_log2FC) > 1 &
                                            data_list[[i]]$p_val_adj < 0.05)
  data_list[[i]] %>%
    dplyr::rename(gene_name = X) %>%
    dplyr::left_join(mouse_ensembl_2_geneName, by = "gene_name") %>%
    dplyr::filter(!is.na(p_val_adj)) %>% 
    dplyr::filter(!duplicated(entrez_id)) %>%
    dplyr::filter(!duplicated(gene_name)) %>% 
    dplyr::filter(p_val_adj < 0.05) -> seu_markers[[i]]
}

## shorten comparison names to fit on Excel spreadsheet tabs
name_vec <- gsub("DE_", "", names(data_list))
# name_vec <- gsub("_resp", "_rs", name_vec)
# name_vec <- gsub("_rel", "_rl", name_vec)
name_vec <- gsub("_vs", "_v", name_vec)
name_vec <- gsub("_in_macro", "", name_vec)
names(seu_markers) <- name_vec

## save file
write_xlsx(path = paste0(cell_type, "_", proj_name, "_DE_sig_markers.xlsx"), seu_markers)
