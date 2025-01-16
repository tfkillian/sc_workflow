## the script takes a list of DE results from the Seurat FindAllMarkers function
## and converts it to an Excel spreadsheet, with each cluster on a new sheet

## load libraries
library("dplyr")
library("ggplot2")
library("writexl")

## read seurat object
file_path <- "~/Documents/tmp/202202_mazzone/data/"
# proj_name <- "gcbc"
proj_name <- "expansion"
read.table(file = paste0(file_path, "seu_markers_", proj_name, ".txt"),
# read.table(paste0(file_path, "unfiltered_endo/seu_markers_", proj_name, ".txt"),
           header = TRUE, sep = "\t") %>%
           dplyr::select(-X) %>%
           dplyr::filter(p_val_adj < 0.05) %>%
           dplyr::select(gene, everything()) -> seu_markers

seu_markers$p_val[seu_markers$p_val == 0 ] <- 1e-300
seu_markers$p_val_adj[seu_markers$p_val_adj == 0 ] <- 1e-300

seu_markers %>%
  dplyr::mutate(cluster = dplyr::case_when(
    cluster == "n/a" ~ 0,
    cluster == "E" ~ 1,
    cluster == "NE" ~ 2)) -> seu_markers

# seu_markers$new_label <- seu_markers$cluster
# seu_markers %>% 
#   dplyr::mutate(cluster = dplyr::case_when(
#     cluster == "B20_Responding" ~ 0,
#     cluster == "B20_Relapsing" ~ 1,
#     cluster == "B20_+_Pi3Ki_Responding" ~ 2,
#     cluster == "Untreated_Untreated" ~ 3,        
#     cluster == "B20_+_aPDL1_Relapsing" ~ 4,
#     cluster == "B20_+_aPD1_+_aCTLA4_Responding" ~ 5,
#     cluster == "B20_+_aPDL1_+_Pi3Ki_Responding" ~ 6,
#     cluster == "B20_+_PDL1_+_LTbR_Responding" ~ 7)) -> seu_markers

# seu_markers$new_label <- seu_markers$cluster
# seu_markers %>%
#   dplyr::mutate(cluster = dplyr::case_when(
#     cluster == "IT-like" ~ 0,
#     cluster == "IT-like B20 PiK3i" ~ 1,
#     cluster == "Unknown B20 PiK3i Rel" ~ 2,
#     cluster == "Ins-hi MLP-like" ~ 3,
#     cluster == "Unknown B20 Resp" ~ 4,
#     cluster == "IT-like B20 Resp" ~ 5,
#     cluster == "Ins-lo MLP-like B20" ~ 6,
#     cluster == "IT-like B20 aPDL1 CTLA4 Resp" ~ 7,
#     cluster == "Ins-lo MLP-like B20 Resp" ~ 8,
#     cluster == "Ins-lo MLP-like" ~ 9,
#     cluster == "Acinar" ~ 10,
#     cluster == "Ductal/Stem" ~ 11,
#     cluster == "Ins-hi MLP-like B20 PiK3i Rel" ~ 12,
#     cluster == "MLP-like CCB104" ~ 13
#     )) -> seu_markers

# seu_markers$new_label <- seu_markers$cluster
# seu_markers %>%
#   dplyr::mutate(cluster = dplyr::case_when(
#     cluster == "IT-like" ~ 0,
#     cluster == "IT-like Ins-lo B20 Responding" ~ 1,
#     cluster == "IT-like B20 aPDL1 CTLA4 Responding" ~ 2,
#     cluster == "IT-like Ins-lo B20 + PiK3i Relapsing" ~ 3,
#     cluster == "IT-like Ins-lo B20" ~ 4,
#     cluster == "MLP-like Ins-lo" ~ 5
#     )) -> seu_markers

# seu_markers$cluster <- seu_markers$new_label 
# seu_markers %>%
#   dplyr::mutate(cluster = dplyr::case_when(
#     cluster == "M1 Il1a+" ~ 0,
#     cluster == "M2 Cd63+" ~ 1,
#     cluster == "MAST" ~ 2,
#     cluster == "M2 Folr2+" ~ 3,
#     cluster == "Monocyte Chil3+" ~ 4,
#     cluster == "M2 Malat1-" ~ 5,
#     cluster == "M1 Cd72+" ~ 6,
#     cluster == "cDCs" ~ 7,
#     cluster == "M2 Birc5+" ~ 8,
#     cluster == "M2 Sepp1+" ~ 9,
#     cluster == "pDCs" ~ 10,
#     cluster == "Unknown pDCs" ~ 11,
#     cluster == "Monocyte Ifit2bl1+" ~ 12,
#     cluster == "Neutrophils" ~ 13,
#     cluster == "Proliferating pDCs" ~ 14
#     )) -> seu_markers

## add column in metadata with the "Cluster" 

# seu_markers %>%
#   dplyr::mutate(cluster = dplyr::case_when(
#     cluster == "Cd8+ Exhausted" ~ 0,
#     cluster == "Unknown T-cells Th17+" ~ 1,
#     cluster == "Tregs" ~ 2,
#     cluster == "Unknown T-cells" ~ 3,
#     cluster == "Cd4+ Th2+" ~ 4,
#     cluster == "Cd8+ Memory" ~ 5,
#     cluster == "NK cells" ~ 6,
#     cluster == "Proliferating T-cells" ~ 7,
#     cluster == "Gamma Delta T-cells" ~ 8
#     )) -> seu_markers

# seu_markers %>%
#   dplyr::mutate(cluster = dplyr::case_when(
#     cluster == "M2 Macrophages" ~ 0,
#     cluster == "M1 Macrophages" ~ 1,
#     cluster == "Proliferating  M2 Macrophages" ~ 2,
#     cluster == "Monocytes" ~ 3,
#     cluster == "Mast cells" ~ 4,
#     cluster == "Unknown M2 Macrophages" ~ 5,
#     cluster == "pDCs" ~ 6,
#     cluster == "cDCs" ~ 7,
#     cluster == "Neutrophils" ~ 8
#     )) -> seu_markers

# seu_markers %>%
#   dplyr::mutate(cluster = dplyr::case_when(
#     cluster == "M2 Macrophages" ~ 0,
#     cluster == "M1 Macrophages" ~ 1,
#     cluster == "Unknown Myeloids" ~ 2,
#     cluster == "Prol. M2 Macrophages" ~ 3,
#     cluster == "Mast cells/pDCs" ~ 4,
#     cluster == "F13a1+ M2 Macrophages" ~ 5,
#     cluster == "cDCs" ~ 6,
#     cluster == "Monocytes" ~ 7,
#     cluster == "Neutrophils" ~ 8
#     )) -> seu_markers

marker_list <- list()
for (i in 1:length(unique(seu_markers$cluster))) {
  j <- (i - 1)
  seu_markers %>% 
    dplyr::filter(cluster == j) -> marker_list[[i]]
}

name_vec <- vector()
for (i in 1:length(unique(seu_markers$cluster))) {
  name_vec[i] <- paste0("cluster_", i - 1)
}

names(marker_list) <- name_vec

## save file
write_xlsx(path = paste0(proj_name, "_cluster_sig_markers.xlsx"), marker_list)
