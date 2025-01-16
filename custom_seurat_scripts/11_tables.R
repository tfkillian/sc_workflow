## this script generates tables of a particular cell type to display in a report

## load libraries
library("dplyr")
file_path <- "/Users/u0112671/Documents/tmp/pancreatic_melanie/"

## these are the samples that are from the whole low filter data
proj_name_0 <- "pancreatic_mice"

## per sample
read.table(paste0(file_path, "cell_per_c_per_O_", proj_name_0, ".txt"),
           header = TRUE, sep = "\t") %>% t() %>% as.data.frame() %>%
  dplyr::slice(-1) %>%
  dplyr::mutate_all(~ as.numeric(.)) -> cell_per_c_0
cell_per_c_0$total <- rowSums(cell_per_c_0)

## per response
read.table(paste0(file_path, "cell_per_c_per_R_", proj_name_0, ".txt"),
           header = TRUE, sep = "\t") %>% t() %>% as.data.frame() %>%
  dplyr::slice(-1) %>%
  dplyr::mutate_all(~ as.numeric(.))-> cell_per_r_0
cell_per_r_0$total <- rowSums(cell_per_r_0)

## per treatment
read.table(paste0(file_path, "cell_per_c_per_T_", proj_name_0, ".txt"),
           header = TRUE, sep = "\t") %>% t() %>% as.data.frame() %>%
  dplyr::slice(-1) %>%
  dplyr::mutate_all(~ as.numeric(.)) -> cell_per_t_0
cell_per_t_0$total <- rowSums(cell_per_t_0)


## these are the samples that are ones taken from combining the MACRO + MONO +
## NEUTRO clusters in the low filter data, then the neutrophil subcluster is
## isolated and has 271 cells
# proj_name_1 <- "neutro1"
# proj_name_1 <- "endo"
# 
# ## per sample
# read.table(paste0(file_path, "cell_per_c_per_O_", proj_name_1, ".txt"),
#            header = TRUE, sep = "\t") %>% t() %>% as.data.frame() -> cell_per_c_1 # %>%
#   dplyr::rename(Macrophages = V1,
#                 Monocytes = V2) %>%
#   dplyr::slice(-1) %>%
#   dplyr::mutate(Macrophages = as.numeric(Macrophages),
#                 Monocytes = as.numeric(Monocytes)
#                 ) -> cell_per_c_1
# cell_per_c_1$cluster_11 <- rowSums(cell_per_c_1)
# 
# ## per response
# read.table(paste0(file_path, "cell_per_c_per_R_", proj_name_1, ".txt"),
#            header = TRUE, sep = "\t") %>% t() %>% as.data.frame() %>%
#   dplyr::rename(Macrophages = V1,
#                 Monocytes = V2) %>%
#   dplyr::slice(-1) %>%
#   dplyr::mutate(Macrophages = as.numeric(Macrophages),
#                 Monocytes = as.numeric(Monocytes)) -> cell_per_r_1
# cell_per_r_1$cluster_11 <- rowSums(cell_per_r_1)
# 
# ## per treatment
# read.table(paste0(file_path, "cell_per_c_per_T_", proj_name_1, ".txt"),
#            header = TRUE, sep = "\t") %>% t() %>% as.data.frame() %>%
#   dplyr::rename(Macrophages = V1,
#                 Monocytes = V2) %>%
#   dplyr::slice(-1) %>%
#   dplyr::mutate(Macrophages = as.numeric(Macrophages),
#                 Monocytes = as.numeric(Monocytes)) -> cell_per_t_1
# cell_per_t_1$cluster_11 <- rowSums(cell_per_t_1)


## these are the samples that are ones taken from subsetting cluster 50 from the
## low filter data, which is neutrophils and has 284 cells
# proj_name_2 <- "neutrophils3"
# 
# ## per sample
# read.table(paste0(file_path, "cell_per_c_per_O_", proj_name_2, ".txt"),
#            header = TRUE, sep = "\t") %>% t() %>% as.data.frame() %>%
#   dplyr::rename(cluster_50 = V1) %>%
#   dplyr::slice(-1) %>%
#   dplyr::mutate(cluster_50 = as.numeric(cluster_50)) -> cell_per_c_2
# 
# ## per response
# read.table(paste0(file_path, "cell_per_c_per_R_", proj_name_2, ".txt"),
#            header = TRUE, sep = "\t") %>% t() %>% as.data.frame() %>%
#   dplyr::rename(cluster_50 = V1) %>%
#   dplyr::slice(-1) %>%
#   dplyr::mutate(cluster_50 = as.numeric(cluster_50)) -> cell_per_r_2
# 
# ## per treatment
# read.table(paste0(file_path, "cell_per_c_per_T_", proj_name_2, ".txt"),
#            header = TRUE, sep = "\t") %>% t() %>% as.data.frame() %>%
#   dplyr::rename(cluster_50 = V1) %>%
#   dplyr::slice(-1) %>%
#   dplyr::mutate(cluster_50 = as.numeric(cluster_50)) -> cell_per_t_2

# proj_name_1 <- "endo"
# proj_name_1 <- "cancer6"
proj_name_1 <- "macro_mono_no_ins2"

## per sample
read.table(paste0(file_path, "cell_per_c_per_O_", proj_name_1, ".txt"),
           header = TRUE, sep = "\t") %>% t() %>% as.data.frame() %>%
  # dplyr::rename(cluster_0 = V1,
  #               cluster_1 = V2,
  #               cluster_2 = V3,
  #               cluster_3 = V4,
  #               cluster_4 = V5,
  #               cluster_5 = V6,
  #               cluster_6 = V7,
  #               cluster_7 = V8,
  #               cluster_8 = V9,
  #               cluster_9 = V10,
  #               cluster_10 = V11,
  #               cluster_11 = V12) %>%
  dplyr::slice(-1) %>%
  dplyr::mutate_all(~ as.numeric(.)) -> cell_per_c_1
cell_per_c_1$endo <- rowSums(cell_per_c_1)

## per response
read.table(paste0(file_path, "cell_per_c_per_R_", proj_name_1, ".txt"),
           header = TRUE, sep = "\t") %>% t() %>% as.data.frame() %>%
  # dplyr::rename(cluster_0 = V1,
  #               cluster_1 = V2,
  #               cluster_2 = V3,
  #               cluster_3 = V4,
  #               cluster_4 = V5,
  #               cluster_5 = V6,
  #               cluster_6 = V7,
  #               cluster_7 = V8,
  #               cluster_8 = V9,
  #               cluster_9 = V10,
  #               cluster_10 = V11,
  #               cluster_11 = V12) %>%
  dplyr::slice(-1) %>%
  dplyr::mutate_all(~ as.numeric(.)) -> cell_per_r_1
cell_per_r_1$endo <- rowSums(cell_per_r_1)

## per treatment
read.table(paste0(file_path, "cell_per_c_per_T_", proj_name_1, ".txt"),
           header = TRUE, sep = "\t") %>% t() %>% as.data.frame() %>%
  # dplyr::rename(cluster_0 = V1,
  #               cluster_1 = V2,
  #               cluster_2 = V3,
  #               cluster_3 = V4,
  #               cluster_4 = V5,
  #               cluster_5 = V6,
  #               cluster_6 = V7,
  #               cluster_7 = V8,
  #               cluster_8 = V9,
  #               cluster_9 = V10,
  #               cluster_10 = V11,
  #               cluster_11 = V12) %>%
  dplyr::slice(-1) %>%
   dplyr::mutate_all(~ as.numeric(.)) -> cell_per_t_1
cell_per_t_1$endo<- rowSums(cell_per_t_1)

### putting it all together
c0 <- sum(cell_per_c_0$total)
# c1 <- sum(cell_per_c_1$cluster_11)
c1 <- sum(cell_per_c_1$endo)
# c2 <- sum(cell_per_c_2$cluster_50)
r0 <- sum(cell_per_r_0$total)
# r1 <- sum(cell_per_r_1$cluster_11)
r1 <- sum(cell_per_r_1$endo)
# r2 <- sum(cell_per_r_2$cluster_50)
t0 <- sum(cell_per_t_0$total)
# t1 <- sum(cell_per_t_1$cluster_11)
t1 <- sum(cell_per_t_1$endo)
# t2 <- sum(cell_per_t_2$cluster_50)

## per sample
cell_per_c_0 %>% 
  tibble::rownames_to_column() %>%
  dplyr::select(rowname, total) %>% 
  dplyr::left_join((cell_per_c_1 %>%
                    tibble::rownames_to_column() %>%
                    # dplyr::select(rowname, cluster_11)), by = "rowname") %>%
                    dplyr::select(rowname, endo)), by = "rowname") %>% 
  # dplyr::left_join((cell_per_c_2 %>%
  #                   tibble::rownames_to_column() %>%
  #                   dplyr::select(rowname, cluster_50)), by = "rowname") %>% 
  dplyr::mutate(total_percent = round((total / c0) * 100, 3),
                # sub1_percent = (cluster_11 / c1) * 100,
                sub1_percent = round((endo / c1) * 100, 3)
                # sub2_percent = (cluster_50 / c2) * 100
                ) %>% 
  dplyr::mutate_all(~ replace(., is.na(.), 0)) %>%
  # mutate_at(5:7, round, 3) %>% 
  dplyr::rename(sample = rowname,
                # neutrophils = cluster_11,
                endo_cells = endo,
                cells_per_sample = total,
                #neutrophil_percent = sub1_percent
                endo_percent = sub1_percent) %>% 
  dplyr::select(sample, cells_per_sample, total_percent,
                #neutrophils, neutrophil_percent,
                endo_cells, endo_percent) -> samples

## per response
cell_per_r_0 %>% 
  tibble::rownames_to_column() %>%
  dplyr::select(rowname, total) %>% 
  dplyr::left_join((cell_per_r_1 %>%
                    tibble::rownames_to_column() %>%
                    #dplyr::select(rowname, cluster_11)), by = "rowname") %>%
                    dplyr::select(rowname, endo)), by = "rowname") %>% 
  # dplyr::left_join((cell_per_r_2 %>%
  #                   tibble::rownames_to_column() %>%
  #                   dplyr::select(rowname, cluster_50)), by = "rowname") %>% 
  dplyr::mutate(total_percent = round((total / r0) * 100, 3),
                #sub1_percent = (cluster_11 / r1) * 100,
                sub1_percent = round((endo / r1) * 100, 3)
                # sub2_percent = (cluster_50 / r2) * 100
                ) %>% 
  dplyr::mutate_all(~ replace(., is.na(.), 0)) %>%
  # mutate_at(5:7, round, 3) %>% 
  dplyr::rename(response = rowname,
                # neutrophils = cluster_11,
                endo_cells = endo,
                cells_per_response = total,
                #neutrophil_percent = sub1_percent
                endo_percent = sub1_percent) %>% 
  dplyr::select(response, cells_per_response, total_percent,
                #neutrophils, neutrophil_percent,
                endo_cells, endo_percent) -> responses

## per treatment
cell_per_t_0 %>% 
  tibble::rownames_to_column() %>%
  dplyr::select(rowname, total) %>% 
  dplyr::left_join((cell_per_t_1 %>%
                    tibble::rownames_to_column() %>%
                    #dplyr::select(rowname, cluster_11)), by = "rowname") %>%
                    dplyr::select(rowname, endo)), by = "rowname") %>% 
  # dplyr::left_join((cell_per_t_2 %>%
  #                   tibble::rownames_to_column() %>%
  #                   dplyr::select(rowname, cluster_50)), by = "rowname") %>% 
  dplyr::mutate(total_percent = round((total / t0) * 100, 3),
                #sub1_percent = (cluster_11 / t1) * 100,
                sub1_percent = round((endo / t1) * 100, 3),
                # sub2_percent = (cluster_50 / t2) * 100,
                rowname = gsub("_._", "_", rowname)) %>% 
  dplyr::mutate_all(~ replace(., is.na(.), 0)) %>%
  # mutate_at(5:7, round, 3) %>% 
  dplyr::rename(treatment = rowname,
                # neutrophils = cluster_11,
                endo_cells = endo,
                cells_per_treatment = total,
                #neutrophil_percent = sub1_percent
                endo_percent = sub1_percent) %>% 
  dplyr::select(treatment, cells_per_treatment, total_percent,
                #neutrophils, neutrophil_percent,
                endo_cells, endo_percent) -> treatments

### contamination ###########

## load libraries
library("dplyr")
library("ggplot2")
library("writexl")

## read seurat object
data_list <- list()
file_path <- "~/Documents/tmp/pancreatic_melanie/contamination_dbl/"
read.table(paste0(file_path, "markers_bcells_cancer_cont_dbl.txt"),
           header = TRUE, sep = "\t") %>% dplyr::rename(gene = X) %>%
           dplyr::filter(p_val_adj < 0.05) %>% 
           dplyr::arrange(desc(avg_log2FC)) -> data_list[[1]]
read.table(paste0(file_path, "markers_endo_cancer_cont_dbl.txt"),
           header = TRUE, sep = "\t") %>% dplyr::rename(gene = X) %>%
           dplyr::filter(p_val_adj < 0.05) %>% 
           dplyr::arrange(desc(avg_log2FC)) -> data_list[[3]]
read.table(paste0(file_path, "markers_macro_cancer_cont_dbl.txt"),
           header = TRUE, sep = "\t") %>% dplyr::rename(gene = X) %>%
           dplyr::filter(p_val_adj < 0.05) %>% 
           dplyr::arrange(desc(avg_log2FC)) -> data_list[[5]]
read.table(paste0(file_path, "markers_tcells_cancer_cont_dbl.txt"),
           header = TRUE, sep = "\t") %>% dplyr::rename(gene = X) %>%
           dplyr::filter(p_val_adj < 0.05) %>% 
           dplyr::arrange(desc(avg_log2FC)) -> data_list[[7]]

read.table(paste0(file_path, "markers_bcells_cancer_cont_dbl.txt"),
           header = TRUE, sep = "\t") %>% dplyr::rename(gene = X) %>%
           dplyr::filter(p_val_adj < 0.05) %>% 
           dplyr::arrange(avg_log2FC) -> data_list[[2]]
read.table(paste0(file_path, "markers_endo_cancer_cont_dbl.txt"),
           header = TRUE, sep = "\t") %>% dplyr::rename(gene = X) %>%
           dplyr::filter(p_val_adj < 0.05) %>% 
           dplyr::arrange(avg_log2FC) -> data_list[[4]]
read.table(paste0(file_path, "markers_macro_cancer_cont_dbl.txt"),
           header = TRUE, sep = "\t") %>% dplyr::rename(gene = X) %>%
           dplyr::filter(p_val_adj < 0.05) %>% 
           dplyr::arrange(avg_log2FC) -> data_list[[6]]
read.table(paste0(file_path, "markers_tcells_cancer_cont_dbl.txt"),
           header = TRUE, sep = "\t") %>% dplyr::rename(gene = X) %>%
           dplyr::filter(p_val_adj < 0.05) %>% 
           dplyr::arrange(avg_log2FC) -> data_list[[8]]

data_list[[1]]$p_val_adj[data_list[[1]]$p_val_adj == 0] <- 1.0e-300
data_list[[2]]$p_val_adj[data_list[[2]]$p_val_adj == 0] <- 1.0e-300
data_list[[3]]$p_val_adj[data_list[[3]]$p_val_adj == 0] <- 1.0e-300
data_list[[4]]$p_val_adj[data_list[[4]]$p_val_adj == 0] <- 1.0e-300
data_list[[5]]$p_val_adj[data_list[[5]]$p_val_adj == 0] <- 1.0e-300
data_list[[6]]$p_val_adj[data_list[[6]]$p_val_adj == 0] <- 1.0e-300
data_list[[7]]$p_val_adj[data_list[[7]]$p_val_adj == 0] <- 1.0e-300
data_list[[8]]$p_val_adj[data_list[[8]]$p_val_adj == 0] <- 1.0e-300
data_list[[1]]$p_val[data_list[[1]]$p_val == 0] <- 1.0e-300
data_list[[2]]$p_val[data_list[[2]]$p_val == 0] <- 1.0e-300
data_list[[3]]$p_val[data_list[[3]]$p_val == 0] <- 1.0e-300
data_list[[4]]$p_val[data_list[[4]]$p_val == 0] <- 1.0e-300
data_list[[5]]$p_val[data_list[[5]]$p_val == 0] <- 1.0e-300
data_list[[6]]$p_val[data_list[[6]]$p_val == 0] <- 1.0e-300
data_list[[7]]$p_val[data_list[[7]]$p_val == 0] <- 1.0e-300
data_list[[8]]$p_val[data_list[[8]]$p_val == 0] <- 1.0e-300

# names(data_list) <- c("B-Cell contamination", "Endothelial contamination",
#                       "Myeloid contamination", "T-cell contamination")

names(data_list) <- c("B-Cell contamination pos", "B-Cell contamination neg",
                      "Endothelial contamination pos", "Endothelial contamination neg",
                      "Myeloid contamination pos", "Myeloid contamination neg",
                      "T-cell contamination pos", "T-cell contamination neg")
## save file
write_xlsx(path = paste0("contamination_cluster_sig_markers_logFC_posneg.xlsx"), data_list)

