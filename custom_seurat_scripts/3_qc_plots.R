################################# QC plots #####################################
## look at doublet states to ensure that there are no doublet clusters
# read.table(file = paste0(local_path, "pann_per_C_Doublet_Singlet_", proj_name,
#                          ".txt"), header = TRUE) %>% dplyr::select(-X) -> pANN

## seurat metadata file
# read.table(paste0(local_path, "cells_metadata_", proj_name, ".txt"),
#            header = TRUE) %>% dplyr::rename(UMI = X) -> seu_metadata

library("dplyr")
library("ggplot2")
library("tidyselect")

## evergreen variables
proj_name <- "endo_final"
proj_version <- "filtered"

pdf(paste0("QC_histograms_", proj_name, "_", proj_version, ".pdf")) 
## relative distribution cluster (x = pct_sample, y = cluster) SAMPLE
# read.table(paste0(local_path, "cell_per_c_per_O_", proj_name, ".txt"), header = TRUE) %>%
read.table(paste0("cell_per_c_per_O_", proj_name, ".txt"), header = TRUE) %>%
  # dplyr::mutate(num_cells = rowSums(across(CCB023:CCB104))) %>%
  dplyr::mutate(num_cells = rowSums(across(CCB023:CCB103))) %>%
  dplyr::rename(cluster = X) %>%
  dplyr::mutate(cluster = as.factor(cluster)) %>% 
  tidyr::gather(-cluster, -num_cells, value = cells, key = sample) %>%
  dplyr::mutate(pct_sample = (cells/num_cells)) %>%
  ggplot2::ggplot(aes(x = cluster, y = pct_sample, fill = sample)) +
  geom_bar(stat = "identity") + coord_flip() + theme_classic() +
  ggtitle("Relative Distribution of Samples by Cluster")

## absolute distribution cluster (x = num_cells, y = cluster) SAMPLE
# read.table(paste0(local_path, "cell_per_c_per_O_", proj_name, ".txt"), header = TRUE) %>%
read.table(paste0("cell_per_c_per_O_", proj_name, ".txt"), header = TRUE) %>%
  # dplyr::mutate(num_cells = rowSums(across(CCB023:CCB104))) %>%
  dplyr::mutate(num_cells = rowSums(across(CCB023:CCB103))) %>%
  dplyr::rename(cluster = X) %>%
  dplyr::mutate(cluster = as.factor(cluster)) %>% 
  tidyr::gather(-cluster, -num_cells, value = cells, key = sample) %>%
  dplyr::mutate(pct_sample = (cells/num_cells)) %>%
  ggplot2::ggplot(aes(x = cluster, y = cells, fill = sample)) +
  geom_bar(stat = "identity") + coord_flip() + theme_classic() +
  ggtitle("Absolute Distribution of Samples by Cluster") 

## relative distribution sample (x = pct_cluster, y = sample) SAMPLE
# read.table(paste0(local_path, "cell_per_c_per_O_", proj_name, ".txt"), header = TRUE) %>%
read.table(paste0("cell_per_c_per_O_", proj_name, ".txt"), header = TRUE) %>%
  dplyr::mutate(X = paste0("x", X)) %>% dplyr::rename(cluster = X) -> rel1
col1 <- rel1$cluster
rel1 %>% dplyr::select(-cluster) %>% t() %>% as.data.frame() -> rel1
names(rel1) <- col1
#names(rel1) <- gsub("x", "", col1) 
rel1 %>%
  tibble::rownames_to_column() %>% 
  dplyr::mutate(num_cells = rowSums(across(x0:tidyselect::last_col()))) %>% 
  # dplyr::mutate(num_cells = rowSums(across(`xM2 Macrophages`:tidyselect::last_col()))) %>% 
  tidyr::gather(-rowname, -num_cells, value = cells, key = cluster) %>% 
  dplyr::rename(sample = rowname) %>% 
  dplyr::mutate(pct_cluster = (cells/num_cells)) %>% 
  dplyr::mutate(cluster = as.factor(gsub("x", "", cluster))) %>% 
  ggplot2::ggplot(aes(x = sample, y = pct_cluster, fill = cluster)) +
  geom_bar(stat = "identity") + coord_flip() + theme_classic() +
  ggtitle("Relative Distribution of Clusters by Samples")

## absolute distribution sample (x = num_cells, y = sample) SAMPLE
# read.table(paste0(local_path, "cell_per_c_per_O_", proj_name, ".txt"), header = TRUE) %>%
read.table(paste0("cell_per_c_per_O_", proj_name, ".txt"), header = TRUE) %>%
  dplyr::mutate(X = paste0("x", X)) %>% dplyr::rename(cluster = X) -> rel1
col1 <- rel1$cluster
rel1 %>% dplyr::select(-cluster) %>% t() %>% as.data.frame() -> rel1
names(rel1) <- col1
rel1 %>% tibble::rownames_to_column() %>% 
  dplyr::mutate(num_cells = rowSums(across(x0:tidyselect::last_col()))) %>%
  # dplyr::mutate(num_cells = rowSums(across(`xM2 Macrophages`:tidyselect::last_col()))) %>%
  tidyr::gather(-rowname, -num_cells, value = cells, key = cluster) %>% 
  dplyr::mutate(pct_cluster = (cells/num_cells)) %>% 
  dplyr::mutate(cluster = as.factor(gsub("x", "", cluster))) %>% 
  dplyr::rename(sample = rowname) %>% 
  ggplot2::ggplot(aes(x = sample, y = cells, fill = cluster)) +
  geom_bar(stat = "identity") + coord_flip() + theme_classic() +
  ggtitle("Absolute Distribution of Clusters by Samples")

## relative distribution cluster (x = pct_sample, y = cluster) RESPONSE
# read.table(paste0(local_path, "cell_per_c_per_R_", proj_name, ".txt"), header = TRUE) %>%
read.table(paste0("cell_per_c_per_R_", proj_name, ".txt"), header = TRUE) %>%
  dplyr::mutate(num_cells = rowSums(across(Relapsing:Untreated))) %>%
  dplyr::rename(cluster = X) %>%
  dplyr::mutate(cluster = as.factor(cluster)) %>% 
  tidyr::gather(-cluster, -num_cells, value = cells, key = response) %>%
  dplyr::mutate(pct_response = (cells/num_cells)) %>%
  ggplot2::ggplot(aes(x = cluster, y = pct_response, fill = response)) +
  geom_bar(stat = "identity") + coord_flip() + theme_classic() +
  ggtitle("Relative Distribution of Response by Cluster")

## absolute distribution cluster (x = num_cells, y = cluster) RESPONSE
# read.table(paste0(local_path, "cell_per_c_per_R_", proj_name, ".txt"), header = TRUE) %>%
read.table(paste0("cell_per_c_per_R_", proj_name, ".txt"), header = TRUE) %>%
  dplyr::mutate(num_cells = rowSums(across(Relapsing:Untreated))) %>%
  dplyr::rename(cluster = X) %>%
  dplyr::mutate(cluster = as.factor(cluster)) %>% 
  tidyr::gather(-cluster, -num_cells, value = cells, key = response) %>%
  dplyr::mutate(pct_response = (cells/num_cells)) %>%
  ggplot2::ggplot(aes(x = cluster, y = cells, fill = response)) +
  geom_bar(stat = "identity") + coord_flip() + theme_classic() +
  ggtitle("Absolute Distribution of Response by Cluster")

## relative distribution sample (x = pct_cluster, y = sample) RESPONSE
#read.table(paste0(local_path, "cell_per_c_per_R_", proj_name, ".txt"), header = TRUE) %>%
read.table(paste0("cell_per_c_per_R_", proj_name, ".txt"), header = TRUE) %>%
  dplyr::mutate(X = paste0("x", X)) %>% dplyr::rename(cluster = X) -> rel1
col1 <- rel1$cluster
rel1 %>% dplyr::select(-cluster) %>% t() %>% as.data.frame() -> rel1
names(rel1) <- col1
rel1 %>% tibble::rownames_to_column() %>% 
  dplyr::mutate(num_cells = rowSums(across(x0:tidyselect::last_col()))) %>%
  # dplyr::mutate(num_cells = rowSums(across(`xM2 Macrophages`:tidyselect::last_col()))) %>%
  tidyr::gather(-rowname, -num_cells, value = cells, key = cluster) %>% 
  dplyr::mutate(pct_cluster = (cells/num_cells)) %>% 
  dplyr::mutate(cluster = as.factor(gsub("x", "", cluster))) %>% 
  dplyr::rename(response = rowname) %>% 
  ggplot2::ggplot(aes(x = response, y = pct_cluster, fill = cluster)) +
  geom_bar(stat = "identity") + coord_flip() + theme_classic() +
  ggtitle("Relative Distribution of Clusters by Response")

## absolute distribution sample (x = num_cells, y = sample) RESPONSE
# read.table(paste0(local_path, "cell_per_c_per_R_", proj_name, ".txt"), header = TRUE) %>%
read.table(paste0("cell_per_c_per_R_", proj_name, ".txt"), header = TRUE) %>%
  dplyr::mutate(X = paste0("x", X)) %>% dplyr::rename(cluster = X) -> rel1
col1 <- rel1$cluster
rel1 %>% dplyr::select(-cluster) %>% t() %>% as.data.frame() -> rel1
names(rel1) <- col1
rel1 %>% tibble::rownames_to_column() %>% 
  dplyr::mutate(num_cells = rowSums(across(x0:tidyselect::last_col()))) %>% 
  # dplyr::mutate(num_cells = rowSums(across(`xM2 Macrophages`:tidyselect::last_col()))) %>%
  tidyr::gather(-rowname, -num_cells, value = cells, key = cluster) %>% 
  dplyr::mutate(pct_cluster = (cells/num_cells)) %>% 
  dplyr::mutate(cluster = as.factor(gsub("x", "", cluster))) %>% 
  dplyr::rename(response = rowname) %>% 
  ggplot2::ggplot(aes(x = response, y = cells, fill = cluster)) +
  geom_bar(stat = "identity") + coord_flip() + theme_classic() +
  ggtitle("Absolute Distribution of Clusters by Response")

## relative distribution cluster (x = pct_sample, y = cluster) TREATMENT
# read.table(paste0(local_path, "cell_per_c_per_T_", proj_name, ".txt"), header = TRUE) %>%
read.table(paste0("cell_per_c_per_T_", proj_name, ".txt"), header = TRUE) %>%
  dplyr::mutate(num_cells = rowSums(across(B20:Untreated))) %>%
  dplyr::rename(cluster = X) %>%
  dplyr::mutate(cluster = as.factor(cluster)) %>% 
  tidyr::gather(-cluster, -num_cells, value = cells, key = treatment) %>%
  dplyr::mutate(pct_treatment = (cells/num_cells)) %>%
  ggplot2::ggplot(aes(x = cluster, y = pct_treatment, fill = treatment)) +
  geom_bar(stat = "identity") + coord_flip() + theme_classic() +
  ggtitle("Relative Distribution of Treatment by Cluster")

## absolute distribution cluster (x = num_cells, y = cluster) TREATMENT
# read.table(paste0(local_path, "cell_per_c_per_T_", proj_name, ".txt"), header = TRUE) %>%
read.table(paste0("cell_per_c_per_T_", proj_name, ".txt"), header = TRUE) %>%
  dplyr::mutate(num_cells = rowSums(across(B20:Untreated))) %>%
  dplyr::rename(cluster = X) %>%
  dplyr::mutate(cluster = as.factor(cluster)) %>% 
  tidyr::gather(-cluster, -num_cells, value = cells, key = treatment) %>%
  dplyr::mutate(pct_treatment = (cells/num_cells)) %>%
  ggplot2::ggplot(aes(x = cluster, y = cells, fill = treatment)) +
  geom_bar(stat = "identity") + coord_flip() + theme_classic() +
  ggtitle("Absolute Distribution of Treatment by Cluster")

## relative distribution sample (x = pct_cluster, y = sample) TREATMENT
# read.table(paste0(local_path, "cell_per_c_per_T_", proj_name, ".txt"), header = TRUE) %>%
read.table(paste0("cell_per_c_per_T_", proj_name, ".txt"), header = TRUE) %>%
  dplyr::mutate(X = paste0("x", X)) %>% dplyr::rename(cluster = X) -> rel1
col1 <- rel1$cluster
rel1 %>% dplyr::select(-cluster) %>% t() %>% as.data.frame() -> rel1
names(rel1) <- col1
rel1 %>% tibble::rownames_to_column() %>% 
  dplyr::mutate(num_cells = rowSums(across(x0:tidyselect::last_col()))) %>% 
  # dplyr::mutate(num_cells = rowSums(across(`xM2 Macrophages`:tidyselect::last_col()))) %>%
  tidyr::gather(-rowname, -num_cells, value = cells, key = cluster) %>% 
  dplyr::mutate(pct_cluster = (cells/num_cells)) %>% 
  dplyr::mutate(cluster = as.factor(gsub("x", "", cluster))) %>% 
  dplyr::rename(treatment = rowname) %>% 
  ggplot2::ggplot(aes(x = treatment, y = pct_cluster, fill = cluster)) +
  geom_bar(stat = "identity") + coord_flip() + theme_classic() +
  ggtitle("Relative Distribution of Clusters by Treatment")

## absolute distribution sample (x = num_cells, y = sample) TREATMENT
# read.table(paste0(local_path, "cell_per_c_per_T_", proj_name, ".txt"), header = TRUE) %>%
read.table(paste0("cell_per_c_per_T_", proj_name, ".txt"), header = TRUE) %>%
  dplyr::mutate(X = paste0("x", X)) %>% dplyr::rename(cluster = X) -> rel1
col1 <- rel1$cluster
rel1 %>% dplyr::select(-cluster) %>% t() %>% as.data.frame() -> rel1
names(rel1) <- col1
rel1 %>% tibble::rownames_to_column() %>% 
  dplyr::mutate(num_cells = rowSums(across(x0:tidyselect::last_col()))) %>% 
  # dplyr::mutate(num_cells = rowSums(across(`xM2 Macrophages`:tidyselect::last_col()))) %>%
  tidyr::gather(-rowname, -num_cells, value = cells, key = cluster) %>% 
  dplyr::mutate(pct_cluster = (cells/num_cells)) %>% 
  dplyr::mutate(cluster = as.factor(gsub("x", "", cluster))) %>% 
  dplyr::rename(treatment = rowname) %>% 
  ggplot2::ggplot(aes(x = treatment, y = cells, fill = cluster)) +
  geom_bar(stat = "identity") + coord_flip() + theme_classic() +
  ggtitle("Absolute Distribution of Clusters by Treatment")
dev.off()

############################# annotated QC plots  ##############################
proj_version <- "filtered"

pdf(paste0("QC_histograms_", proj_name, "_", proj_version, ".pdf")) 
## relative distribution cluster (x = pct_sample, y = cluster) SAMPLE
# read.table(paste0(local_path, "cell_per_c_per_O_", proj_name, ".txt"), header = TRUE) %>%
read.table(paste0("cell_per_c_per_O_", proj_name, ".txt"), header = TRUE) %>%
  # dplyr::mutate(num_cells = rowSums(across(CCB023:CCB104))) %>%
  dplyr::mutate(num_cells = rowSums(across(CCB023:CCB103))) %>%
  dplyr::rename(cluster = X) %>%
  dplyr::mutate(cluster = as.factor(cluster)) %>% 
  tidyr::gather(-cluster, -num_cells, value = cells, key = sample) %>%
  dplyr::mutate(pct_sample = (cells/num_cells)) %>%
  ggplot2::ggplot(aes(x = cluster, y = pct_sample, fill = sample)) +
  geom_bar(stat = "identity") + coord_flip() + theme_classic() +
  ggtitle("Relative Distribution of Samples by Cluster")

## absolute distribution cluster (x = num_cells, y = cluster) SAMPLE
# read.table(paste0(local_path, "cell_per_c_per_O_", proj_name, ".txt"), header = TRUE) %>%
read.table(paste0("cell_per_c_per_O_", proj_name, ".txt"), header = TRUE) %>%
  # dplyr::mutate(num_cells = rowSums(across(CCB023:CCB104))) %>%
  dplyr::mutate(num_cells = rowSums(across(CCB023:CCB103))) %>%
  dplyr::rename(cluster = X) %>%
  dplyr::mutate(cluster = as.factor(cluster)) %>% 
  tidyr::gather(-cluster, -num_cells, value = cells, key = sample) %>%
  dplyr::mutate(pct_sample = (cells/num_cells)) %>%
  ggplot2::ggplot(aes(x = cluster, y = cells, fill = sample)) +
  geom_bar(stat = "identity") + coord_flip() + theme_classic() +
  ggtitle("Absolute Distribution of Samples by Cluster") 

## relative distribution sample (x = pct_cluster, y = sample) SAMPLE
# read.table(paste0(local_path, "cell_per_c_per_O_", proj_name, ".txt"), header = TRUE) %>%
read.table(paste0("cell_per_c_per_O_", proj_name, ".txt"), header = TRUE) %>%
  dplyr::mutate(X = paste0("x", X)) %>% dplyr::rename(cluster = X) -> rel1
col1 <- rel1$cluster
rel1 %>% dplyr::select(-cluster) %>% t() %>% as.data.frame() -> rel1
names(rel1) <- col1
#names(rel1) <- gsub("x", "", col1) 
rel1 %>%
  tibble::rownames_to_column() %>% 
  #dplyr::mutate(num_cells = rowSums(across(x0:tidyselect::last_col()))) %>% 
  dplyr::mutate(num_cells = rowSums(across(`xM2 Macrophages`:tidyselect::last_col()))) %>% 
  tidyr::gather(-rowname, -num_cells, value = cells, key = cluster) %>% 
  dplyr::rename(sample = rowname) %>% 
  dplyr::mutate(pct_cluster = (cells/num_cells)) %>% 
  dplyr::mutate(cluster = as.factor(gsub("x", "", cluster))) %>% 
  ggplot2::ggplot(aes(x = sample, y = pct_cluster, fill = cluster)) +
  geom_bar(stat = "identity") + coord_flip() + theme_classic() +
  ggtitle("Relative Distribution of Clusters by Samples")

## absolute distribution sample (x = num_cells, y = sample) SAMPLE
# read.table(paste0(local_path, "cell_per_c_per_O_", proj_name, ".txt"), header = TRUE) %>%
read.table(paste0("cell_per_c_per_O_", proj_name, ".txt"), header = TRUE) %>%
  dplyr::mutate(X = paste0("x", X)) %>% dplyr::rename(cluster = X) -> rel1
col1 <- rel1$cluster
rel1 %>% dplyr::select(-cluster) %>% t() %>% as.data.frame() -> rel1
names(rel1) <- col1
rel1 %>% tibble::rownames_to_column() %>% 
  # dplyr::mutate(num_cells = rowSums(across(x0:tidyselect::last_col()))) %>%
  dplyr::mutate(num_cells = rowSums(across(`xM2 Macrophages`:tidyselect::last_col()))) %>%
  tidyr::gather(-rowname, -num_cells, value = cells, key = cluster) %>% 
  dplyr::mutate(pct_cluster = (cells/num_cells)) %>% 
  dplyr::mutate(cluster = as.factor(gsub("x", "", cluster))) %>% 
  dplyr::rename(sample = rowname) %>% 
  ggplot2::ggplot(aes(x = sample, y = cells, fill = cluster)) +
  geom_bar(stat = "identity") + coord_flip() + theme_classic() +
  ggtitle("Absolute Distribution of Clusters by Samples")

## relative distribution cluster (x = pct_sample, y = cluster) RESPONSE
# read.table(paste0(local_path, "cell_per_c_per_R_", proj_name, ".txt"), header = TRUE) %>%
read.table(paste0("cell_per_c_per_R_", proj_name, ".txt"), header = TRUE) %>%
  dplyr::mutate(num_cells = rowSums(across(Relapsing:Untreated))) %>%
  dplyr::rename(cluster = X) %>%
  dplyr::mutate(cluster = as.factor(cluster)) %>% 
  tidyr::gather(-cluster, -num_cells, value = cells, key = response) %>%
  dplyr::mutate(pct_response = (cells/num_cells)) %>%
  ggplot2::ggplot(aes(x = cluster, y = pct_response, fill = response)) +
  geom_bar(stat = "identity") + coord_flip() + theme_classic() +
  ggtitle("Relative Distribution of Response by Cluster")

## absolute distribution cluster (x = num_cells, y = cluster) RESPONSE
# read.table(paste0(local_path, "cell_per_c_per_R_", proj_name, ".txt"), header = TRUE) %>%
read.table(paste0("cell_per_c_per_R_", proj_name, ".txt"), header = TRUE) %>%
  dplyr::mutate(num_cells = rowSums(across(Relapsing:Untreated))) %>%
  dplyr::rename(cluster = X) %>%
  dplyr::mutate(cluster = as.factor(cluster)) %>% 
  tidyr::gather(-cluster, -num_cells, value = cells, key = response) %>%
  dplyr::mutate(pct_response = (cells/num_cells)) %>%
  ggplot2::ggplot(aes(x = cluster, y = cells, fill = response)) +
  geom_bar(stat = "identity") + coord_flip() + theme_classic() +
  ggtitle("Absolute Distribution of Response by Cluster")

## relative distribution sample (x = pct_cluster, y = sample) RESPONSE
#read.table(paste0(local_path, "cell_per_c_per_R_", proj_name, ".txt"), header = TRUE) %>%
read.table(paste0("cell_per_c_per_R_", proj_name, ".txt"), header = TRUE) %>%
  dplyr::mutate(X = paste0("x", X)) %>% dplyr::rename(cluster = X) -> rel1
col1 <- rel1$cluster
rel1 %>% dplyr::select(-cluster) %>% t() %>% as.data.frame() -> rel1
names(rel1) <- col1
rel1 %>% tibble::rownames_to_column() %>% 
  # dplyr::mutate(num_cells = rowSums(across(x0:tidyselect::last_col()))) %>%
  dplyr::mutate(num_cells = rowSums(across(`xM2 Macrophages`:tidyselect::last_col()))) %>%
  tidyr::gather(-rowname, -num_cells, value = cells, key = cluster) %>% 
  dplyr::mutate(pct_cluster = (cells/num_cells)) %>% 
  dplyr::mutate(cluster = as.factor(gsub("x", "", cluster))) %>% 
  dplyr::rename(response = rowname) %>% 
  ggplot2::ggplot(aes(x = response, y = pct_cluster, fill = cluster)) +
  geom_bar(stat = "identity") + coord_flip() + theme_classic() +
  ggtitle("Relative Distribution of Clusters by Response")

## absolute distribution sample (x = num_cells, y = sample) RESPONSE
# read.table(paste0(local_path, "cell_per_c_per_R_", proj_name, ".txt"), header = TRUE) %>%
read.table(paste0("cell_per_c_per_R_", proj_name, ".txt"), header = TRUE) %>%
  dplyr::mutate(X = paste0("x", X)) %>% dplyr::rename(cluster = X) -> rel1
col1 <- rel1$cluster
rel1 %>% dplyr::select(-cluster) %>% t() %>% as.data.frame() -> rel1
names(rel1) <- col1
rel1 %>% tibble::rownames_to_column() %>% 
  # dplyr::mutate(num_cells = rowSums(across(x0:tidyselect::last_col()))) %>% 
  dplyr::mutate(num_cells = rowSums(across(`xM2 Macrophages`:tidyselect::last_col()))) %>%
  tidyr::gather(-rowname, -num_cells, value = cells, key = cluster) %>% 
  dplyr::mutate(pct_cluster = (cells/num_cells)) %>% 
  dplyr::mutate(cluster = as.factor(gsub("x", "", cluster))) %>% 
  dplyr::rename(response = rowname) %>% 
  ggplot2::ggplot(aes(x = response, y = cells, fill = cluster)) +
  geom_bar(stat = "identity") + coord_flip() + theme_classic() +
  ggtitle("Absolute Distribution of Clusters by Response")

## relative distribution cluster (x = pct_sample, y = cluster) TREATMENT
# read.table(paste0(local_path, "cell_per_c_per_T_", proj_name, ".txt"), header = TRUE) %>%
read.table(paste0("cell_per_c_per_T_", proj_name, ".txt"), header = TRUE) %>%
  dplyr::mutate(num_cells = rowSums(across(B20:Untreated))) %>%
  dplyr::rename(cluster = X) %>%
  dplyr::mutate(cluster = as.factor(cluster)) %>% 
  tidyr::gather(-cluster, -num_cells, value = cells, key = treatment) %>%
  dplyr::mutate(pct_treatment = (cells/num_cells)) %>%
  ggplot2::ggplot(aes(x = cluster, y = pct_treatment, fill = treatment)) +
  geom_bar(stat = "identity") + coord_flip() + theme_classic() +
  ggtitle("Relative Distribution of Treatment by Cluster")

## absolute distribution cluster (x = num_cells, y = cluster) TREATMENT
# read.table(paste0(local_path, "cell_per_c_per_T_", proj_name, ".txt"), header = TRUE) %>%
read.table(paste0("cell_per_c_per_T_", proj_name, ".txt"), header = TRUE) %>%
  dplyr::mutate(num_cells = rowSums(across(B20:Untreated))) %>%
  dplyr::rename(cluster = X) %>%
  dplyr::mutate(cluster = as.factor(cluster)) %>% 
  tidyr::gather(-cluster, -num_cells, value = cells, key = treatment) %>%
  dplyr::mutate(pct_treatment = (cells/num_cells)) %>%
  ggplot2::ggplot(aes(x = cluster, y = cells, fill = treatment)) +
  geom_bar(stat = "identity") + coord_flip() + theme_classic() +
  ggtitle("Absolute Distribution of Treatment by Cluster")

## relative distribution sample (x = pct_cluster, y = sample) TREATMENT
# read.table(paste0(local_path, "cell_per_c_per_T_", proj_name, ".txt"), header = TRUE) %>%
read.table(paste0("cell_per_c_per_T_", proj_name, ".txt"), header = TRUE) %>%
  dplyr::mutate(X = paste0("x", X)) %>% dplyr::rename(cluster = X) -> rel1
col1 <- rel1$cluster
rel1 %>% dplyr::select(-cluster) %>% t() %>% as.data.frame() -> rel1
names(rel1) <- col1
rel1 %>% tibble::rownames_to_column() %>% 
  # dplyr::mutate(num_cells = rowSums(across(x0:tidyselect::last_col()))) %>% 
  dplyr::mutate(num_cells = rowSums(across(`xM2 Macrophages`:tidyselect::last_col()))) %>%
  tidyr::gather(-rowname, -num_cells, value = cells, key = cluster) %>% 
  dplyr::mutate(pct_cluster = (cells/num_cells)) %>% 
  dplyr::mutate(cluster = as.factor(gsub("x", "", cluster))) %>% 
  dplyr::rename(treatment = rowname) %>% 
  ggplot2::ggplot(aes(x = treatment, y = pct_cluster, fill = cluster)) +
  geom_bar(stat = "identity") + coord_flip() + theme_classic() +
  ggtitle("Relative Distribution of Clusters by Treatment")

## absolute distribution sample (x = num_cells, y = sample) TREATMENT
# read.table(paste0(local_path, "cell_per_c_per_T_", proj_name, ".txt"), header = TRUE) %>%
read.table(paste0("cell_per_c_per_T_", proj_name, ".txt"), header = TRUE) %>%
  dplyr::mutate(X = paste0("x", X)) %>% dplyr::rename(cluster = X) -> rel1
col1 <- rel1$cluster
rel1 %>% dplyr::select(-cluster) %>% t() %>% as.data.frame() -> rel1
names(rel1) <- col1
rel1 %>% tibble::rownames_to_column() %>% 
  # dplyr::mutate(num_cells = rowSums(across(x0:tidyselect::last_col()))) %>% 
  dplyr::mutate(num_cells = rowSums(across(`xM2 Macrophages`:tidyselect::last_col()))) %>%
  tidyr::gather(-rowname, -num_cells, value = cells, key = cluster) %>% 
  dplyr::mutate(pct_cluster = (cells/num_cells)) %>% 
  dplyr::mutate(cluster = as.factor(gsub("x", "", cluster))) %>% 
  dplyr::rename(treatment = rowname) %>% 
  ggplot2::ggplot(aes(x = treatment, y = cells, fill = cluster)) +
  geom_bar(stat = "identity") + coord_flip() + theme_classic() +
  ggtitle("Absolute Distribution of Clusters by Treatment")
dev.off()

############################# QC plots for filter ##############################
## read metadata file
readRDS(file = "filter_metadata.rds") %>% 
  rename(nBarcode = nUMI) -> metadata

pdf(paste0("QC_plots_", proj_name, "_", proj_version, ".pdf")) 
# Visualize the number of cell counts per TREATMENT
metadata %>% 
  	ggplot(aes(x = Treatment, fill = Treatment)) + 
  	geom_bar() + theme_classic() +
  	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  	theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  	ggtitle("Number of Barcodes per treatment")

# Visualize the number of cell counts per SAMPLE
metadata %>%
  	ggplot(aes(x = sample, fill = sample)) +
  	geom_bar() + theme_classic() +
    geom_hline(yintercept = 10000, linetype = "dashed", color = "red") +
  	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  	theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  	ggtitle("Number of Barcodes per sample")

# Visualize the number Barcodes/transcripts per cell by TREATMENT
metadata %>% 
  	ggplot(aes(color = Treatment, x = nBarcode, fill = Treatment)) + 
  	geom_density(alpha = 0.2) + scale_x_log10() + theme_classic() +
  	# ylab("Cell density") +
    xlab("nCount_RNA") +
    geom_vline(xintercept = 1000, linetype = "dashed", color = "red") +
  	ggtitle("Number nCount_RNA per Barcode by Treatment") +
    theme(plot.title = element_text(hjust = 0.5)) +
    annotate(geom = "text", x = 5000, y = 2, color = "red", 
             label = "Filter threshold: \n nCount_RNA = 1000")

# Visualize the number Barcodes/transcripts per cell by SAMPLE
metadata %>% 
  	ggplot(aes(color = sample, x = nBarcode, fill = sample)) + 
  	geom_density(alpha = 0.2) + scale_x_log10() + theme_classic() +
  	# ylab("Cell density") +
    xlab("nCount_RNA") +
    geom_vline(xintercept = 1000, linetype = "dashed", color = "red") +
  	ggtitle("Number nCount_RNA per Barcode by Sample") +
    theme(plot.title = element_text(hjust = 0.5)) +
    annotate(geom = "text", x = 5000, y = 2.5, color = "red",
             label = "Filter threshold: \n nCount_RNA = 1000")

# Visualize the distribution of genes detected per cell via histogram by treatment
metadata %>%
  	ggplot(aes(color = Treatment, x = nGene, fill = Treatment)) + 
  	geom_density(alpha = 0.2) + theme_classic() + scale_x_log10() + 
    xlab("nFeature_RNA") +
  	geom_vline(xintercept = 1000, linetype = "dashed", color = "red") +
  	ggtitle("Number nFeature_RNA per Barcode by Treatment") +
    theme(plot.title = element_text(hjust = 0.5)) +
    annotate(geom = "text", x = 2000, y = 3, color = "red",
             label = "Filter threshold: \n nFeature_RNA \n = 1000")

# Visualize the distribution of genes detected per cell via histogram by SAMPLE
metadata %>%
  	ggplot(aes(color = sample, x = nGene, fill = sample)) + 
  	geom_density(alpha = 0.2) + theme_classic() + scale_x_log10() + 
    xlab("nFeature_RNA") +
  	geom_vline(xintercept = 1000, linetype = "dashed", color = "red") +
  	ggtitle("Number nFeature_RNA per Barcode by Sample") +
    theme(plot.title = element_text(hjust = 0.5)) +
    annotate(geom = "text", x = 2000, y = 4, color = "red", 
             label = "Filter threshold: \n nFeature_RNA \n = 1000")

# Visualize the distribution of genes detected per cell via boxplot by treatment
metadata %>% 
  	ggplot(aes(x = Treatment, y = log10(nGene), fill = Treatment)) + 
  	geom_boxplot() + theme_classic() + ylab("log10(nFeature_RNA)") +
  	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  	theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  	ggtitle("Barplot log10 nFeature_RNA by treatment")

# Visualize the distribution of genes detected per cell via boxplot by sample
metadata %>% 
  	ggplot(aes(x = sample, y = log10(nGene), fill = sample)) + 
  	geom_boxplot() + theme_classic() +
    ylab("log10(nFeature_RNA)") +
  	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  	theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  	ggtitle("Barplot log10 nFeature_RNA by sample")
# 
# # Visualize the correlation between genes detected and number of Barcodes and
# # determine whether strong presence of cells with low numbers of genes/Barcodes
# ## by treatment
# metadata %>% 
#   	ggplot(aes(x = nBarcode, y = nGene, color = mito_ratio)) + 
#   	geom_point() + 
# 	  scale_colour_gradient(low = "gray90", high = "black") +
#   	stat_smooth(method = lm) +
#   	xlab("log10(nCount_RNA)") + ylab("log10(nFeature_RNA)") +
#   	scale_x_log10() + scale_y_log10() + theme_classic() +
#   	geom_vline(xintercept = 500) + geom_hline(yintercept = 250) +
#   	facet_wrap(~Treatment) +
#   	ggtitle("Scatterplots of log10 nCount_RNA vs log10 nFeature_RNA by treatment") +
#     theme(plot.title = element_text(hjust = 0.5))
# 
# # Visualize the correlation between genes detected and number of Barcodes and
# # determine whether strong presence of cells with low numbers of genes/Barcodes
# # by sample
# metadata %>% 
#   	ggplot(aes(x = nBarcode, y = nGene, color = mito_ratio)) + 
#   	geom_point(size = 0.5) + 
# 	  scale_colour_gradient(low = "gray90", high = "black") +
#   	stat_smooth(method = lm) +
#   	xlab("log10(nCount_RNA)") + ylab("log10(nFeature_RNA)") +
#   	scale_x_log10() + scale_y_log10() + theme_classic() +
#   	geom_vline(xintercept = 500) + geom_hline(yintercept = 250) +
#   	facet_wrap(~sample) +
#   	ggtitle("Scatterplots of log10 nCount_RNA vs log10 nFeature_RNA by sample") +
#     theme(plot.title = element_text(hjust = 0.5))
# 
# # Visualize the distribution of mitochondrial gene expression detected per cell
# # by treatment
# metadata %>% 
#   	ggplot(aes(color = Treatment, x = mito_ratio, fill = Treatment)) + 
#   	geom_density(alpha = 0.2) + 
#   	# scale_x_log10() +
#     xlab("Percent Mitochondrial Genes") + theme_classic() +
#     geom_vline(xintercept = 0.2, linetype = "dashed", color = "red") +
#     theme(plot.title = element_text(hjust = 0.5)) +
#     annotate(geom = "text", x = 0.15, y = 25, color = "red",
#              label = "Filter threshold: \n percent.mt = 20") +
#   	ggtitle("Distribution of percent mitochondrial genes by treatment")
# 
# # Visualize the distribution of mitochondrial gene expression detected per cell
# # by sample
# metadata %>% 
#   	ggplot(aes(color = sample, x = mito_ratio, fill = sample)) + 
#   	geom_density(alpha = 0.2) + 
#   	# scale_x_log10() +
#     xlab("Percent Mitochondrial Genes") + theme_classic() +
#     geom_vline(xintercept = 0.2, linetype = "dashed", color = "red") +
#     theme(plot.title = element_text(hjust = 0.5)) +
#     annotate(geom = "text", x = 0.15, y = 35, color = "red",
#              label = "Filter threshold: \n percent.mt = 20") +
#   	ggtitle("Distribution of percent mitochondrial genes by sample")
# 
# # Visualize the overall complexity of the gene expression by visualizing the
# # genes detected per Barcode by treatment
# metadata %>%
#   	ggplot(aes(x = log10GenesPerUMI, color = Treatment, fill = Treatment)) +
#   	geom_density(alpha = 0.2) + theme_classic() +
#   	geom_vline(xintercept = 0.8, linetype = "dashed", color = "red") +
#     theme(plot.title = element_text(hjust = 0.5)) +
#     annotate(geom = "text", x = 0.7, y = 12, color = "red",
#              label = "Complexity \n threshold = 80%") +
#   	ggtitle("Distribution of log10GenesPerUMI by treatment")
# 
# # Visualize the overall complexity of the gene expression by visualizing the
# # genes detected per UMI by sample
# metadata %>%
#   	ggplot(aes(x = log10GenesPerUMI, color = sample, fill = sample)) +
#   	geom_density(alpha = 0.2) + theme_classic() +
#   	geom_vline(xintercept = 0.8, linetype = "dashed", color = "red") +
#     theme(plot.title = element_text(hjust = 0.5)) +
#     annotate(geom = "text", x = 0.7, y = 15, color = "red",
#              label = "Complexity \n threshold = 80%") +
#   	ggtitle("Distribution of log10GenesPerUMI by sample")
dev.off()
