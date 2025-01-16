################################ relative cell populations #####################
## this script takes the numbers of cells by sample, treatment and response for
## a given Seurat object, and creates bar plots of the relative amount of cells
## out of the total numbers of cells by sample, treatment and response
library("dplyr")
library("ggplot2")
library("tidyselect")

### load dataset
# seu <- readRDS("pancreatic_melanie/final_annotations_recluster_v4_correct.rds")
# proj_name <- "Whole High Filter Data"
# seu <- readRDS("pancreatic_melanie/tcell1/after_ann/tcell_recluster_actual_2.rds")
# proj_name <- "T-cells"
# seu <- readRDS("pancreatic_melanie/macro2/macro_recluster_actual_2.rds")
# proj_name <- "Myeloids"
#seu <- readRDS("pancreatic_melanie/new_data/endo1/endo_recluster_actual.rds")
# proj_name <- "Endothelial Cells"
# proj_name <- "Cancer"
proj_name <- "Neutrophils"

## ensure that the active idents are set to the correct annotations
#seu <- SetIdent(object = seu, value = "final_primary_annotations")
# seu <- SetIdent(object = seu, value = "final_subtype_annotations")

## make treatment_response variable in metadata
seu@meta.data$treatment_response <- paste0(seu@meta.data$Treatment, "_",
                                           seu@meta.data$Response)
seu@meta.data$treatment_response <- gsub("Untreated_Untreated", "Untreated",
                                         seu@meta.data$treatment_response)

seu@meta.data$treatment_response_samp <- paste0(seu@meta.data$Treatment, "_",
                                           seu@meta.data$Response, "_",
                                           seu@meta.data$orig.ident)
seu@meta.data$treatment_response_samp <- gsub("Untreated_Untreated", "Untreated",
                                         seu@meta.data$treatment_response_samp)

## how many cells in each cluster per sample
# pdf(paste0("Rel_Cell_Dist_Samp_", proj_name, ".pdf"), width = 12, height = 12)
# as.data.frame(table(seu@active.ident, seu@meta.data$orig.ident)) %>%
#   dplyr::rename(cell_type = Var1, sample = Var2, count = Freq) %>%
#   tidyr::spread(key = cell_type, value = count) %>% 
#   dplyr::mutate(total_cells = rowSums(across(2:tidyselect::last_col()))) %>%
#   tidyr::gather(key = cell_type, value = count, - c(total_cells, sample)) %>%  
#   dplyr::mutate(relative_cells = ((count / total_cells) * 100)) %>% 
#   ggplot2::ggplot(aes(x = sample, y = relative_cells, fill = sample)) +
#   geom_bar(stat = "identity", show.legend = FALSE) +
#   ylab("Relative numbers of cells in %") + # theme_classic() +
#   theme(axis.text.x = element_text(angle = 90)) + facet_grid(cell_type ~ .) +
#   theme(strip.text.y.right = element_text(angle = 0)) +
#   ggtitle(paste0("Relative Cell Distributions by Sample for ", proj_name)) -> p1
# print(p1)
# dev.off()
# 
# ## how many cells in each cluster per treatment
# pdf(paste0("Rel_Cell_Dist_Treat_", proj_name, ".pdf"), width = 12, height = 12)
# as.data.frame(table(seu@active.ident, seu@meta.data$Treatment)) %>%
#   dplyr::rename(cell_type = Var1, treatment = Var2, count = Freq) %>%
#   tidyr::spread(key = cell_type, value = count) %>% 
#   dplyr::mutate(total_cells = rowSums(across(2:tidyselect::last_col()))) %>%
#   tidyr::gather(key = cell_type, value = count, - c(total_cells, treatment)) %>% 
#   dplyr::mutate(relative_cells = ((count / total_cells) * 100)) %>% 
#   ggplot2::ggplot(aes(x = treatment, y = relative_cells, fill = treatment)) +
#   geom_bar(stat = "identity", show.legend = FALSE) +
#   ylab("Relative numbers of cells in %") +  # theme_classic() +
#   theme(axis.text.x = element_text(angle = 90)) + facet_grid(cell_type ~ .) +
#   theme(strip.text.y.right = element_text(angle = 0)) +
#   ggtitle(paste0("Relative Cell Distributions by Treatment for ", proj_name)) -> p2
# print(p2)
# dev.off()
# 
# ## how many cells in each cluster per response
# pdf(paste0("Rel_Cell_Dist_Resp_", proj_name, ".pdf"), width = 12, height = 12)
# as.data.frame(table(seu@active.ident, seu@meta.data$Response)) %>%
#   dplyr::rename(cell_type = Var1, response = Var2, count = Freq) %>%
#   tidyr::spread(key = cell_type, value = count) %>% 
#   dplyr::mutate(total_cells = rowSums(across(2:tidyselect::last_col()))) %>%
#   tidyr::gather(key = cell_type, value = count, - c(total_cells, response)) %>% 
#   dplyr::mutate(relative_cells = ((count / total_cells) * 100)) %>% 
#   ggplot2::ggplot(aes(x = response, y = relative_cells, fill = response)) +
#   geom_bar(stat = "identity", show.legend = FALSE) +
#   ylab("Relative numbers of cells in %") + # theme_classic() +
#   theme(axis.text.x = element_text(angle = 90)) + facet_grid(cell_type ~ .) +
#   theme(strip.text.y.right = element_text(angle = 0)) +
#   ggtitle(paste0("Relative Cell Distributions by Response for ", proj_name)) -> p3
# print(p3)
# dev.off()

## how many cells in each cluster per treatment-response
pdf(paste0("Rel_Cell_Dist_TreatResp_", proj_name, ".pdf"), width = 12, height = 12)
as.data.frame(table(seu@active.ident, seu@meta.data$treatment_response)) %>%
  dplyr::rename(cell_type = Var1, treatment_response = Var2, count = Freq) %>%
  # dplyr::filter(cell_type == c("M2 Macrophages")) %>%
  tidyr::spread(key = cell_type, value = count) %>% 
  dplyr::mutate(total_cells = rowSums(across(2:tidyselect::last_col()))) %>%
  tidyr::gather(key = cell_type, value = count, - c(total_cells, treatment_response)) %>% 
  dplyr::mutate(relative_cells = ((count / total_cells) * 100))  %>% 
  ggplot2::ggplot(aes(x = treatment_response, y = relative_cells, fill = treatment_response)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  ylab("Relative numbers of cells in %") + 
  xlab("Treatment + Response") + # theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) + facet_grid(cell_type ~ .) +
  theme(strip.text.y.right = element_text(angle = 0)) +
  ggtitle(paste0("Relative Cell Distributions by Treatment + Response for ", proj_name)) -> p4
print(p4)
dev.off()

## the same as above but with sample names 
pdf(paste0("Rel_Cell_Dist_TreatRespSamp_", proj_name, ".pdf"), width = 12, height = 12)
as.data.frame(table(seu@active.ident, seu@meta.data$treatment_response_samp)) %>%
  dplyr::rename(cell_type = Var1, treatment_response_samp = Var2, count = Freq) %>%
  # dplyr::filter(cell_type != c("M2 Macrophages", "M1 Macrophages")) %>%
  tidyr::spread(key = cell_type, value = count) %>% 
  dplyr::mutate(total_cells = rowSums(across(2:tidyselect::last_col()))) %>%
  tidyr::gather(key = cell_type, value = count, - c(total_cells, treatment_response_samp)) %>% 
  dplyr::mutate(relative_cells = ((count / total_cells) * 100))  %>% 
  ggplot2::ggplot(aes(x = treatment_response_samp, y = relative_cells, fill = treatment_response_samp)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  ylab("Relative numbers of cells in %") + 
  xlab("Treatment + Response") + # theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) + facet_grid(cell_type ~ .) +
  theme(strip.text.y.right = element_text(angle = 0)) +
  ggtitle(paste0("Relative Cell Distributions by Sample Displaying Treatment + Response for ",
                 proj_name)) -> p5
print(p5)
dev.off()

# seu_c <- subset(seu, subset = final_primary_annotations == "Cancer")
# # s1 <- subset(seu_c, subset = Response == "Untreated")
# # saveRDS(s1, file = "Untreated_Cancer.rds")
# s2 <- subset(seu_c, subset = Response == "Relapsing")
# saveRDS(s2, file = "Relapsing_Cancer.rds")
# s3 <- subset(seu_c, subset = Response == "Responding")
# saveRDS(s3, file = "Responding_Cancer.rds")
##