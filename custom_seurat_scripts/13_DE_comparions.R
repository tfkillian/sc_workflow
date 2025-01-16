############################# DE - major cell types ############################
## this script performs DE on the major cell types on treatment and condition

## load libraries
library("R.utils")
library("dplyr")
library("Seurat")
library("ggplot2")
library("future")
library("cowplot")

## parallelize workflow and set variables
plan("multiprocess", workers = 50) # uses 50 CPU
options(future.globals.maxSize = 14000 * 1024^2) ## 5GB per worker
proj_name <- "rerun" ## name of project
proj_version <- "DE" ## name of version of project
server_path <- "/media/seq-srv-05/vrc/Project/Project_Theo/pancreatic_mice/rerun/"
setwd(server_path)

#################################### CANCER ####################################
cell_type <- "cancer"
readRDS(paste0("/media/seq-srv-05/vrc/Project/Project_Theo/pancreatic_mice/",
               "rerun/cancer1/cancer_actual2.rds")) %>%
  SetIdent(value = "orig.ident") -> c1

## EXPERIMENTAL GROUPS
GRP1 <- grep(paste(pattern = c("CCB081", "CCB082"), ## UNTREATED
                   collapse = "|"), colnames(x = c1), value = TRUE) 
GRP2 <- grep(paste(pattern = c("CCB023", "CCB027", "CCB028"), ## B20_RESPONDING
                   collapse = "|"), colnames(x = c1), value = TRUE)
GRP3 <- grep(paste(pattern = c("CCB024", "CCB026"), ## B20_RELAPSING
                   collapse = "|"), colnames(x = c1), value = TRUE)
GRP4 <- grep(paste(pattern = c("CCB025", "CCB029", "CCB030"), ## B20_aPDL1_RESPONDING
                   collapse = "|"), colnames(x = c1), value = TRUE)
GRP5 <- grep(paste(pattern = c("CCB088", "CCB102"), ## B20_aPDL1_RELAPSING
                   collapse = "|"), colnames(x = c1), value = TRUE)
GRP6 <- grep(paste(pattern = c("CCB055", "CCB080", "CCB089"), ## B20_PiK3i_RESPONDING
                   collapse = "|"), colnames(x = c1), value = TRUE)
GRP7 <- grep(paste(pattern = c("CCB079", "CCB087", "CCB090"), ## B20_PiK3i_RELAPSING
                   collapse = "|"), colnames(x = c1), value = TRUE)
GRP8 <- grep(paste(pattern = c("CCB097", "CCB098", "CCB0101"), ## B20_+_aPDL1_+_Pi3Ki_RESPONDING
                   collapse = "|"), colnames(x = c1), value = TRUE)
GRP9 <- grep(paste(pattern = c("CCB095", "CCB096"), ## B20_+_aPD1_+_aCTLA4_RESPONDING
                   collapse = "|"), colnames(x = c1), value = TRUE)
GRP10 <- grep(paste(pattern = c("CCB103" #, "CCB104"
                                ), ## B20_+_PDL1_+_LTbR_RESPONDING
                    collapse = "|"), colnames(x = c1), value = TRUE)
GRP11 <- grep(paste(pattern = c("CCB023", "CCB027", "CCB028", "CCB025", "CCB029",
                                "CCB030", "CCB055", "CCB080", "CCB089", "CCB097",
                                "CCB098", "CCB0101", "CCB095", "CCB096", "CCB103"
                                ), ## RESPONDING
                    collapse = "|"), colnames(x = c1), value = TRUE)
GRP12 <- grep(paste(pattern = c("CCB024", "CCB026", "CCB088", "CCB102", "CCB079",
                                "CCB087", "CCB090"), ## RELAPSING
                    collapse = "|"), colnames(x = c1), value = TRUE)
GRP13 <- grep(paste(pattern = c("CCB097", "CCB098", "CCB0101", ## B20+PD1+CTLA4 responding 
                                "CCB095", "CCB096", ## B20+PDL1+Pi3Ki responding 
                                "CCB103" ## B20_+_PDL1_+_LTbR_RESPONDING
                                ), collapse = "|"), colnames(x = c1), value = TRUE)

## 1. UT vs B20 responding
# DE <- FindMarkers(c1, ident.1 = GRP1, ident.2 = GRP2, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_untr_vs_B20_resp_in_", cell_type, "_", proj_name,
#                        ".txt"), col.names = NA, sep = "\t")
# 
# ## 2. B20 responding vs B20 relapsing
# DE <- FindMarkers(c1, ident.1 = GRP2, ident.2 = GRP3, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_B20_resp_vs_B20_rel_in_", cell_type, "_", proj_name,
#                        ".txt"), col.names = NA, sep = "\t")
# 
# ## 3. B20 relapsing vs B20+PD1+CTLA4 responding
# DE <- FindMarkers(c1, ident.1 = GRP3, ident.2 = GRP9, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_B20_rel_vs_B20_PDL_CTLA4_resp_in_", cell_type, "_",
#                        proj_name, ".txt"), col.names = NA, sep = "\t")
# 
# ## 4.  B20 relapsing vs B20+PDL1 responding
# DE <- FindMarkers(c1, ident.1 = GRP3, ident.2 = GRP4, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_B20_rel_vs_B20_PDL1_resp_in_", cell_type, "_",
#                        proj_name, ".txt"), col.names = NA, sep = "\t")
# 
# ## 5. B20 relapsing vs B20+Pi3Ki responding
# DE <- FindMarkers(c1, ident.1 = GRP3, ident.2 = GRP6, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_B20_rel_vs_B20_Pi3Ki_resp_in_", cell_type, "_",
#                        proj_name, ".txt"), col.names = NA, sep = "\t")
# 
# ## 6. B20 relapsing vs B20+PDL1+Pi3K responding
# DE <- FindMarkers(c1, ident.1 = GRP3, ident.2 = GRP8, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_B20_rel_vs_B20_PDL1_Pi3Ki_resp_in_", cell_type, "_",
#                        proj_name, ".txt"), col.names = NA, sep = "\t")
# 
# ## 7. B20+PDL1 responding vs B20+PDL1 relapsing
# DE <- FindMarkers(c1, ident.1 = GRP4, ident.2 = GRP5, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_B20_PDL1_resp_vs_B20_PDL1_rel_in_", cell_type, "_",
#                        proj_name, ".txt"), col.names = NA, sep = "\t")
# 
# ## 8. B20+Pi3Ki responding vs B20+Pi3Ki relapsing
# DE <- FindMarkers(c1, ident.1 = GRP6, ident.2 = GRP7, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_B20_PiK3i_resp_vs_B20_PiK3i_rel_in_", cell_type, "_",
#                        proj_name, ".txt"), col.names = NA, sep = "\t")
# 
# ## 9. untreated vs responding 
# DE <- FindMarkers(c1, ident.1 = GRP1, ident.2 = GRP11, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_untreated_vs_responding_in_", cell_type, "_",
#                        proj_name, ".txt"), col.names = NA, sep = "\t")
# 
# ## 10. untreated vs relapsing
# DE <- FindMarkers(c1, ident.1 = GRP1, ident.2 = GRP12, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_untreated_vs_relapsing_in_", cell_type, "_",
#                        proj_name, ".txt"), col.names = NA, sep = "\t")
# 
# ## 11. responding vs relapsing
# DE <- FindMarkers(c1, ident.1 = GRP11, ident.2 = GRP12, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_responding_vs_relapsing_in_", cell_type, "_",
#                        proj_name, ".txt"), col.names = NA, sep = "\t")

## 13. UT vs B20 relapsing
# DE <- FindMarkers(c1, ident.1 = GRP1, ident.2 = GRP3, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_untr_vs_B20_rel_in_", cell_type, "_", proj_name,
#                        ".txt"), col.names = NA, sep = "\t")
# 
# ## 14. UT vs B20+PDL1 responding
# DE <- FindMarkers(c1, ident.1 = GRP1, ident.2 = GRP4, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_untr_vs_B20_PDL1_resp_in_", cell_type, "_", proj_name,
#                        ".txt"), col.names = NA, sep = "\t")
# 
# ## 15. UT vs B20+PDL1 relapsing
# DE <- FindMarkers(c1, ident.1 = GRP1, ident.2 = GRP5, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_untr_vs_B20_PDL1_rel_in_", cell_type, "_", proj_name,
#                        ".txt"), col.names = NA, sep = "\t")
# 
# ## 16. UT vs B20+Pi3Ki responding
# DE <- FindMarkers(c1, ident.1 = GRP1, ident.2 = GRP6, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_untr_vs_B20_Pi3Ki_resp_in_", cell_type, "_", proj_name,
#                        ".txt"), col.names = NA, sep = "\t")
# 
# ## 17. UT vs B20+Pi3Ki relapsing
# DE <- FindMarkers(c1, ident.1 = GRP1, ident.2 = GRP7, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_untr_vs_B20_Pi3Ki_rel_in_", cell_type, "_", proj_name,
#                        ".txt"), col.names = NA, sep = "\t")
# 
# ## 18. UT vs B20+PDL1+Pi3Ki responding
# DE <- FindMarkers(c1, ident.1 = GRP1, ident.2 = GRP8, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_untr_vs_B20_PD1_Pi3Ki_resp_in_", cell_type, "_",
#                        proj_name, ".txt"), col.names = NA, sep = "\t")
# 
# ## 19. UT vs B20+PD1+CTLA4 responding
# DE <- FindMarkers(c1, ident.1 = GRP1, ident.2 = GRP9, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_untr_vs_B20_PD1_CTLA4_resp_in_", cell_type, "_",
#                        proj_name, ".txt"), col.names = NA, sep = "\t")
# 
# ## 20. UT vs B20_+_PDL1_+_LTbR_RESPONDING
# DE <- FindMarkers(c1, ident.1 = GRP1, ident.2 = GRP10, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_untr_vs_B20_PDL1_LTbR_resp_in_", cell_type, "_",
#                        proj_name, ".txt"), col.names = NA, sep = "\t")
# 
# ## 21. UT vs B20+PD1+CTLA4 responding & B20+PDL1+Pi3Ki responding & B20_+_PDL1_+_LTbR_RESPONDING
# DE <- FindMarkers(c1, ident.1 = GRP1, ident.2 = GRP13, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_untr_vs_B20_PDL1_Pi3Ki&CTLA4&LTbR_resp_in_", cell_type,
#                        "_", proj_name, ".txt"), col.names = NA, sep = "\t")

## 22. B20 RELAPSING vs UT 
DE <- FindMarkers(c1, ident.1 = GRP3, ident.2 = GRP1, test.use = "wilcox",
                  only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
                  max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
write.table(DE, paste0("DE_B20_rel_vs_UT_in_", cell_type,
                       "_", proj_name, ".txt"), col.names = NA, sep = "\t")

## 23. B20 aPDL1 RELAPSING vs UT
DE <- FindMarkers(c1, ident.1 = GRP5, ident.2 = GRP1, test.use = "wilcox",
                  only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
                  max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
write.table(DE, paste0("DE_B20_PDL1_rel_vs_UT_in_", cell_type,
                       "_", proj_name, ".txt"), col.names = NA, sep = "\t")

## 23. B20 Pi3KI RELAPSING vs UT
DE <- FindMarkers(c1, ident.1 = GRP7, ident.2 = GRP1, test.use = "wilcox",
                  only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
                  max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
write.table(DE, paste0("DE_B20_Pi3KI_rel_vs_UT_in_", cell_type,
                       "_", proj_name, ".txt"), col.names = NA, sep = "\t")

## 24. B20 RESPONDING vs UT
DE <- FindMarkers(c1, ident.1 = GRP2, ident.2 = GRP1, test.use = "wilcox",
                  only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
                  max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
write.table(DE, paste0("DE_B20_resp_vs_UT_in_", cell_type,
                       "_", proj_name, ".txt"), col.names = NA, sep = "\t")

## 25. B20 aPDL1 RESPONDING vs UT
DE <- FindMarkers(c1, ident.1 = GRP4, ident.2 = GRP1, test.use = "wilcox",
                  only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
                  max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
write.table(DE, paste0("DE_B20_PDL1_resp_vs_UTin_", cell_type,
                       "_", proj_name, ".txt"), col.names = NA, sep = "\t")

## 26. B20_PiK3i_RESPONDING vs UT
DE <- FindMarkers(c1, ident.1 = GRP6, ident.2 = GRP1, test.use = "wilcox",
                  only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
                  max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
write.table(DE, paste0("DE_B20_Pi3KI_resp_vs_UT_in_", cell_type,
                       "_", proj_name, ".txt"), col.names = NA, sep = "\t")

## 27. B20_+_aPDL1_+_Pi3Ki_RESPONDING vs UT
DE <- FindMarkers(c1, ident.1 = GRP8, ident.2 = GRP1, test.use = "wilcox",
                  only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
                  max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
write.table(DE, paste0("DE_B20_PDL1_Pi3KI_resp_vs_UT_in_", cell_type,
                       "_", proj_name, ".txt"), col.names = NA, sep = "\t")

## 28. B20_+_aPD1_+_aCTLA4_RESPONDING vs UT
DE <- FindMarkers(c1, ident.1 = GRP9, ident.2 = GRP1, test.use = "wilcox",
                  only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
                  max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
write.table(DE, paste0("DE_B20_PD1_CTLA4_resp_vs_UT_in_", cell_type,
                       "_", proj_name, ".txt"), col.names = NA, sep = "\t")

## 29. B20_+_PDL1_+_LTbR_RESPONDING vs UT
DE <- FindMarkers(c1, ident.1 = GRP10, ident.2 = GRP1, test.use = "wilcox",
                  only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
                  max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
write.table(DE, paste0("DE_B20_PD1_LTbR_resp_vs_UT_in_", cell_type,
                       "_", proj_name, ".txt"), col.names = NA, sep = "\t")

#################################### MACRO #################################### 
cell_type <- "macro"
readRDS(paste0("/media/seq-srv-05/vrc/Project/Project_Theo/pancreatic_mice/",
               "rerun/macro1/macro_actual2.rds")) %>%
  SetIdent(value = "orig.ident") -> m1

## EXPERIMENTAL GROUPS
GRP1 <- grep(paste(pattern = c("CCB081", "CCB082"), ## UNTREATED
                   collapse = "|"), colnames(x = m1), value = TRUE) 
GRP2 <- grep(paste(pattern = c("CCB023", "CCB027", "CCB028"), ## B20_RESPONDING
                   collapse = "|"), colnames(x = m1), value = TRUE)
GRP3 <- grep(paste(pattern = c("CCB024", "CCB026"), ## B20_RELAPSING
                   collapse = "|"), colnames(x = m1), value = TRUE)
GRP4 <- grep(paste(pattern = c("CCB025", "CCB029", "CCB030"), ## B20_aPDL1_RESPONDING
                   collapse = "|"), colnames(x = m1), value = TRUE)
GRP5 <- grep(paste(pattern = c("CCB088", "CCB102"), ## B20_aPDL1_RELAPSING
                   collapse = "|"), colnames(x = m1), value = TRUE)
GRP6 <- grep(paste(pattern = c("CCB055", "CCB080", "CCB089"), ## B20_PiK3i_RESPONDING
                   collapse = "|"), colnames(x = m1), value = TRUE)
GRP7 <- grep(paste(pattern = c("CCB079", "CCB087", "CCB090"), ## B20_PiK3i_RELAPSING
                   collapse = "|"), colnames(x = m1), value = TRUE)
GRP8 <- grep(paste(pattern = c("CCB097", "CCB098", "CCB0101"), ## B20_+_aPDL1_+_Pi3Ki_RESPONDING
                   collapse = "|"), colnames(x = m1), value = TRUE)
GRP9 <- grep(paste(pattern = c("CCB095", "CCB096"), ## B20_+_aPD1_+_aCTLA4_RESPONDING
                   collapse = "|"), colnames(x = m1), value = TRUE)
GRP10 <- grep(paste(pattern = c("CCB103" #, "CCB104"
                                ), ## B20_+_PDL1_+_LTbR_RESPONDING
                    collapse = "|"), colnames(x = m1), value = TRUE)
GRP11 <- grep(paste(pattern = c("CCB023", "CCB027", "CCB028", "CCB025", "CCB029",
                                "CCB030", "CCB055", "CCB080", "CCB089", "CCB097",
                                "CCB098", "CCB0101", "CCB095", "CCB096", "CCB103"
                                ), ## RESPONDING
                    collapse = "|"), colnames(x = m1), value = TRUE)
GRP12 <- grep(paste(pattern = c("CCB024", "CCB026", "CCB088", "CCB102", "CCB079",
                                "CCB087", "CCB090"), ## RELAPSING
                    collapse = "|"), colnames(x = m1), value = TRUE)
GRP13 <- grep(paste(pattern = c("CCB097", "CCB098", "CCB0101", ## B20+PD1+CTLA4 responding 
                                "CCB095", "CCB096", ## B20+PDL1+Pi3Ki responding 
                                "CCB103" ## B20_+_PDL1_+_LTbR_RESPONDING
                                ), collapse = "|"), colnames(x = m1), value = TRUE)

## 1. UT vs B20 responding
# DE <- FindMarkers(m1, ident.1 = GRP1, ident.2 = GRP2, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_untr_vs_B20_resp_in_", cell_type, "_", proj_name,
#                        ".txt"), col.names = NA, sep = "\t")
# 
# ## 2. B20 responding vs B20 relapsing
# DE <- FindMarkers(m1, ident.1 = GRP2, ident.2 = GRP3, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_B20_resp_vs_B20_rel_in_", cell_type, "_", proj_name,
#                        ".txt"), col.names = NA, sep = "\t")
# 
# ## 3. B20 relapsing vs B20+PD1+CTLA4 responding
# DE <- FindMarkers(m1, ident.1 = GRP3, ident.2 = GRP9, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_B20_rel_vs_B20_PDL_CTLA4_resp_in_", cell_type, "_",
#                        proj_name, ".txt"), col.names = NA, sep = "\t")
# 
# ## 4.  B20 relapsing vs B20+PDL1 responding
# DE <- FindMarkers(m1, ident.1 = GRP3, ident.2 = GRP4, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_B20_rel_vs_B20_PDL1_resp_in_", cell_type, "_",
#                        proj_name, ".txt"), col.names = NA, sep = "\t")
# 
# ## 5. B20 relapsing vs B20+Pi3Ki responding
# DE <- FindMarkers(m1, ident.1 = GRP3, ident.2 = GRP6, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_B20_rel_vs_B20_Pi3Ki_resp_in_", cell_type, "_",
#                        proj_name, ".txt"), col.names = NA, sep = "\t")
# 
# ## 6. B20 relapsing vs B20+PDL1+Pi3K responding
# DE <- FindMarkers(m1, ident.1 = GRP3, ident.2 = GRP8, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_B20_rel_vs_B20_PDL1_Pi3Ki_resp_in_", cell_type, "_",
#                        proj_name, ".txt"), col.names = NA, sep = "\t")
# 
# ## 7. B20+PDL1 responding vs B20+PDL1 relapsing
# DE <- FindMarkers(m1, ident.1 = GRP4, ident.2 = GRP5, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_B20_PDL1_resp_vs_B20_PDL1_rel_in_", cell_type, "_",
#                        proj_name, ".txt"), col.names = NA, sep = "\t")
# 
# ## 8. B20+Pi3Ki responding vs B20+Pi3Ki relapsing
# DE <- FindMarkers(m1, ident.1 = GRP6, ident.2 = GRP7, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_B20_PiK3i_resp_vs_B20_PiK3i_rel_in_", cell_type, "_",
#                        proj_name, ".txt"), col.names = NA, sep = "\t")
# 
# ## 9. untreated vs responding 
# DE <- FindMarkers(m1, ident.1 = GRP1, ident.2 = GRP11, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_untreated_vs_responding_in_", cell_type, "_",
#                        proj_name, ".txt"), col.names = NA, sep = "\t")
# 
# ## 10. untreated vs relapsing
# DE <- FindMarkers(m1, ident.1 = GRP1, ident.2 = GRP12, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_untreated_vs_relapsing_in_", cell_type, "_",
#                        proj_name, ".txt"), col.names = NA, sep = "\t")
# 
# ## 11. responding vs relapsing
# DE <- FindMarkers(m1, ident.1 = GRP11, ident.2 = GRP12, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_responding_vs_relapsing_in_", cell_type, "_",
#                        proj_name, ".txt"), col.names = NA, sep = "\t")
# 
# ## 13. UT vs B20 relapsing
# DE <- FindMarkers(m1, ident.1 = GRP1, ident.2 = GRP3, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_untr_vs_B20_rel_in_", cell_type, "_", proj_name,
#                        ".txt"), col.names = NA, sep = "\t")
# 
# ## 14. UT vs B20+PDL1 responding
# DE <- FindMarkers(m1, ident.1 = GRP1, ident.2 = GRP4, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_untr_vs_B20_PDL1_resp_in_", cell_type, "_", proj_name,
#                        ".txt"), col.names = NA, sep = "\t")
# 
# ## 15. UT vs B20+PDL1 relapsing
# DE <- FindMarkers(m1, ident.1 = GRP1, ident.2 = GRP5, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_untr_vs_B20_PDL1_rel_in_", cell_type, "_", proj_name,
#                        ".txt"), col.names = NA, sep = "\t")
# 
# ## 16. UT vs B20+Pi3Ki responding
# DE <- FindMarkers(m1, ident.1 = GRP1, ident.2 = GRP6, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_untr_vs_B20_Pi3Ki_resp_in_", cell_type, "_", proj_name,
#                        ".txt"), col.names = NA, sep = "\t")
# 
# ## 17. UT vs B20+Pi3Ki relapsing
# DE <- FindMarkers(m1, ident.1 = GRP1, ident.2 = GRP7, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_untr_vs_B20_Pi3Ki_rel_in_", cell_type, "_", proj_name,
#                        ".txt"), col.names = NA, sep = "\t")
# 
# ## 18. UT vs B20+PDL1+Pi3Ki responding
# DE <- FindMarkers(m1, ident.1 = GRP1, ident.2 = GRP8, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_untr_vs_B20_PD1_Pi3Ki_resp_in_", cell_type, "_", proj_name,
#                        ".txt"), col.names = NA, sep = "\t")
# 
# ## 19. UT vs B20+PD1+CTLA4 responding
# DE <- FindMarkers(m1, ident.1 = GRP1, ident.2 = GRP9, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_untr_vs_B20_PD1_CTLA4_resp_in_", cell_type, "_",
#                        proj_name, ".txt"), col.names = NA, sep = "\t")
# 
# ## 20. UT vs B20_+_PDL1_+_LTbR_RESPONDING
# DE <- FindMarkers(m1, ident.1 = GRP1, ident.2 = GRP10, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_untr_vs_B20_PDL1_LTbR_resp_in_", cell_type, "_", proj_name,
#                        ".txt"), col.names = NA, sep = "\t")
# 
# ## 21. UT vs B20+PD1+CTLA4 responding & B20+PDL1+Pi3Ki responding & B20_+_PDL1_+_LTbR_RESPONDING
# DE <- FindMarkers(m1, ident.1 = GRP1, ident.2 = GRP13, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_untr_vs_B20_PDL1_Pi3Ki&CTLA4&LTbR_resp_in_", cell_type,
#                        "_", proj_name, ".txt"), col.names = NA, sep = "\t")

## 22. B20 RELAPSING vs UT 
DE <- FindMarkers(m1, ident.1 = GRP3, ident.2 = GRP1, test.use = "wilcox",
                  only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
                  max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
write.table(DE, paste0("DE_B20_rel_vs_UT_in_", cell_type,
                       "_", proj_name, ".txt"), col.names = NA, sep = "\t")

## 23. B20 aPDL1 RELAPSING vs UT
DE <- FindMarkers(m1, ident.1 = GRP5, ident.2 = GRP1, test.use = "wilcox",
                  only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
                  max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
write.table(DE, paste0("DE_B20_PDL1_rel_vs_UT_in_", cell_type,
                       "_", proj_name, ".txt"), col.names = NA, sep = "\t")

## 23. B20 Pi3KI RELAPSING vs UT
DE <- FindMarkers(m1, ident.1 = GRP7, ident.2 = GRP1, test.use = "wilcox",
                  only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
                  max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
write.table(DE, paste0("DE_B20_Pi3KI_rel_vs_UT_in_", cell_type,
                       "_", proj_name, ".txt"), col.names = NA, sep = "\t")

## 24. B20 RESPONDING vs UT
DE <- FindMarkers(m1, ident.1 = GRP2, ident.2 = GRP1, test.use = "wilcox",
                  only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
                  max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
write.table(DE, paste0("DE_B20_resp_vs_UT_in_", cell_type,
                       "_", proj_name, ".txt"), col.names = NA, sep = "\t")

## 25. B20 aPDL1 RESPONDING vs UT
DE <- FindMarkers(m1, ident.1 = GRP4, ident.2 = GRP1, test.use = "wilcox",
                  only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
                  max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
write.table(DE, paste0("DE_B20_PDL1_resp_vs_UTin_", cell_type,
                       "_", proj_name, ".txt"), col.names = NA, sep = "\t")

## 26. B20_PiK3i_RESPONDING vs UT
DE <- FindMarkers(m1, ident.1 = GRP6, ident.2 = GRP1, test.use = "wilcox",
                  only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
                  max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
write.table(DE, paste0("DE_B20_Pi3KI_resp_vs_UT_in_", cell_type,
                       "_", proj_name, ".txt"), col.names = NA, sep = "\t")

## 27. B20_+_aPDL1_+_Pi3Ki_RESPONDING vs UT
DE <- FindMarkers(m1, ident.1 = GRP8, ident.2 = GRP1, test.use = "wilcox",
                  only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
                  max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
write.table(DE, paste0("DE_B20_PDL1_Pi3KI_resp_vs_UT_in_", cell_type,
                       "_", proj_name, ".txt"), col.names = NA, sep = "\t")

## 28. B20_+_aPD1_+_aCTLA4_RESPONDING vs UT
DE <- FindMarkers(m1, ident.1 = GRP9, ident.2 = GRP1, test.use = "wilcox",
                  only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
                  max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
write.table(DE, paste0("DE_B20_PD1_CTLA4_resp_vs_UT_in_", cell_type,
                       "_", proj_name, ".txt"), col.names = NA, sep = "\t")

## 29. B20_+_PDL1_+_LTbR_RESPONDING vs UT
DE <- FindMarkers(m1, ident.1 = GRP10, ident.2 = GRP1, test.use = "wilcox",
                  only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
                  max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
write.table(DE, paste0("DE_B20_PD1_LTbR_resp_vs_UT_in_", cell_type,
                       "_", proj_name, ".txt"), col.names = NA, sep = "\t")

################################### T-CELLS #################################### 
cell_type <- "t-cells"
readRDS(paste0("/media/seq-srv-05/vrc/Project/Project_Theo/pancreatic_mice/",
               "rerun/tcell1/tcell_recluster_actual3.rds")) %>%
  SetIdent(value = "orig.ident") -> t1

## EXPERIMENTAL GROUPS
GRP1 <- grep(paste(pattern = c("CCB081", "CCB082"), ## UNTREATED
                   collapse = "|"), colnames(x = t1), value = TRUE) 
GRP2 <- grep(paste(pattern = c("CCB023", "CCB027", "CCB028"), ## B20_RESPONDING
                   collapse = "|"), colnames(x = t1), value = TRUE)
GRP3 <- grep(paste(pattern = c("CCB024", "CCB026"), ## B20_RELAPSING
                   collapse = "|"), colnames(x = t1), value = TRUE)
GRP4 <- grep(paste(pattern = c("CCB025", "CCB029", "CCB030"), ## B20_aPDL1_RESPONDING
                   collapse = "|"), colnames(x = t1), value = TRUE)
GRP5 <- grep(paste(pattern = c("CCB088", "CCB102"), ## B20_aPDL1_RELAPSING
                   collapse = "|"), colnames(x = t1), value = TRUE)
GRP6 <- grep(paste(pattern = c("CCB055", "CCB080", "CCB089"), ## B20_PiK3i_RESPONDING
                   collapse = "|"), colnames(x = t1), value = TRUE)
GRP7 <- grep(paste(pattern = c("CCB079", "CCB087", "CCB090"), ## B20_PiK3i_RELAPSING
                   collapse = "|"), colnames(x = t1), value = TRUE)
GRP8 <- grep(paste(pattern = c("CCB097", "CCB098", "CCB0101"), ## B20_+_aPDL1_+_Pi3Ki_RESPONDING
                   collapse = "|"), colnames(x = t1), value = TRUE)
GRP9 <- grep(paste(pattern = c("CCB095", "CCB096"), ## B20_+_aPD1_+_aCTLA4_RESPONDING
                   collapse = "|"), colnames(x = t1), value = TRUE)
GRP10 <- grep(paste(pattern = c("CCB103" #, "CCB104"
                                ), ## B20_+_PDL1_+_LTbR_RESPONDING
                    collapse = "|"), colnames(x = t1), value = TRUE)
GRP11 <- grep(paste(pattern = c("CCB023", "CCB027", "CCB028", "CCB025", "CCB029",
                                "CCB030", "CCB055", "CCB080", "CCB089", "CCB097",
                                "CCB098", "CCB0101", "CCB095", "CCB096", "CCB103"
                                ), ## RESPONDING
                    collapse = "|"), colnames(x = t1), value = TRUE)
GRP12 <- grep(paste(pattern = c("CCB024", "CCB026", "CCB088", "CCB102", "CCB079",
                                "CCB087", "CCB090"), ## RELAPSING
                    collapse = "|"), colnames(x = t1), value = TRUE)
GRP13 <- grep(paste(pattern = c("CCB097", "CCB098", "CCB0101", ## B20+PD1+CTLA4 responding 
                                "CCB095", "CCB096", ## B20+PDL1+Pi3Ki responding 
                                "CCB103" ## B20_+_PDL1_+_LTbR_RESPONDING
                                ), collapse = "|"), colnames(x = t1), value = TRUE)

## 1. UT vs B20 responding
# DE <- FindMarkers(t1, ident.1 = GRP1, ident.2 = GRP2, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_untr_vs_B20_resp_in_", cell_type, "_", proj_name,
#                        ".txt"), col.names = NA, sep = "\t")
# 
# ## 2. B20 responding vs B20 relapsing
# DE <- FindMarkers(t1, ident.1 = GRP2, ident.2 = GRP3, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_B20_resp_vs_B20_rel_in_", cell_type, "_", proj_name,
#                        ".txt"), col.names = NA, sep = "\t")
# 
# ## 3. B20 relapsing vs B20+PD1+CTLA4 responding
# DE <- FindMarkers(t1, ident.1 = GRP3, ident.2 = GRP9, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_B20_rel_vs_B20_PDL_CTLA4_resp_in_", cell_type, "_",
#                        proj_name, ".txt"), col.names = NA, sep = "\t")
# 
# ## 4.  B20 relapsing vs B20+PDL1 responding
# DE <- FindMarkers(t1, ident.1 = GRP3, ident.2 = GRP4, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_B20_rel_vs_B20_PDL1_resp_in_", cell_type, "_",
#                        proj_name, ".txt"), col.names = NA, sep = "\t")
# 
# ## 5. B20 relapsing vs B20+Pi3Ki responding
# DE <- FindMarkers(t1, ident.1 = GRP3, ident.2 = GRP6, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_B20_rel_vs_B20_Pi3Ki_resp_in_", cell_type, "_",
#                        proj_name, ".txt"), col.names = NA, sep = "\t")
# 
# ## 6. B20 relapsing vs B20+PDL1+Pi3K responding
# DE <- FindMarkers(t1, ident.1 = GRP3, ident.2 = GRP8, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_B20_rel_vs_B20_PDL1_Pi3Ki_resp_in_", cell_type, "_",
#                        proj_name, ".txt"), col.names = NA, sep = "\t")
# 
# ## 7. B20+PDL1 responding vs B20+PDL1 relapsing
# DE <- FindMarkers(t1, ident.1 = GRP4, ident.2 = GRP5, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_B20_PDL1_resp_vs_B20_PDL1_rel_in_", cell_type, "_",
#                        proj_name, ".txt"), col.names = NA, sep = "\t")
# 
# ## 8. B20+Pi3Ki responding vs B20+Pi3Ki relapsing
# DE <- FindMarkers(t1, ident.1 = GRP6, ident.2 = GRP7, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_B20_PiK3i_resp_vs_B20_PiK3i_rel_in_", cell_type, "_",
#                        proj_name, ".txt"), col.names = NA, sep = "\t")
# 
# ## 9. untreated vs responding 
# DE <- FindMarkers(t1, ident.1 = GRP1, ident.2 = GRP11, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_untreated_vs_responding_in_", cell_type, "_",
#                        proj_name, ".txt"), col.names = NA, sep = "\t")
# 
# ## 10. untreated vs relapsing
# DE <- FindMarkers(t1, ident.1 = GRP1, ident.2 = GRP12, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_untreated_vs_relapsing_in_", cell_type, "_",
#                        proj_name, ".txt"), col.names = NA, sep = "\t")
# 
# ## 11. responding vs relapsing
# DE <- FindMarkers(t1, ident.1 = GRP11, ident.2 = GRP12, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_responding_vs_relapsing_in_", cell_type, "_",
#                        proj_name, ".txt"), col.names = NA, sep = "\t")
# 
# ## 13. UT vs B20 relapsing
# DE <- FindMarkers(t1, ident.1 = GRP1, ident.2 = GRP3, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_untr_vs_B20_rel_in_", cell_type, "_", proj_name,
#                        ".txt"), col.names = NA, sep = "\t")
# 
# ## 14. UT vs B20+PDL1 responding
# DE <- FindMarkers(t1, ident.1 = GRP1, ident.2 = GRP4, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_untr_vs_B20_PDL1_resp_in_", cell_type, "_", proj_name,
#                        ".txt"), col.names = NA, sep = "\t")
# 
# ## 15. UT vs B20+PDL1 relapsing
# DE <- FindMarkers(t1, ident.1 = GRP1, ident.2 = GRP5, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_untr_vs_B20_PDL1_rel_in_", cell_type, "_", proj_name,
#                        ".txt"), col.names = NA, sep = "\t")
# 
# ## 16. UT vs B20+Pi3Ki responding
# DE <- FindMarkers(t1, ident.1 = GRP1, ident.2 = GRP6, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_untr_vs_B20_Pi3Ki_resp_in_", cell_type, "_", proj_name,
#                        ".txt"), col.names = NA, sep = "\t")
# 
# ## 17. UT vs B20+Pi3Ki relapsing
# DE <- FindMarkers(t1, ident.1 = GRP1, ident.2 = GRP7, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_untr_vs_B20_Pi3Ki_rel_in_", cell_type, "_", proj_name,
#                        ".txt"), col.names = NA, sep = "\t")
# 
# ## 18. UT vs B20+PDL1+Pi3Ki responding
# DE <- FindMarkers(t1, ident.1 = GRP1, ident.2 = GRP8, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_untr_vs_B20_PD1_Pi3Ki_resp_in_", cell_type, "_", proj_name,
#                        ".txt"), col.names = NA, sep = "\t")
# 
# ## 19. UT vs B20+PD1+CTLA4 responding
# DE <- FindMarkers(t1, ident.1 = GRP1, ident.2 = GRP9, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_untr_vs_B20_PD1_CTLA4_resp_in_", cell_type, "_",
#                        proj_name, ".txt"), col.names = NA, sep = "\t")
# 
# ## 20. UT vs B20_+_PDL1_+_LTbR_RESPONDING
# DE <- FindMarkers(t1, ident.1 = GRP1, ident.2 = GRP10, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_untr_vs_B20_PDL1_LTbR_resp_in_", cell_type, "_", proj_name,
#                        ".txt"), col.names = NA, sep = "\t")
# 
# ## 21. UT vs B20+PD1+CTLA4 responding & B20+PDL1+Pi3Ki responding & B20_+_PDL1_+_LTbR_RESPONDING
# DE <- FindMarkers(t1, ident.1 = GRP1, ident.2 = GRP13, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_untr_vs_B20_PDL1_Pi3Ki&CTLA4&LTbR_resp_in_", cell_type,
#                        "_", proj_name, ".txt"), col.names = NA, sep = "\t")

## 22. B20 RELAPSING vs UT 
DE <- FindMarkers(t1, ident.1 = GRP3, ident.2 = GRP1, test.use = "wilcox",
                  only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
                  max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
write.table(DE, paste0("DE_B20_rel_vs_UT_in_", cell_type,
                       "_", proj_name, ".txt"), col.names = NA, sep = "\t")

## 23. B20 aPDL1 RELAPSING vs UT
DE <- FindMarkers(t1, ident.1 = GRP5, ident.2 = GRP1, test.use = "wilcox",
                  only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
                  max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
write.table(DE, paste0("DE_B20_PDL1_rel_vs_UT_in_", cell_type,
                       "_", proj_name, ".txt"), col.names = NA, sep = "\t")

## 23. B20 Pi3KI RELAPSING vs UT
DE <- FindMarkers(t1, ident.1 = GRP7, ident.2 = GRP1, test.use = "wilcox",
                  only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
                  max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
write.table(DE, paste0("DE_B20_Pi3KI_rel_vs_UT_in_", cell_type,
                       "_", proj_name, ".txt"), col.names = NA, sep = "\t")

## 24. B20 RESPONDING vs UT
DE <- FindMarkers(t1, ident.1 = GRP2, ident.2 = GRP1, test.use = "wilcox",
                  only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
                  max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
write.table(DE, paste0("DE_B20_resp_vs_UT_in_", cell_type,
                       "_", proj_name, ".txt"), col.names = NA, sep = "\t")

## 25. B20 aPDL1 RESPONDING vs UT
DE <- FindMarkers(t1, ident.1 = GRP4, ident.2 = GRP1, test.use = "wilcox",
                  only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
                  max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
write.table(DE, paste0("DE_B20_PDL1_resp_vs_UTin_", cell_type,
                       "_", proj_name, ".txt"), col.names = NA, sep = "\t")

## 26. B20_PiK3i_RESPONDING vs UT
DE <- FindMarkers(t1, ident.1 = GRP6, ident.2 = GRP1, test.use = "wilcox",
                  only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
                  max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
write.table(DE, paste0("DE_B20_Pi3KI_resp_vs_UT_in_", cell_type,
                       "_", proj_name, ".txt"), col.names = NA, sep = "\t")

## 27. B20_+_aPDL1_+_Pi3Ki_RESPONDING vs UT
DE <- FindMarkers(t1, ident.1 = GRP8, ident.2 = GRP1, test.use = "wilcox",
                  only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
                  max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
write.table(DE, paste0("DE_B20_PDL1_Pi3KI_resp_vs_UT_in_", cell_type,
                       "_", proj_name, ".txt"), col.names = NA, sep = "\t")

## 28. B20_+_aPD1_+_aCTLA4_RESPONDING vs UT
DE <- FindMarkers(t1, ident.1 = GRP9, ident.2 = GRP1, test.use = "wilcox",
                  only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
                  max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
write.table(DE, paste0("DE_B20_PD1_CTLA4_resp_vs_UT_in_", cell_type,
                       "_", proj_name, ".txt"), col.names = NA, sep = "\t")

## 29. B20_+_PDL1_+_LTbR_RESPONDING vs UT
DE <- FindMarkers(t1, ident.1 = GRP10, ident.2 = GRP1, test.use = "wilcox",
                  only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
                  max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
write.table(DE, paste0("DE_B20_PD1_LTbR_resp_vs_UT_in_", cell_type,
                       "_", proj_name, ".txt"), col.names = NA, sep = "\t")

################################ ENDO CELLS #################################### 
cell_type <- "EC"
readRDS(paste0("/media/seq-srv-05/vrc/Project/Project_Theo/pancreatic_mice/",
               "rerun/endo1/endo_recluster4_actual.rds")) %>%
  SetIdent(value = "orig.ident") -> e1

## EXPERIMENTAL GROUPS
GRP1 <- grep(paste(pattern = c("CCB081", "CCB082"), ## UNTREATED
                   collapse = "|"), colnames(x = e1), value = TRUE) 
GRP2 <- grep(paste(pattern = c("CCB023", "CCB027", "CCB028"), ## B20_RESPONDING
                   collapse = "|"), colnames(x = e1), value = TRUE)
GRP3 <- grep(paste(pattern = c("CCB024", "CCB026"), ## B20_RELAPSING
                   collapse = "|"), colnames(x = e1), value = TRUE)
GRP4 <- grep(paste(pattern = c("CCB025", "CCB029", "CCB030"), ## B20_aPDL1_RESPONDING
                   collapse = "|"), colnames(x = e1), value = TRUE)
GRP5 <- grep(paste(pattern = c("CCB088", "CCB102"), ## B20_aPDL1_RELAPSING
                   collapse = "|"), colnames(x = e1), value = TRUE)
GRP6 <- grep(paste(pattern = c("CCB055", "CCB080", "CCB089"), ## B20_PiK3i_RESPONDING
                   collapse = "|"), colnames(x = e1), value = TRUE)
GRP7 <- grep(paste(pattern = c("CCB079", "CCB087", "CCB090"), ## B20_PiK3i_RELAPSING
                   collapse = "|"), colnames(x = e1), value = TRUE)
GRP8 <- grep(paste(pattern = c("CCB097", "CCB098", "CCB0101"), ## B20_+_aPDL1_+_Pi3Ki_RESPONDING
                   collapse = "|"), colnames(x = e1), value = TRUE)
GRP9 <- grep(paste(pattern = c("CCB095", "CCB096"), ## B20_+_aPD1_+_aCTLA4_RESPONDING
                   collapse = "|"), colnames(x = e1), value = TRUE)
GRP10 <- grep(paste(pattern = c("CCB103" #, "CCB104"
                                ), ## B20_+_PDL1_+_LTbR_RESPONDING
                    collapse = "|"), colnames(x = e1), value = TRUE)
GRP11 <- grep(paste(pattern = c("CCB023", "CCB027", "CCB028", "CCB025", "CCB029",
                                "CCB030", "CCB055", "CCB080", "CCB089", "CCB097",
                                "CCB098", "CCB0101", "CCB095", "CCB096", "CCB103"
                                ), ## RESPONDING
                    collapse = "|"), colnames(x = e1), value = TRUE)
GRP12 <- grep(paste(pattern = c("CCB024", "CCB026", "CCB088", "CCB102", "CCB079",
                                "CCB087", "CCB090"), ## RELAPSING
                    collapse = "|"), colnames(x = e1), value = TRUE)
GRP13 <- grep(paste(pattern = c("CCB097", "CCB098", "CCB0101", ## B20+PD1+CTLA4 responding 
                                "CCB095", "CCB096", ## B20+PDL1+Pi3Ki responding 
                                "CCB103" ## B20_+_PDL1_+_LTbR_RESPONDING
                                ), collapse = "|"), colnames(x = e1), value = TRUE)

## 1. UT vs B20 responding
# DE <- FindMarkers(e1, ident.1 = GRP1, ident.2 = GRP2, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_untr_vs_B20_resp_in_", cell_type, "_", proj_name,
#                        ".txt"), col.names = NA, sep = "\t")
# 
# ## 2. B20 responding vs B20 relapsing
# DE <- FindMarkers(e1, ident.1 = GRP2, ident.2 = GRP3, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_B20_resp_vs_B20_rel_in_", cell_type, "_", proj_name,
#                        ".txt"), col.names = NA, sep = "\t")
# 
# ## 3. B20 relapsing vs B20+PD1+CTLA4 responding
# DE <- FindMarkers(e1, ident.1 = GRP3, ident.2 = GRP9, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_B20_rel_vs_B20_PDL_CTLA4_resp_in_", cell_type, "_",
#                        proj_name, ".txt"), col.names = NA, sep = "\t")
# 
# ## 4.  B20 relapsing vs B20+PDL1 responding
# DE <- FindMarkers(e1, ident.1 = GRP3, ident.2 = GRP4, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_B20_rel_vs_B20_PDL1_resp_in_", cell_type, "_",
#                        proj_name, ".txt"), col.names = NA, sep = "\t")
# 
# ## 5. B20 relapsing vs B20+Pi3Ki responding
# DE <- FindMarkers(e1, ident.1 = GRP3, ident.2 = GRP6, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_B20_rel_vs_B20_Pi3Ki_resp_in_", cell_type, "_",
#                        proj_name, ".txt"), col.names = NA, sep = "\t")
# 
# ## 6. B20 relapsing vs B20+PDL1+Pi3K responding
# DE <- FindMarkers(e1, ident.1 = GRP3, ident.2 = GRP8, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_B20_rel_vs_B20_PDL1_Pi3Ki_resp_in_", cell_type, "_",
#                        proj_name, ".txt"), col.names = NA, sep = "\t")
# 
# ## 7. B20+PDL1 responding vs B20+PDL1 relapsing
# DE <- FindMarkers(e1, ident.1 = GRP4, ident.2 = GRP5, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_B20_PDL1_resp_vs_B20_PDL1_rel_in_", cell_type, "_",
#                        proj_name, ".txt"), col.names = NA, sep = "\t")
# 
# ## 8. B20+Pi3Ki responding vs B20+Pi3Ki relapsing
# DE <- FindMarkers(e1, ident.1 = GRP6, ident.2 = GRP7, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_B20_PiK3i_resp_vs_B20_PiK3i_rel_in_", cell_type, "_",
#                        proj_name, ".txt"), col.names = NA, sep = "\t")
# 
# ## 9. untreated vs responding 
# DE <- FindMarkers(e1, ident.1 = GRP1, ident.2 = GRP11, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_untreated_vs_responding_in_", cell_type, "_",
#                        proj_name, ".txt"), col.names = NA, sep = "\t")
# 
# ## 10. untreated vs relapsing
# DE <- FindMarkers(e1, ident.1 = GRP1, ident.2 = GRP12, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_untreated_vs_relapsing_in_", cell_type, "_",
#                        proj_name, ".txt"), col.names = NA, sep = "\t")
# 
# ## 11. responding vs relapsing
# DE <- FindMarkers(e1, ident.1 = GRP11, ident.2 = GRP12, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_responding_vs_relapsing_in_", cell_type, "_",
#                        proj_name, ".txt"), col.names = NA, sep = "\t")
# 
# ## 13. UT vs B20 relapsing
# DE <- FindMarkers(e1, ident.1 = GRP1, ident.2 = GRP3, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_untr_vs_B20_rel_in_", cell_type, "_", proj_name,
#                        ".txt"), col.names = NA, sep = "\t")
# 
# ## 14. UT vs B20+PDL1 responding
# DE <- FindMarkers(e1, ident.1 = GRP1, ident.2 = GRP4, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_untr_vs_B20_PDL1_resp_in_", cell_type, "_", proj_name,
#                        ".txt"), col.names = NA, sep = "\t")
# 
# ## 15. UT vs B20+PDL1 relapsing
# DE <- FindMarkers(e1, ident.1 = GRP1, ident.2 = GRP5, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_untr_vs_B20_PDL1_rel_in_", cell_type, "_", proj_name,
#                        ".txt"), col.names = NA, sep = "\t")
# 
# ## 16. UT vs B20+Pi3Ki responding
# DE <- FindMarkers(e1, ident.1 = GRP1, ident.2 = GRP6, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_untr_vs_B20_Pi3Ki_resp_in_", cell_type, "_", proj_name,
#                        ".txt"), col.names = NA, sep = "\t")
# 
# ## 17. UT vs B20+Pi3Ki relapsing
# DE <- FindMarkers(e1, ident.1 = GRP1, ident.2 = GRP7, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_untr_vs_B20_Pi3Ki_rel_in_", cell_type, "_", proj_name,
#                        ".txt"), col.names = NA, sep = "\t")
# 
# ## 18. UT vs B20+PDL1+Pi3Ki responding
# DE <- FindMarkers(e1, ident.1 = GRP1, ident.2 = GRP8, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_untr_vs_B20_PD1_Pi3Ki_resp_in_", cell_type, "_",
#                        proj_name, ".txt"), col.names = NA, sep = "\t")
# 
# ## 19. UT vs B20+PD1+CTLA4 responding
# DE <- FindMarkers(e1, ident.1 = GRP1, ident.2 = GRP9, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_untr_vs_B20_PD1_CTLA4_resp_in_", cell_type, "_",
#                        proj_name, ".txt"), col.names = NA, sep = "\t")
# 
# ## 20. UT vs B20_+_PDL1_+_LTbR_RESPONDING
# DE <- FindMarkers(e1, ident.1 = GRP1, ident.2 = GRP10, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_untr_vs_B20_PDL1_LTbR_resp_in_", cell_type, "_", proj_name,
#                        ".txt"), col.names = NA, sep = "\t")
# 
# ## 21. UT vs B20+PD1+CTLA4 responding & B20+PDL1+Pi3Ki responding & B20_+_PDL1_+_LTbR_RESPONDING
# DE <- FindMarkers(e1, ident.1 = GRP1, ident.2 = GRP13, test.use = "wilcox",
#                   only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
#                   max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
# write.table(DE, paste0("DE_untr_vs_B20_PDL1_Pi3Ki&CTLA4&LTbR_resp_in_", cell_type,
#                        "_", proj_name, ".txt"), col.names = NA, sep = "\t")

## 22. B20 RELAPSING vs UT 
DE <- FindMarkers(e1, ident.1 = GRP3, ident.2 = GRP1, test.use = "wilcox",
                  only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
                  max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
write.table(DE, paste0("DE_B20_rel_vs_UT_in_", cell_type,
                       "_", proj_name, ".txt"), col.names = NA, sep = "\t")

## 23. B20 aPDL1 RELAPSING vs UT
DE <- FindMarkers(e1, ident.1 = GRP5, ident.2 = GRP1, test.use = "wilcox",
                  only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
                  max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
write.table(DE, paste0("DE_B20_PDL1_rel_vs_UT_in_", cell_type,
                       "_", proj_name, ".txt"), col.names = NA, sep = "\t")

## 23. B20 Pi3KI RELAPSING vs UT
DE <- FindMarkers(e1, ident.1 = GRP7, ident.2 = GRP1, test.use = "wilcox",
                  only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
                  max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
write.table(DE, paste0("DE_B20_Pi3KI_rel_vs_UT_in_", cell_type,
                       "_", proj_name, ".txt"), col.names = NA, sep = "\t")

## 24. B20 RESPONDING vs UT
DE <- FindMarkers(e1, ident.1 = GRP2, ident.2 = GRP1, test.use = "wilcox",
                  only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
                  max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
write.table(DE, paste0("DE_B20_resp_vs_UT_in_", cell_type,
                       "_", proj_name, ".txt"), col.names = NA, sep = "\t")

## 25. B20 aPDL1 RESPONDING vs UT
DE <- FindMarkers(e1, ident.1 = GRP4, ident.2 = GRP1, test.use = "wilcox",
                  only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
                  max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
write.table(DE, paste0("DE_B20_PDL1_resp_vs_UTin_", cell_type,
                       "_", proj_name, ".txt"), col.names = NA, sep = "\t")

## 26. B20_PiK3i_RESPONDING vs UT
DE <- FindMarkers(e1, ident.1 = GRP6, ident.2 = GRP1, test.use = "wilcox",
                  only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
                  max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
write.table(DE, paste0("DE_B20_Pi3KI_resp_vs_UT_in_", cell_type,
                       "_", proj_name, ".txt"), col.names = NA, sep = "\t")

## 27. B20_+_aPDL1_+_Pi3Ki_RESPONDING vs UT
DE <- FindMarkers(e1, ident.1 = GRP8, ident.2 = GRP1, test.use = "wilcox",
                  only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
                  max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
write.table(DE, paste0("DE_B20_PDL1_Pi3KI_resp_vs_UT_in_", cell_type,
                       "_", proj_name, ".txt"), col.names = NA, sep = "\t")

## 28. B20_+_aPD1_+_aCTLA4_RESPONDING vs UT
DE <- FindMarkers(e1, ident.1 = GRP9, ident.2 = GRP1, test.use = "wilcox",
                  only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
                  max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
write.table(DE, paste0("DE_B20_PD1_CTLA4_resp_vs_UT_in_", cell_type,
                       "_", proj_name, ".txt"), col.names = NA, sep = "\t")

## 29. B20_+_PDL1_+_LTbR_RESPONDING vs UT
DE <- FindMarkers(e1, ident.1 = GRP10, ident.2 = GRP1, test.use = "wilcox",
                  only.pos = FALSE, logfc.threshold = 0, min.cells.group = 0,
                  max.cells.per.ident = Inf, min.cells.feature = 1, min.pct = 0)
write.table(DE, paste0("DE_B20_PD1_LTbR_resp_vs_UT_in_", cell_type,
                       "_", proj_name, ".txt"), col.names = NA, sep = "\t")