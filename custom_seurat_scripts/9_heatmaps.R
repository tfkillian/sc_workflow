## the script generates heatmaps from Seurat objects

## load libraries
library("dplyr")
library("Seurat")
library("ggplot2")

## marker genes to make a heatmap of

## list of marker genes to plot on y axis
# c(#"Epcam", "Krt19", "Ins1", "Ins2", "Cela2a",	"Cpa1",	"Lars2", "Ppy", #cancer
#   "Mki67", ## proliferating marker
#   "Prox1", "Lyve1", "Mmrn1", "Flt4", #lymphatic # "Dcn", "Pdpn",
#   "Cxcl10",	"Cxcl9", "Gbp2",	"Gbp5",	"Igtp", # Interferon
#   "Apoe",	"Sele", # Postcapillary venule
#   "Ackr1", "Selp", "Madcam1", "Lrg1",  # "Chst4", # HEV markers
#   "Emcn", "Vwf", "Bgn", "Nr2f2", "Mgp",	"Ctsh",	"Lbp", "Mmp2", # vein
#   "Rgs5", "Batf3",	"Tbx2", "Bcl6", "Des", "Ng2", # pericytes
#   "Fbln5",	"Gja5", "Bmx", "Serpinf1", "Sox17", "Stmn2", #artery
#   "Col4a1",	"Col4a2",	"Sema6d", # Stalk cell
#   "Esm1",	"Nid1", "Nid2", "Trp53i11",	"Kcne3", "Apln", "Cd34", "Pdgfb", "Kcne3",
#   "Dll4", ## Tip cells
#   "Ca4", "Plvap",	"Igfbp7",	"Gja5",	"Cxcl12", "Rgcc", "Kdr", "Cd300lg",
#   "Ramp3", #capillaries "Cd36",
#   "Gja4",	"Cxcl12",	"Efnb2" # Cap Arterial
# 
# #   Endothelial tip cells	ESM1	NID1	Trp53i11	Kcne3	Apln	Pgf	Adm	Dll4	Cd34	Pdgfb
# # HEVs	ACKR1	SELP	Madcam1	Lrg1	Chst4
# # Endothelial venous	ACKR1	SELP	Nr2f2	Emcn	Vwf	Bgn	Mgp	Ctsh	Lbp	Mmp2
# # Endothelial arterials	FBLN5	GJA5	Sox17	Stmn2	Bmx	Serpinf1
# # Endothelial capillaries	CD36	CA4	PLVAP	IGFBP7	GJA5	CXCL12
# # Cap Arterial	Gja4	Cxcl12	Efnb2
# # Stalk cell	Col4a1	Col4a2	Sema6d
# # Interferon	Cxcl10	Cxcl9	Gbp2	Gbp5	Igtp
# # Endothelial lymphatics	PROX1	PDPN	Mmrn1	Dcn	Lyve1
# # Postcapillary venule	Selp	Ackr1	Lrg1	Apoe	Sele
# ) -> markers_e

# c(# "Epcam", "Krt19", "Ins1", "Ins2", "Cela2a",	"Cpa1",	"Lars2", "Ppy", #cancer
#   "Mki67", ## proliferating marker
#   "Itga2", "Cd226", "Fcgr4", "Klrk1", "Ptprc", "Nkg2d",
#   "Cd6", "Cd3d", "Cd3e", "Sh2d1a", "Trat1", "Cd3g", # TILs
#   "Ncr1",	"Ncam1", "Cd3", "Klrd1", "Klrg1", # NK
#   "Il17a", "Il17f",	"Ccr6", "Rora", #"Klrb1", #"Il22", # Th17 #
#   "Trgc4", "Trgc", "Trdc", "Trgc1", "Tcrg-c",
#   "Thy1",
#   "Il12rb1", "Stat4", # Th1
#   "Gata3",	"Il4", # Th2
#   "Cxcr5", "Tox",	"Slamf6", #Tfh
#   "Foxp3", "Cd25", "Ikzf2",	"Il2ra", "Batf", "Maf", #Tregs
#   "Gnly",	"Ifng",	"Nkg7", "Prf1",	"Gzma", "Gzmb",	"Gzmh",	"Gzmk", # Cytotoxic CD8+
#   "Havcr2",	"Lag3",	"Pdcd1", "Ctla4", "Tigit", #"Cxcl13",
#   "Btla", "Klrc1", # Inhibitory marker
#   "Anxa1", "Ankrd28",	"Il7r",	"Cd69",	"Cd40lg", # #CD4 effector/Memory effector
#   "Pd1", "Tim3", "Cd244", # exhausted CD8
#   "Ccr7",	"Lef1",	"Sell",	"Tcf1", # naive marker
#   "Cd4", # CD4
#   "Cd8a",	"Cd8b" # CD8
#   ) -> markers_t

# c(# "Epcam", "Krt19", "Ins1", "Ins2", "Cela2a",	"Cpa1",	"Lars2", "Ppy", #cancer
# "Mki67",
# "Pdgfrb", "Gata2", "Ms4a2", "Itgae", "Arg1","Cd84", "Cd226", "F13a1", "F5",
# "Stat1",
# "Tnf", "Flt3",  "Ccl18", "Ifitm1", "Ccl3",
# "Siglech", "Ccr9",  "Gsr", "Bst2",	"Ly6g",	#Neutrophils
# "Clec9a", "Xcr1", "Ido9", "Cd1c", "Clec10a", "Sirpa", "Fcer1a",
# "Ccr7", "Ccl17", "Ccl19", #(cDC1)
# "Pira2", "Cxcr3", "Irf7", #Plasmacytoid denderocytes (pDCs)
# "Axl",  "Cx3cr", "Cd2", "Cd5", "Cd81", "Siglec1", ## tDCs
# "Tcf4", "Bcl11a", "Irf8", "Spib", "Runx2",
# "Gpnmb", "Kit", "Cpa3", "Cma1", "Ilrl1",	"Il4", #MAST cells
# "Selenop", "Sepp1", "Stab1","CCL13", "Ccl13", "Cd163", #M2 markers
# "Il1b", "Cxcl9", "Cdcl10", "Socs3", "Cd86", "Cd80", #M1 markers
# "Cd14",	"Msr1",	"Mrc1", "Mx1", "Cd40", # Monocytes (in general)
# "S100a8", "S100a9", "S100a12",
# "Fcgr1", "Cxcl10", "C1qa", "C1qb", "C1qc", "Ear2", "Lyz1", "Cd83",
# "Vcan",	"Lilrb5", "Marco", "Adgre1", "Fcer1g", #"Cd40" # Macrophages
# "Ccr2", "Ly6c1", "Cd68", "Fcgr1a" # Myeloids (in general)
# ) -> markers_m

## list of marker genes to plot on y axis
# c(# "Ins1", "Ins2", "Cela2a",	"Cpa1",	"Lars2", "Ppy", ## CANCER
#   "Mki67",
#   "Gpr183", "Ccr7", "Rgs13",
#   "Mzb1",	"Cd79a", "Cd79b", "Sdc1", "Ms4a1", 
#   "Prdm1", "Malat1",
#   "Ighm", "Ighd", "Ighg1", "Igha1",
#   "Cd38", "Cd27", "Cd19"
#   ) -> markers_b

#list of marker genes to plot on y axis
c(# "Ins1", "Ins2",  "Iapp", "Insm1",	"Pax6",	"Sez6l", "Cfc1", "Slc8a2",
  # "Ppp1r1a", "Wnt4", "Pdx1", "Gck", "Glut2", ## IT
  # "Enpp2", "Sema3a", "Epha3", "Gpm6a", "Hk1", "Mct1", "Hnf1b", "Gata6", "Lin28b",
  # "Chga", "Scgn", ## MLP
  # "Cela2a", "Cpa1", "Lars2", "Prss2", "Try5", "Ctrb1", ## acinar  # "Cela3a",
  # "Lgr5", "Sox9", "Foxj1", "Ehf", ## stem
  # "Clu", "Krt18", "Krt8", "Krt19", "Krt7", "Epcam", ## ductal # "Tff1",
  # "Pecam1",
  # "Ins1", "Ins2", "Iapp", "Glut2", "Gck", "Isl1", "Insm1", "Ppp1r1a", "Slc2a2",
  # "Slc8a2", "Cfc1", "Sez6l", "Nkx6-1", "Nkx6-2", "Pdx1", "Wnt4", "Pax6", ## IT
  # "Enpp2", "Chga", "Scgn", "Slc16a1", "Slc2a3", "Hk1", "Gmp6a", "Mct1", "Hnf1b",
  # "Gata6", "Lin28b", "Sema3e", "Epha3", "Pou3f4", ## MLP?
  # "Mki67", "Vim", "Alt", "Lamc2", "Sna", "Twist1", "Licam", ## invasion
  # "Cela2a", "Cpa1", "Lars2", "Prss2", "Try5", "Ctrb1", ## acinar
  # "Lgr5", "Sox9", "Foxj1", "Ehf", ## stem
  # "Clu", "Krt18", "Krt8", "Krt19", "Krt7", "Epcam" ## ductal
"MARCO", "TLR2", "TLR4", "CD80", "CD86", "TNF", "IL1B", "IL6",
                   "CSF2", "CXCL2", "IFNG", "IL1R1", ##M1 NOS1",
                   "PPARG", "CLEC10A", "CLEC7A", "PDCD1LG2",  ##M2 "MCR1", "ARG1", "MHCII",
                   "CCL22", "CD40", "IL10",  "IRF4", "PDGFB", "STAT6",
                   "CD68", "CCR5", "TFRC", "ITGAM", "FCGR1A", "CSF1R", "MRC1", ## generally
                   "CD163", "PTPRC", "CD14", "FCN1", "LYZ", "HLA-DRA", "HLA-DRB1",
                   "S100A4", "CD74", "FTH1", "FTL", "CD302", "FCGR3A", "ITGAX"
  ) -> markers_c

# c(# "Epcam", "Krt19", "Ins1", "Ins2", "Cela2a",	"Cpa1",	"Lars2", "Ppy", #cancer
#   # "Pdgfrb", "Gata2", "Ms4a2", "Itgae", "Arg1","Cd84", "Cd226", "F13a1", "F5",
#   "Ly6g", "Mki67",
#   "Tnf", "Ccl18", "Ifitm1", "Ccl3",
#   "Fcgr1", "Cxcl10", "C1qa", "C1qb", "C1qc", "Cd83",
#   "Cxcl1", "Cxcl4", "Cxcl10", "Cxcl11", "Ccl11",
#   "Il1a", "Il23a",  #Should be up in responding tumors:
#   "Tgfb1", "Nsfl1c", #Should be down in responding tumors:
#   "Gsr", "Bst2", #	"Ly6g",	#Neutrophils
#   "Ido9", "Cd1c", "Sirpa", "Fcer1a",
#   "Ccl19", #(cDC1)
#   "Lilra4", "Cxcr3", "Irf7", #Plasmacytoid denderocytes (pDCs)
#   "Kit", "Ilrl1",	"Il4", #MAST cells
#   "Sepp1", "Ccl13", #M2 markers
#   "Il1b", "Cdcl10", "Socs3", "Cd86", "Cd80", #M1 markers
#   "Cd14",	"Mrc1", # Monocytes (in general)
#   "S100a8", "S100a9", "S100a12",
#   "Lilrb5", "Adgre1", "Fcer1g", #"Cd40" # Macrophages
#   "Ly6c1", "Cd68", "Fcgr1a" # Myeloids (in general)
#   ) -> markers_n

# c(# "Ighm", "Cd27", "Cd19", "Ccr7", "Ms4a1", ## b-cells
#   "Cd3d", "Cd3e",	"Ncr1",	"Xcl1", "Gzma", ## T-cells in general
#   "Cd4", # CD4 T-cell
#   "Cd8a",	"Cd8b", # CD8 T-cell
#   "Myl9", "Mylk", "Tagln", "Des", "Acta2", ## smooth muscle
#   "Itgam", "Cd68", "Csf1R", "Cx3cr1",
#   "Msr1",	"Mrc1", "Mx1", "Cd40", "Ly6c1", ## Monocytes
#   "Cd163", "Lilrb5", "Marco", "Adgre1", "Fcer1g", ## Macrophages
#   "Pdgfra", "Col1a1", "Col1a2", "Col3a1", "Col6a2", "Lum", "Mmp2", "Col5a2", # fibro
#   "Cldn5", "Ramp2", "Plvap", "Flt1", "Pecam1", ## endo
#   "Lgr5", "Sox9", "Foxj1", "Ehf", ## stem
#   "Clu", "Krt18", "Krt8", "Krt19", ## ductal
#   "Insm1","Pax6", "Ins1", "Ins2", "Iapp", ## cancer
#   "Prss2", "Try5", "Cela2a", "Cpa1", "Cela3a" ## acinar
# ) -> markers_a

## read seurat object
file_path <- "~/Documents/tmp/pancreatic_melanie/"
heat_type <- "Whole_data"
# heat_type <- "Neutrophils"
# heat_type <- "Cancer"
# heat_type <- "Myeloids"
#heat_type <- "Mast cells & pDCs"
# heat_type <- "DCs"
# heat_type <- "Endothelial"
# heat_type <- "T-Cells"
# file_name <- "merged_mm.rds"
# seu <- readRDS(file = paste0(file_path, file_name))
# seu <- readRDS(file = paste0(file_path, "FINAL_latest_ann_endo.rds"))

## filter Seurat object to contain genes of interest
# seu@meta.data$Cluster
# seu <- SetIdent(object = seu, value = 'RNA_snn_res.1.4')
# seu <- StashIdent(seu, save.name = "Cluster")
seu <- StashIdent(seu, save.name = "cluster")
#seu_counts <- GetAssay(object = seu, assay = "RNA")

seu_genes <- subset(seu, features = markers_c)
seu_counts <- GetAssay(object = seu_genes, assay = "RNA")
seu_mat <- as.matrix(seu_counts@counts)
seu_df <- as.data.frame(seu_mat)
seu@meta.data %>%
  #dplyr::select(orig.ident, primary_annotations_v3) %>%
  dplyr::select(orig.ident, cluster) %>%
  tibble::rownames_to_column() %>%
  dplyr::rename(cells = rowname,
                #Cluster = primary_annotations_v3) -> seu_metadata ## %>% ########
                Cluster = cluster) -> seu_metadata # %>% # 
  # dplyr::filter(Cluster %in% c(0, 1, 2, 3, 4, 5, 6, 10)) -> seu_metadata
  # Cluster = cluster) -> seu_metadata
  # seu_metadata %>%
  #   dplyr::filter(Cluster %in% c(
  #   "Cap arterial", "Cap Esm1+", "Cap Cd36+", "Tip cells Apln+",
  #   "Cap Plpp3+", "Cap Plpp1+", "Tip cells Prcp+")) %>%
  #   dplyr::filter(!is.na(Cluster)) -> seu_metadata

## turn the filtered Seurat counts into a long format tibble and make heatmap
seu_df %>%
  tibble::rownames_to_column() %>%
  dplyr::rename(genes = rowname) %>%
  tidyr::gather(key = "cells", value = "counts", -genes) %>%
  dplyr::left_join(seu_metadata, by = "cells") %>%
  #dplyr::filter(!is.na(Cluster)) %>%
  dplyr::mutate(genes = factor(genes, levels = unique(genes))) %>%
  # dplyr::mutate(Cluster = as.character(Cluster)) %>%
  # dplyr::mutate(Cluster = factor(Cluster, levels = unique(Cluster))) %>%
  dplyr::group_by(Cluster, genes) %>%
  # dplyr::group_by(cluster, genes) %>%
  dplyr::summarise(avg_counts = mean(counts, na.rm = TRUE))  %>%
  dplyr::group_by(genes) %>%
  dplyr::mutate(zscore = (avg_counts - mean(avg_counts)) / sd(avg_counts)) %>%
  ggplot(aes(x = Cluster, y = genes, fill = zscore)) +
  # ggplot(aes(x = cluster, y = genes, fill = zscore)) +
         geom_tile() +
         theme_classic() +
         theme(axis.text.x = element_text(angle = 90),
               axis.text.y = element_text(size = 8)) +
         scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
         ggtitle(paste0("Heatmap of Marker Genes for ",
                        gsub("_", " ", heat_type))) -> p1

## print pdf
pdf(paste0("Heatmap_", heat_type, "_", proj_name, ".pdf"), width = 8.5, height = 11)
print(p1)
dev.off()

############################### heatmap3 #######################################
# library("heatmap3")
# 
# seu <- SetIdent(object = seu, value = RESOLUTION)
# seu <- StashIdent(seu, save.name = "cluster")
# seu_genes <- subset(seu, features = markers_c)
# seu_counts <- GetAssay(object = seu_genes, assay = "RNA")
# seu_mat <- as.matrix(seu_counts@counts)
# seu_df <- as.data.frame(seu_mat)
# seu@meta.data %>%
#   dplyr::select(orig.ident, Cluster) %>%
#   # dplyr::select(orig.ident, cluster) %>%
#   tibble::rownames_to_column() %>%
#   dplyr::rename(cells = rowname) -> seu_metadata
# 
# ## turn the filtered Seurat counts into a long format tibble and make heatmap
# seu_df %>%
#   tibble::rownames_to_column() %>%
#   dplyr::rename(genes = rowname) %>%
#   tidyr::gather(key = "cells", value = "counts", -genes) %>%
#   dplyr::left_join(seu_metadata, by = "cells") %>%
#   dplyr::mutate(genes = factor(genes, levels = unique(genes))) %>%
#   dplyr::mutate(Cluster = as.character(Cluster)) %>%
#   dplyr::group_by(Cluster, genes) %>%
#   # dplyr::group_by(cluster, genes) %>%
#   dplyr::summarise(avg_counts = mean(counts, na.rm = TRUE)) %>% ungroup() %>% 
#   tidyr::spread(key = "cluster", value = "avg_counts", fill = NA) -> g1
# 
# pdf(paste0("Heatmap_dendro_", heat_type, "_", proj_name, ".pdf"),
#     width = 8.5, height = 11)
# g1 %>% dplyr::select(-genes) %>% as.data.frame() -> h1
# rownames(h1) <- g1$genes
# h1 %>% as.matrix() %>%
#   heatmap3::heatmap3(keep.dendro = TRUE, cexRow = 0.5, cexCol = 0.8, Rowv = NA,
#                      main = paste0(gsub("_", " ", heat_type)))
# dev.off()
# 
# ############################### dotplot ########################################
# seu %>% 
#   DotPlot(features = c("Epcam", "Krt19", "Ins1", "Ins2", "Cela2a",	"Cpa1",
#                        "Lars2")) -> p2
# 
# ## print pdf
# pdf(paste0("dotplot_", proj_name, ".pdf"))
# print(p2)
# dev.off()
# 
# ############################### ridgeplot ######################################
# seu %>% 
#   RidgePlot(features = c("Ins1")) +
#   RotatedAxis() +
#   theme(legend.position = "none") -> p3
# 
# ## print pdf
# pdf(paste0("ridgeplot_", proj_name, ".pdf"))
# print(p3)
# dev.off()
# 
# ################################################################################
# 
# ## filter Seurat object to contain genes of interest
# # seu@meta.data$Cluster
# # seu <- StashIdent(seu, save.name = "Cluster")
# seu_ins <- subset(seu, features = c("Epcam", "Krt19", "Ins1", "Ins2", "Cela2a",
#                                     "Cpa1",	"Lars2"))
# seu_counts <- GetAssay(object = seu_ins, assay = "RNA")
# seu_mat <- as.matrix(seu_counts@counts)
# seu_df <- as.data.frame(seu_mat)
# seu@meta.data %>%
#   #dplyr::select(orig.ident, Cluster) %>%
#   dplyr::select(orig.ident, RNA_snn_res.0.4) %>%
#   tibble::rownames_to_column() %>%
#   dplyr::rename(cells = rowname,
#                 Cluster = RNA_snn_res.0.4) -> seu_metadata
# 
# # c("Ins1") -> markers
# 
# ## turn the filtered Seurat counts into a long format tibble and make heatmap
# seu_df %>%
#   # dplyr::filter(row.names(seu_df) %in% markers) %>%
#   tibble::rownames_to_column() %>%
#   dplyr::rename(genes = rowname) %>%
#   tidyr::gather(key = "cells", value = "counts", -genes) %>%
#   dplyr::left_join(seu_metadata, by = "cells") %>%
#   # dplyr::mutate(genes = factor(genes, levels = unique(genes)))  %>%
#   dplyr::group_by(Cluster, genes) %>%
#   dplyr::summarise(avg_counts = mean(counts, na.rm = TRUE)) %>%
#   dplyr::arrange(desc(avg_counts)) %>%
#   # dplyr::mutate(Cluster = factor(Cluster, levels = unique(Cluster))) %>%
#   ggplot(mapping = aes(x = reorder(Cluster, avg_counts), y = avg_counts,
#                        fill = Cluster)) +
#          geom_bar(stat = "identity") +
#          theme_classic() +
#          theme(axis.text.x = element_text(angle = 90)) +
#          ggtitle(paste0("Barplot Avg Expression of Cancer Genes by Cluster")) -> p4
# 
# ## print pdf
# pdf(paste0("barplot_Ins1_", proj_name, ".pdf"))
# print(p4)
# dev.off()
# 
# #    Cluster                    genes avg_counts
# #    <fct>                      <fct>      <dbl>
# #  1 Cancer                     Ins1     1046.  
# #  2 Macrophages                Ins1       87.3 
# #  3 Acinar                     Ins1      576.  
# #  4 Cancer_B20_PDL1_CTLA4_Resp Ins1     1088.  
# #  5 Endothelial                Ins1       62.6 
# #  6 T-cells/NK                 Ins1       15.9 
# #  7 B-Cell                     Ins1        5.04
# #  8 Cancer_B20_PDL1_Resp1      Ins1     1267.  
# #  9 Cancer_B20_PDL1_Resp2      Ins1      672.  
# # 10 Monocytes                  Ins1       49.3 
# # 11 Cancer_B20_Pi3Ki_Resp      Ins1      841.  
# # 12 Fibroblasts                Ins1      105.  
# # 13 Ductal/Stem                Ins1      176.  
# # 14 Smooth Muscle              Ins1       36.2 
# 
# 
# ################################# nCount ######################################
# seu_ins <- subset(seu, features = c("Ins1"))
# seu_counts <- GetAssay(object = seu_ins, assay = "RNA")
# seu_mat <- as.matrix(seu_counts@counts)
# seu_df <- as.data.frame(seu_mat)
# seu@meta.data %>%
#   dplyr::select(orig.ident, nCount_RNA, nFeature_RNA, Cluster) %>%
#   tibble::rownames_to_column() %>%
#   dplyr::rename(cells = rowname) -> seu_metadata
# 
# ## turn the filtered Seurat counts into a long format tibble and make heatmap
# seu_df %>%
#   tibble::rownames_to_column() %>%
#   dplyr::rename(genes = rowname) %>%
#   tidyr::gather(key = "cells", value = "counts", -c(genes)) %>%
#   dplyr::left_join(seu_metadata, by = "cells") %>% 
#   dplyr::filter(Cluster == "Macrophages" | 
#                 Cluster == "Endothelial" |
#                 Cluster == "T-cells/NK" | 
#                 Cluster == "B-Cell" |
#                 Cluster == "Monocytes" | 
#                 Cluster == "Fibroblasts" |
#                 Cluster == "Smooth Muscle") %>% 
#   dplyr::filter(counts > 100) %>%
#   dplyr::group_by(Cluster, genes) %>% 
#   # dplyr::summarise(n = n()) -> h1
#   dplyr::summarise(avg_nCounts = mean(nCount_RNA, na.rm = TRUE))  %>% 
#   ggplot(mapping = aes(x = reorder(Cluster, avg_nCounts), y = avg_nCounts,
#                        fill = Cluster)) +
#          geom_bar(stat = "identity") +
#          theme_classic() +
#          theme(axis.text.x = element_text(angle = 90)) +
#          ggtitle(paste0("Average nCount_RNA of Cells Expressing Ins1 by Cluster")) -> p5
# 
# ## print pdf
# pdf(paste0("nCount_Ins1_", proj_name, ".pdf"))
# print(p5)
# dev.off()
# 
# # 1 Macrophages   Ins1       10362.
# # 2 Endothelial   Ins1        8359.
# # 3 T-cells/NK    Ins1        6884 
# # 4 B-Cell        Ins1        8364.
# # 5 Monocytes     Ins1        9692.
# # 6 Fibroblasts   Ins1       11765.
# # 7 Smooth Muscle Ins1        8223.
# 
# ################################### nFeature ###################################
# 
# seu_ins <- subset(seu, features = c("Ins1"))
# seu_counts <- GetAssay(object = seu_ins, assay = "RNA")
# seu_mat <- as.matrix(seu_counts@counts)
# seu_df <- as.data.frame(seu_mat)
# seu@meta.data %>%
#   dplyr::select(orig.ident, nCount_RNA, nFeature_RNA, Cluster) %>%
#   tibble::rownames_to_column() %>%
#   dplyr::rename(cells = rowname) -> seu_metadata
# 
# ## turn the filtered Seurat counts into a long format tibble and make heatmap
# seu_df %>%
#   tibble::rownames_to_column() %>%
#   dplyr::rename(genes = rowname) %>%
#   tidyr::gather(key = "cells", value = "counts", -c(genes)) %>%
#   dplyr::left_join(seu_metadata, by = "cells") %>% 
#   dplyr::filter(Cluster == "Macrophages" | 
#                 Cluster == "Endothelial" |
#                 Cluster == "T-cells/NK" | 
#                 Cluster == "B-Cell" |
#                 Cluster == "Monocytes" | 
#                 Cluster == "Fibroblasts" |
#                 Cluster == "Smooth Muscle") %>% 
#   dplyr::filter(counts > 500) %>%
#   dplyr::group_by(Cluster, genes) %>% 
#   dplyr::summarise(avg_nFeature = mean(nFeature_RNA, na.rm = TRUE))  %>% 
#   ggplot(mapping = aes(x = reorder(Cluster, avg_nFeature), y = avg_nFeature,
#                        fill = Cluster)) +
#          geom_bar(stat = "identity") +
#          theme_classic() +
#          theme(axis.text.x = element_text(angle = 90)) +
#          ggtitle(paste0("Average nFeature_RNA of Cells Expressing Ins1 by Cluster")) -> p6
# 
# ## print pdf
# pdf(paste0("nFeature_Ins1_", proj_name, ".pdf"))
# print(p6)
# dev.off()

############################### old code #######################################
# sig_dims <- 10
# var_feat <- 2000 ## variableFeature threshold

# c("Siglech", "Bst2", "Ly6g", "Adgre1", "Kit", "Cd68", "Cd86", "Cd80", "Cd163",
#   "Mrc1", "Cd14", "Hck", "Cxcr4", "Xcr1", "Epcam", "Krt19", "Cela2a") -> markers

# clust_avgs <- AverageExpression(seu, return.seurat = TRUE)
# pdf(paste0("Heatmap.pdf"))
# DoHeatmap(clust_avgs,
#           features = c("Siglech", "Bst2", "Ly6g", "Adgre1", "Kit", "Cd68",
#                             "Cd86", "Cd80", "Cd163", "Mrc1", "Cd14", "Hck",
#                             "Cxcr4", "Xcr1", "Epcam", "Krt19", "Cela2a"))
# dev.off()

  # dplyr::mutate(avg_count = mean(counts, na.rm = TRUE)) -> a1 # %>%
  # dplyr::ungroup() %>%
  # dplyr::group_by(Cluster, genes) %>%
  # dplyr::mutate(zscore = (counts - mean(counts))/ sd(counts)) %>%
  # dplyr::ungroup() %>%
  # dplyr::group_by(Cluster) -> seu_filt
# seu_filt$zscore[is.nan(seu_filt$zscore)] <- 0

## create new subsetted Seurat object
# seu_new <- CreateSeuratObject(counts = as.matrix(seu_filt),
#                               min.cells = 1,
#                               min.features = 1)
# seu_new@meta.data <- seu_metadata
# seu_new %>%
#   FindVariableFeatures(selection.method = "vst", nfeatures = var_feat) %>%
#   ScaleData() %>% 
#   RunPCA(npcs = sig_dims, verbose = FALSE) -> seu_new

## NOTE! your heatmap will have x-axis labels of whatever the Seurat Idents are
## and you may need to change them to the Cluster labels, like so:
# seu_new <- SetIdent(object = seu_new, value = "RNA_snn_res.0.8")

## Expression Heatmap
# clust_avgs <- AverageExpression(seu_new, return.seurat = TRUE)
# DoHeatmap(seu_new, features = )
#pdf(paste0("Exp_Heatmap_5_", proj_name, ".pdf"))
# DoHeatmap(clust_avgs, features = unlist(TopFeatures(seu[["pca"]], balanced = TRUE)),
#           size = 5, draw.lines = FALSE)
#dev.off()