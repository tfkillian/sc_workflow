############################### Marker genes ##################################
## Annotate the cluster based on gene expression. To annotate all clusters, run
## all the following lines (with the correct names and genes) then, exclude
## clusters you are sure about (e.g T-cells, B-cells) for each other cluster,
## with the "FindAllMarker", to see the most highly expressed genes for each
## cluster. With literature and knowledge, it is possible to ID all clusters
## normally based on that info, (e.g. for Cluster 2 is Dendritic cells because
## the top 10 expressed genes are marker genes for Dendritic cells). Repeat as
## many times as needed. Marker genes can be found in literature, or use the top
## genes expressed in FindAllMarker table. NOTE: don't forget: CAPITAL GENES:
## HUMAN, lower case: mouse. e.g.: PECAM or Pecam1 are the same gene, but the
## first is human and the second one is mouse. Make pdfs with 4 genes at a time
## so that they show up better in the pdfs

## load libraries
library("dplyr")
library("ggplot2")
library("tidyselect")
library("Seurat")

## evergreen variables
# DefaultAssay(seu) <- def_assay
# proj_name <- "pancreatic_mice"
# proj_version <- "unfiltered"
ggl_a <- list(); ggl_b <- list(); ggl_c <- list(); ggl_d <- list();
ggl_e <- list(); ggl_f <- list(); ggl_g <- list(); ggl_h <- list();
ggl_i <- list(); ggl_j <- list(); ggl_k <- list(); ggl_l <- list();
ggl_m <- list(); ggl_n <- list(); ggl_o <- list(); ggl_p <- list();
ggl_q <- list(); ggl_r <- list(); ggl_s <- list(); ggl_t <- list();
ggl_u <- list(); ggl_v <- list(); ggl_w <- list(); ggl_x <- list();
ggl_y <- list(); ggl_z <- list(); ggl_aa <- list(); ggl_bb <- list();
ggl_cc <- list(); ggl_dd <- list(); ggl_ee <- list(); ggl_ff <- list();
ggl_gg <- list(); ggl_hh <- list(); ggl_ii <- list(); ggl_jj <- list();
ggl_kk <- list(); ggl_ll <- list(); ggl_mm <- list(); ggl_nn <- list();
ggl_oo <- list(); ggl_pp <- list(); ggl_qq <- list(); ggl_rr <- list();
ggl_ss <- list(); ggl_tt <- list(); ggl_uu <- list(); ggl_vv <- list();
ggl_ww <- list(); ggl_xx <- list(); ggl_yy <- list(); ggl_zz <- list();
ggl_aaa <- list(); ggl_bbb <- list(); ggl_ccc <- list(); ggl_ddd <- list();
ggl_eee <- list(); ggl_fff <- list(); ggl_ggg <- list(); ggl_hhh <- list(); 
ggl_iii <- list(); ggl_jjj <- list(); ggl_kkk <- list(); ggl_lll <- list(); 

############################### Endothelial cells ##############################
# endo_list <- list("Pecam1", "Chst4", "Cldn5", "Ramp2", "Plvap", "Flt1", "Cdh5",
#                   "Cd34", "Hspg2", "Igfbp7", "Plvap", "Cd34") 
# for (i in 1:length(endo_list)){
# FeaturePlot(object = seu, features = endo_list[[i]], reduction = "umap",
#             cols = c("azure3", "deeppink3"), pt.size = 0.1) -> ggl_a[[i]]
# VlnPlot(object = seu, features = endo_list[[i]], pt.size = 0) +
#   theme(axis.text.x = element_text(angle = 45, size = 10)) -> ggl_b[[i]]
# }
# pdf(paste0("endothelial_cells_", proj_name, ".pdf"))
# for (i in 1:length(endo_list)){
# print(ggl_a[[i]])
# print(ggl_b[[i]])
# }
# dev.off()
# 
# ############################## Myeloids ########################################
# # myl_list <- list("Cd68", "Fcer1g", "Ear2", "Lyz1")
# myl_list <- list("CD68", "PARP14", "SEPP1")
# for (i in 1:length(myl_list)){
# FeaturePlot(object = seu, features = myl_list[[i]], reduction = "umap",
#             cols = c("azure3", "deeppink3"), pt.size = 0.1) -> ggl_c[[i]]
# VlnPlot(object = seu, features = myl_list[[i]], pt.size = 0) +
#   theme(axis.text.x = element_text(angle = 45, size = 10)) -> ggl_d[[i]]
# }
# pdf(paste0("myeloids_", proj_name, ".pdf"))
# for (i in 1:length(myl_list)){
# print(ggl_c[[i]])
# print(ggl_d[[i]])
# }
# dev.off()
# 
# ############################## B cells 1 #######################################
# bcell_list <- list("Cd79a", "Cd79b", "Sdc1", "Cd27", "Ms4a1", "Cd19", "Ly6d",
#                    "Ms4a1") # "B220"
# for (i in 1:length(bcell_list)){
# FeaturePlot(object = seu, features = bcell_list[[i]], reduction = "umap",
#             cols = c("azure3", "deeppink3"), pt.size = 0.1) -> ggl_e[[i]]
# VlnPlot(object = seu, features = bcell_list[[i]], pt.size = 0) +
#   theme(axis.text.x = element_text(angle = 45, size = 10)) -> ggl_f[[i]]
# }
# pdf(paste0("bcell_", proj_name, ".pdf"))
# for (i in 1:length(bcell_list)){
# print(ggl_e[[i]])
# print(ggl_f[[i]])
# }
# dev.off()
# 
# ################################ Cancer 1 ######################################
# cancer_list <- list(# "Pecam1", "Epcam",
#                     "Ins1", "Ins2", "Ppy", "Iapp"
#                     ) # "Krt6" , "Scgb2a2" "Krt19", "Krt14", "Krt7", 
# for (i in 1:length(cancer_list)){
# FeaturePlot(object = seu, features = cancer_list[[i]], reduction = "umap",
#             cols = c("azure3", "deeppink3"), pt.size = 0.1) -> ggl_g[[i]]
# VlnPlot(object = seu, features = cancer_list[[i]], pt.size = 0) +
#   theme(axis.text.x = element_text(angle = 45, size = 10)) -> ggl_h[[i]]
# }
# pdf(paste0("cancer_", proj_name, ".pdf"))
# for (i in 1:length(cancer_list)){
# print(ggl_g[[i]])
# print(ggl_h[[i]])
# }
# dev.off()

############################### Macrophages ####################################
# macro_list <- list("Cd68", "Adgre1", "Cd163", "Fcer1g", "Cd40", "Itgam", "Fcgr1",
#                    "Cxcl10", "Arg1", "Mrc1", "Cd74", "C1qb", "Tgfb1", "Il10",
#                    "C1qa", "C1qc", "Nos2", "Fpr2", "Il1a", "Cxcl9", "Apoe",
#                    "Saa3")
# macro_list <- list("VCAN", "MACROD2", "CD163", "LILRB5", "S100A8", "S100A12",
#                    "MARCO", "CD40")
macro_list <- list(#"MARCO", "TLR2", "TLR4", "CD80", "CD86", "TNF", "IL1B", "IL6",
                   #"CSF2", "CXCL2", "IFNG", "IL1R1", "NOS1", ##M1
                   #"PPARG", "CLEC10A", "CLEC7A", "PDCD1LG2",  ##M2  "ARG1", 
                   #"CCL22", "CD40", "IL10",  "IRF4", "PDGFB", "STAT6",
                   "CD68", "CCR5", "TFRC", "ITGAM", "FCGR1A", "CSF1R", ## generally
                   "CD163", "PTPRC", "CD14", "FCN1", "LYZ", "HLA-DRA", "HLA-DRB1",
                   "S100A4", "CD74", "FTH1", "FTL", "CD302", "FCGR3A", "ITGAX"
                   ) ## "IL4",
for (i in 1:length(macro_list)){
FeaturePlot(object = seu, features = macro_list[[i]], reduction = "umap",
            cols = c("azure3", "deeppink3"), pt.size = 0.1) -> ggl_i[[i]]
VlnPlot(object = seu, features = macro_list[[i]], pt.size = 0) +
  theme(axis.text.x = element_text(angle = 45, size = 10)) -> ggl_j[[i]]
}
pdf(paste0("macro_", proj_name, ".pdf"))
for (i in 1:length(macro_list)){
print(ggl_i[[i]])
print(ggl_j[[i]])
}
dev.off()

################################# Monocytes ###################################
# mono_list <- list("Cd40", "Hck", "Cxcr4", "Cd86", "Cd14", "Fcgr2b", "Mx1",
#                   "Il1rn", "Ifit1", "Ifit3", "Cxcl10") # "Ifr7"
mono_list <- list("CD14", "S100A8", "S100A9", "CD68", "MSR1", "APOBEC3A", "CFP",
                  "CD7", "TET2", "CD40", "DYSF", "HCK", "SELE", "CXCR4")
for (i in 1:length(mono_list)){
FeaturePlot(object = seu, features = mono_list[[i]], reduction = "umap",
            cols = c("azure3", "deeppink3"), pt.size = 0.1) -> ggl_k[[i]]
VlnPlot(object = seu, features = mono_list[[i]], pt.size = 0) +
  theme(axis.text.x = element_text(angle = 45, size = 10)) -> ggl_l[[i]]
}
pdf(paste0("mono_", proj_name, ".pdf"))
for (i in 1:length(mono_list)){
print(ggl_k[[i]])
print(ggl_l[[i]])
}
dev.off()

############################### Erythrocytes ###################################
# eryth_list <- list("Hba-a1", "Hba-a2", "Hbb-bs", "Hbb-bt")
# for (i in 1:length(eryth_list)){
# FeaturePlot(object = seu, features = eryth_list[[i]], reduction = "umap",
#             cols = c("azure3", "deeppink3"), pt.size = 0.1) -> ggl_m[[i]]
# VlnPlot(object = seu, features = eryth_list[[i]], pt.size = 0) +
#   theme(axis.text.x = element_text(angle = 45, size = 10)) -> ggl_n[[i]]
# }
# pdf(paste0("eryth_", proj_name, ".pdf"))
# for (i in 1:length(eryth_list)){
# print(ggl_m[[i]])
# print(ggl_n[[i]])
# }
# dev.off()
# 
# ################################# Fibroblasts ##################################
# fibro_list <- list("Pdgfra", "Col1a1", "Col1a2", "Col3a1", "Col6a2", "Lum",
#                    "Mmp2", "Col5a2")
# for (i in 1:length(fibro_list)){
# FeaturePlot(object = seu, features = fibro_list[[i]], reduction = "umap",
#             cols = c("azure3", "deeppink3"), pt.size = 0.1) -> ggl_o[[i]]
# VlnPlot(object = seu, features = fibro_list[[i]], pt.size = 0) +
#   theme(axis.text.x = element_text(angle = 45, size = 10)) -> ggl_p[[i]]
# }
# pdf(paste0("fibro_", proj_name, ".pdf"))
# for (i in 1:length(fibro_list)){
# print(ggl_o[[i]])
# print(ggl_p[[i]])
# }
# dev.off()
# 
# ############################### Smooth muscle 1 ################################
# smooth_list <- list("Des", "Acta2", "Myl9", "Mylk", "Tagln", "Pcp4l1", "Myh11",
#                     "Rgs5", "Mustn1", "Mylk", "Lamb2", "Notch3", "Cnn1", "Gja4",
#                     "Map3k7cl", "Lmod1", "Aaed1")
# for (i in 1:length(smooth_list)){
# FeaturePlot(object = seu, features = smooth_list[[i]], reduction = "umap",
#             cols = c("azure3", "deeppink3"), pt.size = 0.1) -> ggl_q[[i]]
# VlnPlot(object = seu, features = smooth_list[[i]], pt.size = 0) +
#   theme(axis.text.x = element_text(angle = 45, size = 10)) -> ggl_r[[i]]
# }
# pdf(paste0("smooth_", proj_name, ".pdf"))
# for (i in 1:length(smooth_list)){
# print(ggl_q[[i]])
# print(ggl_r[[i]])
# }
# dev.off()
# 
# ################################# Neutrophils ################################## 
# neutro_list <- list("Gsr", "S100a9", "Ly6g", "S100a8", "G0s2", "Ncf1", "Cd177",
#                     "Lrg1", "Camp") # Cd45, "Cd11b",
neutro_list <- list("TREM1", "BST1", "CSF3R", "CTSG", "OSM", "MMP9")
for (i in 1:length(neutro_list)){
FeaturePlot(object = seu, features = neutro_list[[i]], reduction = "umap",
            cols = c("azure3", "deeppink3"), pt.size = 0.1) -> ggl_s[[i]]
VlnPlot(object = seu, features = neutro_list[[i]], pt.size = 0) +
  theme(axis.text.x = element_text(angle = 45, size = 10)) -> ggl_t[[i]]
}
pdf(paste0("neutro_", proj_name, ".pdf"))
for (i in 1:length(neutro_list)){
print(ggl_s[[i]])
print(ggl_t[[i]])
}
dev.off()
# 
# ########################## T-cells and NK-cells ################################
# tcell_list <- list("Cd3d", "Cd3e", "Ncr1", "Xcl1", "Cd8a", "Cd8b1", "Cd4",
#                    "Cd3g", "Gzma", "Nkg7")
# for (i in 1:length(tcell_list)){
# FeaturePlot(object = seu, features = tcell_list[[i]], reduction = "umap",
#             cols = c("azure3", "deeppink3"), pt.size = 0.1) -> ggl_u[[i]]
# VlnPlot(object = seu, features = tcell_list[[i]], pt.size = 0) +
#   theme(axis.text.x = element_text(angle = 45, size = 10)) -> ggl_v[[i]]
# }
# pdf(paste0("tcell_", proj_name, ".pdf"))
# for (i in 1:length(tcell_list)){
# print(ggl_u[[i]])
# print(ggl_v[[i]])
# }
# dev.off()
# 
# ############################# Acinar cells ##################################### 
# acinar_list <- list("Cela2a", "Cpa1", "Cela3a", "Lars2", "Prss2", "Try5", "Ctrb1")
# for (i in 1:length(acinar_list)){
# FeaturePlot(object = seu, features = acinar_list[[i]], reduction = "umap",
#             cols = c("azure3", "deeppink3"), pt.size = 0.1) -> ggl_w[[i]]
# VlnPlot(object = seu, features = acinar_list[[i]], pt.size = 0) +
#   theme(axis.text.x = element_text(angle = 45, size = 10)) -> ggl_x[[i]]
# }
# pdf(paste0("acinar_", proj_name, ".pdf"))
# for (i in 1:length(acinar_list)){
# print(ggl_w[[i]])
# print(ggl_x[[i]])
# }
# dev.off()
# 
# ############################### Perivascular Cells #############################
# peri_list <- list("Rgs5", "Acta2")
# for (i in 1:length(peri_list)){
# FeaturePlot(object = seu, features = peri_list[[i]], reduction = "umap",
#             cols = c("azure3", "deeppink3"), pt.size = 0.1) -> ggl_y[[i]]
# VlnPlot(object = seu, features = peri_list[[i]], pt.size = 0) +
#   theme(axis.text.x = element_text(angle = 45, size = 10)) -> ggl_z[[i]]
# }
# pdf(paste0("peri_", proj_name, ".pdf"))
# for (i in 1:length(peri_list)){
# print(ggl_y[[i]])
# print(ggl_z[[i]])
# }
# dev.off()
# 
# ################ Epithelial-mesenchymal transition (EMT) Cells #################
# emt_list <- list("Cdkn2a", "S100a6", "Igfbp4", "Vim", "Spp1")
# for (i in 1:length(emt_list)){
# FeaturePlot(object = seu, features = emt_list[[i]], reduction = "umap",
#             cols = c("azure3", "deeppink3"), pt.size = 0.1) -> ggl_aa[[i]]
# VlnPlot(object = seu, features = emt_list[[i]], pt.size = 0) +
#   theme(axis.text.x = element_text(angle = 45, size = 10)) -> ggl_bb[[i]]
# }
# pdf(paste0("emt_", proj_name, ".pdf"))
# for (i in 1:length(emt_list)){
# print(ggl_aa[[i]])
# print(ggl_bb[[i]])
# }
# dev.off()
# 
# ############################ Dendritic Cells ###################################
# dend_list <- list("Ccl5", "Ccr7", "Ccl22")
# for (i in 1:length(dend_list)){
# FeaturePlot(object = seu, features = dend_list[[i]], reduction = "umap",
#             cols = c("azure3", "deeppink3"), pt.size = 0.1) -> ggl_cc[[i]]
# VlnPlot(object = seu, features = dend_list[[i]], pt.size = 0) +
#   theme(axis.text.x = element_text(angle = 45, size = 10)) -> ggl_dd[[i]]
# }
# pdf(paste0("dend_", proj_name, ".pdf"))
# for (i in 1:length(dend_list)){
# print(ggl_cc[[i]])
# print(ggl_dd[[i]])
# }
# dev.off()
# 
# ############################## Ductal Cells ####################################
# duct_list <- list("Clu", "Tff1", "Krt18", "Krt8", "Krt19", "Krt7", "Epcam")
# for (i in 1:length(duct_list)){
# FeaturePlot(object = seu, features = duct_list[[i]], reduction = "umap",
#             cols = c("azure3", "deeppink3"), pt.size = 0.1) -> ggl_ee[[i]]
# VlnPlot(object = seu, features = duct_list[[i]], pt.size = 0) +
#   theme(axis.text.x = element_text(angle = 45, size = 10)) -> ggl_ff[[i]]
# }
# pdf(paste0("duct_", proj_name, ".pdf"))
# for (i in 1:length(duct_list)){
# print(ggl_ee[[i]])
# print(ggl_ff[[i]])
# }
# dev.off()
# 
# ################################ Alpha Cells ###################################
# alpha_list <- list("Gcg", #"Pcsk2", "Mafb", "Chga", "Neurod1",
#                    "Ppy" #, "Cryba2"
#                    )
# for (i in 1:length(alpha_list)){
# FeaturePlot(object = seu, features = alpha_list[[i]], reduction = "umap",
#             cols = c("azure3", "deeppink3"), pt.size = 0.1) -> ggl_gg[[i]]
# VlnPlot(object = seu, features = alpha_list[[i]], pt.size = 0) +
#   theme(axis.text.x = element_text(angle = 45, size = 10)) -> ggl_hh[[i]]
# }
# pdf(paste0("alpha_", proj_name, ".pdf"))
# for (i in 1:length(alpha_list)){
# print(ggl_gg[[i]])
# print(ggl_hh[[i]])
# }
# dev.off()
# 
# ################################## Beta Cells ##################################
# beta_list <- list("Ins1", "Ins2", "Ppy" #,"Nkx6-1", "Mafa", "Pdx1", "Syp", "Chga", "Pax6"
#                   )
# for (i in 1:length(beta_list)){
# FeaturePlot(object = seu, features = beta_list[[i]], reduction = "umap",
#             cols = c("azure3", "deeppink3"), pt.size = 0.1) -> ggl_ii[[i]]
# VlnPlot(object = seu, features = beta_list[[i]], pt.size = 0) +
#   theme(axis.text.x = element_text(angle = 45, size = 10)) -> ggl_jj[[i]]
# }
# pdf(paste0("beta_", proj_name, ".pdf"))
# for (i in 1:length(beta_list)){
# print(ggl_ii[[i]])
# print(ggl_jj[[i]])
# }
# dev.off()
# 
# ################################# Delta Cells ##################################
# delta_list <- list("Frzb", "Pcsk1", "Etv1", "Ffar4", "Hhex", "Isl1")
# for (i in 1:length(delta_list)){
# FeaturePlot(object = seu, features = delta_list[[i]], reduction = "umap",
#             cols = c("azure3", "deeppink3"), pt.size = 0.1) -> ggl_kk[[i]]
# VlnPlot(object = seu, features = delta_list[[i]], pt.size = 0) +
#   theme(axis.text.x = element_text(angle = 45, size = 10)) -> ggl_ll[[i]]
# }
# pdf(paste0("delta_", proj_name, ".pdf"))
# for (i in 1:length(delta_list)){
# print(ggl_kk[[i]])
# print(ggl_ll[[i]])
# }
# dev.off()
# 
# ################################## Stellate Cells ##############################
# stell_list <- list("Col1a2", "Mmp14", "Pdgfrb") # "Colga1", "Colga3", "Col3a3", 
# for (i in 1:length(stell_list)){
# FeaturePlot(object = seu, features = stell_list[[i]], reduction = "umap",
#             cols = c("azure3", "deeppink3"), pt.size = 0.1) -> ggl_mm[[i]]
# VlnPlot(object = seu, features = stell_list[[i]], pt.size = 0) +
#   theme(axis.text.x = element_text(angle = 45, size = 10)) -> ggl_nn[[i]]
# }
# pdf(paste0("stellate_", proj_name, ".pdf"))
# for (i in 1:length(stell_list)){
# print(ggl_mm[[i]])
# print(ggl_nn[[i]])
# }
# dev.off()

################################# Platelets ####################################
# plate_list <- list("Gp1ba", "Cd84", "Cd226", "Gp5", "F13a1", "F5", "Mfsd2b",
#                    "Ost4", "Ptgs1")
# for (i in 1:length(plate_list)){
# FeaturePlot(object = seu, features = plate_list[[i]], reduction = "umap",
#             cols = c("azure3", "deeppink3"), pt.size = 0.1) -> ggl_oo[[i]]
# VlnPlot(object = seu, features = plate_list[[i]], pt.size = 0) +
#   theme(axis.text.x = element_text(angle = 45, size = 10)) -> ggl_pp[[i]]
# }
# pdf(paste0("platelet_", proj_name, ".pdf"))
# for (i in 1:length(plate_list)){
# print(ggl_oo[[i]])
# print(ggl_pp[[i]])
# }
# dev.off()

####################### pancreatic neuroendocrine hormones #####################
# hormone_list <- list("Gast", "Gcg", "Ppy", "Sst", "Vip")
# for (i in 1:length(hormone_list)){
# FeaturePlot(object = seu, features = hormone_list[[i]], reduction = "umap",
#             cols = c("azure3", "deeppink3"), pt.size = 0.1) -> ggl_qq[[i]]
# VlnPlot(object = seu, features = hormone_list[[i]], pt.size = 0) +
#   theme(axis.text.x = element_text(angle = 45, size = 10)) -> ggl_rr[[i]]
# }
# pdf(paste0("hormone_", proj_name, ".pdf"))
# for (i in 1:length(hormone_list)){
# print(ggl_qq[[i]])
# print(ggl_rr[[i]])
# }
# dev.off()
# 
# ####################### Pancreatic stem cell markers ###########################
# stem_list <- list("Lgr5", "Sox9", "Foxj1", "Ehf")
# for (i in 1:length(stem_list)){
# FeaturePlot(object = seu, features = stem_list[[i]], reduction = "umap",
#             cols = c("azure3", "deeppink3"), pt.size = 0.1) -> ggl_ss[[i]]
# VlnPlot(object = seu, features = stem_list[[i]], pt.size = 0) +
#   theme(axis.text.x = element_text(angle = 45, size = 10)) -> ggl_tt[[i]]
# }
# pdf(paste0("stem_", proj_name, ".pdf"))
# for (i in 1:length(stem_list)){
# print(ggl_ss[[i]])
# print(ggl_tt[[i]])
# }
# dev.off()
# 
# ############################### MAST cells #####################################
# mast_list <- list("Kit", "Cpa3", "Cma1", "Il1rl1", "Osbpl8", "Adora3", "Csf2rb",
#                   "Cyp11a1", "Hdc", "Il4", "Tcf4", "Bcl11a", "Irf8", "Spib",
#                   "Runx2", "Gpnmb")
# for (i in 1:length(mast_list)){
# FeaturePlot(object = seu, features = mast_list[[i]], reduction = "umap",
#             cols = c("azure3", "deeppink3"), pt.size = 1) -> ggl_uu[[i]]
# VlnPlot(object = seu, features = mast_list[[i]], pt.size = 0) +
#   theme(axis.text.x = element_text(angle = 45, size = 10)) -> ggl_vv[[i]]
# }
# pdf(paste0("mast_", proj_name, ".pdf"))
# for (i in 1:length(mast_list)){
# print(ggl_uu[[i]])
# print(ggl_vv[[i]])
# }
# dev.off()
# 
# ######################### peri-islet schwann cells #############################
# peri_list <- list("Egfl8", "Gulp1", "Gfra3", "Fign", "Ngfr")
# for (i in 1:length(peri_list)){
# FeaturePlot(object = seu, features = peri_list[[i]], reduction = "umap",
#             cols = c("azure3", "deeppink3"), pt.size = 0.1) -> ggl_ww[[i]]
# VlnPlot(object = seu, features = peri_list[[i]], pt.size = 0) +
#   theme(axis.text.x = element_text(angle = 45, size = 10)) -> ggl_xx[[i]]
# }
# pdf(paste0("peri-islet_", proj_name, ".pdf"))
# for (i in 1:length(peri_list)){
# print(ggl_ww[[i]])
# print(ggl_xx[[i]])
# }
# dev.off()
# 
# ############################## NK-cells ########################################
# nk_list <- list("Gzma", "Ms4a4b", "Ctla2a", "Ctsw", "Ncr1", "Klrk1",
#                 "Il2rb", "Dusp2", "Klrb1c", "Klre1", "Ccl5")
# for (i in 1:length(nk_list)){
# FeaturePlot(object = seu, features = nk_list[[i]], reduction = "umap",
#             cols = c("azure3", "deeppink3"), pt.size = 0.1) -> ggl_yy[[i]]
# VlnPlot(object = seu, features = nk_list[[i]], pt.size = 0) +
#   theme(axis.text.x = element_text(angle = 45, size = 10)) -> ggl_zz[[i]]
# }
# pdf(paste0("NK_", proj_name, ".pdf"))
# for (i in 1:length(nk_list)){
# print(ggl_yy[[i]])
# print(ggl_zz[[i]])
# }
# dev.off()
# 
# ############################## pDCs ##########################################
# 
# pdc_list <- list("Pira2", "Cxcr3", "Irf7")
# for (i in 1:length(pdc_list)){
# FeaturePlot(object = seu, features = pdc_list[[i]], reduction = "umap",
#             cols = c("azure3", "deeppink3"), pt.size = 1) -> ggl_yy[[i]]
# VlnPlot(object = seu, features = pdc_list[[i]], pt.size = 0) +
#   theme(axis.text.x = element_text(angle = 45, size = 10)) -> ggl_zz[[i]]
# }
# pdf(paste0("pdcs_", proj_name, ".pdf"))
# for (i in 1:length(pdc_list)){
# print(ggl_yy[[i]])
# print(ggl_zz[[i]])
# }
# dev.off()
# 
# ######################### Endothelial cells subtypes  ##########################
# endo_list2 <- list("Aqp7", "Meox2", "Lpl", "Fabp5", "Cd36", "Tcf15", "Rgcc",
#                   "Kdr", "Prox1", "Pdpn",
#                   "Madcam1", "Lrg1", "Ackr1", # HEV markers
#                   "Cd34", "Pdgfb", "Kcne3", "Nid2", "Dll4", ## tip
#                   "Rgcc", "Kdr", "Cd300lg", "Ramp3", #capillaries
#                   "Mmrn1", "Flt4", #lymphatic
#                   "Emcn", "Vwf", "Bgn", # vein
#                   "Sox17", "Stmn2", #artery
#   "Gja4",	"Cxcl12",	"Efnb2", # Cap Arterial
#   "Col4a1",	"Col4a2",	"Sema6d", # Stalk cell
#   "Esm1",	"Nid1", "Nid2", "Trp53i11",	"Kcne3", "Apln", "Cd34", "Pdgfb", "Kcne3",
#   "Dll4", ## Tip cells
#   "Cd36",	"Plvap",	"Igfbp7",	"Gja5",	"Cxcl12", "Rgcc", "Kdr", "Cd300lg",
#   "Ramp3" #capillaries
# )
# 
# for (i in 1:length(endo_list2)){
# FeaturePlot(object = seu, features = endo_list2[[i]], reduction = "umap",
#             cols = c("azure3", "deeppink3"), pt.size = 0.1) -> ggl_aaa[[i]]
# VlnPlot(object = seu, features = endo_list2[[i]], pt.size = 0) +
#   theme(axis.text.x = element_text(angle = 45, size = 10)) -> ggl_bbb[[i]]
# }
# pdf(paste0("endothelial_cells_sub", proj_name, ".pdf"))
# for (i in 1:length(endo_list2)){
# print(ggl_aaa[[i]])
# print(ggl_bbb[[i]])
# }
# dev.off()
# 
######################### macro cells subtypes  ##########################
macro_list2 <- list(
  "Hgf", "Ccl3", "Hdc", "Ms4a2", "Il4", "Itgb7", "Cpa3",
   "Il3ra", "Ccl9", "Cd69", # basophils
   "Ear2", "Ear2", "Prg3", "Prg2", "Pglyrp1", # eosinophils
   "Lyz1", "Csf1r", "Cd300e", "Mafb", "Batf3", "Krt79", "Msr1", # DC
   "Siglech", "Ccr9", "Bst2", "Pacsin1", "Tcf4", "Xcr1", "Clec9a", "Cd83",
   "S100a9", "S100a8", "Mmp9", "Mmp8", "Csf3r", "Il1rn", "Cxcr2", "Ly6g", #Neutro
   "Il6", "Gata2", "Cpa3", "Ms4a2", "Fcer1a", # basophils (Cd117 neg-)
   "Adgre1", "Kit", "Ccr3", # eosinophils Kit = Cd117, Ccr3 = Cd193
   "Cd86", "Cd80", # M1-like macrophages
   "Cd163", "Mrc1", # M2-like macrophages
   "Arg1", # suppressor cells
   "Cd14", # monocytes "Adgre1" neg-
   "Il1b", "Tnf", ## "Adgre1" is F4/80 marker # MAST (Fcer1a+ Cd117+ )
   "Klf2", "Gata6", "Cd93", "Irf4", "H2-Eb1", "Cd74", "Tgfbi",
   # Lilra is pDC marker gene
   "Xcr1", "Flt3", "Itgae", "Zbtb46", "Itgax", # for cDC
   "Siglech", "Bst2", # for the pDCs
   "Cx3cr1", "Ccr2", "Sell", "Ly6c1"  # "Cd62l", "Cd43", "H2ab" , not present!
   )
for (i in 1:length(macro_list2)){
FeaturePlot(object = seu, features = macro_list2[[i]], reduction = "umap",
            cols = c("azure3", "deeppink3"), pt.size = 0.1) -> ggl_ccc[[i]]
VlnPlot(object = seu, features = macro_list2[[i]], pt.size = 0) +
  theme(axis.text.x = element_text(angle = 45, size = 10)) -> ggl_ddd[[i]]
}
pdf(paste0("macro_marker_", proj_name, ".pdf"))
for (i in 1:length(macro_list2)){
print(ggl_ccc[[i]])
print(ggl_ddd[[i]])
}
dev.off()

# ############################ tcell cells subtypes  #############################
# tcell_list2 <- list("Ikzf2", "Maf", "Ctla4", "Foxp3", "Gzmb", # TREGS
#                     "Ptprc", "Lat", "Arhgap45", "Ms4a6b", "Ltb", "Cd3g", "Cd3d", 
#                     "Arhgef1", "H2-Q7", # T-memory
#                     "Tox", "Slamf6", #Tfh
#                     "Cx3cr1", "Pdcd1", ## CD8 Exhausted
#                     "Ifng", "Klra1", "S1pr1", "Il17a", "Cxcr5",
#                     "Ifngr1" # not present "Ccr7, Tcf7" #resting T-cells/naive
# )
# 
# for (i in 1:length(tcell_list2)){
# FeaturePlot(object = seu, features = tcell_list2[[i]], reduction = "umap",
#             cols = c("azure3", "deeppink3"), pt.size = 0.1) -> ggl_eee[[i]]
# VlnPlot(object = seu, features = tcell_list2[[i]], pt.size = 0) +
#   theme(axis.text.x = element_text(angle = 45, size = 10)) -> ggl_fff[[i]]
# }
# pdf(paste0("tcells_sub", proj_name, ".pdf"))
# for (i in 1:length(tcell_list2)){
# print(ggl_eee[[i]])
# print(ggl_fff[[i]])
# }
# dev.off()
# 
# ############################ proliferating marker  #############################
# prolif <- list("Mki67")
# 
# for (i in 1:length(prolif)){
# FeaturePlot(object = seu, features = prolif[[i]], reduction = "umap",
#             cols = c("azure3", "deeppink3"), pt.size = 0.1) -> ggl_iii[[i]]
# VlnPlot(object = seu, features = prolif[[i]], pt.size = 0) +
#   theme(axis.text.x = element_text(angle = 45, size = 10)) -> ggl_jjj[[i]]
# }
# pdf(paste0("proliferating_marker_", proj_name, ".pdf"))
# for (i in 1:length(prolif)){
# print(ggl_iii[[i]])
# print(ggl_jjj[[i]])
# }
# dev.off()
# 
# ############################ extra T-cell markers  #############################
# pdf(paste0("Foxp3_", proj_name, ".pdf"))
# FeaturePlot(object = seu, features = "Foxp3", reduction = "umap",
#             cols = c("azure3", "deeppink3"), pt.size = 0.2) -> p1
# print(p1)
# dev.off()
# 
# pdf(paste0("Itga2_", proj_name, ".pdf"))
# FeaturePlot(object = seu, features = "Itga2", reduction = "umap",
#             cols = c("azure3", "deeppink3"), pt.size = 0.2) -> p2
# print(p2)
# dev.off()
# 
# pdf(paste0("Lag3_", proj_name, ".pdf"))
# FeaturePlot(object = seu, features = "Lag3", reduction = "umap",
#             cols = c("azure3", "deeppink3"), pt.size = 0.2) -> p3
# print(p3)
# dev.off()
# 
# pdf(paste0("Cd4_", proj_name, ".pdf"))
# FeaturePlot(object = seu, features = "Cd4", reduction = "umap",
#             cols = c("azure3", "deeppink3"), pt.size = 0.2) -> p4
# print(p4)
# dev.off()
# 
# pdf(paste0("Cd226_", proj_name, ".pdf"))
# FeaturePlot(object = seu, features = "Cd226", reduction = "umap",
#             cols = c("azure3", "deeppink3"), pt.size = 0.2) -> p5
# print(p5)
# dev.off()
# 
# ###
# pdf(paste0("Chst4_", proj_name, ".pdf"))
# FeaturePlot(object = seu, features = "Chst4", reduction = "umap",
#             cols = c("azure3", "deeppink3"), pt.size = 0.2) -> p5
# print(p5)
# dev.off()
# ###


############################# custom markers ###################################

# pdc_list <- list("Atg5", "Pdgfb")
# for (i in 1:length(pdc_list)){
# FeaturePlot(object = seu, features = pdc_list[[i]], reduction = "umap",
#             cols = c("azure3", "deeppink3"), pt.size = 1) -> ggl_yy[[i]]
# VlnPlot(object = seu, features = pdc_list[[i]], pt.size = 1) +
#   theme(axis.text.x = element_text(angle = 45, size = 10)) -> ggl_zz[[i]]
# DotPlot(object = seu, features = pdc_list[[i]], group.by = "primary_annotations_v1"
#         ) -> ggl_aa[[i]]
# }
# pdf(paste0("Atg5_", proj_name, ".pdf"))
# for (i in 1:length(pdc_list)){
# print(ggl_yy[[i]])
# print(ggl_zz[[i]])
# print(ggl_aa[[i]])
# }
# dev.off()