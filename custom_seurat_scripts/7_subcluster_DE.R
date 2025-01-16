############################### sc-pipeline ####################################
## this script performs downstream analysis on Seurat objects taken from the
## previous script `2_seurat.R` and filters the objects further, normalizes the
## data and finds variable features. running PCA, UMAP and generating various
## diagnostic plots. Lastly, the clusters resulting from the UMAP are annotated 
## sample and response type, filters the cells, produces QC plots,
## NOTE: if you need to run harmony, you need to use R 3.6.3 libraries
## /media/seq-srv-05/vrc/Project/Project_Gino/software/R/R-3.6.3
## /data/projects/Project_Theo/R/R-4.0.3

################################################################################
# For EC or immune cells will be quite easy, as this is already studied. 
# 
# For exemple, in the T/NK cluster you should be able to identify T CD4+ cells,
# T CD8+ cells , regulatory T cells, NK cells and others sub clusters of
# activated cells, memory cells or exhausted cells etc... For myeloid cells we
# can probably find resident macrophages vs monocytes derived macrophages but
# also M1 vs M2 for exemple. The same for the neutrophils N1 vs N2. M1 or N1 are
# some anti-tumoral cells and N2 or M2 are more pro-angiogenic and pro-tumoral
# cells.... I am sure that we can probably identify some dendritic cells
# somewhere too. 
# 
# For the EC this is well know too. Tip cells, capillary or venous cells, HEVs,
# Lymphatic cells, arterial cells etc...
# 
# For the tumor cells it is a bit less defined and we will have to see.... But
# I can guess that we can find some invasive cell cluster, stem cell like or
# proliferative cells.
# 
# We can always subclustered and see the top express genes to identify the
# subclusters or at least have some hypothesis.
# 
# Other strategy we can also try to identify clusters that are affected by
# treatment or between responding and relapsing samples and try to identify
# them. 
# 
# But indeed, this we can discuss with Gabriele and Yichao as soon as you have
# the definitive clustering. 
################################################################################

## load libraries
library("R.utils")
library("dplyr")
library("Seurat")
library("ggplot2")
library("future")
library("cowplot")

## parallelize workflow
## only run this code chunk if you are on the VSC
plan("multiprocess", workers = 30) # uses 50 CPU
options(future.globals.maxSize = 14000 * 1024^2) ## 5GB per worker

## evergreen variables
proj_name <- "pancreatic_mice" ## name of project
proj_version <- "filter"       ## name of version of project
# server_path <- "/media/seq-srv-05/vrc/Project/Project_Theo/pancreatic_mice/unfiltered/macro_sub2/"
server_path <- "/media/seq-srv-05/vrc/Project/Project_Theo/pancreatic_mice/raw_1000dr/"
mito <- "^mt-" ## "^mt-" if mouse, "^MT-" if human
up_feat <- 6000 ## upper nFeature threshold
low_feat <- 1000 ## lower nFeature threshold
var_feat <- 2000 ## variableFeature threshold
low_count <- 1000 ## lower nCount threshold
p_mt <- 20 ## percent mito threshold
NPC <- 35 ## number of PCs of first PCA
sig_dims <- 10 ## number of PCs of final PCA
resolution <- 0.4 ## resolution
RESOLUTION <- 'RNA_snn_res.0.4' ## resolution string

########################### subclustering CANCER ##11111111111111111111111111111

### original cancer subclustering
# seu <- readRDS(file = paste0("FINAL_pancreatic_mice.rds"))
# proj_name <- "cancer" ## name of project
# c1 <- subset(seu, subset = Cluster == c("Cancer"))
# c1 <- SetIdent(c1, value = "orig.ident")

### second round of subclustering looking for IT and MLP subtypes
# seu <- readRDS("FINAL_cancer_ann_pancreatic_mice.rds")
seu <- readRDS("../FINAL_latest_ann_actual.rds")
c1 <- subset(seu, subset = Cluster == c("Cancer"))
c2 <- subset(seu, subset = Cluster == c("Cancer_B20_PDL1_CTLA4_Resp"))
c3 <- subset(seu, subset = Cluster == c("Cancer_B20_PDL1_Resp1"))
c4 <- subset(seu, subset = Cluster == c("Cancer_B20_PDL1_Resp2"))
c5 <- subset(seu, subset = Cluster == c("Cancer_B20_Pi3Ki_Resp"))
c6 <- subset(seu, subset = Cluster == c("Ductal/Stem"))
c7 <- subset(seu, subset = Cluster == c("Acinar"))

#  "Cancer"                        "Macrophages"                  
#  [3] "Acinar"                        "Cancer_B20_PDL1_CTLA4_Resp"   
#  [5] "Endothelial"                   "T-cells/NK_B20_PDL1_LTbR_Resp"
#  [7] "B-Cell_B20_PDL1_LTbR_Resp"     "Cancer_B20_PDL1_Resp1"        
#  [9] "T-cells/NK"                    "Cancer_B20_PDL1_Resp2"        
# [11] "Monocytes"                     "Cancer_B20_Pi3Ki_Resp"        
# [13] "Fibroblasts"                   "Ductal/Stem"

# subset(seu, subset = Cluster == c("Cancer", "Cancer_B20_PDL1_CTLA4_Resp",
#                                   "Cancer_B20_PDL1_Resp1", "Cancer_B20_PDL1_Resp2",
#                                   "Cancer_B20_Pi3Ki_Resp", "Ductal/Stem")) -> c_main
## merge the cancer clusters and recluster
# c1 %>% merge(y = c(c2, c3, c4, c5, c6), project = proj_name) -> c_g
c1 %>% merge(y = c(c2, c3, c4, c5, c6, c7), project = proj_name) -> c_0
# saveRDS(c_main, file = "cancer_main.rds")
# saveRDS(c_g, file = "cancer_grouped.rds")

### add contamination data to see where it clusters
## B-cells
# seu1 <- readRDS("sub_b_cells.rds")
# seu_ins <- subset(seu1, subset = Ins1 > 1 & Ins2 > 1)
# saveRDS(seu_ins, file = "cancer_bcells.rds")
# ## endo
# seu2 <- readRDS("endo_sub2/FINAL_latest_ann_endo.rds") ## 3590 samples
# seu3 <- readRDS("endo_sub/final_ann_endo.rds") ## 2871 samples
# seu_ins <- subset(seu2, cells = colnames(seu2)[!(Cells(seu2) %in% Cells(seu3))])
# saveRDS(seu_ins, file = "cancer_endo.rds")
# ## macro
# seu1 <- readRDS("sub_macro.rds") ## 3529 samples
# seu2 <- readRDS("macro_recluster_4.rds") ## 2832 samples
# seu_ins <- subset(seu1, cells = colnames(seu1)[!(Cells(seu1) %in% Cells(seu2))])
# saveRDS(seu_ins, file = "cancer_macro.rds") # 697 samples
# ## t-cells
# seu1 <- readRDS("sub_tcell.rds") ## 5531 samples
# seu2 <- readRDS("tcells_actual_FINAL.rds") ## 4679 samples
# seu_ins <- subset(seu1, cells = colnames(seu1)[!(Cells(seu1) %in% Cells(seu2))])
# saveRDS(seu_ins, file = "cancer_tcell.rds") # 852 samples

## read objects and change IDs
c_b <- readRDS("cancer_bcells.rds")
c_e <- readRDS("cancer_endo.rds")
c_m <- readRDS("cancer_macro.rds")
c_t <- readRDS("cancer_tcell.rds")
c_b@meta.data$Cluster <- "Contamination_B-Cells"
c_e@meta.data$Cluster <- "Contamination_EC"
c_m@meta.data$Cluster <- "Contamination_Macro"
c_t@meta.data$Cluster <- "Contamination_T-Cells"
#c_0 <- readRDS("cancer_grouped.rds")

seu@meta.data$new_labels <- as.character(seu@meta.data$Cluster)
seu@meta.data$new_labels <- gsub("Cancer_B20_PDL1_CTLA4_Resp", "Cancer", seu@meta.data$new_labels)
seu@meta.data$new_labels <- gsub("Cancer_B20_PDL1_Resp1", "Cancer", seu@meta.data$new_labels)
seu@meta.data$new_labels <- gsub("Cancer_B20_PDL1_Resp2", "Cancer", seu@meta.data$new_labels)
seu@meta.data$new_labels <- gsub("Cancer_B20_Pi3Ki_Resp", "Cancer", seu@meta.data$new_labels)
seu@meta.data$new_labels <- gsub("_", " ", seu@meta.data$new_labels)

## merge the cancer clusters and recluster
c_0 %>% merge(y = c(c_b, c_e, c_m, c_t), project = proj_name) -> c_z
# c_z <- SetIdent(object = c_z, value = "Cluster")
# c_z <- SetIdent(object = c_0, value = "Cluster")

c_z %>% ## trying with only the cancer therapies and ductal
# c1 %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = var_feat) %>%
  # ScaleData(vars.to.regress = c("nCount_RNA", "percent.mt", "orig.ident")) %>% 
  ScaleData() %>% 
  RunPCA(npcs = sig_dims, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:sig_dims, k.param = 30,
                compute.SNN = TRUE, nn.method = "rann", nn.eps = 0) %>% 
  FindClusters(resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                              1.1, 1.2, 1.3, 1.4)) %>% 
  RunUMAP(reduction = "pca", dims = 1:sig_dims) -> FC

## multiple resolution plots
pdf(file = paste0("clusters_diff_res_", proj_name, ".pdf"), width = 15, height = 8)
Seurat_02 <- SetIdent(object = FC, value = 'RNA_snn_res.0.2')
Seurat02 <- RunUMAP(Seurat_02, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat02, reduction = "umap", pt.size = 0.25, label = TRUE) +
  labs(title = "UMAP with res 0.2") + theme(plot.title = element_text(hjust = 0.5))
Seurat_04 <- SetIdent(object = FC, value = 'RNA_snn_res.0.4')
Seurat04 <- RunUMAP(Seurat_04, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat04, reduction = "umap", pt.size = 0.25, label = TRUE) +
  labs(title = "UMAP with res 0.4") + theme(plot.title = element_text(hjust = 0.5))
Seurat_06 <- SetIdent(object = FC, value = 'RNA_snn_res.0.6')
Seurat06 <- RunUMAP(Seurat_06, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat06, reduction = "umap", pt.size = 0.25, label = TRUE) +
  labs(title = "UMAP with res 0.6") + theme(plot.title = element_text(hjust = 0.5))
Seurat_08 <- SetIdent(object = FC, value = 'RNA_snn_res.0.8')
Seurat08 <- RunUMAP(Seurat_08, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat08, reduction = "umap", pt.size = 0.25, label = TRUE) +
  labs(title = "UMAP with res 0.8") + theme(plot.title = element_text(hjust = 0.5))
Seurat_1 <- SetIdent(object = FC, value = 'RNA_snn_res.1.0')
Seurat1 <- RunUMAP(Seurat_1, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat1, reduction = "umap", pt.size = 0.25, label = TRUE) + 
  labs(title = "UMAP with res 1.0") + theme(plot.title = element_text(hjust = 0.5))
Seurat_12 <- SetIdent(object = FC, value = 'RNA_snn_res.1.2')
Seurat12 <- RunUMAP(Seurat_12, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat12, reduction = "umap", pt.size = 0.25, label = TRUE) + 
  labs(title = "UMAP with res 1.2") + theme(plot.title = element_text(hjust = 0.5))
Seurat_14 <- SetIdent(object = FC, value = 'RNA_snn_res.1.4')
Seurat14 <- RunUMAP(Seurat_14, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat14, reduction = "umap", pt.size = 0.25, label = TRUE) +
  labs(title = "UMAP with res 1.4") + theme(plot.title = element_text(hjust = 0.5))
dev.off() 

## write the column name of the chosen resolution with the SetIdent function
## Seurat will set the object with this column
seu <- SetIdent(object = FC, value = RESOLUTION)
# seu <- SetIdent(object = FC, value = "Cluster")
seu <- SetIdent(object = seu, value = "new_labels")
seu <- RunUMAP(seu, reduction = "pca", dims = 1:sig_dims)
pdf(file = paste0("res_", RESOLUTION, "_clusters_", proj_name, ".pdf"),
    width = 15, height = 8)
UMAPPlot(seu, reduction = "umap", pt.size = 0.1, label = TRUE,
         label.size = 5) + NoLegend() +
  labs(title = paste0("UMAP with ", RESOLUTION)) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

## produce dimplots
pdf(paste0("umap_", proj_name, ".pdf"), width = 15, height = 8)
p1 <- DimPlot(seu, reduction = "umap", group.by = "orig.ident",
              pt.size = 0.1, label = FALSE, label.size = 2.0)
p2 <- DimPlot(seu, reduction = "umap", group.by = "Response",
              pt.size = 0.1, label = FALSE, label.size = 3.5,  repel = TRUE)
p3 <- DimPlot(seu, reduction = "umap", group.by = "Treatment",
              pt.size = 0.1, label = FALSE, label.size = 3.5,  repel = TRUE)
p4 <- DimPlot(seu, reduction = "umap",
              pt.size = 0.1, label = TRUE, label.size = 3.5,  repel = FALSE)
plot_grid(p1, p2, p3, p4, ncol = 2)
dev.off()

## umap split on response, sample and treatment 
# pdf(paste0("Split_Response_", proj_name, ".pdf"), width = 15, height = 8)
# DimPlot(seu, reduction = "umap", split.by = "Response")
# dev.off()
# pdf(paste0("Split_Sample_", proj_name, ".pdf"), width = 15, height = 8)
# DimPlot(seu, reduction = "umap", split.by = "orig.ident")
# dev.off()
# pdf(paste0("Split_Treatment_", proj_name, ".pdf"), width = 15, height = 8)
# DimPlot(seu, reduction = "umap", split.by = "Treatment")
# dev.off()

########################## tables for QC plots ################################
## make sure you see the correct number of clusters. if not, SetIdent to the
## correct column in  metadata
table(seu@active.ident) 

## save select metadata
u <- seu@meta.data

## table of number of cells 
write.table(u, paste0("cells_metadata_", proj_name, ".txt"), col.names = NA,
            sep = "\t")

## how many cells in each cluster per sample
cell_perC_perO <- table(seu@active.ident, seu@meta.data$orig.ident)
write.table(cell_perC_perO, paste0("cell_per_c_per_O_", proj_name, ".txt"),
            col.names = NA, sep = "\t")

## how many cells in each cluster per treatment
cell_perC_perT <- table(seu@active.ident, seu@meta.data$Treatment)
write.table(cell_perC_perT, paste0("cell_per_c_per_T_", proj_name, ".txt"),
            col.names = NA, sep = "\t")

## how many cells in each cluster per response
cell_perC_perR <- table(seu@active.ident, seu@meta.data$Response)
write.table(cell_perC_perR, paste0("cell_per_c_per_R_", proj_name, ".txt"),
            col.names = NA, sep = "\t")

c("IT-like", "IT-like", "IT-like", "IT-like", "IT-like B20 PiK3i",
  "Unknown B20 PiK3i Rel", "Ins-hi MLP-like", "IT-like",
  "Unknown B20 Resp", "IT-like", "IT-like B20 Resp", "Ins-lo MLP-like B20",
  "IT-like B20 aPDL1 CTLA4 Resp",  "IT-like B20 aPDL1 CTLA4 Resp",
  "Ins-lo MLP-like B20 Resp", "IT-like", "Unknown B20 PiK3i Rel",
  "Ins-lo MLP-like B20", "Unknown B20 PiK3i Rel", "IT-like", "Ins-lo MLP-like",
  "Acinar", "Ductal/Stem", "IT-like B20 aPDL1 CTLA4 Resp",
  "Ins-hi MLP-like B20 PiK3i Rel", "MLP-like CCB104"
  ) -> new.cluster.ids

# c("Invasive front A CCB081",
#   "Invasive front B CCB081",
#   "IT-like CCB082",
#   "Ins-hi MLP", "IT-like CCB082", "Acinar cells CCB082", "Ins-lo MLP",
#   "IT-like") -> new.cluster.ids

# ## add column in metadata with the "Cluster" title and the cluster name
names(new.cluster.ids) <- levels(seu)
seu <- RenameIdents(seu, new.cluster.ids)
seu <- StashIdent(seu, save.name = "new_labels")
head(seu@meta.data) ## verify new column
#
# pdf(paste0("Annotated_", proj_name, ".pdf"))
# DimPlot(seu, reduction = "umap", label = FALSE, pt.size = 0.1) +
#   labs(title = "Annotated UMAP") + theme(plot.title = element_text(hjust = 0.5))
# dev.off()

pdf(paste0("Annotated_labels_", proj_name, ".pdf"), width = 11, height = 8.5) 
DimPlot(seu, reduction = "umap", label = FALSE, pt.size = 0.1, repel = TRUE) +
  labs(title = "Annotated UMAP") + theme(plot.title = element_text(hjust = 0.5))
dev.off()
# 
# ## save the results


########################

## Plot number of genes, nUMIs and %mito on a tSNE - 3 plots per sample
## These allows you to see if part of your umap that has a high myeloid content
## (dying cells) or high number of genes (normally, no 2 cells in a droplets
## because you still filtered these bad cells!)
pdf(paste0("graph_nfeature_", proj_name, ".pdf"), width = 15, height = 8)
FeaturePlot(seu, features = "nFeature_RNA", cols = c("cadetblue2", "darkred"),
            pt.size = 0.1, reduction = "umap", label = FALSE)
dev.off()
pdf(paste0("graph_nCount_", proj_name, ".pdf"), width = 15, height = 8)
FeaturePlot(seu, features = "nCount_RNA", cols = c("cadetblue2", "darkred"),
            pt.size = 0.1, reduction = "umap", label = FALSE)
dev.off()
pdf(paste0("graph_percentmt_", proj_name, ".pdf"), width = 15, height = 8)
FeaturePlot(seu, features = "percent.mt", cols = c("cadetblue2", "darkred"),
            pt.size = 0.1, reduction = "umap", label = FALSE)
dev.off()
# 
# ## export nFeature, nCOunt and %mito per cluster 
# features <- FetchData(seu, vars = c("nFeature_RNA", "nCount_RNA",
#                                     "percent.mt", RESOLUTION), cells = NULL)
# write.table(features, paste0("table_features_", proj_name, ".txt"),
#             col.names = NA, sep = "\t")

## enter as many number as you have clusters and picked the resolution you chose
# f <- ordered(features$RNA_snn_res.0.4, levels = c(0:clust_num))
# f <- ordered(features$RNA_snn_res.1.2,
#              levels = c(0:(max(as.integer(seu@meta.data$RNA_snn_res.1.2)) - 1)))
# pdf(paste0("boxplot_nFeature_v2_", proj_name, ".pdf"), width = 15, height = 8)
# box = boxplot(nFeature_RNA ~ f , data = features, xlab = "# genes",
#               ylab = "Cluster", main = "Number of Genes per cluster",
#               col = "dodgerblue3",  horizontal = TRUE)
# dev.off()
# 
# pdf(paste0("boxplot_nCount_v2_", proj_name, ".pdf"), width = 15, height = 8)
# box = boxplot(nCount_RNA ~ f , data = features, xlab = "# Barcodes", ylab = "Cluster",
#               main = "Number of Barcodes per cluster", col = "dodgerblue3",
#               horizontal = TRUE, ylim = c(0, 60000)) # , ylim = c(0, 35)
# dev.off()
# 
# pdf(paste0("boxplot_percentmito_v2_", proj_name, ".pdf"), width = 15, height = 8)
# box = boxplot(percent.mt ~ f , data = features, xlab = "% mito", ylab = "Cluster",
#               main = "Percentage of mito genes per cluster", col = "dodgerblue3",
#               horizontal = TRUE)
# dev.off()
# 
# ## save the results
# saveRDS(seu, file = paste0("sub_", proj_name, ".rds"))
#seu <- readRDS(file = paste0("sub_cancer", ".rds"))

# sample.cluster <- subset(seu, idents = sample)
 # print(paste0("sample:", sample))
 length(rownames(seu@meta.data))
 ## ceiling value is a default - can be adapted for specific projects
 expected.doublets <- ceiling(0.039 * length(rownames(seu@meta.data)))
 ## apply doubletFinder function
 doubletFinder_v3(seu,
                  PCs = 1:20, # number of significant PCs
                  pN = 0.25, # defines number of artificial doublets
                  pK = 0.01, # needs to be estimated
                  nExp = expected.doublets,
                  reuse.pANN = FALSE, sct = TRUE) -> seu
 seu@meta.data[colnames(seu),
     paste("pANN_0.25_0.01", expected.doublets, sep = "_")] -> seu$pANN
 seu@meta.data[colnames(seu),
     paste("DF.classifications_0.25_0.01",
           expected.doublets, sep = "_")] -> seu$pANNPredictions
 seu$pANN[colnames(seu)] <- seu$pANN[colnames(seu)]
 seu$pANNPredictions[colnames(seu)] <- seu$pANNPredictions[colnames(seu)]

## new doublets can be found here:
# unique(seu@meta.data$DF.classifications_0.25_0.01_2916) ## doublets YES!
seu <- SetIdent(seu, value = "DF.classifications_0.25_0.01_2916")
seu <- SetIdent(seu, value = "Cluster")

## produce dimplots
pdf(paste0("feature_umap_", proj_name, ".pdf"), width = 15, height = 15)
p1 <- DimPlot(seu, reduction = "umap", group.by = "orig.ident",
              pt.size = 0.1, label = FALSE, label.size = 2.0)
p2 <- DimPlot(seu, reduction = "umap", group.by = "Response",
              pt.size = 0.1, label = FALSE, label.size = 3.5,  repel = TRUE)
p3 <- DimPlot(seu, reduction = "umap", group.by = "Treatment",
              pt.size = 0.1, label = FALSE, label.size = 3.5,  repel = TRUE)
p4 <- DimPlot(seu, reduction = "umap", group.by = "Cluster",
              pt.size = 0.1, label = FALSE, label.size = 3.5,  repel = FALSE)
p5 <- DimPlot(seu, reduction = "umap", group.by = "DF.classifications_0.25_0.01_2916",
              pt.size = 0.1, label = FALSE, label.size = 3.5,  repel = FALSE)
p6 <- FeaturePlot(seu, features = "nFeature_RNA", cols = c("cadetblue2", "darkred"),
            pt.size = 0.1, reduction = "umap", label = FALSE)
p7 <- FeaturePlot(seu, features = "nCount_RNA", cols = c("cadetblue2", "darkred"),
            pt.size = 0.1, reduction = "umap", label = FALSE)
p8 <- FeaturePlot(seu, features = "percent.mt", cols = c("cadetblue2", "darkred"),
            pt.size = 0.1, reduction = "umap", label = FALSE)
plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, ncol = 2)
dev.off()

########################### subclustering MACRO/MONO ##2222222222222222222222222
## cluster with monocytes
# m1 <- readRDS("macro_recluster_4.rds")
# m2 <- readRDS("monocytes_2.rds")
# m1 %>% merge(y = m2, project = proj_name) -> merged_mm
# merged_mm <- SetIdent(merged_mm, value = "orig.ident")
# seu_no_ins <- subset(seu, subset = Cluster != 14)
# seu_ins <- subset(seu, subset = Cluster == 14)
# saveRDS(seu_ins, "monocyte_contamination1.rds")
# saveRDS(seu_no_ins, file = "macro_mono2.rds")

seu_no_ins <- readRDS("macro_mono2.rds")

# seu_ins <- subset(m2, subset = Ins1 > 5 & Ins2 > 5)
# seu_no_ins <- subset(seu, subset = Ins1 < 1 | Ins2 < 1 | Ppy < 1)

## original round
# ## recluster
# seu_ins <- subset(seu, subset = Ins1 > 10 & Ins2 > 10)
# seu_no_ins <- subset(seu, subset = Ins1 < 1 | Ins2 < 1 | Ppy < 1)
# 
# 
# seu_no_ins <- subset(seu, subset = Cluster != 4)
# 
# seu_no_ins <- subset(seu, subset = Cluster != 7)
# 
# # seu_7 <- subset(seu, subset = Cluster == 7)
# 
# ########## macro subset of unfiltered 2
# # seu <- readRDS(file = paste0("FINAL_pancreatic_mice.rds"))
# # levels(seu)
# #  [1] "Acinar"                    "Neurons"                  
# #  [3] "Macrophages/Myeloids"      "Endothelial cells"        
# #  [5] "T cells"                   "B cells"                  
# #  [7] "Acinar/Cancer?"            "Fibroblasts/Smooth Muscle"
# #  [9] "Fibroblasts/Cancer?"       "Monocytes/Myeloids"       
# # [11] "Granulocytes/Neutrophils"
# 
# #    0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
# # 5034 2469 1457 1325 1307  900  898  853  849  352  339  271  243  209  132  108 
# #   16 
# #   61 
# 
# proj_name <- "macro4" ## name of project
# # m1 <- subset(seu, subset = Cluster == c("Macrophages", "Monocytes"))
# m1 <- subset(seu, subset = Cluster == c("Macrophages"))
# m2 <- subset(seu, subset = Cluster == c("Monocytes"))
# # m3 <- subset(seu, subset = Cluster == c("Granulocytes/Neutrophils"))
# m1 %>% merge(y = c(m2), project = proj_name) -> merged_mm
# merged_mm <- SetIdent(merged_mm, value = "orig.ident")
#saveRDS(merged_mm, file = paste0("merged_mm.rds"))
## merged_mm <- readRDS(file = paste0("merged_mm.rds")) ### this is low filter
# m1
# An object of class Seurat 
# 29043 features across 6842 samples within 1 assay 
# Active assay: RNA (29043 features, 2000 variable features)
#  2 dimensional reductions calculated: pca, umap
# > m2
# An object of class Seurat 
# 29043 features across 11117 samples within 1 assay 
# Active assay: RNA (29043 features, 2000 variable features)
#  2 dimensional reductions calculated: pca, umap
# > m3
# An object of class Seurat 
# 29043 features across 4998 samples within 1 assay 
# Active assay: RNA (29043 features, 2000 variable features)
#  2 dimensional reductions calculated: pca, umap
# > merged_mm
# An object of class Seurat 
# 29043 features across 22957 samples within 1 assay 
# Active assay: RNA (29043 features, 0 variable features)

######################## old macro subclustering ##############################
# seu <- readRDS(file = paste0("FINAL_pancreatic_mice.rds"))
# seu <- readRDS(file = paste0("FINAL_cancer_ann_pancreatic_mice.rds")) ### cstc
# proj_name <- "macro" ## name of project
# # m1 <- subset(seu, subset = Cluster == c("Macrophages", "Monocytes"))
# m1 <- subset(seu, subset = Cluster == c("Macrophages"))
# m2 <- subset(seu, subset = Cluster == c("Monocytes"))
# m1 %>% merge(y = m2) -> merged_mm
# merged_mm <- SetIdent(merged_mm, value = "orig.ident")

#merged_mm %>% 
seu_no_ins %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = var_feat) %>%
  ScaleData() %>% 
  RunPCA(npcs = NPC, verbose = FALSE) -> SC

## whatever cutoff for significant PCs, you will use to set sig_dims
## rerun PCA with significant PCs
SC %>% RunPCA(npcs = sig_dims, verbose = FALSE) -> SC

## save scaled
# saveRDS(SC, file = paste0("scaled_", proj_name, ".rds"))
# SC <- readRDS(file = paste0("scaled_", proj_name, ".rds"))

## Clustering and UMAP
## usually NPC = between 15-30 (around 20 is plenty)
## we calculate the clustering using multiple resolutions, the higher, the more
## number of clusters you will have. The point is that it makes biological sense
SC %>%
  FindNeighbors(reduction = "pca", dims = 1:sig_dims, k.param = 30,
                compute.SNN = TRUE, nn.method = "rann", nn.eps = 0) %>% 
  FindClusters(resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                              1.1, 1.2, 1.3, 1.4)) %>% 
  RunUMAP(reduction = "pca", dims = 1:sig_dims) -> FC

## multiple resolution plots
pdf(file = paste0("clusters_diff_res_", proj_name, ".pdf"), width = 15, height = 8)
Seurat_02 <- SetIdent(object = FC, value = 'RNA_snn_res.0.2')
Seurat02 <- RunUMAP(Seurat_02, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat02, reduction = "umap", pt.size = 0.25, label = TRUE) +
  labs(title = "UMAP with res 0.2") + theme(plot.title = element_text(hjust = 0.5))
Seurat_04 <- SetIdent(object = FC, value = 'RNA_snn_res.0.4')
Seurat04 <- RunUMAP(Seurat_04, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat04, reduction = "umap", pt.size = 0.25, label = TRUE) +
  labs(title = "UMAP with res 0.4") + theme(plot.title = element_text(hjust = 0.5))
Seurat_06 <- SetIdent(object = FC, value = 'RNA_snn_res.0.6')
Seurat06 <- RunUMAP(Seurat_06, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat06, reduction = "umap", pt.size = 0.25, label = TRUE) +
  labs(title = "UMAP with res 0.6") + theme(plot.title = element_text(hjust = 0.5))
Seurat_08 <- SetIdent(object = FC, value = 'RNA_snn_res.0.8')
Seurat08 <- RunUMAP(Seurat_08, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat08, reduction = "umap", pt.size = 0.25, label = TRUE) +
  labs(title = "UMAP with res 0.8") + theme(plot.title = element_text(hjust = 0.5))
Seurat_1 <- SetIdent(object = FC, value = 'RNA_snn_res.1.0')
Seurat1 <- RunUMAP(Seurat_1, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat1, reduction = "umap", pt.size = 0.25, label = TRUE) +
  labs(title = "UMAP with res 1.0") + theme(plot.title = element_text(hjust = 0.5))
Seurat_12 <- SetIdent(object = FC, value = 'RNA_snn_res.1.2')
Seurat12 <- RunUMAP(Seurat_12, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat12, reduction = "umap", pt.size = 0.25, label = TRUE) +
  labs(title = "UMAP with res 1.2") + theme(plot.title = element_text(hjust = 0.5))
Seurat_14 <- SetIdent(object = FC, value = 'RNA_snn_res.1.4')
Seurat14 <- RunUMAP(Seurat_14, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat14, reduction = "umap", pt.size = 0.25, label = TRUE) +
  labs(title = "UMAP with res 1.4") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

# seu <- SetIdent(object = FC, value = 'RNA_snn_res.0.6')

## write the column name of the chosen resolution with the SetIdent function
## Seurat will set the object with this column
seu <- SetIdent(object = FC, value = RESOLUTION)
seu <- RunUMAP(seu, reduction = "pca", dims = 1:sig_dims)
pdf(file = paste0("res_", RESOLUTION, "_clusters_", proj_name, ".pdf"),
    width = 15, height = 8)
UMAPPlot(seu, reduction = "umap", pt.size = 0.1, label = TRUE,
         label.size = 5) + NoLegend() +
  labs(title = paste0("UMAP with ", RESOLUTION)) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

## produce dimplots
pdf(paste0("umap_", proj_name, ".pdf"), width = 15, height = 8)
p1 <- DimPlot(seu, reduction = "umap", group.by = "orig.ident",
              pt.size = 0.1, label = FALSE, label.size = 2.0)
p2 <- DimPlot(seu, reduction = "umap", group.by = "Response",
              pt.size = 0.1, label = FALSE, label.size = 3.5,  repel = TRUE)
p3 <- DimPlot(seu, reduction = "umap", group.by = "Treatment",
              pt.size = 0.1, label = FALSE, label.size = 3.5,  repel = TRUE)
p4 <- DimPlot(seu, reduction = "umap",
              pt.size = 0.1, label = TRUE, label.size = 3.5,  repel = FALSE)
plot_grid(p1, p2, p3, p4, ncol = 2)
dev.off()

## umap split on response, sample and treatment 
# pdf(paste0("Split_Response_", proj_name, ".pdf"), width = 15, height = 8)
# DimPlot(seu, reduction = "umap", split.by = "Response")
# dev.off()
# pdf(paste0("Split_Sample_", proj_name, ".pdf"), width = 15, height = 8)
# DimPlot(seu, reduction = "umap", split.by = "orig.ident")
# dev.off()
# pdf(paste0("Split_Treatment_", proj_name, ".pdf"), width = 15, height = 8)
# DimPlot(seu, reduction = "umap", split.by = "Treatment")
# dev.off()

########################## tables for QC plots ################################
## make sure you see the correct number of clusters. if not, SetIdent to the
## correct column in  metadata
table(seu@active.ident) 

## save select metadata
u <- seu@meta.data

## table of number of cells 
write.table(u, paste0("cells_metadata_", proj_name, ".txt"), col.names = NA,
            sep = "\t")

## how many cells in each cluster per sample
cell_perC_perO <- table(seu@active.ident, seu@meta.data$orig.ident)
write.table(cell_perC_perO, paste0("cell_per_c_per_O_", proj_name, ".txt"),
            col.names = NA, sep = "\t")

## how many cells in each cluster per treatment
cell_perC_perT <- table(seu@active.ident, seu@meta.data$Treatment)
write.table(cell_perC_perT, paste0("cell_per_c_per_T_", proj_name, ".txt"),
            col.names = NA, sep = "\t")

## how many cells in each cluster per response
cell_perC_perR <- table(seu@active.ident, seu@meta.data$Response)
write.table(cell_perC_perR, paste0("cell_per_c_per_R_", proj_name, ".txt"),
            col.names = NA, sep = "\t")

## Plot number of genes, nUMIs and %mito on a tSNE - 3 plots per sample
## These allows you to see if part of your umap that has a high myeloid content
## (dying cells) or high number of genes (normally, no 2 cells in a droplets
## because you still filtered these bad cells!)
# pdf(paste0("graph_nfeature_", proj_name, ".pdf"), width = 15, height = 8)
# FeaturePlot(seu, features = "nFeature_RNA", cols = c("cadetblue2", "darkred"),
#             pt.size = 0.1, reduction = "umap", label = FALSE)
# dev.off()
# pdf(paste0("graph_nCount_", proj_name, ".pdf"), width = 15, height = 8)
# FeaturePlot(seu, features = "nCount_RNA", cols = c("cadetblue2", "darkred"),
#             pt.size = 0.1, reduction = "umap", label = FALSE)
# dev.off()
# pdf(paste0("graph_percentmt_", proj_name, ".pdf"), width = 15, height = 8)
# FeaturePlot(seu, features = "percent.mt", cols = c("cadetblue2", "darkred"),
#             pt.size = 0.1, reduction = "umap", label = FALSE)
# dev.off()

## export nFeature, nCOunt and %mito per cluster 
# features <- FetchData(seu, vars = c("nFeature_RNA", "nCount_RNA",
#                                     "percent.mt", RESOLUTION), cells = NULL)
# write.table(features, paste0("table_features_", proj_name, ".txt"),
#             col.names = NA, sep = "\t")

## enter as many number as you have clusters and picked the resolution you chose
# f <- ordered(features$RNA_snn_res.0.4, levels = c(0:clust_num))
# f <- ordered(features$RNA_snn_res.0.4,
#              levels = c(0:(max(as.integer(seu@meta.data$RNA_snn_res.0.4)) - 1)))
# pdf(paste0("boxplot_nFeature_v2_", proj_name, ".pdf"), width = 15, height = 8)
# box = boxplot(nFeature_RNA ~ f , data = features, xlab = "# genes",
#               ylab = "Cluster", main = "Number of Genes per cluster",
#               col = "dodgerblue3",  horizontal = TRUE)
# dev.off()
# 
# pdf(paste0("boxplot_nCount_v2_", proj_name, ".pdf"), width = 15, height = 8)
# box = boxplot(nCount_RNA ~ f , data = features, xlab = "# Barcodes", ylab = "Cluster",
#               main = "Number of Barcodes per cluster", col = "dodgerblue3",
#               horizontal = TRUE, ylim = c(0, 60000)) # , ylim = c(0, 35)
# dev.off()
# 
# pdf(paste0("boxplot_percentmito_v2_", proj_name, ".pdf"), width = 15, height = 8)
# box = boxplot(percent.mt ~ f , data = features, xlab = "% mito", ylab = "Cluster",
#               main = "Percentage of mito genes per cluster", col = "dodgerblue3",
#               horizontal = TRUE)
# dev.off()

## save the results
# saveRDS(seu, file = paste0("sub_", proj_name, ".rds"))

#DefaultAssay(seu) <- def_assay
seu_markers <- FindAllMarkers(seu, test.use = "wilcox", min.pct = 0, ###########
                              logfc.threshold = 0, max.cells.per.ident = Inf,
                              latent.vars = NULL, min.cells.feature = 0,
                              min.cells.group = 0)

# seu_markers %>%
#   dplyr::filter(gene != "Ins1",
#                 gene != "Ins2",
#                 gene != "Iapp") -> seu_markers_1

write.table(seu_markers_1, paste0("seu_markers_", proj_name, ".txt"),
            col.names = NA, sep = "\t")

# read.table(paste0("tables/seu_markers_endo4.txt"),
#            header = TRUE, sep = "\t") -> seu_markers

## heatmap on avg log FC
seu_markers_1 %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 5, wt = avg_log2FC) -> topn_fc
pdf(paste0("Heatmap_logFC_", proj_name, ".pdf"), width = 30, height = 10)
DoHeatmap(object = seu, slot = "scale.data", features = topn_fc$gene) +
  #theme(axis.title.x = element_text(angle = 90, hjust = 1)) +
  NoLegend()
dev.off()


saveRDS(seu, file = paste0("macro_recluster_4.rds"))
## save the results
# saveRDS(seu, file = paste0("sub_", proj_name, ".rds"))
# saveRDS(seu, file = paste0("tcell_recluster.rds")) 
## cluster 12 was removed from tcell_recluster to make tcell_recluster_2
# saveRDS(seu, file = paste0("tcell_recluster_10.rds")) ## final file you want

# c("M1 Jun+",
#   "M2 S100a4-",
#   "M2 Zfp36-",
#   "M1 H2+",
#   "Prol. Macrophages", #"M2 Birc5+", ##
#   "Neutrophils", # formerly  "Monocytes",
#   "MAST",
#   "M1 F13a1+", # Formerly "Macrophages",
#   "M2 Sepp1+", ## formerly TAMs
#   "pDCs", "cDCs",
#   "Unknown pDCs" # Formerly, "Neutrophils"
#   ) -> new.cluster.ids

c("M1 Il1a+",
  "M2 Cd63+",
  "MAST",
  "M2 Folr2+",
  "Monocyte Chil3+",
  "M2 Malat1-",
  "M1 Cd72+",
  "cDCs",
  "M2 Birc5+",
  "M2 Sepp1+",
  "pDCs",
  "Unknown pDCs",
  "Monocyte Ifit2bl1+",
  "Neutrophils",
  "Proliferating pDCs"
  ) -> new.cluster.ids

# ## add column in metadata with the "Cluster" title and the cluster name
names(new.cluster.ids) <- levels(seu)
seu <- RenameIdents(seu, new.cluster.ids)
seu <- StashIdent(seu, save.name = "Cluster")
head(seu@meta.data) ## verify new column
#
# pdf(paste0("Annotated_", proj_name, ".pdf"))
# DimPlot(seu, reduction = "umap", label = FALSE, pt.size = 0.1) +
#   labs(title = "Annotated UMAP") + theme(plot.title = element_text(hjust = 0.5))
# dev.off()

pdf(paste0("Annotated_labels_", proj_name, ".pdf")) # the same with labels
DimPlot(seu, reduction = "umap", label = TRUE, pt.size = 0.1, repel = TRUE) +
  labs(title = "Annotated UMAP") + theme(plot.title = element_text(hjust = 0.5))
dev.off()
# 
# ## save the results

########################### subclustering T/NK cells ##3333333333333333333333333
# seu <- readRDS("sub_tcell.rds")
seu <- readRDS(file = paste0("FINAL_latest_ann_tcell.rds"))

seu_can <- subset(seu, subset = Cluster == "T cells CD8+ Naive")
seu_no_can <- subset(seu, subset = Cluster != "T cells CD8+ Naive")

saveRDS(seu_can, file = paste0("sub_tcell_cancer_.rds"))
saveRDS(seu_no_can, file = paste0("sub_tcell_cancer_free.rds"))

seu_can <- subset(seu, subset = Cluster == 13) ## second round at res = 0.6
saveRDS(seu_can, file = paste0("sub_tcell9_cancer_free.rds"))

seu_no_can <- subset(seu, subset = Cluster != 13) ## second round at res = 0.6

seu_can <- subset(seu, subset = Cluster == 12) ## second round at res = 0.6
# saveRDS(seu_can, file = paste0("sub_tcell10_cancer_free.rds"))

seu_no_can <- subset(seu, subset = Cluster != 12) ## second round at res = 0.6

seu_no_can <- subset(seu, subset = Cluster != 6) ## second round at res = 0.6

seu_no_can <- subset(seu, subset = Cluster != 4) ## second round at res = 0.6
# 
# 
# seu_ins <- subset(seu, subset = Ins1 > 1 | Ins2 > 1)
# saveRDS(seu_ins, file = paste0("sub_tcell_ins_cancer_.rds"))
# 
# seu_no_ins <- subset(seu, subset = Ins1 < 3 & Ins2 < 3)
# seu_ins <- subset(seu, subset = Ins1 > 5 | Ins2 > 5)
# 
# saveRDS(seu_ins, file = paste0("sub_tcell6_cancer.rds"))
# saveRDS(seu_no_ins, file = paste0("sub_tcell6_no_cancer.rds"))

# 
# saveRDS(seu_no_ins, file = paste0("sub_tcell_new_cancerfree_.rds"))

# seu_ins <- subset(seu, subset = Ins1 > 1 & Ins2 > 1)

# seu_can <- subset(seu, subset = Cluster == 3)
# seu_can <- subset(seu, subset = Cluster == 12) ## second round at res = 0.6

#seu_no_can <- subset(seu, subset = Cluster != 12)

#seu <- readRDS(file = paste0("FINAL_pancreatic_mice.rds"))
proj_name <- "tcell12" ## name of project
# t1 <- subset(seu, subset = Cluster == "T-cells/NK")
# t1 <- SetIdent(t1, value = "orig.ident")

t1 <- seu_no_can
# t1 <- seu

t1 %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = var_feat) %>%
  # ScaleData(vars.to.regress = c("nCount_RNA", "percent.mt", "orig.ident")) %>%
  ScaleData() %>% 
  RunPCA(npcs = NPC, verbose = FALSE) -> SC

## whatever cutoff for significant PCs, you will use to set sig_dims
## rerun PCA with significant PCs
SC %>% RunPCA(npcs = sig_dims, verbose = FALSE) -> SC

## save scaled
# saveRDS(SC, file = paste0("scaled_", proj_name, ".rds"))
# SC <- readRDS(file = paste0("scaled_", proj_name, ".rds"))

## Clustering and UMAP
## usually NPC = between 15-30 (around 20 is plenty)
## we calculate the clustering using multiple resolutions, the higher, the more
## number of clusters you will have. The point is that it makes biological sense
SC %>%
  FindNeighbors(reduction = "pca", dims = 1:sig_dims, k.param = 30,
                compute.SNN = TRUE, nn.method = "rann", nn.eps = 0) %>% 
  FindClusters(resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                              1.1, 1.2, 1.3, 1.4)) %>% 
  RunUMAP(reduction = "pca", dims = 1:sig_dims) -> FC

## multiple resolution plots
pdf(file = paste0("clusters_diff_res_", proj_name, ".pdf"), width = 15, height = 8)
Seurat_02 <- SetIdent(object = FC, value = 'RNA_snn_res.0.2')
Seurat02 <- RunUMAP(Seurat_02, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat02, reduction = "umap", pt.size = 0.25, label = TRUE) +
  labs(title = "UMAP with res 0.2") + theme(plot.title = element_text(hjust = 0.5))
Seurat_04 <- SetIdent(object = FC, value = 'RNA_snn_res.0.4')
Seurat04 <- RunUMAP(Seurat_04, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat04, reduction = "umap", pt.size = 0.25, label = TRUE) +
  labs(title = "UMAP with res 0.4") + theme(plot.title = element_text(hjust = 0.5))
Seurat_06 <- SetIdent(object = FC, value = 'RNA_snn_res.0.6')
Seurat06 <- RunUMAP(Seurat_06, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat06, reduction = "umap", pt.size = 0.25, label = TRUE) +
  labs(title = "UMAP with res 0.6") + theme(plot.title = element_text(hjust = 0.5))
Seurat_08 <- SetIdent(object = FC, value = 'RNA_snn_res.0.8')
Seurat08 <- RunUMAP(Seurat_08, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat08, reduction = "umap", pt.size = 0.25, label = TRUE) +
  labs(title = "UMAP with res 0.8") + theme(plot.title = element_text(hjust = 0.5))
Seurat_1 <- SetIdent(object = FC, value = 'RNA_snn_res.1.0')
Seurat1 <- RunUMAP(Seurat_1, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat1, reduction = "umap", pt.size = 0.25, label = TRUE) + 
  labs(title = "UMAP with res 1.0") + theme(plot.title = element_text(hjust = 0.5))
Seurat_12 <- SetIdent(object = FC, value = 'RNA_snn_res.1.2')
Seurat12 <- RunUMAP(Seurat_12, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat12, reduction = "umap", pt.size = 0.25, label = TRUE) + 
  labs(title = "UMAP with res 1.2") + theme(plot.title = element_text(hjust = 0.5))
Seurat_14 <- SetIdent(object = FC, value = 'RNA_snn_res.1.4')
Seurat14 <- RunUMAP(Seurat_14, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat14, reduction = "umap", pt.size = 0.25, label = TRUE) +
  labs(title = "UMAP with res 1.4") + theme(plot.title = element_text(hjust = 0.5))
dev.off() 

## write the column name of the chosen resolution with the SetIdent function
## Seurat will set the object with this column
seu <- SetIdent(object = FC, value = RESOLUTION)
# seu <- SetIdent(object = FC, value = 'RNA_snn_res.0.6') ##0.8
seu <- RunUMAP(seu, reduction = "pca", dims = 1:sig_dims)
pdf(file = paste0("res_", RESOLUTION, "_clusters_", proj_name, ".pdf"),
    width = 15, height = 8)
UMAPPlot(seu, reduction = "umap", pt.size = 0.1, label = TRUE,
         label.size = 5) + NoLegend() +
  labs(title = paste0("UMAP with ", RESOLUTION)) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

## produce dimplots
pdf(paste0("umap_", proj_name, ".pdf"), width = 15, height = 8)
p1 <- DimPlot(seu, reduction = "umap", group.by = "orig.ident",
              pt.size = 0.1, label = FALSE, label.size = 2.0)
p2 <- DimPlot(seu, reduction = "umap", group.by = "Response",
              pt.size = 0.1, label = FALSE, label.size = 3.5,  repel = TRUE)
p3 <- DimPlot(seu, reduction = "umap", group.by = "Treatment",
              pt.size = 0.1, label = FALSE, label.size = 3.5,  repel = TRUE)
p4 <- DimPlot(seu, reduction = "umap",
              pt.size = 0.1, label = TRUE, label.size = 3.5,  repel = FALSE)
plot_grid(p1, p2, p3, p4, ncol = 2)
dev.off()

## umap split on response, sample and treatment 
pdf(paste0("Split_Response_", proj_name, ".pdf"), width = 15, height = 8)
DimPlot(seu, reduction = "umap", split.by = "Response")
dev.off()
pdf(paste0("Split_Sample_", proj_name, ".pdf"), width = 15, height = 8)
DimPlot(seu, reduction = "umap", split.by = "orig.ident")
dev.off()
pdf(paste0("Split_Treatment_", proj_name, ".pdf"), width = 15, height = 8)
DimPlot(seu, reduction = "umap", split.by = "Treatment")
dev.off()

########################## tables for QC plots ################################
## make sure you see the correct number of clusters. if not, SetIdent to the
## correct column in  metadata
table(seu@active.ident) 

## save select metadata
u <- seu@meta.data

## table of number of cells 
write.table(u, paste0("cells_metadata_", proj_name, ".txt"), col.names = NA,
            sep = "\t")

## how many cells in each cluster per sample
cell_perC_perO <- table(seu@active.ident, seu@meta.data$orig.ident)
write.table(cell_perC_perO, paste0("cell_per_c_per_O_", proj_name, ".txt"),
            col.names = NA, sep = "\t")

## how many cells in each cluster per treatment
cell_perC_perT <- table(seu@active.ident, seu@meta.data$Treatment)
write.table(cell_perC_perT, paste0("cell_per_c_per_T_", proj_name, ".txt"),
            col.names = NA, sep = "\t")

## how many cells in each cluster per response
cell_perC_perR <- table(seu@active.ident, seu@meta.data$Response)
write.table(cell_perC_perR, paste0("cell_per_c_per_R_", proj_name, ".txt"),
            col.names = NA, sep = "\t")

## Plot number of genes, nUMIs and %mito on a tSNE - 3 plots per sample
## These allows you to see if part of your umap that has a high myeloid content
## (dying cells) or high number of genes (normally, no 2 cells in a droplets
## because you still filtered these bad cells!)
# pdf(paste0("graph_nfeature_", proj_name, ".pdf"), width = 15, height = 8)
# FeaturePlot(seu, features = "nFeature_RNA", cols = c("cadetblue2", "darkred"),
#             pt.size = 0.1, reduction = "umap", label = FALSE)
# dev.off()
# pdf(paste0("graph_nCount_", proj_name, ".pdf"), width = 15, height = 8)
# FeaturePlot(seu, features = "nCount_RNA", cols = c("cadetblue2", "darkred"),
#             pt.size = 0.1, reduction = "umap", label = FALSE)
# dev.off()
# pdf(paste0("graph_percentmt_", proj_name, ".pdf"), width = 15, height = 8)
# FeaturePlot(seu, features = "percent.mt", cols = c("cadetblue2", "darkred"),
#             pt.size = 0.1, reduction = "umap", label = FALSE)
# dev.off()
# 
# ## export nFeature, nCOunt and %mito per cluster 
# features <- FetchData(seu, vars = c("nFeature_RNA", "nCount_RNA",
#                                     "percent.mt", RESOLUTION), cells = NULL)
# write.table(features, paste0("table_features_", proj_name, ".txt"),
#             col.names = NA, sep = "\t")
# 
# ## enter as many number as you have clusters and picked the resolution you chose
# # f <- ordered(features$RNA_snn_res.0.4, levels = c(0:clust_num))
# f <- ordered(features$RNA_snn_res.0.4,
#              levels = c(0:(max(as.integer(seu@meta.data$RNA_snn_res.0.4)) - 1)))
# pdf(paste0("boxplot_nFeature_v2_", proj_name, ".pdf"), width = 15, height = 8)
# box = boxplot(nFeature_RNA ~ f , data = features, xlab = "# genes",
#               ylab = "Cluster", main = "Number of Genes per cluster",
#               col = "dodgerblue3",  horizontal = TRUE)
# dev.off()
# 
# pdf(paste0("boxplot_nCount_v2_", proj_name, ".pdf"), width = 15, height = 8)
# box = boxplot(nCount_RNA ~ f , data = features, xlab = "# Barcodes", ylab = "Cluster",
#               main = "Number of Barcodes per cluster", col = "dodgerblue3",
#               horizontal = TRUE, ylim = c(0, 60000)) # , ylim = c(0, 35)
# dev.off()
# 
# pdf(paste0("boxplot_percentmito_v2_", proj_name, ".pdf"), width = 15, height = 8)
# box = boxplot(percent.mt ~ f , data = features, xlab = "% mito", ylab = "Cluster",
#               main = "Percentage of mito genes per cluster", col = "dodgerblue3",
#               horizontal = TRUE)
# dev.off()

#DefaultAssay(seu) <- def_assay
seu_markers <- FindAllMarkers(seu, test.use = "wilcox", min.pct = 0, ###########
                              logfc.threshold = 0, max.cells.per.ident = Inf,
                              latent.vars = NULL, min.cells.feature = 0,
                              min.cells.group = 0)
## 
seu_markers %>%
  dplyr::filter(gene != "Ins1",
                gene != "Ins2",
                gene != "Iapp") -> seu_markers_1

write.table(seu_markers_1, paste0("seu_markers_", proj_name, ".txt"),
            col.names = NA, sep = "\t")

## heatmap on avg log FC
seu_markers_1 %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 5, wt = avg_log2FC) -> topn_fc
pdf(paste0("Heatmap_logFC_", proj_name, ".pdf"), width = 30, height = 10)
DoHeatmap(object = seu, slot = "scale.data", features = topn_fc$gene) +
  #theme(axis.title.x = element_text(angle = 90, hjust = 1)) +
  NoLegend()
dev.off()
# 

saveRDS(seu, file = paste0("tcell_recluster_12.rds"))
## save the results
# saveRDS(seu, file = paste0("sub_", proj_name, ".rds"))
# saveRDS(seu, file = paste0("tcell_recluster.rds")) 
## cluster 12 was removed from tcell_recluster to make tcell_recluster_2
# saveRDS(seu, file = paste0("tcell_recluster_10.rds")) ## final file you want


c("Cd8+ Ex Ccl5+", "Cd8+ Naive Dapl1+", "Cd4+ Naive Igfbp4+", "Cd4+ Naive Satb1+",
  "Tregs", "Cd8+ Ex Ifit3+", "Gamma Delta T", "Cd8+ Naive Npm1+",
  "Unknown T-cells", # Formerly "Th1",
  "Proliferating T-cells", # Formerly "TILs", 
  "Cd8+ Activated", # Formerly "Th2",
  "NK", "NKT",
  "Cd4+ Th17+" # Formerly "Th17"
  ) -> new.cluster.ids

# ## add column in metadata with the "Cluster" title and the cluster name
names(new.cluster.ids) <- levels(seu)
seu <- RenameIdents(seu, new.cluster.ids)
seu <- StashIdent(seu, save.name = "Cluster")
head(seu@meta.data) ## verify new column
#
# pdf(paste0("Annotated_", proj_name, ".pdf"))
# DimPlot(seu, reduction = "umap", label = FALSE, pt.size = 0.1) +
#   labs(title = "Annotated UMAP") + theme(plot.title = element_text(hjust = 0.5))
# dev.off()

pdf(paste0("Annotated_labels_", proj_name, ".pdf")) # the same with labels
DimPlot(seu, reduction = "umap", label = TRUE, pt.size = 0.1, repel = TRUE) +
  labs(title = "Annotated UMAP") + theme(plot.title = element_text(hjust = 0.5))
dev.off()
# 
# ## save the results
# saveRDS(seu, file = paste0("final_ann_endo.rds"))

########################### subclustering Endothelial #4444444444444444444444444
# seu <- readRDS(file = paste0("FINAL_pancreatic_mice.rds"))
# seu <- readRDS(file = paste0("FINAL_latest_ann_2_pancreatic_mice.rds"))

seu <- readRDS(paste0("/media/seq-srv-05/vrc/Project/Project_Theo/",
                      "pancreatic_mice/raw_1000dr/endo_sub/sub_endo.rds"))

seu_ins <- subset(seu, subset = Ins1 > 1 & Ins2 > 1 & Ppy > 1)
seu_no_ins <- seu[, !colnames(seu) %in% colnames(seu_ins)]
seu_aci <- subset(seu_no_ins, subset =  Cela2a > 1 & Cpa1 > 1)
seu_no_can <- seu_no_ins[, !colnames(seu_no_ins) %in% colnames(seu_aci)]

saveRDS(seu_ins, file = paste0("sub_endo_cancer.rds"))
saveRDS(seu_aci, file = paste0("sub_endo_acinar.rds"))
saveRDS(seu_no_can, file = paste0("sub_endo_cancer_removed.rds"))

## we removed cluster #3 because it was highly expressing cancer markers
# seu_3 <- subset(seu, subset = RNA_snn_res.0.4 == 3)
# saveRDS(seu_3, file = paste0("sub_endo_cluster3.rds"))

# saveRDS(seu_ins, file = paste0("sub_endo_leaky_cancer_cells.rds"))
# saveRDS(seu_aci, file = paste0("sub_endo_leaky_acinar_cells.rds"))
# saveRDS(seu_no_can, file = paste0("sub_endo_cancer_removed.rds"))
# saveRDS(seu, file = paste0("sub_endo_some_cancer_removed.rds"))

# seu_not_3 <- subset(seu, subset = RNA_snn_res.0.4 != 3)
# e1 <- seu_not_3
e1 <- seu_no_can
proj_name <- "endo4" ## name of project
# e1 <- subset(seu, subset = Cluster == "Endothelial cells")
# e1 <- SetIdent(e1, value = "orig.ident")

e1 %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = var_feat) %>%
  # ScaleData(vars.to.regress = c("nCount_RNA", "percent.mt", "orig.ident")) %>%
  ScaleData() %>% 
  RunPCA(npcs = NPC, verbose = FALSE) -> SC

## whatever cutoff for significant PCs, you will use to set sig_dims
## rerun PCA with significant PCs
SC %>% RunPCA(npcs = sig_dims, verbose = FALSE) -> SC

## save scaled
# saveRDS(SC, file = paste0("scaled_", proj_name, ".rds"))
# SC <- readRDS(file = paste0("scaled_", proj_name, ".rds"))

## Clustering and UMAP
## usually NPC = between 15-30 (around 20 is plenty)
## we calculate the clustering using multiple resolutions, the higher, the more
## number of clusters you will have. The point is that it makes biological sense
SC %>%
  FindNeighbors(reduction = "pca", dims = 1:sig_dims, k.param = 30,
                compute.SNN = TRUE, nn.method = "rann", nn.eps = 0) %>% 
  FindClusters(resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                              1.1, 1.2, 1.3, 1.4)) %>% 
  RunUMAP(reduction = "pca", dims = 1:sig_dims) -> FC

# multiple resolution plots
pdf(file = paste0("clusters_diff_res_", proj_name, ".pdf"), width = 15, height = 8)
Seurat_02 <- SetIdent(object = FC, value = 'RNA_snn_res.0.2')
Seurat02 <- RunUMAP(Seurat_02, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat02, reduction = "umap", pt.size = 0.25, label = TRUE) +
  labs(title = "UMAP with res 0.2") + theme(plot.title = element_text(hjust = 0.5))
Seurat_04 <- SetIdent(object = FC, value = 'RNA_snn_res.0.4')
Seurat04 <- RunUMAP(Seurat_04, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat04, reduction = "umap", pt.size = 0.25, label = TRUE) +
  labs(title = "UMAP with res 0.4") + theme(plot.title = element_text(hjust = 0.5))
Seurat_06 <- SetIdent(object = FC, value = 'RNA_snn_res.0.6')
Seurat06 <- RunUMAP(Seurat_06, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat06, reduction = "umap", pt.size = 0.25, label = TRUE) +
  labs(title = "UMAP with res 0.6") + theme(plot.title = element_text(hjust = 0.5))
Seurat_08 <- SetIdent(object = FC, value = 'RNA_snn_res.0.8')
Seurat08 <- RunUMAP(Seurat_08, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat08, reduction = "umap", pt.size = 0.25, label = TRUE) +
  labs(title = "UMAP with res 0.8") + theme(plot.title = element_text(hjust = 0.5))
Seurat_1 <- SetIdent(object = FC, value = 'RNA_snn_res.1.0')
Seurat1 <- RunUMAP(Seurat_1, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat1, reduction = "umap", pt.size = 0.25, label = TRUE) +
  labs(title = "UMAP with res 1.0") + theme(plot.title = element_text(hjust = 0.5))
Seurat_12 <- SetIdent(object = FC, value = 'RNA_snn_res.1.2')
Seurat12 <- RunUMAP(Seurat_12, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat12, reduction = "umap", pt.size = 0.25, label = TRUE) +
  labs(title = "UMAP with res 1.2") + theme(plot.title = element_text(hjust = 0.5))
Seurat_14 <- SetIdent(object = FC, value = 'RNA_snn_res.1.4')
Seurat14 <- RunUMAP(Seurat_14, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat14, reduction = "umap", pt.size = 0.25, label = TRUE) +
  labs(title = "UMAP with res 1.4") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

## write the column name of the chosen resolution with the SetIdent function
## Seurat will set the object with this column
# seu <- SetIdent(object = FC, value = RESOLUTION)
seu <- SetIdent(object = FC, value = 'RNA_snn_res.0.8')
seu <- RunUMAP(seu, reduction = "pca", dims = 1:sig_dims)
pdf(file = paste0("res_", RESOLUTION, "_clusters_", proj_name, ".pdf"),
    width = 15, height = 8)
UMAPPlot(seu, reduction = "umap", pt.size = 0.1, label = TRUE,
         label.size = 5) + NoLegend() +
  labs(title = paste0("UMAP with ", RESOLUTION)) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

## produce dimplots
pdf(paste0("umap_", proj_name, ".pdf"), width = 15, height = 8)
p1 <- DimPlot(seu, reduction = "umap", group.by = "orig.ident",
              pt.size = 0.1, label = FALSE, label.size = 2.0)
p2 <- DimPlot(seu, reduction = "umap", group.by = "Response",
              pt.size = 0.1, label = FALSE, label.size = 3.5,  repel = TRUE)
p3 <- DimPlot(seu, reduction = "umap", group.by = "Treatment",
              pt.size = 0.1, label = FALSE, label.size = 3.5,  repel = TRUE)
p4 <- DimPlot(seu, reduction = "umap",
              pt.size = 0.1, label = TRUE, label.size = 3.5,  repel = FALSE)
plot_grid(p1, p2, p3, p4, ncol = 2)
dev.off()

## umap split on response, sample and treatment 
pdf(paste0("Split_Response_", proj_name, ".pdf"), width = 15, height = 8)
DimPlot(seu, reduction = "umap", split.by = "Response")
dev.off()
pdf(paste0("Split_Sample_", proj_name, ".pdf"), width = 15, height = 8)
DimPlot(seu, reduction = "umap", split.by = "orig.ident")
dev.off()
pdf(paste0("Split_Treatment_", proj_name, ".pdf"), width = 15, height = 8)
DimPlot(seu, reduction = "umap", split.by = "Treatment")
dev.off()

########################## tables for QC plots ################################
## make sure you see the correct number of clusters. if not, SetIdent to the
## correct column in  metadata
table(seu@active.ident) 

## save select metadata
u <- seu@meta.data

## table of number of cells 
write.table(u, paste0("cells_metadata_", proj_name, ".txt"), col.names = NA,
            sep = "\t")

## how many cells in each cluster per sample
cell_perC_perO <- table(seu@active.ident, seu@meta.data$orig.ident)
write.table(cell_perC_perO, paste0("cell_per_c_per_O_", proj_name, ".txt"),
            col.names = NA, sep = "\t")

## how many cells in each cluster per treatment
cell_perC_perT <- table(seu@active.ident, seu@meta.data$Treatment)
write.table(cell_perC_perT, paste0("cell_per_c_per_T_", proj_name, ".txt"),
            col.names = NA, sep = "\t")

## how many cells in each cluster per response
cell_perC_perR <- table(seu@active.ident, seu@meta.data$Response)
write.table(cell_perC_perR, paste0("cell_per_c_per_R_", proj_name, ".txt"),
            col.names = NA, sep = "\t")

## Plot number of genes, nUMIs and %mito on a tSNE - 3 plots per sample
## These allows you to see if part of your umap that has a high myeloid content
## (dying cells) or high number of genes (normally, no 2 cells in a droplets
## because you still filtered these bad cells!)
pdf(paste0("graph_nfeature_", proj_name, ".pdf"), width = 15, height = 8)
FeaturePlot(seu, features = "nFeature_RNA", cols = c("cadetblue2", "darkred"),
            pt.size = 0.1, reduction = "umap", label = FALSE)
dev.off()
pdf(paste0("graph_nCount_", proj_name, ".pdf"), width = 15, height = 8)
FeaturePlot(seu, features = "nCount_RNA", cols = c("cadetblue2", "darkred"),
            pt.size = 0.1, reduction = "umap", label = FALSE)
dev.off()
pdf(paste0("graph_percentmt_", proj_name, ".pdf"), width = 15, height = 8)
FeaturePlot(seu, features = "percent.mt", cols = c("cadetblue2", "darkred"),
            pt.size = 0.1, reduction = "umap", label = FALSE)
dev.off()

## export nFeature, nCOunt and %mito per cluster 
features <- FetchData(seu, vars = c("nFeature_RNA", "nCount_RNA",
                                    "percent.mt", RESOLUTION), cells = NULL)
write.table(features, paste0("table_features_", proj_name, ".txt"),
            col.names = NA, sep = "\t")

## enter as many number as you have clusters and picked the resolution you chose
# f <- ordered(features$RNA_snn_res.0.4, levels = c(0:clust_num))
f <- ordered(features$RNA_snn_res.0.4,
             levels = c(0:(max(as.integer(seu@meta.data$RNA_snn_res.0.4)) - 1)))
pdf(paste0("boxplot_nFeature_v2_", proj_name, ".pdf"), width = 15, height = 8)
box = boxplot(nFeature_RNA ~ f , data = features, xlab = "# genes",
              ylab = "Cluster", main = "Number of Genes per cluster",
              col = "dodgerblue3",  horizontal = TRUE)
dev.off()

pdf(paste0("boxplot_nCount_v2_", proj_name, ".pdf"), width = 15, height = 8)
box = boxplot(nCount_RNA ~ f , data = features, xlab = "# Barcodes", ylab = "Cluster",
              main = "Number of Barcodes per cluster", col = "dodgerblue3",
              horizontal = TRUE, ylim = c(0, 60000)) # , ylim = c(0, 35)
dev.off()

pdf(paste0("boxplot_percentmito_v2_", proj_name, ".pdf"), width = 15, height = 8)
box = boxplot(percent.mt ~ f , data = features, xlab = "% mito", ylab = "Cluster",
              main = "Percentage of mito genes per cluster", col = "dodgerblue3",
              horizontal = TRUE)
dev.off()

## save the results
#saveRDS(seu, file = paste0("sub_endo.rds"))

# seu <- readRDS(file = paste0("endo_sub/sub_endo.rds"))

#    0    1    2    3    4    5    6    7    8    9   10 
# 1495 1359  900  857  712  554  503  341  302   93   42 

c("Cap arterial", "Capillaries Plpp3+", "Tip cells",
  "Capillaries Fmo2+",
  "Capillaries Xist+",
  "Capillaries Prkcdbp+",
  "Pericytes", "Venous", "Cap Interferon", "Artery",
  "Postcap venule",  # "Capillaries Rac2+",
  "Lymphatic") -> new.cluster.ids

## add column in metadata with the "Cluster" title and the cluster name
names(new.cluster.ids) <- levels(seu)
seu <- RenameIdents(seu, new.cluster.ids)
seu <- StashIdent(seu, save.name = "Cluster")
head(seu@meta.data) ## verify new column
# 
# pdf(paste0("Annotated_v2", proj_name, ".pdf"))
# DimPlot(seu, reduction = "umap", label = FALSE, pt.size = 0.1) +
#   labs(title = "Annotated UMAP") + theme(plot.title = element_text(hjust = 0.5))
# dev.off()

pdf(paste0("Annotated_labels_", proj_name, ".pdf")) # the same with labels
DimPlot(seu, reduction = "umap", label = TRUE, pt.size = 0.1, repel = TRUE) +
  labs(title = "Annotated UMAP") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

## save the results
saveRDS(seu, file = paste0("final_ann_endo.rds"))

############################### sub_endo processing ############################
# seu <- readRDS(paste0("/media/seq-srv-05/vrc/Project/Project_Theo/",
#                       "pancreatic_mice/raw_1000dr/endo_sub/sub_endo.rds"))

# sub_endo <- readRDS(paste0("~/Documents/tmp/pancreatic_melanie/results/raw1000dr/",
#                            "subclustering/endo_sub/sub_endo.rds"))
#levels(sub_endo)
# [1] "0"  "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11"
# tcells_actual_FINAL <- readRDS("~/Documents/tmp/pancreatic_melanie/tcells_actual_FINAL.rds")
# macro_recluster_4 <- readRDS("~/Documents/tmp/pancreatic_melanie/macro_recluster_4.rds")
# seu <- macro_recluster_4
seu <- readRDS(file = paste0("FINAL_cancer_ann_pancreatic_mice.rds"))
b_cells <- subset(seu, subset = Cluster == "B-Cell") 

## clean up metadata
# sub <- seu
sub <- b_cells
sub <- StashIdent(sub, save.name = "Cluster")
k1 <- sub@meta.data$Cluster
md1 <- as.data.frame(sub@meta.data)
md1 %>%
  tibble::rownames_to_column() %>%
  dplyr::rename(cells = rowname) %>%
  dplyr::select(cells, Cluster, orig.ident, nCount_RNA, nFeature_RNA, Response,
                Treatment, percent.mt) %>%
  as.data.frame() -> a1
rownames(a1) <- a1$cells
sub@meta.data <- a1
#View(sub@meta.data)
seu0 <- sub
seu0 <- SetIdent(object = seu0, value = "Cluster")
# saveRDS(seu, file = "final_endo_clean_metadata.rds")
# saveRDS(seu, file = "final_tcell_clean_metadata.rds")
# saveRDS(seu, file = "final_macro_clean_metadata.rds")
saveRDS(seu0, file = "final_bcell_clean_metadata.rds")

############################### subcluster neutrophils #########################
merged_mm <- readRDS("macro_sub2/sub_macro.rds")
proj_name <- "neutrophils3"
# seu <- SetIdent(object = seu, value = 'RNA_snn_res.0.4')
# seu <- SetIdent(object = seu, value = "Cluster")
# neutro1 <- subset(seu, subset = RNA_snn_res.0.4 == 11) ## 271 samples
seu <- SetIdent(object = seu, value = 'RNA_snn_res.1.4')
seu <- SetIdent(object = seu, value = "Cluster")
neutro2 <- subset(seu, subset = RNA_snn_res.1.4 == 50) ## 284 samples

# file_path <- "/media/seq-srv-05/vrc/Project/Project_Theo/pancreatic_mice/unfiltered"
# seu1 <- readRDS(paste0(file_path, "FINAL_pancreatic_mice.rds"))
# m1 <- subset(seu1, subset = Cluster == "Macrophages/Myeloids")
# m2 <- subset(seu1, subset = Cluster == "Monocytes/Myeloids")
# m3 <- subset(seu1, subset = Cluster == "Granulocytes/Neutrophils")
# m1 %>% merge(y = c(m2, m3), project = proj_name) -> merged_myeloids
# proj_name <- "neutrophils"
# 
# merged_myeloids %>% 
#   FindVariableFeatures(selection.method = "vst", nfeatures = var_feat) %>%
#   ScaleData() %>% 
#   RunPCA(npcs = sig_dims, verbose = FALSE) %>%
#   FindNeighbors(reduction = "pca", dims = 1:sig_dims, k.param = 30,
#                 compute.SNN = TRUE, nn.method = "rann", nn.eps = 0) %>% 
#   FindClusters(resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
#                               1.1, 1.2, 1.3, 1.4)) %>% 
#   RunUMAP(reduction = "pca", dims = 1:sig_dims) -> FC
# 
## multiple resolution plots
# pdf(file = paste0("clusters_diff_res_", proj_name, ".pdf"), width = 15, height = 8)
# Seurat_02 <- SetIdent(object = seu, value = 'RNA_snn_res.0.2')
# Seurat02 <- RunUMAP(Seurat_02, reduction = "pca", dims = 1:sig_dims)
# UMAPPlot(Seurat02, reduction = "umap", pt.size = 0.25, label = TRUE) +
#   labs(title = "UMAP with res 0.2") + theme(plot.title = element_text(hjust = 0.5))
# Seurat_04 <- SetIdent(object = seu, value = 'RNA_snn_res.0.4')
# Seurat04 <- RunUMAP(Seurat_04, reduction = "pca", dims = 1:sig_dims)
# UMAPPlot(Seurat04, reduction = "umap", pt.size = 0.25, label = TRUE) +
#   labs(title = "UMAP with res 0.4") + theme(plot.title = element_text(hjust = 0.5))
# Seurat_06 <- SetIdent(object = seu, value = 'RNA_snn_res.0.6')
# Seurat06 <- RunUMAP(Seurat_06, reduction = "pca", dims = 1:sig_dims)
# UMAPPlot(Seurat06, reduction = "umap", pt.size = 0.25, label = TRUE) +
#   labs(title = "UMAP with res 0.6") + theme(plot.title = element_text(hjust = 0.5))
# Seurat_08 <- SetIdent(object = seu, value = 'RNA_snn_res.0.8')
# Seurat08 <- RunUMAP(Seurat_08, reduction = "pca", dims = 1:sig_dims)
# UMAPPlot(Seurat08, reduction = "umap", pt.size = 0.25, label = TRUE) +
#   labs(title = "UMAP with res 0.8") + theme(plot.title = element_text(hjust = 0.5))
# Seurat_1 <- SetIdent(object = seu, value = 'RNA_snn_res.1.0')
# Seurat1 <- RunUMAP(Seurat_1, reduction = "pca", dims = 1:sig_dims)
# UMAPPlot(Seurat1, reduction = "umap", pt.size = 0.25, label = TRUE) +
#   labs(title = "UMAP with res 1.0") + theme(plot.title = element_text(hjust = 0.5))
# Seurat_12 <- SetIdent(object = seu, value = 'RNA_snn_res.1.2')
# Seurat12 <- RunUMAP(Seurat_12, reduction = "pca", dims = 1:sig_dims)
# UMAPPlot(Seurat12, reduction = "umap", pt.size = 0.25, label = TRUE) +
#   labs(title = "UMAP with res 1.2") + theme(plot.title = element_text(hjust = 0.5))
# Seurat_14 <- SetIdent(object = seu, value = 'RNA_snn_res.1.4')
# Seurat14 <- RunUMAP(Seurat_14, reduction = "pca", dims = 1:sig_dims)
# UMAPPlot(Seurat14, reduction = "umap", pt.size = 0.25, label = TRUE) +
#   labs(title = "UMAP with res 1.4") + theme(plot.title = element_text(hjust = 0.5))
# dev.off()
# 
# ## write the column name of the chosen resolution with the SetIdent function
# ## Seurat will set the object with this column
# seu <- SetIdent(object = FC, value = RESOLUTION)
# # seu <- SetIdent(object = FC, value = 'RNA_snn_res.0.6')
# seu <- RunUMAP(seu, reduction = "pca", dims = 1:sig_dims)
# pdf(file = paste0("res_", RESOLUTION, "_clusters_", proj_name, ".pdf"),
#     width = 15, height = 8)
# UMAPPlot(seu, reduction = "umap", pt.size = 0.1, label = TRUE,
#          label.size = 5) + NoLegend() +
#   labs(title = paste0("UMAP with ", RESOLUTION)) +
#   theme(plot.title = element_text(hjust = 0.5))
# dev.off()
# 
# produce dimplots
# pdf(paste0("umap_", proj_name, ".pdf"), width = 15, height = 8)
# p1 <- DimPlot(seu, reduction = "umap", group.by = "orig.ident",
#               pt.size = 0.1, label = FALSE, label.size = 2.0)
# p2 <- DimPlot(seu, reduction = "umap", group.by = "Response",
#               pt.size = 0.1, label = FALSE, label.size = 3.5,  repel = TRUE)
# p3 <- DimPlot(seu, reduction = "umap", group.by = "Treatment",
#               pt.size = 0.1, label = FALSE, label.size = 3.5,  repel = TRUE)
# p4 <- DimPlot(seu, reduction = "umap",
#               pt.size = 0.1, label = TRUE, label.size = 3.5,  repel = FALSE)
# plot_grid(p1, p2, p3, p4, ncol = 2)
# dev.off()
# 
# ## umap split on response, sample and treatment 
# pdf(paste0("Split_Response_", proj_name, ".pdf"), width = 15, height = 8)
# DimPlot(seu, reduction = "umap", split.by = "Response")
# dev.off()
# pdf(paste0("Split_Sample_", proj_name, ".pdf"), width = 15, height = 8)
# DimPlot(seu, reduction = "umap", split.by = "orig.ident")
# dev.off()
# pdf(paste0("Split_Treatment_", proj_name, ".pdf"), width = 15, height = 8)
# DimPlot(seu, reduction = "umap", split.by = "Treatment")
# dev.off()
# 
# ########################## tables for QC plots ################################
# ## make sure you see the correct number of clusters. if not, SetIdent to the
# ## correct column in  metadata
# # table(seu@active.ident) 
# 
# ## save select metadata
# u <- seu@meta.data
# 
# ## table of number of cells 
# write.table(u, paste0("cells_metadata_", proj_name, ".txt"), col.names = NA,
#             sep = "\t")
# 
# ## how many cells in each cluster per sample
# cell_perC_perO <- table(seu@active.ident, seu@meta.data$orig.ident)
# write.table(cell_perC_perO, paste0("cell_per_c_per_O_", proj_name, ".txt"),
#             col.names = NA, sep = "\t")
# 
# ## how many cells in each cluster per treatment
# cell_perC_perT <- table(seu@active.ident, seu@meta.data$Treatment)
# write.table(cell_perC_perT, paste0("cell_per_c_per_T_", proj_name, ".txt"),
#             col.names = NA, sep = "\t")
# 
# ## how many cells in each cluster per response
# cell_perC_perR <- table(seu@active.ident, seu@meta.data$Response)
# write.table(cell_perC_perR, paste0("cell_per_c_per_R_", proj_name, ".txt"),
#             col.names = NA, sep = "\t")
# 
# ## Plot number of genes, nUMIs and %mito on a tSNE - 3 plots per sample
# ## These allows you to see if part of your umap that has a high myeloid content
# ## (dying cells) or high number of genes (normally, no 2 cells in a droplets
# ## because you still filtered these bad cells!)
# pdf(paste0("graph_nfeature_", proj_name, ".pdf"), width = 15, height = 8)
# FeaturePlot(seu, features = "nFeature_RNA", cols = c("cadetblue2", "darkred"),
#             pt.size = 0.1, reduction = "umap", label = FALSE)
# dev.off()
# pdf(paste0("graph_nCount_", proj_name, ".pdf"), width = 15, height = 8)
# FeaturePlot(seu, features = "nCount_RNA", cols = c("cadetblue2", "darkred"),
#             pt.size = 0.1, reduction = "umap", label = FALSE)
# dev.off()
# pdf(paste0("graph_percentmt_", proj_name, ".pdf"), width = 15, height = 8)
# FeaturePlot(seu, features = "percent.mt", cols = c("cadetblue2", "darkred"),
#             pt.size = 0.1, reduction = "umap", label = FALSE)
# dev.off()

proj_name <- "neutro_final2"
seu <- readRDS("neutro1.rds")
seu@meta.data$treatment_response <- paste0(seu@meta.data$Treatment, "_", seu@meta.data$Response)
seu <- SetIdent(object = seu, value = "treatment_response")


################################## B-cell ######################################
# seu <- readRDS("FINAL_cancer_ann_pancreatic_mice.rds")
# b_cells <- subset(seu, subset = Cluster == "B-Cell")
seu_ins <- subset(b_cells, subset = Ins1 > 1 & Ins2 > 1)
seu_no_ins <- subset(b_cells, subset = Ins1 < 1 | Ins2 < 1)
# seu <- b_cells
seu <- seu_no_ins
seu %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = var_feat) %>%
  ScaleData() %>%
  RunPCA(npcs = sig_dims, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:sig_dims, k.param = 30,
                compute.SNN = TRUE, nn.method = "rann", nn.eps = 0) %>%
  FindClusters(resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                              1.1, 1.2, 1.3, 1.4)) %>%
  RunUMAP(reduction = "pca", dims = 1:sig_dims) -> FC

# multiple resolution plots
pdf(file = paste0("clusters_diff_res_", proj_name, ".pdf"), width = 15, height = 8)
Seurat_02 <- SetIdent(object = seu, value = 'RNA_snn_res.0.2')
Seurat02 <- RunUMAP(Seurat_02, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat02, reduction = "umap", pt.size = 0.25, label = TRUE) +
  labs(title = "UMAP with res 0.2") + theme(plot.title = element_text(hjust = 0.5))
Seurat_04 <- SetIdent(object = seu, value = 'RNA_snn_res.0.4')
Seurat04 <- RunUMAP(Seurat_04, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat04, reduction = "umap", pt.size = 0.25, label = TRUE) +
  labs(title = "UMAP with res 0.4") + theme(plot.title = element_text(hjust = 0.5))
Seurat_06 <- SetIdent(object = seu, value = 'RNA_snn_res.0.6')
Seurat06 <- RunUMAP(Seurat_06, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat06, reduction = "umap", pt.size = 0.25, label = TRUE) +
  labs(title = "UMAP with res 0.6") + theme(plot.title = element_text(hjust = 0.5))
Seurat_08 <- SetIdent(object = seu, value = 'RNA_snn_res.0.8')
Seurat08 <- RunUMAP(Seurat_08, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat08, reduction = "umap", pt.size = 0.25, label = TRUE) +
  labs(title = "UMAP with res 0.8") + theme(plot.title = element_text(hjust = 0.5))
Seurat_1 <- SetIdent(object = seu, value = 'RNA_snn_res.1.0')
Seurat1 <- RunUMAP(Seurat_1, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat1, reduction = "umap", pt.size = 0.25, label = TRUE) +
  labs(title = "UMAP with res 1.0") + theme(plot.title = element_text(hjust = 0.5))
Seurat_12 <- SetIdent(object = seu, value = 'RNA_snn_res.1.2')
Seurat12 <- RunUMAP(Seurat_12, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat12, reduction = "umap", pt.size = 0.25, label = TRUE) +
  labs(title = "UMAP with res 1.2") + theme(plot.title = element_text(hjust = 0.5))
Seurat_14 <- SetIdent(object = seu, value = 'RNA_snn_res.1.4')
Seurat14 <- RunUMAP(Seurat_14, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat14, reduction = "umap", pt.size = 0.25, label = TRUE) +
  labs(title = "UMAP with res 1.4") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

## write the column name of the chosen resolution with the SetIdent function
## Seurat will set the object with this column
seu <- SetIdent(object = FC, value = RESOLUTION)
# seu <- SetIdent(object = FC, value = 'RNA_snn_res.0.6')
seu <- RunUMAP(seu, reduction = "pca", dims = 1:sig_dims)
pdf(file = paste0("res_", RESOLUTION, "_clusters_", proj_name, ".pdf"),
    width = 15, height = 8)
UMAPPlot(seu, reduction = "umap", pt.size = 0.1, label = TRUE,
         label.size = 5) + NoLegend() +
  labs(title = paste0("UMAP with ", RESOLUTION)) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

## produce dimplots
pdf(paste0("umap_", proj_name, ".pdf"), width = 15, height = 8)
p1 <- DimPlot(seu, reduction = "umap", group.by = "orig.ident",
              pt.size = 0.1, label = FALSE, label.size = 2.0)
p2 <- DimPlot(seu, reduction = "umap", group.by = "Response",
              pt.size = 0.1, label = FALSE, label.size = 3.5,  repel = TRUE)
p3 <- DimPlot(seu, reduction = "umap", group.by = "Treatment",
              pt.size = 0.1, label = FALSE, label.size = 3.5,  repel = TRUE)
p4 <- DimPlot(seu, reduction = "umap",
              pt.size = 0.1, label = TRUE, label.size = 3.5,  repel = FALSE)
plot_grid(p1, p2, p3, p4, ncol = 2)
dev.off()

## umap split on response, sample and treatment
pdf(paste0("Split_Response_", proj_name, ".pdf"), width = 15, height = 8)
DimPlot(seu, reduction = "umap", split.by = "Response")
dev.off()
pdf(paste0("Split_Sample_", proj_name, ".pdf"), width = 15, height = 8)
DimPlot(seu, reduction = "umap", split.by = "orig.ident")
dev.off()
pdf(paste0("Split_Treatment_", proj_name, ".pdf"), width = 15, height = 8)
DimPlot(seu, reduction = "umap", split.by = "Treatment")
dev.off()

########################## tables for QC plots ################################
## make sure you see the correct number of clusters. if not, SetIdent to the
## correct column in  metadata
# table(seu@active.ident)

## save select metadata
u <- seu@meta.data

## table of number of cells
write.table(u, paste0("cells_metadata_", proj_name, ".txt"), col.names = NA,
            sep = "\t")

## how many cells in each cluster per sample
cell_perC_perO <- table(seu@active.ident, seu@meta.data$orig.ident)
write.table(cell_perC_perO, paste0("cell_per_c_per_O_", proj_name, ".txt"),
            col.names = NA, sep = "\t")

## how many cells in each cluster per treatment
cell_perC_perT <- table(seu@active.ident, seu@meta.data$Treatment)
write.table(cell_perC_perT, paste0("cell_per_c_per_T_", proj_name, ".txt"),
            col.names = NA, sep = "\t")

## how many cells in each cluster per response
cell_perC_perR <- table(seu@active.ident, seu@meta.data$Response)
write.table(cell_perC_perR, paste0("cell_per_c_per_R_", proj_name, ".txt"),
            col.names = NA, sep = "\t")

## Plot number of genes, nUMIs and %mito on a tSNE - 3 plots per sample
## These allows you to see if part of your umap that has a high myeloid content
## (dying cells) or high number of genes (normally, no 2 cells in a droplets
## because you still filtered these bad cells!)
pdf(paste0("graph_nfeature_", proj_name, ".pdf"), width = 15, height = 8)
FeaturePlot(seu, features = "nFeature_RNA", cols = c("cadetblue2", "darkred"),
            pt.size = 0.1, reduction = "umap", label = FALSE)
dev.off()
pdf(paste0("graph_nCount_", proj_name, ".pdf"), width = 15, height = 8)
FeaturePlot(seu, features = "nCount_RNA", cols = c("cadetblue2", "darkred"),
            pt.size = 0.1, reduction = "umap", label = FALSE)
dev.off()
pdf(paste0("graph_percentmt_", proj_name, ".pdf"), width = 15, height = 8)
FeaturePlot(seu, features = "percent.mt", cols = c("cadetblue2", "darkred"),
            pt.size = 0.1, reduction = "umap", label = FALSE)
dev.off()

# saveRDS(seu, file = "sub_b_cells.rds")

seu@meta.data$Cluster <- as.character(seu@meta.data$Cluster)
seu@meta.data$Cluster <- as.factor(seu@meta.data$Cluster)

#seu$meta.data$Cluster <- gsub("TILs", "Proliferating T-cells", seu$meta.data$Cluster)

c("Mature Naïve B-cells Gapdh-",
  "Mature Naïve B-cells Hsp90ab1+",
  "Mature Naïve B-cells Fau+",
  "Proliferating B-cells",
  "Memory B-cells Class Switched") -> new.cluster.ids

## add column in metadata with the "Cluster" title and the cluster name
names(new.cluster.ids) <- levels(seu)
seu <- RenameIdents(seu, new.cluster.ids)
seu <- StashIdent(seu, save.name = "Cluster")
head(seu@meta.data) ## verify new column
# 
# pdf(paste0("Annotated_v2", proj_name, ".pdf"))
# DimPlot(seu, reduction = "umap", label = FALSE, pt.size = 0.1) +
#   labs(title = "Annotated UMAP") + theme(plot.title = element_text(hjust = 0.5))
# dev.off()

pdf(paste0("Annotated_labels_", proj_name, ".pdf"), width = 11, height = 8.5)
DimPlot(seu, reduction = "umap", label = TRUE, pt.size = 0.1, repel = TRUE) +
  labs(title = "Annotated UMAP") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

################################## pericytes ###################################

seu@meta.data$new_labels <- paste0(seu@meta.data$Treatment, "_",
                                   seu@meta.data$Response)
seu <- SetIdent(object = seu, value = 'new_labels')

############################### neutrophils ####################################
# seu <- readRDS(file = "neutro_no104.rds")
# seu <- SetIdent(object = seu, value = 'orig.ident')
seu0 <- readRDS(file = paste0("neutro1.rds"))
proj_name <- "neutro_no104"
seu <- subset(seu0, subset = orig.ident != c("CCB104"))

seu %>%  
  FindVariableFeatures(selection.method = "vst", nfeatures = var_feat) %>%
  ScaleData() %>% 
  RunPCA(npcs = sig_dims, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:sig_dims, k.param = 30,
                compute.SNN = TRUE, nn.method = "rann", nn.eps = 0) %>% 
  FindClusters(resolution = c(0.1, 0.2, 0.3, 0.4, 0.5)) %>% 
  RunUMAP(reduction = "pca", dims = 1:sig_dims) -> FC

## multiple resolution plots
pdf(file = paste0("clusters_diff_res_", proj_name, ".pdf"), width = 15, height = 8)
Seurat_02 <- SetIdent(object = FC, value = 'RNA_snn_res.0.2')
Seurat02 <- RunUMAP(Seurat_02, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat02, reduction = "umap", pt.size = 0.25, label = TRUE) +
  labs(title = "UMAP with res 0.2") + theme(plot.title = element_text(hjust = 0.5))
Seurat_04 <- SetIdent(object = FC, value = 'RNA_snn_res.0.4')
Seurat04 <- RunUMAP(Seurat_04, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat04, reduction = "umap", pt.size = 0.25, label = TRUE) +
  labs(title = "UMAP with res 0.4") + theme(plot.title = element_text(hjust = 0.5))
dev.off() 

## write the column name of the chosen resolution with the SetIdent function
## Seurat will set the object with this column
seu <- SetIdent(object = FC, value = RESOLUTION)
seu <- RunUMAP(seu, reduction = "pca", dims = 1:sig_dims)
pdf(file = paste0("res_", RESOLUTION, "_clusters_", proj_name, ".pdf"),
    width = 15, height = 8)
UMAPPlot(seu, reduction = "umap", pt.size = 0.1, label = TRUE,
         label.size = 5) + NoLegend() +
  labs(title = paste0("UMAP with ", RESOLUTION)) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

## produce dimplots
pdf(paste0("umap_", proj_name, ".pdf"), width = 15, height = 8)
p1 <- DimPlot(seu, reduction = "umap", group.by = "orig.ident",
              pt.size = 0.1, label = FALSE, label.size = 2.0)
p2 <- DimPlot(seu, reduction = "umap", group.by = "Response",
              pt.size = 0.1, label = FALSE, label.size = 3.5,  repel = TRUE)
p3 <- DimPlot(seu, reduction = "umap", group.by = "Treatment",
              pt.size = 0.1, label = FALSE, label.size = 3.5,  repel = TRUE)
p4 <- DimPlot(seu, reduction = "umap",
              pt.size = 0.1, label = TRUE, label.size = 3.5,  repel = FALSE)
plot_grid(p1, p2, p3, p4, ncol = 2)
dev.off()

