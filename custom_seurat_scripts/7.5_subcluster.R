############################### sc-pipeline ####################################
## this script performs downstream analysis on Seurat objects taken from the
## previous script `2_seurat.R` and filters the objects further, normalizes the
## data and finds variable features. running PCA, UMAP and generating various
## diagnostic plots. Lastly, the clusters resulting from the UMAP are annotated 
## sample and response type, filters the cells, produces QC plots

## load libraries
library("R.utils")
library("dplyr")
library("Seurat")
library("ggplot2")
library("future")
library("cowplot")

## parallelize workflow
## only run this code chunk if you are on the VSC
plan("multiprocess", workers = 45) # uses 50 CPU
options(future.globals.maxSize = 15000 * 1024^2) ## 5GB per worker

## evergreen variables
proj_name <- "rerun" ## name of project
proj_version <- "whole_data"       ## name of version of project
# server_path <- "/media/seq-srv-05/vrc/Project/Project_Theo/pancreatic_mice/unfiltered/macro_sub2/"
server_path <- "/media/seq-srv-05/vrc/Project/Project_Theo/pancreatic_mice/rerun/"
mito <- "^mt-" ## "^mt-" if mouse, "^MT-" if human
up_feat <- 6000 ## upper nFeature threshold
low_feat <- 1000 ## lower nFeature threshold
var_feat <- 2000 ## variableFeature threshold
low_count <- 1000 ## lower nCount threshold
p_mt <- 20 ## percent mito threshold
NPC <- 35 ## number of PCs of first PCA
k_para <- 50 ## number of k paramaters in FindNeighbors default is 30
sig_dims <- 10 ## number of PCs of final PCA
resolution <- 0.8 ## resolution
RESOLUTION <- 'RNA_snn_res.0.8' ## resolution string

########################### subclustering contamination ##0000000000000000000000
# seu0 <- readRDS("final_annotations_recluster.rds")
proj_name <- "contamination1"
seu0 <- readRDS("hf_rerun_complete_metadata.rds")
t1 <- subset(seu0, subset = subtype_annotations_v3 == "Unknown Lymphoid")
m1 <- subset(seu0, subset = subtype_annotations_v3 == "Unknown Myeloid")
e1 <- subset(seu0, subset = subtype_annotations_v3 == "Unknown Endothelial Cell")
b1 <- subset(seu0, subset = subtype_annotations_v3 == "Unknown B-cell")
m1 %>% merge(y = c(t1, e1, b1)) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = var_feat) %>%
  ScaleData() %>% 
  RunPCA(npcs = sig_dims, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:sig_dims, k.param = k_para,
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
#seu <- SetIdent(object = FC, value = 'RNA_snn_res.1.4')
seu <- RunUMAP(seu, reduction = "pca", dims = 1:sig_dims)
pdf(file = paste0("res_", RESOLUTION, "_clusters_", proj_name, ".pdf"),
    width = 15, height = 8)
UMAPPlot(seu, reduction = "umap", pt.size = 0.1, label = TRUE,
         label.size = 5) + NoLegend() +
  labs(title = paste0("UMAP with ", RESOLUTION)) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

seu <- SetIdent(object = seu, value = 'subtype_annotations_v3')

#### umap old labels
pdf(file = paste0("umap_doublets_", proj_name, ".pdf"), width = 8, height = 8)
DimPlot(seu, reduction = "umap", group.by = "pANNPredictions", pt.size = 0.1,
        label = FALSE, label.size = 2.0)
dev.off()

## produce dimplots
pdf(paste0("umap_", proj_name, ".pdf"), width = 15, height = 8)
p1 <- DimPlot(seu, reduction = "umap", group.by = "orig.ident",
              pt.size = 0.1, label = FALSE, label.size = 2.0)
p2 <- DimPlot(seu, reduction = "umap", group.by = "Response",
              pt.size = 0.1, label = FALSE, label.size = 3.5,  repel = TRUE)
p3 <- DimPlot(seu, reduction = "umap", group.by = "Treatment",
              pt.size = 0.1, label = FALSE, label.size = 3.5,  repel = TRUE)
p4 <- DimPlot(seu, reduction = "umap", group.by = "subtype_annotations_v3",
              pt.size = 0.1, label = TRUE, label.size = 3.5,  repel = TRUE)
plot_grid(p1, p2, p3, p4, ncol = 2)
dev.off()

pdf(paste0("markers_cancer_IT.pdf"), width = 15, height = 8)
FeaturePlot(seu, features = c(
  "Ins1", "Ins2", "Iapp", "Ppy", "Insm1", "Ppp1r1a" ## IT
  ), cols = c("lightblue", "red"), order = TRUE, pt.size = 0.1)
dev.off()

pdf(paste0("markers_cancer_acinar.pdf"), width = 15, height = 8)
FeaturePlot(seu, features = c(
  "Cela2a", "Cpa1", "Lars2", "Prss2", "Try5", "Ctrb1" ## acinar
  ), cols = c("lightblue", "red"), order = TRUE, pt.size = 0.1)
dev.off()

saveRDS(seu, file = "contam_cells.rds")

########################### subclustering CANCER ##11111111111111111111111111111
# seu0 <- readRDS("final_annotations_recluster.rds")
proj_name <- "cancer1"
# c1 <- subset(seu0, subset = final_primary_annotations == "Cancer")
c1 <- readRDS("hf_rerun_cancer.rds")
# m1 <- readRDS(file = "../macro1/macro_contam_all_groups1.rds")
# t1 <- readRDS(file = "../tcell1/tcell_contam_actual.rds")
# t2 <- readRDS(file = "../tcell1/tcell_also_contam_actual1.rds")
# e1 <- readRDS(file = "../endo1/contamination/endo_cancer_contam_endo.rds")
# c1 %>% merge(y = c(m1, t1, t2, e1)) -> merged_can

## run pipeline
# merged_can %>% 
c1 %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = var_feat) %>%
  ScaleData() %>% 
  RunPCA(npcs = sig_dims, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:sig_dims, k.param = k_para,
                compute.SNN = TRUE, nn.method = "rann", nn.eps = 0) %>% 
  FindClusters(resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                              1.1, 1.2, 1.3, 1.4, 1.6)) %>% 
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
#seu <- SetIdent(object = FC, value = 'RNA_snn_res.1.4')
seu <- RunUMAP(seu, reduction = "pca", dims = 1:sig_dims)
pdf(file = paste0("res_", RESOLUTION, "_clusters_", proj_name, ".pdf"),
    width = 15, height = 8)
UMAPPlot(seu, reduction = "umap", pt.size = 0.1, label = TRUE,
         label.size = 5) + NoLegend() +
  labs(title = paste0("UMAP with ", RESOLUTION)) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

#### umap old labels
pdf(file = paste0("umap_doublets_", proj_name, ".pdf"), width = 8, height = 8)
DimPlot(seu, reduction = "umap", group.by = "pANNPredictions", pt.size = 0.1,
        label = FALSE, label.size = 2.0)
dev.off()

pdf(file = paste0("umap_old_labels_", proj_name, ".pdf"), width = 12, height = 8)
DimPlot(seu, reduction = "umap", group.by = "final_subtype_annotations",
        pt.size = 0.1, label = TRUE, label.size = 2.0)
dev.off()

## produce dimplots
pdf(paste0("umap_", proj_name, ".pdf"), width = 15, height = 8)
p1 <- DimPlot(seu, reduction = "umap", group.by = "orig.ident",
              pt.size = 0.1, label = FALSE, label.size = 2.0)
p2 <- DimPlot(seu, reduction = "umap", group.by = "Response",
              pt.size = 0.1, label = FALSE, label.size = 3.5,  repel = TRUE)
p3 <- DimPlot(seu, reduction = "umap", group.by = "Treatment",
              pt.size = 0.1, label = FALSE, label.size = 3.5,  repel = TRUE)
p4 <- DimPlot(seu, reduction = "umap", #group.by = "initial_subtype_annotations",
              pt.size = 0.1, label = TRUE, label.size = 3.5,  repel = TRUE)
plot_grid(p1, p2, p3, p4, ncol = 2)
dev.off()

# saveRDS(seu, file = "cancer_w_contam.rds")
saveRDS(seu, file = "cancer_no_contam.rds")

seu <- SetIdent(object = seu, value = "RNA_snn_res.0.2")
#seu <- SetIdent(object = seu, value = "current_subtype_annotations")
# cl0 <- subset(seu, subset = RNA_snn_res.0.8 == 0)
# cl1 <- subset(seu, subset = RNA_snn_res.0.8 == 1)
# length(rownames(cl0@assays$RNA@counts))
# length(rownames(cl1@assays$RNA@counts))

c("IT-like", "IT-like Ins-lo B20 Responding", "IT-like B20 aPDL1 CTLA4 Responding",
"IT-like Ins-lo B20 + PiK3i Relapsing", "Proliferating IT-like", "IT-like",
"IT-like Ins-lo B20", "MLP-like Ins-lo") -> new.cluster.ids

## add column in metadata with the "Cluster" title and the cluster name
names(new.cluster.ids) <- levels(seu)
seu <- RenameIdents(seu, new.cluster.ids)
seu <- StashIdent(seu, save.name = "current_subtype_annotations")
head(seu@meta.data) ## verify new column

pdf(paste0("Annotated_", proj_name, ".pdf"), width = 11, height = 8.5)
DimPlot(seu, reduction = "umap", label = FALSE, pt.size = 0.1) +
  labs(title = "Annotated UMAP") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf(paste0("Annotated_labels_", proj_name, ".pdf")) 
DimPlot(seu, reduction = "umap", label = TRUE, pt.size = 0.1, repel = TRUE) +
  labs(title = "Annotated UMAP") + theme(plot.title = element_text(hjust = 0.5))
dev.off() # the same with labels

saveRDS(seu, file = "cancer_actual1.rds")

## gradients
pdf(paste0("cancer_invasion.pdf"), width = 15, height = 8)
FeaturePlot(seu, features = c(
  "Pecam1", "Mki67", "Vim", "Alt", "Lamc2", "Sna", "Twist1", "Licam" ## invasion
  ), cols = c("lightblue", "red"), order = TRUE, pt.size = 0.1)
dev.off()

pdf(paste0("cancer_IT1.pdf"), width = 15, height = 8)
FeaturePlot(seu, features = c(
  "Ins1", "Ins2", "Iapp", "Glut2", "Gck", "Isl1", "Insm1", "Ppp1r1a", "Slc2a2" ## IT
  ), cols = c("lightblue", "red"), order = TRUE, pt.size = 0.1)
dev.off()

pdf(paste0("cancer_IT2.pdf"), width = 15, height = 8)
FeaturePlot(seu, features = c(
  "Slc8a2", "Cfc1", "Sez6l", "Nkx6-1", "Nkx6-2", "Pdx1", "Wnt4", "Pax6" ## IT
  ), cols = c("lightblue", "red"), order = TRUE, pt.size = 0.1)
dev.off()

pdf(paste0("cancer_MLP1.pdf"), width = 15, height = 8)
FeaturePlot(seu, features = c(
  "Enpp2", "Chga", "Scgn", "Slc16a1", "Slc2a3", "Hk1", "Gmp6a", "Mct1" ## MLP
  ), cols = c("lightblue", "red"), order = TRUE, pt.size = 0.1)
dev.off()

pdf(paste0("cancer_MLP2.pdf"), width = 15, height = 8)
FeaturePlot(seu, features = c(
  "Hnf1b","Gata6", "Lin28b", "Sema3e", "Epha3", "Pou3f4" ## MLP?
  ), cols = c("lightblue", "red"), order = TRUE, pt.size = 0.1)
dev.off()

pdf(paste0("cancer_acinar.pdf"), width = 15, height = 8)
FeaturePlot(seu, features = c(
  "Cela2a", "Cpa1", "Lars2", "Prss2", "Try5", "Ctrb1" ## acinar
  ), cols = c("lightblue", "red"), order = TRUE, pt.size = 0.1)
dev.off()

pdf(paste0("cancer_ductalstem.pdf"), width = 15, height = 8)
FeaturePlot(seu, features = c(
  "Lgr5", "Sox9", "Foxj1", # "Ehf", ## stem
  "Clu", "Krt18", "Krt8", "Krt19", "Krt7", "Epcam" ## ductal
  ), cols = c("lightblue", "red"), order = TRUE, pt.size = 0.1)
dev.off()

########################### subclustering MACRO/MONO ##2222222222222222222222222
proj_name <- "myeloid_rerun"
# seu0 <- readRDS("../final_annotations_recluster.rds")
# ma1 <- subset(seu0, subset = final_primary_annotations == "Macrophages")
# mo1 <- subset(seu0, subset = final_primary_annotations == "Monocytes")
# ma1 %>% merge(y = mo1, project = proj_name) -> merged_mm
# merged_mm <- readRDS("myeloid_recluster_v3_correct.rds")
# saveRDS(seu, file = "myeloid_recluster_v3_correct_no_cont.rds")
# seu <- readRDS(seu, file = "myeloid_recluster_v3_correct_no_cont.rds")
# mm_contam <- subset(seu, subset = RNA_snn_res.0.8 == 5)
# saveRDS(mm_contam, file = "myeloid_recluster_v3_correct_CONT.rds")
# merged_mm <- subset(seu, subset = RNA_snn_res.0.8 != 5) ## remove contamination
# saveRDS(merged_mm, file = "myeloid_recluster_v3_correct_NO_CONT.rds")
# merged_mm <- SetIdent(object = merged_mm, value = 'orig.ident')
#m0 <- readRDS("hf_rerun_myeloids.rds") ## 7753 cells
m0 <- readRDS(file = "macro_temp1.rds") ## 7753 cells
m3 <- subset(m0, subset = cluster == 3) # 774
m5 <- subset(m0, subset = cluster == 5) # 441
m8 <- subset(m0, subset = cluster == 8) # 370
m13 <- subset(m0, subset = cluster == 13) # 238
M1 <- subset(m0, subset = cluster != 3)
M2 <- subset(M1, subset = cluster != 5)
M3 <- subset(M2, subset = cluster != 8)
M4 <- subset(M3, subset = cluster != 13) ## 5930
cont_m <- subset(M4, subset = final_subtype_annotations == "Macrophage cancer contamination") 
M5 <- subset(M4, subset = final_subtype_annotations != "Macrophage cancer contamination") 
seu_ins <- subset(M5, subset = Ins1 > 4 & Ins2 > 4) # 356
M6 <- subset(M5, subset = Ins1 < 4 | Ins2 < 4)
m3 %>% merge(y = c(m5, m8, m13, cont_m, seu_ins), project = proj_name) -> merged_mm_contam # 1823
saveRDS(merged_mm_contam, file = "macro_contam_all_groups1.rds")
b0 <- subset(M5, subset = cluster == 14)
saveRDS(b0, file = "bcells.rds")
M6 <- subset(M5, subset = cluster != 14)
saveRDS(M6, file = "macro_no_contam_actual1.rds") # 5531

## run pipeline
M6 %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = var_feat) %>%
  ScaleData() %>% 
  RunPCA(npcs = sig_dims, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:sig_dims, k.param = k_para,
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

#### umap old labels
# pdf(file = paste0("umap_", RESOLUTION, "_oldlabels_", proj_name, ".pdf"),
#     width = 15, height = 8)
# DimPlot(seu, reduction = "umap", group.by = "initial_subtype_annotations",
#               pt.size = 0.1, label = TRUE, label.size = 3.5,  repel = TRUE)
# dev.off()

pdf(file = paste0("umap_doublets_", proj_name, ".pdf"), width = 8, height = 8)
DimPlot(seu, reduction = "umap", group.by = "pANNPredictions", pt.size = 0.1,
        label = FALSE, label.size = 2.0)
dev.off()

pdf(file = paste0("umap_old_labels_", proj_name, ".pdf"), width = 12, height = 8)
DimPlot(seu, reduction = "umap", group.by = "final_subtype_annotations",
        pt.size = 0.1, label = TRUE, label.size = 2.0)
dev.off()

## produce dimplots
pdf(paste0("umap_", proj_name, ".pdf"), width = 15, height = 8)
p1 <- DimPlot(seu, reduction = "umap", group.by = "orig.ident",
              pt.size = 0.1, label = FALSE, label.size = 2.0)
p2 <- DimPlot(seu, reduction = "umap", group.by = "Response",
              pt.size = 0.1, label = FALSE, label.size = 3.5,  repel = TRUE)
p3 <- DimPlot(seu, reduction = "umap", group.by = "Treatment",
              pt.size = 0.1, label = FALSE, label.size = 3.5,  repel = TRUE)
p4 <- DimPlot(seu, reduction = "umap", #group.by = "initial_subtype_annotations",
              pt.size = 0.1, label = TRUE, label.size = 3.5,  repel = TRUE)
plot_grid(p1, p2, p3, p4, ncol = 2)
dev.off()

# saveRDS(seu, file = "macro_temp1.rds")
saveRDS(seu, file = "macro_temp2.rds")

########################## tables for QC plots ################################
## table of number of cells 
# write.table(seu@meta.data, paste0("cells_metadata_", proj_name, ".txt"),
#             col.names = NA, sep = "\t")
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

################################ SingleR #######################################
# library("SingleR")
# library("celldex")
# mouse <- celldex::MouseRNAseqData()
# seuratObj <- seu
# singlerObj <- as.SingleCellExperiment(seuratObj)
# common <- intersect(rownames(singlerObj), rownames(mouse))
# mouse <- mouse[common,]
# singlerObj <- singlerObj[common,]
# singler <- SingleR(test = singlerObj, ref = mouse, labels = mouse$label.main,
#                    method = "cluster", ## change the resolution you chose
#                    clusters = seuratObj@meta.data[,"RNA_snn_res.1.4"],
#                    genes = "de") 

## save singlr results 
# saveRDS(singler, file = paste0("singler_object_", proj_name, ".rds"))
# write.table(singler, file = paste0( "singler_", proj_name, ".txt"),
#             col.names = NA, sep = "\t")

############# Assigning cell type identity to clusters (annotation) ############
## the cluster names are out of order
# cluster_ids <- singler@rownames
# singler_labels <- singler@listData$labels
# as.data.frame(cbind(cluster_ids, singler_labels)) %>%
#   dplyr::mutate(real_cluster_ids = as.integer(as.character(cluster_ids))) %>% 
#   dplyr::arrange(real_cluster_ids) -> singler_annotations
## unfiltered2 @RES = 1.4
# c("M2 Macrophages", "M2 Macrophages", "M2 Macrophages", "M1 Macrophages",
#   "M2 Macrophages", "M1 Macrophages", "M1 Macrophages", "Prol. M2 Macrophages",
#   "M2 Macrophages", "Monocytes", "Mast cells", "Clec10a+ M2 Macrophages",
#   "Monocytes", "M2 Macrophages", "pDCs", "cDCs", "Neutrophils") -> new.cluster.ids

c("M2 Macrophages", "M2 Macrophages", "M1 Macrophages", "Unknown Myeloids",
  "M1 Macrophages", "M1 Macrophages", "Prol. M2 Macrophages", "Monocytes/Neutrophils",
  "Mast cells/pDCs", "F13a1+ M2 Macrophages", "M2 Macrophages", "cDCs") -> new.cluster.ids

## add column in metadata with the "Cluster" title and the cluster name
names(new.cluster.ids) <- levels(seu)
seu <- RenameIdents(seu, new.cluster.ids)
seu <- StashIdent(seu, save.name = "current_subtype_annotations")
head(seu@meta.data) ## verify new column

pdf(paste0("Annotated_", proj_name, ".pdf"), width = 11, height = 8.5)
DimPlot(seu, reduction = "umap", label = FALSE, pt.size = 0.1) +
  labs(title = "Annotated UMAP") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf(paste0("Annotated_labels_", proj_name, ".pdf")) 
DimPlot(seu, reduction = "umap", label = TRUE, pt.size = 0.1, repel = TRUE) +
  labs(title = "Annotated UMAP") + theme(plot.title = element_text(hjust = 0.5))
dev.off() # the same with labels

saveRDS(seu, file = "macro_temp3.rds")
# seu <- readRDS(file = "macro_temp3.rds")
MAST <- subset(seu, subset = current_subtype_annotations == "Mast cells/pDCs")
saveRDS(MAST, file = "macro_mast_pdc.rds")
no_mast <- subset(seu, subset = current_subtype_annotations != "Mast cells/pDCs")
mono <- subset(no_mast, subset = current_subtype_annotations == "Monocytes/Neutrophils")
saveRDS(mono, file = "macro_mono_neutro.rds")
no_subs <- subset(no_mast, subset = current_subtype_annotations != "Monocytes/Neutrophils")
saveRDS(no_subs, file = "macro_temp4_missing2.rds")

########################### subclustering mast/pdc cells
proj_name <- "mast_pdc_rerun"
MAST <- readRDS("macro_mast_pdc.rds")
MAST  %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = var_feat) %>%
  ScaleData() %>% 
  RunPCA(npcs = sig_dims, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:sig_dims, k.param = k_para,
                compute.SNN = TRUE, nn.method = "rann", nn.eps = 0) %>% 
  FindClusters(resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                              1.1, 1.2, 1.3, 1.4)) %>% 
  RunUMAP(reduction = "pca", dims = 1:sig_dims) -> FC

## multiple resolution plots
pdf(file = paste0("clusters_diff_res_", proj_name, ".pdf"), width = 15, height = 8)
Seurat_02 <- SetIdent(object = FC, value = 'RNA_snn_res.0.2')
Seurat02 <- RunUMAP(Seurat_02, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat02, reduction = "umap", pt.size = 1, label = TRUE) +
  labs(title = "UMAP with res 0.2") + theme(plot.title = element_text(hjust = 0.5))
Seurat_04 <- SetIdent(object = FC, value = 'RNA_snn_res.0.4')
Seurat04 <- RunUMAP(Seurat_04, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat04, reduction = "umap", pt.size = 1, label = TRUE) +
  labs(title = "UMAP with res 0.4") + theme(plot.title = element_text(hjust = 0.5))
Seurat_06 <- SetIdent(object = FC, value = 'RNA_snn_res.0.6')
Seurat06 <- RunUMAP(Seurat_06, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat06, reduction = "umap", pt.size = 1, label = TRUE) +
  labs(title = "UMAP with res 0.6") + theme(plot.title = element_text(hjust = 0.5))
Seurat_08 <- SetIdent(object = FC, value = 'RNA_snn_res.0.8')
Seurat08 <- RunUMAP(Seurat_08, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat08, reduction = "umap", pt.size = 1, label = TRUE) +
  labs(title = "UMAP with res 0.8") + theme(plot.title = element_text(hjust = 0.5))
Seurat_1 <- SetIdent(object = FC, value = 'RNA_snn_res.1.0')
Seurat1 <- RunUMAP(Seurat_1, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat1, reduction = "umap", pt.size = 1, label = TRUE) + 
  labs(title = "UMAP with res 1.0") + theme(plot.title = element_text(hjust = 0.5))
Seurat_12 <- SetIdent(object = FC, value = 'RNA_snn_res.1.2')
Seurat12 <- RunUMAP(Seurat_12, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat12, reduction = "umap", pt.size = 1, label = TRUE) + 
  labs(title = "UMAP with res 1.2") + theme(plot.title = element_text(hjust = 0.5))
Seurat_14 <- SetIdent(object = FC, value = 'RNA_snn_res.1.4')
Seurat14 <- RunUMAP(Seurat_14, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat14, reduction = "umap", pt.size = 1, label = TRUE) +
  labs(title = "UMAP with res 1.4") + theme(plot.title = element_text(hjust = 0.5))
dev.off() 

## write the column name of the chosen resolution with the SetIdent function

## write the column name of the chosen resolution with the SetIdent function
## Seurat will set the object with this column
# seu <- SetIdent(object = FC, value = RESOLUTION)
seu <- SetIdent(object = FC, value = 'RNA_snn_res.1.4')
seu <- RunUMAP(seu, reduction = "pca", dims = 1:sig_dims)
pdf(file = paste0("res_", RESOLUTION, "_clusters_", proj_name, ".pdf"),
    width = 15, height = 8)
UMAPPlot(seu, reduction = "umap", pt.size = 2, label = TRUE,
         label.size = 5) + NoLegend() +
  labs(title = paste0("UMAP with ", RESOLUTION)) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf(file = paste0("umap_old_labels_", proj_name, ".pdf"), width = 12, height = 8)
DimPlot(seu, reduction = "umap", group.by = "final_subtype_annotations",
        pt.size = 1, label = TRUE, label.size = 2.0)
dev.off()

## produce dimplots
pdf(paste0("umap_", proj_name, ".pdf"), width = 15, height = 8)
p1 <- DimPlot(seu, reduction = "umap", group.by = "orig.ident",
              pt.size = 1, label = FALSE, label.size = 2.0)
p2 <- DimPlot(seu, reduction = "umap", group.by = "Response",
              pt.size = 1, label = FALSE, label.size = 3.5,  repel = TRUE)
p3 <- DimPlot(seu, reduction = "umap", group.by = "Treatment",
              pt.size = 1, label = FALSE, label.size = 3.5,  repel = TRUE)
p4 <- DimPlot(seu, reduction = "umap", #group.by = "initial_subtype_annotations",
              pt.size = 1, label = TRUE, label.size = 3.5,  repel = TRUE)
plot_grid(p1, p2, p3, p4, ncol = 2)
dev.off()

saveRDS(seu, file = "macro_mast_pdc_temp1.rds")

pdf(paste0("mast_markers.pdf"), width = 15, height = 8)
FeaturePlot(seu, features = c(
  "Kit", "Cpa3", "Cma1", "Il1rl1", "Osbpl8", "Adora3", "Csf2rb","Cyp11a1",
  "Hdc", "Il4", "Tcf4", "Bcl11a", "Irf8", "Spib", "Runx2", "Gpnmb"),
            cols = c("lightblue", "red"), order = TRUE, pt.size = 1)
dev.off()

pdf(paste0("tpdc_markers.pdf"), width = 15, height = 8)
FeaturePlot(seu, features = c(
  "Axl", "Cx3cr", "Cd2", "Cd5", "Cd81", "Siglec1", "Pira2", "Cxcr3", "Irf7"),
            cols = c("lightblue", "red"), order = TRUE, pt.size = 1)
dev.off()

########################### subclustering mast/pdc cells
proj_name <- "mono_neutro_rerun"
MONO <- readRDS("macro_mono_neutro.rds")
MONO  %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = var_feat) %>%
  ScaleData() %>% 
  RunPCA(npcs = sig_dims, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:sig_dims, k.param = k_para,
                compute.SNN = TRUE, nn.method = "rann", nn.eps = 0) %>% 
  FindClusters(resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                              1.1, 1.2, 1.3, 1.4)) %>% 
  RunUMAP(reduction = "pca", dims = 1:sig_dims) -> FC

## multiple resolution plots
pdf(file = paste0("clusters_diff_res_", proj_name, ".pdf"), width = 15, height = 8)
Seurat_02 <- SetIdent(object = FC, value = 'RNA_snn_res.0.2')
Seurat02 <- RunUMAP(Seurat_02, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat02, reduction = "umap", pt.size = 1, label = TRUE) +
  labs(title = "UMAP with res 0.2") + theme(plot.title = element_text(hjust = 0.5))
Seurat_04 <- SetIdent(object = FC, value = 'RNA_snn_res.0.4')
Seurat04 <- RunUMAP(Seurat_04, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat04, reduction = "umap", pt.size = 1, label = TRUE) +
  labs(title = "UMAP with res 0.4") + theme(plot.title = element_text(hjust = 0.5))
Seurat_06 <- SetIdent(object = FC, value = 'RNA_snn_res.0.6')
Seurat06 <- RunUMAP(Seurat_06, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat06, reduction = "umap", pt.size = 1, label = TRUE) +
  labs(title = "UMAP with res 0.6") + theme(plot.title = element_text(hjust = 0.5))
Seurat_08 <- SetIdent(object = FC, value = 'RNA_snn_res.0.8')
Seurat08 <- RunUMAP(Seurat_08, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat08, reduction = "umap", pt.size = 1, label = TRUE) +
  labs(title = "UMAP with res 0.8") + theme(plot.title = element_text(hjust = 0.5))
Seurat_1 <- SetIdent(object = FC, value = 'RNA_snn_res.1.0')
Seurat1 <- RunUMAP(Seurat_1, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat1, reduction = "umap", pt.size = 1, label = TRUE) + 
  labs(title = "UMAP with res 1.0") + theme(plot.title = element_text(hjust = 0.5))
Seurat_12 <- SetIdent(object = FC, value = 'RNA_snn_res.1.2')
Seurat12 <- RunUMAP(Seurat_12, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat12, reduction = "umap", pt.size = 1, label = TRUE) + 
  labs(title = "UMAP with res 1.2") + theme(plot.title = element_text(hjust = 0.5))
Seurat_14 <- SetIdent(object = FC, value = 'RNA_snn_res.1.4')
Seurat14 <- RunUMAP(Seurat_14, reduction = "pca", dims = 1:sig_dims)
UMAPPlot(Seurat14, reduction = "umap", pt.size = 1, label = TRUE) +
  labs(title = "UMAP with res 1.4") + theme(plot.title = element_text(hjust = 0.5))
dev.off() 

## write the column name of the chosen resolution with the SetIdent function

## write the column name of the chosen resolution with the SetIdent function
## Seurat will set the object with this column
# seu <- SetIdent(object = FC, value = RESOLUTION)
seu <- SetIdent(object = FC, value = 'RNA_snn_res.0.8')
seu <- RunUMAP(seu, reduction = "pca", dims = 1:sig_dims)
pdf(file = paste0("res_", RESOLUTION, "_clusters_", proj_name, ".pdf"),
    width = 15, height = 8)
UMAPPlot(seu, reduction = "umap", pt.size = 2, label = TRUE,
         label.size = 5) + NoLegend() +
  labs(title = paste0("UMAP with ", RESOLUTION)) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf(file = paste0("umap_old_labels_", proj_name, ".pdf"), width = 12, height = 8)
DimPlot(seu, reduction = "umap", group.by = "final_subtype_annotations",
        pt.size = 1, label = TRUE, label.size = 2.0)
dev.off()

## produce dimplots
pdf(paste0("umap_", proj_name, ".pdf"), width = 15, height = 8)
p1 <- DimPlot(seu, reduction = "umap", group.by = "orig.ident",
              pt.size = 1, label = FALSE, label.size = 2.0)
p2 <- DimPlot(seu, reduction = "umap", group.by = "Response",
              pt.size = 1, label = FALSE, label.size = 3.5,  repel = TRUE)
p3 <- DimPlot(seu, reduction = "umap", group.by = "Treatment",
              pt.size = 1, label = FALSE, label.size = 3.5,  repel = TRUE)
p4 <- DimPlot(seu, reduction = "umap", #group.by = "initial_subtype_annotations",
              pt.size = 1, label = TRUE, label.size = 3.5,  repel = TRUE)
plot_grid(p1, p2, p3, p4, ncol = 2)
dev.off()

pdf(paste0("mono_markers.pdf"), width = 15, height = 8)
FeaturePlot(seu, features = c(
  "Cd40", "Hck", "Cxcr4", "Cd86", "Cd14", "Fcgr2b", "Mx1", "Il1rn", "Ifit1",
  "Ifit3", "Cxcl10"),
            cols = c("lightblue", "red"), order = TRUE, pt.size = 1)
dev.off()

pdf(paste0("neutro_markers.pdf"), width = 15, height = 8)
FeaturePlot(seu, features = c(
  "Gsr", "S100a9", "Ly6g", "S100a8", "G0s2", "Ncf1", "Cd177", "Lrg1", "Camp"),
            cols = c("lightblue", "red"), order = TRUE, pt.size = 1)
dev.off()

c("Monocytes", "Monocytes", "Monocytes", "Neutrophils") -> new.cluster.ids

## add column in metadata with the "Cluster" title and the cluster name
names(new.cluster.ids) <- levels(seu)
seu <- RenameIdents(seu, new.cluster.ids)
seu <- StashIdent(seu, save.name = "current_subtype_annotations")
head(seu@meta.data) ## verify new column
saveRDS(seu, file = "mono_neutro_actual.rds")

## need to clean up metadata
seu <- readRDS(file = "macro_temp3.rds")
mono1 <- readRDS(file = "mono_neutro_actual.rds")
no_subs <- readRDS(file = "macro_temp4_missing2.rds")
mast1 <- readRDS(file = "macro_mast_pdc.rds")

seu@meta.data$cells <- rownames(seu@meta.data)
seu@meta.data %>% 
  dplyr::select(-current_subtype_annotations) -> meta_seu

mono1@meta.data$cells <- rownames(mono1@meta.data)
mono1@meta.data %>% 
  dplyr::select(cells, current_subtype_annotations) -> meta_mono1

mast1@meta.data$cells <- rownames(mast1@meta.data)
mast1@meta.data %>% 
  dplyr::select(cells, current_subtype_annotations) -> meta_mast1

no_subs@meta.data$cells <- rownames(no_subs@meta.data)
no_subs@meta.data %>% 
  dplyr::select(cells, current_subtype_annotations) -> meta_no_subs1

dplyr::bind_rows(meta_no_subs1, meta_mono1, meta_mast1) -> combined_meta
meta_seu %>%
  dplyr::left_join(combined_meta, by = "cells") %>% 
  as.data.frame() -> meta_clean
rownames(meta_clean) <- meta_clean$cells
seu@meta.data <- meta_clean
seu <- SetIdent(object = seu, value = 'current_subtype_annotations')

seu@meta.data$current_subtype_annotations <- gsub("Unknown Myeloids", "Unknown Macrophages",
                                                  seu@meta.data$current_subtype_annotations)
seu <- SetIdent(object = seu, value = 'current_subtype_annotations')

pdf(paste0("Annotated_", proj_name, ".pdf"), width = 11, height = 8)
DimPlot(seu, reduction = "umap", label = FALSE, pt.size = 0.1) +
  labs(title = "Annotated UMAP") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf(paste0("Annotated_labels_", proj_name, ".pdf"), width = 11, height = 8) 
DimPlot(seu, reduction = "umap", label = TRUE, pt.size = 0.1, repel = TRUE) +
  labs(title = "Annotated UMAP") + theme(plot.title = element_text(hjust = 0.5))
dev.off() # the same with labels

# saveRDS(seu, file = "macro_actual1.rds")
saveRDS(seu, file = "macro_actual2.rds")

########################### subclustering T/NK cells ##3333333333333333333333333
proj_name <- "tcell_subcluster"
# seu0 <- readRDS("rds2/final_annotations_recluster_v2_correct.rds")
# t1 <- subset(seu0, subset = final_primary_annotations == "T-cells")
#t1 <- readRDS("hf_rerun_tcells.rds") # 2965 cells round 1
t1 <- readRDS(file = "tcell_temp1.rds")
t5 <- subset(t1, subset = cluster == 5) # 210 cells
saveRDS(t5, file = "tcell_c5_cancer.rds")
t2 <- subset(t1, subset = cluster != 5)
cont_t2 <- subset(t2, subset = final_subtype_annotations == "T-cell cancer contamination") ## 706 cells
saveRDS(cont_t2, file = "tcell_contam_actual.rds")
t3 <- subset(t2, subset = final_subtype_annotations != "T-cell cancer contamination") ## 2049
saveRDS(t3, file = "tcell_no_contam_actual1.rds")
seu_ins <- subset(t3, subset = Ins1 > 4 & Ins2 > 4) ## 109 cells
saveRDS(seu_ins, file = "tcell_also_contam_actual1.rds")
seu_no_ins <- subset(t3, subset = Ins1 < 4 | Ins2 < 4)
saveRDS(seu_no_ins, file = "tcell_supposedly_free_of_cancer1.rds")
t4 <- seu_no_ins

## run pipeline
t4 %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = var_feat) %>%
  ScaleData() %>% 
  RunPCA(npcs = sig_dims, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:sig_dims, k.param = k_para,
                compute.SNN = TRUE, nn.method = "rann", nn.eps = 0) %>% 
  FindClusters(resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                              1.1, 1.2, 1.3, 1.4, 1.6)) %>% 
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
# seu <- SetIdent(object = FC, value = RESOLUTION)
# seu <- SetIdent(object = FC, value = 'RNA_snn_res.0.8')
seu <- SetIdent(object = FC, value = 'RNA_snn_res.1.2')
seu <- RunUMAP(seu, reduction = "pca", dims = 1:sig_dims)
pdf(file = paste0("res_", RESOLUTION, "_clusters_", proj_name, ".pdf"),
    width = 15, height = 8)
UMAPPlot(seu, reduction = "umap", pt.size = 0.1, label = TRUE,
         label.size = 5) + NoLegend() +
  labs(title = paste0("UMAP with ", RESOLUTION)) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

#### umap old labels
# pdf(file = paste0("umap_", RESOLUTION, "_oldlabels_", proj_name, ".pdf"),
#     width = 15, height = 8)
# DimPlot(seu, reduction = "umap", group.by = "initial_subtype_annotations",
#               pt.size = 0.1, label = TRUE, label.size = 3.5,  repel = TRUE)
# dev.off()

pdf(file = paste0("umap_doublets_", proj_name, ".pdf"), width = 8, height = 8)
DimPlot(seu, reduction = "umap", group.by = "pANNPredictions", pt.size = 0.1,
        label = FALSE, label.size = 2.0)
dev.off()

pdf(file = paste0("umap_old_labels_", proj_name, ".pdf"), width = 12, height = 8)
DimPlot(seu, reduction = "umap", group.by = "final_subtype_annotations",
        pt.size = 0.1, label = TRUE, label.size = 2.0)
dev.off()

## produce dimplots
pdf(paste0("umap_", proj_name, ".pdf"), width = 15, height = 8)
p1 <- DimPlot(seu, reduction = "umap", group.by = "orig.ident",
              pt.size = 0.1, label = FALSE, label.size = 2.0)
p2 <- DimPlot(seu, reduction = "umap", group.by = "Response",
              pt.size = 0.1, label = FALSE, label.size = 3.5,  repel = TRUE)
p3 <- DimPlot(seu, reduction = "umap", group.by = "Treatment",
              pt.size = 0.1, label = FALSE, label.size = 3.5,  repel = TRUE)
p4 <- DimPlot(seu, reduction = "umap", #group.by = "initial_subtype_annotations",
              pt.size = 0.1, label = TRUE, label.size = 3.5,  repel = TRUE)
plot_grid(p1, p2, p3, p4, ncol = 2)
dev.off()

# saveRDS(seu, file = "tcell_temp1.rds")
saveRDS(seu, file = "tcell_temp2.rds")

########################## tables for QC plots ################################
## table of number of cells 
# write.table(seu@meta.data, paste0("cells_metadata_", proj_name, ".txt"),
#             col.names = NA, sep = "\t")
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

################################ SingleR #######################################
# library("SingleR")
# library("celldex")
# mouse <- celldex::MouseRNAseqData()
# seuratObj <- seu
# singlerObj <- as.SingleCellExperiment(seuratObj)
# common <- intersect(rownames(singlerObj), rownames(mouse))
# mouse <- mouse[common,]
# singlerObj <- singlerObj[common,]
# singler <- SingleR(test = singlerObj, ref = mouse, labels = mouse$label.main,
#                    method = "cluster", ## change the resolution you chose
#                    clusters = seuratObj@meta.data[,"RNA_snn_res.1.8"],
#                    genes = "de")

## save singlr results 
# saveRDS(singler, file = paste0("singler_object_", proj_name, ".rds"))
# write.table(singler, file = paste0( "singler_", proj_name, ".txt"),
#             col.names = NA, sep = "\t")

############# Assigning cell type identity to clusters (annotation) ############
# seu <- SetIdent(object = seu, value = "initial_subtype_annotations")
## the cluster names are out of order
# cluster_ids <- singler@rownames
# singler_labels <- singler@listData$labels
# as.data.frame(cbind(cluster_ids, singler_labels)) %>%
#   dplyr::mutate(real_cluster_ids = as.integer(as.character(cluster_ids))) %>%
#   dplyr::arrange(real_cluster_ids) -> singler_annotations

# c("Cd8+ Exhausted", "Cd8+ Exhausted", "Cd8+ Naive", "Cd8+ Exhausted Th17+",
#   "Tregs", "Central Memory T-cells", "NK", "Proliferating T-cells",
#   "Cd4+ Th2+", "Cd8+ Exhausted", "Cd8+ Exhausted", "Gamma Delta T-cells",
#   "Proliferating T-cells") -> new.cluster.ids

c("Cd8+ Exhausted", "Unknown T-cells Th17+", "Tregs", "Unknown T-cells",
  "Cd4+ Th2+", "Cd8+ Memory", "NK cells", "Cd8+ Exhausted",
  "Proliferating T-cells", "Gamma Delta T-cells", "Cd8+ Exhausted"
  ) -> new.cluster.ids

## add column in metadata with the "Cluster" title and the cluster name
names(new.cluster.ids) <- levels(seu)
seu <- RenameIdents(seu, new.cluster.ids)
seu <- StashIdent(seu, save.name = "current_subtype_annotations")
head(seu@meta.data) ## verify new column

seu <- SetIdent(object = seu, value = "current_subtype_annotations")

pdf(paste0("Annotated_", proj_name, ".pdf"), width = 11, height = 8.5)
DimPlot(seu, reduction = "umap", label = FALSE, pt.size = 0.1) +
  labs(title = "Annotated UMAP") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf(paste0("Annotated_labels_", proj_name, ".pdf"), width = 11, height = 8.5) 
DimPlot(seu, reduction = "umap", label = TRUE, pt.size = 0.1, repel = TRUE) +
  labs(title = "Annotated UMAP") + theme(plot.title = element_text(hjust = 0.5))
dev.off() # the same with labels

saveRDS(seu, file = "tcell_recluster_actual3.rds")

########################### subclustering Endothelial #4444444444444444444444444
proj_name <- "endo_subcluster"
# file_path <- "/vsc-hard-mounts/leuven-data/341/vsc34181/RIPTAG_dorothyjr/high_filter_recluster/"
# endo1 <- readRDS(file = paste0(file_path, "endo_recluster_v3_correct.rds"))
# seu0 <- readRDS("final_annotations_recluster_v4_correct.rds")
# e1 <- subset(seu0, subset = final_primary_annotations == c("Endothelial cells"))
# e2 <- subset(seu0, subset = final_primary_annotations == c("Endothelial cells B20 PiK3i Resp"))
# e1 %>% merge(y = e2, project = proj_name) -> endo1
#endo1 <- e1
#endo1 <- readRDS("hf_rerun_endo.rds") ## round one
#endo1 <- readRDS(file = "endo_rerun_temp.rds")  ## round two
## clusters 4 and 5 are cancer and will be removed
# ce4 <- subset(endo1, subset = cluster == 4)
# ce5 <- subset(endo1, subset = cluster == 5)
# saveRDS(ce4, file = "endo_cancer_contam4.rds")
# saveRDS(ce5, file = "endo_cancer_contam5.rds")
# endo2 <- subset(endo1, subset = cluster != 4)
# endo3 <- subset(endo2, subset = cluster != 5)
# saveRDS(endo3, file = "endo_clean_temp.rds") ## 3224 cells
endo1 <- readRDS(file = "endo_rerun_temp2.rds")
cont_e <- subset(endo1, subset = final_subtype_annotations == "Endothelial cancer contamination")
saveRDS(cont_e, file = "endo_cancer_contam_endo.rds")
endo2 <- subset(endo1, subset = final_subtype_annotations != "Endothelial cancer contamination")
ce10 <- subset(endo2, subset = cluster == 10)
saveRDS(ce10, file = "endo_c10_doublets.rds")
endo3 <- subset(endo2, subset = cluster != 10)
saveRDS(endo3, file = "endo_clean_temp2.rds") ## 2810

## run pipeline
endo3 %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = var_feat) %>%
  ScaleData() %>% 
  RunPCA(npcs = sig_dims, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:sig_dims, k.param = k_para,
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
seu <- RunUMAP(seu, reduction = "pca", dims = 1:sig_dims)
pdf(file = paste0("res_", RESOLUTION, "_clusters_", proj_name, ".pdf"),
    width = 15, height = 8)
UMAPPlot(seu, reduction = "umap", pt.size = 0.1, label = TRUE,
         label.size = 5) + NoLegend() +
  labs(title = paste0("UMAP with ", RESOLUTION)) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

## produce dimplots
pdf(paste0("umap_", proj_name, "_old_labels2.pdf"), width = 15, height = 8)
p1 <- DimPlot(seu, reduction = "umap", group.by = "orig.ident",
              pt.size = 0.1, label = FALSE, label.size = 2.0)
p2 <- DimPlot(seu, reduction = "umap", group.by = "Response",
              pt.size = 0.1, label = FALSE, label.size = 3.5,  repel = TRUE)
p3 <- DimPlot(seu, reduction = "umap", group.by = "Treatment",
              pt.size = 0.1, label = FALSE, label.size = 3.5,  repel = TRUE)
p4 <- DimPlot(seu, reduction = "umap", #group.by = "final_subtype_annotations",
              pt.size = 0.1, label = TRUE, label.size = 3.5,  repel = TRUE)
plot_grid(p1, p2, p3, p4, ncol = 2)
dev.off()

pdf(file = paste0("umap_doublets_", proj_name, ".pdf"), width = 8, height = 8)
DimPlot(seu, reduction = "umap", group.by = "pANNPredictions", pt.size = 0.1,
        label = FALSE, label.size = 2.0)
dev.off()

pdf(file = paste0("umap_old_labels_", proj_name, ".pdf"), width = 12, height = 8)
DimPlot(seu, reduction = "umap", group.by = "final_subtype_annotations",
        pt.size = 0.1, label = TRUE, label.size = 2.0)
dev.off()

# saveRDS(seu, file = "endo_rerun_temp2.rds")
saveRDS(seu, file = "endo_rerun_temp2.rds")

########################## tables for QC plots ################################
## table of number of cells 
# write.table(seu@meta.data, paste0("cells_metadata_", proj_name, ".txt"),
#             col.names = NA, sep = "\t")
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

################################ annotate ######################################
# seu <- SetIdent(object = seu, value = "RNA_snn_res.1.4")
# 
# c("Cap arterial", "Cap Esm1+", "Cap Cd36+", "Tip cells Apln+",
#   "Cap Plpp3+", "Cap Plpp1+", "Tip cells Prcp+", "Pericytes",
#   "Cap interferon", "Cap Gm8817+", "Venous/HVC", "Artery",
#   "Post-cap Venule", "Lymphatic") -> new.cluster.ids

seu <- SetIdent(object = seu, value = "RNA_snn_res.0.8")

c("Tip cells", "Capillaries arterial", "Capillaries Oaz2+", "Capillaries Cd36+",
  "Arterial", "Pericytes", "Venous", "Capillaries interferon",  "Lymphatic"
  ) -> new.cluster.ids

## add column in metadata with the "Cluster" title and the cluster name
names(new.cluster.ids) <- levels(seu)
seu <- RenameIdents(seu, new.cluster.ids)
seu <- StashIdent(seu, save.name = "current_subtype_annotations")
head(seu@meta.data) ## verify new column

pdf(paste0("Annotated_", proj_name, ".pdf"), width = 11, height = 8.5)
DimPlot(seu, reduction = "umap", label = FALSE, pt.size = 0.1) +
  labs(title = "Annotated UMAP") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf(paste0("Annotated_labels_", proj_name, ".pdf"), width = 11, height = 8.5) 
DimPlot(seu, reduction = "umap", label = TRUE, pt.size = 0.1, repel = TRUE) +
  labs(title = "Annotated UMAP") + theme(plot.title = element_text(hjust = 0.5))
dev.off() # the same with labels

# saveRDS(seu, file = "endo_recluster3_actual.rds")
saveRDS(seu, file = "endo_recluster4_actual.rds")

############################ subclustering DCs #################################
# seu <- readRDS("macro_recluster_actual_2")
# pdcs <- subset(seu, subset = final_subtype_annotations == "pDCs")
# cdcs <- subset(seu, subset = final_subtype_annotations == "cDCs")
# pdcs %>% merge(y = cdcs, project = proj_name) -> DC
# 
# ## run pipeline
# DC %>% 
#   FindVariableFeatures(selection.method = "vst", nfeatures = var_feat) %>%
#   ScaleData(vars.to.regress = c("nCount_RNA", "percent.mt", "orig.ident")) %>% 
#   RunPCA(npcs = sig_dims, verbose = FALSE) %>%
#   FindNeighbors(reduction = "pca", dims = 1:sig_dims, k.param = k_para,
#                 compute.SNN = TRUE, nn.method = "rann", nn.eps = 0) %>% 
#   FindClusters(resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)) %>% 
#   RunUMAP(reduction = "pca", dims = 1:sig_dims) -> FC
# 
# #### umap old labels
# pdf(file = paste0("umap_", RESOLUTION, "_oldlabels_", proj_name, ".pdf"),
#     width = 15, height = 8)
# DimPlot(seu, reduction = "umap", group.by = "initial_subtype_annotations",
#               pt.size = 0.1, label = TRUE, label.size = 3.5,  repel = TRUE)
# dev.off()
# 
# ## write the column name of the chosen resolution with the SetIdent function
# ## Seurat will set the object with this column
# seu <- SetIdent(object = FC, value = RESOLUTION)
# seu <- RunUMAP(seu, reduction = "pca", dims = 1:sig_dims)
# pdf(file = paste0("res_", RESOLUTION, "_clusters_", proj_name, ".pdf"),
#     width = 15, height = 8)
# UMAPPlot(seu, reduction = "umap", pt.size = 0.1, label = TRUE,
#          label.size = 5) + NoLegend() +
#   labs(title = paste0("UMAP with ", RESOLUTION)) +
#   theme(plot.title = element_text(hjust = 0.5))
# dev.off()
# 
# ## produce dimplots
# pdf(paste0("umap_", proj_name, "_old_labels2.pdf"), width = 15, height = 8)
# p1 <- DimPlot(seu, reduction = "umap", group.by = "orig.ident",
#               pt.size = 0.1, label = FALSE, label.size = 2.0)
# p2 <- DimPlot(seu, reduction = "umap", group.by = "Response",
#               pt.size = 0.1, label = FALSE, label.size = 3.5,  repel = TRUE)
# p3 <- DimPlot(seu, reduction = "umap", group.by = "Treatment",
#               pt.size = 0.1, label = FALSE, label.size = 3.5,  repel = TRUE)
# p4 <- DimPlot(seu, reduction = "umap", group.by = "final_subtype_annotations",
#               pt.size = 0.1, label = TRUE, label.size = 3.5,  repel = TRUE)
# plot_grid(p1, p2, p3, p4, ncol = 2)
# dev.off()
# 
# ####################### subclustering neutrophils ##############################
# seu <- readRDS("~/Documents/tmp/pancreatic_melanie/lf_neutro_rerun/lf_myeloid_neutros.rds")
# 
# ## run pipeline
# seu %>% 
#   FindVariableFeatures(selection.method = "vst", nfeatures = var_feat) %>%
#   ScaleData() %>% 
#   RunPCA(npcs = sig_dims, verbose = FALSE) %>%
#   FindNeighbors(reduction = "pca", dims = 1:sig_dims, k.param = k_para,
#                 compute.SNN = TRUE, nn.method = "rann", nn.eps = 0) %>% 
#   FindClusters(resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
#                               1.1, 1.2)) %>% 
#   RunUMAP(reduction = "pca", dims = 1:sig_dims) -> FC
# 
# ## multiple resolution plots
# pdf(file = paste0("clusters_diff_res_", proj_name, ".pdf"), width = 15, height = 8)
# Seurat_02 <- SetIdent(object = FC, value = 'RNA_snn_res.0.2')
# Seurat02 <- RunUMAP(Seurat_02, reduction = "pca", dims = 1:sig_dims)
# UMAPPlot(Seurat02, reduction = "umap", pt.size = 0.25, label = TRUE) +
#   labs(title = "UMAP with res 0.2") + theme(plot.title = element_text(hjust = 0.5))
# Seurat_04 <- SetIdent(object = FC, value = 'RNA_snn_res.0.4')
# Seurat04 <- RunUMAP(Seurat_04, reduction = "pca", dims = 1:sig_dims)
# UMAPPlot(Seurat04, reduction = "umap", pt.size = 0.25, label = TRUE) +
#   labs(title = "UMAP with res 0.4") + theme(plot.title = element_text(hjust = 0.5))
# Seurat_06 <- SetIdent(object = FC, value = 'RNA_snn_res.0.6')
# Seurat06 <- RunUMAP(Seurat_06, reduction = "pca", dims = 1:sig_dims)
# UMAPPlot(Seurat06, reduction = "umap", pt.size = 0.25, label = TRUE) +
#   labs(title = "UMAP with res 0.6") + theme(plot.title = element_text(hjust = 0.5))
# Seurat_08 <- SetIdent(object = FC, value = 'RNA_snn_res.0.8')
# Seurat08 <- RunUMAP(Seurat_08, reduction = "pca", dims = 1:sig_dims)
# UMAPPlot(Seurat08, reduction = "umap", pt.size = 0.25, label = TRUE) +
#   labs(title = "UMAP with res 0.8") + theme(plot.title = element_text(hjust = 0.5))
# Seurat_1 <- SetIdent(object = FC, value = 'RNA_snn_res.1.0')
# Seurat1 <- RunUMAP(Seurat_1, reduction = "pca", dims = 1:sig_dims)
# UMAPPlot(Seurat1, reduction = "umap", pt.size = 0.25, label = TRUE) + 
#   labs(title = "UMAP with res 1.0") + theme(plot.title = element_text(hjust = 0.5))
# Seurat_12 <- SetIdent(object = FC, value = 'RNA_snn_res.1.2')
# Seurat12 <- RunUMAP(Seurat_12, reduction = "pca", dims = 1:sig_dims)
# UMAPPlot(Seurat12, reduction = "umap", pt.size = 0.25, label = TRUE) + 
#   labs(title = "UMAP with res 1.2") + theme(plot.title = element_text(hjust = 0.5))
# dev.off() 
# 
# ## write the column name of the chosen resolution with the SetIdent function
# ## Seurat will set the object with this column
# # seu <- SetIdent(object = FC, value = RESOLUTION)
# seu <- SetIdent(object = FC, value = RESOLUTION)
# seu <- RunUMAP(seu, reduction = "pca", dims = 1:sig_dims)
# pdf(file = paste0("res_", RESOLUTION, "_clusters_", proj_name, ".pdf"),
#     width = 15, height = 8)
# UMAPPlot(seu, reduction = "umap", pt.size = 0.1, label = TRUE,
#          label.size = 5) + NoLegend() +
#   labs(title = paste0("UMAP with ", RESOLUTION)) +
#   theme(plot.title = element_text(hjust = 0.5))
# dev.off()
# 
# #### umap old labels
# # pdf(file = paste0("umap_", RESOLUTION, "_oldlabels_", proj_name, ".pdf"),
# #     width = 15, height = 8)
# # DimPlot(seu, reduction = "umap", group.by = "final_subtype_annotations",
# #               pt.size = 0.1, label = TRUE, label.size = 3.5,  repel = TRUE)
# # dev.off()
# 
# pdf(file = paste0("umap_doublets_", proj_name, ".pdf"), width = 8, height = 8)
# DimPlot(seu, reduction = "umap", group.by = "pANNPredictions", pt.size = 0.1,
#         label = FALSE, label.size = 2.0)
# dev.off()
# 
# ## produce dimplots
# pdf(paste0("umap_", proj_name, ".pdf"), width = 15, height = 8)
# p1 <- DimPlot(seu, reduction = "umap", group.by = "orig.ident",
#               pt.size = 0.1, label = FALSE, label.size = 2.0)
# p2 <- DimPlot(seu, reduction = "umap", group.by = "Response",
#               pt.size = 0.1, label = FALSE, label.size = 3.5,  repel = TRUE)
# p3 <- DimPlot(seu, reduction = "umap", group.by = "Treatment",
#               pt.size = 0.1, label = FALSE, label.size = 3.5,  repel = TRUE)
# p4 <- DimPlot(seu, reduction = "umap", #group.by = "initial_subtype_annotations",
#               pt.size = 0.1, label = TRUE, label.size = 3.5,  repel = TRUE)
# plot_grid(p1, p2, p3, p4, ncol = 2)
# dev.off()
# 
# ### save the clustered myeloid object
# saveRDS(seu, file = "lf_myeloid_neutro_subcluster.rds")
# n1 <- subset(seu, subset = RNA_snn_res.0.6 == 2)
# saveRDS(n1, file = "lf_neutro_actual.rds")
# proj_name <- "neutro"
# 
# n1 <- SetIdent(object = n1, value = "Response")
# 
# n1 %>% 
#   FindVariableFeatures(selection.method = "vst", nfeatures = var_feat) %>%
#   ScaleData() %>%
#   RunPCA(npcs = sig_dims, verbose = FALSE) %>%
#   FindNeighbors(reduction = "pca", dims = 1:sig_dims, k.param = k_para,
#                 compute.SNN = TRUE, nn.method = "rann", nn.eps = 0) %>% 
#   FindClusters(resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)) %>% 
#   RunUMAP(reduction = "pca", dims = 1:sig_dims) -> FC
# 
# seu <- SetIdent(object = FC, value = RESOLUTION)
# 
# pdf(paste0("umap_", proj_name, ".pdf"), width = 12, height = 10)
# p1 <- DimPlot(seu, reduction = "umap", group.by = "orig.ident",
#               pt.size = 0.1, label = FALSE, label.size = 2.0)
# p2 <- DimPlot(seu, reduction = "umap", group.by = "Response",
#               pt.size = 0.1, label = FALSE, label.size = 3.5,  repel = TRUE)
# p3 <- DimPlot(seu, reduction = "umap", group.by = "Treatment",
#               pt.size = 0.1, label = FALSE, label.size = 3.5,  repel = TRUE)
# p4 <- DimPlot(seu, reduction = "umap", #group.by = "initial_subtype_annotations",
#               pt.size = 0.1, label = TRUE, label.size = 3.5,  repel = TRUE)
# plot_grid(p1, p2, p3, p4, ncol = 2)
# dev.off()
# 
# seu <- SetIdent(object = seu, value = "Response")
# 
# seu@meta.data$Treatment_Response <- paste0(seu@meta.data$Treatment, "_", seu@meta.data$Response)
# seu@meta.data$Treatment_Response <- gsub("Untreated_", "", seu@meta.data$Treatment_Response)
# 
# seu <- SetIdent(object = seu, value = "Treatment_Response")
# 
# saveRDS(seu, file = "lf_neutro_actual_metadata.rds")
# seu <- readRDS("lf_neutro_actual_metadata.rds")
# seu <- SetIdent(object = seu, value = "Response")
# 
# ## neutrophil volcano plots
# neutro_comp <- FindMarkers(seu, ident.1 = "Responding", ident.2 = "Relapsing")
# 
# neutro_comp %>%
#   tibble::rownames_to_column() %>%
#   dplyr::rename(gene = rowname) %>%
#   dplyr::filter(gene != "Ins1",
#                 gene != "Ins2",
#                 gene != "Iapp",
#                 gene != "Cela2a",
#                 gene != "Cela3b",
#                 gene != "Cela1",
#                 gene != "Ppy",
#                 gene != "Try4",
#                 gene != "Try5",
#                 gene != "Prss2",
#                 gene != "Reg1",
#                 gene != "Cpb1",
#                 gene != "Pnlip",
#                 gene != "Clps",
#                 gene != "Ctrb1",
#                 gene != "Rnase1",
#                 gene != "2210010C04Rik",
#                 gene != "Nnat",
#                 gene != "Chga",
#                 gene != "Pcsk1n",
#                 gene != "Sycn") -> neutro_clean
# 
# res_final1 <- neutro_clean[!is.na(neutro_clean$p_val_adj), ]
# res_final1$threshold <- as.factor(abs(res_final1$avg_log2FC) > 1 &
#                                       res_final1$p_val_adj < 0.05)
# 
# pdf(paste0("Vol_plot_neutrophils_resp_rel_", proj_name, ".pdf"))
# ggplot2::ggplot(data = res_final1,
#                 aes(x = avg_log2FC, y = -log10(p_val_adj), color = threshold)) +
#                 geom_text(label = res_final1$gene, nudge_x = 0.1,
#                           nudge_y = 0.1, check_overlap = TRUE) +
#                 theme(legend.position = "none") +
#                 geom_point(alpha = 0.4, size = 0.5) +
#                 xlab("log2 fold change") + ylab("-log10 padj-value") +
#                 theme(plot.title = element_text(hjust = 0.5)) +
#                 scale_color_manual(values = c("#000000", "#FF0000")) -> p1
# print(p1)
# dev.off()
# #
# res_final1 <- neutro_clean[!is.na(neutro_clean$p_val), ]
# res_final1$threshold <- as.factor(abs(res_final1$avg_log2FC) > 1 &
#                                       res_final1$p_val < 0.05)
# 
# pdf(paste0("Vol_plot_neutrophils_resp_rel_pval_", proj_name, ".pdf"))
# ggplot2::ggplot(data = res_final1,
#                 aes(x = avg_log2FC, y = -log10(p_val), color = threshold)) +
#                 geom_text(label = res_final1$gene, nudge_x = 0.1,
#                           nudge_y = 0.1, check_overlap = TRUE) +
#                 theme(legend.position = "none") +
#                 geom_point(alpha = 0.4, size = 0.5) +
#                 xlab("log2 fold change") + ylab("-log10 p-value") +
#                 theme(plot.title = element_text(hjust = 0.5)) +
#                 scale_color_manual(values = c("#000000", "#FF0000")) -> p1
# print(p1)
# dev.off()

# library("Seurat")
# file_path <- "/data/leuven/341/vsc34181/RIPTAG_dorothyjr/sc_data/dorothy_files/new_rds/"
# seu <- readRDS(file = paste0(file_path, "cancer_actual2.rds"))
# c1 <- subset(seu, subset = Treatment == "Untreated") ## only untreated
# c2 <- subset(seu, subset = Treatment != "Untreated") ## all treatments EXCEPT untreated
# saveRDS(c1, file = "cancer_untreated.rds")
# saveRDS(c2, file = "cancer_all_non_untreated.rds")
