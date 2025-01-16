############################### sc-pipeline ####################################
## this script reads the pre-processed Seurat objects taken from the previous
## script `1_preprocess.R` and filters the objects further, normalizes the data
## and finds variable features. running PCA, UMAP and generating various
## diagnostic plots. Lastly, the clusters resulting from the UMAP are annotated 
## sample and response type, filters the cells, produces QC plots,
## NOTE: if you need to run harmony, you need to use R 3.6.3 libraries
## /media/seq-srv-05/vrc/Project/Project_Gino/software/R/R-3.6.3
## /data/projects/Project_Theo/R/R-4.0.3

## load libraries
library("R.utils")
library("dplyr")
library("Seurat")
library("ggplot2")
library("future")
library("cowplot")
## library("harmony")

## parallelize workflow
plan("multiprocess", workers = 50) # uses 50 CPU
options(future.globals.maxSize = 14000 * 1024^2) ## 5GB per worker

## evergreen variables
proj_name <- "pancreatic_mice" ## name of project
proj_version <- "unf2"       ## name of version of project
# server_path <- "/media/seq-srv-05/vrc/Project/Project_Theo/pancreatic_mice/unfiltered/"
server_path <- "/media/seq-srv-05/vrc/Project/Project_Theo/pancreatic_mice/raw_1000dr/"
mito <- "^mt-" ## "^mt-" if mouse, "^MT-" if human
up_feat <- 6000 ## upper nFeature threshold
low_feat <- 1000 ## lower nFeature threshold default 400
var_feat <- 2000 ## variableFeature threshold
low_count <- 1000 ## lower nCount threshold default 400
p_mt <- 20 ## percent mito threshold
NPC <- 35 ## number of PCs of first PCA
k_para <- 50 ## number of k paramaters in FindNeighbors default is 30
sig_dims <- 12 ## number of PCs of final PCA
resolution <- 1.4 ## resolution
RESOLUTION <- 'RNA_snn_res.1.4'; ## resolution string

##############################################################################
## this chunk of code reads the processed Seurat files for each sample 

## read files into list (if running the saved files on the server)
file_list <- list.files(server_path, pattern = "processed_")
fls <- lapply(file_list, function(i) readRDS(file = i))
names(fls) <- gsub("processed_", "", gsub(".rds", "", unlist(file_list)))

##############################################################################
## this chunk of code only needs to be ran to prepare unfiltered metadata to
## determine filter thresholds in the next step. The saved file from this code
## chunk can be used to generate QC plots in the `3_qc_plots.R` script.

# fls[[1]] %>% merge(y = fls[2:length(fls)], project = proj_name,
#                    add.cell.ids = names(fls)) -> merged_seurat
# seu <- readRDS(file = paste0("FINAL_pancreatic_mice.rds"))
# seu <- readRDS(file = paste0("seu_", proj_name, ".rds"))
# merged_seurat <- seu ## or if you already have this object

## Compute number of genes per UMI for each cell and percent mito ratio
# merged_seurat@meta.data$log10GenesPerUMI <- log10(merged_seurat@meta.data$nFeature_RNA) /
#                                             log10(merged_seurat@meta.data$nCount_RNA)
# merged_seurat@meta.data$percent.mt <- PercentageFeatureSet(object = merged_seurat,
#                                                            pattern = mito)
# 
# ## Create metadata dataframe from merged Seurat object and save .rds for QC plots
# merged_seurat@meta.data %>%
#   dplyr::select(orig.ident, nCount_RNA, nFeature_RNA, Response, Treatment,
#                 percent.mt, pANN, pANNPredictions, log10GenesPerUMI) %>%
#   dplyr::mutate(mito_ratio = percent.mt / 100) %>%
#   dplyr::rename(sample = orig.ident, nUMI = nCount_RNA, nGene = nFeature_RNA) %>%
#   dplyr::mutate(mito_ratio = mito_ratio$nCount_RNA,
#                 percent_mt = percent.mt$nCount_RNA) %>%
#   dplyr::select(-percent.mt) -> metadata
# metadata$cells <- rownames(metadata) ## add cell IDs  metadata
# saveRDS(metadata, file = paste0(proj_version, "_metadata.rds"))

##############################################################################
## this chunk of code reads the processed Seurat files for each sample and 
## prepares the metadata for QC analysis

## merge dataset, and add cell-ids and subset on parameters
fls[[1]] %>%
merge(y = fls[2:length(fls)], project = proj_name, add.cell.ids = names(fls)) %>%
  subset(subset = nFeature_RNA < up_feat &
                  nFeature_RNA > low_feat &
                  nCount_RNA > low_count &
                  percent.mt < p_mt) -> merged_seu

## remove doublets
filtered_seurat <- subset(merged_seu, subset = pANNPredictions == c("Singlet"))

## whole high filter object
# An object of class Seurat 
# 29043 features across 92346 samples within 1 assay 
# Active assay: RNA (29043 features, 0 variable features)

## no doublets high filter object
# An object of class Seurat 
# 29043 features across 96338 samples within 1 assay 
# Active assay: RNA (29043 features, 0 variable features)

## get counts, output a logical vector for every gene on whether the more than
## zero counts per UMI. Sums all TRUE values and returns TRUE if more than 10
## TRUE values per gene, only keeping those genes expressed in more than 10,
## cells. Then assign to filtered Seurat object
# counts <- GetAssayData(object = merged_seu, slot = "counts")
# nonzero <- counts > 0
# keep_genes <- Matrix::rowSums(nonzero) >= 10
# filtered_counts <- counts[keep_genes, ]
# filtered_seurat <- CreateSeuratObject(filtered_counts,
#                                       meta.data = merged_seu@meta.data)

# filtered_seurat <- merged_seu
filtered_seurat %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = var_feat) %>%
  ScaleData(vars.to.regress = c("nCount_RNA", "percent.mt", "orig.ident")) %>% 
  RunPCA(npcs = NPC, verbose = FALSE) -> SC

## Jackstraw
# SC <- JackStraw(SC, num.replicate = 100, dims = NPC)
# SC <- ScoreJackStraw(SC, dims = 1:NPC)
# pdf(paste0("jackstraw_", proj_name, ".pdf"), width = 15, height = 8)
# JackStrawPlot(SC, dims = 1:NPC)
# dev.off()
# 
# ## elbow plot visualizing meaningful PCs
# pdf(paste0("elbowplot_", proj_name, ".pdf"), width = 15, height = 8)
# ElbowPlot(SC, ndims = NPC, reduction = "pca")
# dev.off()

## how to find the elbow point
# library("PCAtools")
# dims <- Stdev(SC, reduction = "pca")
# elbow <- PCAtools::findElbowPoint(dims)

## Heatmap
# pdf(paste0("PC_heatmap_", proj_name, ".pdf"), width = 15, height = 8)
# # DimHeatmap(SC, dims = 1:NPC, cells = 2000, balanced = TRUE)
# DimHeatmap(SC, dims = 1, cells = 2000, balanced = TRUE)
# DimHeatmap(SC, dims = 2, cells = 2000, balanced = TRUE)
# DimHeatmap(SC, dims = 3, cells = 2000, balanced = TRUE)
# DimHeatmap(SC, dims = 4, cells = 2000, balanced = TRUE)
# DimHeatmap(SC, dims = 5, cells = 2000, balanced = TRUE)
# DimHeatmap(SC, dims = 6, cells = 2000, balanced = TRUE)
# DimHeatmap(SC, dims = 7, cells = 2000, balanced = TRUE)
# DimHeatmap(SC, dims = 8, cells = 2000, balanced = TRUE)
# DimHeatmap(SC, dims = 9, cells = 2000, balanced = TRUE)
# DimHeatmap(SC, dims = 10, cells = 2000, balanced = TRUE)
# DimHeatmap(SC, dims = 11, cells = 2000, balanced = TRUE)
# DimHeatmap(SC, dims = 12, cells = 2000, balanced = TRUE)
# DimHeatmap(SC, dims = 13, cells = 2000, balanced = TRUE)
# DimHeatmap(SC, dims = 14, cells = 2000, balanced = TRUE)
# DimHeatmap(SC, dims = 15, cells = 2000, balanced = TRUE)
# DimHeatmap(SC, dims = 16, cells = 2000, balanced = TRUE)
# DimHeatmap(SC, dims = 17, cells = 2000, balanced = TRUE)
# DimHeatmap(SC, dims = 18, cells = 2000, balanced = TRUE)
# DimHeatmap(SC, dims = 19, cells = 2000, balanced = TRUE)
# DimHeatmap(SC, dims = 20, cells = 2000, balanced = TRUE)
# dev.off()

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
  FindNeighbors(reduction = "pca", dims = 1:sig_dims, k.param = k_para,
                compute.SNN = TRUE, nn.method = "rann", nn.eps = 0) %>% 
  FindClusters(resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                              1.1, 1.2, 1.3, 1.4, 1.6, 1.8, 2.0, 3.0)) %>% 
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

## save the results
# saveRDS(seu, file = paste0("res_", proj_name, ".rds"))
# seu <- readRDS(file = paste0("res_", proj_name, ".rds"))

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
# f <- ordered(features$RNA_snn_res.0.8, levels = c(0:clust_num))
f <- ordered(features$RNA_snn_res.1.4,
             levels = c(0:(max(as.integer(seu@meta.data$RNA_snn_res.1.4)) - 1)))

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

################################ SingleR #######################################
## SingleR is Unsupervised clustering, which affords differential expression
## between 2 condition in a cluster - for each cluster. we usually use marker
## genes to annotate our clusters, but what if the marker genes are unknown,
## or an automated annotation is desired to compare for confirmation, SingleR is
## used. Basically, there are 2 type of annotations: cluster annotation or cell
## annotation
library("SingleR")
library("celldex")
# DefaultAssay(seu) <- def_assay
# head(seu@meta.data)

## select desired resolution. this is easier after the FindCluster function,
## where the headers of the resolutions are renamed in the metadata into
## "RNA...."you should have "Integrated_..."
# levels(seu@meta.data$RNA_snn_res.0.4) 
# levels(seu@meta.data$RNA_snn_res.1.2)
# levels(seu@meta.data$RNA_snn_res.0.8) 

## human
# encode = BlueprintEncodeData()
## mouse
# mouse = MouseRNAseqData()
mouse <- celldex::MouseRNAseqData()

## HUMAN
# seuratObj = seu
# singlerObj <- as.SingleCellExperiment(seuratObj)
# common <- intersect(rownames(singlerObj), rownames(encode))
# encode <- encode[common,]
# singlerObj <- singlerObj[common,]
## here, only change the resolution you chose
# singler <- SingleR(test = singlerObj, ref = encode, labels = encode$label.main,
#                   method = "cluster",
#                   clusters = seuratObj@meta.data[,"RNA_snn_res.1.4"], genes = "de")

## MOUSE
seuratObj <- seu
singlerObj <- as.SingleCellExperiment(seuratObj)
common <- intersect(rownames(singlerObj), rownames(mouse))
mouse <- mouse[common,]
singlerObj <- singlerObj[common,]
singler <- SingleR(test = singlerObj, ref = mouse, labels = mouse$label.main,
                   method = "cluster", ## change the resolution you chose
                   clusters = seuratObj@meta.data[,"RNA_snn_res.0.8"],
                   # clusters = seuratObj@meta.data[,"RNA_snn_res.1.4"],
                   genes = "de") 

## These results table can be opened in excel. it displays the list of marker
## genes per cluster, labels of cluster (like T-cells, macrophages, etc...)
## Labels are shown before fine-tuning (first.labels), after fine-tuning
## (labels) and after pruning (pruned.labels) but it is not super detailed, so
## don't use these annotations for subclustering, it may yield undesirable
## results. always use your own marker genes to really annotate the clusters,
## but it gives a nice confirmation (if you already plotted your markers) or
## idea (if you didn't plot your markers yet) of what to expect
# save(singler, file = 'singler_object.RData')
saveRDS(singler, file = paste0("singler_object_", proj_name, ".rds"))
write.table(singler, file = paste0( "singler_", proj_name, ".txt"),
            col.names = NA, sep = "\t")

################################ singler object ################################ 
# singler <- readRDS(file = paste0("singler_object_", proj_name, ".rds"))

## the resolution you chose - column of your metadata - note! with integrated
## analysis, the name will be "integrated_snn_res.1.4" for example
table(seu@meta.data$RNA_snn_res.1.4)
## seu <- SetIdent(object = seu, value = RESOLUTION)

## generate a nice table for your presentation ################################
library("dplyr")
library("plyr")
library("ggplot2")

c = ddply(seu@meta.data, ~ RNA_snn_res.1.4, summarise,
# c = ddply(seu@meta.data, ~ RNA_snn_res.0.4, summarise,
          mean = mean(as.numeric(pANN)))
c

y = table(seu@active.ident, seu@meta.data$pANNPredictions)
y

y = as.data.frame.matrix(y)
w = cbind(c, Doublet = y$Doublet)
w = cbind(w, Singlet = y$Singlet)

w$TOTAL <- w$Doublet + w$Singlet
w$pct_Doublet <- w$Doublet / w$TOTAL
w$pct_Singlet <- w$Singlet / w$TOTAL
w

library("tidyverse")
w = w %>% dplyr::rename(mean_pANN_score = mean)
w = w[, c(1, 2, 4, 7, 3, 6, 5)]
w

write.table(w, file = paste0("pann_per_C_doublet", proj_name, ".txt"),
            col.names = NA, sep = "\t")

# read.table(file = paste0("pann_per_C_doublet", proj_name, ".txt"),
#            header = TRUE, sep = "\t") -> a1

# saveRDS(seu, file = paste0("seu_", proj_name, ".rds"))

## make a plot - 1 color for doublets and 1 for singlets
# pdf(paste0("doublets_", proj_name, ".pdf"))
# DimPlot(seu, pt.size = 0.1, label = FALSE, label.size = 0,
#         reduction = "umap", group.by = "pANNPredictions") + theme(aspect.ratio = 1)
# dev.off()

############# Assigning cell type identity to clusters (annotation) ############
#seu <- readRDS(file = paste0("seu_", proj_name, ".rds"))
seu <- SetIdent(object = seu, value = RESOLUTION)

## the cluster names are out of order
cluster_ids <- singler@rownames
singler_labels <- singler@listData$labels
as.data.frame(cbind(cluster_ids, singler_labels)) %>%
  dplyr::mutate(real_cluster_ids = as.integer(as.character(cluster_ids))) %>% 
  dplyr::arrange(real_cluster_ids) -> singler_annotations

## clusters at this resolution, here we are manually annotating them
# c("Cancer", "Cancer", "Cancer", "Cancer", "Cancer", "Cancer", "Macrophages",
#   "Cancer", "Cancer", "Acinar", "Cancer", "Endothelial", "T-cells/NK", "B-Cell",
#   "Cancer", "Cancer", "T-cells/NK", "Cancer", "Cancer", "Cancer", "Cancer",
#   "Cancer", "Monocytes", "Macrophages", "Cancer", "Cancer", "Acinar",
#   "T-cells/NK", "Fibroblasts", "Ductal/Stem", "Endothelial", "Cancer", "Cancer",
#   "Cancer", "Cancer", "Macrophages", "Cancer", "B-Cell", "Endothelial",
#   "Monocytes", "Cancer", "Smooth Muscle", "B-Cell", "Cancer",
#   "T-cells/NK") -> manual_labels

#   cluster_ids    singler_labels real_cluster_ids
# 1            0           Neurons                0
# 2            1           Neurons                1
# 3            2           Neurons                2
# 4            3           Neurons                3
# 5            4           Neurons                4
# 6            5           Neurons                5
# 7            6           Neurons                6
# 8            7           Neurons                7
# 9            8           Neurons                8
# 10           9       Macrophages                9
# 11          10           Neurons               10
# 12          11           Neurons               11
# 13          12           Neurons               12
# 14          13           Neurons               13
# 15          14       Macrophages               14
# 16          15           T cells               15
# 17          16       Macrophages               16
# 18          17 Endothelial cells               17
# 19          18           Neurons               18
# 20          19           Neurons               19
# 21          20           Neurons               20
# 22          21           B cells               21
# 23          22           Neurons               22
# 24          23           Neurons               23
# 25          24           T cells               24
# 26          25           Neurons               25
# 27          26           Neurons               26
# 28          27           Neurons               27
# 29          28           Neurons               28
# 30          29           Neurons               29
# 31          30 Endothelial cells               30
# 32          31           Neurons               31
# 33          32           Neurons               32
# 34          33           Neurons               33
# 35          34           Neurons               34
# 36          35           Neurons               35
# 37          36       Fibroblasts               36
# 38          37           Neurons               37
# 39          38           Neurons               38
# 40          39           Neurons               39
# 41          40         Monocytes               40
# 42          41       Fibroblasts               41
# 43          42           Neurons               42
# 44          43           B cells               43
# 45          44           Neurons               44
# 46          45           Neurons               45
# 47          46 Endothelial cells               46
# 48          47           Neurons               47
# 49          48         Monocytes               48
# 50          49           Neurons               49
# 51          50 Endothelial cells               50
# 52          51           Neurons               51
# 53          52          NK cells               52
# 54          53           Neurons               53
# 55          54       Macrophages               54
# 56          55       Macrophages               55

c("Neurons", "Neurons", "Neurons", "Neurons", "Neurons", "Neurons", "Neurons",
  "Neurons", "Neurons", "Macrophages", "Neurons", "Neurons", "Neurons",
  "Neurons", "Macrophages", "T cells", "Macrophages", "Endothelial cells",
  "Neurons", "Neurons", "Neurons", "B cells", "Neurons", "Neurons", "T cells",
  "Neurons", "Neurons", "Neurons", "Neurons", "Neurons", "Endothelial cells",
  "Neurons", "Neurons", "Neurons", "Neurons", "Neurons", "Fibroblasts",
  "Neurons", "Neurons", "Neurons", "Monocytes", "Fibroblasts", "Neurons",
  "B cells", "Neurons", "Neurons", "Endothelial cells", "Neurons", "Monocytes",
  "Neurons", "Endothelial cells", "Neurons", "NK cells", "Neurons",
  "Macrophages", "Macrophages") -> manual_labels

# cluster_ids <- as.character(0:(max(as.integer(seu@meta.data$RNA_snn_res.0.4)) - 1))
cluster_ids <- as.character(0:(max(as.integer(seu@meta.data$RNA_snn_res.1.4)) - 1))
manual_annotations <- as.data.frame(cbind(cluster_ids, manual_labels))

# manual_annotations %>%
#   dplyr::full_join(singler_annotations, by = "cluster_ids") -> new_cluster_ids
singler_annotations %>%
  dplyr::full_join(manual_annotations, by = "cluster_ids") -> new_cluster_ids
new_cluster_ids %>% 
  write.table(file = paste0("cluster_ids_", proj_name, ".txt"), col.names = NA,
              sep = "\t")
new.cluster.ids <- new_cluster_ids$manual_labels

## unfiltered2 @RES = 1.4
c("Cancer", "Cancer", "Cancer", "Cancer", "Cancer", "Cancer", "Cancer",
  "Cancer", "Cancer", "Macrophages", "Cancer", "Cancer", "Cancer", "Cancer",
  "Macrophages", "T cells", "Macrophages", "Endothelial cells", "Cancer",
  "Cancer", "Cancer", "B cells", "Cancer", "Cancer", "T cells", "Cancer",
  "Cancer", "Cancer", "Cancer", "Cancer", "Endothelial cells", "Cancer",
  "Cancer", "Cancer", "Cancer", "Cancer", "Fibroblasts", "Cancer", "Cancer",
  "Cancer", "Monocytes", "Fibroblasts", "Cancer", "B cells", "Cancer", "Cancer",
  "Endothelial cells", "Cancer", "Monocytes", "Cancer", "Endothelial cells",
  "Cancer", "NK cells", "Cancer", "Macrophages", "Macrophages") -> new.cluster.ids

## add column in metadata with the "Cluster" title and the cluster name
names(new.cluster.ids) <- levels(seu)
seu <- RenameIdents(seu, new.cluster.ids)
seu <- StashIdent(seu, save.name = "Cluster")
head(seu@meta.data) ## verify new column

pdf(paste0("Annotated_v2", proj_name, ".pdf"), width = 11, height = 8.5)
DimPlot(seu, reduction = "umap", label = FALSE, pt.size = 0.1) +
  labs(title = "Annotated UMAP") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf(paste0("Annotated_labels_v2", proj_name, ".pdf")) 
DimPlot(seu, reduction = "umap", label = TRUE, pt.size = 0.1, repel = TRUE) +
  labs(title = "Annotated UMAP") + theme(plot.title = element_text(hjust = 0.5))
dev.off() # the same with labels

## Expression Heatmap
clust_avgs <- AverageExpression(seu, return.seurat = TRUE)
pdf(paste0("Exp_Heatmap_5_", proj_name, ".pdf"))
DoHeatmap(clust_avgs, features = unlist(TopFeatures(seu[["pca"]], balanced = TRUE)),
          size = 5, draw.lines = FALSE)
dev.off()

saveRDS(seu, file = paste0("FINAL_latest_ann_2_", proj_name, ".rds"))
# saveRDS(seu, file = paste0("FINAL_cancer_ann_", proj_name, ".rds"))
# saveRDS(seu, file = paste0("FINAL_", proj_name, ".rds"))
## seu <- readRDS(file = paste0("FINAL_cancer_ann_", proj_name, ".rds")) ### cstc
# seu <- readRDS(file = paste0("FINAL_latest_ann_2_pancreatic_mice.rds")) ### cstc
# seu <- readRDS(file = paste0("FINAL_", proj_name, ".rds"))

# seu1 <- subset(seu, subset = Cluster != "Fibroblasts")
# seu2 <- subset(seu1, subset = Cluster != "Smooth Muscle")
# 
seu@meta.data$new_labels <- as.character(seu@meta.data$Cluster)
seu@meta.data$new_labels <- gsub("Cancer_B20_PDL1_CTLA4_Resp", "Cancer", seu@meta.data$new_labels)
seu@meta.data$new_labels <- gsub("Cancer_B20_PDL1_Resp1", "Cancer", seu@meta.data$new_labels)
seu@meta.data$new_labels <- gsub("Cancer_B20_PDL1_Resp2", "Cancer", seu@meta.data$new_labels)
seu@meta.data$new_labels <- gsub("Cancer_B20_Pi3Ki_Resp", "Cancer", seu@meta.data$new_labels)
seu@meta.data$new_labels <- gsub("_", " ", seu@meta.data$new_labels)
seu <- SetIdent(object = seu, value = 'new_labels')
#  [1] "Cancer"                     "Macrophages"               
#  [3] "Acinar"                     "Cancer_B20_PDL1_CTLA4_Resp"
#  [5] "Endothelial"                "T-cells/NK"                
#  [7] "B-Cell"                     "Cancer_B20_PDL1_Resp1"     
#  [9] "Cancer_B20_PDL1_Resp2"      "Monocytes"                 
# [11] "Cancer_B20_Pi3Ki_Resp"      "Fibroblasts"               
# [13] "Ductal/Stem"                "Smooth Muscle"    


## raw_400
# c("Cancer", "Cancer", "Macro/Myeloid/Neutro", "Acinar/Erythrocytes", "Acinar",
#   "Cancer", "Endothelial cells", "Cancer", "Cancer", "T cells", "Cancer",
#   "B cells", "T cells", "Cancer", "Cancer", "Cancer", "Cancer",
#   "Fibroblasts/Smooth Muscle", "Cancer", "Cancer", "Endothelial cells", 
#   "Monocytes", "B cells", "Cancer") -> new.cluster.ids

## raw_1000
# c("Cancer", "Cancer", "Cancer", "Macro/Myeloids/Neutro", "Acinar", "Cancer",
#   "Endothelial cells", "Cancer", "T cells", "B cells", "T cells", "Cancer",
#   "Cancer", "Cancer/Erythrocytes", "EC/Fibro/Smooth Muscle",
#   "Fibroblasts", "Cancer", "Monocytes", "Cancer", "B cells") -> new.cluster.ids

## raw_1000dr @1.4
# Cluster 10 – CCB095, Cancer B20+aPDL1+aCTLA4 – Responding (3182 cells)
# Cluster 15 – CCB030, Cancer B20+aPDL1– Responding (2492 cells)
# Cluster 17 – CCB095, Cancer B20+aPDL1– Responding (1918 cells)
# Cluster 25 – CCB095, Cancer B20+Pi3Ki – Responding (836 cells)
# Cluster 12 – CCB104, T-cells/NK B20_+_PDL1_+_LTbR_RESPONDING
# Cluster 44 – CCB104, T-cells/NK B20_+_PDL1_+_LTbR_RESPONDING
# Cluster 13 – CCB104, B-cells B20_+_PDL1_+_LTbR_RESPONDING

c("Cancer", "Cancer", "Cancer", "Cancer", "Cancer", "Cancer", "Macrophages",
  "Cancer", "Cancer", "Acinar", "Cancer_B20_PDL1_CTLA4_Resp", "Endothelial",
  "T-cells/NK_B20_PDL1_LTbR_Resp", "B-Cell_B20_PDL1_LTbR_Resp", "Cancer",
  "Cancer_B20_PDL1_Resp1", "T-cells/NK", "Cancer_B20_PDL1_Resp2", "Cancer",
  "Cancer", "Cancer", "Cancer", "Monocytes", "Macrophages", "Cancer",
  "Cancer_B20_Pi3Ki_Resp", "Acinar", "T-cells/NK", "Fibroblasts", "Ductal/Stem",
  "Endothelial", "Cancer", "Cancer", "Cancer", "Cancer", "Macrophages",
  "Cancer", "B-Cell", "Endothelial", "Monocytes", "Cancer", "Smooth Muscle",
  "B-Cell", "Cancer", "T-cells/NK_B20_PDL1_LTbR_Resp") -> new.cluster.ids

## unfiltered
# c("Acinar", "Acinar", "Cancer", "Cancer", "Cancer",
# "Macrophages/Myeloids", "Cancer", "Cancer", "Endothelial cells", "Cancer",
# "Macrophages/Myeloids", "T cells", "Cancer", "Cancer", "Cancer", "B cells",
# "Cancer", "T cells", "Acinar", "Cancer", "Cancer", "Cancer",
# "Cancer", "Acinar", "Fibroblasts/Smooth Muscle", "Cancer", "Cancer",
# "Fibroblasts", "Monocytes", "Cancer", "Cancer", "B cells",
# "Neutrophils") -> new.cluster.ids
