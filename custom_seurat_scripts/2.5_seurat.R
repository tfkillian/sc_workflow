############################### sc-pipeline ####################################
## this script reads the pre-processed Seurat objects taken from the previous
## script `1_preprocess.R` and filters the objects further, normalizes the data
## and finds variable features. running PCA, UMAP and generating various
## diagnostic plots. Lastly, the clusters resulting from the UMAP are annotated 
## sample and response type, filters the cells, produces QC plots,

## load libraries
library("R.utils")
library("dplyr")
library("Seurat")
library("ggplot2")
library("future")
library("cowplot")

## parallelize workflow
plan("multiprocess", workers = 20) # uses 50 CPU
options(future.globals.maxSize = 14000 * 1024^2) ## 5GB per worker

## evergreen variables
proj_name <- "RIPTAG" ## name of project
proj_version <- "rerun_hf" ## name of version of project
server_path <- "/media/seq-srv-05/vrc/Project/Project_Theo/processed"
mito <- "^mt-" ## "^mt-" if mouse, "^MT-" if human
up_feat <- 6000 ## upper nFeature threshold
low_feat <- 200 ## lower nFeature threshold default 200 #####
var_feat <- 2000 ## variableFeature threshold
low_count <- 400 ## lower nCount threshold default 400 ######
p_mt <- 20 ## percent mito threshold
NPC <- 40 ## number of PCs of first PCA
k_para <- 50 ## number of k paramaters in FindNeighbors default is 30
sig_dims <- 12 ## number of PCs of final PCA
resolution <- 1.4 ## resolution
RESOLUTION <- 'RNA_snn_res.1.4'; ## resolution string

################# read Seurat object from "voltron" script ####################
# filtered_seurat <- readRDS(paste0(server_path , "merged_seurat_no_104_no_contam.rds"))
# seu <- readRDS("seurat_rerun.rds") ## with regressed variables

##############################################################################
## this chunk of code reads the processed Seurat files for each sample 
## read files into list (if running the saved files on the server)
file_list <- list.files(server_path, pattern = "processed_")
fls <- lapply(file_list, function(i) readRDS(file = i))
names(fls) <- gsub("processed_", "", gsub(".rds", "", unlist(file_list)))

## merge dataset, and add cell-ids and subset on parameters
fls[[1]] %>%
merge(y = fls[2:length(fls)], project = proj_name, add.cell.ids = names(fls)) %>%
  subset(subset = nFeature_RNA < up_feat &
                  nFeature_RNA > low_feat &
                  nCount_RNA > low_count &
                  percent.mt < p_mt) -> merged_seu

saveRDS(merged_seu, file = "filtered_merged_seu.rds")

## remove doublets
#filtered_seurat <- subset(merged_seu, subset = pANNPredictions == c("Singlet"))
filtered_seurat <- merged_seu

## run pipeline
filtered_seurat %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = var_feat) %>%
  # ScaleData(vars.to.regress = c("nCount_RNA", "percent.mt", "orig.ident")) %>% 
  ScaleData() %>% 
  RunPCA(npcs = sig_dims, verbose = FALSE) %>%
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
pdf(file = paste0("umap_", proj_name, "_ann.pdf"), width = 12, height = 8)
p1 <- DimPlot(seu, reduction = "umap", group.by = "orig.ident",
              pt.size = 0.1, label = FALSE, label.size = 2.0)
p2 <- DimPlot(seu, reduction = "umap", group.by = "Response",
              pt.size = 0.1, label = FALSE, label.size = 3.5,  repel = TRUE)
p3 <- DimPlot(seu, reduction = "umap", group.by = "Treatment",
              pt.size = 0.1, label = FALSE, label.size = 3.5,  repel = TRUE)
p4 <- DimPlot(seu, reduction = "umap",
              pt.size = 0.1, label = TRUE, label.size = 3.5,  repel = TRUE)
plot_grid(p1, p2, p3, p4, ncol = 2)
dev.off()

pdf(file = paste0("umap_doublets_", proj_name, ".pdf"), width = 8, height = 8)
DimPlot(seu, reduction = "umap", group.by = "pANNPredictions", pt.size = 0.1,
        label = FALSE, label.size = 2.0)
dev.off()

########################## tables for QC plots ################################
## table of number of cells 
write.table(seu@meta.data, paste0("cells_metadata_", proj_name, ".txt"),
            col.names = NA, sep = "\t")

## how many cells in each cluster per sample
cell_perC_perO <- table(seu@meta.data$ctive.ident, seu@meta.data$orig.ident)
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
# f <- ordered(features$RNA_snn_res.0.8, levels = c(0:clust_num))
# f <- ordered(features$RNA_snn_res.1.4,
#              levels = c(0:(max(as.integer(seu@meta.data$RNA_snn_res.1.4)) - 1)))
# 
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

############################### checkpont ######################################
seu <- StashIdent(seu, save.name = "seurat_cluster_res_1.4")
# saveRDS(seu, file = "seurat_rerun.rds")

################################ SingleR #######################################
## SingleR is Unsupervised clustering, which affords differential expression
## between 2 condition in a cluster - for each cluster. we usually use marker
## genes to annotate our clusters, but what if the marker genes are unknown,
## or an automated annotation is desired to compare for confirmation, SingleR is
## used. Basically, there are 2 type of annotations: cluster annotation or cell
## annotation
library("SingleR")
library("celldex")

# mouse = MouseRNAseqData()
mouse <- celldex::MouseRNAseqData()

## MOUSE
seuratObj <- seu
singlerObj <- as.SingleCellExperiment(seuratObj)
common <- intersect(rownames(singlerObj), rownames(mouse))
mouse <- mouse[common,]
singlerObj <- singlerObj[common,]
singler <- SingleR(test = singlerObj, ref = mouse, labels = mouse$label.main,
                   method = "cluster", ## change the resolution you chose
                   clusters = seuratObj@meta.data[,"RNA_snn_res.1.4"],
                   genes = "de") 

## save singlr results 
saveRDS(singler, file = paste0("singler_object_", proj_name, ".rds"))
write.table(singler, file = paste0( "singler_", proj_name, ".txt"),
            col.names = NA, sep = "\t")

## save checkpoint
# seu <- StashIdent(seu, save.name = "seurat_cluster_res_1.4")
# saveRDS(seu, file = "seurat_rerun.rds")

############# Assigning cell type identity to clusters (annotation) ############
#seu <- readRDS(file = paste0("seurat_recluster_no_104_no_contam.rds"))
seu <- SetIdent(object = seu, value = RESOLUTION)

## the cluster names are occasional put out of order, especially if you save and
## re-load a Seurat object. this code chunk ensures that misordering does not
## occur
cluster_ids <- singler@rownames
singler_labels <- singler@listData$labels
as.data.frame(cbind(cluster_ids, singler_labels)) %>%
  dplyr::mutate(real_cluster_ids = as.integer(as.character(cluster_ids))) %>% 
  dplyr::arrange(real_cluster_ids) -> singler_annotations

## SingleR annotations
c("Cancer", "Cancer", "Cancer", "Cancer", "Cancer",
  "Cancer", "Cancer", "Cancer", "Cancer", "Acinar",
  "Cancer", "Macrophages", "Cancer", "Acinar", "Cancer",
  "Endothelial cells", "T cells", "Cancer", "Cancer", "Macrophages",
  "Monocytes", "Macrophages", "Endothelial cells", "Acinar", "Fibroblasts",
  "Cancer", "Fibroblasts", "Cancer", "T cells", "Fibroblasts"
  ) -> manual_labels

cluster_ids <- as.character(0:(max(as.integer(seu@meta.data$RNA_snn_res.1.4)) - 1))
singler_annotations %>%
  dplyr::full_join(as.data.frame(cbind(cluster_ids, manual_labels)),
                   by = "cluster_ids") -> new_cluster_ids
new_cluster_ids %>% 
  write.table(file = paste0("cluster_ids_", proj_name, ".txt"), col.names = NA,
              sep = "\t")
# new.cluster.ids <- new_cluster_ids$manual_labels

## recluster @RES = 1.4 primary annotations
c("Cancer", "Cancer", "Cancer", "Cancer", "Cancer",
  "Cancer", "Cancer", "Cancer", "Cancer", "Acinar",
  "Cancer", "Myeloids", "Cancer", "Acinar", "Cancer",
  "Endothelial cells", "T/NK cells", "Cancer", "Cancer", "Myeloids",
  "Myeloids", "Myeloids", "Endothelial cells", "Acinar", "Ductal/Stem",
  "Cancer", "Fibroblasts", "Cancer", "T/NK cells", "Smooth Muscle"
  ) -> new.cluster.ids

## add column in metadata with the "Cluster" title and the cluster name
names(new.cluster.ids) <- levels(seu) ## "seurat_cluster_res_1.4"
seu <- RenameIdents(seu, new.cluster.ids)
# seu <- StashIdent(seu, save.name = "final_primary_annotations")
seu <- StashIdent(seu, save.name = "rerun_primary_annotations")
head(seu@meta.data) ## verify new column

seu <- SetIdent(object = seu, value = "rerun_primary_annotations")

## make annotated UMAPs
pdf(paste0("Annotated_primary_", proj_name, ".pdf"), width = 11, height = 8.5)
DimPlot(seu, reduction = "umap", label = FALSE, pt.size = 0.1) +
  labs(title = "Annotated UMAP") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf(paste0("Annotated_labels_primary_", proj_name, ".pdf"), width = 11, height = 8.5) 
DimPlot(seu, reduction = "umap", label = TRUE, pt.size = 0.1, repel = TRUE) +
  labs(title = "Annotated UMAP") + theme(plot.title = element_text(hjust = 0.5))
dev.off() # the same with labels

## recluster @RES = 1.4 secondary annotations
# c("Cancer", "Cancer", "Cancer", "Cancer", "Cancer", "Cancer", "Cancer",
#   "Acinar", "Cancer", "Cancer", "Cancer B20 aPDL1 CTLA4 Resp", "Macrophages",
#   "Cancer B20 aPDL1 Resp", "Cancer", "Endothelial cells", "Cancer",
#   "Cancer B20 PiK3i", "Cancer", "T-cells", "Cancer", "Cancer B20 aPDL1 Resp",
#   "Macrophages", "Cancer B20 aPDL1 Rel", "Cancer B20 Pik3i Rel", "Monocytes",
#   "Cancer B20 PiK3i Resp", "Acinar", "Fibroblasts", "Ductal/Stem", "Cancer",
#   "Cancer", "Endothelial cells B20 PiK3i Resp", "Smooth Muscle",
#   "Endothelial cells", "Cancer B20 aPDL1 CTLA4 Resp", "Cancer B20 PiK3i"
#   ) -> new.cluster.ids
# 
## add column in metadata with the "Cluster" title and the cluster name
seu <- SetIdent(object = seu, value = "seurat_cluster_res_1.4")
head(seu@meta.data) ## verify new column
names(new.cluster.ids) <- levels(seu)
seu <- RenameIdents(seu, new.cluster.ids)
seu <- StashIdent(seu, save.name = "final_subtype_annotations")
head(seu@meta.data) ## verify new column

seu <- SetIdent(object = seu, value = "old_merged_primary_annotations")

## make annotated UMAPs
pdf(paste0("Annotated_secondary_", proj_name, ".pdf"), width = 11, height = 8.5)
DimPlot(seu, reduction = "umap", label = FALSE, pt.size = 0.1) +
  labs(title = "Annotated UMAP") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf(paste0("Annotated_labels_secondary_", proj_name, ".pdf"), width = 11, height = 8.5)
DimPlot(seu, reduction = "umap", label = TRUE, pt.size = 0.1, repel = TRUE) +
  labs(title = "Annotated UMAP") + theme(plot.title = element_text(hjust = 0.5))
dev.off() # the same with labels

### save tables again

## table of number of cells 
write.table(seu@meta.data,
            paste0("cells_metadata_", proj_name, "_final_secondary_annotations.txt"),
            col.names = NA, sep = "\t")

## how many cells in each cluster per sample
cell_perC_perO <- table(seu@meta.data$final_secondary_annotations,
                        seu@meta.data$orig.ident)
write.table(cell_perC_perO,
            paste0("cell_per_c_per_O_", proj_name,  "_final_secondary_annotations.txt"),
            col.names = NA, sep = "\t")

## how many cells in each cluster per treatment
cell_perC_perT <- table(seu@meta.data$final_secondary_annotations,
                        seu@meta.data$Treatment)
write.table(cell_perC_perT,
            paste0("cell_per_c_per_T_", proj_name,  "_final_secondary_annotations.txt"),
            col.names = NA, sep = "\t")

## how many cells in each cluster per response
cell_perC_perR <- table(seu@meta.data$final_secondary_annotations,
                        seu@meta.data$Response)
write.table(cell_perC_perR,
            paste0("cell_per_c_per_R_", proj_name, "_final_secondary_annotations.txt"),
            col.names = NA, sep = "\t")

## save RDS file
saveRDS(seu, file = paste0("final_annotations_recluster.rds"))
##

# sub_e <- subset(seu, cells = colnames(seu)[seu$subtype_annotations_v3 %in% c(
#   "Lymphatic", "Capillaries Oaz2+", "Capillaries Cd36+", "Pericytes",
#   "Tip cells", "Venous", "Capillaries arterial", "Arterial",
#   "Capillaries interferon")])


# sub_m <- subset(seu, cells = colnames(seu)[seu$subtype_annotations_v3 %in% c(
#   "cDCs", "Mast cells/pDCs", "Neutrophils", "F13a1+ M2 Macrophages",
#   "M2 Macrophages")])
# 
# sub_t <- subset(seu, cells = colnames(seu)[seu$subtype_annotations_v3 %in% c(
#   "NK cells", "Cd8+ Exhausted", "Cd4+ Th2+", "Proliferating T-cells",
#   "Unknown T-cells Th17+", "Unknown T-cells", "Gamma Delta T-cells",
#   "Cd8+ Memory", "Tregs")])

# meta_data <- seu@meta.data
# 
# meta_data$is_contamination <- meta_data$subtype_annotations_v3
# meta_data %>%
#   dplyr::mutate(is_contamination = dplyr::case_when(
#     "Unknown Cancer" ~ TRUE,
#     "Unknown Myeloids" ~ TRUE,             
#     "Unknown T-cells" ~ TRUE,                   
#     "Unknown B-cell" ~ TRUE,
#     "Unknown Endothelial cells" ~ TRUE,
#     "Unknown Acinar" ~ TRUE
#     )) -> meta_data
