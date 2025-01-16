############################# FAM - Find All Markers ###########################
## This script uses Seurats FindAllMarkers function to find the percentage of
## cells in one cluster (in the last column of output) expressing a gene 
## (in column A) (and give a measure of its strength with logFC) and a
## significance with p-value (adjusted with Bonferroni correction)
## compared to the percentage of remaining cells (so all the other ones)
## expressing this gene as well can help you with the annotation of clusters.
## use the MAST test (performs better) and ensure the default assay integrated!

## load libraries
library("R.utils")
library("dplyr")
library("Seurat")
library("ggplot2")
library("future")
library("cowplot")
## library("harmony")

## parallelize workflow
plan("multiprocess", workers = 40) # uses 50 CPU
options(future.globals.maxSize = 14000 * 1024^2) ## 5GB per worker

## evergreen variables
proj_name <- "pancreatic_mice" ## name of project
proj_version <- "filter"       ## name of version of project
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

############################# FAM - Find All Markers ###########################
# proj_name <- "cancer" ## name of project
# seu <- readRDS(file = paste0("cancer_sub_pancreatic_mice.rds"))

proj_name <- "macro" ## name of project
seu <- readRDS(file = paste0("macro_sub/FINAL_latest_ann_macro.rds"))
# seu <- readRDS(file = paste0("macro_sub/sub_macro.rds"))
# seu <- readRDS(file = paste0("sub_macro.rds"))

proj_name <- "tcell" ## name of project
# seu <- readRDS(file = paste0("tcell_sub/FINAL_latest_ann_tcell.rds"))
# seu <- readRDS(file = paste0("tcell_sub/sub_tcell.rds"))
seu <- readRDS(file = paste0("tcell_sub/FINAL_latest_ann_tcell.rds"))

proj_name <- "endo" ## name of project
# seu <- readRDS(file = paste0("endo_sub/sub_endo.rds"))
seu <- readRDS(file = paste0("endo_sub/FINAL_latest_ann_endo.rds"))

saveRDS(seu, file = "hf_rerun_complete_metadata3.rds")

#DefaultAssay(seu) <- def_assay
seu_markers <- FindAllMarkers(seu, test.use = "wilcox", min.pct = 0, ###########
                              logfc.threshold = 0, max.cells.per.ident = Inf,
                              latent.vars = NULL, min.cells.feature = 0,
                              min.cells.group = 0)

# seu_markers %>%
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
#                 gene != "Lars2",
#                 gene != "Sycn",
#                 gene != "Ctrl",
#                 gene != "Reg2") -> seu_markers

write.table(seu_markers, paste0("seu_markers_", proj_name, "_global",  ".txt"),
            col.names = NA, sep = "\t")

# file_path <- "~/Documents/tmp/pancreatic_melanie/"
# read.table(paste0(file_path, "seu_markers_cancer_markers6.txt"),
#            header = TRUE, sep = "\t") -> seu_markers

## heatmap on avg log FC
seu_markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 5, wt = avg_log2FC) -> topn_fc
pdf(paste0("Heatmap_logFC_", proj_name, ".pdf"), width = 30, height = 16)
DoHeatmap(object = seu, slot = "scale.data", features = topn_fc$gene, angle = 90) +
  #theme(axis.title.x = element_text(angle = 90, hjust = 1)) +
  NoLegend()
dev.off()
# 

# saveRDS(seu, file = "cancer_recluster_actual.rds")
# heatmap on p-val
topn_padj <- seu_markers %>% dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 5, wt = -p_val_adj)
pdf(paste0("Heatmap_pval_", proj_name, ".pdf"), width = 30, height = 20)
DoHeatmap(seu, slot = "scale.data", features = topn_padj$gene) +
  NoLegend()
dev.off()


### markers_40 <- FindMarkers(seu, ident.1 = 40) 

# data.frame(cluster_0 = topn_fc$gene[1:5],
#            cluster_1 = topn_fc$gene[6:10],
#            cluster_2 = topn_fc$gene[11:15],
#            cluster_3 = topn_fc$gene[16:20],
#            cluster_4 = topn_fc$gene[21:25],
#            cluster_5 = topn_fc$gene[26:30],
#            cluster_6 = topn_fc$gene[31:35],
#            cluster_7 = topn_fc$gene[36:40],
#            cluster_8 = topn_fc$gene[41:45],
#            cluster_9 = topn_fc$gene[46:50],
#            cluster_10 = topn_fc$gene[51:55],
#            cluster_11 = topn_fc$gene[56:60]#,
#            #cluster_12 = topn_fc$gene[61:65]
#            ) %>%
#            t() %>% as.data.frame() %>% 
#            tibble::rownames_to_column() -> tcell_df
# names(tcell_df) <- c("cluster", "top_gene1", "top_gene2", "top_gene3",
#                     "top_gene4", "top_gene5")

# write.table(tcell_df, paste0("tcell_cluster_markers.txt"),
#             col.names = NA, sep = "\t")
# write.table(endo_df, paste0("endo_cluster_markers_.txt"),
#             col.names = NA, sep = "\t")


#c("basophil", "macro_M1", "macro_M1", "macro_M1", "macro_M1",
  # "macro_M1_B20_Pi3Ki_resp", "neutrophils", "macro_M1", "MAST", "macro_M2",
  # "macro_M1", "monocytes", "DCs") -> new.cluster.ids

# c("T cells CD8+", "T cells CD8+", "T cells CD4+", "Tregs", "T cells CD8+",
#   "T cells CD8+ Early", "T cells CD8+", "NK cells", "T cells CD4+ Naive",
#   "T cells CD8+ Naive", "T cells CD8+ Ex", "Tfh") -> manual_labels

# c("artery Esm1+", "capillary Aplnr+", "capillary Cd36+", "capillary Chga+",
#   "artery Rgs5+", "vein Madcam1+", "tip cell Cxcl10+", "vein C1qb+",
#   "tip cell Gm42418+", "lymphatic", "vein Fn1+",
#   "artery CCB023 Cxx1a+") -> new.cluster.ids


########################### Differential expression - DE #######################

########################### DE - neutrophol clusters ###########################
seu0 <- readRDS(file = paste0("neutro1.rds"))
proj_name <- "neutro_no104"
seu <- subset(seu0, subset = orig.ident != c("CCB104"))
neu_resp <- subset(seu, subset = Response == c("Responding"))
neu_rel <- subset(seu, subset = Response == c("Relapsing"))
neu_resp %>% merge(y = c(neu_rel), project = proj_name) -> neu1
saveRDS(neu1, file = "neutro_no104.rds")
neu1 <- SetIdent(object = neu1, value = "Response")
neutro_comp <- FindMarkers(neu1, ident.1 = "Responding", ident.2 = "Relapsing")
##
neutro_comp %>%
  tibble::rownames_to_column() %>%
  dplyr::rename(gene = rowname) %>%
  dplyr::filter(gene != "Ins1",
                gene != "Ins2",
                gene != "Iapp",
                gene != "Cela2a",
                gene != "Cela3b",
                gene != "Cela1",
                gene != "Ppy",
                gene != "Try4",
                gene != "Try5",
                gene != "Prss2",
                gene != "Reg1",
                gene != "Cpb1",
                gene != "Pnlip",
                gene != "Clps",
                gene != "Ctrb1",
                gene != "Rnase1",
                gene != "2210010C04Rik",
                gene != "Nnat",
                gene != "Chga",
                gene != "Pcsk1n",
                gene != "Sycn") -> neutro_clean

write.table(neutro_clean, paste0("neutrophils_resp_vs_rel_", proj_name, ".txt"),
            col.names = NA, sep = "\t")

## heatmap on avg log FC
# neutro_clean %>%
#   #dplyr::group_by(Response) %>%
#   dplyr::slice(1:10) -> topn_fc
# pdf(paste0("Heatmap_logFC_", proj_name, ".pdf"), width = 30, height = 10)
# DoHeatmap(object = neu1, slot = "scale.data",
#           features = c("Cebpb", "Pfn1", "Malat1", "Tmsb4x", "Fos", "Fau",
#                        "Ypel3", "Ftl1", "Ube2d3", "Eif1")) +
#   #theme(axis.title.x = element_text(angle = 90, hjust = 1)) +
#   NoLegend()
# dev.off()

res_final1 <- neutro_clean[!is.na(neutro_clean$p_val_adj), ]
res_final1$threshold <- as.factor(abs(res_final1$avg_log2FC) > 1 &
                                      res_final1$p_val_adj < 0.05)

pdf(paste0("Vol_plot_neutrophils_resp_rel", proj_name, ".pdf"))
ggplot2::ggplot(data = res_final1,
                aes(x = avg_log2FC, y = -log10(p_val_adj), color = threshold)) +
                geom_text(label = res_final1$gene, nudge_x = 0.1,
                          nudge_y = 0.1, check_overlap = TRUE) +
                theme(legend.position = "none") +
                geom_point(alpha = 0.4, size = 0.5) +
                xlab("log2 fold change") + ylab("-log10 padj-value") +
                theme(plot.title = element_text(hjust = 0.5)) +
                scale_color_manual(values = c("#000000", "#FF0000")) -> p1
print(p1)
dev.off()

##################### DE - cancer/contamination clusters #######################
seu <- readRDS("cancer_contam.rds")
seu@meta.data$new_labels <- as.character(seu@meta.data$Cluster)
seu@meta.data$new_labels <- gsub("Cancer_B20_PDL1_CTLA4_Resp", "Cancer", seu@meta.data$new_labels)
seu@meta.data$new_labels <- gsub("Cancer_B20_PDL1_Resp1", "Cancer", seu@meta.data$new_labels)
seu@meta.data$new_labels <- gsub("Cancer_B20_PDL1_Resp2", "Cancer", seu@meta.data$new_labels)
seu@meta.data$new_labels <- gsub("Cancer_B20_Pi3Ki_Resp", "Cancer", seu@meta.data$new_labels)
seu@meta.data$new_labels <- gsub("_", " ", seu@meta.data$new_labels)
seu <- SetIdent(object = seu, value = 'new_labels')

mark_con1 <- FindMarkers(seu, ident.1 = "Cancer", ident.2 = "Contamination B-Cells")
mark_con2 <- FindMarkers(seu, ident.1 = "Cancer", ident.2 = "Contamination EC")
mark_con3 <- FindMarkers(seu, ident.1 = "Cancer", ident.2 = "Contamination Macro")
mark_con4 <- FindMarkers(seu, ident.1 = "Cancer", ident.2 = "Contamination T-Cells")

cell_type <- "cancer"
write.table(mark_con1, paste0("markers_bcells_", cell_type, "_", proj_name, ".txt"),
            col.names = NA, sep = "\t")
write.table(mark_con2, paste0("markers_endo_", cell_type, "_", proj_name, ".txt"),
            col.names = NA, sep = "\t")
write.table(mark_con3, paste0("markers_macro_", cell_type, "_", proj_name, ".txt"),
            col.names = NA, sep = "\t")
write.table(mark_con4, paste0("markers_tcells_", cell_type, "_", proj_name, ".txt"),
            col.names = NA, sep = "\t")

########################### DE - cancer clusters ###############################
## we will subset the cancer clusters to determine what are interesting markers
seu <- readRDS(file = paste0("FINAL_cancer_ann_", proj_name, ".rds"))
# c("0", "1", "2", "3", "4", "5", "7", "8", "10", "14", "15", "17", "18", "19",
#   "20", "21", "24", "25", "31", "32", "33", "34", "36", "40", "43")
# c("Cancer_B20_PDL1_CTLA4_Resp", "Cancer_B20_PDL1_Resp1",
#   "Cancer_B20_PDL1_Resp2",  "Cancer_B20_Pi3Ki_Resp")

# Cluster 10 – CCB095, B20+aPDL1+aCTLA4 – Responding (3182 cells)
can1 <- subset(seu, subset = Cluster == c("Cancer", "Cancer_B20_PDL1_CTLA4_Resp")) 
mark_can1 <- FindMarkers(can1, ident.1 = "Cancer_B20_PDL1_CTLA4_Resp")

## Cluster 15 – CCB030, B20+aPDL1– Responding (2492 cells)
can2 <- subset(seu, subset = Cluster == c("Cancer", "Cancer_B20_PDL1_Resp1"))
mark_can2 <- FindMarkers(can2, ident.1 = "Cancer_B20_PDL1_Resp1") 

## Cluster 17 – CCB095, B20+aPDL1– Responding (1918 cells)
can3 <- subset(seu, subset = Cluster == c("Cancer", "Cancer_B20_PDL1_Resp2"))
mark_can3 <- FindMarkers(can3, ident.1 = "Cancer_B20_PDL1_Resp2") 

## Cluster 25 – CCB095, B20+Pi3Ki – Responding (836 cells)
can4 <- subset(seu, subset = Cluster == c("Cancer", "Cancer_B20_Pi3Ki_Resp"))
mark_can4 <- FindMarkers(can4, ident.1 = "Cancer_B20_Pi3Ki_Resp")

## Cluster 15 vs 17
can5 <- subset(seu, subset = Cluster == c("Cancer_B20_PDL1_Resp1",
                                          "Cancer_B20_PDL1_Resp2"))
mark_can5 <- FindMarkers(can5, ident.1 = "Cancer_B20_PDL1_Resp1")

cell_type <- "cancer"
write.table(mark_can1, paste0("markers_B20_PDL1_CTLA4_Resp_", cell_type, "_",
                              proj_name, ".txt"), col.names = NA, sep = "\t")
write.table(mark_can2, paste0("markers_B20_PDL1_Resp1_", cell_type, "_",
                              proj_name, ".txt"), col.names = NA, sep = "\t")
write.table(mark_can3, paste0("markers_B20_PDL1_Resp2_", cell_type, "_",
                              proj_name, ".txt"), col.names = NA, sep = "\t")
write.table(mark_can4, paste0("markers_B20_Pi3Ki_Resp_", cell_type, "_",
                              proj_name, ".txt"), col.names = NA, sep = "\t")
write.table(mark_can5, paste0("markers_B20_PDL1_Resp_1v2_correct", cell_type, "_",
                              proj_name, ".txt"), col.names = NA, sep = "\t")
## read data
# Melanie.... I was looking at the DE genes in the TC clusters that you send to
# me and I realised that the txt files of the comparison between cluster 15vs17
# (last email) and B20+PDL1+Pi3Ki (cluster 25) are the same. This is normal? 
# setwd(paste0("/Users/u0112671/Documents/tmp/pancreatic_melanie/results/",
#              "raw1000dr/cancer_treatment_markers/marker_genes_cancer_clusters/"))
# read.table("./markers_B20_Pi3Ki_Resp_cancer_pancreatic_mice.txt",
#            header = TRUE) -> a1
# read.table("./markers_B20_PDL1_Resp_1v2_cancer_pancreatic_mice.txt",
#            header = TRUE) -> a2

## Expression Heatmap
# clust_avgs <- AverageExpression(seu, return.seurat = TRUE)
# pdf(paste0("Exp_Heatmap_5_", proj_name, ".pdf"))
# DoHeatmap(clust_avgs, features = unlist(TopFeatures(seu[["pca"]], balanced = TRUE)),
#           size = 5, draw.lines = FALSE)
# dev.off()
