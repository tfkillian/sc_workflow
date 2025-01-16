############################### pre-processing #################################
## this script reads the 10X data files into merged Seurat objects for each
## sample and response type, filters the cells, produces QC plots, normalizes
## the data and finds variable features.

## /media/seq-srv-05/vrc/Project/Project_Gino/software/R/R-3.6.3
## /data/projects/Project_Theo/R/R-4.0.3

## load libraries
library("R.utils")
library("dplyr")
library("Seurat")
library("ggplot2")
library("future")
library("cowplot")
library("DoubletFinder")

## parallelize workflow
plan("multiprocess", workers = 24) # uses 50 CPU
options(future.globals.maxSize = 14000 * 1024^2) ## 5GB per worker

## evergreen variables
proj_name <- "pancreatic_mice"
proj_version <- "filter"
server_path <- "/media/seq-srv-05/vrc/Project/Project_Theo/"
local_file_path <- "~/Documents/tmp/202008_guyot/pancreatic_melanie/"  ## path to analysis output
list_path <- paste0(local_file_path, "sample_data/")      ## path to metadata file
data_path <- paste0(local_file_path, "sc_data/")          ## path to raw data 
img_path <- paste0(local_file_path, "plots/")             ## path to image output
res_path <- paste0(local_file_path)           ## path to results output
mito <- "^mt-" ## "^mt-" if mouse, "^MT-" if human
up_feat <- 6000; low_feat <- 200; var_feat <- 2000; percent_mito <- 30;
low_count <- 400 ## lower nCount threshold default 400
NPC <- 35; sig_dims <- 15; resolution <- 0.2; RESOLUTION <- 'RNA_snn_res.0.2';
ggl_a <- list(); ggl_b <- list(); ggl_c <- list(); ggl_d <- list(); ggl_e <- list();

## csv file with the sample metadata
read.csv(paste0(server_path, "sample_df.csv"), sep = ",") %>%
  # dplyr::filter(treatment == "B20") %>% # done
  # dplyr::filter(treatment == "B20_+_aPDL1") %>% # done
  # dplyr::filter(treatment == "B20_+_Pi3Ki") %>% # done
  # dplyr::filter(treatment == "Untreated") %>% # done
  # dplyr::filter(treatment == "B20_+_aPD1_+_aCTLA4") %>% # done
  # dplyr::filter(treatment == "B20_+_aPDL1_+_Pi3Ki") %>% # done
  # dplyr::filter(treatment == "B20_+_PDL1_+_LTbR") %>% # done
  # dplyr::filter(sample_name == "CCB090") %>%
  #as.data.frame() -> samp_df0
  as.data.frame() -> samp_df

# for (j in 1:nrow(samp_df0)) {
# samp_df0 %>% dplyr::filter(sample_number == j) -> samp_df

## determine number of samples to loop through, a vector which corresponds to
## the number of rows in the csv file
samp_num <- nrow(samp_df)

## initialize empty lists of 10X and Seurat files to be populated below:
fl10x <- list() ## list of 10X objects
fls <- list() ## list of seurat objects
#data_list <- as.list(samp_df$local_file_path)
# data_list <- as.list(samp_df$vsc_file_path)
data_list <- as.list(samp_df$filter_file_path)

## for loop reads all barcode, matrix and metadata files as 10X samples and puts
## those objects in a list
fl10x <- lapply(data_list, function(i) Seurat::Read10X(data.dir = i))
names(fl10x) <- samp_df$sample_name

## for loop reads all 10X samples as Seurat objects and puts those objects in a
## list. the name of each Seurat object is the name of that sample 
fls <- lapply(fl10x, function(i) CreateSeuratObject(counts = i,
                                                    min.cells = 1,
                                                    min.features = 1))

## assign experiment conditions, orig.ident and idents 
for (i in 1:length(fls)) {
fls[[i]]@meta.data$Response <- as.factor(samp_df$response[i])
fls[[i]]@meta.data$Treatment <- as.factor(samp_df$treatment[i])
fls[[i]]@meta.data$orig.ident <- as.factor(samp_df$sample_name[i])
Idents(x = fls[[i]]) <- as.factor(paste0(samp_df$sample_name[i]))
}

### ensure Seurat metadata is factor and not character, calculate percent.mt
for (i in 1:length(fls)) {
fls[[i]]@meta.data[, 1] <- as.factor(fls[[i]]@meta.data[, 1])
fls[[i]]@meta.data[, 4] <- as.factor(fls[[i]]@meta.data[, 4])
fls[[i]]@meta.data[, 5] <- as.factor(fls[[i]]@meta.data[, 5])

### percent mitochondrial NOTE: # ! mouse: lower case, human: CAPITALS: "^MT-"
PercentageFeatureSet(fls[[i]], pattern = mito) -> fls[[i]][["percent.mt"]]
}

## pre-filtered  violin plot and feature scatter plot
# for (i in 1:length(fls)) {
# VlnPlot(fls[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
#         ncol = 3, pt.size = 0.2, fill.by = "orig.ident") -> ggl_a[[i]]
# 
# CombinePlots(plots = list(
#   FeatureScatter(fls[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt"),
#   FeatureScatter(fls[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#   )) -> ggl_b[[i]]
# }

## filter samples
for (i in 1:length(fls)) {
# fls[[i]]@meta.data %>%
#     group_by(orig.ident) %>%
#     filter(nFeature_RNA > low_feat,
#            percent.mt > -Inf,
#            percent.mt < percent_mito,
#            nCount_RNA > low_count,
#            nCount_RNA < +Inf) %>%
#     summarise(n = n()) %>%
# write.table(file = paste0(res_path, "filt_", names(fls)[i], ".txt"), sep = "\t")
  
subset(fls[[i]],
       subset = nFeature_RNA > low_feat &
                percent.mt < percent_mito &
                nCount_RNA > low_count) -> fls[[i]]
}

## post-filter violin plot and feature scatter plot
# for (i in 1:length(fls)) {
# VlnPlot(fls[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
#         ncol = 3, pt.size = 0.2) -> ggl_c[[i]]
# 
# CombinePlots(plots = list(
#   FeatureScatter(fls[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt"),
#   FeatureScatter(fls[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#   )) -> ggl_d[[i]]
# }

## Normalizing the data and find variable features
for (i in 1:length(fls)) {
NormalizeData(fls[[i]], normalization.method = "LogNormalize",
              scale.factor = 10000) -> fls[[i]]
FindVariableFeatures(fls[[i]], selection.method = "vst",
                     nfeatures = var_feat) -> fls[[i]]
}

## save pre-filter and post-filter QC plots
# for (i in 1:length(fls)) {
# png(filename = paste0(img_path, "pre_filt_vlnplt_qc_", names(fls)[i], ".png"),
#     width = 1500, height = 700)
# print(ggl_a[[i]])
# dev.off()
# png(filename = paste0(img_path, "pre_filt_vlnplt_", names(fls)[i], ".png"),
#     width = 1500, height = 700)
# print(ggl_b[[i]])
# dev.off()
# png(filename = paste0(img_path, "filt_vlnplt_qc_", names(fls)[i], ".png"),
#     width = 1500, height = 700)
# print(ggl_c[[i]])
# dev.off()
# png(filename = paste0(img_path, "filt_vlnplt", names(fls)[i], ".png"),
#     width = 1500, height = 700)
# print(ggl_d[[i]])
# dev.off()
# }

## save files as rds
# for (i in 1:length(fls)) {
# saveRDS(fls[[i]], file = paste0(res_path, names(fls)[i], ".rds"))
# }

## read files into list (if running the saved files on the server)
# file_list <- list.files(server_path, pattern = "103|104")
# file_list <- list.files(server_path)
# fls <- lapply(file_list, function(i) readRDS(file = i))
# names(fls) <- gsub(".rds", "", unlist(file_list))

## scale data, run PCA and UMAP
for (i in 1:length(fls)) {
fls[[i]] %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = NPC, verbose = FALSE) %>% 
  RunUMAP(reduction = "pca", dims = 1:NPC) %>% 
  FindNeighbors(reduction = "pca", dims = 1:NPC, k.param = 30,
                compute.SNN = TRUE, nn.method = "rann", nn.eps = 0) %>%
  FindClusters(resolution = resolution) -> fls[[i]]
}

############################# doublet finder ###################################
## This detects doublet clusters (a doublet is created when  multiple cells are
## captured in a single gel bead). Artificial doublet cells are generated for
## each sample and Euclidian distances between each real cell and artificial
## cells is calculated. If a real cell is close to artificial doublets, it is
## considered as a doublet the goal is to end up with a table like this
# RNA_snn_res.0.8 | mean_pANN_score | Singlet | pct_Singlet | Doublet | pct_Doublet | TOTAL
# 0               | 0.27            | 5,873   | 92.0%       | 512     | 8.0%        | 6,385
# ...
## we observe the mean pann score (there is no precise threshold above which you
## consider a cluster to be a doublet cluster, but the more high it is, the more
## chance it is a doubet cluster) and the percentage of doublet cells. you can
## also have a plot that shows the doublet cells - depending on the colors, you
## can also suspect a doublet cluster. this step is performed a after picking
## the resolution, so extra time is not spent trying to ID a cluster with marker
## genes and singleR that clearly is a doublet cluster. Doublet finder is not a
## package but a function. Don't adjust parameters in the function definition,
## only the arguments after
# doubletFinder_v3 <- function(seu, PCs, pN = 0.25, pK, nExp, reuse.pANN = FALSE,
#                              sct = FALSE) {
#   require("Seurat"); require("fields"); require("KernSmooth")
#   ## Generate new list of doublet classificatons from existing pANN vector
#   if (reuse.pANN != FALSE ) {
#     pANN.old <- seu@meta.data[ , reuse.pANN]
#     classifications <- rep("Singlet", length(pANN.old))
#     classifications[order(pANN.old, decreasing = TRUE)[1:nExp]] <- "Doublet"
#     seu@meta.data[, paste("DF.classifications", pN, pK, nExp, sep = "_")] <- classifications
#     return(seu)
#   }
#   if (reuse.pANN == FALSE) {
#     ## Make merged real-artifical data
#     real.cells <- rownames(seu@meta.data)
#     data <- seu@assays$RNA@counts[, real.cells]
#     n_real.cells <- length(real.cells)
#     n_doublets <- round(n_real.cells/(1 - pN) - n_real.cells)
#     print(paste("Creating", n_doublets, "artificial doublets...", sep = " "))
#     real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
#     real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
#     doublets <- (data[, real.cells1] + data[, real.cells2])/2
#     colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
#     data_wdoublets <- cbind(data, doublets)
# 
#     ## Store important pre-processing information
#     orig.commands <- seu@commands
# 
#     ## Pre-process Seurat object
#     if (sct == FALSE) {
#       print("Creating Seurat object...")
#       seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
#       print("Normalizing Seurat object...")
#       NormalizeData(seu_wdoublets,
#                     normalization.method = orig.commands$NormalizeData.RNA@params$normalization.method,
#                     scale.factor = orig.commands$NormalizeData.RNA@params$scale.factor,
#                     margin = orig.commands$NormalizeData.RNA@params$margin
#                     ) -> seu_wdoublets
#       print("Finding variable genes...")
#       FindVariableFeatures(seu_wdoublets,
#                            selection.method = orig.commands$FindVariableFeatures.RNA$selection.method,
#                            loess.span = orig.commands$FindVariableFeatures.RNA$loess.span,
#                            clip.max = orig.commands$FindVariableFeatures.RNA$clip.max,
#                            mean.function = orig.commands$FindVariableFeatures.RNA$mean.function,
#                            dispersion.function = orig.commands$FindVariableFeatures.RNA$dispersion.function,
#                            num.bin = orig.commands$FindVariableFeatures.RNA$num.bin,
#                            binning.method = orig.commands$FindVariableFeatures.RNA$binning.method,
#                            nfeatures = orig.commands$FindVariableFeatures.RNA$nfeatures,
#                            mean.cutoff = orig.commands$FindVariableFeatures.RNA$mean.cutoff,
#                            dispersion.cutoff = orig.commands$FindVariableFeatures.RNA$dispersion.cutoff
#                            ) -> seu_wdoublets
#       print("Scaling data...")
#       ScaleData(seu_wdoublets,
#                 features = orig.commands$ScaleData.RNA$features,
#                 model.use = orig.commands$ScaleData.RNA$model.use,
#                 do.scale = orig.commands$ScaleData.RNA$do.scale,
#                 do.center = orig.commands$ScaleData.RNA$do.center,
#                 scale.max = orig.commands$ScaleData.RNA$scale.max,
#                 block.size = orig.commands$ScaleData.RNA$block.size,
#                 min.cells.to.block = orig.commands$ScaleData.RNA$min.cells.to.block
#                 ) -> seu_wdoublets
#       print("Running PCA...")
#       RunPCA(seu_wdoublets,
#              features = orig.commands$ScaleData.RNA$features,
#              npcs = length(PCs),
#              rev.pca =  orig.commands$RunPCA.RNA$rev.pca,
#              weight.by.var = orig.commands$RunPCA.RNA$weight.by.var,
#              verbose = FALSE) -> seu_wdoublets
#       pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[, PCs]
#       cell.names <- rownames(seu_wdoublets@meta.data)
#       nCells <- length(cell.names)
#       rm(seu_wdoublets); gc() # Free up memory
#     }
#     if (sct == TRUE) {
#       require("sctransform")
#       print("Creating Seurat object...")
#       seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
# 
#       print("Running SCTransform...")
#       seu_wdoublets <- SCTransform(seu_wdoublets)
# 
#       print("Running PCA...")
#       seu_wdoublets <- RunPCA(seu_wdoublets, npcs = length(PCs))
#       pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[ , PCs]
#       cell.names <- rownames(seu_wdoublets@meta.data)
#       nCells <- length(cell.names)
#       rm(seu_wdoublets); gc()
#     }
#     ## Compute PC distance matrix
#     print("Calculating PC distance matrix...")
#     dist.mat <- fields::rdist(pca.coord)
# 
#     ## Compute pANN
#     print("Computing pANN...")
#     pANN <- as.data.frame(matrix(0L, nrow = n_real.cells, ncol = 1))
#     rownames(pANN) <- real.cells
#     colnames(pANN) <- "pANN"
#     k <- round(nCells * pK)
#     for (i in 1:n_real.cells) {
#       neighbors <- order(dist.mat[, i])
#       neighbors <- neighbors[2:(k + 1)]
#       neighbor.names <- rownames(dist.mat)[neighbors]
#       pANN$pANN[i] <- length(which(neighbors > n_real.cells)) / k
#     }
#     print("Classifying doublets..")
#     classifications <- rep("Singlet",n_real.cells)
#     classifications[order(pANN$pANN[1:n_real.cells], decreasing = TRUE)[1:nExp]] <- "Doublet"
#     seu@meta.data[, paste("pANN", pN, pK, nExp, sep = "_")] <- pANN[rownames(seu@meta.data), 1]
#     seu@meta.data[, paste("DF.classifications", pN, pK, nExp, sep = "_")] <- classifications
#     return(seu)
#   }
# }

## format seurat metadata
for (i in 1:length(fls)) {
fls[[i]] <- SetIdent(object = fls[[i]], value = "orig.ident")
fls[[i]]$pANN <- "NA"
fls[[i]]$pANNPredictions <- "NA"
}

## for each sample, it will run doublet finder basically
# for (i in 1:length(fls)) {
#  # sample.cluster <- subset(fls[[i]], idents = sample)
#  # print(paste0("sample:", sample))
#  # length(rownames(fls[[i]]@meta.data))
#  ## ceiling value is a default - can be adapted for specific projects
#  expected.doublets <- ceiling(0.039 * length(rownames(fls[[i]]@meta.data)))
#  ## apply doubletFinder function
#  doubletFinder_v3(fls[[i]],
#                   PCs = 1:20, # number of significant PCs
#                   pN = 0.25, # defines number of artificial doublets
#                   pK = 0.01, # needs to be estimated
#                   nExp = expected.doublets,
#                   reuse.pANN = FALSE, sct = TRUE) -> fls[[i]]
#  fls[[i]]@meta.data[colnames(fls[[i]]),
#      paste("pANN_0.25_0.01", expected.doublets, sep = "_")] -> fls[[i]]$pANN
#  fls[[i]]@meta.data[colnames(fls[[i]]),
#      paste("DF.classifications_0.25_0.01",
#            expected.doublets, sep = "_")] -> fls[[i]]$pANNPredictions
#  fls[[i]]$pANN[colnames(fls[[i]])] <- fls[[i]]$pANN[colnames(fls[[i]])]
#  fls[[i]]$pANNPredictions[colnames(fls[[i]])] <- fls[[i]]$pANNPredictions[colnames(fls[[i]])]
#  # fls[[i]] <- NULL
# }

######
for (i in 1:length(fls)) {
## pK Identification (no ground-truth) -----------------------------------------
sweep_res_list <- DoubletFinder::paramSweep_v3(fls[[i]], PCs = 1:10, sct = FALSE)
sweep_stats <- DoubletFinder::summarizeSweep(sweep_res_list, GT = FALSE)
# bcmvn_95 <- find.pK(sweep.stats_95)

## apply doubletFinder function
expected_doublets <- ceiling(0.039 * length(rownames(fls[[i]]@meta.data)))
DoubletFinder::doubletFinder_v3(fls[[i]],
                 PCs = 1:10, # number of significant PCs
                 pN = 0.25, # defines number of artificial doublets
                 pK = max(as.numeric(levels(sweep_stats$pK))),
                 nExp = expected_doublets,
                 reuse.pANN = FALSE,
                 sct = FALSE) -> fls[[i]]

## add the doublet annotations to the metadata
fls[[i]]@meta.data[colnames(fls[[i]]),
     paste("pANN_0.25", max(as.numeric(levels(sweep_stats$pK))),
           expected_doublets, sep = "_")] -> fls[[i]]$pANN
fls[[i]]@meta.data[colnames(fls[[i]]),
     paste("DF.classifications_0.25", max(as.numeric(levels(sweep_stats$pK))),
           expected_doublets, sep = "_")] -> fls[[i]]$pANNPredictions
fls[[i]]$pANN[colnames(fls[[i]])] <- fls[[i]]$pANN[colnames(fls[[i]])]
fls[[i]]$pANNPredictions[colnames(fls[[i]])] <- fls[[i]]$pANNPredictions[colnames(fls[[i]])]
}

## write metadata
## check the mean PANN per cluster --> if high: doublet cluster, make a table
## with the number of doublet / cluster, the number of singlets. From there,
## calculate the % of doublets and singlets then calculate an average of "pANN"
# for (i in 1:length(fls)) {
# fls[[i]]@meta.data %>%
#   write.table(file = paste0("DF_metadata_", names(fls[i]), ".txt"),
#               col.names = NA, sep = "\t")
# }

## make a plot - 1 colour for doublets and 1 for singlets
# for (i in 1:length(fls)) {
# DimPlot(fls[[i]], pt.size = 0.1, label = FALSE, label.size = 0,
#         reduction = "umap", group.by = "pANNPredictions") +
#         theme(aspect.ratio = 1) -> ggl_e[[i]]
# }

## save doublet plots
# for (i in 1:length(fls)) {
# pdf(paste0("doublets_", names(fls[i]), ".pdf"))
# print(ggl_e[[i]])
# dev.off()
# }

## generate tables
#library("plyr")
## plyr needs to be replaced by dplyr, probably like this:
## fls[[i]]@meta.data %>%
##  group_by(RNA_snn_res.0.2) %>%
##  summarize(mean = mean(as.numeric(pANN)))

# for (i in 1:length(fls)) {
# c <- ddply(fls[[i]]@meta.data, ~ RNA_snn_res.0.2, summarise,
#            mean = mean(as.numeric(pANN)))
# y <- table(fls[[i]]@active.ident,
#            fls[[i]]@meta.data$pANNPredictions)
# y <- as.data.frame.matrix(y)
# w <- cbind(c, Doublet = y$Doublet)
# w <- cbind(w, Singlet = y$Singlet)
# w$TOTAL <- w$Doublet + w$Singlet
# w$pct_Doublet <- w$Doublet / w$TOTAL
# w$pct_Singlet <- w$Singlet / w$TOTAL
# w <- w %>% dplyr::rename(mean_pANN_score = mean)
# w <- w[, c(1, 2, 4, 7, 3, 6, 5)]
# w %>%
#   write.table(file = paste0("pann_per_C_Doublet", names(fls[i]), ".txt"),
#               col.names = NA, sep = "\t")
# }

## save processed files
for (i in 1:length(fls)) {
fls[[i]] %>%
  saveRDS(file = paste0("processed_", names(fls[i]), ".rds"))
}

#}
#detach("package:plyr", unload = TRUE)