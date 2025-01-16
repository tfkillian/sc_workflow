The files in this directory describe an example R Markdown file that demonstrates how to:

- Load the filtered feature-barcode matrix (i.e., barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz) from Cell Ranger output.
- Create a Seurat object.
- Perform quality control (QC) and filtering.
- Run PCA, UMAP, and clustering.
- Visualize results with Seurat’s plotting functions.

Notes/Assumptions

- This example assumes you have already run the Snakemake pipeline that generates Cell Ranger outputs in directories such as results/cellranger_count/{sample}/outs/filtered_feature_bc_matrix/.
- Adjust paths, filter cutoffs, and clustering parameters to suit your data and biological context.
- We use a standard Seurat pipeline (using NormalizeData(), FindVariableFeatures(), etc.). For some datasets, you may choose SCTransform or other normalization methods.
