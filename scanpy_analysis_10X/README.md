The files in this directory describe an example Python + Scanpy workflow in a Jupyter notebook format
for single-cell RNA-seq analysis. It mirrors the logic of the previous R + Seurat pipeline but uses Scanpy to:

- Load 10x Genomics Cell Ranger outputs (the filtered matrix files).
- Perform QC (filter low-quality cells, remove cells with excessive mitochondrial gene expression).
- Normalize and log-transform the data.
- Identify highly variable genes.
- Scale the data and run PCA, UMAP, and clustering.
- Visualize key QC metrics and clustering results.

Notes / Assumptions

- Your Cell Ranger–filtered matrix directory contains barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz, typically in a folder like:
- results/cellranger_count/SAMPLE_1/outs/filtered_feature_bc_matrix/
- Adjust file paths, filtering thresholds, and parameter settings to your data and biological needs.
