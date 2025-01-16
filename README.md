# sc_workflow

This repository contains reproducible workflows for processing and analyzing 10x Genomics single-cell RNA-seq data. The pipelines include:

- **Snakemake Workflow**: Automates the end-to-end processing of raw FASTQ files through Cell Ranger, quality control (FastQC, MultiQC), and additional QC/doublet detection using Scrublet.
- **Seurat Workflow (RMarkdown)**: Loads Cell Ranger output data into a Seurat object, performs quality control (QC), filtering, normalization, dimensionality reduction (PCA/UMAP), clustering, and marker gene identification.
- **Scanpy Workflow (Jupyter Notebook)**: Implements a similar analysis pipeline in Python. It loads data from Cell Ranger output, performs QC, filtering, normalization, variable gene detection, PCA, UMAP, clustering, and optional marker gene analysis.

## Workflows

### 1. Snakemake Workflow

- **Purpose**:  
  Automates the processing of raw 10x FASTQ files through:
  - **Cell Ranger**: Aligns reads and generates filtered feature-barcode matrix output.
  - **FastQC & MultiQC**: Performs preliminary quality control on FASTQ files.
  - **Scrublet QC**: Runs an additional QC step using Scrublet for doublet detection and filtering.
  
- **Files**:
  - `Snakefile`: Defines the rules and steps to run the analysis.
  - `config.yaml`: Contains sample lists, reference genome paths, and resource definitions.
  - `scripts/scrublet_qc.py`: A helper script invoked by the Snakemake workflow for performing Scrublet-based QC.

- **Usage**:  
  Run the pipeline via command line. For example, with SLURM-based cluster execution:
  ```bash
  snakemake --jobs 50 --cluster-config cluster.json --cluster "sbatch --cpus-per-task={cluster.cores} --mem={cluster.mem} --time={cluster.time}" --latency-wait 60
  ```
  
### 2. Seurat Workflow (RMarkdown)

- Purpose: Provides an interactive analysis using R and Seurat. The workflow:
* Loads the Cell Ranger generated feature-barcode matrix files.
* Creates a Seurat object.
* Performs quality control (mitochondrial content, gene counts).
* Normalizes the data, identifies variable features, runs PCA, performs clustering, and visualizes results via UMAP.
* Optionally identifies marker genes for clusters.
* Files: R/single_cell_analysis.Rmd: An RMarkdown document that steps through the analysis.
* Usage: Open the RMarkdown file in RStudio or your preferred IDE and knit the document to generate an HTML report.

### 3. Scanpy Workflow (Jupyter Notebook)

- Purpose: Implements a similar analysis pipeline in Python using Scanpy. The notebook:
* Reads 10x Genomics filtered data into an AnnData object.
* Calculates QC metrics and filters out low-quality cells.
* Normalizes and log-transforms the data.
* Identifies highly variable genes and scales the data.
* Runs dimensionality reduction (PCA, UMAP) and clusters cells (using the Leiden algorithm).
* Provides visualizations for quality metrics and clustering outcomes.
* Files: Python/single_cell_analysis.ipynb: A Jupyter Notebook demonstrating the entire Scanpy-based workflow.
* Usage: Open the notebook with Jupyter Notebook or JupyterLab and run the cells sequentially to execute the analysis. The notebook includes code to save the processed AnnData object for further analysis.