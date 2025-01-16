This directory contains an extended example of a Snakemake workflow that processes 10x Genomics
single-cell RNA-seq data—from raw FASTQs, through Cell Ranger alignment and basic QC, to an
additional QC/doublet detection step with Scrublet. This version also shows how you might set
up rules with resource declarations so that you can run the workflow in cluster execution mode.

Assumptions and Setup:

- Your raw FASTQs are organized by sample under data/{sample}/.
- A 10x Genomics–compatible reference is available at a path given by config["cellranger_ref"].
- You have installed Cell Ranger, FastQC, MultiQC, and Scrublet (or a custom Python script that wraps Scrublet) accessible from your environment.
- For Scrublet QC, we assume a helper script (e.g., scripts/scrublet_qc.py) that takes as input the path to a Cell Ranger–filtered matrix folder and produces a QC report and/or filtered barcode list.
- We provide resource settings for each rule and outline a simple cluster configuration example in the comments.
- This workflow is intended as a starting point; adjust file names, paths, and resource settings as needed.

Directory Structure (Example):

```bash
.
├── Snakefile
├── config.yaml
├── scripts/
│   └── scrublet_qc.py     # Custom wrapper around scrublet for QC & doublet detection.
├── data/
│   ├── SAMPLE_1/
│   │   ├── SAMPLE_1_S1_L001_R1_001.fastq.gz
│   │   ├── SAMPLE_1_S1_L001_R2_001.fastq.gz
│   ├── SAMPLE_2/
│   │   ├── SAMPLE_2_S1_L001_R1_001.fastq.gz
│   │   ├── SAMPLE_2_S1_L001_R2_001.fastq.gz
└── results/
    ├── fastqc/
    ├── cellranger_count/
    ├── multiqc/
    └── scrublet/
```

Then, run Snakemake with a command similar to:

```bash
snakemake --jobs 50 --cluster-config cluster.json --cluster "sbatch --cpus-per-task={cluster.cores} --mem={cluster.mem} --time={cluster.time}" --latency-wait 60
```

Adjust the cluster submission command (shown here for SLURM with sbatch) as needed for your computing environment.

Summary

This extended workflow:

    Processes raw 10x FASTQ files with FastQC and aggregates QC via MultiQC.
    Aligns and quantifies reads with Cell Ranger, generating a feature–barcode matrix.
    Performs additional QC including doublet detection and filtering using a custom Scrublet step.
    Supports cluster execution via rule-specific resource requests and a sample cluster configuration file.

Customize paths, resource parameters, and the Scrublet wrapper script to best suit your experimental setup and computing environment. Enjoy your analysis!
