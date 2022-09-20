This directory contains analysis scripts used to analyze the single-cell RNA-seq data and the single-cell Multiome (ATAC + GEX) data. 

## System requirements
### OS and Software dependencies
- Red Hat Enterprise Linux Server 7.9 (Maipo)
- R v4.1.2
- Python v3.7.6
- Bedtools v2.27.1
- Cell Ranger v6.0.1 (10x Genomics)
- Scrublet v0.2.1
- Seurat v4.0.4
- SingleR v1.10.0
- celldex v1.6.0
- Cell Ranger ARC v2.0.0 (10x Genomics)
- Signac v1.4.0
- GenomicRanges v1.46.1
- Azimuth v0.4.3
- ChromVAR v1.16.0
- Cicero v1.3.5
- Other R libraries are listed in sessionInfo.txt

### Hardware
Analysis was performed on National Center for Supercomputing Applications (NCSA) high-performance computing cluster - mforge. 

## Installation
git clone https://github.com/LabFunEpi/CC_multiomics.git

## Running the scripts
Scripts included in this repository were run in the following order: 
- cellranger_jobs.sh (This file contains commands that call upon scRNA_cr.sh and scMULT_cr.sh, to run Cell Ranger on each sample). 
- scRNA_cr.sh (Script that can be submitted to the job scheduler to run Cell Ranger)
- scMULT_cr.sh (Script that can be submitted to the job scheduler to run Cell Ranger ARC)
- scRNA.R (R script to combine output from Cell Ranger and perform downstream analysis)
- scMULT.R (R script to combine output from Cell Ranger ARC and perform downstream analysis)
- figurePlots.R (Downstream analysis and plot generation for figures)
