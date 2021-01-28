# 0. Introduction
This pipeline allows performing rapid differential expression (DE) analysis in non-model organisms, without the need for computationally expensive de-novo transcriptome assembly. 

For non-model organisms, the conventional strategy for DE analysis has been to begin by constructing a de-novo transcriptome assembly and annotating it against a high-confidence protein database -- the assembly serving as a reference for read mapping and the annotation allowing for functional analysis of genes found to have DE. 

This pipeline takes in paired-end RNA-seq reads and a reference proteome (which would have been used for annotation in the assembly-based approach) and does the following:
1. Uses [LAST](http://last.cbrc.jp) to learn the alignment scoring parameters suitable for the input data, and to estimate the paired-end fragment size distribution
2. Uses [LAST](http://last.cbrc.jp) to directly align RNA-seq reads to the high-confidence proteome that would have been otherwise used for annotation, and generates counts
3. Estimates the counts of reads aligning to each protein entry, using a rescue strategy for multi-mapping reads.
4. Runs [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) for differential expression analysis based on the count data.

# 1. Installation
This pipeline requires the package manager **Conda** and the workflow management system **Snakemake**.
All other dependencies are handled automatically by Snakemake.
## 1.1 Install Conda 
Download Miniconda3  installer for Linux from  [here](https://docs.conda.io/en/latest/miniconda.html#linux-installers).
Installation instructions are [here](https://conda.io/projects/conda/en/latest/user-guide/install/linux.html).
Once installation is complete, you can test your Miniconda installation by running:
```
$ conda list
```
## 1.2 Install Snakemake
Snakemake recommends installation via Conda:
```
$ conda install -c conda-forge mamba
$ mamba create -c conda-forge -c bioconda -n snakemake snakemake
```
This creates an isolated enviroment containing the latest Snakemake. To activate it:
```
$ conda activate snakemake
```
To test snakemake installation 
```
$ snakemake --help
```
# 1.3 Download the pipeline
Download our pipeline from the online [repository](https://bitbucket.org/refeless_rnaseq/pipeline), or using the command line:
```
git clone https://bitbucket.org/refeless_rnaseq/pipeline.git
```
# 2. Quickstart Usage Guide

# 2.1 Input
The pipeline requires, at the very least: (1) paired-end RNA-seq reads in fastq format and (2) a reference protein database in fasta format. 
These are specified via a config file, a template for which is provided in the config folder.

# 2.2 Example
As an example, suppose we wish to compute counts for the sets of paired-end reads in the *test_data* folder, using the reference proteome in the same folder. 
The folder also contains a config file *config_test.yaml* file specifying the following: paths to the input data, the experiment design for differential expression analysis, and the output folder.
To run the pipeline, from the pipeline root directory and with the snakemake conda environment activated, run:
```
$ snakemake --configfile testdata/config_test.yaml --use-conda --cores all 
```
This might take some time since for the first run, Snakemake needs to install all dependencies prior to running the actual computations.
# 2.3 Ouput
The results of the differential expression analysis can be found inside the *DEanalysis* folder. 

Other intermediate output data can also be found in the output folder.  For example, the count data can be found inside the *counts* folder, in case you wish to perform your own differential gene expression analysis.
The count file contains five columns, of which the last column (NumReads) should be used for differential expression analysis. 


