# 0. Introduction
This pipeline allows performing rapid differential expression (DE) analysis in non-model organisms, without the need for computationally expensive de-novo transcriptome assembly. 

For non-model organisms, the conventional strategy for DE analysis has been to begin by constructing a de-novo transcriptome assembly and annotating it against a high-confidence protein database -- the assembly serving as a reference for read mapping and the annotation allowing for functional analysis of genes found to have DE. 

This pipeline uses [LAST](http://last.cbrc.jp) to directly align RNA-seq reads to the high-confidence proteome that would have been otherwise used for annotation, and generates counts that can be fed to a counts-based differential expression analysis tool (e.g. [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)).
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
# 2.2 Example
As an example, suppose we wish to compute counts for the sets of paired-end reads in the *test_data* folder, using the reference proteome in the same folder. 
The folder also contains a config file *config_test.yaml* file specifying the paths to the input data, and where the output should be stored.
To run the pipeline, from the pipeline root directory and with the snakemake conda environment activated:
```
$ snakemake --configfile test_data/config_test.yaml --use-conda --cores all 
```
This might take some time since Snakemake needs to install all dependencies prior to running the actual computations.
# 2.3 Ouput
The count data required for differential gene expression analysis can be found inside the *counts* folder in the output folder. If the output folder is not specified in the config file, the default output location is the root of the project directory. 
The count file contains five columns, of which the last column (NumReads) should be used for differential expression analysis.

