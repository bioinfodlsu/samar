# 0. Introduction
This pipeline allows performing rapid differential expression (DE) analysis in non-model organisms, without the need for computationally expensive de-novo transcriptome assembly. 

For non-model organisms, the conventional strategy for DE analysis has been to first construct a de-novo transcriptome assembly and annotate it against a high-confidence protein database -- the assembly serving as a reference for read mapping and the annotation allowing for functional analysis of genes found to have DE. 

This pipeline uses [LAST](www.last.cbrc.jp) to directly RNA-seq reads to the high-confidence proteome that would have been otherwise used for annotation, and generates counts that can be fed to a counts-based differential expression analysis tool (e.g. [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)).
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
# 2. Usage

# 2.1 Input
The pipeline requires, at the very least: (1) paired-end RNA-seq reads in fastq format and (2) a reference protein database in fasta format. Additionally, you may also provide the path to an output folder where the results will be placed. The file/folder paths need to be specified in a config file (YAML) format. A template config file is provided in the *config* folder of the repository. 

# 2.2 Example
As an example, suppose we wish to compute counts for the sets of paired-end reads in the *sample_data* folder, using the reference proteome in the same folder. The folder also contains a config.yaml file. From the pipeline root directory, with the snakemake conda environment activated, we run the pipeline as follows:

```
$ snakemake --configfile sample_data/config.yaml --use-conda --cores all 
```
This might take some time -- especially if it is the first ever run, Snakemake needs to install all dependencies prior to running the actual computations. This is a good time for a short break.

# 2.3 Ouput
The count data required for differential gene expression analysis is inside the *last_counts* folder inside the output folder (if not specified in the config file, the default output location is the root directory of the repository). To adhere to conventional data formats, the count file contains five columns, of which the last column (NumReads) should be used for differential expression analysis.





