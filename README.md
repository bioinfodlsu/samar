# 0. Introduction
Welcome to SAMAR (Speedy, Assembly-free Method to Analyze RNA-seq expression data) -- a quick-and-easy way to perform differential expression (DE) analysis in non-model organisms.

For non-model organisms, conventional DE analysis begins by constructing a de-novo transcriptome assembly and annotating it against a high-confidence protein database -- the assembly serving as a reference for read mapping and the annotation allowing for functional analysis of DE genes. 
This pipeline provides a shortcut in which RNA-seq reads are directly aligned to the reference proteome which would have been otherwise used for annotation in the assembly-based approach. This saves you a LOT of time and other computational resources.

Specifically, SAMAR performs the following steps:

1. Uses [LAST](https://gitlab.com/mcfrith/last) to learn the alignment scoring parameters suitable for the input data, and to estimate the paired-end fragment size distribution of paired-end reads,
2. Uses LAST to directly align RNA-seq reads to the high-confidence proteome that would have been otherwise used for annotation, and generates counts,
3. Estimates the counts of reads aligning to each protein entry, using a simple rescue strategy for multi-mapping reads,
4. Uses [tximport](https://bioconductor.org/packages/release/bioc/html/tximport.html) to read in the counts (and optionally to aggregate isoform-level counts into gene-level counts),
5. Runs [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) for differential expression analysis based on the count data.

# 1. Installation
This pipeline requires the package manager **Conda** and the workflow management system **Snakemake**.
All other dependencies are handled automatically by Snakemake.

## 1.1. Install Conda 
Download Miniconda3  installer for Linux from  [here](https://docs.conda.io/en/latest/miniconda.html#linux-installers).
Installation instructions are [here](https://conda.io/projects/conda/en/latest/user-guide/install/linux.html).
Once installation is complete, you can test your Miniconda installation by running:
```
$ conda list
```

## 1.2. Install Snakemake
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

# 1.3. Download the pipeline
Download SAMAR from the online  [repository](https://bitbucket.org/refeless_rnaseq/pipeline), or using the command line:
```
git clone https://bitbucket.org/project_samar/samar.git
```
# 2. Quickstart Usage Guide

# 2.1. Input
The pipeline requires, at the very least: (1) paired-end RNA-seq fastq (or fasta) reads and (2) a reference protein database in fasta format.  These and other input parameters are specified via a YAML-format config file -- a template *config_template.yam* is provided in the config folder. Section~3 below provides a list of all parameters that can be defined in the config file.

# 2.2. Running the pipeline
After constructing a config.yaml file and with the snakemake conda environment you created earlier activated, you can call the pipeline from the top-level directory of SAMAR:
```
cd samar 
snakemake --configfile <my_config.yaml> --use-conda --cores all 
```
Please see 2.4 below for a simple example.

# 2.3. Ouput
Output is stored inside the folder specified in the *out_dir* entry of the config file. Inside it, the alignments can be found in the *last_alignments* folder, the count data in the *counts* folder, and the results of differential expression analysis in the *DEanalysis* folder. The contents of the *DEanalysis* folder are explained through an example below.

# 2.4. Example
A toy example is provided in the *test_data* folder. In this example, there are 6 RNA-seq samples, of which 3 are the "control" group replicates and 3 are the "treated" group replicates. The paths to the files and information about the experiment design is specified in  *config_test.yaml*.
From the top directory of the pipeline and with the snakemake conda environment activated, run:
```
$ snakemake --configfile testdata/config_test.yaml --use-conda --cores all 
```
Warning: This might take some time, since for the first run, Snakemake needs to install all dependencies prior to running the actual computations. 

Once the run is complete, you can check the output in the folder *test_output* which is created in the top-level directory of SAMAR. Inside the *DEanalysis* folder in *test_output*, you will find the following files:
- *DE_count_summary.txt*, which contains the count of up- and downregulated genes
- *Test_result-Group_treated_vs_control.csv*, which contains the table of the results of hypothesis testing for each gene in the reference. 
- *DESeq2_fit.RDS* which can be loaded using R, if further analysis is required.

For the example run, the first file reveals that only 1 gene is differentially expressed; 
while the second file allows the identification of this *gene* to be *UPI0000000AE2* which is upregulated (*log2FoldChange* > 0 and *padj* < 0.01).

In general, the number of hypothesis tests depends on the number of factors, the number of levels per factor, and the design (i.e. the factorial model that is hypothesized to cause a difference in gene expression). These have to be specified in the config file using the keywords: *factors* , *sample_info*, and *design*, respectively.  

Please note, especially for complicated designs, that the first sample in *sample_info* will be the reference group, and that significance testing is only done for the difference of the other levels (or interaction of these levels) vs this reference group. For advanced DE analysis like testing contrasts, a DESeq2 object is available in the loadable R data to perform such functions. 

In case you wish to perform your own differential gene expression analysis, the *counts* folder contains a *.counts* table for each sample, where the last column (NumReads) should be used for differential expression analysis.

## 3. Config
Apart from the paths to the read file and the reference proteome, the following parameters can be set in the YAML config file.

| Keyword       |   Possible values         | Default  |  Description  |
| ------------- |------------------------| ------ |  ------------|
| filetype | fasta, fastq  | fastq | -|
| training | yes, no | yes | Train the alignment scoring parameter for the input data. The default scoring scheme is BLOSUM62 |
|deanalysis | yes, no | yes |  If set to "no", the pipeline does not proceed to DE anlaysis and halts after counting. If keyword is not provided, value defaults to "yes" |
|factor| [factor1, factor2, .. ] | - | Name of factors in the study |
| sample_info: sample_i: |  [level of factor1, level of factor2, ..] | -| Describes the factor levels which the sample corresponds to |
