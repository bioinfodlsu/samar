# Installation

## Install Miniconda 
Donwload Miniconda3 installer for Linux 
```
https://docs.conda.io/en/latest/miniconda.html#linux-installers
```
Verify the downloaded file in your teminal window 
```
$ sha256sum <filename>
```
Run installer in your terminal
```
$ bash Miniconda3-latest-Linux-x86_64.sh
```
To test Miniconda installation run
```
$ conda list
```
## Install Snakemake
### Installation through Miniconda (Preferred)
Install conda package manager "mamba"
```
$ conda install -c conda-forge mamba
```
Install snakemake using mamba
```
$ mamba create -c conda-forge -c bioconda -n snakemake snakemake
```
To run snakemake environment
```
$ conda activate snakemake
```
To test snakemake installation 
```
$ snakemake --help
```
To deactivate snakemake environment 
```
$ conda deactivate
```
### Installation through Pip
Install snakemake 
```
$ pip3 install snakemake
```
To test snakemake installation
```
$ snakemake --help
```
# Usage
## Update Config File
A sample config file with sample data is present in the repository 
## Running the pipeline
Within the snakemake environment (for conda) or directly in the terminal (for pip)
```
$ snakemake --use-conda --cores all 
```
To dry run the pipeline 
```
$ snakemake --use-conda --cores all -n
```

