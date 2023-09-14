# About
This repository contains an easy-to-use pipeline for the assembly and analysis of bacterial genomes using ONT long-read or Illumina short-read technology.

# Introduction
Advances in sequencing technology during the COVID-19 pandemic has led to massive increases in the generation of sequencing data. Many bioinformatics tools have been developed to analyze this data, but very few tools can be utilized by individuals without prior bioinformatics training.

This pipeline was designed to encapsulate pre-existing tools to automate analysis of whole genome sequencing of bacteria. Installation is fast and straightfoward. The pipeline is easy to setup and contains rationale defaults, but is highly modular and configurable by more advance users.
A successful run generates consensus sequences, typing information, phylogenetic tree, and quality control report.

# Features
We anticipate the pipeline will be able to perform the following functions:
- [x] Reference-based assembly of Illumina paired-end reads
- [ ] *De novo* assembly of Illumina paired-end reads
- [ ] *De novo* assembly of ONT long reads
- [x] Run quality control checks
- [x] Variant calling using [bcftools](https://github.com/samtools/bcftools)
- [x] Maximum-likelihood phylogenetic inference of processed samples and background dataset using [iqtree](https://github.com/iqtree/iqtree2) 
- [x] MLST profiling and virulence factor detection
- [ ] Antimicrobial resistance genes and plasmid detection

# Installation
1. Install `miniconda` by running the following two command:
```commandline
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
bash Mambaforge-$(uname)-$(uname -m).sh
```

2. Clone the repository:
```commandline
git clone https://github.com/CholGen/bacpage.git
```

3. Install and activate the pipeline's conda environment:
```commandline
mamba install -y -n bacpage -f environment.yaml
mamba activate bacpage
```

4. Test the installation:
```commandline
snakemake --configfile test/test.yaml --all-temp --cores 8
```
This command should run to completion without a problem. Please create an issue if this is not the case.

# Usage
0. Navigate to the pipeline's directory.
1. Copy the `example/` directory to create a directory specifically for each batch of samples.
```commandline
cp example/ <your-project-directory-name>
```
2. Place raw sequencing reads in the `input/` directory of your project directory.
3. Record the name and absolute path of raw sequencing reads in the `sample_data.csv` found within your project directory.
4. Replace the values `<your-project-directory-name>` and `<sequencing-directory>` in `config.yaml` found within your project directory, with the absolute path of your project directory and pipeline directory, respectively.
5. Determine how many cores are available on your computer:
```commandline
cat /proc/cpuinfo | grep processor
```
6. From the pipeline's directory, run the entire pipeline on your samples using the following command:
```commandline
snakemake --configfile <your-project-directory-name>/config.yaml --cores <cores>
```
This will generate a consensus sequence in FASTA format for each of your samples and place them in `<your-project-directory-name>/results/consensus_sequences/<sample>.masked.fasta`. An HTML report containing alignment and quality metrics for your samples can be found at `<your-project-directory-name>/results/reports/qc_report.html`. A phylogeny comparing your sequences to the background dataset can be found at `<your-project-directory-name>/results/phylogeny/phylogeny.tree`
