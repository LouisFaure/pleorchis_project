# pleorchis_project

This repo contains the code to reproduce the analysis of the unkown trematoda DNA.

## Setup environment

It is recommended to set the environment using `conda` or `mamba`.

```bash
mamba create -n trematoda -c bioconda -c conda-forge seqkit=2.3.1 \
             trimmomatic=0.36 spades=3.11.1 bwa=0.7.17 blast=2.12.0 \
             barrnap=0.9 python=3.9 --yes

mamba activate trematoda
pip install snakemake ffq gget
```

## Download data
Data will be made available on SRA, to obtain it, we recommned to use `ffq`:

```
pip install ffq
ffq SRR23342065
```

## Run pipeline

### Extract contigs

```bash
snakemake -j 80 assembly_Sailors_S2_L002/contigs_500.fasta assembly_Undetermined_S0_L001/contigs_500.fasta
```

### Obtain the sequences of the new 18S, 5.8S and 28S ribosomal DNA

```bash
snakemake -j 1 main_contig_rDNA.fasta
```

#### Compare SSU/LSU sequences from "passenger" sample to main contig rDNA

```bash
snakemake -j 1 blast-main_contig-SSU-Undetermined_S0_L001.txt
snakemake -j 1 blast-main_contig-LSU-Undetermined_S0_L001.txt
```

