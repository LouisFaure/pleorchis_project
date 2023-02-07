# pleorchis_project

This repo contains the code to reproduce the analysis of the unkown trematoda DNA.

## Setup environment

It is recommended to set the environment using `conda` or `mamba`.

```bash
mamba create -n de_novo -c bioconda -c conda-forge fastqc=0.11.5 seqkit=2.3.1 \
             trimmomatic=0.36 spades=3.11.1 megahit=1.1.2 quast=5.0.2 \
             bowtie2=2.2.5 anvio=5.5.0 centrifuge=1.0.4 java-jdk=8.0.112 \
             blast gget barrnap --yes
mamba activate de_novo
pip install snakemake
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
snakemake -j 80 \ 
	assembly_Sailors_S2_L002/contigs_500.fasta \
	assembly_Undetermined_S0_L001/contigs_500.fasta
```

### Obtain the sequences of the new 18S, 5.8S and 28S ribosomal DNA

```bash
snakemake main_contig_rDNA.fasta
```

