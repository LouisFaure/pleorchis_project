# pleorchis_project

This repo contains the code to reproduce the analysis of the unkown trematoda DNA, for the following study:

```
paper in publication process, citation coming soon!
```

## Setup environment

It is recommended to set the environment using `conda` or `mamba`.

```bash
mamba create -n trematoda -c bioconda -c conda-forge seqkit=2.3.1 \
             trimmomatic=0.36 spades=3.11.1 bwa=0.7.17 blast=2.12.0 \
             barrnap=0.9 python=3.9 --yes

mamba activate trematoda
pip install snakemake ffq gget
git clone https://github.com/LouisFaure/pleorchis_project
cd pleorchis_project
```

## Download data

Raw data will be made available on SRA, to download the fastq files, run the following:

```
pip install ffq
ffq SRR23342065
```


## Run pipeline

### Extract contigs

In the study, we used high quality reads from samples "Sailors" from Lane 2 of the flowcell. While not enough reads could be determined for the "Passengers" sample, we leveraged the NovaSeq Xp Workflow, which allowed us to use the Undetermined reads from Lane 1.

```bash
snakemake -j 80 assembly_Sailors_S2_L002/contigs_500.fasta assembly_Undetermined_S0_L001/contigs_500.fasta
```

The outputs are fasta files containing all contigs that are at least 500bp long.

### Extract the main contig of interest

All the contigs of the high quality "Sailors" sample are BLASTed against SSU and LSU sequences of SILVA database, restricted to Holozoa clade, only the ones with hit are kept.

```bash
snakemake -j 40 main_contig.fasta 
```
The output is a single contig in fasta format, as we have found out that only one (the longest) was having multiple hits on Holozoa clade.

### Obtain the sequences of the new 18S, 5.8S and 28S ribosomal DNA

From that main contig, we predict ribosomal DNA seqs using HMMER 3.1 (using barrnap wrapper).

```bash
snakemake -j 1 main_contig_rDNA.fasta
```

#### Compare SSU/LSU sequences from "passenger" sample to main contig rDNA

To confirm the identity of Passengers, we BLAST holozoa SSU or LSU sequences identified in Passengers samples against the identified main contig from Sailors

```bash
snakemake -j 1 blast-main_contig-SSU-Undetermined_S0_L001.txt
snakemake -j 1 blast-main_contig-LSU-Undetermined_S0_L001.txt
```

### Various explorations

#### Check read alignments on the main contig

Here, we re-map raw corrected reads against newly found main contig. 

```bash
snakemake -j 40 main_contig_mapping_Sailors.bam
```

The output bam and bai files can be visualised in [IGV browser](http://igv.org/app/).

#### Blast the main contig against Nucleotide database

```bash
snakemake -j 1 main_contig_blast_nt.csv
```

#### Similarity with Polyorchis and Brachycladium goliath

Here we use sequences included in fasta files `pleorchis_rDNA.fasta` and `goliath_rDNA.fasta`

```bash
snakemake -j 1 main_contig_blast_goliath.txt main_contig_blast_pleorchis.txt
```
