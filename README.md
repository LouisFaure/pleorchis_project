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
             barrnap=0.9 python=3.9 bamtocov=2.7.0 seqtk=1.4 \
             make=4.3 gxx=12.2.0 --yes

mamba activate trematoda
pip install snakemake ffq gget
git clone https://github.com/LouisFaure/pleorchis_project
cd pleorchis_project
```

## Download data

Raw data will be made available on SRA, to download the fastq files, run the following:

```bash
ffq SRR23342065
```

### Saving undertermined reads from lane 1 as passengers

Data from passengers was of poor quality (not much reads associated to illumina indexes recovered), 
hinting on a possible issue while prepary libraries.
As the 2 libraries were separately loaded into each lane of the Illumina SP flow cell according to Xp workflow, we could 
recover the data from Undetermined reads from lane 1 as data containing sequences from passengers sample.

```bash
cd data
mv Undetermined_S0_L001_R1_001.fastq.gz Passengers_S0_L001_R1_001.fastq.gz
mv Undetermined_S0_L001_R2_001.fastq.gz Passengers_S0_L001_R2_001.fastq.gz
cd ..
```

## Run pipeline

### Extract contigs

This is the main part of the analysis, and also the most time consuming and compute intensive, the following steps are applied:
1. Trim raw reads with `trimmomatic`.
2. Do a first pass of `spades`, for error correction of the trimmed reads.
  - For Undetermined, add one step of `trimmomatic`, trimming reads using *illumina.fa* provided.
3. Do a second pass of `spades`, this time for assembly of contigs.
4. Filter contigs by removing short ones.

```bash
snakemake -j 80 contigs_500-Sailors_S2_L002.fasta \
                contigs_500-Passengers_S0_L001.fasta
```

The outputs are fasta files containing all contigs that are at least 500bp long.

### Extract the main contig of interest

All the contigs of the high quality "Sailors" sample are BLASTed against SSU and LSU sequences of SILVA database, restricted to Holozoa clade, only the ones with hit are kept.

```bash
snakemake -j 40 main_contig-Sailors_S2_L002-seq.fasta \ 
                main_contig-Passengers_S0_L001-seq.fasta 
```
The output is a single contig in fasta format, as we have found out that only one (the longest) was having multiple hits with holozoan LSU/SSU.

### Obtain the sequences of the new 18S, 5.8S and 28S ribosomal DNA

From that main contig, we predict ribosomal DNA seqs using HMMER 3.1 (using barrnap wrapper).

```bash
snakemake -j 1 main_contig-Sailors_S2_L002-rRNA.fasta \
               main_contig-Passengers_S0_L001-rRNA.fasta 
```
The output sequences are the ones that have been submitted to GenBank.


### Various explorations

#### Check read alignments on the main contig

Here, we re-map raw corrected reads against newly found main contig. 

```bash
snakemake -j 40 main_contig-Sailors_S2_L002-mapping.bam
```

The output bam and bai files can be visualised in [IGV browser](http://igv.org/app/).

#### Blast the main contig against Nucleotide database

```bash
snakemake -j 1 main_contig-Sailors_S2_L002-blast_nt.csv
```

#### Similarity with Polyorchis and Brachycladium goliath

Here we use sequences included in fasta files `pleorchis_rDNA.fasta` and `goliath_rDNA.fasta`

```bash
snakemake -j 1 main_contig-Sailors_S2_L002-blast-goliath.txt main_contig-Sailors_S2_L002-blast-pleorchis.txt
```

#### Generate read stats

```bash
snakemake -j 40 contigs_500-Sailors_S2_L002-fcounts_classified.csv \
                contigs_500-Passengers_S0_L001-fcounts_classified.csv
```

#### Study contig overlap between the two samples

```bash
snakemake -j 1 common_hits.csv
```