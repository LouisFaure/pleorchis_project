from snakemake.shell import shell


rule all:
    input:
        "assembly_{sample}/contigs_500.fasta"

rule trimmomatic:
    # Trim raw reads
    input:
        r1="data/{sample}_R1_001.fastq.gz",
        r2="data/{sample}_R2_001.fastq.gz",
    output:
        pr1="trimmed/{sample}_R1_paired.fastq.gz",
        pr2="trimmed/{sample}_R2_paired.fastq.gz",
        ur1="trimmed/{sample}_R1_unpaired.fastq.gz",
        ur2="trimmed/{sample}_R2_unpaired.fastq.gz",
    threads: 8
    shell: 
        "trimmomatic PE {input.r1} {input.r2} {output.pr1} {output.ur1} {output.pr2} {output.ur2} LEADING:10 TRAILING:10 SLIDINGWINDOW:5:20 MINLEN:50 -threads 8"


rule spades_error_correct:
    # First run of spades before assembly, generate error corrected reads
    input:
        pr1="trimmed/{sample}_R1_paired.fastq.gz",
        pr2="trimmed/{sample}_R2_paired.fastq.gz",
    output:
        "spades_{sample}/corrected/{sample}_R1_paired.fastq.00.0_0.cor.fastq.gz",
        "spades_{sample}/corrected/{sample}_R2_paired.fastq.00.0_0.cor.fastq.gz",
    threads: 40
    shell:
        "spades.py -1 {input.pr1} -2 {input.pr2} -o spades_{wildcards.sample} -t {threads} -m 500 --only-error-correction"


rule spades_assembly:
    # Do actual assembly to get contigs
    input:
        ecr1="spades_{sample}/corrected/{sample}_R1_paired.fastq.00.0_0.cor.fastq.gz",
        ecr2="spades_{sample}/corrected/{sample}_R2_paired.fastq.00.0_0.cor.fastq.gz",
    output:
        "assembly_{sample}/contigs.fasta"
    threads: 40
    shell:
        "spades.py -1 {input.ecr1} -2 {input.ecr2} -o assembly_{wildcards.sample} -t 40 --only-assembler"

rule filter_contigs:
    # Filter out sequence shorter than 500bp from the contigs
    input:
        "assembly_{sample}/contigs.fasta"
    output:
        "assembly_{sample}/contigs_500.fasta"
    shell:
        "seqkit seq -g -m 300 {input}  > {output}"

rule generate_blastdb:
    # Prepare blastdb with Silva databases containing SSU and LSU sequences
    output:
        "db/silva_{type}.nsq",
    shell:
        """
        mkdir -p db
        wget -O - https://www.arb-silva.de/fileadmin/silva_databases/release_138.1/Exports/SILVA_138.1_{wildcards.type}Ref_NR99_tax_silva.fasta.gz | zcat | tr -d "[:blank:]" > db/silva_{wildcards.type}.fasta
        makeblastdb -in db/silva_{wildcards.type}.fasta -dbtype nucl -out db/silva_{wildcards.type}
        """


rule blast_contigs:
    # BLAST contigs against SSU and LSU sequences, keep only from Holozoa clade
    input:
        db="db/silva_{type}.nsq",
        contigs="assembly_{sample}/contigs_500.fasta",
    output:
        "blast_silva_holozoa-{sample}-{type}.txt",
    threads: 20
    shell:
        """
        blastn -query {input.contigs} -db db/silva_{wildcards.type} -num_threads 20 -outfmt 6 -max_target_seqs 5 -out blast_temp.txt 
        cat blast_temp.txt | grep "Holozoa" | grep -v "metagenome"> {output}
        rm blast_temp.txt
        """

rule get_hit_contigs:
    # Extract contigs that succesfully matched to LSU/SSU for inspection
    input:
        blast="blast_silva_holozoa-{sample}-{subunit}.txt",
        contigs="assembly_{sample}/contigs_500.fasta"
    output:
        "holozoa-{subunit}_contigs-{sample}.fasta",
    shell:
        """
        array=($(cat {input.blast} | awk -F'\t' '{{print $1}}' | uniq))
        for seqname in "${{array[@]}}"
        do
            bash get_seq.sh $seqname {input.contigs} >> temp
        done
        awk 'NF' temp > {output} && rm temp
        """

rule extract_main_contig:
    # From the Sailors sample, extract the main contig that match to both LSU and SSU seqs from trematoda
    input:
        blast="blast_silva_holozoa-Sailors_S2_L002-SSU.txt",
        contigs="assembly_Sailors_S2_L002/contigs_500.fasta"
    output:
        fasta="main_contig.fasta",
        fai="main_contig.fasta.fai",
    shell:
        """
        seqname=$(head -n 1 {input.blast} | awk '{{print $1}}')
        bash get_seq.sh $seqname {input.contigs} > {output.fasta} 
        samtools faidx {output.fasta}
        """

rule check_alignment:
    # Re-map raw corrected reads against newly found main contig. Generate files needed to inspect aligment in
    # IGV browser http://igv.org/app/ (Genome: output from extract_main_contig, Tracks: present ouputs)
    input:
        contig="main_contig.fasta",
        r1="spades_{sample}/corrected/{sample}_R1_paired.fastq.00.0_0.cor.fastq.gz",
        r2="spades_{sample}/corrected/{sample}_R2_paired.fastq.00.0_0.cor.fastq.gz",
    output:
        bam="main_contig_mapping_{sample}.bam",
        bai="main_contig_mapping_{sample}.bam.bai",
    shell:
        """
        mkdir -p bwa
        bwa index -p bwa/main_contig {input.contig}
        bwa mem -t 40 bwa/main_contig {input.r1} {input.r2} > alignment.sam
        samtools view -@ 10 -F 4 alignment.sam -o mapped.bam
        samtools sort mapped.bam > {output.bam}
        rm mapped.bam alignment.sam
        samtools index {output.bam}
        """

rule blast_nt:
    # BLAST the contig against whole NCBI nt databse (using gget tool)
    input:
        "main_contig.fasta",
    output:
        "main_contig_blast_nt.csv",
    shell: "bash blast_nt.sh {input} {output}"

rule trematoda_rDNA_blast:
    # Explore ribosomal DNA similarites with Polyorchis and Brachycladium goliath
    input:
        seqs="{specie}_rDNA.fasta",
        contig="main_contig.fasta",
    output:
        "main_contig_blast_{specie}.txt"
    shell:
        """
        makeblastdb -in {input.seqs} -dbtype nucl -out db/{wildcards.specie}_rDNA
        blastn -query {input.contig} -db db/{wildcards.specie}_rDNA -num_threads 1 -outfmt 6 -max_target_seqs 2 > {output}
        """

rule predict_rDNA:
    # Predict ribosomal DNA seqs using HMMER 3.1 (using barrnap wrapper)
    input:
        "main_contig.fasta",
    output:
        "main_contig_rDNA.fasta",
    shell:
        """
        barrnap --quiet --kingdom euk {input} > main_contig.gff3
        # if 5.8S is present, barrnap wrongly assign the start of 28S at the same location
        # see https://github.com/tseemann/barrnap/issues/47 for more info
        # so we need to truncate the sequence to ingore 5.8S to detect the start of 28S
        end=$(cat main_contig.gff3 | grep "Name=5_8S" | awk '{{print $5}}')
        name_contig=$(cat {input} | head -n 1 | awk -F '>' '{{print $2}}')

        samtools faidx {input} "$name_contig":-"$end" > main_contig_part1.fasta
        samtools faidx {input} "$name_contig":"$end"- > main_contig_part2.fasta

        barrnap --outseq main_contig_rDNA_part1.fasta --quiet --kingdom euk main_contig_part1.fasta > main_contig_part1.gff3
        barrnap --outseq main_contig_rDNA_part2.fasta --quiet --kingdom euk main_contig_part2.fasta > main_contig_part2.gff3

        cat main_contig_rDNA_part1.fasta main_contig_rDNA_part2.fasta > {output}

        rm main_contig_part*.fasta* main_contig*gff3 main_contig_rDNA_*
        """

rule match_passengers_to_sailors:
    # Confirm identity of Passengers by blasting holozoa SSU or LSU sequences identified in Passengers samples
    # against the identified main contig from Sailors
    input:
        seqs="holozoa-{type}_contigs-Undetermined_S0_L001.fasta",
        contig="main_contig.fasta",
    output:
        "blast-main_contig-{type}-Undetermined_S0_L001.txt",
    shell:
        """
        makeblastdb -in {input.contig} -dbtype nucl -out db/main_contig
        blastn -query {input.seqs} -db db/main_contig -num_threads 1 -outfmt 6 -max_target_seqs 2 > {output}
        """

