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
        """
        spades.py -1 {input.pr1} -2 {input.pr2} -o spades_{wildcards.sample} -t {threads} -m 500 --only-error-correction
        ## remove illumina index contamination from passenger sample
        if [ {wildcards.sample} = 'Passengers_S0_L001' ] 
        then
            cd spades_Passengers_S0_L001/corrected/
            trimmomatic PE Passengers_S0_L001_R1_paired.fastq.00.0_0.cor.fastq.gz Passengers_S0_L001_R2_paired.fastq.00.0_0.cor.fastq.gz u_R1_f.fastq.gz u_R1_fu.fastq.gz u_R2_f.fastq.gz u_R2_fu.fastq.gz ILLUMINACLIP:../../illumina.fa:2:30:10 -threads 8
            mv Passengers_S0_L001_R1_paired.fastq.00.0_0.cor.fastq.gz u_R1.fastq.gz # temp file
            mv Passengers_S0_L001_R2_paired.fastq.00.0_0.cor.fastq.gz u_R2.fastq.gz
            mv u_R1_f.fastq.gz Passengers_S0_L001_R1_paired.fastq.00.0_0.cor.fastq.gz
            mv u_R2_f.fastq.gz Passengers_S0_L001_R2_paired.fastq.00.0_0.cor.fastq.gz
            rm u_R1_fu.fastq.gz u_R2_fu.fastq.gz
        fi
        """


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
        "contigs_500-{sample}-seq.fasta"
    shell:
        "seqkit seq -g -m 500 {input}  > {output}"

rule generate_silva_blastdb:
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
        contigs="contigs_500-{sample}-seq.fasta",
    output:
        "blast_silva_holozoa-{sample}-{type}.txt",
    threads: 20
    shell:
        """
        blastn -query {input.contigs} -db db/silva_{wildcards.type} -num_threads 20 -outfmt 6 -max_target_seqs 5 -out blast_temp-{wildcards.sample}-{wildcards.type}.txt 
        cat blast_temp-{wildcards.sample}-{wildcards.type}.txt | grep "Holozoa" | grep -v "metagenome"> {output}
        rm blast_temp-{wildcards.sample}-{wildcards.type}.txt
        """

rule get_hit_contigs:
    # Extract contigs that succesfully matched to LSU/SSU for inspection
    input:
        blast="blast_silva_holozoa-{sample}-{subunit}.txt",
        contigs="contigs_500-{sample}-seq.fasta"
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
    # Extract the main contig that match to both LSU and SSU seqs from trematoda
    input:
        blast_ssu="blast_silva_holozoa-{sample}-SSU.txt",
        blast_lsu="blast_silva_holozoa-{sample}-LSU.txt",
        contigs="contigs_500-{sample}-seq.fasta"
    output:
        fasta="main_contig-{sample}-seq.fasta",
        fai="main_contig-{sample}-seq.fasta.fai",
    shell:
        """
        cat blast_silva_holozoa-{wildcards.sample}-SSU.txt | cut -f1 | uniq | sort > {wildcards.sample}-SSU.txt
        cat blast_silva_holozoa-{wildcards.sample}-LSU.txt | cut -f1 | uniq | sort > {wildcards.sample}-LSU.txt
        seqname=$(comm -12 {wildcards.sample}-LSU.txt {wildcards.sample}-SSU.txt)
        bash get_seq.sh $seqname {input.contigs} > {output.fasta} 
        samtools faidx {output.fasta}
        rm {wildcards.sample}-SSU.txt {wildcards.sample}-LSU.txt
        """


rule predict_rRNA:
    # Predict ribosomal DNA seqs using HMMER 3.1 (using barrnap wrapper)
    input:
        "{seq}-{sample}-seq.fasta",
    output:
        "{seq}-{sample}-rRNA.fasta",
    shell:
        """
        # When 5.8S is present, barrnap wrongly assign the start of 28S at the same location
        # see https://github.com/tseemann/barrnap/issues/47 for more info
        # so here is a hack where we first detect 5_8S, then mask it to detect the start of 28S
        barrnap --quiet --kingdom euk {input} > {wildcards.seq}-{wildcards.sample}-temp.gff3
        seq_=$(cat {wildcards.seq}-{wildcards.sample}-temp.gff3 | grep "Name=5_8S")
        start=$(echo "$seq_" | awk '{{print $4}}')
        end=$(echo "$seq_" | grep "Name=5_8S" | awk '{{print $5}}')

        seq=$(grep -v ">" {input} | tr -d '\\n')
        subseq=${{seq:($start - 1):($end - $start + 1)}}
        subseq_length=$(($end - $start + 1))
        # Generate the 'N' replacement string
        replacement=$(printf '%*s' "$subseq_length" | tr ' ' 'N')

        grep ">" {input} > {wildcards.seq}-{wildcards.sample}-masked.fasta

        echo $seq | sed "s/$subseq/$replacement/g" >> {wildcards.seq}-{wildcards.sample}-masked.fasta

        barrnap --quiet --kingdom euk {wildcards.seq}-{wildcards.sample}-masked.fasta > {wildcards.seq}-{wildcards.sample}.gff3
        grep "5_8S" {wildcards.seq}-{wildcards.sample}-temp.gff3 >> {wildcards.seq}-{wildcards.sample}.gff3
        rm {wildcards.seq}-{wildcards.sample}-masked.fasta {wildcards.seq}-{wildcards.sample}-temp.gff3
        bedtools getfasta -fi {input} -bed {wildcards.seq}-{wildcards.sample}.gff3 -fo {output}
        """


rule check_alignment:
    # Re-map raw corrected reads against newly found main contig. Generate files needed to inspect aligment in
    # IGV browser http://igv.org/app/ (Genome: output from extract_main_contig, Tracks: present ouputs)
    input:
        contig="{seq}-{sample}-seq.fasta",
        r1="spades_{sample}/corrected/{sample}_R1_paired.fastq.00.0_0.cor.fastq.gz",
        r2="spades_{sample}/corrected/{sample}_R2_paired.fastq.00.0_0.cor.fastq.gz",
    output:
        bam="{seq}-{sample}-mapping.bam",
        bai="{seq}-{sample}-mapping.bam.bai",
    shell:
        """
        mkdir -p bwa
        bwa index -p bwa/{wildcards.seq}-{wildcards.sample} {input.contig}
        bwa mem -t 40 bwa/{wildcards.seq}-{wildcards.sample} {input.r1} {input.r2} > {wildcards.seq}-{wildcards.sample}-alignment.sam
        samtools view -@ 10 -F 4 {wildcards.seq}-{wildcards.sample}-alignment.sam -o {wildcards.seq}-{wildcards.sample}-mapped.bam
        samtools sort {wildcards.seq}-{wildcards.sample}-mapped.bam > {output.bam}
        rm {wildcards.seq}-{wildcards.sample}-mapped.bam {wildcards.seq}-{wildcards.sample}-alignment.sam
        samtools index {output.bam}
        """

rule setup_taxonomy:
    output:
        names="taxo/names.dmp",
        nodes="taxo/nodes.dmp",
        taxtrace="tax_trace.pl",
    shell:
        """
        wget https://raw.githubusercontent.com/theo-allnutt-bioinformatics/scripts/master/tax_trace.pl
        wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
        mkdir -p taxo && tar xvf taxdump.tar.gz -C taxo
        """

rule feature_counts:
    input:
        bam="{seq}-{sample}-mapping.bam",
        contig="{seq}-{sample}-seq.fasta",
    output:
        fcounts="{seq}-{sample}-fcounts.tsv"
    shell:
        """
        samtools view -@ 10 -F 2432 {input.bam} -o {wildcards.seq}-{wildcards.sample}-unique.bam
        samtools index {wildcards.seq}-{wildcards.sample}-unique.bam
        grep "^>" {input.contig} | awk 'sub(/^>/, "")' > {wildcards.seq}-{wildcards.sample}.txt
        cat /dev/null >| {output.fcounts}
        while read s; do
            count=$(samtools view {wildcards.seq}-{wildcards.sample}-unique.bam "$s" -c)
            echo -e "$s\t$count" >> {output.fcounts}
        done <{wildcards.seq}-{wildcards.sample}.txt
        rm {wildcards.seq}-{wildcards.sample}.txt {wildcards.seq}-{wildcards.sample}-unique.bam*
        """

rule summarize_stats:
    input:
        contig="{seq}-{sample}-seq.fasta",
        bam="{seq}-{sample}-mapping.bam",
        fcounts="{seq}-{sample}-fcounts.tsv",
        names="taxo/names.dmp",
        nodes="taxo/nodes.dmp",
        taxtrace="tax_trace.pl",
        db="db/nt.nto",
    output:
        fcounts_class="{seq}-{sample}-fcounts_classified.csv",
        summary="{seq}-{sample}-summary.csv",
    shell:
        """
        ## Get taxonomy
        blastn -num_threads 40 -db db/nt -query {input.contig} -out {wildcards.seq}-{wildcards.sample}-blast.txt -outfmt '6 qseqid sseqid evalue pident stitle staxids' -max_target_seqs 1 -max_hsps 1
        cut -f1,6 {wildcards.seq}-{wildcards.sample}-blast.txt  > {wildcards.seq}-{wildcards.sample}-taxo_temp.txt
        perl tax_trace.pl taxo/nodes.dmp taxo/names.dmp {wildcards.seq}-{wildcards.sample}-taxo_temp.txt {wildcards.seq}-{wildcards.sample}-taxo_id.txt
        cat {wildcards.seq}-{wildcards.sample}-taxo_id.txt | cut -f1 > {wildcards.seq}-{wildcards.sample}-matched.txt
        grep 'artificial' {wildcards.seq}-{wildcards.sample}-taxo_id.txt | cut -f1 > {wildcards.seq}-{wildcards.sample}-artificial.txt
        grep 'Viruses' {wildcards.seq}-{wildcards.sample}-taxo_id.txt | cut -f1 > {wildcards.seq}-{wildcards.sample}-viruses.txt
        grep 'Bacteria' {wildcards.seq}-{wildcards.sample}-taxo_id.txt | cut -f1 > {wildcards.seq}-{wildcards.sample}-bacteria.txt
        grep 'Trematoda' {wildcards.seq}-{wildcards.sample}-taxo_id.txt | cut -f1 > {wildcards.seq}-{wildcards.sample}-trematoda.txt
        grep -v 'Trematoda' {wildcards.seq}-{wildcards.sample}-taxo_id.txt| grep Eukaryota | cut -f1 > {wildcards.seq}-{wildcards.sample}-other_euk.txt
        
        ## Get mean coverage
        average-coverage.py {input.bam} > {wildcards.seq}-{wildcards.sample}-coverage.txt

        ## Summarize
        python generate_stats.py {wildcards.seq}-{wildcards.sample}
        rm {wildcards.seq}-{wildcards.sample}-*.txt
        """

rule get_coverage:
    # Compute mean read coverage for each sequences (each contigs)
    input:
        bam="{seq}-{sample}-mapping.bam",
        fcounts_class="{seq}-{sample}-fcounts_classified.csv",
    output:
        tsv="{seq}-{sample}-coverage.tsv"
    shell:
        """
        average-coverage.py {input.bam} > {output.tsv}

        """

rule setup_nt_blastdb:
    # WARNING: this step will take a lot of time, as it will download the whole eukaryote blast database
    output:
        db="db/nt.nto",
    shell:
        """
        cd db
        wget https://raw.githubusercontent.com/jrherr/bioinformatics_scripts/master/perl_scripts/update_blastdb.pl
        perl update_blastdb.pl --decompress nt
        """

rule blast_nt:
    # BLAST the contig against whole NCBI nt databse (using gget tool)
    input:
        "main_contig-{sample}.fasta",
    output:
        "main_contig-{sample}-blast_nt.csv",
    shell: "bash blast_nt.sh {input} {output}"

rule trematoda_rRNA_blast:
    # Explore ribosomal DNA similarites with Polyorchis and Brachycladium goliath
    input:
        seqs="{specie}_rRNA.fasta",
        contig="main_contig-{sample}.fasta",
    output:
        "main_contig-{sample}-blast-{specie}.txt"
    shell:
        """
        makeblastdb -in {input.seqs} -dbtype nucl -out db/{wildcards.specie}_rRNA
        blastn -query {input.contig} -db db/{wildcards.specie}_rRNA -num_threads 1 -outfmt 6 -max_target_seqs 2 > {output}
        """
        
        
rule common_hits:
    # Extract overlapping Trematoda contigs between the two samples
    input:
        sailors="contigs_500-Sailors_S2_L002-seq.fasta",
        passeng="contigs_500-Passengers_S0_L001-seq.fasta",
        s_stat="contigs_500-Passengers_S0_L001-fcounts_classified.csv",
        p_stat="contigs_500-Sailors_S2_L002-fcounts_classified.csv",
    output:
        "common_hits.csv",
    shell:
        """
        makeblastdb -in contigs_500-Sailors_S2_L002-seq.fasta -dbtype nucl -out db/contigs_500-Sailors_S2_L002
        python get_overlap.py
        """
