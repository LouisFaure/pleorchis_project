import pandas as pd
from Bio import SeqIO
from io import StringIO
import subprocess

command="blastn -query contigs_500-Passengers_S0_L001-seq.fasta -db db/contigs_500-Sailors_S2_L002 -num_threads 20 -outfmt 6 -max_target_seqs 1 -perc_identity 98 -max_hsps 1"
result = subprocess.run(command,shell=True, capture_output=True, text=True)
blst=pd.read_table(StringIO(result.stdout),header=None)

blst=blst.loc[:,:1]
blst.columns=["passengers","sailors"]


fcounts=pd.read_csv("contigs_500-Passengers_S0_L001-fcounts_classified.csv")
fcounts["common"]=fcounts.contig.isin(blst["passengers"])
seq_names=fcounts.contig.loc[((fcounts["class"]=="Trematoda") & fcounts.common)].values

with open("contigs_500-Passengers_S0_L001.fasta", "r") as original, open("contigs_500-Passengers_S0_L001-common.fasta", "w") as output:
    sequences = SeqIO.parse(original, "fasta")
    for sequence in sequences:
        if sequence.id in seq_names:
            SeqIO.write(sequence, output, "fasta")

command="blastn -query contigs_500-Passengers_S0_L001-common.fasta -db db/nt -num_threads 60 -outfmt 6 -max_target_seqs 1 -max_hsps 1"
result = subprocess.run(command,shell=True, capture_output=True, text=True)
from Bio import Entrez
Entrez.email = "louis.faure@meduniwien.ac.at"
acc=pd.read_table(StringIO(result.stdout),header=None)[1].values
desc=[]
for accession in acc:
	handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
	record = SeqIO.read(handle, "genbank")
	desc.append(record.description)
	handle.close()

blst=blst.set_index("passengers").loc[seq_names].reset_index()
blst["description"]=desc
blst.to_csv("common_hit.csv", index=False)
