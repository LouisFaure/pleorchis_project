import sys
import pandas as pd
prefix = sys.argv[1]+"-"

fcounts=pd.read_csv(prefix+"fcounts.tsv",sep="\t",header=None,index_col=0,names=["counts"])

fcounts.index.name="contig"

fcounts["coverage"]=pd.read_table(prefix+"coverage.txt").iloc[:,1].values

matched=pd.read_table(prefix+"matched.txt",header=None).values.ravel()
artificial=pd.read_table(prefix+"artificial.txt",header=None).values.ravel()
trematoda=pd.read_table(prefix+"trematoda.txt",header=None).values.ravel()
euk=pd.read_table(prefix+"other_euk.txt",header=None).values.ravel()
viruses=pd.read_table(prefix+"viruses.txt",header=None).values.ravel()
bacteria=pd.read_table(prefix+"bacteria.txt",header=None).values.ravel()

fcounts["class"]="unkown"
fcounts.loc[fcounts.index.isin(matched),"class"]="matched"
fcounts.loc[fcounts.index.isin(artificial),"class"]="artificial"
fcounts.loc[fcounts.index.isin(trematoda),"class"]="Trematoda"
fcounts.loc[fcounts.index.isin(euk),"class"]="Other Eukaryotes"
fcounts.loc[fcounts.index.isin(viruses),"class"]="Virus"
fcounts.loc[fcounts.index.isin(bacteria),"class"]="bacteria"

summary=fcounts.groupby("class")["counts"].sum()

fcounts.groupby("class")["counts"].sum().to_frame().to_csv(prefix+"summary.csv")
fcounts.to_csv(prefix+"fcounts_classified.csv")
