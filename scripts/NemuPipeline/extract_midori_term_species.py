import pandas as pd
from Bio import SeqIO

PATH_TO_TABLE = '/mnt/data/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/MIDORI/midori_CO1_table.csv'
PATH_TO_MIDORI = '/home/gabs/Documents/lab/TermitesAndCockroaches/MIDORI/CO1_longest.faa'

df = pd.read_csv(PATH_TO_TABLE)

df = df[df.Orders == 'Blattodea_85823']
df = df[df.CO1 >= 20] #should use 10, but we need to get rid of a tad more species

species_list = list(df['Species'])
ter_fams = ['Termitidae_46569', 'Rhinotermitidae_36985', 'Kalotermitidae_46562', 'Hodotermitidae_70920'] #Relevant only for midori_CO1_table.csv

for entry in SeqIO.parse(PATH_TO_MIDORI, 'fasta'):
    taxonomy = entry.id.split(";")
    sp_name = taxonomy[7]
    id = taxonomy[0].split('.')
    fam = taxonomy[5]
    if sp_name in species_list:
        if fam in ter_fams:
            with open(f"/mnt/data/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/NemuPipeline/TermAndCock/midori_separate_ids/termite/{sp_name}_CO1.fasta", "w") as fout:
                fout.write(f">{id[0]}.{id[1]} [{sp_name}]\n{entry.seq}")
        else:
            with open(f"/mnt/data/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/NemuPipeline/TermAndCock/midori_separate_ids/non-termite/{sp_name}_CO1.fasta", "w") as fout:
                fout.write(f">{id[0]}.{id[1]} [{sp_name}]\n{entry.seq}")
            