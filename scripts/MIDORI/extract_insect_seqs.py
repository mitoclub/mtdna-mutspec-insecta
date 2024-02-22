import pandas as pd
from Bio import SeqIO


NUCL_PROT_CHECK = False #true if nucl, false if prots


PATH_TO_TABLE = '/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/MIDORI/midori_CO1_table.csv'
if NUCL_PROT_CHECK:
    PATH_TO_MIDORI = '/home/gabs/Documents/lab/TermitesAndCockroaches/MIDORI/MIDORI2_LONGEST_NUC_GB255_CO1_BLAST.fasta'
    PATH_TO_OUTPUT = '/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/MIDORI/all_insects_CO1seqs.fa'
else:
    PATH_TO_MIDORI = '/home/gabs/Documents/lab/TermitesAndCockroaches/MIDORI/CO1_longest.faa'
    PATH_TO_OUTPUT = '/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/MIDORI/all_insects_CO1seqs.faa'

table_df = pd.read_csv(PATH_TO_TABLE)

table_df = table_df[table_df.CO1 >= 10]


species_list = list(table_df['Species'])

insecta_df = []
taxonomy_df = []
genuses_df = []
seq_mapping = []
for entry in SeqIO.parse(PATH_TO_MIDORI, 'fasta'):
    taxonomy = entry.id.split(";")
    sp_name = taxonomy[7]
    id = taxonomy[0].split('.')
    fam = taxonomy[5]
    if sp_name in species_list:
        insecta_df.append(f">{id[0]}.{id[1]} [{sp_name}]\n{entry.seq}\n")
        taxonomy = reversed(re.sub(r'_[0-9]+', '', entry.id).split(";")[1:])
        taxonomy = ";".join(taxonomy)
        genuses_df.append(f'{taxonomy.split(";")[1]}\n')
        seq_mapping.append(f'{taxonomy.split(";")[0]}\t{taxonomy.split(";")[0]}\n')
        taxonomy_df.append(f'{taxonomy}\n')

with open('/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/MIDORI/all_insects_CO1.tax', 'w') as fout:
    for entry in taxonomy_df:
        fout.write(entry)

genuses_df = set(genuses_df)
with open('/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/MIDORI/all_insects_CO1.genus', 'w') as fout:
    for entry in genuses_df:
        fout.write(entry)

with open('/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/MIDORI/all_insects_CO1.sequence_mapping', 'w') as fout:
    for entry in seq_mapping:
        fout.write(entry)


with open(PATH_TO_OUTPUT, 'w') as fout:
    for entry in insecta_df:
        fout.write(entry)


        
