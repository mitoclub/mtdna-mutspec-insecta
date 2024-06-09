from Bio import SeqIO
import numpy as np
import pandas as pd
from Bio.Data import CodonTable
from tqdm import tqdm

import warnings 
warnings.filterwarnings('ignore') 

translation = CodonTable.unambiguous_dna_by_id[5].forward_table
for g in CodonTable.unambiguous_dna_by_id[5].stop_codons:
    translation[g] = '_'

DROP_WRONG_AMINO_GENES = False #true to remove genes with any amount of wrong amino. false to keep as is
TAXA = 'Orthoptera'
#set to True to use gb with CO1 only
JUST_CO1 = True
if JUST_CO1 == True:
    PATH_TO_GB = f'/mnt/data/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/MIDORI/mergedAllGenes{TAXA}_CO1.gb'
    PATH_TO_CODON_USAGE_TABLE = f'/mnt/data/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/DescriptiveStat/codontable_midori_{TAXA}_CO1.csv'
else:
    PATH_TO_GB = f'/mnt/data/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/MIDORI/mergedAllGenes{TAXA}.gb'
    PATH_TO_CODON_USAGE_TABLE = f'/mnt/data/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/DescriptiveStat/codontable_midori_{TAXA}.csv'
records_num = sum(1 for rec in SeqIO.parse(PATH_TO_GB, format='genbank').records)

cl = list(CodonTable.unambiguous_dna_by_id[5].forward_table.keys())
item_table = ['Species_name','GenbankID', 'Taxonomy', 'Gene_name','Gene_start_end_and_trend', 'GeneID', 'Aminoacids_from_genbank',
             'Translated_aminoacids_by_Python', 'Not_standart_codons', 'Wrong_amino_num', 'Wrong_nucl_num','wrong_amino_%','Sequence','mtDNA_length',
             'nA','nT','nC','nG','nNA','%A','%T','%C','%G','%NA','neutralA','neutralG','neutralC','neutralT', 'Neutral_count']
full_items = item_table + cl
btable = pd.DataFrame(columns=full_items)

for bc in tqdm(SeqIO.parse(PATH_TO_GB, format='genbank'), total=records_num, unit='sequences', desc='Sequences calculated'):
    item_table = ['Species_name','GenbankID', 'Taxonomy', 'Gene_name','Gene_start_end_and_trend', 'GeneID', 'Aminoacids_from_genbank',
                 'Translated_aminoacids_by_Python', 'Not_standart_codons', 'Wrong_amino_num', 'Wrong_nucl_num','wrong_amino_%','Sequence','mtDNA_length',
                 'nA','nT','nC','nG','nNA','%A','%T','%C','%G','%NA','wrong_amino_%','neutralA','neutralG','neutralC','neutralT', 'Neutral_count']

    full_items = item_table + cl
    items_manage = {}
    for item in full_items:
        items_manage[item] = 0

    counter = 0
    start_c = 0
    proc_c = 0
    codon_list = []
    triplet = ''
    codons = []
    trnsl_amino = ''
    kl = 0
    proc_c = 0
    nucl_c = 0
    wrong_c = ''
    ndc = 1
    neuc = 0
    items_manage = {}
    for item in full_items:
        items_manage[item] = 0
    ndc = 1
    for i in bc.features:
        if i.type == 'CDS':
            items_manage['Species_name'] = bc.annotations['organism']
            items_manage['GenbankID'] = bc.annotations['accessions']
            items_manage['Taxonomy'] = bc.annotations['taxonomy']
            if 'gene' not in i.qualifiers:
                items_manage['Gene_name'] = 'NA'
            else:
                items_manage['Gene_name'] = ''.join(i.qualifiers['gene'])
            items_manage['Gene_start_end_and_trend'] = i.location
            if 'db_xref' not in i.qualifiers:
                items_manage['GeneID'] = 'NA'
            else:
                items_manage['GeneID'] = i.qualifiers['db_xref'][-1]
            items_manage['Aminoacids_from_genbank'] = i.qualifiers['translation'][0]
            ###!!! REVERSE COMPLEMENTING NEGATIVE GENES, REMOVE IT IF THESE GENES AREN'T ON A "MINUS STRAND"
            if i.qualifiers['gene'] == ['ND1'] or i.qualifiers['gene'] == ['ND4'] or i.qualifiers['gene'] == ['ND4L'] or i.qualifiers['gene'] == ['ND5']:
                items_manage['Sequence'] = i.location.extract(bc).seq.reverse_complement()
            else:
                items_manage['Sequence'] = i.location.extract(bc).seq
            ###
            items_manage['mtDNA_length'] = len(i.location.extract(bc).seq)
            for nucl in i.location.extract(bc).seq:
                if nucl == 'A':
                    items_manage['nA'] += 1
                if nucl == 'T':
                    items_manage['nT'] += 1
                if nucl == 'C':
                    items_manage['nC'] += 1
                if nucl == 'G':
                    items_manage['nG'] += 1
                if nucl !='A' and nucl !='G' and nucl !='C'and nucl !='T':
                    items_manage['nNA'] += 1
            items_manage['%A'] = (items_manage['nA'])/(len(i.location.extract(bc).seq))
            items_manage['%T'] = (items_manage['nT'])/(len(i.location.extract(bc).seq))
            items_manage['%C'] = (items_manage['nC'])/(len(i.location.extract(bc).seq))
            items_manage['%G'] = (items_manage['nG'])/(len(i.location.extract(bc).seq))
            items_manage['%NA'] = (items_manage['nNA'])/(len(i.location.extract(bc).seq))
            ###!!! REVERSE COMPLEMENTING NEGATIVE GENES, REMOVE IT IF THESE GENES AREN'T ON A "MINUS STRAND"
            if i.qualifiers['gene'] == ['ND1'] or i.qualifiers['gene'] == ['ND4'] or i.qualifiers['gene'] == ['ND4L'] or i.qualifiers['gene'] == ['ND5']:
                a = list(i.location.extract(bc).seq.reverse_complement())
            else:
                a = i.location.extract(bc).seq
            ###
            for j in a:
                triplet = triplet + j
                counter += 1
                if counter == 3:
                    codon_list.append(triplet)
                    counter = 0
                    triplet = ''
            triplet = ''
            counter = 0
            codons.append(codon_list)
            codon_list = []
            a = 0
            for lst in codons:
                for trpl in lst:
                    start_c +=1
                    if trpl not in cl:
                        wrong_c += trpl
                        items_manage['Not_standart_codons'] += 1
                        trnsl_amino += 'X'
                    else:
                        if start_c == 1 and translation[trpl] != 'M':
                            trnsl_amino += 'M'
                        else:
                            trnsl_amino += translation[trpl]
                            items_manage[trpl] += 1
                        if trpl in ['CTT','CTC','CTA','CTG','GTT','GTC','GTA','GTG','TCT','TCA','TCG','TCC','CCT','CCG','CCA','CCC',
                                   'ACT','ACG','ACA','ACC','GCT','GCA','GCG','GCC','CGT','CGA','CGG','CGC','GGT','GGC','GGG','GGA']:
                            neuc += 1
                            if trpl[2] == 'A':
                                items_manage['neutralA'] += 1
                            if trpl[2] == 'G':
                                items_manage['neutralG'] += 1
                            if trpl[2] == 'C':
                                items_manage['neutralC'] += 1
                            if trpl[2] == 'T':
                                items_manage['neutralT'] += 1
            items_manage['Translated_aminoacids_by_Python'] = trnsl_amino
            items_manage['Neutral_count'] = neuc
            for nuc in (items_manage['Aminoacids_from_genbank']):
                if len(items_manage['Translated_aminoacids_by_Python']) < len(items_manage['Aminoacids_from_genbank']):
                    diff = len(items_manage['Aminoacids_from_genbank']) - len(items_manage['Translated_aminoacids_by_Python'])
                    items_manage['Translated_aminoacids_by_Python'] = items_manage['Translated_aminoacids_by_Python'] + ('_' * diff)
                if nuc != items_manage['Translated_aminoacids_by_Python'][nucl_c]:
                    proc_c += 1
                nucl_c += 1
            items_manage['wrong_amino_%'] = proc_c * 100 / len(items_manage['Translated_aminoacids_by_Python'])
            items_manage['Wrong_amino_num'] = proc_c
            items_manage['Wrong_nucl_num'] = wrong_c



            proc_c = 0
            nucl_c = 0
            start_c = 0
            trnsl_amino = ''
            codons = []
            wrong_c = ''
            neuc = 0

            if ndc <= 13:
                items_manage = pd.DataFrame.from_dict(items_manage, orient='index')
                items_manage = items_manage.transpose()
                btable = pd.concat([btable,items_manage], ignore_index=True)
                ndc += 1
            else:
                ndc = 1
            for k in items_manage:
                items_manage[k] = 0
btable.sort_values(['Species_name', 'Gene_name'], ascending=[True, True], inplace=True)

if TAXA == 'Blattodea':
    btable.insert(3, "Workers", pd.Series(dtype='int'))
    workers = {'Termitidae_46569': 1.0,
    'Anaplectidae_2163898': np.nan,
    'Blattidae_6974': np.nan,
    'Termopsidae_7501': 0.0,
    'Kalotermitidae_46562': 0.0,
    'Blaberidae_6979': np.nan,
    'Ectobiidae_1049651': np.nan,
    'Rhinotermitidae_36985': 0.0,
    'Cryptocercidae_36982': np.nan,
    'Corydiidae_30007': np.nan,
    'Lamproblattidae_1080998': np.nan,
    'Tryonicidae_1560744': np.nan,
    'Nocticolidae_85826':np.nan,
    'Mastotermitidae_37434': 1.0,
    'Hodotermitidae_70920': 1.0,
    'Serritermitidae_119664': 0.0}
    for i in btable.index:
        btable['Workers'][i] = workers[btable['Taxonomy'][i][4]]

if DROP_WRONG_AMINO_GENES == True:
    btable = btable.drop(btable[btable["Wrong_amino_num"] != 0].index)
    btable.drop(btable[btable["Wrong_amino_num"] != 0].index)
    gene_count = btable['Species_name'].value_counts().to_dict()
    for sp, count in gene_count.items():
        if count != 13:
            btable = btable.drop(btable[btable['Species_name'] == sp].index)
btable.to_csv(PATH_TO_CODON_USAGE_TABLE, sep=',')


