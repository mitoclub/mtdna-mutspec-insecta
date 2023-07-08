from Bio import SeqIO
import numpy as np
import pandas as pd
from Bio.Data import CodonTable

translation = CodonTable.standard_dna_table.forward_table
for g in CodonTable.standard_dna_table.stop_codons:
    translation[g] = '_'
translation['TGA'] = 'W'
translation['ATA'] = 'M'
translation['AGA'] = '_'
translation['AGG'] = '_'

PATH_TO_GB = '/home/gabs/Documents/lab/TermitesAndCockroaches/MutSpec-Redone/interim/MIDORI/mergedAllGenesBlattodea.gb'
#PATH_TO_GB = '../../raw/blattodea_genbank/TermiteRefSeqs.gb'
PATH_TO_CODON_USAGE_TABLE = '/home/gabs/Documents/lab/TermitesAndCockroaches/MutSpec-Redone/interim/DescriptiveStat/codontable_midori_blattodea.csv'
#bc = SeqIO.parse(PATH_TO_GB, 'genbank')
item_table = ['Species_name','GenbankID', 'Taxonomy', 'Gene_name','Gene_start_end_and_trend', 'GeneID', 'Aminoacids_from_genbank',
             'Translated_aminoacids_by_Python', 'Not_standart_codons', 'Wrong_amino_num', 'Wrong_nucl_num','wrong_amino_%','Sequence','mtDNA_length',
             'nA','nT','nC','nG','nNA','%A','%T','%C','%G','%NA','wrong_amino_%','neutralA','neutralG','neutralC','neutralT', 'Neutral_count']
cl = list(CodonTable.standard_dna_table.forward_table.keys())
full_items = item_table + cl
items_manage = {}
for item in full_items:
    items_manage[item] = 0

btable = pd.DataFrame(columns=full_items)

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


for bc in SeqIO.parse(PATH_TO_GB, format='genbank'):
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
            items_manage['Sequence'] = i.location.extract(bc).seq
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
            a = list(i.location.extract(bc).seq)
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
            '''
            #blocked it, cuz it was causing some stupid problems I don't know how to solve. Just fuck it for now.
            #TODO: 
            Solve it
            if nuc != items_manage['Translated_aminoacids_by_Python'][nucl_c]:
            IndexError: string index out of range
            '''
            #for nuc in (items_manage['Aminoacids_from_genbank']):
            #    if nuc != items_manage['Translated_aminoacids_by_Python'][nucl_c]:
            #        proc_c += 1
            #    nucl_c += 1
            #items_manage['wrong_amino_%'] = proc_c * 100 / len(items_manage['Translated_aminoacids_by_Python'])
            #items_manage['Wrong_amino_num'] = proc_c
            #items_manage['Wrong_nucl_num'] = wrong_c



            proc_c = 0
            nucl_c = 0
            start_c = 0
            trnsl_amino = ''
            codons = []
            wrong_c = ''
            neuc = 0
            if ndc <= 13:
                btable = btable.append(items_manage, ignore_index=True)
                ndc += 1
            else:
                ndc = 1
            for k in items_manage:
                items_manage[k] = 0
btable.sort_values(['Species_name', 'Gene_name'], ascending=[True, True], inplace=True)
btable.to_csv(PATH_TO_CODON_USAGE_TABLE, sep=',')
