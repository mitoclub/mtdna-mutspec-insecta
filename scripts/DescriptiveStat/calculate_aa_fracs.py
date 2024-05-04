from Bio.Seq import Seq
import pandas as pd
import numpy as np
import re


FAMILY = 'Blattodea'
PATH_TO_CODONTABLE = f'/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/DescriptiveStat/codontable_midori_{FAMILY}.csv'
codontable = pd.read_csv(PATH_TO_CODONTABLE, index_col=0)
#reverse complementig. probably unreasonable or is it?
codons = ['TTT', 'TTC', 'TTA', 'TTG', 'TCT', 'TCC',
       'TCA', 'TCG', 'TAT', 'TAC', 'TGT', 'TGC', 'TGG', 'CTT', 'CTC', 'CTA',
       'CTG', 'CCT', 'CCC', 'CCA', 'CCG', 'CAT', 'CAC', 'CAA', 'CAG', 'CGT',
       'CGC', 'CGA', 'CGG', 'ATT', 'ATC', 'ATA', 'ATG', 'ACT', 'ACC', 'ACA',
       'ACG', 'AAT', 'AAC', 'AAA', 'AAG', 'AGT', 'AGC', 'AGA', 'AGG', 'GTT',
       'GTC', 'GTA', 'GTG', 'GCT', 'GCC', 'GCA', 'GCG', 'GAT', 'GAC', 'GAA',
       'GAG', 'GGT', 'GGC', 'GGA', 'GGG', 'TAA', 'TAG', 'TGA']
revcomp_codons = {}
for cod in codons:
    revcomp = Seq(cod).reverse_complement()
    revcomp_codons[cod] = str(revcomp)
codontable.rename(columns=revcomp_codons, inplace=True)
total_codons_sum = codontable['TTT'] + codontable['TTC'] + codontable['TTA'] + codontable['TTG'] + codontable['TCT'] + codontable['TCC'] + codontable['TCA'] + codontable['TCG'] + codontable['TAT'] + codontable['TAC'] + codontable['TGT'] + codontable['TGC'] + codontable['TGG'] + codontable['CTT'] + codontable['CTC'] + codontable['CTA'] + codontable['CTG'] + codontable['CCT'] + codontable['CCC'] + codontable['CCA'] + codontable['CCG'] + codontable['CAT'] + codontable['CAC'] + codontable['CAA'] + codontable['CAG'] + codontable['CGT'] + codontable['CGC'] + codontable['CGA'] + codontable['CGG'] + codontable['ATT'] + codontable['ATC'] + codontable['ATA'] + codontable['ATG'] + codontable['ACT'] + codontable['ACC'] + codontable['ACA'] + codontable['ACG'] + codontable['AAT'] + codontable['AAC'] + codontable['AAA'] + codontable['AAG'] + codontable['AGT'] + codontable['AGC'] + codontable['AGA'] + codontable['AGG'] + codontable['GTT'] + codontable['GTC'] + codontable['GTA'] + codontable['GTG'] + codontable['GCT'] + codontable['GCC'] + codontable['GCA'] + codontable['GCG'] + codontable['GAT'] + codontable['GAC'] + codontable['GAA'] + codontable['GAG'] + codontable['GGT'] + codontable['GGC'] + codontable['GGA'] + codontable['GGG'] + codontable['TAA'] + codontable['TAG'] + codontable['TGA']
phe_frac = (codontable['TTT'] + codontable['TTC'])/total_codons_sum
leu_frac = (codontable['TTA'] + codontable['TTG'])/total_codons_sum
pro_frac = (codontable['CCT'] + codontable['CCC'] + codontable['CCA'] + codontable['CCG'])/total_codons_sum

if FAMILY == 'Blattodea':
        workers = codontable['Workers']
    #such a mess :(
elif FAMILY == 'Diptera':
    codontable['Organism'] = ""
    nematocera_families = ['Anisopodidae_52748','Bibionidae_52729','Cecidomyiidae_33406','Ceratopogonidae_41819','Chaoboridae_41811','Chironomidae_7149','Culicidae_7157','Keroplatidae_58254',
        'Limoniidae_43823','Mycetophilidae_29035','Psychodidae_7197','Ptychopteridae_79304','Sciaridae_7184','Simuliidae_7190','Tipulidae_41042']
    
    for fam in nematocera_families:
        for i in codontable.index:
            real_fam = re.sub("['\[\]]",'',codontable['Taxonomy'][i].split(',')[4].strip())
            if real_fam == fam:
                codontable.at[i, 'Organism'] = 'Nematocera'

    brachycera_families = ['Agromyzidae_127399','Anthomyiidae_30062','Asilidae_50673','Aulacigastridae_286480','Calliphoridae_7371','Chamaemyiidae_189958','Chloropidae_29032',
        'Clusiidae_286472','Conopidae_115263','Dolichopodidae_92558','Drosophilidae_7214','Dryomyzidae_169441','Empididae_92557','Ephydridae_48991',
        'Fanniidae_27471','Fergusoninidae_156410','Glossinidae_7392','Heleomyzidae_219548','Hippoboscidae_81710','Hybotidae_1446258','Lauxaniidae_189929',
        'Milichiidae_305559','Muscidae_7366','Mydidae_50677','Nemestrinidae_92615','Nycteribiidae_81707','Oestridae_7387','Opomyzidae_286476','Phoridae_36164',
        'Piophilidae_28629','Pipunculidae_43835','Platypezidae_43827','Platystomatidae_28632','Polleniidae_54279','Rhagionidae_92609','Sarcophagidae_7381',
        'Scathophagidae_43756','Sciomyzidae_169447','Sepsidae_137503','Sphaeroceridae_114620','Stratiomyidae_34687','Streblidae_81697','Syrphidae_34680',
        'Tachinidae_27474','Tephritidae_7211','Tabanidae_7205','Xylophagaidae_92613']
    
    for fam in brachycera_families:
        for i in codontable.index:
            real_fam = re.sub("['\[\]]",'',codontable['Taxonomy'][i].split(',')[4].strip())
            if real_fam == fam:
                codontable.at[i, 'Organism'] = 'Brachycera'

    #Removing species without set organism type
    codontable = codontable.drop(codontable[codontable['Organism'] == ""].index)
    organism = codontable['Organism']

IDs = codontable['Species_name']
genes = codontable['Gene_name']
aa_fracs = IDs.to_frame().merge(phe_frac.rename('Phe_frac'), left_index=True, right_index=True)
aa_fracs = aa_fracs.merge(leu_frac.rename('Leu_frac'), left_index=True, right_index=True)
aa_fracs = aa_fracs.merge(pro_frac.rename('Pro_frac'), left_index=True, right_index=True)
aa_fracs = aa_fracs.merge(genes.rename('Gene_name'), left_index=True, right_index=True)

if FAMILY == 'Blattodea':
    aa_fracs = aa_fracs.merge(workers.rename('Organism'), left_index=True, right_index=True)
else:
    aa_fracs = aa_fracs.merge(organism.rename('Organism'), left_index=True, right_index=True)

aa_fracs.sort_values(by=['Organism', 'Species_name'], ascending=True)

if FAMILY == 'Blattodea':
    aa_fracs['Organism'] = aa_fracs['Organism'].map({1.0 : 'Termites', 0.0 : 'Termites'}) 
    aa_fracs['Organism'].fillna('Cockroaches', inplace=True)

aa_fracs.to_csv(f"/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/DescriptiveStat/midori_{FAMILY}_aa_fracs.csv", na_rep='NA')