from Bio.Seq import Seq
import pandas as pd
import numpy as np
import re

# IMPORTANT! SET THIS VAL TO TRUE IF YOU WANT FRACTIONS, ELSE YOU'LL GET ABSOLUTE VALS
FRAC=True

# Set to true to reverse complement codons
REV_COMP = True
TAXA = 'Orthoptera'
#set to True to use gb with CO1 only
JUST_CO1 = True
if JUST_CO1 == True:
    PATH_TO_CODONTABLE = f'/mnt/data/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/DescriptiveStat/codontable_midori_{TAXA}_CO1.csv'
else:
    PATH_TO_CODONTABLE = f'/mnt/data/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/DescriptiveStat/codontable_midori_{TAXA}.csv'
codontable = pd.read_csv(PATH_TO_CODONTABLE, index_col=0)

if REV_COMP == True:
    #reverse complementing
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

if FRAC == True:
    total_codons_sum = codontable['TTT'] + codontable['TTC'] + codontable['TTA'] + codontable['TTG'] + codontable['TCT'] + codontable['TCC'] + codontable['TCA'] + codontable['TCG'] + codontable['TAT'] + codontable['TAC'] + codontable['TGT'] + codontable['TGC'] + codontable['TGG'] + codontable['CTT'] + codontable['CTC'] + codontable['CTA'] + codontable['CTG'] + codontable['CCT'] + codontable['CCC'] + codontable['CCA'] + codontable['CCG'] + codontable['CAT'] + codontable['CAC'] + codontable['CAA'] + codontable['CAG'] + codontable['CGT'] + codontable['CGC'] + codontable['CGA'] + codontable['CGG'] + codontable['ATT'] + codontable['ATC'] + codontable['ATA'] + codontable['ATG'] + codontable['ACT'] + codontable['ACC'] + codontable['ACA'] + codontable['ACG'] + codontable['AAT'] + codontable['AAC'] + codontable['AAA'] + codontable['AAG'] + codontable['AGT'] + codontable['AGC'] + codontable['AGA'] + codontable['AGG'] + codontable['GTT'] + codontable['GTC'] + codontable['GTA'] + codontable['GTG'] + codontable['GCT'] + codontable['GCC'] + codontable['GCA'] + codontable['GCG'] + codontable['GAT'] + codontable['GAC'] + codontable['GAA'] + codontable['GAG'] + codontable['GGT'] + codontable['GGC'] + codontable['GGA'] + codontable['GGG'] + codontable['TAA'] + codontable['TAG'] + codontable['TGA']
else:
    total_codons_sum = 1

phe_frac = (codontable['TTT'] + codontable['TTC'])/total_codons_sum
leu_frac = (codontable['TTA'] + codontable['TTG'])/total_codons_sum
leu2_frac = (codontable['CTT'] + codontable['CTC'] + codontable['CTA'] + codontable['CTG'])/total_codons_sum
ser_frac = (codontable['TCT'] + codontable['TCC'] + codontable['TCA'] + codontable['TCG'])/total_codons_sum
pro_frac = (codontable['CCT'] + codontable['CCC'] + codontable['CCA'] + codontable['CCG'])/total_codons_sum

cys_frac = (codontable['TGT'] + codontable['TGC'])/total_codons_sum
trp_frac = (codontable['TGA'] + codontable['TGG'])/total_codons_sum
tyr_frac = (codontable['TAT'] + codontable['TAC'])/total_codons_sum
arg_frac = (codontable['CGT'] + codontable['CGC'] + codontable['CGA'] + codontable['CGG'])/total_codons_sum
his_frac = (codontable['CAT'] + codontable['CAC'])/total_codons_sum
gln_frac = (codontable['CAA'] + codontable['CAG'])/total_codons_sum

gly_frac = (codontable['GGT'] + codontable['GGC'] + codontable['GGA'] + codontable['GGG'])/total_codons_sum
ser2_frac = (codontable['AGA'] + codontable['AGG'] + codontable['AGT'] + codontable['AGC'])/total_codons_sum
asp_frac = (codontable['GAT']+ codontable['GAC'])/total_codons_sum
glu_frac = (codontable['GAA'] + codontable['GAG'])/total_codons_sum
asn_frac = (codontable['AAT'] + codontable['AAC'])/total_codons_sum
lys_frac = (codontable['AAA'] + codontable['AAG'])/total_codons_sum

val_frac = (codontable['GTC'] + codontable['GTA'] + codontable['GTG'] + codontable['GTT'])/total_codons_sum
ile_frac = (codontable['ATT'] + codontable['ATC'])/total_codons_sum
met_frac = (codontable['ATA'] + codontable['ATG'])/total_codons_sum
ala_frac = (codontable['GCT'] + codontable['GCC'] + codontable['GCA'] + codontable['GCG'])/total_codons_sum
thr_frac = (codontable['ACT'] + codontable['ACC'] + codontable['ACA'] + codontable['ACG'])/total_codons_sum

if TAXA == 'Blattodea':
        workers = codontable['Workers']
    #such a mess :(
elif TAXA == 'Diptera':
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
elif TAXA == 'Apidae':
    codontable['Organism'] = ""
    social_genuses = ['Bombus', 'Apis', 'Exoneurella', 'Ceratina', 'Melipona', 'Partamona', 'Euglossa', 'Tetragonula']
    nonsocial_genuses = ['Xylocopa', 'Amegilla', 'Anthophora', 'Nomada', 'Eulaema', 'Eucera', 'Epeolus']
    for genus in social_genuses:
        for i in codontable.index:
            real_genus = codontable['Species_name'][i].split('_')[0]
            if real_genus == genus:
                codontable.at[i, 'Organism'] = 'Social'
    for genus in nonsocial_genuses:
        for i in codontable.index:
            real_genus = codontable['Species_name'][i].split('_')[0]
            if real_genus == genus:
                codontable.at[i, 'Organism'] = 'Nonsocial'
    #Removing species without set organism type
    codontable = codontable.drop(codontable[codontable['Organism'] == ""].index)
    organism = codontable['Organism']
elif TAXA == 'Hymenoptera':
    codontable['Organism'] = ""
    long_tongue_fams = ['Apidae_7458', 'Megachilidae_124286']
    for fam in long_tongue_fams:
        for i in codontable.index:
            real_fam = re.sub("['\[\]]",'',codontable['Taxonomy'][i].split(',')[4].strip())
            if real_fam == fam:
                codontable.at[i, 'Organism'] = 'Long tongued'
    short_tongue_fams = ['Andrenidae_48719', 'Halictidae_77572', 'Colletidae_156309']
    for fam in short_tongue_fams:
        for i in codontable.index:
            real_fam = re.sub("['\[\]]",'',codontable['Taxonomy'][i].split(',')[4].strip())
            if real_fam == fam:
                codontable.at[i, 'Organism'] = 'Short tongued'
    #Removing species without set organism type
    codontable = codontable.drop(codontable[codontable['Organism'] == ""].index)
    organism = codontable['Organism']
else:
    codontable['Organism'] = ""
    organism = codontable['Organism']

IDs = codontable['Species_name']
genes = codontable['Gene_name']

if TAXA == 'Blattodea':
    aa_fracs = IDs.to_frame().merge(workers.rename('Organism'), left_index=True, right_index=True)
else:
    aa_fracs = IDs.to_frame().merge(organism.rename('Organism'), left_index=True, right_index=True)

aa_fracs = aa_fracs.merge(genes.rename('Gene_name'), left_index=True, right_index=True)
aa_fracs = aa_fracs.merge(phe_frac.rename('Phe_frac'), left_index=True, right_index=True)
aa_fracs = aa_fracs.merge(leu_frac.rename('Leu_frac'), left_index=True, right_index=True)
aa_fracs = aa_fracs.merge(leu2_frac.rename('Leu2_frac'), left_index=True, right_index=True)
aa_fracs = aa_fracs.merge(ser_frac.rename('Ser_frac'), left_index=True, right_index=True)
aa_fracs = aa_fracs.merge(pro_frac.rename('Pro_frac'), left_index=True, right_index=True)
aa_fracs = aa_fracs.merge(cys_frac.rename('Cys_frac'), left_index=True, right_index=True)
aa_fracs = aa_fracs.merge(trp_frac.rename('Trp_frac'), left_index=True, right_index=True)
aa_fracs = aa_fracs.merge(tyr_frac.rename('Tyr_frac'), left_index=True, right_index=True)
aa_fracs = aa_fracs.merge(arg_frac.rename('Arg_frac'), left_index=True, right_index=True)
aa_fracs = aa_fracs.merge(his_frac.rename('His_frac'), left_index=True, right_index=True)
aa_fracs = aa_fracs.merge(gln_frac.rename('Gln_frac'), left_index=True, right_index=True)
aa_fracs = aa_fracs.merge(gly_frac.rename('Gly_frac'), left_index=True, right_index=True)
aa_fracs = aa_fracs.merge(ser2_frac.rename('Ser2_frac'), left_index=True, right_index=True)
aa_fracs = aa_fracs.merge(asp_frac.rename('Asp_frac'), left_index=True, right_index=True)
aa_fracs = aa_fracs.merge(glu_frac.rename('Glu_frac'), left_index=True, right_index=True)
aa_fracs = aa_fracs.merge(asn_frac.rename('Asn_frac'), left_index=True, right_index=True)
aa_fracs = aa_fracs.merge(lys_frac.rename('Lys_frac'), left_index=True, right_index=True)
aa_fracs = aa_fracs.merge(val_frac.rename('Val_frac'), left_index=True, right_index=True)
aa_fracs = aa_fracs.merge(ile_frac.rename('Ile_frac'), left_index=True, right_index=True)
aa_fracs = aa_fracs.merge(met_frac.rename('Met_frac'), left_index=True, right_index=True)
aa_fracs = aa_fracs.merge(ala_frac.rename('Ala_frac'), left_index=True, right_index=True)
aa_fracs = aa_fracs.merge(thr_frac.rename('Thr_frac'), left_index=True, right_index=True)


aa_fracs.sort_values(by=['Organism', 'Species_name'], ascending=True)

if TAXA == 'Blattodea':
    aa_fracs['Organism'] = aa_fracs['Organism'].map({1.0 : 'Termites', 0.0 : 'Termites'}) 
    aa_fracs['Organism'].fillna('Cockroaches', inplace=True)

if FRAC == True:
    if JUST_CO1 == True:
        aa_fracs.to_csv(f"/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/DescriptiveStat/midori_{TAXA}_CO1_aa_fracs.tsv", sep='\t', na_rep='NA', index=False)
    else:
        aa_fracs.to_csv(f"/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/DescriptiveStat/midori_{TAXA}_aa_fracs.tsv", sep='\t', na_rep='NA', index=False)
else:
    aa_fracs.to_csv(f"/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/DescriptiveStat/midori_{TAXA}_aa_abs_vals.tsv", sep='\t', na_rep='NA', index=False)