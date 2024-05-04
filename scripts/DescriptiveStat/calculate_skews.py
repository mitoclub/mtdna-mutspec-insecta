import pandas as pd
import numpy as np
import re


FAMILY = 'Blattodea'
PATH_TO_CODONTABLE = f'/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/DescriptiveStat/codontable_midori_{FAMILY}.csv'

def get_skew_df(codon_table):
    codon_table = pd.read_csv(codon_table)
    
    ###ADD IT TO SEPARATE SUBSOCIAL CRIPTOCERCUS
    if FAMILY == 'Blattodea':
        codon_table.loc[codon_table["Taxonomy"] == "['Eukaryota_2759', 'Arthropoda_6656', 'Insecta_50557', 'Blattodea_85823', 'Cryptocercidae_36982']",
                                        "Workers"] = 'Sub'
    
    ###Assigning organism type
    if FAMILY == 'Blattodea':
        workers = codon_table['Workers']
    #such a mess :(
    elif FAMILY == 'Diptera':
        codon_table['Organism'] = ""
        nematocera_families = ['Anisopodidae_52748','Bibionidae_52729','Cecidomyiidae_33406','Ceratopogonidae_41819','Chaoboridae_41811','Chironomidae_7149','Culicidae_7157','Keroplatidae_58254',
            'Limoniidae_43823','Mycetophilidae_29035','Psychodidae_7197','Ptychopteridae_79304','Sciaridae_7184','Simuliidae_7190','Tipulidae_41042']

        for fam in nematocera_families:
            for i in codon_table.index:
                real_fam = re.sub("['\[\]]",'',codon_table['Taxonomy'][i].split(',')[4].strip())
                if real_fam == fam:
                    codon_table.at[i, 'Organism'] = 'Nematocera'

        brachycera_families = ['Agromyzidae_127399','Anthomyiidae_30062','Asilidae_50673','Aulacigastridae_286480','Calliphoridae_7371','Chamaemyiidae_189958','Chloropidae_29032',
            'Clusiidae_286472','Conopidae_115263','Dolichopodidae_92558','Drosophilidae_7214','Dryomyzidae_169441','Empididae_92557','Ephydridae_48991',
            'Fanniidae_27471','Fergusoninidae_156410','Glossinidae_7392','Heleomyzidae_219548','Hippoboscidae_81710','Hybotidae_1446258','Lauxaniidae_189929',
            'Milichiidae_305559','Muscidae_7366','Mydidae_50677','Nemestrinidae_92615','Nycteribiidae_81707','Oestridae_7387','Opomyzidae_286476','Phoridae_36164',
            'Piophilidae_28629','Pipunculidae_43835','Platypezidae_43827','Platystomatidae_28632','Polleniidae_54279','Rhagionidae_92609','Sarcophagidae_7381',
            'Scathophagidae_43756','Sciomyzidae_169447','Sepsidae_137503','Sphaeroceridae_114620','Stratiomyidae_34687','Streblidae_81697','Syrphidae_34680',
            'Tachinidae_27474','Tephritidae_7211','Tabanidae_7205','Xylophagaidae_92613']
        for fam in brachycera_families:
            for i in codon_table.index:
                real_fam = re.sub("['\[\]]",'',codon_table['Taxonomy'][i].split(',')[4].strip())
                if real_fam == fam:
                    codon_table.at[i, 'Organism'] = 'Brachycera'
        #Removing species without set organism type
        codon_table = codon_table.drop(codon_table[codon_table['Organism'] == ""].index)
        organism = codon_table['Organism']

    #!!! NUCL VALUES ARE INVERTED FOR NOW
    GAskew = (codon_table['nC'] - codon_table['nT'])/(codon_table['nC'] + codon_table['nT'])
    TCskew = (codon_table['nA'] - codon_table['nG'])/(codon_table['nA'] + codon_table['nG'])
    #Stg_Sac = (codon_table['neutralT'] + codon_table['neutralG']) - (codon_table['neutralA'] + codon_table['neutralC'])
    IDs = codon_table['Species_name']
    genes = codon_table['Gene_name']
    
            
    skews = IDs.to_frame().merge(GAskew.rename('GAskew'), left_index=True, right_index=True)
    skews = skews.merge(genes.rename('Gene_name'), left_index=True, right_index=True)
    skews = skews.merge(TCskew.rename('TCskew'), left_index=True, right_index=True)
    #skews = skews.merge(Stg_Sac.rename('Stg-Sac'), left_index=True, right_index=True)
    if FAMILY == 'Blattodea':
        skews = skews.merge(workers.rename('Organism'), left_index=True, right_index=True)
    else:
        skews = skews.merge(organism.rename('Organism'), left_index=True, right_index=True)

    return (skews)

derrived_skew = get_skew_df(PATH_TO_CODONTABLE)
if FAMILY == 'Blattodea' or FAMILY == 'Diptera':
    derrived_skew = derrived_skew.sort_values(by=['Organism', 'Species_name'], ascending=True)
else:
    derrived_skew = derrived_skew.sort_values(by=['Species_name'], ascending=True)
if FAMILY == 'Blattodea':
    derrived_skew['Organism'] = derrived_skew['Organism'].map({1.0 : 'Termites w/ workers', 0.0 : 'Termites w/o workers', "Sub" : 'Sub-social Cryptocercus'}) ### LAST ONE SEPARATES SUBSOCIAL SPECIES
    derrived_skew['Organism'].fillna('Cockroaches', inplace=True)

derrived_skew.to_csv(f"/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/DescriptiveStat/midori_{FAMILY}_skew.csv", na_rep='NA')