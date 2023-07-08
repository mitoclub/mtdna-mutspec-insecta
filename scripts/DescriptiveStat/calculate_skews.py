import pandas as pd

PATH_TO_BLATTODEA = '/home/gabs/Documents/lab/TermitesAndCockroaches/MutSpec-Redone/interim/DescriptiveStat/codontable_midori_blattodea.csv'

def get_skew_df(codon_table):
    codon_table = pd.read_csv(codon_table)
    
    ###ADD IT TO SEPARATE SUBSOCIAL CRIPTOCERCUS
    codon_table.loc[codon_table["Taxonomy"] == "['Eukaryota_2759', 'Arthropoda_6656', 'Insecta_50557', 'Blattodea_85823', 'Cryptocercidae_36982']",
                                    "Workers"] = 'Sub'
    ###

    GAskew = (codon_table['nG'] - codon_table['nA'])/(codon_table['nG'] + codon_table['nA'])
    TCskew = (codon_table['nT'] - codon_table['nC'])/(codon_table['nT'] + codon_table['nC'])
    Stg_Sac = (codon_table['neutralT'] + codon_table['neutralG']) - (codon_table['neutralA'] + codon_table['neutralC'])
    IDs = codon_table['Species_name']
    genes = codon_table['Gene_name']
    workers = codon_table['Workers']

    skews = IDs.to_frame().merge(GAskew.rename('GAskew'), left_index=True, right_index=True)
    skews = skews.merge(genes.rename('Gene_name'), left_index=True, right_index=True)
    skews = skews.merge(TCskew.rename('TCskew'), left_index=True, right_index=True)
    skews = skews.merge(Stg_Sac.rename('Stg-Sac'), left_index=True, right_index=True)
    skews = skews.merge(workers.rename('Organism'), left_index=True, right_index=True)

    return (skews)



blattodea_skew = get_skew_df(PATH_TO_BLATTODEA)
blattodea_skew = blattodea_skew.sort_values(by=['Organism', 'Species_name'], ascending=True)

blattodea_skew['Organism'] = blattodea_skew['Organism'].map({1.0 : 'Termites w/ workers', 0.0 : 'Termites w/o workers', "Sub" : 'Sub-social Cryptocercus'}) ### LAST ONE SEPARATES SUBSOCIAL SPECIES
blattodea_skew['Organism'].fillna('Cockroaches', inplace=True)


blattodea_skew.to_csv("/home/gabs/Documents/lab/TermitesAndCockroaches/MutSpec-Redone/interim/DescriptiveStat/midori_blattodea_skew.csv", na_rep='NA')


#cock_skew.to_csv("../../interim/DescriptiveStat/cock_skew.csv", na_rep='NA')
#term_skew.to_csv("../../interim/DescriptiveStat/term_skew.csv", na_rep='NA')
