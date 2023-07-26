#there's an API, so fuck scrapper for now
import wikipediaapi
import pandas as pd
import numpy as np
import re

PATH_TO_META = '/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/NemuPipeline/AllInsects/midori_all_insects_mutspec/msMetaData.tsv'
PATH_TO_METAWINGS = '/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/NemuPipeline/AllInsects/midori_all_insects_mutspec/meta_wiki_wings.tsv'

meta_df = pd.read_csv(PATH_TO_META, sep='\t')
meta_df['Wings'] = np.nan

#insects = ['Periplaneta_americana', 'Papilio_helenus', 'afslalf', 'Heart and dart'] # test data
insects = meta_df['Species'].to_list()

wiki = wikipediaapi.Wikipedia('insecta_wings_db_bot (tr4sh.cannot@gmail.com)', 'en')
for sp in insects:
    try:
        page = wiki.page(sp)
        if page.exists():
            if re.search('wing$|wings|wing | wing', page.text) and len(page.title.split(' ')) > 1 and 'list' not in page.title.lower():
                print(f'{sp} - {page.title}')
                with open('/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/scripts/WikipediaScrapper/scrapper.log', 'a') as log:
                    log.write(f'{sp} - {page.title}\n')
                meta_df.loc[meta_df['Species'] == sp, 'Wings'] = '+'
            else:
                meta_df.loc[meta_df['Species'] == sp, 'Wings'] = '-'
    except:
        with open('/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/scripts/WikipediaScrapper/scrapper.log', 'a') as log:
                    log.write(f'{sp} - {page.title} - CONNECTION ERROR\n')
        continue

#meta_df['Wings'] = meta_df['Wings'].fillna('No_Data') 

meta_df.to_csv(PATH_TO_METAWINGS, sep='\t', index=False)



