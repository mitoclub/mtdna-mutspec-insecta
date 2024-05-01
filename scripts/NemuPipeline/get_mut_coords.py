import dask.dataframe as dd
import pandas as pd
import plotly.express as px
import plotly.io as pio

PATH_TO_MUTSPEC_META = '/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/NemuPipeline/AllInsects/midori_all_insects_mutspec/msMetaData.tsv'
PATH_TO_COORDS = '/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/map_insecta/test_coords.tsv'

coord_ddf = dd.read_csv(PATH_TO_COORDS, sep='\t', dtype={'gbifID' : 'string', 'datasetKey' : 'string','occurrenceID' : 'string','kingdom' : 'string','phylum' : 'string','class' : 'string','order' : 'string','family' : 'string','genus' : 'string','species' : 'string','infraspecificEpithet' : 'string','taxonRank' : 'string','scientificName' : 'string','verbatimScientificName' : 'string','verbatimScientificNameAuthorship' : 'string','countryCode' : 'string','locality' : 'string','stateProvince' : 'string','occurrenceStatus' : 'string','individualCount' : 'string','publishingOrgKey' : 'string','decimalLatitude' : 'float','decimalLongitude' : 'float','coordinateUncertaintyInMeters' : 'string','coordinatePrecision' : 'string','elevation' : 'string','elevationAccuracy' : 'string','depth' : 'string','depthAccuracy' : 'string','eventDate' : 'string','day' : 'string','month' : 'string','year' : 'string','taxonKey' : 'string','speciesKey' : 'string','basisOfRecord' : 'string','institutionCode' : 'string','collectionCode' : 'string','catalogNumber' : 'string','recordNumber' : 'string','identifiedBy' : 'string','dateIdentified' : 'string','license' : 'string','rightsHolder' : 'string','recordedBy' : 'string','typeStatus' : 'string','establishmentMeans' : 'string','lastInterpreted' : 'string','mediaType' : 'string','issue' : 'string'}, engine='python', encoding='utf8', quotechar='"', on_bad_lines='skip')
coord_ddf = coord_ddf[['species', 'order', 'decimalLatitude', 'decimalLongitude']]
coord_ddf = coord_ddf.dropna(subset=['species', 'decimalLatitude', 'decimalLongitude'])
meta_df = pd.read_csv(PATH_TO_MUTSPEC_META, sep='\t') 
meta_df.drop(meta_df[meta_df['Species'] == 'Stenopelmatus'].index, inplace=True) # this fucker is not species!
species = list(map(lambda sp: f'{sp.split("_")[0]} {sp.split("_")[1]}', meta_df['Species']))
print('###COMPUTE SPECIES LIST###')
df = coord_ddf[coord_ddf['species'].isin(species)].compute()
print('###SAVING DATA###')
df.to_csv('relevant_coords.tsv', sep='\t', index=False)
print('###PLOT###')
fig = px.scatter_geo(df, lat='decimalLatitude', lon='decimalLongitude', color='species',
                     hover_name='species')
fig.update_layout(showlegend=False)
pio.write_html(fig, file='muts_on_map.html')
