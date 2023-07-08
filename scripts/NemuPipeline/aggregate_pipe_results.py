import pandas as pd
import numpy as np
from ete3 import Tree
import os

'''
Task:
- aggregate all info into one (or several) table
- get number of collected mutations from observed_mutations_iqtree.tsv file (column ProbaFull, sum values in this column to get estimation of nmuts, 
  don't include values(probabilities) less then 0.3 OR you can get ObsNum columns from ms12syn_internal_iqtree.tsv and ms12syn_iqtree.tsv - these files 
  contains all info about used observed muts
- get spectra values from ms12syn_internal_iqtree.tsv and ms192syn_internal_iqtree.tsv
- get number of tree nodes from iqtree_anc_tree.nwk or iqtree_anc.state

Final columns in the table:
- spectra values (12 + 192 (if exists) columns)
- number of collected mutations (2 columns: mutations from all branches (ms12syn_iqtree.tsv) and from non-terminal (internal)
- number of seqs in alignment (leaves in tree)
- number of nodes in tree (optional)
'''

PATH_TO_FOLDER = '/home/gabs/Documents/lab/TermitesAndCockroaches/MutSpec-Redone/results/Dolphin/test/'
PATH_TO_SAVE = '/home/gabs/Documents/lab/TermitesAndCockroaches/MutSpec-Redone/results/Dolphin/midori_all_insects_mutspecs/'

filler_muts12 = ['A>C', 'A>G', 'A>T', 'C>A', 'C>G', 'C>T', 'G>A', 'G>C', 'G>T', 'T>A', 'T>C', 'T>G']
filler_muts192 = ['A[A>C]A', 'A[A>C]C', 'A[A>C]G', 'A[A>C]T', 'A[A>G]A', 'A[A>G]C', 'A[A>G]G', 'A[A>G]T', 'A[A>T]A', 'A[A>T]C', 'A[A>T]G', 'A[A>T]T', 'A[C>A]A', 'A[C>A]C', 'A[C>A]G', 'A[C>A]T', 'A[C>G]A', 'A[C>G]C', 'A[C>G]G', 'A[C>G]T', 'A[C>T]A', 'A[C>T]C', 'A[C>T]G', 'A[C>T]T', 'A[G>A]A', 'A[G>A]C', 'A[G>A]G', 'A[G>A]T', 'A[G>C]A', 'A[G>C]C', 'A[G>C]G', 'A[G>C]T', 'A[G>T]A', 'A[G>T]C', 'A[G>T]G', 'A[G>T]T', 'A[T>A]A', 'A[T>A]C', 'A[T>A]G', 'A[T>A]T', 'A[T>C]A', 'A[T>C]C', 'A[T>C]G', 'A[T>C]T', 'A[T>G]A', 'A[T>G]C', 'A[T>G]G', 'A[T>G]T', 'C[A>C]A', 'C[A>C]C', 'C[A>C]G', 'C[A>C]T', 'C[A>G]A', 'C[A>G]C', 'C[A>G]G', 'C[A>G]T', 'C[A>T]A', 'C[A>T]C', 'C[A>T]G', 'C[A>T]T', 'C[C>A]A', 'C[C>A]C', 'C[C>A]G', 'C[C>A]T', 'C[C>G]A', 'C[C>G]C', 'C[C>G]G', 'C[C>G]T', 'C[C>T]A', 'C[C>T]C', 'C[C>T]G', 'C[C>T]T', 'C[G>A]A', 'C[G>A]C', 'C[G>A]G', 'C[G>A]T', 'C[G>C]A', 'C[G>C]C', 'C[G>C]G', 'C[G>C]T', 'C[G>T]A', 'C[G>T]C', 'C[G>T]G', 'C[G>T]T', 'C[T>A]A', 'C[T>A]C', 'C[T>A]G', 'C[T>A]T', 'C[T>C]A', 'C[T>C]C', 'C[T>C]G', 'C[T>C]T', 'C[T>G]A', 'C[T>G]C', 'C[T>G]G', 'C[T>G]T', 'G[A>C]A', 'G[A>C]C', 'G[A>C]G', 'G[A>C]T', 'G[A>G]A', 'G[A>G]C', 'G[A>G]G', 'G[A>G]T', 'G[A>T]A', 'G[A>T]C', 'G[A>T]G', 'G[A>T]T', 'G[C>A]A', 'G[C>A]C', 'G[C>A]G', 'G[C>A]T', 'G[C>G]A', 'G[C>G]C', 'G[C>G]G', 'G[C>G]T', 'G[C>T]A', 'G[C>T]C', 'G[C>T]G', 'G[C>T]T', 'G[G>A]A', 'G[G>A]C', 'G[G>A]G', 'G[G>A]T', 'G[G>C]A', 'G[G>C]C', 'G[G>C]G', 'G[G>C]T', 'G[G>T]A', 'G[G>T]C', 'G[G>T]G', 'G[G>T]T', 'G[T>A]A', 'G[T>A]C', 'G[T>A]G', 'G[T>A]T', 'G[T>C]A', 'G[T>C]C', 'G[T>C]G', 'G[T>C]T', 'G[T>G]A', 'G[T>G]C', 'G[T>G]G', 'G[T>G]T', 'T[A>C]A', 'T[A>C]C', 'T[A>C]G', 'T[A>C]T', 'T[A>G]A', 'T[A>G]C', 'T[A>G]G', 'T[A>G]T', 'T[A>T]A', 'T[A>T]C', 'T[A>T]G', 'T[A>T]T', 'T[C>A]A', 'T[C>A]C', 'T[C>A]G', 'T[C>A]T', 'T[C>G]A', 'T[C>G]C', 'T[C>G]G', 'T[C>G]T', 'T[C>T]A', 'T[C>T]C', 'T[C>T]G', 'T[C>T]T', 'T[G>A]A', 'T[G>A]C', 'T[G>A]G', 'T[G>A]T', 'T[G>C]A', 'T[G>C]C', 'T[G>C]G', 'T[G>C]T', 'T[G>T]A', 'T[G>T]C', 'T[G>T]G', 'T[G>T]T', 'T[T>A]A', 'T[T>A]C', 'T[T>A]G', 'T[T>A]T', 'T[T>C]A', 'T[T>C]C', 'T[T>C]G', 'T[T>C]T', 'T[T>G]A', 'T[T>G]C', 'T[T>G]G', 'T[T>G]T']

df_meta = []
df_all192 = []
df_all12 = []
for sp in os.listdir(PATH_TO_FOLDER):
    if '.' in sp:
        continue
    print(sp)
    PATH_TO_ms12syn = f'{PATH_TO_FOLDER}{sp}/mutspec_tables/ms12syn_iqtree.tsv'
    PATH_TO_ms12syn_inter = f'{PATH_TO_FOLDER}{sp}/mutspec_tables/ms12syn_internal_iqtree.tsv'
    PATH_TO_ms192syn =f'{PATH_TO_FOLDER}{sp}/mutspec_tables/ms192syn_iqtree.tsv'
    PATH_TO_ms192syn_inter =f'{PATH_TO_FOLDER}{sp}/mutspec_tables/ms192syn_internal_iqtree.tsv'
    
    if os.path.isfile(PATH_TO_ms12syn):
        df12 = pd.read_csv(PATH_TO_ms12syn, sep='\t')
    else:
        df12 = pd.DataFrame({'Mut' : filler_muts12, 'ObsNum' : np.nan, 'ExpNum' : np.nan, 'MutSpec' : np.nan})
    if os.path.isfile(PATH_TO_ms12syn_inter):
        df_inter12 = pd.read_csv(PATH_TO_ms12syn_inter, sep='\t')
    else:
        df_inter12 = pd.DataFrame({'Mut' : filler_muts12, 'ObsNum' : np.nan, 'ExpNum' : np.nan, 'MutSpec' : np.nan})
    if os.path.isfile(PATH_TO_ms192syn):
        df192 = pd.read_csv(PATH_TO_ms192syn, sep='\t')
    else:
        df192 = pd.DataFrame({'Mut' : filler_muts192, 'ObsNum' : np.nan, 'ExpNum' : np.nan, 'MutSpec' : np.nan})
    if os.path.isfile(PATH_TO_ms192syn_inter):
        df_inter192 = pd.read_csv(PATH_TO_ms192syn_inter, sep='\t')
    else:
        df_inter192 = pd.DataFrame({'Mut' : filler_muts192, 'ObsNum' : np.nan, 'ExpNum' : np.nan, 'MutSpec' : np.nan})

    df12['Species'] = sp.split('__')[1]
    df_inter12['Species'] = sp.split('__')[1]
    df192['Species'] = sp.split('__')[1]
    df_inter192['Species'] = sp.split('__')[1]

    df12['ObsNum_inter'] = df_inter12['ObsNum']
    df12['ExpNum_inter'] = df_inter12['ExpNum']
    df12['MutSpec12inter'] = df_inter12['MutSpec']
    df192['ObsNum_inter'] = df_inter192['ObsNum']
    df192['ExpNum_inter'] = df_inter192['ExpNum']
    df192['MutSpec192inter'] = df_inter192['MutSpec']
    df_all12.append(df12[['Species','Mut', 'ObsNum', 'ObsNum_inter', 'ExpNum', 'ExpNum_inter', 'MutSpec12inter']])
    df_all192.append(df192[['Species','Mut', 'ObsNum', 'ObsNum_inter', 'ExpNum', 'ExpNum_inter', 'MutSpec192inter']])
    
    if os.path.isfile(f'{PATH_TO_FOLDER}{sp}/IQTREE/iqtree_anc_tree.nwk'):
        tree = Tree(f'{PATH_TO_FOLDER}{sp}/IQTREE/iqtree_anc_tree.nwk', format=1)
        num_of_nodes = len(tree) - 1 # -1 to remove the outgroup
    else:
        num_of_nodes = np.nan
    with open (f'{PATH_TO_FOLDER}{sp}/{sp}.fasta', 'r') as fasta_seqs:
        taxonomy = fasta_seqs.readline().split('###')[1]
    taxonomy = taxonomy.split(';')
    df_meta.append(pd.DataFrame({'Species' : [sp.split('__')[1]], 
                                'Class' : [taxonomy[1]], 
                                'Order' : [taxonomy[2]], 
                                'Family' : [taxonomy[3]], 
                                'Genus' : [taxonomy[4]],
                                'Nodes_in_tree' : [num_of_nodes]}))


df_all12 = pd.concat(df_all12, ignore_index=True)
df_all192 = pd.concat(df_all192, ignore_index=True)
df_meta = pd.concat(df_meta, ignore_index=True)


print(df_all12)
#df_all12.to_csv(f'{PATH_TO_SAVE}mutspec12.tsv', sep='\t')
#df_all192.to_csv(f'{PATH_TO_SAVE}mutspec192.tsv', sep='\t')
#df_meta.to_csv(f'{PATH_TO_SAVE}mutspecMetaData.tsv', sep='\t')