# mtdna-mutstep-insecta

# How everything works

1. filter_and_extract_prots.ipynb - allows to filter MIDORI NUCL UNIQ db to extract species with at least n number of unicue seqs present in the db (n is currently 5), also filters seqs by length (currently >= 600nuc). Finaly, it extracts relevant AA seqs from MIDORI for further work with pipeline
2. NeMuPipeline - does mutspec magic
3. aggregate_pipe_results.py - merge mutspec tables from different species into one table
4. obsnum_qc.ipynb - merge observed_mutations.tsv tables from different species (no saving) and plot relevant qc data (muts per sp, muts per branch, transitions fraction). For now, it automatically filters out the trees' branches with more than 6 mutations and then updates obsnum and therefore recalculates mutspecs in the merged table from the previous step (change paths accordingly).
5.  midori_parser_codonusage.ipynb - creates mega .gb file of specified taxa from MIDORI sequences
6.  codonusage_table.py - calculates nucleotide content (general and neutral) and number of codons. automaticaly reverse complements genes on negative strand (<ins>set manually</ins> to ND1, ND4L, ND4 and ND5), can drop species with wrong amino translations, which significantly cuts number of CockTer species for phylo PCA (from 32 to 11)
7.  calculate_skews.py - calculates skews, can invert codonusage data (including changing skew for GAskew to CTskew). Can only get stats for inverted data
8.  skews_figures.r - plots skews, can use inverted and not inverted data
9.  calculate_aa_fracs.py - calculates amino acid fractions or absolute values (set via trigger), can also reverse complement codons
10.  aa_fracs.r - plots amino acid fraction or absolute values
11.  plot_mutspec_and_nuc_content.ipynb - makes nice plots of mutspecs and nucl content for all the insects, exclusively CockVsTerm, Diptera (and Syrphidae) and Hymenoptera (badly implemented, be carefull). This code also can invert nucls for any included plot. Also makes mutspec_as_header_12internal.csv for cocks_phylo_stats.r
12.  cocks_phylo_stats.r - Makes PGLS, simple PCA and phyloPCA for skews or mutspecs. 
13.  aa_frac_comparison.ipynb plots AA fractions side by side (hue by organism type), <ins>in this case run it on reverse complemented codons,</ins> see different "butterfly" codonusage tables (pic's below). Also makes codonusage heatmap
14. cock_term_linalg.ipynb - calculate some linalg stats on ms with relation to nucleotide and aminoacid changes (WIP)

It is generally advised to run the scripts in that order as some of them might create files for the scripts bellow 

# Code related picks

### Codonusage tables for different strands (and strand explanation)
![wepik-butterfly-and-strands-20240510155245VQcE](https://github.com/user-attachments/assets/4ee8b431-e9e1-4f0d-a594-51ad16e33d2c)