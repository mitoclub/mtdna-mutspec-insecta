# mtdna-mutstep-insecta

# How everything works

&nbsp; 0. NeMuPipeline - does mutspec magic

&nbsp; 0.5. aggregate_pipe_results.py - merge mutspec tables from different species into one table

1.  midori_parser_codonusage.ipynb - creates mega .gb file of specified taxa from MIDORI sequences
2.  codonusage_table.py - calculates nucleotide content (general and neutral) and number of codons. automaticaly reverse complements genes on negative strand (<ins>set manually</ins> to ND1, ND4L, ND4 and ND5), can drop species with wrong amino translations, which significantly cuts number of CockTer species for phylo PCA (from 32 to 11)
3.  calculate_skews.py - calculates skews, can invert codonusage data (including changing skew for GAskew to CTskew). Can only get stats for inverted data
4.  skews_figures.r - plots skews, can use inverted and not inverted data
5.  calculate_aa_fracs.py - calculates amino acid fractions or absolute values (set via trigger), can also reverse complement codons
6.  aa_fracs.r - plots amino acid fraction or absolute values
7.  plot_mutspec_and_nuc_content.ipynb - makes nice plots of mutspecs and nucl content for all the insects, exclusively CockVsTerm, Diptera (and Syrphidae) and Hymenoptera (badly implemented, be carefull). This code also can invert nucls for any included plot. Also makes mutspec_as_header_12internal.csv for cocks_phylo_stats.r
8.  cocks_phylo_stats.r - Makes PGLS, simple PCA and phyloPCA for skews or mutspecs. 
9.  aa_frac_comparison.ipynb plots AA fractions side by side (hue by organism type), <ins>in this case run it on reverse complemented codons,</ins> see different "butterfly" codonusage tables (pick's below). Also makes codonusage heatmap

# Code related picks

### Codonusage tables for different strands (and strand explanation)
![wepik-butterfly-and-strands-20240510155245VQcE](https://github.com/user-attachments/assets/4ee8b431-e9e1-4f0d-a594-51ad16e33d2c)
