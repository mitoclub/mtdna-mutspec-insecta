library(magrittr)
library(dplyr)

PATH_TO_CODON_USAGE = "/home/gab/Documents/lab/TermitesAndCockroaches/MutSpec-Redone/interim/DescriptiveStat/codonusage_table_cock.csv"
#PATH_TO_CODON_USAGE = "/home/gab/Documents/lab/TermitesAndCockroaches/MutSpec-Redone/interim/DescriptiveStat/codonusage_table_term.csv"
df_codon_usage_table <- read.csv(file = PATH_TO_CODON_USAGE)

vec_all = c('TTC','TTT','TCC','TCT','TAC','TAT','TGC','TGT',
            'TTA','TTG','TCA','TCG','TAA','TAG','TGA','TGG',
            'CTC','CTT','CCC','CCT','CAC','CAT','CGC','CGT',
            'CTA','CTG','CCA','CCG','CAA','CAG','CGA','CGG',
            'ATC','ATT','ACC','ACT','AAC','AAT','AGC','AGT',
            'ATA','ATG','ACA','ACG','AAA','AAG','AGA','AGG',
            'GTC','GTT','GCC','GCT','GAC','GAT','GGC','GGT',
            'GTA','GTG','GCA','GCG','GAA','GAG','GGA','GGG')

needed_codons = c('TTC','TCC','TAC','TGC',
                  'TTA','TCA','TAA','TGA',
                  'CTC','CCC','CAC','CGC',
                  'CTA','CCA','CAA','CGA',
                  'ATC','ACC','AAC','AGC',
                  'ATA','ACA','AAA','AGA',
                  'GTC','GCC','GAC','GGC',
                  'GTA','GCA','GAA','GGA')

df_codons = df_codon_usage_table %>% select(Species_name, Gene_name,all_of(vec_all))

sp_sum_gen = data.frame(unique(df_codons$Species_name))

for ( codon in vec_all)
{
  
  sum_of_codon = aggregate(df_codons[ ,codon], by = list(df_codons$Species_name), FUN = 'sum')[2]
  sp_sum_gen = cbind(sp_sum_gen, sum_of_codon)
  
}
names(sp_sum_gen) = c('Species', vec_all)


codon_norm = data.frame()

for (i in 1:nrow(sp_sum_gen))
{
  org_gen = sp_sum_gen[i,]
  org_gen = as.vector(org_gen)
  df_out= data.frame(sp_sum_gen[i,]$Species) 
  for (codon in seq(from = 2, to = 65, by = 2))
  {
    if (org_gen[1,codon] == 0) {df_out = cbind(df_out, 0)}
    else {
      norm_cod = org_gen[1,codon] / (org_gen[1,codon+1] + org_gen[1,codon])
      df_out = cbind(df_out, norm_cod)
    }
  }
  names(df_out) = c('Species', needed_codons)
  codon_norm = rbind(codon_norm,df_out)
}


names(codon_norm) = c('species_name', needed_codons)
codon_norm = codon_norm %>% select(-c('TAA','AGA'))


df_try = data.frame(unique(codon_norm))


final = data.frame()
for (org in 1:nrow(df_try))
{
  sp_r = df_try[org,]
  
  vec_of_c = sp_r %>% select(TTC, TCC, TAC, TGC, CTC, CCC, CAC, CGC,
                             ATC, ACC, AAC, AGC, GTC, GCC, GAC, GGC)
  vec_of_a = sp_r %>% select(TTA, TCA, TGA, CTA, CCA, CAA, CGA,
                             ATA, ACA, AAA, GTA, GCA, GAA, GGA)
  
  med_c = median(as.numeric(vec_of_c), na.rm = TRUE)
  med_a = median(as.numeric(vec_of_a), na.rm = TRUE)
  sp_out = data.frame(sp_r$species_name, med_c, med_a) 
  final = rbind(final,sp_out)
}
names(final) = c('species_name', 'med_c', 'med_a')
write.csv(final, file = 'GTasymmetry_cock.csv')
#write.csv(final, file = 'GTasymmetry_term.csv')