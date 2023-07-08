library(ggplot2)
library(ggpubr)

cock_term_skew <- read.csv("/home/gabs/Documents/lab/TermitesAndCockroaches/MutSpec-Redone/interim/DescriptiveStat/midori_blattodea_skew.csv")
genes <- c( "CO1", "CO2", "A8", "A6",  "CO3", "ND3", "ND4L", "ND4", "ND5", "Cytb")
cock_term_skew$Gene_name <- factor(cock_term_skew$Gene_name, levels = genes)

#74B72E - cock color, #466D1D -term color (didn't use)

# Inverting skew data for these genes, because they are located on a different strand! NOT SUITABLE FOR MIDORI!!!
#cock_term_skew <- transform(cock_term_skew, TCskew= ifelse(Gene_name == 'ND4L', GAskew, TCskew), GAskew = ifelse(Gene_name == 'ND4L', TCskew, GAskew))
#cock_term_skew <- transform(cock_term_skew, TCskew= ifelse(Gene_name == 'ND4', GAskew, TCskew), GAskew = ifelse(Gene_name == 'ND4', TCskew, GAskew))
#cock_term_skew <- transform(cock_term_skew, TCskew= ifelse(Gene_name == 'ND5', GAskew, TCskew), GAskew = ifelse(Gene_name == 'ND5', TCskew, GAskew))

#cock_term_skew[cock_term_skew$Gene_name == "ND4L", c("GAskew", "TCskew")] <- cock_term_skew[cock_term_skew$Gene_name == "ND4L", c("TCskew", "GAskew")]

######GA#########
cock_termGAskew <- ggplot(data=subset(cock_term_skew, !is.na(cock_term_skew$Gene_name)), aes(x=Gene_name, y=GAskew, fill=Organism)) + ylab('(G-A)/(G+A)') +
  geom_boxplot() + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

cock_termGAskew


#######TC#######
cock_termTCskew <- ggplot(data=subset(cock_term_skew, !is.na(cock_term_skew$Gene_name)), aes(x=Gene_name, y=TCskew, fill=Organism)) + ylab('(T-C)/(T+C)') +
  geom_boxplot() + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

cock_termTCskew

ggarrange(cock_termGAskew + ylim(-0.9, 0.7), cock_termTCskew + ylim(-0.9, 0.7),
          common.legend = TRUE,
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

#######Stg-Sac######
#cock_term_StgSac_skew <- ggplot(data=subset(cock_term_skew, !is.na(cock_term_skew$Gene_name)), aes(x=Gene_name, y=Stg.Sac, fill=Organism)) + 
#  geom_boxplot() + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#cock_term_StgSac_skew

#######STATISTICS######
cocks <- subset(subset(cock_term_skew, !is.na(cock_term_skew$Gene_name)), Organism == 'Cockroaches')
lower <- subset(subset(cock_term_skew, !is.na(cock_term_skew$Gene_name)), Organism == 'Termites w/o workers')
higher <- subset(subset(cock_term_skew, !is.na(cock_term_skew$Gene_name)), Organism == 'Termites w/ workers')
subsocial <- subset(subset(cock_term_skew, !is.na(cock_term_skew$Gene_name)), Organism == 'Sub-social Cryptocercus')

#TTEST
for (gene in genes)
{
  cocksVShigher <- t.test(subset(cocks, Gene_name == gene)$GAskew, subset(higher, Gene_name == gene)$GAskew)
  cocksVSlower <- t.test(subset(cocks, Gene_name == gene)$GAskew, subset(lower, Gene_name == gene)$GAskew)
  higherVSlower <- t.test(subset(higher, Gene_name == gene)$GAskew, subset(lower, Gene_name == gene)$GAskew)
  if (cocksVShigher$p.value < 0.05){
  print(sprintf('Cockroaches Vs Termites w/ workers: %s (<0.05) - %s', cocksVShigher$p.value, gene))
  }
  else
  {
    print(sprintf('Cockroaches Vs Termites w/ workers: %s - %s', cocksVShigher$p.value, gene))
  }
  if (cocksVSlower$p.value < 0.05){
  print(sprintf('Cockroaches Vs Termites w/o workers: %s (<0.05) - %s', cocksVSlower$p.value, gene))
  }
  else
  {
    print(sprintf('Cockroaches Vs Termites w/o workers: %s - %s', cocksVSlower$p.value, gene))
  }
  if (higherVSlower$p.value < 0.05){
  print(sprintf('Termites w/ workers Vs Termites w/o workers: %s (<0.05) - %s', higherVSlower$p.value, gene))
  }
  else
  {
    print(sprintf('Termites w/ workers Vs Termites w/o workers: %s - %s', higherVSlower$p.value, gene))
  }
  print('--------------------------------------------------------------------------')
}

#WILCOXONTEST
for (gene in genes)
{
  cocksVShigher <- wilcox.test(subset(cocks, Gene_name == gene)$GAskew, subset(higher, Gene_name == gene)$GAskew)
  cocksVSlower <- wilcox.test(subset(cocks, Gene_name == gene)$GAskew, subset(lower, Gene_name == gene)$GAskew)
  higherVSlower <- wilcox.test(subset(higher, Gene_name == gene)$GAskew, subset(lower, Gene_name == gene)$GAskew)
  if (cocksVShigher$p.value < 0.05){
    print(sprintf('Cockroaches Vs Termites w/ workers: %s (<0.05) - %s', cocksVShigher$p.value, gene))
  }
  else
  {
    print(sprintf('Cockroaches Vs Termites w/ workers: %s - %s', cocksVShigher$p.value, gene))
  }
  if (cocksVSlower$p.value < 0.05){
    print(sprintf('Cockroaches Vs Termites w/o workers: %s (<0.05) - %s', cocksVSlower$p.value, gene))
  }
  else
  {
    print(sprintf('Cockroaches Vs Termites w/o workers: %s - %s', cocksVSlower$p.value, gene))
  }
  if (higherVSlower$p.value < 0.05){
    print(sprintf('Termites w/ workers Vs Termites w/o workers: %s (<0.05) - %s', higherVSlower$p.value, gene))
  }
  else
  {
    print(sprintf('Termites w/ workers Vs Termites w/o workers: %s - %s', higherVSlower$p.value, gene))
  }
  print('--------------------------------------------------------------------------')
}
cocksVShigher
compare_pairs = list(c('Cockroaches', 'Termites w/ workers'),c('Cockroaches', 'Termites w/o workers'),
                     c('Termites w/o workers', 'Termites w/ workers'), c('Cockroaches', 'Sub-social Cryptocercus'))

cock_termGAskew + stat_compare_means(method = 'anova')
cock_termGAskew + compare_means(GAskew ~ Organism, data = subset(cock_term_skew, Gene_name == 'CO1'))

####















#####GTassymetry####

# I merged them manually and also manually added organism names
#gt_assymetry <- read.csv("/home/gab/Documents/lab/TermitesAndCockroaches/MutSpec-Redone/interim/DescriptiveStat/GTasymmetry_merged.csv")
#g_assymetry <- ggplot(data=gt_assymetry, aes(x=organism, y=med_c)) + ylab("Median G") +
#  geom_boxplot() + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 
#g_assymetry

#t_assymetry <- ggplot(data=gt_assymetry, aes(x=organism, y=med_a)) + ylab("Median T") +
#  geom_boxplot() + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 
#t_assymetry
