library(ggplot2)
library(ggpubr)
library(glue)

#set to true to use inverted skews (don't forget to first calculate them)
INVERT = TRUE
#set to true to use only CO1 skews (calculate them first)
JUST_CO1 = FALSE

TAXA <- 'Blattodea'

if (JUST_CO1 == TRUE){
derrived_skew <- read.csv(glue("/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/DescriptiveStat/midori_{TAXA}_CO1_skew.csv"))
}else{
derrived_skew <- read.csv(glue("/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/DescriptiveStat/midori_{TAXA}_skew.csv"))
} 
genes <- c("ND2", "CO1", "CO2", "A8", "A6",  "CO3", "ND3", "ND5", "ND4", "ND4L", "ND6", "Cytb", "ND1")
derrived_skew$Gene_name <- factor(derrived_skew$Gene_name, levels = genes)

#derrived_skew[derrived_skew$Gene_name == "ND4L", c("GAskew", "TCskew")] <- derrived_skew[derrived_skew$Gene_name == "ND4L", c("TCskew", "GAskew")]

if(INVERT == TRUE){
derrivedGAskew <- ggplot(data=subset(derrived_skew, !is.na(derrived_skew$Gene_name)), aes(x=Gene_name, y=GAskew, fill=Organism)) + ylab('(G-A)/(G+A)') +
  geom_boxplot() + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

derrivedTCskew <- ggplot(data=subset(derrived_skew, !is.na(derrived_skew$Gene_name)), aes(x=Gene_name, y=TCskew, fill=Organism)) + ylab('(T-C)/(T+C)') +
  geom_boxplot() + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


ggarrange(derrivedGAskew + ylim(-1.5, 1.5), derrivedTCskew + ylim(-1.5, 1.5),
          common.legend = TRUE,
          labels = c("A", "B"),
          ncol = 2, nrow = 1)
}else{
  derrivedAGskew <- ggplot(data=subset(derrived_skew, !is.na(derrived_skew$Gene_name)), aes(x=Gene_name, y=AGskew, fill=Organism)) + ylab('(A-G)/(A+G)') +
    geom_boxplot() + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  derrivedCTskew <- ggplot(data=subset(derrived_skew, !is.na(derrived_skew$Gene_name)), aes(x=Gene_name, y=CTskew, fill=Organism)) + ylab('(C-T)/(C+T)') +
    geom_boxplot() + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  
  ggarrange(derrivedCTskew + ylim(-1.5, 1.5), derrivedAGskew + ylim(-1.5, 1.5),
            common.legend = TRUE,
            labels = c("A", "B"),
            ncol = 2, nrow = 1)  
}
derrivedGAskew
derrivedTCskew

#######STATISTICS######
#TODO: Add INVERT trigger
cocks <- subset(subset(derrived_skew, !is.na(derrived_skew$Gene_name)), Organism == 'Cockroaches')
lower <- subset(subset(derrived_skew, !is.na(derrived_skew$Gene_name)), Organism == 'Termites w/o workers')
higher <- subset(subset(derrived_skew, !is.na(derrived_skew$Gene_name)), Organism == 'Termites w/ workers')
subsocial <- subset(subset(derrived_skew, !is.na(derrived_skew$Gene_name)), Organism == 'Sub-social Cryptocercus')

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

derrivedGAskew + stat_compare_means(method = 'anova')
derrivedGAskew + compare_means(GAskew ~ Organism, data = subset(derrived_skew, Gene_name == 'CO1'))

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
