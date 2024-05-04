library(ggplot2)
library(ggpubr)
library(glue)

FAMILY <- 'Blattodea'

aa_fracs <- read.csv(glue("/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/DescriptiveStat/midori_{FAMILY}_aa_fracs.csv"))
genes <- c( "CO1", "CO2", "A8", "A6",  "CO3", "ND3", "ND4L", "ND4", "ND5", "Cytb")
aa_fracs$Gene_name <- factor(aa_fracs$Gene_name, levels = genes)

phe_fracs <- ggplot(data=subset(aa_fracs, !is.na(aa_fracs$Gene_name)), aes(x=Gene_name, y=Phe_frac, fill=Organism)) + ylab('Phe fraction') +
  geom_boxplot() + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
phe_fracs

leu_frac <- ggplot(data=subset(aa_fracs, !is.na(aa_fracs$Gene_name)), aes(x=Gene_name, y=Leu_frac, fill=Organism)) + ylab('Leu fraction') +
  geom_boxplot() + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
leu_frac

pro_frac <- ggplot(data=subset(aa_fracs, !is.na(aa_fracs$Gene_name)), aes(x=Gene_name, y=Pro_frac, fill=Organism)) + ylab('Pro fraction') +
  geom_boxplot() + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
pro_frac


ggarrange(phe_fracs + ylim(0, 0.2), leu_frac + ylim(0, 0.2), 
          pro_frac + ylim(0, 0.2),
          common.legend = TRUE,
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)
#ylim(-0.05, 0.2)
#ylim(-0.05, 130)

