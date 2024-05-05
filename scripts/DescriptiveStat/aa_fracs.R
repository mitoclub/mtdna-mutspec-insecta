library(ggplot2)
library(ggpubr)
library(glue)
library(readr)

FAMILY <- 'Blattodea'
FRACS <- TRUE
if(FRACS == TRUE){
aa_fracs <- read_tsv(glue("/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/DescriptiveStat/midori_{FAMILY}_aa_fracs.tsv"))
}else{
  aa_fracs <- read_tsv(glue("/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/DescriptiveStat/midori_{FAMILY}_aa_abs_vals.tsv"))  
}
genes <- c( "CO1", "CO2", "A8", "A6",  "CO3", "ND3", "ND4L", "ND4", "ND5", "Cytb")
aa_fracs$Gene_name <- factor(aa_fracs$Gene_name, levels = genes)

Phe_fracs <- ggplot(data=subset(aa_fracs, !is.na(aa_fracs$Gene_name)), aes(x=Gene_name, y=Phe_frac, fill=Organism)) + ylab('Phe fraction') +
  geom_boxplot() + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

Leu_frac <- ggplot(data=subset(aa_fracs, !is.na(aa_fracs$Gene_name)), aes(x=Gene_name, y=Leu_frac, fill=Organism)) + ylab('Leu fraction') +
  geom_boxplot() + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

Leu2_frac <- ggplot(data=subset(aa_fracs, !is.na(aa_fracs$Gene_name)), aes(x=Gene_name, y=Leu2_frac, fill=Organism)) + ylab('Leu2 fraction') +
  geom_boxplot() + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

Ser_frac <- ggplot(data=subset(aa_fracs, !is.na(aa_fracs$Gene_name)), aes(x=Gene_name, y=Ser_frac, fill=Organism)) + ylab('Ser fraction') +
  geom_boxplot() + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

Pro_frac <- ggplot(data=subset(aa_fracs, !is.na(aa_fracs$Gene_name)), aes(x=Gene_name, y=Pro_frac, fill=Organism)) + ylab('Pro fraction') +
  geom_boxplot() + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

Cys_frac <- ggplot(data=subset(aa_fracs, !is.na(aa_fracs$Gene_name)), aes(x=Gene_name, y=Cys_frac, fill=Organism)) + ylab('Cys fraction') +
  geom_boxplot() + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

Trp_frac <- ggplot(data=subset(aa_fracs, !is.na(aa_fracs$Gene_name)), aes(x=Gene_name, y=Trp_frac, fill=Organism)) + ylab('Trp fraction') +
  geom_boxplot() + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

Tyr_frac <- ggplot(data=subset(aa_fracs, !is.na(aa_fracs$Gene_name)), aes(x=Gene_name, y=Tyr_frac, fill=Organism)) + ylab('Tyr fraction') +
  geom_boxplot() + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

Arg_frac <- ggplot(data=subset(aa_fracs, !is.na(aa_fracs$Gene_name)), aes(x=Gene_name, y=Arg_frac, fill=Organism)) + ylab('Arg fraction') +
  geom_boxplot() + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

His_frac <- ggplot(data=subset(aa_fracs, !is.na(aa_fracs$Gene_name)), aes(x=Gene_name, y=His_frac, fill=Organism)) + ylab('His fraction') +
  geom_boxplot() + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

Gln_frac <- ggplot(data=subset(aa_fracs, !is.na(aa_fracs$Gene_name)), aes(x=Gene_name, y=Gln_frac, fill=Organism)) + ylab('Gln fraction') +
  geom_boxplot() + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

Gly_frac <- ggplot(data=subset(aa_fracs, !is.na(aa_fracs$Gene_name)), aes(x=Gene_name, y=Gly_frac, fill=Organism)) + ylab('Gly fraction') +
  geom_boxplot() + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

Ser2_frac <- ggplot(data=subset(aa_fracs, !is.na(aa_fracs$Gene_name)), aes(x=Gene_name, y=Ser2_frac, fill=Organism)) + ylab('Ser2 fraction') +
  geom_boxplot() + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

Asp_frac <- ggplot(data=subset(aa_fracs, !is.na(aa_fracs$Gene_name)), aes(x=Gene_name, y=Asp_frac, fill=Organism)) + ylab('Asp fraction') +
  geom_boxplot() + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

Glu_frac <- ggplot(data=subset(aa_fracs, !is.na(aa_fracs$Gene_name)), aes(x=Gene_name, y=Glu_frac, fill=Organism)) + ylab('Glu fraction') +
  geom_boxplot() + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

Asn_frac <- ggplot(data=subset(aa_fracs, !is.na(aa_fracs$Gene_name)), aes(x=Gene_name, y=Asn_frac, fill=Organism)) + ylab('Asn fraction') +
  geom_boxplot() + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

Lys_frac <- ggplot(data=subset(aa_fracs, !is.na(aa_fracs$Gene_name)), aes(x=Gene_name, y=Lys_frac, fill=Organism)) + ylab('Lys fraction') +
  geom_boxplot() + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

Val_frac <- ggplot(data=subset(aa_fracs, !is.na(aa_fracs$Gene_name)), aes(x=Gene_name, y=Val_frac, fill=Organism)) + ylab('Val fraction') +
  geom_boxplot() + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

Ile_frac <- ggplot(data=subset(aa_fracs, !is.na(aa_fracs$Gene_name)), aes(x=Gene_name, y=Ile_frac, fill=Organism)) + ylab('Ile fraction') +
  geom_boxplot() + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

Met_frac <- ggplot(data=subset(aa_fracs, !is.na(aa_fracs$Gene_name)), aes(x=Gene_name, y=Met_frac, fill=Organism)) + ylab('Met fraction') +
  geom_boxplot() + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

Ala_frac <- ggplot(data=subset(aa_fracs, !is.na(aa_fracs$Gene_name)), aes(x=Gene_name, y=Ala_frac, fill=Organism)) + ylab('Ala fraction') +
  geom_boxplot() + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

Thr_frac <- ggplot(data=subset(aa_fracs, !is.na(aa_fracs$Gene_name)), aes(x=Gene_name, y=Thr_frac, fill=Organism)) + ylab('Thr fraction') +
  geom_boxplot() + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

ggarrange(Phe_fracs, Leu_frac, Leu2_frac, Ser_frac, Pro_frac,
          Cys_frac, Trp_frac, Tyr_frac , Arg_frac , His_frac,
          Gln_frac, Gly_frac, Ser2_frac, Asp_frac, Glu_frac,
          Asn_frac, Lys_frac, Val_frac, Ile_frac, Met_frac,
          Ala_frac, Thr_frac,
          common.legend = TRUE,
          labels = c('L', 'L', 'L', 'L','G',
                     'L', 'L', 'L', 'L', 'G',
                     'G', 'L', 'L', 'L', 'L', 
                     'G', 'G', 'L', 'L', 'L',
                     'L', 'G'),
          ncol = 5, nrow = 5)

#ylim(0, 0.2)
#ylim(0, 130)

####Testing####
Val_frac <- ggplot(data=subset(aa_fracs, !is.na(aa_fracs$Gene_name)), aes(x=Organism, y=Val_frac)) + ylab('Val fraction') +
  geom_boxplot() + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
Val_frac
