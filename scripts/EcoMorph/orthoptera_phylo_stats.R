#courtesy of "Phylogenetic comparative methods in R"
library(ape)
library(geiger)
library(nlme)
library(plyr)
library(dplyr)
library(ggplot2)
library(ggfortify)
library(stringr)
library('corrr')
library(ggcorrplot)
library("FactoMineR")
library(factoextra)
library(tidyr)
library(readr)

#Set to true to use mutspecs, else use skews (they are more plentiful)
MUTSPEC = TRUE

nwk_path <- "/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/MIDORI/all_insects_co1seqs_sp_names.aln.treefile"
tree <- read.tree(file=nwk_path)

####NUCLS####
ortho_traits <- read.csv('/mnt/data/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/EcoMorphDBs/orthoptera_lifetraits.tsv', header=TRUE, sep = '\t')
ortho_traits <- ortho_traits %>% drop_na(c('Length', 'Clutch_S'))
ortho_traits <- ortho_traits[c(1, 2, 11)]

if (MUTSPEC == TRUE){
 ms12_internal <- read.csv(file = "/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/NemuPipeline/AllInsects/midori_all_insects_mutspec/mutspec_as_header_12internal.csv", header = TRUE)
 colnames(ms12_internal)[1] = "Species"
 ms12_internal <- ms12_internal[(ms12_internal$Species %in% ortho_traits$Species),]
 ortho_traits <- ortho_traits[(ortho_traits$Species %in% ms12_internal$Species),]
 ortho_data <- cbind(ortho_traits, ms12_internal[, -which(names(ms12_internal) %in% c("Species"))])
 
 #meta_data <- read.csv(file='/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/NemuPipeline/AllInsects/midori_all_insects_mutspec/msMetaData.tsv', sep='\t')
 #meta_data <- meta_data[(meta_data$Species %in% ortho_data$Species),]
 #ortho_data <- cbind(ortho_data, meta_data[, -which(names(meta_data) %in% c("Species"))])
 
 chk <- name.check(tree, ortho_data$Species, data.names = ortho_data$Species)
 summary(chk)
 ecomorph_tree <- drop.tip(tree, chk$tree_not_data)
 ortho_data <- ortho_data[(ortho_data$Species %in% ecomorph_tree$tip.label),]
 name.check(ecomorph_tree, ortho_data$Species, data.names = ortho_data$Species)
 
 #extracting genuses for PCA coloring
 taxa_coloring <- sub('_.*','',ortho_data$Species)
 
 rownames(ortho_data) <- ortho_data$Species
 ortho_data[1] <- NULL
 
}else{
 skews <- read.csv('/mnt/data/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/DescriptiveStat/midori_Insecta_CO1_skew.csv')
 skews <- skews %>% mutate(Species_name = str_replace(Species_name, "_[0-9]*$", ""))
 skews <- skews[(skews$Species %in% ortho_traits$Species),]
 #Drop gene name, X column (extra index bug) and Organism type (all NA)
 skews[c(1, 4, 6)] <- list(NULL)
 ortho_traits <- ortho_traits[(ortho_traits$Species %in% skews$Species_name),]
 ortho_data <- cbind(ortho_traits, skews[, -which(names(skews) %in% c("Species_name"))])
 chk <- name.check(tree, ortho_data$Species, data.names = ortho_data$Species)
 summary(chk)
 
 ecomorph_tree <- drop.tip(tree, chk$tree_not_data)
 ortho_data <- ortho_data[(ortho_data$Species %in% ecomorph_tree$tip.label),]
 name.check(ecomorph_tree, ortho_data$Species, data.names = ortho_data$Species)
 
 #extracting genuses for PCA coloring
 taxa_coloring <- sub('_.*','',ortho_data$Species)
 
 rownames(ortho_data) <- ortho_data$Species
 ortho_data[1] <- NULL
}

morph_pca <- phyl.pca(ecomorph_tree, ortho_data)
morph_pca
par(mar=c(4.1,4.1,2.1,1.1),las=1) ## set margins
plot(morph_pca,main="")

par(cex.axis=0.8,mar=c(5.1,5.1,1.1,1.1))
phylomorphospace(ecomorph_tree,
                 scores(morph_pca)[,1:2],
                 ftype="off",node.size=c(0,1),bty="n",las=1,
                 xlab="PC1",
                 ylab=expression(paste("PC2")))
eco<-setNames(taxa_coloring,rownames(ortho_data))
ECO<-to.matrix(eco,unique(taxa_coloring))
tiplabels(pie=ECO[ecomorph_tree$tip.label,],cex=0.5)
legend(x="topleft",legend=unique(taxa_coloring),cex=0.8,pch=21,
       pt.bg=rainbow(n=length(unique(taxa_coloring))),pt.cex=1.5, text.width = 12)

####AA####

nwk_path <- "/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/MIDORI/all_insects_co1seqs_sp_names.aln.treefile"
tree <- read.tree(file=nwk_path)

ortho_traits <- read.csv('/mnt/data/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/EcoMorphDBs/orthoptera_lifetraits.tsv', header=TRUE, sep = '\t')
ortho_traits <- ortho_traits %>% drop_na(c('Length', 'Clutch_S'))
ortho_traits <- ortho_traits[c(1, 2, 11)]

fracs <- read.csv('/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/DescriptiveStat/midori_Orthoptera_CO1_aa_fracs.tsv', sep='\t')
fracs <- fracs %>% mutate(Species_name = str_replace(Species_name, "_[0-9]*$", ""))
fracs['Gainer/Loser'] = (fracs['Phe_frac'] + fracs['Leu_frac'] + fracs['Cys_frac'] + fracs['Trp_frac'] + fracs['Ser2_frac'] + fracs['Gly_frac'] + fracs['Val_frac'])/(fracs['Pro_frac'] + fracs['Thr_frac'] + fracs['Asn_frac'] + fracs['Lys_frac'] + fracs['His_frac'] + fracs['Gln_frac'])

fracs <- fracs[(fracs$Species_name %in% ortho_traits$Species),]
ortho_traits <- ortho_traits[(ortho_traits$Species %in% fracs$Species_name),]
ortho_data <- cbind(ortho_traits, fracs[, -which(names(fracs) %in% c("Species_name"))])

ortho_data$Organism <- NULL
ortho_data$Gene_name <- NULL

chk <- name.check(tree, ortho_data$Species, data.names = ortho_data$Species)
summary(chk)
ecomorph_tree <- drop.tip(tree, chk$tree_not_data)
ortho_data <- ortho_data[(ortho_data$Species %in% ecomorph_tree$tip.label),]
name.check(ecomorph_tree, ortho_data$Species, data.names = ortho_data$Species)

#extracting genuses for PCA coloring
taxa_coloring <- sub('_.*','',ortho_data$Species)

rownames(ortho_data) <- ortho_data$Species
ortho_data[1] <- NULL

morph_pca <- phyl.pca(ecomorph_tree, ortho_data[c(1,2,3,4,8,9,14,20,7,12,13,18,19,24)])
morph_pca
par(mar=c(4.1,4.1,2.1,1.1),las=1) ## set margins
plot(morph_pca,main="")

par(cex.axis=0.8,mar=c(5.1,5.1,1.1,1.1))
phylomorphospace(ecomorph_tree,
                 scores(morph_pca)[,1:2],
                 ftype="off",node.size=c(0,1),bty="n",las=1,
                 xlab="PC1",
                 ylab=expression(paste("PC2")))
eco<-setNames(taxa_coloring,rownames(ortho_data))
ECO<-to.matrix(eco,unique(taxa_coloring))
tiplabels(pie=ECO[ecomorph_tree$tip.label,],cex=0.5, )
legend(x="topleft",legend=unique(taxa_coloring),cex=0.8,pch=21,
      pt.bg=rainbow(n=length(unique(taxa_coloring))),pt.cex=1.5, text.width = 9)

