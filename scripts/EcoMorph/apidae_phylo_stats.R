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


nwk_path <- "/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/MIDORI/all_insects_co1seqs_sp_names.aln.treefile"
tree <- read.tree(file=nwk_path)

skews <- read.csv('/mnt/data/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/DescriptiveStat/midori_Apidae_CO1_skew.csv')
skews$X <- NULL
skews$Gene_name <- NULL
skews <- skews %>% mutate(Species_name = str_replace(Species_name, "_[0-9]*$", ""))

ms12_internal <- read.csv(file = "/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/NemuPipeline/AllInsects/midori_all_insects_mutspec/mutspec_as_header_12internal.csv", header = TRUE)
colnames(ms12_internal)[1] = "Species_name"
ms12_internal <- ms12_internal[(ms12_internal$Species_name %in% skews$Species_name),]
skews <- skews[(skews$Species_name %in% ms12_internal$Species_name),]
skews <- cbind(skews, ms12_internal[, -which(names(ms12_internal) %in% c("Species_name"))])

chk <- name.check(tree, skews$Species_name, data.names = skews$Species_name)
summary(chk)
ecomorph_tree <- drop.tip(tree, chk$tree_not_data)
skews <- skews[(skews$Species_name %in% ecomorph_tree$tip.label),]
name.check(ecomorph_tree, skews$Species_name, data.names = skews$Species_name)

rownames(skews) <- skews$Species_name
skews[1] <- NULL

organism_types <- skews
skews$Organism <- NULL

morph_pca <- phyl.pca(ecomorph_tree, skews[c(1,2,4,8)])
morph_pca
par(mar=c(4.1,4.1,2.1,1.1),las=1) ## set margins
plot(morph_pca,main="")

par(cex.axis=0.8,mar=c(5.1,5.1,1.1,1.1))
phylomorphospace(ecomorph_tree,
                 scores(morph_pca)[,1:2],
                 ftype="off",node.size=c(0,1),bty="n",las=1,
                 xlab="PC1",
                 ylab=expression(paste("PC2")))
eco<-setNames(organism_types$Organism,rownames(organism_types))
ECO<-to.matrix(eco,unique(organism_types$Organism))
tiplabels(pie=ECO[ecomorph_tree$tip.label,],cex=0.5)
legend(x="topleft",legend=unique(organism_types$Organism),cex=0.8,pch=21,
       pt.bg=rainbow(n=length(unique(organism_types$Organism))),pt.cex=1.5, text.width = 0.2)
