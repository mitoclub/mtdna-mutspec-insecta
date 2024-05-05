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

#set to true to use inverted skews and nucls (first calculate this all)
INVERT = FALSE
nwk_path <- "/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/MIDORI/all_insects_co1seqs_sp_names.aln.treefile"
tree <- read.tree(file=nwk_path)

#edit Species_name col to fit our data
cocks_terms <- read.csv(file="/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/DescriptiveStat/midori_Blattodea_skew_onehot_encoded.csv", header = TRUE)
cocks_terms <- cocks_terms %>% mutate(Species_name = str_replace(Species_name, "_[0-9]*$", ""))
codon_table <- read.csv(file='/mnt/data/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/DescriptiveStat/codontable_midori_Blattodea.csv')
codon_table <- codon_table %>% mutate(Species_name = str_replace(Species_name, "_[0-9]*$", ""))
cocks_terms <- cocks_terms[order(cocks_terms$Species_name), ]
codon_table <- codon_table[order(codon_table$Species_name), ]

####SIMPLE PCA MS12####
ms12_internal <- read.csv(file = "/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/NemuPipeline/AllInsects/midori_all_insects_mutspec/mutspec_as_header_12internal.csv", header = TRUE)
ms_meta <- read_tsv(file = '/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/NemuPipeline/AllInsects/midori_all_insects_mutspec/msMetaData.tsv')
blattodea_meta <- filter(ms_meta, Order == 'Blattodea_85823')
colnames(ms12_internal)[1] = "Species_name"
ms12_internal <- ms12_internal[(ms12_internal$Species_name %in% blattodea_meta$Species),]
blattodea_meta <- blattodea_meta[(blattodea_meta$Species %in% ms12_internal$Species_name),]
rownames(ms12_internal) <- ms12_internal$Species_name
ms12_internal$Species_name <- NULL

corr_matrix <- cor(ms12_internal)
ggcorrplot(corr_matrix)
data_pca <- prcomp(ms12_internal, scale. =TRUE)
summary(data_pca)
fviz_eig(data_pca, addlabels = TRUE)
#data.pca$loadings[, 1:2]
fviz_pca_var(data_pca, col.var = "black")
autoplot(data_pca, data=blattodea_meta, colour='Family')

####TREE STUFF####
#DO NOT EXECUTE SIMPLE PCA!!!
#pruning tree and db
chk <- name.check(tree, cocks_terms$Species_name, data.names = cocks_terms$Species_name)
summary(chk)

ecomorph_tree <- drop.tip(tree, chk$tree_not_data)
cocks_terms <- cocks_terms[(cocks_terms$Species_name %in% ecomorph_tree$tip.label),]
codon_table <- codon_table[(codon_table$Species_name %in% ecomorph_tree$tip.label),]
name.check(ecomorph_tree, cocks_terms$Species_name, data.names = cocks_terms$Species_name)
length(ecomorph_tree$tip.label)
ecomorph_tree



####PGLS####
#droping all genes but co1
morph_data <- cocks_terms[cocks_terms$Gene_name %in% c('CO1'), ]
codon_table <- codon_table[codon_table$Gene_name %in% c('CO1'), ]

morph_tree <- ecomorph_tree
rownames(morph_data) <- morph_data$Species_name
morph_data$Species_name <- NULL
rownames(codon_table) <- codon_table$Species_name
codon_table$Species_name <- NULL

spp<-rownames(morph_data)
corBM<-corBrownian(phy = morph_tree, form = ~spp)
if (INVERT == TRUE){
pgls <- gls(GAskew~TCskew, data = morph_data, correlation = corBM)
summary(pgls)
pgls <- gls(GAskew~Cockroaches, data = morph_data, correlation = corBM)
summary(pgls)
pgls <- gls(GAskew~Termites.w..workers + Cockroaches, data = morph_data, correlation = corBM)
summary(pgls)
}else{
  pgls <- gls(CTskew~AGskew, data = morph_data, correlation = corBM)
  summary(pgls)
  pgls <- gls(CTskew~Cockroaches, data = morph_data, correlation = corBM)
  summary(pgls)
  pgls <- gls(CTskew~Termites.w..workers + Cockroaches, data = morph_data, correlation = corBM)
  summary(pgls)
}

####TESTING####
ms12_internal <- read.csv(file = "/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/NemuPipeline/AllInsects/midori_all_insects_mutspec/mutspec_as_header_12internal.csv", header = TRUE)
colnames(ms12_internal)[1] = "Species_name"
ms12_internal <- ms12_internal[(ms12_internal$Species_name %in% morph_tree$tip.label),]
chk <- name.check(morph_tree, ms12_internal$Species_name, data.names = ms12_internal$Species_name)
morph_tree <- drop.tip(morph_tree, chk$tree_not_data)
name.check(morph_tree, ms12_internal$Species_name, data.names = ms12_internal$Species_name)
morph_data <- morph_data[(rownames(morph_data) %in% morph_tree$tip.label),]
t <- cbind(morph_data, ms12_internal[, -which(names(ms12_internal) %in% c("Species_name"))])

spp<-rownames(morph_data)
corBM<-corBrownian(phy = morph_tree, form = ~spp)
if(INVERT == TRUE){
pgls <- gls(G.A~Cockroaches, data = t, correlation = corBM)
summary(pgls)
pgls <- gls(G.A~Termites.w..workers, data = t, correlation = corBM)
summary(pgls)
}else{
  pgls <- gls(T.C~Cockroaches, data = t, correlation = corBM)
  summary(pgls)
  pgls <- gls(T.C~Termites.w..workers, data = t, correlation = corBM)
  summary(pgls)
}
####PCA####

#DO NOT EXECUTE SIMPLE PCA MS12, PGLS AND TESTING SNIPPETS!!!

#getting means for all the genes
morph_data <- ddply(cocks_terms,"Species_name",numcolwise(mean))
codon_table <- ddply(codon_table,"Species_name",numcolwise(mean))

#merging nucleotide freqs and skews
codon_table <- subset(codon_table, select = c('Species_name','X.A', 'X.T', 'X.G', 'X.C'))
#swapping codon table's nucleotides, required for Blattodea
if (INVERT == TRUE){
colnames(codon_table) <- c ('Species_name','T', 'A', 'C', 'G')
}
morph_tree <- ecomorph_tree
rownames(morph_data) <- morph_data$Species_name
morph_data$Species_name <- NULL
rownames(codon_table) <- codon_table$Species_name
codon_table$Species_name <- NULL
if (INVERT == TRUE){
pca_data <- subset(morph_data, select = c('TCskew', 'GAskew'))
}else{
  pca_data <- subset(morph_data, select = c('CTskew', 'AGskew'))
}
pca_data <- pca_data[order(row.names(pca_data)),]
codon_table <- codon_table[order(row.names(codon_table)),]
pca_data <- merge(pca_data, codon_table, by = "row.names")
rownames(pca_data) <- pca_data$Row.names
pca_data$Row.names <- NULL

morph_pca <- phyl.pca(morph_tree, pca_data)
morph_pca

#Separating organisms for visualization
organism_types <- morph_data
organism_types['Sub.social.Cryptocercus'][organism_types['Sub.social.Cryptocercus'] == 1] <- 2
organism_types['Termites.w..workers'][organism_types['Termites.w..workers'] == 1] <- 3
organism_types['Termites.w.o.workers'][organism_types['Termites.w.o.workers'] == 1] <- 4
organism_types <- unite(organism_types, Organism, c("Cockroaches", "Sub.social.Cryptocercus", "Termites.w..workers", "Termites.w.o.workers"), remove = TRUE)

par(mar=c(4.1,4.1,2.1,1.1),las=1) ## set margins
plot(morph_pca,main="")

par(cex.axis=0.8,mar=c(5.1,5.1,1.1,1.1))
phylomorphospace(morph_tree,
                 scores(morph_pca)[,1:2],
                 ftype="off",node.size=c(0,1),bty="n",las=1,
                 xlab="PC1",
                 ylab=expression(paste("PC2")))
#get tip labels (species, basically)
#text(scores(morph_pca)[,1], scores(morph_pca)[,2], rownames(morph_data), cex=1, adj=c(NA,1))

eco<-setNames(organism_types$Organism,rownames(organism_types))
ECO<-to.matrix(eco,unique(organism_types$Organism))
tiplabels(pie=ECO[morph_tree$tip.label,],cex=0.5)
#ALWAYS CHECK IF LEGEN CORRESPONDS WITH POINTS, DAMNIT!
legend(x="topright",legend=c('Cockroaches', 'Termites w/ workers', 'Termites w/o workers', 'Sub social Cryptocercus'),cex=0.8,pch=21,
       pt.bg=rainbow(n=length(unique(organism_types$Organism))),pt.cex=1.5, text.width = 0.09)

#simplePCA (non phylo)
corr_matrix <- cor(pca_data)
ggcorrplot(corr_matrix)
data.pca <- princomp(corr_matrix)
summary(data.pca)
fviz_eig(data.pca, addlabels = TRUE)
data.pca$loadings[, 1:2]
fviz_pca_var(data.pca, col.var = "black")


# mutspecPCA, sucks butt. MAYBE NOT ANYMORE :)
pca_mutspec <- t[,-c(2,4:14,16:18,20)]
morph_pca <- phyl.pca(morph_tree, pca_mutspec)
morph_pca

par(mar=c(4.1,4.1,2.1,1.1),las=1) ## set margins
plot(morph_pca,main="")

par(cex.axis=0.8,mar=c(5.1,5.1,1.1,1.1))
phylomorphospace(morph_tree,
                 scores(morph_pca)[,1:2],
                 ftype="off",node.size=c(0,1),bty="n",las=1,
                 xlab="PC1",
                 ylab=expression(paste("PC2")))
#get tip labels (species, basically)
text(scores(morph_pca)[,1], scores(morph_pca)[,2], rownames(morph_data), cex=1, adj=c(NA,1))
eco<-setNames(t$Cockroaches,rownames(t))
ECO<-to.matrix(eco,unique(t$Cockroaches))
tiplabels(pie=ECO[morph_tree$tip.label,],cex=0.5)
legend(x="topright",legend=c('Cockroaches', 'Termites'),cex=0.8,pch=21,
       pt.bg=rainbow(n=length(unique(t$Cockroaches))),pt.cex=1.5, text.width = 0.1)
