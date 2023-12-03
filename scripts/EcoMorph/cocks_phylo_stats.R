#courtesy of "Phylogenetic comparative methods in R"
library(ape)
library(geiger)
library(nlme)
library(plyr)
library(dplyr)
library(stringr)

nwk_path <- "/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/MIDORI/all_insects_co1seqs_sp_names.aln.treefile"
tree <- read.tree(file=nwk_path)

#edit Species_name col to fit our data
cocks_terms <- read.csv(file="/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/DescriptiveStat/midori_blattodea_skew_onehot_encoded.csv", header = TRUE)
cocks_terms <- cocks_terms %>% mutate(Species_name = str_replace(Species_name, "_[0-9]*$", ""))
codon_table <- read.csv(file='/mnt/data/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/DescriptiveStat/codontable_midori_blattodea.csv')
codon_table <- codon_table %>% mutate(Species_name = str_replace(Species_name, "_[0-9]*$", ""))
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
pgls <- gls(GAskew~TCskew, data = morph_data, correlation = corBM)
summary(pgls)
pgls <- gls(GAskew~Cockroaches, data = morph_data, correlation = corBM)
summary(pgls)
pgls <- gls(GAskew~Termites.w..workers, data = morph_data, correlation = corBM)
summary(pgls)

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
pgls <- gls(GAskew~G.A + Cockroaches, data = t, correlation = corBM)
summary(pgls)
pgls <- gls(GAskew~G.A + Termites.w..workers, data = t, correlation = corBM)
summary(pgls)

####PCA####

#DO NOT EXECUTE PGLS AND TESTING SNIPPETS!!!

#getting means for all the genes
morph_data <- ddply(cocks_terms,"Species_name",numcolwise(mean))
codon_table <- ddply(codon_table,"Species_name",numcolwise(mean))

#merging nucleotide freqs and skews
codon_table <- subset(codon_table, select = c('X.A', 'X.T', 'X.G', 'X.C'))
morph_tree <- ecomorph_tree
rownames(morph_data) <- morph_data$Species_name
morph_data$Species_name <- NULL
rownames(codon_table) <- codon_table$Species_name
codon_table$Species_name <- NULL
pca_data <- subset(morph_data, select = c('TCskew', 'GAskew'))
pca_data <- pca_data[order(row.names(pca_data)),]
codon_table <- codon_table[order(row.names(codon_table)),]
pca_data <- merge(pca_data, codon_table, by = "row.names")
rownames(pca_data) <- pca_data$Row.names
pca_data$Row.names <- NULL

morph_pca <- phyl.pca(morph_tree, pca_data)
morph_pca

par(mar=c(4.1,4.1,2.1,1.1),las=1) ## set margins
plot(morph_pca,main="")

par(cex.axis=0.8,mar=c(5.1,5.1,1.1,1.1))
phylomorphospace(morph_tree,
                 scores(morph_pca)[,1:2],
                 ftype="off",node.size=c(0,1),bty="n",las=1,
                 xlab="PC1",
                 ylab=expression(paste("PC2")))
eco<-setNames(morph_data$Cockroaches,rownames(morph_data))
ECO<-to.matrix(eco,unique(morph_data$Cockroaches))
tiplabels(pie=ECO[morph_tree$tip.label,],cex=0.5)
legend(x="topright",legend=c('Termites', 'Cockroaches'),cex=0.8,pch=21,
       pt.bg=rainbow(n=length(unique(morph_data$Cockroaches))),pt.cex=1.5, text.width = 0.09)

# mutspecPCA, sucks butt
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
eco<-setNames(t$Cockroaches,rownames(t))
ECO<-to.matrix(eco,unique(t$Cockroaches))
tiplabels(pie=ECO[morph_tree$tip.label,],cex=0.5)
legend(x="topright",legend=c('Termites', 'Cockroaches'),cex=0.8,pch=21,
       pt.bg=rainbow(n=length(unique(t$Cockroaches))),pt.cex=1.5, text.width = 0.1)
