#courtesy of "Phylogenetic comparative methods in R"
library(ape)
library(geiger)
library(nlme)
library(plyr)
library(dplyr)
library(stringr)

nwk_path <- "/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/MIDORI/all_insects_co1seqs_sp_names.aln.treefile"
tree <- read.tree(file=nwk_path)

#edit species col to fit our data
diptera_traits <- read.csv(file="/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/EcoMorphDBs/DipteraTrait_label_encoded.csv", header = TRUE)
diptera_traits <- diptera_traits %>% mutate(species = str_replace(species, " ", "_"))
diptera_skew <- read.csv(file="/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/DescriptiveStat/midori_Diptera_skew_onehot_encoded.csv", header = TRUE)
diptera_skew <- diptera_skew %>% mutate(Species_name = str_replace(Species_name, "_[0-9]*$", ""))
codon_table <- read.csv(file='/mnt/data/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/DescriptiveStat/codontable_midori_Diptera.csv')
codon_table <- codon_table %>% mutate(Species_name = str_replace(Species_name, "_[0-9]*$", ""))
diptera_traits <- diptera_traits[order(diptera_traits$species), ]
diptera_skew <- diptera_skew[order(diptera_skew$Species_name), ]
codon_table <- codon_table[order(codon_table$Species_name), ]
#pruning tree and db
chk <- name.check(tree, diptera_traits$species, data.names = diptera_traits$species)
summary(chk)

ecomorph_tree <- drop.tip(tree, chk$tree_not_data)
diptera_traits <- diptera_traits[(diptera_traits$species %in% ecomorph_tree$tip.label),]
diptera_skew <- diptera_skew[(diptera_skew$Species_name %in% diptera_traits$species),]
codon_table <- codon_table[(codon_table$Species_name %in% diptera_traits$species),]
diptera_traits <- diptera_traits[(diptera_traits$species %in% diptera_skew$Species_name),]
chk <- name.check(ecomorph_tree, diptera_traits$species, data.names = diptera_traits$species)
ecomorph_tree <- drop.tip(ecomorph_tree, chk$tree_not_data)
length(ecomorph_tree$tip.label)
ecomorph_tree

#choosing between larvae and adult diets and removing non-morphological data
morph_data <- diptera_traits[diptera_traits$traitName == 'adultDiet',]
morph_data$traitName <- NULL
diptera_skew <- diptera_skew[diptera_skew$Gene_name %in% c('CO1'), ]
codon_table <- codon_table[codon_table$Gene_name %in% c('CO1'), ]
#merging diptera traits with skews
colnames(codon_table)[1] <- "species"
colnames(diptera_skew)[1] <- "species"
df_list <- list(diptera_skew, morph_data)
morph_data <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)
morph_data <- subset(morph_data, select = -c(3, 5, 6, 7))
#plus extracting genuses for PCA coloring
morph_genuses <- sub('_.*','',morph_data$species)


####PCA####
chk <- name.check(ecomorph_tree, morph_data$species, data.names = morph_data$species)

#removing duplicate species from tree, how's it even possible? Might be usefull, might not be
morph_tree <- drop.tip(ecomorph_tree, 'Culicoides_obsoletus')
morph_data <- morph_data[(morph_data$species %in% morph_tree$tip.label),]

name.check(morph_tree, morph_data$species, data.names = morph_data$species)
rownames(morph_data) <- morph_data$species
morph_data$species <- NULL
morph_pca <- phyl.pca(morph_tree, morph_data)
morph_pca

par(mar=c(4.1,4.1,2.1,1.1),las=1) ## set margins
plot(morph_pca,main="")

par(cex.axis=0.8,mar=c(5.1,5.1,1.1,1.1))
phylomorphospace(morph_tree,
                 scores(morph_pca)[,1:2],
                 ftype="off",node.size=c(0,1),bty="n",las=1,
                 xlab="PC1",
                 ylab=expression(paste("PC2")))
#text(scores(morph_pca)[,1], scores(morph_pca)[,2], rownames(morph_data), cex=1, adj=c(NA,1))
eco<-setNames(diptera_skew$Brachycera,rownames(morph_data))
ECO<-to.matrix(eco,unique(diptera_skew$Brachycera))
tiplabels(pie=ECO[morph_tree$tip.label,],cex=0.5)
legend(x="topleft",legend=c('Nematocera','Brachycera'),cex=0.8,pch=21,
       pt.bg=rainbow(n=length(unique(diptera_skew$Brachycera))),pt.cex=1.5, text.width = 2)


####PIC####
habitat <- setNames(morph_data[,"habitat"], rownames(morph_data))
traitValue <- setNames(morph_data[,"traitValue"], rownames(morph_data))

pic_habitat <- pic(habitat, morph_tree)
pic_traitValue <- pic(traitValue, morph_tree)

fit_pic <- lm(pic_habitat~pic_traitValue+0)

## set margins
par(mar=c(5.1,5.1,1.1,1.1))
## graph scatterplot of contrasts
plot(pic_habitat~pic_traitValue,
     xlab="PICs for log(traitValue)",
     ylab="PICs for log(habitat)",
     pch=21,bg="gray",cex=1.2,las=1,
     cex.axis=0.7,cex.lab=0.9,bty="n")
## add gridlines to the plot
abline(h=0,lty="dotted")
abline(v=0,lty="dotted")
## reset graphing limits of the plot to the
## x/y range of our PICs
clip(min(pic_traitValue),max(pic_traitValue),
     min(pic_habitat),max(pic_habitat))
## graph our fitted line
abline(fit_pic,lwd=2,col="darkgray")

####PGLS####

spp<-rownames(morph_data)
corBM<-corBrownian(phy = morph_tree, form = ~spp)
pgls <- gls(habitat~traitValue, data = morph_data, correlation = corBM)
anova(pgls)

####TESTING####
ms12_internal <- read.csv(file = "/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/NemuPipeline/AllInsects/midori_all_insects_mutspec/mutspec_as_header_12internal.csv", header = TRUE)
colnames(ms12_internal)[1] = "species"
ms12_internal <- ms12_internal[(ms12_internal$species %in% morph_tree$tip.label),]
chk <- name.check(morph_tree, ms12_internal$species, data.names = ms12_internal$species)
morph_tree <- drop.tip(morph_tree, chk$tree_not_data)
name.check(morph_tree, ms12_internal$species, data.names = ms12_internal$species)
morph_data <- morph_data[(rownames(morph_data) %in% morph_tree$tip.label),]
t <- cbind(morph_data, ms12_internal[, -which(names(ms12_internal) %in% c("species"))])

spp<-rownames(morph_data)
corBM<-corBrownian(phy = morph_tree, form = ~spp)
pgls <- gls(T.C~G.A + habitat + traitValue, data = t, correlation = corBM)
anova(pgls)




####ONEHOT_ENCODED####
diptera_traits <- read.csv(file="/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/EcoMorphDBs/DipteraTrait_onehot_encoded.csv", header = TRUE)
diptera_traits <- diptera_traits %>% mutate(species = str_replace(species, " ", "_"))
diptera_skew <- read.csv(file="/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/DescriptiveStat/midori_Diptera_skew_onehot_encoded.csv", header = TRUE)
diptera_skew <- diptera_skew %>% mutate(Species_name = str_replace(Species_name, "_[0-9]*$", ""))
codon_table <- read.csv(file='/mnt/data/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/DescriptiveStat/codontable_midori_Diptera.csv')
codon_table <- codon_table %>% mutate(Species_name = str_replace(Species_name, "_[0-9]*$", ""))
diptera_traits <- diptera_traits[order(diptera_traits$species), ]
diptera_skew <- diptera_skew[order(diptera_skew$Species_name), ]
codon_table <- codon_table[order(codon_table$Species_name), ]

#pruning tree and db
chk <- name.check(tree, diptera_traits$species, data.names = diptera_traits$species)
summary(chk)

ecomorph_tree <- drop.tip(tree, chk$tree_not_data)
diptera_traits <- diptera_traits[(diptera_traits$species %in% ecomorph_tree$tip.label),]
diptera_skew <- diptera_skew[(diptera_skew$Species_name %in% diptera_traits$species),]
codon_table <- codon_table[(codon_table$Species_name %in% diptera_traits$species),]
diptera_traits <- diptera_traits[(diptera_traits$species %in% diptera_skew$Species_name),]
chk <- name.check(ecomorph_tree, diptera_traits$species, data.names = diptera_traits$species)
ecomorph_tree <- drop.tip(ecomorph_tree, chk$tree_not_data)
length(ecomorph_tree$tip.label)
ecomorph_tree

#chosing between larvae and adult diets and removing non-morphological data
morph_data <- diptera_traits[diptera_traits$traitName == 'adultDiet',]
morph_data$traitName <- NULL
diptera_skew <- diptera_skew[diptera_skew$Gene_name %in% c('CO1'), ]
codon_table <- codon_table[codon_table$Gene_name %in% c('CO1'), ]
#merging diptera traits with skews
colnames(codon_table)[1] <- "species"
colnames(diptera_skew)[1] <- "species"
df_list <- list(diptera_skew, morph_data)
morph_data <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)
#removing excessive cols and cols with only zero values
morph_data <- subset(morph_data, select = -c(3, 5, 6, 7))
morph_data <- morph_data[, colSums(morph_data != 0) > 0]
#PCA
chk <- name.check(ecomorph_tree, morph_data$species, data.names = morph_data$species)

#removing duplicate species from tree, how's it even possible? Might be usefull, might not be
morph_tree <- drop.tip(ecomorph_tree, 'Culicoides_obsoletus')
morph_data <- morph_data[(morph_data$species %in% morph_tree$tip.label),]

name.check(morph_tree, morph_data$species, data.names = morph_data$species)
rownames(morph_data) <- morph_data$species
morph_data$species <- NULL
morph_pca <- phyl.pca(morph_tree, morph_data)
morph_pca

par(mar=c(4.1,4.1,2.1,1.1),las=1) ## set margins
plot(morph_pca,main="")

par(cex.axis=0.8,mar=c(5.1,5.1,1.1,1.1))
phylomorphospace(morph_tree,
                 scores(morph_pca)[,1:2],
                 ftype="off",node.size=c(0,1),bty="n",las=1,
                 xlab="PC1",
                 ylab=expression(paste("PC2")))
#text(scores(morph_pca)[,1], scores(morph_pca)[,2], rownames(morph_data), cex=1, adj=c(NA,1))
eco<-setNames(diptera_skew$Brachycera,rownames(morph_data))
ECO<-to.matrix(eco,unique(diptera_skew$Brachycera))
tiplabels(pie=ECO[morph_tree$tip.label,],cex=0.5)
legend(x="topleft",legend=c('Nematocera','Brachycera'),cex=0.8,pch=21,
       pt.bg=rainbow(n=length(unique(diptera_skew$Brachycera))),pt.cex=1.5, text.width = 0.5)

#PGLS
ms12_internal <- ms12_internal[(ms12_internal$species %in% morph_tree$tip.label),]
name.check(morph_tree, ms12_internal$species, data.names = ms12_internal$species)
morph_data <- morph_data[(morph_data$species %in% morph_tree$tip.label),]
t <- cbind(morph_data, ms12_internal[, -which(names(ms12_internal) %in% c("species"))])

spp<-morph_data$species
corBM<-corBrownian(phy = morph_tree, form = ~spp)
pgls <- gls(T.C~G.A + NonFeeding + Nectivore + Phytophagous + Predaceous + Saprophagous + Aquatic + Terrestrial, data = t, correlation = corBM)
summary(pgls)


####GENUS_COLAPSE####
nwk_path <- "/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/MIDORI/genus_collapsed_co1seqs_sp_names.aln.treefile"
tree <- read.tree(file=nwk_path)

#edit cols to fit our data
diptera_traits <- read.csv(file="/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/EcoMorphDBs/DipteraTrait_label_encoded.csv", header = TRUE)
morph_data <- diptera_traits[diptera_traits$traitName == 'adultDiet',]
morph_data$species <- sub(' .*','',morph_data$species)
colnames(morph_data)[1] = "genus"
morph_data <- morph_data %>% select(c('genus', 'habitat', 'traitValue'))
morph_data <- ddply(morph_data,"genus",numcolwise(mean))
morph_data <- morph_data[complete.cases(morph_data),]

chk <- name.check(tree, morph_data$genus, data.names = morph_data$genus)
summary(chk)

morph_tree <- drop.tip(tree, chk$tree_not_data)
morph_data <- morph_data[(morph_data$genus %in% morph_tree$tip.label),]
name.check(morph_tree, morph_data$genus, data.names = morph_data$genus)

#remove dublicates
n_occur <- data.frame(table(morph_tree$tip.label))
dubs <- unique(morph_tree$tip.label[morph_tree$tip.label %in% n_occur$Var1[n_occur$Freq > 1]])
morph_tree <- drop.tip(morph_tree, dubs)
morph_data <- morph_data[(morph_data$genus %in% morph_tree$tip.label),]

spp<-morph_data$genus
corBM<-corBrownian(phy = morph_tree, form = ~spp)
pgls <- gls(habitat~traitValue, data = morph_data, correlation = corBM)
anova(pgls)

#mutspec
ms12_internal <- read.csv(file = "/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/NemuPipeline/AllInsects/midori_all_insects_mutspec/mutspec_as_header_12internal.csv", header = TRUE)
colnames(ms12_internal)[1] = "species"
ms12_internal$species <- sub('_.*','',ms12_internal$species)
ms12_internal <- ddply(ms12_internal,"species",numcolwise(mean))
chk <- name.check(morph_tree, ms12_internal$species, data.names = ms12_internal$species)
morph_tree <- drop.tip(morph_tree, chk$tree_not_data)
ms12_internal <- ms12_internal[(ms12_internal$species %in% morph_tree$tip.label),]
name.check(morph_tree, ms12_internal$species, data.names = ms12_internal$species)
morph_data <- morph_data[(morph_data$genus %in% morph_tree$tip.label),]
t <- cbind(morph_data, ms12_internal[, -which(names(ms12_internal) %in% c("species"))])

spp<-morph_data$genus
corBM<-corBrownian(phy = morph_tree, form = ~spp)
pgls <- gls(T.C~G.A + habitat + traitValue, data = t, correlation = corBM)
summary(pgls)


