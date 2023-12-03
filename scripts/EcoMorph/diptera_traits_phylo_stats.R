#courtesy of "Phylogenetic comparative methods in R"
library(ape)
library(geiger)
library(nlme)
library(plyr)
library(dplyr)
library(stringr)

nwk_path <- "/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/MIDORI/all_insects_co1seqs(sp_names).aln.treefile"
tree <- read.tree(file=nwk_path)

#edit species col to fit our data
diptera_traits <- read.csv(file="/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/EcoMorphDBs/DipteraTrait_label_encoded.csv", header = TRUE)
diptera_traits <- diptera_traits %>% mutate(species = str_replace(species, " ", "_"))

#pruning tree and db
chk <- name.check(tree, diptera_traits$species, data.names = diptera_traits$species)
summary(chk)

ecomorph_tree <- drop.tip(tree, chk$tree_not_data)
diptera_traits <- diptera_traits[(diptera_traits$species %in% ecomorph_tree$tip.label),]
name.check(ecomorph_tree, diptera_traits$species, data.names = diptera_traits$species)
length(ecomorph_tree$tip.label)
ecomorph_tree

#chosing between larvae and adult diets and removing non-morphological data
morph_data <- diptera_traits[diptera_traits$traitName == 'adultDiet',]
morph_data$traitName <- NULL

#plus extracting genuses for PCA coloring
morph_genuses <- sub('_.*','',morph_data$species)


####PCA####
chk <- name.check(ecomorph_tree, morph_data$species, data.names = morph_data$species)

#removing duplicate species from tree, how's it even possible?
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
eco<-setNames(morph_genuses,rownames(morph_data))
ECO<-to.matrix(eco,unique(morph_genuses))
tiplabels(pie=ECO[morph_tree$tip.label,],cex=0.5)
legend(x="topright",legend=unique(morph_genuses),cex=0.8,pch=21,
       pt.bg=rainbow(n=length(unique(morph_genuses))),pt.cex=1.5, text.width = 2)


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


#pruning tree and db
chk <- name.check(tree, diptera_traits$species, data.names = diptera_traits$species)
summary(chk)

ecomorph_tree <- drop.tip(tree, chk$tree_not_data)
diptera_traits <- diptera_traits[(diptera_traits$species %in% ecomorph_tree$tip.label),]
name.check(ecomorph_tree, diptera_traits$species, data.names = diptera_traits$species)
length(ecomorph_tree$tip.label)
ecomorph_tree

#chosing between larvae and adult diets and removing non-morphological data
morph_data <- diptera_traits[diptera_traits$traitName == 'adultDiet',]
morph_data$traitName <- NULL

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


