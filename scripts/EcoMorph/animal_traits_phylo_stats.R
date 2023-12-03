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
animal_traits <- read.csv(file="/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/EcoMorphDBs/animal_traits.csv", header = TRUE)
animal_traits <- animal_traits %>% mutate(species = str_replace(species, " ", "_"))

#pruning tree and db
chk <- name.check(tree, animal_traits$species, data.names = animal_traits$species)
summary(chk)

ecomorph_tree <- drop.tip(tree, chk$tree_not_data)
animal_traits <- animal_traits[(animal_traits$species %in% ecomorph_tree$tip.label),]
name.check(ecomorph_tree, animal_traits$species, data.names = animal_traits$species)
length(ecomorph_tree$tip.label)
ecomorph_tree

#merging duplicates in db and removing non-morphological data
morph_data <- animal_traits %>% select(c('species', 'body.mass', 'metabolic.rate', 'mass.specific.metabolic.rate'))
morph_data <- ddply(morph_data,"species",numcolwise(mean))




morph_data <- morph_data[complete.cases(morph_data),]


#plus extracting genuses for PCA coloring
morph_genuses <- sub('_.*','',morph_data$species)



####PCA####
chk <- name.check(ecomorph_tree, morph_data$species, data.names = morph_data$species)
morph_tree <- drop.tip(ecomorph_tree, chk$tree_not_data)

#removing duplicate species from tree, how's it even possible?
morph_tree <- drop.tip(morph_tree, "Tetramorium_caespitum")
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
body_mass <- setNames(morph_data[,"body.mass"], rownames(morph_data))
metabolic_rate <- setNames(morph_data[,"metabolic.rate"], rownames(morph_data))

pic_body_mass <- pic(log(body_mass), morph_tree)
pic_metabolic_rate <- pic(log(metabolic_rate), morph_tree)

fit_pic <- lm(pic_body_mass~pic_metabolic_rate+0)

## set margins
par(mar=c(5.1,5.1,1.1,1.1))
## graph scatterplot of contrasts
plot(pic_body_mass~pic_metabolic_rate,
     xlab="PICs for log(metabolic_rate)",
     ylab="PICs for log(body_mass)",
     pch=21,bg="gray",cex=1.2,las=1,
     cex.axis=0.7,cex.lab=0.9,bty="n")
## add gridlines to the plot
abline(h=0,lty="dotted")
abline(v=0,lty="dotted")
## reset graphing limits of the plot to the
## x/y range of our PICs
clip(min(pic_metabolic_rate),max(pic_metabolic_rate),
     min(pic_body_mass),max(pic_body_mass))
## graph our fitted line
abline(fit_pic,lwd=2,col="darkgray")

####PGLS####

spp<-rownames(morph_data)
corBM<-corBrownian(phy = morph_tree, form = ~spp)
pgls <- gls(log(body.mass)~log(metabolic.rate), data = morph_data, correlation = corBM)
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
pgls <- gls(T.C~G.A + log(mass.specific.metabolic.rate) + log(body.mass) + log (metabolic.rate), data = t, correlation = corBM)
summary(pgls)


####GENUS_COLAPSE####
#body_mass, metabolic_rate, mass_specific_metabolic_rate
nwk_path <- "/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/MIDORI/genus_collapsed_co1seqs_sp_names.aln.treefile"
tree <- read.tree(file=nwk_path)

#edit cols to fit our data
animal_traits <- read.csv(file="/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/EcoMorphDBs/animal_traits.csv", header = TRUE)
morph_data <- animal_traits %>% select(c('genus', 'body.mass', 'metabolic.rate', 'mass.specific.metabolic.rate'))
morph_data <- ddply(morph_data,"genus",numcolwise(mean))
morph_data <- morph_data[complete.cases(morph_data),]

chk <- name.check(tree, morph_data$genus, data.names = morph_data$genus)
summary(chk)

morph_tree <- drop.tip(tree, chk$tree_not_data)
morph_data <- morph_data[(morph_data$genus %in% morph_tree$tip.label),]
name.check(morph_tree, morph_data$genus, data.names = morph_data$genus)

#remove dublicates
morph_tree <- drop.tip(morph_tree, c("Tetramorium", "Lasius", "Pimelia", "Aphaenogaster"))
morph_data <- morph_data[(morph_data$genus %in% morph_tree$tip.label),]

spp<-morph_data$genus
corBM<-corBrownian(phy = morph_tree, form = ~spp)
pgls <- gls(log(body.mass)~log(metabolic.rate), data = morph_data, correlation = corBM)
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
pgls <- gls(T.C~G.A + log(mass.specific.metabolic.rate) + log(body.mass) + log (metabolic.rate), data = t, correlation = corBM)
summary(pgls)


####GENUS_COLAPSE####
#brain_size
nwk_path <- "/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/MIDORI/genus_collapsed_co1seqs_sp_names.aln.treefile"
tree <- read.tree(file=nwk_path)

#edit cols to fit our data
animal_traits <- read.csv(file="/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/EcoMorphDBs/animal_traits.csv", header = TRUE)
morph_data <- animal_traits %>% select(c('genus', 'brain.size'))
morph_data <- ddply(morph_data,"genus",numcolwise(mean))
morph_data <- morph_data[complete.cases(morph_data),]

chk <- name.check(tree, morph_data$genus, data.names = morph_data$genus)
summary(chk)

morph_tree <- drop.tip(tree, chk$tree_not_data)
morph_data <- morph_data[(morph_data$genus %in% morph_tree$tip.label),]
name.check(morph_tree, morph_data$genus, data.names = morph_data$genus)

#remove dublicates
morph_tree <- drop.tip(morph_tree, c("Tetramorium", "Lasius", "Pimelia", "Aphaenogaster", "Cyphomyrmex", "Pheidole"))
morph_data <- morph_data[(morph_data$genus %in% morph_tree$tip.label),]

spp<-morph_data$genus
corBM<-corBrownian(phy = morph_tree, form = ~spp)
pgls <- gls(log(brain.size)~log(metabolic.rate), data = morph_data, correlation = corBM)
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
pgls <- gls(T.C~G.A + log(brain.size), data = t, correlation = corBM)
summary(pgls)
