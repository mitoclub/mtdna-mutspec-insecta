#courtesy of "Phylogenetic comparative methods in R"
library(ape)
library(geiger)
library(nlme)
library(plyr)
library(dplyr)

nwk_path <- "/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/MIDORI/all_insects_co1seqs(sp_names).aln.treefile"
tree <- read.tree(file=nwk_path)

#edit species col to fit our data
global_ants <- read.csv(file="/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/EcoMorphDBs/global_ants.csv", header = TRUE)
global_ants$Species <- paste(global_ants$Genus, global_ants$Species, sep = "_")
global_ants[global_ants == "Not in paper"] <- NA

#pruning tree and db
chk <- name.check(tree, global_ants$Species, data.names = global_ants$Species)
summary(chk)

ecomorph_tree <- drop.tip(tree, chk$tree_not_data)
global_ants <- global_ants[(global_ants$Species %in% ecomorph_tree$tip.label),]
name.check(ecomorph_tree, global_ants$Species, data.names = global_ants$Species)
length(ecomorph_tree$tip.label)
ecomorph_tree

#merging duplicates in db and removing non-morphological data
morph_data <- global_ants[, -which(names(global_ants) %in% c("Visibility","Traits.ID","Invasive.To.Political.Region1", "Morphospecies","Latitude","Longitude","Caste","Polymorphism","Queen.number","Colony.type","Colony.founding","Nest.Site","Activity","Diet","Political.Region1","Continent"))]

morph_data$Genus <- NULL
morph_data[,2:13] <- apply(morph_data[,2:13], 2, function(x) as.numeric(as.character(x)))
morph_data <- ddply(morph_data,"Species",numcolwise(mean))



#removing some cols to minimize loss due to NAs, still got only 34 species. Can we do some BOOTSTRAP?
morph_data <- morph_data[, -which(names(morph_data) %in% c("Whole.body.length..mm.", "Number.of.Spines...."))]
morph_data <- morph_data[complete.cases(morph_data),]


#plus extracting genuses for PCA coloring. JUST GET GENUSES FROM SP NAMES. How?
morph_genuses <- sub('_.*','',morph_data$Species)



####PCA####
chk <- name.check(ecomorph_tree, morph_data$Species, data.names = morph_data$Species)
morph_tree <- drop.tip(ecomorph_tree, chk$tree_not_data)
name.check(morph_tree, morph_data$Species, data.names = morph_data$Species)
rownames(morph_data) <- morph_data$Species
morph_data$Species <- NULL
morph_pca <- phyl.pca(morph_tree, morph_data)
morph_pca

par(mar=c(4.1,4.1,2.1,1.1),las=1) ## set margins
plot(morph_pca,main="")

morph_pca$Evec[,1]<--morph_pca$Evec[,1]
morph_pca$L[,1]<--morph_pca$L[,1]
morph_pca$S<-scores(morph_pca,
                       newdata=morph_data)

par(cex.axis=0.8,mar=c(5.1,5.1,1.1,1.1))
phylomorphospace(morph_tree,
                 scores(morph_pca)[,1:2],
                 ftype="off",node.size=c(0,1),bty="n",las=1,
                 xlab="PC1",
                 ylab=expression(paste("PC2")))
eco<-setNames(morph_genuses,rownames(morph_data))
ECO<-to.matrix(eco,unique(morph_genuses))
tiplabels(pie=ECO[morph_tree$tip.label,],cex=0.5)
legend(x="bottomright",legend=unique(morph_genuses),cex=0.8,pch=21,
       pt.bg=rainbow(n=length(unique(morph_genuses))),pt.cex=1.5, text.width = .7)


####PIC####
clypeus_length <- setNames(morph_data[,"Clypeus.length..mm."], rownames(morph_data))
max_eye_width <- setNames(morph_data[,"Max.eye.width..mm."], rownames(morph_data))

pic_clypeus <- pic(log(clypeus_length), morph_tree)
pic_eye <- pic(log(max_eye_width), morph_tree)

fit_pic <- lm(pic_clypeus~pic_eye+0)

## set margins
par(mar=c(5.1,5.1,1.1,1.1))
## graph scatterplot of contrasts
plot(pic_clypeus~pic_eye,
     xlab="PICs for log(max_eye_width)",
     ylab="PICs for log(clypeus_length)",
     pch=21,bg="gray",cex=1.2,las=1,
     cex.axis=0.7,cex.lab=0.9,bty="n")
## add gridlines to the plot
abline(h=0,lty="dotted")
abline(v=0,lty="dotted")
## reset graphing limits of the plot to the
## x/y range of our PICs
clip(min(pic_eye),max(pic_eye),
     min(pic_clypeus),max(pic_clypeus))
## graph our fitted line
abline(fit_pic,lwd=2,col="darkgray")

####PGLS####

spp<-rownames(morph_data)
corBM<-corBrownian(phy = morph_tree, form = ~spp)
pgls <- gls(log(Pronotum.width..mm.)~log(Max.eye.width..mm.), data = morph_data, correlation = corBM)
anova(pgls)

####TESTING####
ms12_internal <- read.csv(file = "/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/NemuPipeline/AllInsects/midori_all_insects_mutspec/mutspec_as_header_12internal.csv", header = TRUE)
colnames(ms12_internal)[1] = "Species"
ms12_internal <- ms12_internal[(ms12_internal$Species %in% morph_tree$tip.label),]
chk <- name.check(morph_tree, ms12_internal$Species, data.names = ms12_internal$Species)
morph_tree <- drop.tip(morph_tree, chk$tree_not_data)
name.check(morph_tree, ms12_internal$Species, data.names = ms12_internal$Species)
morph_data <- morph_data[(rownames(morph_data) %in% morph_tree$tip.label),]
t <- cbind(morph_data, ms12_internal[, -which(names(ms12_internal) %in% c("Species"))])

spp<-rownames(morph_data)
corBM<-corBrownian(phy = morph_tree, form = ~spp)
pgls <- gls(T.C~G.A + log(Pronotum.width..mm.), data = t, correlation = corBM)
anova(pgls)




                                     