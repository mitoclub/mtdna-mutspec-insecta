rm(list=ls(all=TRUE))
library(dplyr)
library(stringr)


############ Syn mut
SynNuc = read.table("/mnt/data/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/DescriptiveStat/codontable_midori_Blattodea.csv", header = TRUE, sep = ',')
#keeping code as is
colnames(SynNuc)[colnames(SynNuc) == 'Species_name'] <- 'Species'
SynNuc <- SynNuc %>% mutate(Species = str_replace(Species, "_[0-9]*$", ""))

names(SynNuc)

### make neg genes complementary:
neg_genes = subset(SynNuc, (Gene_name %in% c('ND1', 'ND4L', 'ND4', 'ND5')))
pos_genes = subset(SynNuc, !(Gene_name %in% c('ND1', 'ND4L', 'ND4', 'ND5')))

A = neg_genes$neutralT
T = neg_genes$neutralA
G = neg_genes$neutralC
C = neg_genes$neutralG
neg_genes$neutralA = A
neg_genes$neutralT = T
neg_genes$neutralG = G
neg_genes$neutralC = C
SynNuc = rbind(pos_genes,neg_genes)

#Temp changing taxonomy to class to keep code as is, but technicaly it's not a Class but a Family
colnames(SynNuc)[colnames(SynNuc) == 'Taxonomy'] <- 'Class'

VecOfTaxa = unique(SynNuc$Class)
table(SynNuc$Class)/13

########################################## GENOME WIDE SKEW FOR EACH SPECIES

AGG = aggregate(list(SynNuc$neutralA,SynNuc$neutralT,SynNuc$neutralG,SynNuc$neutralC), by = list(SynNuc$Species,SynNuc$Class), FUN = sum)
names(AGG) = c('Species','Class','neutralA','neutralT','neutralG','neutralC')

### count fraction of nucleotides
AGG$FrA = AGG$neutralA / (AGG$neutralA + AGG$neutralT + AGG$neutralG + AGG$neutralC)
AGG$FrT = AGG$neutralT / (AGG$neutralA + AGG$neutralT + AGG$neutralG + AGG$neutralC) 
AGG$FrG = AGG$neutralG / (AGG$neutralA + AGG$neutralT + AGG$neutralG + AGG$neutralC) 
AGG$FrC = AGG$neutralC / (AGG$neutralA + AGG$neutralT + AGG$neutralG + AGG$neutralC) 

## all six different skews
AGG$CTSkew = (AGG$neutralC - AGG$neutralT)/(AGG$neutralC + AGG$neutralT); summary(AGG$CTSkew) # GA on heavy
AGG$CGSkew = (AGG$neutralC - AGG$neutralG)/(AGG$neutralC + AGG$neutralG); summary(AGG$CGSkew) # 
AGG$CASkew = (AGG$neutralC - AGG$neutralA)/(AGG$neutralC + AGG$neutralA); summary(AGG$CASkew) # 
AGG$TGSkew = (AGG$neutralT - AGG$neutralG)/(AGG$neutralT + AGG$neutralG); summary(AGG$TGSkew) # 
AGG$TASkew = (AGG$neutralT - AGG$neutralA)/(AGG$neutralT + AGG$neutralA); summary(AGG$TASkew) # 
AGG$GASkew = (AGG$neutralG - AGG$neutralA)/(AGG$neutralG + AGG$neutralA); summary(AGG$GASkew) # 

AGG$TCSkew = (AGG$neutralT - AGG$neutralC)/(AGG$neutralT + AGG$neutralC); summary(AGG$CTSkew) # AG on heavy. Added it for fun, just to be sure, that it is opposite to CT (GA on heavy) 

#######opening MutSpec data
MUT = read.table('/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/NemuPipeline/AllInsects/midori_all_insects_mutspec/mutspec_as_header_12internal.csv', header = TRUE, sep=',')
colnames(MUT) <- c("Species", "A_C", "A_G", "A_T", "C_A", "C_G", "C_T", "G_A", "G_C", "G_T", "T_A", "T_C", "T_G")
MutCTskew = merge(MUT, AGG)
######mut spec with nuc content analyses
cor.test(MutCTskew$A_G, MutCTskew$GASkew, method = "spearman") # p > 0.05
cor.test(MutCTskew$C_T, MutCTskew$TCSkew, method = "spearman") # p > 0.05


#add workers designation
SynNuc = read.table("/mnt/data/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/DescriptiveStat/codontable_midori_Blattodea.csv", header = TRUE, sep = ',')
SynNuc = subset(SynNuc, !(Gene_name %in% c('ND1', 'ND4L', 'ND4', 'ND5')))
#fill 0 (terms w/o workers) to be 1
SynNuc$Workers[SynNuc$Workers == 0] <- 1
#fill na (cocks) to be 0
SynNuc <- SynNuc %>% mutate(Workers = ifelse(is.na(Workers), 0, Workers))


colnames(SynNuc)[colnames(SynNuc) == 'Species_name'] <- 'Species'
SynNuc <- SynNuc %>% mutate(Species = str_replace(Species, "_[0-9]*$", ""))
workers = SynNuc[!duplicated(SynNuc[c('Species')]), ]

Mam = merge(AGG,workers[c("Species", "Workers")], by ='Species')
MutCTskew = merge(MutCTskew, workers[c("Species", "Workers")])
summary(lm(GASkew ~ A_G, data=MutCTskew))
summary(lm(A_G ~ FrA + FrG, data=MutCTskew))
cor.test(MutCTskew$A_G, MutCTskew$FrA, method = "spearman")
cor.test(MutCTskew$A_G, MutCTskew$FrA, method = "spearman")
MutCTskew$FrGFrA = MutCTskew$FrG/MutCTskew$FrA
cor.test(MutCTskew$A_G, MutCTskew$FrGFrA, method = "spearman")

############### MIGHT BE WHAT I NEED! organism type versus all six skews: 6 rank correlations and one multiple model:
###### pairwise rank corr:
cor.test((Mam$Workers),Mam$CTSkew, method = 'spearman') # p < 0.05, rho = 0.6391
cor.test((Mam$Workers),Mam$CGSkew, method = 'spearman')  
cor.test((Mam$Workers),Mam$CASkew, method = 'spearman') 
cor.test((Mam$Workers),Mam$TGSkew, method = 'spearman')  
cor.test((Mam$Workers),Mam$TASkew, method = 'spearman')  
cor.test((Mam$Workers),Mam$GASkew, method = 'spearman') # p < 0.05, rho = 0.4816

######### multiple Linear Model (it is reasonable to add only 4 significant parameters - without CG and TA. But If I add everything, results are unusual and close to opposite...). Interactions?
A<-lm((Mam$Workers) ~ Mam$CTSkew + Mam$CGSkew + Mam$CASkew + Mam$TGSkew + Mam$TASkew + Mam$GASkew); summary(A)
#(Intercept)   0.7460     0.2600   2.870  0.00433 ** 
#  Mam$CTSkew    2.0415     0.3694   5.526 6.00e-08 *** --- good
#  Mam$CGSkew   -2.0368     0.2036 -10.002  < 2e-16 ***
#  Mam$CASkew    0.6124     0.7925   0.773  0.44014    
#  Mam$TGSkew    1.8636     0.2694   6.918 1.89e-11 ***
#  Mam$TASkew   -0.7394     0.3280  -2.255  0.02471 *  
#  Mam$GASkew   -0.1538     0.8269  -0.186  0.85252     --- trash

A<-lm((Mam$Workers) ~ Mam$CTSkew + Mam$GASkew); summary(A)
A<-lm((Mam$Workers) ~ Mam$CTSkew +  Mam$CASkew + Mam$TGSkew + Mam$GASkew); summary(A)
#this one is fucking greate
#  (Intercept)   1.9987     0.2431   8.221 2.98e-15 ***
#  Mam$CTSkew    2.3047     0.1871  12.318  < 2e-16 ***
#  Mam$CASkew   -3.7783     0.3633 -10.400  < 2e-16 ***
#  Mam$TGSkew    1.5271     0.2328   6.560 1.70e-10 ***
#  Mam$GASkew    5.1434     0.6058   8.490 4.33e-16 ***  


####WIP####
##################### PICs

library(ape)
library(geiger)
library(caper)

nwk_path <- "/home/gabs/Documents/lab/TermitesAndCockroaches/mtdna-mutspec-insecta/data/MIDORI/all_insects_co1seqs_sp_names.aln.treefile"

tree <- read.tree(file=nwk_path)

row.names(Mam) = Mam$Species

#Error in if (warnings) { : the condition has lenWorkersh > 1 --- Solved by setting sort and warnings to False, not sure if that's good
tree_w = treedata(tree, Mam, sort=F, warnings=F)$phy
#Error in if (warnings) { : the condition has lenWorkersh > 1 --- Solved by setting sort and warnings to False, not sure if that's good
data<-as.data.frame(treedata(tree_w, Mam, sort=F, warnings=F)$data)

data$Species = as.character(data$Species)

data$CTSkew = as.numeric(as.character(data$CTSkew))
data$Workers = as.numeric(as.character(data$Workers))

cor.test(pic(data$CTSkew, tree_w), pic(data$Workers, tree_w), method = 'spearman')

# rho 
# 0.926
# p-value = 8.197e-14

MutComp = comparative.data(tree_w, data, Species, vcv=TRUE)

model = pgls(scale(Workers) ~ scale(CTSkew), MutComp, lambda="ML")
summary(model)

#  lambda [ ML]  : 1.000
#  Estimate Std. Error t value  Pr(>|t|)    
#  (Intercept)   -0.45315    0.16638 -2.7236 0.0106609 *  
#  scale(CTSkew)  0.46127    0.11919  3.8699 0.0005452 ***

crunch(scale(Workers) ~ scale(CTSkew), MutComp)

# p-value: 0.0005452

################################### more PICs
library(ape)
library(geiger)
library(caper)

#MutCTskew = MutCTskew[-156,]
MutCTskew$FrCFrT = MutCTskew$FrC/MutCTskew$FrT
tree <- read.tree(nwk_path)
row.names(MutCTskew) = MutCTskew$Species
#Error in if (warnings) { : the condition has lenWorkersh > 1 --- Solved by setting sort and warnings to False, not sure if that's good
tree_w = treedata(tree, MutCTskew, sort=F, warnings=F)$phy

#Error in if (warnings) { : the condition has lenWorkersh > 1 --- Solved by setting sort and warnings to False, not sure if that's good
data<-as.data.frame(treedata(tree_w, MutCTskew, sort=F, warnings=F)$data)
data$Species = as.character(data$Species)
data$T_C = as.numeric(as.character(data$T_C))
data$FrT = as.numeric(as.character(data$FrT))
data$FrC = as.numeric(as.character(data$FrC))
data$FrCFrT = as.numeric(as.character(data$FrCFrT))
data$CTSkew = as.numeric(as.character(data$CTSkew))

MutComp = comparative.data(tree_w, data, Species, vcv=TRUE)
summary(pgls(T_C ~ FrT + FrC, MutComp, lambda="ML"))
summary(pgls(T_C ~ FrCFrT, MutComp, lambda="ML"))
summary(pgls(T_C ~ FrC, MutComp, lambda="ML"))
summary(pgls(T_C ~ FrT, MutComp, lambda="ML"))
summary(pgls(T_C ~ FrT * FrC, MutComp, lambda="ML"))
summary(pgls(T_C ~ CTSkew, MutComp, lambda="ML"))

str(MutComp)

#################################### GENE_SPECIFIC SKEW FOR EACH SPECIES
colnames(SynNuc)[colnames(SynNuc) == 'Taxonomy'] <- 'Class'
AGG = aggregate(list(SynNuc$neutralA,SynNuc$neutralT,SynNuc$neutralG,SynNuc$neutralC), by = list(SynNuc$Species,SynNuc$Class,SynNuc$Gene_name, SynNuc$Workers), FUN = sum)
names(AGG) = c('Species','Class','Gene_name','Workers', 'neutralA','neutralT','neutralG','neutralC')
AGG$CTSkew = (AGG$neutralC - AGG$neutralT)/(AGG$neutralC + AGG$neutralT); summary(AGG$TCSkew) # GA on heavy

M =AGG

M$Gene =  ordered(M$Gene, levels = c('ND2', 'CO1','CO2','A8','A6','CO3','ND3','ND6','Cytb'))
# M$Gene =  ordered(M$Gene, levels = c('COX1','COX2','ATP8','ATP6','COX3','ND3','ND4L','ND4','ND5','CytB'))
M = M[order(M$Gene),]

# par(mfrow=c(2,1), oma = c(3, 1, 1, 1), cex = 2)
boxplot(CTSkew ~ Workers*Gene, data = M,  notch = TRUE, outline = FALSE, las = 2, col = c('red','green'), main = 'Blattodea, GA skew')

####WHAT I NEED######
####### naive multiple linear model, assigning numbers to order of genes. We can improve it substituting integers by real time or using more correct stat
FromGenesToNumbers = data.frame(c('ND2', 'CO1','CO2','A8','A6','COX3','ND3', 'ND6','Cytb'),seq(1:9)); names(FromGenesToNumbers)=c('Gene','TSSS') # time spend single stranded
names(FromGenesToNumbers)=c('Gene','TBSS') # time spend single stranded
gl_and_tbss = merge(M,FromGenesToNumbers) 

summary(lm(gl_and_tbss$CTSkew ~ (gl_and_tbss$Workers) + gl_and_tbss$TBSS))
summary(lm(gl_and_tbss$CTSkew ~ (gl_and_tbss$Workers) * gl_and_tbss$TBSS))



####beautiful boxplots###################

library("ggpubr")

MM = M
MM$Workers<- factor(MM$Workers, levels = c(0, 1))
ggboxplot(MM, "Gene", "CTSkew",
          fill = "Workers", palette = c("#22AF1D", "#0B350B"), xlab="Genes", ylab="GA skew", title = "GA skew in cock vs ter", legend.title = "Blattodea", width = 0.7, notch = TRUE)
dev.off()