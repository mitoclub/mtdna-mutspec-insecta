library(ggplot2)
library(ggpubr)
library(tidyverse)

library(GGally)
library(ggExtra)
library(ggalluvial)
library(plotly)


df <- read.csv("/home/gab/Documents/lab/TermitesAndCockroaches/MutSpec-Redone/interim/MIDORI/midori_sp_table_Fam.csv")

at_least_5 <- subset(df, CO1>=5)

p <- ggplot(at_least_5, aes(x=CO1)) + 
  geom_histogram(color="black", fill="white") + labs(y='')

p + scale_y_continuous(trans='log2') 
