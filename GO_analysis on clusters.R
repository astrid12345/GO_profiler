#################################
# Name: A. Alsema
# Date: 2021 Sept
# Purpose: run clusterProfiler GO analysis on enriched cluster markers
#################################

# clear workspace
rm(list = ls())
# configuration
library(enrichplot)
library(ggplot2)
library(org.Hs.eg.db)
library(clusterProfiler)
setwd("C:/Users/astri/Documents/UMCG/Schizo/for_analysis/")
# load data
mylist <- readRDS("subclusters_upmarkers_list.rds")

# Perform Gene Ontology (GO) enrichment analysis
compareGO <- compareCluster(geneClusters =   mylist, fun = "enrichGO", 
                                   OrgDb =   "org.Hs.eg.db",
                                   keyType = "SYMBOL",
                                   ont   =   "BP")
p <- dotplot(compareGO,showCategory = 5, font.size = 5) 
p <- p + scale_colour_gradientn(colours = c("darkred","#fee0d2")) + 
         xlab("") +
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Save the dot plot as a TIFF file
tiff("Routput/GO_subclustermarkers.tiff", width = 7, height = 7, units = 'in', res = 300)
p
dev.off()

# export GO table
GOinfo <- as.data.frame(compareGO)
write.csv(file = "Routput/GO_Tables_neuronalsubclusters.csv", GOinfo)

# select the top 5 most significant per cluster
GOtop <- by(GOinfo, GOinfo["Cluster"], head, n = 5) 
GOtop <- Reduce(rbind, GOtop)
write.csv(GOtop, file = "./Routput/Top5_GO_Table_neuronalsubclusters.csv" )