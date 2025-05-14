## Export variables and load libraries
rm(list=ls())
library(tidyverse)

options(scipen=999)

args <- strsplit(commandArgs(trailingOnly=T), "=", fixed=T)
for(i in 1:length(args)) {
  assign(args[[i]][1], args[[i]][2])
}

### Import data
data <- read.delim(GENES, head=T, skipNul = T)
data_pos <- read.delim(paste(WORKDIR1, "/D7_positive_selection/metadata/ABSREL_results.tsv", sep=""), head=T, skipNul = T)
data_rel1 <- read.delim(paste(WORKDIR1, "/D8_relaxed_selection1/metadata/relax_results.tsv", sep=""), head=T, skipNul = T)
data_rel2 <- read.delim(paste(WORKDIR1, "/D8_relaxed_selection2/metadata/relax_results.tsv", sep=""), head=T, skipNul = T)

### Merge data

## aBSREL
colnames(data_pos) <- c("geneID", "SexLinked_aBSREL_p.Z1", "SexLinked_aBSREL_p.W1", "SexLinked_aBSREL_p.A1",
                        "SexLinked_aBSREL_p.Z2", "SexLinked_aBSREL_p.W2", "SexLinked_aBSREL_p.A2",
                        "SexLinked_aBSREL_p.Z3", "SexLinked_aBSREL_p.W3", "SexLinked_aBSREL_p.A3")
data_pos <- data_pos[,c("geneID", "SexLinked_aBSREL_p.A1", "SexLinked_aBSREL_p.Z1", "SexLinked_aBSREL_p.W1",
                        "SexLinked_aBSREL_p.A2", "SexLinked_aBSREL_p.Z2", "SexLinked_aBSREL_p.W2",
                        "SexLinked_aBSREL_p.A3", "SexLinked_aBSREL_p.Z3", "SexLinked_aBSREL_p.W3")]
data_pos[,2:ncol(data_pos)] <- lapply(data_pos[,2:ncol(data_pos)], as.numeric)
data_pos$Species <- "Skylark"                  
                        
data1 <- merge(data, data_pos, by = c("geneID", "Species"), all.x = TRUE, sort = FALSE)
data <- data1[c(colnames(data), setdiff(colnames(data1), names(data)))]

## RELAX1
colnames(data_rel1) <- c("geneID", "SexLinked_RELAX1_K.WvsZ", "SexLinked_RELAX1_p.WvsZ1", "SexLinked_RELAX1_p.WvsZ2")    
data_rel1[,2:ncol(data_rel1)] <- lapply(data_rel1[,2:ncol(data_rel1)] , as.numeric)
data_rel1$Species <- "Skylark"

data1 <- merge(data, data_rel1, by = c("geneID", "Species"), all.x = TRUE,  sort = FALSE)
data <- data1[c(colnames(data), setdiff(colnames(data1), names(data)))]

## RELAX2
colnames(data_rel2) <- c("geneID", "SexLinked_RELAX2_K.ZWvsA", "SexLinked_RELAX2_p.ZWvsA1", "SexLinked_RELAX2_p.ZWvsA2")       
data_rel2[,2:ncol(data_rel2)] <- lapply(data_rel2[,2:ncol(data_rel2)] , as.numeric)
data_rel2$Species <- "Skylark"

data1 <- merge(data, data_rel2, by = c("geneID", "Species"), all.x = TRUE, sort = FALSE)
data <- data1[c(colnames(data), setdiff(colnames(data1), names(data)))]

## Write organised table
write.table(data, paste(OUTDIR, "/", PROJECT1, "_", PROJECT2, "_organised_data2.tsv", sep=""), quote=F, sep='\t', row.names = F, col.names = T)
q()
