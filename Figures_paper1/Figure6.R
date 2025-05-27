## Export variables and load libraries
rm(list=ls())
library(tidyverse)
library(ggtree)
library(tidytree)
library(ape)
library(jsonlite)
library(patchwork)
library(viridis)
sessionInfo()
#library(phytools)
#library(treeio)

options(scipen=999)
setwd("C:/Users/Simon JE/OneDrive - Lund University/Dokument/Simon/PhD/Projects/Skylark_2021/Results/")

# Import data
data <- read.delim("Genes/Skylark_2021_Rasolark_2021_organised_data2.tsv", sep="\t", head=T)
data1 <- data[which(data$Filter1=="OK" & data$Filter3=="OK" & data$Filter4=="OK" & data$Filter5=="OK"),]
data1 <- data1[which(data1$Strata != "S0" & data1$Strata != "S1" & data1$Strata != "S2" & data1$Strata != "S3" & data1$Strata != "Ancestral unknown"),]
data1 <- data1[which(data1$Region == "sex_phase"),]
genes <- names(which(table(data1$geneID) == 2))
data1 <- data1[data1$geneID %in% genes,]
data1 <- data1[which(data1$Species == "Skylark"),]

pos_test <- data[which(data1$SexLinked_aBSREL_p.A2 < 0.05 | data1$SexLinked_aBSREL_p.Z2 < 0.05 |data1$SexLinked_aBSREL_p.W2 < 0.05),c("geneID", "Strata", "SexLinked_aBSREL_p.Z2", "SexLinked_aBSREL_p.W2", "SexLinked_aBSREL_p.A2" )]

for(i in 1:nrow(pos_test)) {
  gene <- pos_test$geneID[i]
  data_json <- fromJSON(paste("Gene_trees/D7_positive_selection/genes/", gene, "/Skylark_2021_Rasolark_2021_", gene, "_aBSREL.json", sep=""), simplifyVector = FALSE)
  branch_data <- data_json$`branch attributes`$`0`
  print(c(gene, unlist(branch_data$BranchW$`Rate Distributions`)))
}


# Import tree
tree <- read.tree("Gene_trees/Skylark_2021_Rasolark_2021_sex_phase_sex_phase_females.nwk")
tree <- as_tibble(tree)
tree$label <- gsub("-", "_", tree$label)
tree$branch.length <- 0
tree$label[43] <- "Node1"

## All genes
data2 <- data1
tree_all <- tree
genes <- data1$geneID
Branch_length_matrix <- matrix(NA, length(tree_all$label), length(genes))
colnames(Branch_length_matrix) <- genes
rownames(Branch_length_matrix) <- tree_all$label
MLE_matrix <- matrix(NA, length(tree_all$label), length(genes))
colnames(MLE_matrix) <- genes
rownames(MLE_matrix) <- tree_all$label
gene_length <- rep(NA, length(genes))
names(gene_length) <- genes

for(i in 1:length(genes)) {
  data_json <- fromJSON(paste("Gene_trees/sex_phase_sex_phase/", genes[i], "/Skylark_2021_Rasolark_2021_", genes[i], "_MG94.json", sep=""), simplifyVector = FALSE)
  branch_data <- data_json$`branch attributes`$`0`
  Branch_length <- sapply(branch_data, function(x) x$nonsynonymous +  x$synonymous)
  
  Branch_length_matrix[match(names(Branch_length), rownames(Branch_length_matrix)), which(colnames(Branch_length_matrix) == genes[i])] <- Branch_length
  MLE <- sapply(branch_data, function(x) x$`Confidence Intervals`$MLE)
  MLE_matrix[match(names(MLE), rownames(Branch_length_matrix)), which(colnames(Branch_length_matrix) == genes[i])] <- MLE
  gene_length[i] <- data1$geneLengthTrimmed[i]
}

Branch_length_matrix <- Branch_length_matrix[, colSums(!is.na(Branch_length_matrix)) > 0]
MLE_matrix <- MLE_matrix[, colSums(!is.na(MLE_matrix)) > 0]
gene_length <- as.numeric(gene_length[!is.na(gene_length)])

Branch_length_matrix_length <- sweep(Branch_length_matrix, 2, gene_length, `*`)
MLE_matrix_length <- sweep(MLE_matrix, 2, gene_length, `*`)

Branch_length_median <- rep(NA, length(tree_all$label))
names(Branch_length_median) <- rownames(Branch_length_matrix_length)
for(i in 1:length(Branch_length_median)) {
  Branch_length_median[i] <- median(Branch_length_matrix_length[i,], na.rm=T)/median(gene_length)
}

MLE_median <- rep(NA, length(tree_all$label))
names(MLE_median) <- rownames(MLE_matrix_length)
for(i in 1:length(MLE_median)) {
  MLE_median[i] <- median(MLE_matrix_length[i,], na.rm=T)/median(gene_length)
}

# Map tree data to tree
tree_all$branch.length <- Branch_length_median[match(tree_all$label, names(Branch_length_median))]
tree_all$branch.length[is.na(tree_all$branch.length)] <- 0
tree_all$Omega <- MLE_median[match(tree_all$label, names(MLE_median))]
tree_all$Omega[which(tree_all$label == "")] <- tree_all$Omega[which(tree_all$label == "BranchO")]
tree_all$log10Omega <- log10(tree_all$Omega)
node_data_all <- tree_all[,c("label", "Omega", "log10Omega")]
node_data_all2 <- node_data_all[which(node_data_all$label == "BranchA" | node_data_all$label == "BranchZ" | node_data_all$label == "BranchW"),]
tree_all <- as.phylo(tree_all)


## 4A
data2 <- data1[which(data1$Strata == "4A"),]
tree_4A <- tree
genes <- data2$geneID
Branch_length_matrix <- matrix(NA, length(tree_4A$label), length(genes))
colnames(Branch_length_matrix) <- genes
rownames(Branch_length_matrix) <- tree_4A$label
MLE_matrix <- matrix(NA, length(tree_4A$label), length(genes))
colnames(MLE_matrix) <- genes
rownames(MLE_matrix) <- tree_4A$label
gene_length <- rep(NA, length(genes))
names(gene_length) <- genes


for(i in 1:length(genes)) {
  data_json <- fromJSON(paste("Gene_trees/sex_phase_sex_phase/", genes[i], "/Skylark_2021_Rasolark_2021_", genes[i], "_MG94.json", sep=""), simplifyVector = FALSE)
  branch_data <- data_json$`branch attributes`$`0`
  Branch_length <- sapply(branch_data, function(x) x$nonsynonymous +  x$synonymous)
  
  Branch_length_matrix[match(names(Branch_length), rownames(Branch_length_matrix)), which(colnames(Branch_length_matrix) == genes[i])] <- Branch_length
  MLE <- sapply(branch_data, function(x) x$`Confidence Intervals`$MLE)
  MLE_matrix[match(names(MLE), rownames(Branch_length_matrix)), which(colnames(Branch_length_matrix) == genes[i])] <- MLE
  gene_length[i] <- data2$geneLengthTrimmed[i]
}

Branch_length_matrix <- Branch_length_matrix[, colSums(!is.na(Branch_length_matrix)) > 0]
MLE_matrix <- MLE_matrix[, colSums(!is.na(MLE_matrix)) > 0]
gene_length <- as.numeric(gene_length[!is.na(gene_length)])

Branch_length_matrix_length <- sweep(Branch_length_matrix, 2, gene_length, `*`)
MLE_matrix_length <- sweep(MLE_matrix, 2, gene_length, `*`)

Branch_length_median <- rep(NA, length(tree_4A$label))
names(Branch_length_median) <- rownames(Branch_length_matrix_length)
for(i in 1:length(Branch_length_median)) {
  Branch_length_median[i] <- median(Branch_length_matrix_length[i,], na.rm=T)/median(gene_length)
}

MLE_median <- rep(NA, length(tree_4A$label))
names(MLE_median) <- rownames(MLE_matrix_length)
for(i in 1:length(MLE_median)) {
  MLE_median[i] <- median(MLE_matrix_length[i,], na.rm=T)/median(gene_length)
}

# Map tree data to tree
tree_4A$branch.length <- Branch_length_median[match(tree_4A$label, names(Branch_length_median))]
tree_4A$branch.length[is.na(tree_4A$branch.length)] <- 0
tree_4A$Omega <- MLE_median[match(tree_4A$label, names(MLE_median))]
tree_4A$Omega[which(tree_4A$label == "")] <- tree_4A$Omega[which(tree_4A$label == "BranchO")]
tree_4A$log10Omega <- log10(tree_4A$Omega)
node_data_4A <- tree_4A[,c("label", "Omega", "log10Omega")]
tree_4A <- as.phylo(tree_4A)


## 3a
data2 <- data1[which(data1$Strata == "3a"),]
tree_3a <- tree
genes <- data2$geneID
Branch_length_matrix <- matrix(NA, length(tree_3a$label), length(genes))
colnames(Branch_length_matrix) <- genes
rownames(Branch_length_matrix) <- tree_3a$label
MLE_matrix <- matrix(NA, length(tree_3a$label), length(genes))
colnames(MLE_matrix) <- genes
rownames(MLE_matrix) <- tree_3a$label
gene_length <- rep(NA, length(genes))
names(gene_length) <- genes


for(i in 1:length(genes)) {
  data_json <- fromJSON(paste("Gene_trees/sex_phase_sex_phase/", genes[i], "/Skylark_2021_Rasolark_2021_", genes[i], "_MG94.json", sep=""), simplifyVector = FALSE)
  branch_data <- data_json$`branch attributes`$`0`
  Branch_length <- sapply(branch_data, function(x) x$nonsynonymous +  x$synonymous)
  
  Branch_length_matrix[match(names(Branch_length), rownames(Branch_length_matrix)), which(colnames(Branch_length_matrix) == genes[i])] <- Branch_length
  MLE <- sapply(branch_data, function(x) x$`Confidence Intervals`$MLE)
  MLE_matrix[match(names(MLE), rownames(Branch_length_matrix)), which(colnames(Branch_length_matrix) == genes[i])] <- MLE
  gene_length[i] <- data2$geneLengthTrimmed[i]
}

Branch_length_matrix <- Branch_length_matrix[, colSums(!is.na(Branch_length_matrix)) > 0]
MLE_matrix <- MLE_matrix[, colSums(!is.na(MLE_matrix)) > 0]
gene_length <- as.numeric(gene_length[!is.na(gene_length)])

Branch_length_matrix_length <- sweep(Branch_length_matrix, 2, gene_length, `*`)
MLE_matrix_length <- sweep(MLE_matrix, 2, gene_length, `*`)

Branch_length_median <- rep(NA, length(tree_3a$label))
names(Branch_length_median) <- rownames(Branch_length_matrix_length)
for(i in 1:length(Branch_length_median)) {
  Branch_length_median[i] <- median(Branch_length_matrix_length[i,], na.rm=T)/median(gene_length)
}

MLE_median <- rep(NA, length(tree_3a$label))
names(MLE_median) <- rownames(MLE_matrix_length)
for(i in 1:length(MLE_median)) {
  MLE_median[i] <- median(MLE_matrix_length[i,], na.rm=T)/median(gene_length)
}

# Map tree data to tree
tree_3a$branch.length <- Branch_length_median[match(tree_3a$label, names(Branch_length_median))]
tree_3a$branch.length[is.na(tree_3a$branch.length)] <- 0
tree_3a$Omega <- MLE_median[match(tree_3a$label, names(MLE_median))]
tree_3a$Omega[which(tree_3a$label == "")] <- tree_3a$Omega[which(tree_3a$label == "BranchO")]
tree_3a$log10Omega <- log10(tree_3a$Omega)
node_data_3a <- tree_3a[,c("label", "Omega", "log10Omega")]
tree_3a <- as.phylo(tree_3a)


## 3b
data2 <- data1[which(data1$Strata == "3b"),]
tree_3b <- tree
genes <- data2$geneID
Branch_length_matrix <- matrix(NA, length(tree_3b$label), length(genes))
colnames(Branch_length_matrix) <- genes
rownames(Branch_length_matrix) <- tree_3b$label
MLE_matrix <- matrix(NA, length(tree_3b$label), length(genes))
colnames(MLE_matrix) <- genes
rownames(MLE_matrix) <- tree_3b$label
gene_length <- rep(NA, length(genes))
names(gene_length) <- genes


for(i in 1:length(genes)) {
  data_json <- fromJSON(paste("Gene_trees/sex_phase_sex_phase/", genes[i], "/Skylark_2021_Rasolark_2021_", genes[i], "_MG94.json", sep=""), simplifyVector = FALSE)
  branch_data <- data_json$`branch attributes`$`0`
  Branch_length <- sapply(branch_data, function(x) x$nonsynonymous +  x$synonymous)
  
  Branch_length_matrix[match(names(Branch_length), rownames(Branch_length_matrix)), which(colnames(Branch_length_matrix) == genes[i])] <- Branch_length
  MLE <- sapply(branch_data, function(x) x$`Confidence Intervals`$MLE)
  MLE_matrix[match(names(MLE), rownames(Branch_length_matrix)), which(colnames(Branch_length_matrix) == genes[i])] <- MLE
  gene_length[i] <- data2$geneLengthTrimmed[i]
}

Branch_length_matrix <- Branch_length_matrix[, colSums(!is.na(Branch_length_matrix)) > 0]
MLE_matrix <- MLE_matrix[, colSums(!is.na(MLE_matrix)) > 0]
gene_length <- as.numeric(gene_length[!is.na(gene_length)])

Branch_length_matrix_length <- sweep(Branch_length_matrix, 2, gene_length, `*`)
MLE_matrix_length <- sweep(MLE_matrix, 2, gene_length, `*`)

Branch_length_median <- rep(NA, length(tree_3b$label))
names(Branch_length_median) <- rownames(Branch_length_matrix_length)
for(i in 1:length(Branch_length_median)) {
  Branch_length_median[i] <- median(Branch_length_matrix_length[i,], na.rm=T)/median(gene_length)
}

MLE_median <- rep(NA, length(tree_3b$label))
names(MLE_median) <- rownames(MLE_matrix_length)
for(i in 1:length(MLE_median)) {
  MLE_median[i] <- median(MLE_matrix_length[i,], na.rm=T)/median(gene_length)
}

# Map tree data to tree
tree_3b$branch.length <- Branch_length_median[match(tree_3b$label, names(Branch_length_median))]
tree_3b$branch.length[is.na(tree_3b$branch.length)] <- 0
tree_3b$Omega <- MLE_median[match(tree_3b$label, names(MLE_median))]
tree_3b$Omega[which(tree_3b$label == "")] <- tree_3b$Omega[which(tree_3b$label == "BranchO")]
tree_3b$log10Omega <- log10(tree_3b$Omega)
node_data_3b <- tree_3b[,c("label", "Omega", "log10Omega")]
tree_3b <- as.phylo(tree_3b)


## 5
data2 <- data1[which(data1$Strata == "5"),]
tree_5 <- tree
genes <- data2$geneID
Branch_length_matrix <- matrix(NA, length(tree_5$label), length(genes))
colnames(Branch_length_matrix) <- genes
rownames(Branch_length_matrix) <- tree_5$label
MLE_matrix <- matrix(NA, length(tree_5$label), length(genes))
colnames(MLE_matrix) <- genes
rownames(MLE_matrix) <- tree_5$label
gene_length <- rep(NA, length(genes))
names(gene_length) <- genes


for(i in 1:length(genes)) {
  data_json <- fromJSON(paste("Gene_trees/sex_phase_sex_phase/", genes[i], "/Skylark_2021_Rasolark_2021_", genes[i], "_MG94.json", sep=""), simplifyVector = FALSE)
  branch_data <- data_json$`branch attributes`$`0`
  Branch_length <- sapply(branch_data, function(x) x$nonsynonymous +  x$synonymous)
  
  Branch_length_matrix[match(names(Branch_length), rownames(Branch_length_matrix)), which(colnames(Branch_length_matrix) == genes[i])] <- Branch_length
  MLE <- sapply(branch_data, function(x) x$`Confidence Intervals`$MLE)
  MLE_matrix[match(names(MLE), rownames(Branch_length_matrix)), which(colnames(Branch_length_matrix) == genes[i])] <- MLE
  gene_length[i] <- data2$geneLengthTrimmed[i]
}

Branch_length_matrix <- Branch_length_matrix[, colSums(!is.na(Branch_length_matrix)) > 0]
MLE_matrix <- MLE_matrix[, colSums(!is.na(MLE_matrix)) > 0]
gene_length <- as.numeric(gene_length[!is.na(gene_length)])

Branch_length_matrix_length <- sweep(Branch_length_matrix, 2, gene_length, `*`)
MLE_matrix_length <- sweep(MLE_matrix, 2, gene_length, `*`)

Branch_length_median <- rep(NA, length(tree_5$label))
names(Branch_length_median) <- rownames(Branch_length_matrix_length)
for(i in 1:length(Branch_length_median)) {
  Branch_length_median[i] <- median(Branch_length_matrix_length[i,], na.rm=T)/median(gene_length)
}

MLE_median <- rep(NA, length(tree_5$label))
names(MLE_median) <- rownames(MLE_matrix_length)
for(i in 1:length(MLE_median)) {
  MLE_median[i] <- median(MLE_matrix_length[i,], na.rm=T)/median(gene_length)
}

# Map tree data to tree
tree_5$branch.length <- Branch_length_median[match(tree_5$label, names(Branch_length_median))]
tree_5$branch.length[is.na(tree_5$branch.length)] <- 0
tree_5$Omega <- MLE_median[match(tree_5$label, names(MLE_median))]
tree_5$Omega[which(tree_5$label == "")] <- tree_5$Omega[which(tree_5$label == "BranchO")]
tree_5$log10Omega <- log10(tree_5$Omega)
node_data_5 <- tree_5[,c("label", "Omega", "log10Omega")]
tree_5 <- as.phylo(tree_5)


## 3c
data2 <- data1[which(data1$Strata == "3c"),]
tree_3c <- tree
genes <- data2$geneID
Branch_length_matrix <- matrix(NA, length(tree_3c$label), length(genes))
colnames(Branch_length_matrix) <- genes
rownames(Branch_length_matrix) <- tree_3c$label
MLE_matrix <- matrix(NA, length(tree_3c$label), length(genes))
colnames(MLE_matrix) <- genes
rownames(MLE_matrix) <- tree_3c$label
gene_length <- rep(NA, length(genes))
names(gene_length) <- genes


for(i in 1:length(genes)) {
  data_json <- fromJSON(paste("Gene_trees/sex_phase_sex_phase/", genes[i], "/Skylark_2021_Rasolark_2021_", genes[i], "_MG94.json", sep=""), simplifyVector = FALSE)
  branch_data <- data_json$`branch attributes`$`0`
  Branch_length <- sapply(branch_data, function(x) x$nonsynonymous +  x$synonymous)
  
  Branch_length_matrix[match(names(Branch_length), rownames(Branch_length_matrix)), which(colnames(Branch_length_matrix) == genes[i])] <- Branch_length
  MLE <- sapply(branch_data, function(x) x$`Confidence Intervals`$MLE)
  MLE_matrix[match(names(MLE), rownames(Branch_length_matrix)), which(colnames(Branch_length_matrix) == genes[i])] <- MLE
  gene_length[i] <- data2$geneLengthTrimmed[i]
}

Branch_length_matrix <- Branch_length_matrix[, colSums(!is.na(Branch_length_matrix)) > 0]
MLE_matrix <- MLE_matrix[, colSums(!is.na(MLE_matrix)) > 0]
gene_length <- as.numeric(gene_length[!is.na(gene_length)])

Branch_length_matrix_length <- sweep(Branch_length_matrix, 2, gene_length, `*`)
MLE_matrix_length <- sweep(MLE_matrix, 2, gene_length, `*`)

Branch_length_median <- rep(NA, length(tree_3c$label))
names(Branch_length_median) <- rownames(Branch_length_matrix_length)
for(i in 1:length(Branch_length_median)) {
  Branch_length_median[i] <- median(Branch_length_matrix_length[i,], na.rm=T)/median(gene_length)
}

MLE_median <- rep(NA, length(tree_3c$label))
names(MLE_median) <- rownames(MLE_matrix_length)
for(i in 1:length(MLE_median)) {
  MLE_median[i] <- median(MLE_matrix_length[i,], na.rm=T)/median(gene_length)
}

# Map tree data to tree
tree_3c$branch.length <- Branch_length_median[match(tree_3c$label, names(Branch_length_median))]
tree_3c$branch.length[is.na(tree_3c$branch.length)] <- 0
tree_3c$Omega <- MLE_median[match(tree_3c$label, names(MLE_median))]
tree_3c$Omega[which(tree_3c$label == "")] <- tree_3c$Omega[which(tree_3c$label == "BranchO")]
tree_3c$log10Omega <- log10(tree_3c$Omega)
node_data_3c <- tree_3c[,c("label", "Omega", "log10Omega")]
tree_3c <- as.phylo(tree_3c)

all_branches <- as.numeric(unlist(c(node_data_all$Omega, node_data_4A$Omega, node_data_3a$Omega, node_data_3b$Omega, node_data_5$Omega, node_data_3c$Omega)))
min_branch_val <- min(log10(all_branches[which(all_branches > 0)]))
node_data_all$log10Omega[which(node_data_all$Omega == 0)] <- min_branch_val
node_data_4A$log10Omega[which(node_data_4A$Omega == 0)] <- min_branch_val
node_data_3a$log10Omega[which(node_data_3a$Omega == 0)] <- min_branch_val
node_data_3b$log10Omega[which(node_data_3b$Omega == 0)] <- min_branch_val
node_data_5$log10Omega[which(node_data_5$Omega == 0)] <- min_branch_val
node_data_3c$log10Omega[which(node_data_3c$Omega == 0)] <- min_branch_val






## SLC6A14 Zpos_Apos
data<- fromJSON("Gene_trees/sex_phase_sex_phase/SLC6A14/Skylark_2021_Rasolark_2021_SLC6A14_MG94.json", simplifyVector = FALSE)
branch_data <- data$`branch attributes`$`0`
Branch_length <- sapply(branch_data, function(x) x$nonsynonymous +  x$synonymous)
MLE <- sapply(branch_data, function(x) x$`Confidence Intervals`$MLE)

# Map tree data to tree
tree_SLC6A14 <- tree
tree_SLC6A14$branch.length <- Branch_length[match(tree_SLC6A14$label, names(Branch_length))]
tree_SLC6A14$branch.length[is.na(tree_SLC6A14$branch.length)] <- 0
tree_SLC6A14$Omega <- MLE[match(tree_SLC6A14$label, names(MLE))]
tree_SLC6A14$Omega[which(tree_SLC6A14$label == "")] <- tree_SLC6A14$Omega[which(tree_SLC6A14$label == "BranchO")]
tree_SLC6A14$log10Omega <- log10(tree_SLC6A14$Omega)
node_data <- tree_SLC6A14[,c("label", "Omega", "log10Omega")]
node_data_SLC6A14 <- node_data[which(node_data$label == "BranchA" | node_data$label == "BranchZ" | node_data$label == "BranchW"),]
tree_SLC6A14 <- as.phylo(tree_SLC6A14)


## LOC100223640 (FNDC3A) Zpos_W_pos 
data<- fromJSON("Gene_trees/sex_phase_sex_phase/LOC100223640/Skylark_2021_Rasolark_2021_LOC100223640_MG94.json", simplifyVector = FALSE)
branch_data <- data$`branch attributes`$`0`
Branch_length <- sapply(branch_data, function(x) x$nonsynonymous +  x$synonymous)
MLE <- sapply(branch_data, function(x) x$`Confidence Intervals`$MLE)

# Map tree data to tree
tree_LOC100223640 <- tree
tree_LOC100223640$branch.length <- Branch_length[match(tree_LOC100223640$label, names(Branch_length))]
tree_LOC100223640$branch.length[is.na(tree_LOC100223640$branch.length)] <- 0
tree_LOC100223640$Omega <- MLE[match(tree_LOC100223640$label, names(MLE))]
tree_LOC100223640$Omega[which(tree_LOC100223640$label == "")] <- tree_LOC100223640$Omega[which(tree_LOC100223640$label == "BranchO")]
tree_LOC100223640$log10Omega <- log10(tree_LOC100223640$Omega)
node_data <- tree_LOC100223640[,c("label", "Omega", "log10Omega")]
node_data_LOC100223640 <- node_data[which(node_data$label == "BranchA" | node_data$label == "BranchZ" | node_data$label == "BranchW"),]
tree_LOC100223640 <- as.phylo(tree_LOC100223640)


## SPDYA Wpos_Zrel
data<- fromJSON("Gene_trees/sex_phase_sex_phase/SPDYA/Skylark_2021_Rasolark_2021_SPDYA_MG94.json", simplifyVector = FALSE)
branch_data <- data$`branch attributes`$`0`
Branch_length <- sapply(branch_data, function(x) x$nonsynonymous +  x$synonymous)
MLE <- sapply(branch_data, function(x) x$`Confidence Intervals`$MLE)

# Map tree data to tree
tree_SPDYA <- tree
tree_SPDYA$branch.length <- Branch_length[match(tree_SPDYA$label, names(Branch_length))]
tree_SPDYA$branch.length[is.na(tree_SPDYA$branch.length)] <- 0
tree_SPDYA$Omega <- MLE[match(tree_SPDYA$label, names(MLE))]
tree_SPDYA$Omega[which(tree_SPDYA$label == "")] <- tree_SPDYA$Omega[which(tree_SPDYA$label == "BranchO")]
tree_SPDYA$log10Omega <- log10(tree_SPDYA$Omega)
node_data <- tree_SPDYA[,c("label", "Omega", "log10Omega")]
node_data_SPDYA <- node_data[which(node_data$label == "BranchA" | node_data$label == "BranchZ" | node_data$label == "BranchW"),]
tree_SPDYA <- as.phylo(tree_SPDYA)


## CKAP5 Wpos_Apos
data<- fromJSON("Gene_trees/sex_phase_sex_phase/CKAP5/Skylark_2021_Rasolark_2021_CKAP5_MG94.json", simplifyVector = FALSE)
branch_data <- data$`branch attributes`$`0`
Branch_length <- sapply(branch_data, function(x) x$nonsynonymous +  x$synonymous)
MLE <- sapply(branch_data, function(x) x$`Confidence Intervals`$MLE)

# Map tree data to tree
tree_CKAP5 <- tree
tree_CKAP5$branch.length <- Branch_length[match(tree_CKAP5$label, names(Branch_length))]
tree_CKAP5$branch.length[is.na(tree_CKAP5$branch.length)] <- 0
tree_CKAP5$Omega <- MLE[match(tree_CKAP5$label, names(MLE))]
tree_CKAP5$Omega[which(tree_CKAP5$label == "")] <- tree_CKAP5$Omega[which(tree_CKAP5$label == "BranchO")]
tree_CKAP5$log10Omega <- log10(tree_CKAP5$Omega)
node_data <- tree_CKAP5[,c("label", "Omega", "log10Omega")]
node_data_CKAP5 <- node_data[which(node_data$label == "BranchA" | node_data$label == "BranchZ" | node_data$label == "BranchW"),]
tree_CKAP5 <- as.phylo(tree_CKAP5)


## C5H11orf16 Apos_Wrel
data<- fromJSON("Gene_trees/sex_phase_sex_phase/C5H11orf16/Skylark_2021_Rasolark_2021_C5H11orf16_MG94.json", simplifyVector = FALSE)
branch_data <- data$`branch attributes`$`0`
Branch_length <- sapply(branch_data, function(x) x$nonsynonymous +  x$synonymous)
MLE <- sapply(branch_data, function(x) x$`Confidence Intervals`$MLE)

# Map tree data to tree
tree_C5H11orf16 <- tree
tree_C5H11orf16$branch.length <- Branch_length[match(tree_C5H11orf16$label, names(Branch_length))]
tree_C5H11orf16$branch.length[is.na(tree_C5H11orf16$branch.length)] <- 0
tree_C5H11orf16$Omega <- MLE[match(tree_C5H11orf16$label, names(MLE))]
tree_C5H11orf16$Omega[which(tree_C5H11orf16$label == "")] <- tree_C5H11orf16$Omega[which(tree_C5H11orf16$label == "BranchO")]
tree_C5H11orf16$log10Omega <- log10(tree_C5H11orf16$Omega)
node_data <- tree_C5H11orf16[,c("label", "Omega", "log10Omega")]
node_data_C5H11orf16 <- node_data[which(node_data$label == "BranchA" | node_data$label == "BranchZ" | node_data$label == "BranchW"),]
tree_C5H11orf16 <- as.phylo(tree_C5H11orf16)


## URB2 Zpos_Wrel
data<- fromJSON("Gene_trees/sex_phase_sex_phase/URB2/Skylark_2021_Rasolark_2021_URB2_MG94.json", simplifyVector = FALSE)
branch_data <- data$`branch attributes`$`0`
Branch_length <- sapply(branch_data, function(x) x$nonsynonymous +  x$synonymous)
MLE <- sapply(branch_data, function(x) x$`Confidence Intervals`$MLE)

# Map tree data to tree
tree_URB2 <- tree
tree_URB2$branch.length <- Branch_length[match(tree_URB2$label, names(Branch_length))]
tree_URB2$branch.length[is.na(tree_URB2$branch.length)] <- 0
tree_URB2$Omega <- MLE[match(tree_URB2$label, names(MLE))]
tree_URB2$Omega[which(tree_URB2$label == "")] <- tree_URB2$Omega[which(tree_URB2$label == "BranchO")]
tree_URB2$log10Omega <- log10(tree_URB2$Omega)
node_data <- tree_URB2[,c("label", "Omega", "log10Omega")]
node_data_URB2 <- node_data[which(node_data$label == "BranchA" | node_data$label == "BranchZ" | node_data$label == "BranchW"),]
tree_URB2 <- as.phylo(tree_URB2)


### Plots
plot_tree_all <- ggtree(tree_all) %<+% as.data.frame(node_data_all) +
  geom_tree(aes(color = log10Omega), size = 3) +
  scale_color_viridis(limits=c(-1.5,1.5), breaks=c(-1, 0, 1), name=expression("log"[10]*"(dN/dS)")) +
  xlim(0,0.04) +
  ggtitle("All neo-strata") + 
  theme(legend.position= "none",
        plot.title = element_text(size=30, hjust=0.1, vjust = -10, face = "bold"))

plot_tree_4A <- ggtree(tree_4A) %<+% as.data.frame(node_data_4A) +
  geom_tree(aes(color = log10Omega), size = 3) +
  scale_color_viridis(limits=c(-1.5,1.5), breaks=c(-1, 0, 1)) +
  ggtitle("4A") +
  xlim(0,0.04) +
  theme(legend.position= "none",
        plot.title = element_text(size=30, hjust=0.1, vjust = -10, face = "bold"))

plot_tree_3a <- ggtree(tree_3a) %<+% as.data.frame(node_data_3a) +
  geom_tree(aes(color = log10Omega), size = 3) +
  scale_color_viridis(limits=c(-1.5,1.5), breaks=c(-1, 0, 1)) +
  xlim(0,0.04) +
  ggtitle("3-a") +
  theme(legend.position= "none",
        plot.title = element_text(size=30, hjust=0.1, vjust = -10, face = "bold"))

plot_tree_3b <- ggtree(tree_3b) %<+% as.data.frame(node_data_3b) +
  geom_tree(aes(color = log10Omega), size = 3) +
  scale_color_viridis(limits=c(-1.5,1.5), breaks=c(-1, 0, 1)) +
  xlim(0,0.04) +
  ggtitle("3-b") +
  theme(legend.position= "none",
        plot.title = element_text(size=30, hjust=0.1, vjust = -10, face = "bold"))

plot_tree_5 <- ggtree(tree_5) %<+% as.data.frame(node_data_5) +
  geom_tree(aes(color = log10Omega), size = 3) +
  scale_color_viridis(limits=c(-1.5,1.5), breaks=c(-1, 0, 1)) +
  xlim(0,0.04) +
  ggtitle("5") +
  theme(legend.position= "none",
        plot.title = element_text(size=30, hjust=0.1, vjust = -10, face = "bold"))

plot_tree_3c <- ggtree(tree_3c) %<+% as.data.frame(node_data_3c) +
  geom_tree(aes(color = log10Omega), size = 3) +
  scale_color_viridis(limits=c(-1.5,1.5), breaks=c(-1, 0, 1)) +
  xlim(0,0.04) +
  ggtitle("3-c") +
  theme(legend.position= "none",
        plot.title = element_text(size=30, hjust=0.1, vjust = -10, face = "bold"))

annotation_plot <- ggtree(tree_3c) %<+% as.data.frame(node_data_3c) +
  geom_tree(aes(color = log10Omega), alpha = 0) +
  geom_tree(color="white", size=3) +
  annotate(geom="text", x=c(0, 0, 0, 0, 0, 0, 0) , y=c(36.5, 27.1, 17.5, 8, 3, 2, 1), label=c("SW", "RW", "SZ", "RZ", "GTA", "ZFA", "CFA"), color="Black", size=7, hjust=0) +
  theme_void() +
  theme(legend.position = "none")

legend_plot <- ggtree(tree_3c) %<+% as.data.frame(node_data_3c) +
  geom_tree(aes(color = log10Omega), alpha = 0) +
  geom_tree(color="white", size=3) +
  scale_color_viridis(limits=c(-1.5,1.5), breaks=c(-1, 0, 1), name=expression("log"[10]*"(dN/dS)")) +
  theme_void() +
  theme(legend.position = c(0.5, 0.5),
        legend.key.height = unit(2, 'cm'),
        legend.key.width = unit(1.5, 'cm'),
        legend.title = element_text(size=30),
        legend.text = element_text(size=30, hjust=1.3),
        legend.key.size = unit(3,"line"))


scale_plot <- ggtree(tree_3c) %<+% as.data.frame(node_data_3c) +
  geom_tree(aes(color = log10Omega), alpha = 0) +
  geom_tree(color="white", size=3) +
  geom_treescale(x = 0, y = -2, width = 0.02, linesize=2, fontsize=8)  +
  xlim(0,0.04) +
  theme_void() +
  theme(legend.position = "none")

row1a <- (plot_tree_all | annotation_plot | plot_tree_4A | annotation_plot | plot_spacer()) + plot_layout(widths=c(2,0.2,2,0.2,1)) & theme(plot.margin = margin(10,10,10,10))
row2a <- (plot_tree_3a | annotation_plot | plot_tree_3b | annotation_plot | legend_plot) + plot_layout(widths=c(2,0.2,2,0.2,1)) & theme(plot.margin = margin(10,10,10,10))
row3a <- (plot_tree_5 | annotation_plot | plot_tree_3c | annotation_plot | plot_spacer()) + plot_layout(widths=c(2,0.2,2,0.2,1)) & theme(plot.margin = margin(10,10,10,10))
row4a <- (scale_plot | plot_spacer() | plot_spacer() | plot_spacer() | plot_spacer()) + plot_layout(widths=c(2,0.2,2,0.2,1)) & theme(plot.margin = margin(10,10,10,10))

trees_plot <- row1a/row2a/row3a/row4a

jpeg("Figures/Supplementary/trees_plot.jpg", width=6000, height=12000, res=300)
trees_plot
dev.off()


# Main plot
plot_tree_all2 <- ggtree(tree_all) %<+% as.data.frame(node_data_all2) +
  geom_tree(aes(color = log10Omega), size = 3) +
  scale_color_viridis(limits=c(-1.5,1.5), breaks=c(-1, 0, 1), name=expression("log"[10]*"(dN/dS)")) +
  xlim(0,0.05) +
  ggtitle("All neo-strata") + 
  annotate(geom="text", x=c(0.0125, 0.022, 0.0215) , y=c(23.5, 33, 11.5), label=c("BranchA", "BranchW", "BranchZ"), color="Black", size=5) +
  theme(legend.position= "none",
        plot.title = element_text(size=25, hjust=0.1, vjust = -10, face = "bold"))

plot_tree_LOC100223640 <- ggtree(tree_LOC100223640) %<+% as.data.frame(node_data_LOC100223640) +
  geom_tree(aes(color = log10Omega), size = 3) +
  scale_color_viridis(limits=c(-1.5,1.5), breaks=c(-1, 0, 1), name=expression("log"[10]*"(dN/dS)")) +
  xlim(0,0.05) +
  ggtitle("4A: LOC100223640/FNDC3A") +
  geom_segment(aes(x=0.0235, y=24, xend=0.0235, yend=31), arrow = arrow(), linewidth=1) +
  geom_segment(aes(x=0.0235, y=20.5, xend=0.0235, yend=13.5), arrow = arrow(), linewidth=1) +
  annotate(geom="text", x=c(0.009, 0.009, 0.025, 0.025, 0.022, 0.022, 0.022, 0.0235) , y=c(25, 23.5, 34.5, 33, 11.5, 10, 8.5, 22.25), label=c("BranchA", "[0.254, 1.00]", "BranchW *", "[0.336, 1.00]", "BranchZ ***", "[0.00, 0.998]", "[5.00, 0.00243]", "K=0.457"), color="Black", size=5) +
  theme(legend.position= "none",
        plot.title = element_text(size=25, hjust=0.1, vjust = -10, face = "bold"))

plot_tree_SLC6A14 <- ggtree(tree_SLC6A14) %<+% as.data.frame(node_data_SLC6A14) +
  geom_tree(aes(color = log10Omega), size = 3) +
  scale_color_viridis(limits=c(-1.5,1.5), breaks=c(-1, 0, 1), name=expression("log"[10]*"(dN/dS)")) +
  xlim(0,0.10) +
  ggtitle("4A: SLC6A14") +
  geom_segment(aes(x=0.028, y=24, xend=0.028, yend=31), arrow = arrow(), linewidth=1) +
  geom_segment(aes(x=0.028, y=20.5, xend=0.028, yend=13.5), arrow = arrow(), linewidth=1) +
  annotate(geom="text", x=c(0.012, 0.012, 0.012, 0.028, 0.028, 0.05, 0.05, 0.05, 0.028) , y=c(26.5, 25, 23.5, 34.5, 33, 11.5, 10, 8.5, 22.25), label=c("BranchA ***", "[NaN, 0.938]", "[1.83, 0.0616]", "BranchW", "[0.226, 1.00]", "BranchZ **", "[0.00, 0.978]", "[1.75, 0.0216]", "K=-0.248"), color="Black", size=5) +
  theme(legend.position= "none",
        plot.title = element_text(size=25, hjust=0.1, vjust = -10, face = "bold"))

plot_tree_SPDYA <- ggtree(tree_SPDYA) %<+% as.data.frame(node_data_SPDYA) +
  geom_tree(aes(color = log10Omega), size = 3) +
  scale_color_viridis(limits=c(-1.5,1.5), breaks=c(-1, 0, 1), name=expression("log"[10]*"(dN/dS)")) +
  xlim(0,0.05) +
  ggtitle("3-a: SPDYA") +
  geom_segment(aes(x=0.027, y=24, xend=0.027, yend=31), arrow = arrow(), linewidth=1) +
  geom_segment(aes(x=0.027, y=20.5, xend=0.027, yend=13.5), arrow = arrow(), linewidth=1) +
  annotate(geom="text", x=c(0.01, 0.01, 0.027, 0.027, 0.027, 0.027, 0.027, 0.027) , y=c(25, 23.5, 36, 34.5, 33, 11.5, 10, 22.25), label=c("BranchA", "[-1.01, 1.00]", "BranchW **", "[0.00, 0.867]", "[1.92, 0.133]", "BranchZ", "[-0.414, 1.00]", "K=0.674 ** "), color="Black", size=5) +
  theme(legend.position= "none",
        plot.title = element_text(size=25, hjust=0.1, vjust = -10, face = "bold"))

plot_tree_CKAP5 <- ggtree(tree_CKAP5) %<+% as.data.frame(node_data_CKAP5) +
  geom_tree(aes(color = log10Omega), size = 3) +
  scale_color_viridis(limits=c(-1.5,1.5), breaks=c(-1, 0, 1), name=expression("log"[10]*"(dN/dS)")) +
  xlim(0,0.05) +
  ggtitle("5: CKAP5") +
  geom_segment(aes(x=0.019, y=23.25, xend=0.019, yend=31), arrow = arrow(), linewidth=1) +
  geom_segment(aes(x=0.019, y=19.75, xend=0.019, yend=13.5), arrow = arrow(), linewidth=1) +
  annotate(geom="text", x=c(0.0115, 0.0115, 0.0115, 0.019, 0.019, 0.019, 0.019, 0.019, 0.019) , y=c(26.5, 25, 23.5, 36, 34.5, 33, 11.5, 10, 21.5), label=c("BranchA **", "[-1.16, 0.999]", "[2.00, 0.00149]", "BranchW *", "[-0.789, 0.986]", "[1.83, 0.0138]", "BranchZ", "[-1.27, 1.00]", "K=-5.21"), color="Black", size=5) +
  theme(legend.position= "none",
        plot.title = element_text(size=25, hjust=0.1, vjust = -10, face = "bold"))


plot_tree_C5H11orf16 <- ggtree(tree_C5H11orf16) %<+% as.data.frame(node_data_C5H11orf16) +
  geom_tree(aes(color = log10Omega), size = 3) +
  xlim(0,0.05) +
  ggtitle("5: C5H11orf16") +
  scale_color_viridis(limits=c(-1.5,1.5), breaks=c(-1, 0, 1), name=expression("log"[10]*"(dN/dS)")) +
  geom_segment(aes(x=0.046, y=24, xend=0.046, yend=31), arrow = arrow(), linewidth=1) +
  geom_segment(aes(x=0.046, y=20.5, xend=0.046, yend=13.5), arrow = arrow(), linewidth=1) +
  annotate(geom="text", x=c(0.028, 0.028, 0.028, 0.047, 0.047, 0.045, 0.045, 0.046) , y=c(26.5, 25, 23.5, 34.5, 33, 11.5, 10, 22.5), label=c("BranchA *", "[-0.134, 0.993]", "[2.89, 0.00729]", "BranchW", "[-0.219, 1.00]", "BranchZ", "[10.0, 1.00]", "K=-1.08 *"), color="Black", size=5) +
  theme(legend.position= "none",
        plot.title = element_text(size=25, hjust=0.1, vjust = -10, face = "bold"))


plot_tree_URB2 <- ggtree(tree_URB2) %<+% as.data.frame(node_data_URB2) +
  geom_tree(aes(color = log10Omega), size = 3) +
  xlim(0,0.05) +
  ggtitle("3-c: URB2") +
  scale_color_viridis(limits=c(-1.5,1.5), breaks=c(-1, 0, 1), name=expression("log"[10]*"(dN/dS)")) +
  geom_segment(aes(x=0.033, y=24, xend=0.033, yend=31), arrow = arrow(), linewidth=1) +
  geom_segment(aes(x=0.033, y=20.5, xend=0.033, yend=13.5), arrow = arrow(), linewidth=1) +
  annotate(geom="text", x=c(0.0205, 0.0205, 0.0205, 0.034, 0.034, 0.0325, 0.0325, 0.0325, 0.033) , y=c(26.5, 25, 23.5, 34.5, 33, 11.5, 10, 8.5, 22.5), label=c("BranchA", "[NaN, 0.869]", "[0.471, 0.131]", "BranchW", "[-0.817, 1.00]", "BranchZ *", "[-1.10, 0.988]", "[2.34, 0.0124]", "K=-2.31 *"), color="Black", size=5) +
  theme(legend.position= "none",
        plot.title = element_text(size=25, hjust=0.1, vjust = -10, face = "bold"))

scale_plot <- ggtree(tree_3c) %<+% as.data.frame(node_data_3c) +
  geom_tree(aes(color = log10Omega), alpha = 0) +
  geom_tree(color="white", size=3) +
  geom_treescale(x = 0, y = -2, width = 0.02, linesize=2, fontsize=8)  +
  xlim(0,0.05) +
  theme_void() +
  theme(legend.position = "none")

 

row1b <- (plot_tree_all2 | annotation_plot | plot_tree_LOC100223640 | annotation_plot) + plot_layout(widths=c(1.8,0.35,1.8,0.35)) & theme(plot.margin = margin(10,10,10,10))
row2b <- (plot_tree_SLC6A14 | annotation_plot | plot_spacer() | legend_plot) + plot_layout(widths=c(3.6,0.2,0,0.5)) & theme(plot.margin = margin(10,10,10,10))
row3b <- (plot_tree_SPDYA | annotation_plot | plot_tree_CKAP5 | annotation_plot) + plot_layout(widths=c(1.8,0.35,1.8,0.35)) & theme(plot.margin = margin(10,10,10,10))
row4b <- (scale_plot | plot_spacer() | plot_spacer() | plot_spacer()) + plot_layout(widths=c(1.8,0.35,1.8,0.35)) & theme(plot.margin = margin(10,10,10,10))

fig6_plot <- row1b/row2b/row3b/row4b


jpeg("Figures/Figure6.jpg", width=6000, height=12000, res=300)
fig6_plot
dev.off()
