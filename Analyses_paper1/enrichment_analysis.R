## Export variables and load libraries
rm(list=ls())
library(tidyverse)
library(patchwork)
library(biomaRt)
library(clusterProfiler)
library(enrichplot)
sessionInfo()

options(scipen=999)
setwd("C:/Users/Simon JE/OneDrive - Lund University/Dokument/Simon/PhD/Projects/Skylark_2021/Results/")
data <- read.delim("Genes/Skylark_2021_Rasolark_2021_organised_data2.tsv", sep="\t", head=T)
data1 <- data[which(data$Filter1=="OK" & data$Filter3=="OK" & data$Filter4=="OK" & data$Filter5=="OK"),]



###### Set up data for GO and KEGG analysis ########################################################################################################################

# Set up the Ensembl connection for Zebra Finch
ensembl <- useMart("ensembl", host="https://may2024.archive.ensembl.org") # Version 112
ensembl <- useMart("ensembl")
ensembl_zf <- useDataset("tguttata_gene_ensembl", mart = ensembl)
ensembl_gg <- useDataset("ggallus_gene_ensembl", mart = ensembl)
ensembl_hs <- useDataset("hsapiens_gene_ensembl", mart = ensembl)


# Create list
gene_info <- as.data.frame(unique(data1$geneID))
colnames(gene_info)[1] <- "entrezgene_accession"

# Find alternative Zebra finch names
zf_biomart <- getBM(attributes = c("entrezgene_accession", "ensembl_gene_id", "entrezgene_id"),
                    filters = "entrezgene_accession",
                    values = gene_info$entrezgene_accession,
                    mart = ensembl_zf)

gene_info <- merge(gene_info, zf_biomart, by="entrezgene_accession")
colnames(gene_info)[c(1,3)] <- c("zf_entrezgene_accession", "zf_entrezgene_id") 

# Find chicken ortologs
chicken_ortologs <- getBM(attributes = c("ensembl_gene_id", "ggallus_homolog_ensembl_gene"),
                          filters = "ensembl_gene_id",
                          values = gene_info$ensembl_gene_id,  # Replace with your gene list
                          mart = ensembl_zf)

# Find human ortologs
human_ortologs <- getBM(attributes = c("ensembl_gene_id", "hsapiens_homolog_ensembl_gene"),
                        filters = "ensembl_gene_id",
                        values = gene_info$ensembl_gene_id,  # Replace with your gene list
                        mart = ensembl_zf)


gene_info <- merge(gene_info, chicken_ortologs, by="ensembl_gene_id") 
gene_info <- merge(gene_info, human_ortologs, by="ensembl_gene_id") 
colnames(gene_info)[c(1,4)] <- c("zf_ensembl_gene_id", "ensembl_gene_id") 

# Find alternative chicken names
gg_biomart <- getBM(attributes = c("ensembl_gene_id", "entrezgene_accession", "entrezgene_id", "name_1006"),
                    filters = "ensembl_gene_id",
                    values = gene_info$ensembl_gene_id,  # Replace with your gene list
                    mart = ensembl_gg)

gene_info_gg <- merge(gene_info, gg_biomart, by="ensembl_gene_id") 
colnames(gene_info_gg)[c(1,6,7,8)] <- c("gg_ensembl_gene_id", "gg_entrezgene_accession", "gg_entrezgene_id", "gg_name_1006") 

# Find alternative human names
colnames(gene_info)[c(4,5)] <- c("gg_ensembl_gene_id", "ensembl_gene_id") 
hs_biomart <- getBM(attributes = c("ensembl_gene_id", "entrezgene_accession", "entrezgene_id", "name_1006"),
                    filters = "ensembl_gene_id",
                    values = gene_info$ensembl_gene_id,  # Replace with your gene list
                    mart = ensembl_hs)

gene_info_hs <- merge(gene_info, hs_biomart, by="ensembl_gene_id") 
colnames(gene_info_hs)[c(1,6,7,8)] <- c("hs_ensembl_gene_id", "hs_entrezgene_accession", "hs_entrezgene_id", "hs_name_1006") 



###### Retrieve GO terms ########################################################################################################################
data1$gg_entrezgene_id <- rep(NA, nrow(data1))
data1$hs_entrezgene_id <- rep(NA, nrow(data1))
data1$GO_gg_name_1006 <- rep(NA, nrow(data1))
data1$GO_hs_name_1006 <- rep(NA, nrow(data1))
data1$GO_all_name_1006 <- rep(NA, nrow(data1))

genes <- unique(data1$geneID)

for(i in 1:length(genes)) {
  index <- which(data1$geneID == genes[i])
  
  gg_data <- gene_info_gg[which(gene_info_gg$zf_entrezgene_accession == genes[i]),]
  gg_id <- paste(Filter(nzchar, unique(gg_data$gg_entrezgene_id)), collapse=";")
  if(gg_id != "") {
    data1$gg_entrezgene_id[index] <- gg_id
  }
  gg_GO <- paste(Filter(nzchar, unique(gg_data$gg_name_1006)), collapse=";")
  if(gg_GO != "") {
    data1$GO_gg_name_1006[index] <- gg_GO
  }
  
  hs_data <- gene_info_hs[which(gene_info_hs$zf_entrezgene_accession == genes[i]),]
  hs_id <- paste(Filter(nzchar, unique(hs_data$hs_entrezgene_id)), collapse=";")
  if(hs_id != "") {
    data1$hs_entrezgene_id[index] <- hs_id
  }
  hs_GO <- paste(Filter(nzchar, unique(hs_data$hs_name_1006)), collapse=";")
  if(hs_GO != "") {
    data1$GO_hs_name_1006[index] <- hs_GO
  }
  
  all_GO <- paste(Filter(nzchar, unique(c(gg_data$gg_name_1006, hs_data$hs_name_1006))), collapse=";")
  if(all_GO != "") {
    data1$GO_all_name_1006[index] <- all_GO
  }
}


###### Retrieve KEGG terms ########################################################################################################################
genes_gg <- unique(as.vector(str_split(paste(data1$gg_entrezgene_id, collapse=";"), ";", simplify=T)))
genes_gg <- genes_gg[-which(genes_gg == "NA")]
KEGG_database_gg <- enrichKEGG(gene = genes_gg, organism = "gga", keyType = "ncbi-geneid", pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod="none", minGSSize=0, maxGSSize=100000) 

genes_hs <- unique(as.vector(str_split(paste(data1$hs_entrezgene_id, collapse=";"), ";", simplify=T)))
genes_hs <- genes_hs[-which(genes_hs == "NA")]
KEGG_database_hs <- enrichKEGG(gene = genes_hs, organism = "hsa", keyType = "ncbi-geneid", pvalueCutoff = 1, qvalueCutoff = 1, "none", minGSSize=0, maxGSSize=100000)


data1$KEGG_gg_category <- rep(NA, nrow(data1))
data1$KEGG_hs_category <- rep(NA, nrow(data1))
data1$KEGG_all_category <- rep(NA, nrow(data1))

data1$KEGG_gg_subcategory <- rep(NA, nrow(data1))
data1$KEGG_hs_subcategory <- rep(NA, nrow(data1))
data1$KEGG_all_subcategory <- rep(NA, nrow(data1))

data1$KEGG_gg_pathway <- rep(NA, nrow(data1))
data1$KEGG_hs_pathway <- rep(NA, nrow(data1))
data1$KEGG_all_pathway <- rep(NA, nrow(data1))

for(i in 1:nrow(gene_info)) {
  index <- which(data1$geneID == genes[i])
  gg_category <- ""
  hs_category <- ""
  gg_subcategory <- ""
  hs_subcategory <- ""
  gg_pathway <- ""
  hs_pathway <- ""
  
  gg_id <- str_split(data1$gg_entrezgene_id[i], ";", simplify=T)
  if(length(gg_id) > 1 || (gg_id != "NA" & !is.na(gg_id))) {
    for(j in 1:length(gg_id)) {
      kegg_gg_results <- KEGG_database_gg@result[grep(paste("/", gg_id[j], "/", sep=""), KEGG_database_gg@result$geneID),]
      gg_category <- c(gg_category, kegg_gg_results$category)
      gg_subcategory <- c(gg_subcategory, kegg_gg_results$subcategory)
      gg_pathway <- c(gg_pathway, kegg_gg_results$Description)
    }
    gg_category <- paste(Filter(nzchar, unique(gg_category)), collapse=";")
    if(gg_category != "") {
      data1$KEGG_gg_category[index] <- gg_category
    }
    gg_subcategory <- paste(Filter(nzchar, unique(gg_subcategory)), collapse=";")
    if(gg_subcategory != "") {
      data1$KEGG_gg_subcategory[index] <- gg_subcategory
    }
    gg_pathway <- paste(Filter(nzchar, unique(gg_pathway)), collapse=";")
    if(gg_pathway != "") {
      data1$KEGG_gg_pathway[index] <- gg_pathway
    }
  }
  
  hs_id <- str_split(data1$hs_entrezgene_id[i], ";", simplify=T)
  if(length(hs_id) > 1 || (hs_id != "NA" & !is.na(hs_id))) {
    for(j in 1:length(hs_id)) {
      kegg_hs_results <- KEGG_database_hs@result[grep(paste("/", hs_id[j], "/", sep=""), KEGG_database_hs@result$geneID),]
      hs_category <- c(hs_category, kegg_hs_results$category)
      hs_subcategory <- c(hs_subcategory, kegg_hs_results$subcategory)
      hs_pathway <- c(hs_pathway, kegg_hs_results$Description)
    }
    hs_category <- paste(Filter(nzchar, unique(hs_category)), collapse=";")
    if(hs_category != "") {
      data1$KEGG_hs_category[index] <- hs_category
    }
    hs_subcategory <- paste(Filter(nzchar, unique(hs_subcategory)), collapse=";")
    if(hs_subcategory != "") {
      data1$KEGG_hs_subcategory[index] <- hs_subcategory
    }
    hs_pathway <- paste(Filter(nzchar, unique(hs_pathway)), collapse=";")
    if(hs_pathway != "") {
      data1$KEGG_hs_pathway[index] <- hs_pathway
    }
  }
  
  all_category <- paste(Filter(nzchar, unique(c(str_split(gg_category, ";", simplify=T), str_split(hs_category, ";", simplify=T)))), collapse=";")
  if(all_category != "") {
    data1$KEGG_all_category[index] <- all_category
  }
  all_subcategory <- paste(Filter(nzchar, unique(c(str_split(gg_subcategory, ";", simplify=T), str_split(hs_subcategory, ";", simplify=T)))), collapse=";")
  if(all_subcategory != "") {
    data1$KEGG_all_subcategory[index] <- all_subcategory
  }
  all_pathway <- paste(Filter(nzchar, unique(c(str_split(gg_pathway, ";", simplify=T), str_split(hs_pathway, ";", simplify=T)))), collapse=";")
  if(all_pathway != "") {
    data1$KEGG_all_pathway[index] <- all_pathway
  }
}

### Create data sets
data_background1 <- data1[which(data1$Species == "Skylark"),]

data1 <- data1[which(data1$Strata != "S0" & data1$Strata != "S1" & data1$Strata != "S2" & data1$Strata != "S3" & data1$Strata != "Ancestral unknown"),]
data1 <- data1[which(data1$Region == "sex_phase"),]
genes <- names(which(table(data1$geneID) == 2))
data1 <- data1[data1$geneID %in% genes,]
data1 <- data1[which(data1$Species == "Skylark"),]

data1$SexLinked_DNDS_DOSZ <- data1$SexLinkedBetaZ/(data1$SexLinkedAlphaZ + data1$SexLinkedBetaZ)
data1$SexLinked_DNDS_DOSW <- data1$SexLinkedBetaW/(data1$SexLinkedAlphaW + data1$SexLinkedBetaW)
data1$SexLinked_DNDS_DOSA <- data1$SexLinkedBetaA/(data1$SexLinkedAlphaA + data1$SexLinkedBetaA)

gene_summary <- data1[which(data1$SexLinked_aBSREL_p.A2 < 0.05 | data1$SexLinked_aBSREL_p.Z2 < 0.05 | data1$SexLinked_aBSREL_p.W2 < 0.05 | data1$SexLinked_RELAX1_p.WvsZ1 < 0.05 | data1$SexLinked_RELAX2_p.ZWvsA1 < 0.05), c("geneID", "Strata", "SexLinked_aBSREL_p.Z2", "SexLinked_aBSREL_p.W2", "SexLinked_aBSREL_p.A2", "SexLinked_RELAX1_p.WvsZ1",  "SexLinked_RELAX2_p.ZWvsA1", "SexLinked_RELAX1_K.WvsZ", "SexLinked_RELAX2_K.ZWvsA", "SexLinkedAlphaZ","SexLinkedAlphaW", "SexLinkedAlphaA", "SexLinkedBetaZ","SexLinkedBetaW", "SexLinkedBetaA", "PNZ", "PNW", "PNA", "SexLinkedNIZ", "SexLinkedNIW", "SexLinked_DNDS_DOSZ", "SexLinked_DNDS_DOSW", "SexLinked_DNDS_DOSA", "SexLinkedDoSZ", "SexLinkedDoSW", "pHaplo", "geneLengthDataBase", "geneLengthTrimmed", "Wdegeneration", "Zdegeneration")]
write.table(gene_summary, "gene_summary.tsv", quote=F, sep = "\t")

# Define selection regimes ########################################################################################################################
data1$A_pos <- rep("Other", nrow(data1))
data1$Z_pos <- rep("Other", nrow(data1))
data1$W_pos <- rep("Other", nrow(data1))
data1$Z_rel <- rep("Other", nrow(data1))
data1$W_rel <- rep("Other", nrow(data1))

data1$A_pos[which(!is.na(data1$SexLinked_aBSREL_p.A2) & data1$SexLinked_aBSREL_p.A2 < 0.05)] <- "Positive A selection"
data1$Z_pos[which(!is.na(data1$SexLinked_aBSREL_p.Z2) & data1$SexLinked_aBSREL_p.Z2 < 0.05)] <- "Positive Z selection"
data1$W_pos[which(!is.na(data1$SexLinked_aBSREL_p.W2) & data1$SexLinked_aBSREL_p.W2 < 0.05)] <- "Positive W selection"
data1$Z_rel[which(!is.na(data1$SexLinked_RELAX1_K.WvsZ) & data1$SexLinked_RELAX1_K.WvsZ > 1 & data1$SexLinked_RELAX1_p.WvsZ1 < 0.05)] <- "Relaxed Z selection"
data1$W_rel[which(!is.na(data1$SexLinked_RELAX1_K.WvsZ) & data1$SexLinked_RELAX1_K.WvsZ < 1 & data1$SexLinked_RELAX1_p.WvsZ1 < 0.05)] <- "Relaxed W selection"

data1$Z_pos_W_pos <- rep("Other", nrow(data1))
data1$Z_pos_A_pos <- rep("Other", nrow(data1))
data1$W_pos_A_pos <- rep("Other", nrow(data1))
data1$Z_pos_W_rel <- rep("Other", nrow(data1))
data1$W_pos_Z_rel <- rep("Other", nrow(data1))
data1$A_pos_Z_rel <- rep("Other", nrow(data1))
data1$A_pos_W_rel <- rep("Other", nrow(data1))

data1$Z_pos_A_pos[which(data1$Z_pos != "Other" & data1$A_pos != "Other")] <- "Positive Z positive A selection"
data1$W_pos_A_pos[which(data1$W_pos != "Other" & data1$A_pos != "Other")] <- "Positive W positive A selection"
data1$Z_pos_W_rel[which(data1$Z_pos != "Other" & data1$W_rel != "Other")] <- "Positive Z relaxed W selection"
data1$W_pos_Z_rel[which(data1$W_pos != "Other" & data1$Z_rel != "Other")] <- "Positive W relaxed Z selection"
data1$A_pos_Z_rel[which(data1$A_pos != "Other" & data1$Z_rel != "Other")] <- "Positive A relaxed Z selection"
data1$A_pos_W_rel[which(data1$A_pos != "Other" & data1$W_rel != "Other")] <- "Positive A relaxed W selection"

###### Analyse selection regimes with GO and KEGG enrichment ########################################################################################################################

# Background
back_genes1 <- as.character(unique(as.vector(str_split(paste(data_background1$gg_entrezgene_id, collapse=";"), ";", simplify=T))))
back_genes2 <- as.character(unique(as.vector(str_split(paste(data_background1$hs_entrezgene_id, collapse=";"), ";", simplify=T))))


selection <- c("A_pos", "Z_pos", "W_pos", "Z_rel", "W_rel",
               "Z_pos_W_pos", "Z_pos_A_pos", "W_pos_A_pos", "Z_pos_W_rel", "W_pos_Z_rel",
               "A_pos_Z_rel", "A_pos_W_rel")
GO_list1 <- list()

for(i in 1:length(selection)) {
  
  test_genes1 <-  as.character(unique(as.vector(str_split(paste(data1$gg_entrezgene_id[which(data1[,which(colnames(data1) == selection[i])] != "Other")], collapse=";"), ";", simplify=T))))
  
  GO_list1[[paste("GO_chicken", selection[i], sep="_")]] <- enrichGO(gene = test_genes1,
                                                                     universe = back_genes1, 
                                                                     OrgDb = "org.Gg.eg.db",
                                                                     keyType = "ENTREZID", # Adjust if using gene symbols or other IDs
                                                                     ont = "BP",  # Choose BP, MF, or CC, or use "ALL" for all categories
                                                                     pvalueCutoff = 0.05)
  
  GO_list1[[paste("KEGG_chicken", selection[i], sep="_")]] <- enrichKEGG(gene = test_genes1,
                                                                         organism = "gga",  # For Chicken, use "gga"; for Human, use "hsa"
                                                                         keyType = "ncbi-geneid",   # Use the appropriate keyType (e.g., "ENTREZID")
                                                                         pvalueCutoff = 0.05)  # Set the p-value cutoff
}

for(i in 1:length(selection)) {
  test_genes2 <- as.character(unique(as.vector(str_split(paste(data1$hs_entrezgene_id[which(data1[,which(colnames(data1) == selection[i])] != "Other")], collapse=";"), ";", simplify=T))))
  
  GO_list1[[paste("GO_human", selection[i], sep="_")]] <- enrichGO(gene = test_genes2,
                                                                   universe = back_genes2, 
                                                                   OrgDb = "org.Hs.eg.db",
                                                                   keyType = "ENTREZID", # Adjust if using gene symbols or other IDs
                                                                   ont = "BP",  # Choose BP, MF, or CC, or use "ALL" for all categories
                                                                   pvalueCutoff = 0.05)
  
  GO_list1[[paste("KEGG_human", selection[i], sep="_")]] <- enrichKEGG(gene = test_genes2,
                                                                       organism = "hsa",  # For Chicken, use "gga"; for Human, use "hsa"
                                                                       keyType = "ncbi-geneid",   # Use the appropriate keyType (e.g., "ENTREZID")
                                                                       pvalueCutoff = 0.05)  # Set the p-value cutoff
}

GO_list2 <- GO_list1

for(i in 1:length(selection)) {
  if(length(GO_list2[[paste("GO_chicken",  selection[i], sep="_")]]) > 0) {
    GO_list2[[paste("GO_chicken",  selection[i], sep="_")]] <- simplify(GO_list1[[paste("GO_chicken",  selection[i], sep="_")]])
  }
  if(length(GO_list2[[paste("GO_human",  selection[i], sep="_")]]) > 0) {
    GO_list2[[paste("GO_human",  selection[i], sep="_")]] <- simplify(GO_list1[[paste("GO_human",  selection[i], sep="_")]])
  }
}

for(i in 1:length(selection)) {
  
  print(paste("GO_chicken", selection[i], sep="_"))
  if(is.null(GO_list2[[paste("GO_chicken", selection[i], sep="_")]])) {
    print("NA")
  } else {
    GO <- GO_list2[[paste("GO_chicken", selection[i], sep="_")]]@result
    GO <- GO$Description[which(GO$p.adjust < 0.05)]
    print(paste(unique(GO), collapse="; "))
  }
  
  print(paste("GO_human", selection[i], sep="_"))
  if(is.null(GO_list2[[paste("GO_human", selection[i], sep="_")]])) {
    print("NA")
  } else {
    GO <- GO_list2[[paste("GO_human", selection[i], sep="_")]]@result
    GO <- GO$Description[which(GO$p.adjust < 0.05)]
    print(paste(unique(GO), collapse="; "))
  }
  
  print(paste("KEGG_chicken", selection[i], sep="_"))
  if(is.null(GO_list2[[paste("KEGG_chicken", selection[i], sep="_")]])) {
    print("NA")
  } else {
    KEGG <- GO_list2[[paste("KEGG_chicken", selection[i], sep="_")]]@result
    KEGG <- KEGG$Description[which(KEGG$p.adjust < 0.05)]
    print(paste(unique(KEGG), collapse="; "))
  }
  
  print(paste("KEGG_human", selection[i], sep="_"))
  if(is.null(GO_list2[[paste("KEGG_human", selection[i], sep="_")]])) {
    print("NA")
  } else {
    KEGG <- GO_list2[[paste("KEGG_human", selection[i], sep="_")]]@result
    KEGG <- KEGG$Description[which(KEGG$p.adjust < 0.05)]
    print(paste(unique(KEGG), collapse="; "))
  }
}

