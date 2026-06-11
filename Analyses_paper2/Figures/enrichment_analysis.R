## Export variables and load libraries
rm(list=ls())
library(tidyverse)
library(patchwork)
library(biomaRt)
library(clusterProfiler)
library(enrichplot)
sessionInfo()

options(scipen=999)
setwd("C:/Users/Simon JE/OneDrive - Lund University/Dokument/Simon/PhD/Projects/")

### Import evolutionary data
data_evo <- read.delim("Skylark_2021/Results/Genes/Skylark_2021_Rasolark_2021_organised_data2.tsv", sep="\t", head=T)
GERP <- read.delim("Rasolark_2021/Results/GERP/GERP_scores_genes.tsv", sep="\t", head=T)

### Filter evolutionary data
data_evo <- data_evo[which(data_evo$Filter3=="OK" & data_evo$Filter4=="OK" & data_evo$Filter5=="OK"),]
data_evo$Strata2 <- data_evo$Strata
data_evo$Strata[which(data_evo$Strata == "PAR3" | data_evo$Strata == "PAR5")] <- "PAR"
data_evo$Strata <- factor(data_evo$Strata, order=T, labels=rev(c("S0", "S1", "S2", "S3", "4A", "3-a", "3-b", "5", "3-c", "PAR", "Autosomal")), levels=rev(c("S0", "S1", "S2", "S3", "4A", "3a", "3b", "5", "3c", "PAR", "autosomal")))
data_evo$Species <- factor(data_evo$Species, order=T, labels=c("Skylark", "Raso lark"), levels=c("Skylark", "Rasolark"))
data <- data_evo

### Import genotype data
data_gt <- read.delim("Rasolark_2021/Results/Genotypes/genotype_info.txt", sep="\t", head=T)
gene_strata <- unique(cbind(data_evo$geneID, as.character(data_evo$Strata2)))
gene_strata <- gene_strata[which(gene_strata[,2] != "Autosomal"), ]
for(i in 1:nrow(gene_strata)) {
  data_gt$Strata[which(data_gt$geneID == gene_strata[i,1])] <- gene_strata[i,2] 
}
data_gt$Heterozygosity_score_LOF[which(is.na(data_gt$Heterozygosity_score_LOF))] <- 0
data_gt$Heterozygosity_score_HIGH[which(is.na(data_gt$Heterozygosity_score_HIGH))] <- 0

### Filter genotype data
data_gt <- data_gt %>% filter(geneID %in% data_evo$geneID)
data_evo <- data_evo %>% filter(geneID %in% data_gt$geneID)
data_gt$Strata2 <- data_gt$Strata
data_gt$Strata[which(data_gt$Strata == "PAR3" | data_gt$Strata == "PAR5")] <- "PAR"
data_gt$Strata <- factor(data_gt$Strata, order=T, labels=rev(c("S0", "S1", "S2", "S3", "4A", "3-a", "3-b", "5", "3-c", "PAR", "Autosomal")), levels=rev(c("S0", "S1", "S2", "S3", "4A", "3a", "3b", "5", "3c", "PAR", "autosomal")))
data_gt$Species <- factor(str_split(data_gt$Project, "_", simplify=T)[,1], labels=c("Raso lark", "Skylark"), levels=c("Rasolark", "Skylark"))
data_gt$Sex <-factor(data_gt$Sex, levels=c("Female", "Male"))

# Add gene specific data
data_gt$pHaplo <- rep(NA, nrow(data_gt))
data_gt$geneLengthDataBase <- rep(NA, nrow(data_gt))
data_gt$Strata_Age_Generations <- rep(NA, nrow(data_gt))
data_gt$GERP_n <- rep(NA, nrow(data_gt))
data_gt$GERP_mean <- rep(NA, nrow(data_gt))
data_gt$GERP_sd <- rep(NA, nrow(data_gt))
data_gt$GERP_max <- rep(NA, nrow(data_gt))

unique_genes <- data_evo[which(data_evo$Species == "Skylark"),]

for(i in 1:nrow(unique_genes)) {
  indices <- which(data_gt$geneID == unique_genes$geneID[i])
  data_gt$Strata_Age_Generations[indices] <- unique_genes$Strata_Age_Generations[i]
  data_gt$geneLengthDataBase[indices] <- unique_genes$geneLengthDataBase[i]
  data_gt$pHaplo[indices] <- unique_genes$pHaplo[i]
  indices2 <- which(GERP$geneID == unique_genes$geneID[i])
  data_gt$GERP_n[indices] <- GERP$n[indices2]
  data_gt$GERP_mean[indices] <- GERP$mean[indices2]
  data_gt$GERP_sd[indices] <- GERP$sd[indices2]
  data_gt$GERP_max[indices] <- GERP$max[indices2]
}
rm(unique_genes)
data_gt <- data_gt[which(data_gt$GERP_n > 100),]
data_gt$logGeneLen <- log10(data_gt$geneLengthDataBase)
data_gt <- data_gt[which(!is.na(data_gt$pHaplo)),]

# Set up genotypes
data_gt$func_hom_score_LOF <- rep(0, nrow(data_gt))
data_gt$func_hom_score_LOF[which(data_gt$Heterozygosity_score_LOF == 0 & data_gt$LOF1 == 0 & data_gt$LOF2 == 0 & data_gt$W_Dropout == 0)] <- 1
data_gt$LOF_hom_score_LOF <- rep(0, nrow(data_gt))
data_gt$LOF_hom_score_LOF[which((data_gt$LOF1 == 1 & (data_gt$LOF2 == 1 & data_gt$W_Dropout == 0)) | (data_gt$LOF1 == 1 & data_gt$LOF2 == 0 & data_gt$W_Dropout == 1))] <- 1
data_gt$Mask_Het_score_LOF <- rep(0, nrow(data_gt))
data_gt$Mask_Het_score_LOF[which((data_gt$LOF1 == 0 & (data_gt$LOF2 == 1 | data_gt$W_Dropout == 1)) | (data_gt$LOF1 == 1 & data_gt$LOF2 == 0 & data_gt$W_Dropout == 0))] <- 1
data_gt$Func_Het_score_LOF <- data_gt$Heterozygosity_score_LOF
data_gt$Tot_Het_score_LOF <- rep(0, nrow(data_gt))
data_gt$Tot_Het_score_LOF[which(data_gt$Func_Het_score_LOF == 1 | data_gt$Mask_Het_score_LOF == 1)] <- 1
data_gt2 <- data_gt[which(data_gt$Region != "autosomal" & data_gt$Sex == "Female"),]
data_gt2 <- data_gt2[which(data_gt2$Species == "Raso lark"),]
data_gt2 <- data_gt2[which(data_gt2$Strata != "S0" & data_gt2$Strata != "S1" & data_gt2$Strata != "S2" & data_gt2$Strata != "S3"),]
data1 <- data[which(data$Species == "Raso lark" & data$Region != "autosomal"),]
data1 <- data1[which(data1$Strata != "S0" & data1$Strata != "S1" & data1$Strata != "S2" & data1$Strata != "S3"),]

###### Set up data for GO and KEGG analysis ########################################################################################################################

# Set up the Ensembl connection for Zebra Finch
#ensembl <- useMart("ensembl", host="https://may2024.archive.ensembl.org") # Version 112
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

######### Find common genotypes ################################
prop <- 1

data1$funchom <- rep(NA, nrow(data1))
funchom <- data_gt2 |> group_by(geneID) |> summarise(n_total = n(), n   = sum(func_hom_score_LOF == 1, na.rm = TRUE), prop = mean(func_hom_score_LOF, na.rm = TRUE))
funchom_genes <- funchom$geneID[which(funchom$prop >= prop)]
data1$funchom[which(data1$geneID %in% funchom_genes)] <- 1

data1$funclof <- rep(NA, nrow(data1))
funclof <- data_gt2 |> group_by(geneID) |> summarise(n_total = n(), n   = sum(LOF_hom_score_LOF == 1, na.rm = TRUE), prop = mean(LOF_hom_score_LOF, na.rm = TRUE))
funclof_genes <- funclof$geneID[which(funclof$prop >= prop)]
data1$funclof[which(data1$geneID %in% funclof_genes)] <- 1

data1$maskhet <- rep(NA, nrow(data1))
maskhet <- data_gt2 |> group_by(geneID) |> summarise(n_total = n(), n   = sum(Mask_Het_score_LOF == 1, na.rm = TRUE), prop = mean(Mask_Het_score_LOF, na.rm = TRUE))
maskhet_genes <- maskhet$geneID[which(maskhet$prop >= prop)]
data1$maskhet[which(data1$geneID %in% maskhet_genes)] <- 1

data1$funchet <- rep(NA, nrow(data1))
funchet <- data_gt2 |> group_by(geneID) |> summarise(n_total = n(), n   = sum(Func_Het_score_LOF == 1, na.rm = TRUE), prop = mean(Func_Het_score_LOF, na.rm = TRUE))
funchet_genes <- funchet$geneID[which(funchet$prop >= prop)]
data1$funchet[which(data1$geneID %in% funchet_genes)] <- 1

data1$Zshelt <- rep(NA, nrow(data1))
Zshelt <- data_gt2 |> group_by(geneID) |> summarise(n_total = n(), n_Z = sum(Sheltering_gametolog_LOF == "Z", na.rm = TRUE), prop = n_Z / n_total,, .groups = "drop")
Zshelt_genes <- Zshelt$geneID[which(Zshelt$prop >= prop)]
data1$Zshelt[which(data1$geneID %in% Zshelt_genes)] <- 1

data1$Wshelt <- rep(NA, nrow(data1))
Wshelt <- data_gt2 |> group_by(geneID) |> summarise(n_total = n(), n_W = sum(Sheltering_gametolog_LOF == "W", na.rm = TRUE), prop = n_W / n_total,, .groups = "drop")
Wshelt_genes <- Wshelt$geneID[which(Wshelt$prop >= prop)]
data1$Wshelt[which(data1$geneID %in% Wshelt_genes)] <- 1


###### Analyse gene groups with GO and KEGG enrichment ########################################################################################################################

# Background all genes
data_background1 <- data1
back_genes1 <- as.character(unique(as.vector(str_split(paste(data_background1$gg_entrezgene_id, collapse=";"), ";", simplify=T))))
back_genes2 <- as.character(unique(as.vector(str_split(paste(data_background1$hs_entrezgene_id, collapse=";"), ";", simplify=T))))


genes <- c("funchom", "funclof", "maskhet", "funchet", "Zshelt", "Wshelt")
GO_list1 <- list()

for(i in 1:length(genes)) {
  
  test_genes1 <-  as.character(unique(as.vector(str_split(paste(data1$gg_entrezgene_id[which(data1[,which(colnames(data1) == genes[i])] != "Other")], collapse=";"), ";", simplify=T))))
  
  GO_list1[[paste("GO_chicken", genes[i], sep="_")]] <- enrichGO(gene = test_genes1,
                                                                     universe = back_genes1, 
                                                                     OrgDb = "org.Gg.eg.db",
                                                                     keyType = "ENTREZID", # Adjust if using gene symbols or other IDs
                                                                     ont = "BP",  # Choose BP, MF, or CC, or use "ALL" for all categories
                                                                     pvalueCutoff = 0.05)
  
  GO_list1[[paste("KEGG_chicken", genes[i], sep="_")]] <- enrichKEGG(gene = test_genes1,
                                                                         organism = "gga",  # For Chicken, use "gga"; for Human, use "hsa"
                                                                         keyType = "ncbi-geneid",   # Use the appropriate keyType (e.g., "ENTREZID")
                                                                         pvalueCutoff = 0.05)  # Set the p-value cutoff
}

for(i in 1:length(genes)) {
  test_genes2 <- as.character(unique(as.vector(str_split(paste(data1$hs_entrezgene_id[which(data1[,which(colnames(data1) == genes[i])] != "Other")], collapse=";"), ";", simplify=T))))
  
  GO_list1[[paste("GO_human", genes[i], sep="_")]] <- enrichGO(gene = test_genes2,
                                                                   universe = back_genes2, 
                                                                   OrgDb = "org.Hs.eg.db",
                                                                   keyType = "ENTREZID", # Adjust if using gene symbols or other IDs
                                                                   ont = "BP",  # Choose BP, MF, or CC, or use "ALL" for all categories
                                                                   pvalueCutoff = 0.05)
  
  GO_list1[[paste("KEGG_human", genes[i], sep="_")]] <- enrichKEGG(gene = test_genes2,
                                                                       organism = "hsa",  # For Chicken, use "gga"; for Human, use "hsa"
                                                                       keyType = "ncbi-geneid",   # Use the appropriate keyType (e.g., "ENTREZID")
                                                                       pvalueCutoff = 0.05)  # Set the p-value cutoff
}

GO_list2 <- GO_list1

for(i in 1:length(genes)) {
  if(length(GO_list2[[paste("GO_chicken",  genes[i], sep="_")]]) > 0) {
    GO_list2[[paste("GO_chicken",  genes[i], sep="_")]] <- simplify(GO_list1[[paste("GO_chicken",  genes[i], sep="_")]])
  }
  if(length(GO_list2[[paste("GO_human",  genes[i], sep="_")]]) > 0) {
    GO_list2[[paste("GO_human",  genes[i], sep="_")]] <- simplify(GO_list1[[paste("GO_human",  genes[i], sep="_")]])
  }
}

for(i in 1:length(genes)) {
  
  print(paste("GO_chicken", genes[i], sep="_"))
  if(is.null(GO_list2[[paste("GO_chicken", genes[i], sep="_")]])) {
    print("NA")
  } else {
    GO <- GO_list2[[paste("GO_chicken", genes[i], sep="_")]]@result
    GO <- GO$Description[which(GO$p.adjust < 0.05)]
    print(paste(unique(GO), collapse="; "))
  }
  
  print(paste("GO_human", genes[i], sep="_"))
  if(is.null(GO_list2[[paste("GO_human", genes[i], sep="_")]])) {
    print("NA")
  } else {
    GO <- GO_list2[[paste("GO_human", genes[i], sep="_")]]@result
    GO <- GO$Description[which(GO$p.adjust < 0.05)]
    print(paste(unique(GO), collapse="; "))
  }
  
  print(paste("KEGG_chicken", genes[i], sep="_"))
  if(is.null(GO_list2[[paste("KEGG_chicken", genes[i], sep="_")]])) {
    print("NA")
  } else {
    KEGG <- GO_list2[[paste("KEGG_chicken", genes[i], sep="_")]]@result
    KEGG <- KEGG$Description[which(KEGG$p.adjust < 0.05)]
    print(paste(unique(KEGG), collapse="; "))
  }
  
  print(paste("KEGG_human", genes[i], sep="_"))
  if(is.null(GO_list2[[paste("KEGG_human", genes[i], sep="_")]])) {
    print("NA")
  } else {
    KEGG <- GO_list2[[paste("KEGG_human", genes[i], sep="_")]]@result
    KEGG <- KEGG$Description[which(KEGG$p.adjust < 0.05)]
    print(paste(unique(KEGG), collapse="; "))
  }
}