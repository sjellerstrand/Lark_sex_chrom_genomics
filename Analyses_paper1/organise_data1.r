## Export variables and load libraries
rm(list=ls())
library(tidyverse)

options(scipen=999)

args <- strsplit(commandArgs(trailingOnly=T), "=", fixed=T)
for(i in 1:length(args)) {
  assign(args[[i]][1], args[[i]][2])
}

### Import data
data <- as.matrix(read.delim(GENES, head=F, skipNul = T))

### Seperate data not present in both datasets
data_other <- as.data.frame(data[which(data[,14] != "Skylark_2021_Rasolark_2021"),])
data_other$V15 <- str_split(data_other[,14], "_", simplify=T)[,1]
data <- data[which(data[,14]== "Skylark_2021_Rasolark_2021"),]
data <- as.data.frame(rbind(cbind(data, rep("Skylark", nrow(data))), cbind(data, rep("Rasolark", nrow(data)))))
data <- rbind(data, data_other)

colnames(data) <- c("geneID", "transcriptID", "Region1", "Region2", "Strata1", "Strata2", "geneLengthAlignment", "geneLengthDataBase", "geneLengthSpecies", "geneLengthTrimmed", "ExonsAll", "ExonsCallable", "LiftOverStatus", "Datasets", "Species")
data2 <- as.data.frame(read.delim(paste(WORKDIR1, "/D5_sequence_evolution/metadata/genes_sequence_evolution_info.tsv", sep=""), head=F))
colnames(data2) <- c("geneID", "Branch", "Length", "DN", "DS", "DNDS", "Beta", "Alpha", "Omega", "PN", "PS", "IngroupOmega")
data3 <- as.data.frame(read.delim(paste(WORKDIR1, "/D6_snpEff_substitutions/metadata/features_shared.tsv", sep=""), head=F))
colnames(data3) <- c("DB_Scaffold", "DB_Start", "DB_End", "DB_Strandedness", "Scaffold", "Start", "End", "Strandedness", "Callable1", "Callable2", "Shared1", "Shared2", "Region1", "Region2", "Strata1", "Strata2", "Exon")
data3$Gene <- str_split(data3$Exon, "_", simplify=T)[,1]

# Prepare data
# General information
data$geneLengthAlignment <- as.numeric(data$geneLengthAlignment)
data$geneLengthDataBase <- as.numeric(data$geneLengthDataBase)
data$geneLengthSpecies <- as.numeric(data$geneLengthSpecies)
data$geneLengthTrimmed <- as.numeric(data$geneLengthTrimmed)
data$ExonsAll <- as.integer(data$ExonsAll)
data$ExonsCallable <- as.integer(data$ExonsCallable)
data$geneScaffold <- rep(NA, nrow(data))
data$geneMidPos <- as.numeric(rep(NA, nrow(data)))

# Within lineage polymorphisms
data$PSA <- as.numeric(rep(NA, nrow(data)))
data$PSZ <- as.numeric(rep(NA, nrow(data)))
data$PSW <- as.numeric(rep(NA, nrow(data)))
data$PNA <- as.numeric(rep(NA, nrow(data)))
data$PNZ <- as.numeric(rep(NA, nrow(data)))
data$PNW <- as.numeric(rep(NA, nrow(data)))
data$IngroupOmegaA <- as.numeric(rep(NA, nrow(data)))
data$IngroupOmegaZ <- as.numeric(rep(NA, nrow(data)))
data$IngroupOmegaW <- as.numeric(rep(NA, nrow(data)))

# Substitutions since Pmaj divergence
data$FullDSA <- as.numeric(rep(NA, nrow(data)))
data$FullDSZ <- as.numeric(rep(NA, nrow(data)))
data$FullDSW <- as.numeric(rep(NA, nrow(data)))
data$FullAlphaA <- as.numeric(rep(NA, nrow(data)))
data$FullAlphaZ <- as.numeric(rep(NA, nrow(data)))
data$FullAlphaW <- as.numeric(rep(NA, nrow(data)))
data$FullBetaA <- as.numeric(rep(NA, nrow(data)))
data$FullBetaZ <- as.numeric(rep(NA, nrow(data)))
data$FullBetaW <- as.numeric(rep(NA, nrow(data)))
data$FullOmegaA <- as.numeric(rep(NA, nrow(data)))
data$FullOmegaZ <- as.numeric(rep(NA, nrow(data)))
data$FullOmegaW <- as.numeric(rep(NA, nrow(data)))
data$FullNIA <- as.numeric(rep(NA, nrow(data)))
data$FullNIZ <- as.numeric(rep(NA, nrow(data)))
data$FullNIW <- as.numeric(rep(NA, nrow(data)))
data$FullDoSA <- as.numeric(rep(NA, nrow(data)))
data$FullDoSZ <- as.numeric(rep(NA, nrow(data)))
data$FullDoSW <- as.numeric(rep(NA, nrow(data)))

# Substitutions since recombination cessation
data$SexLinkedAlphaA <- as.numeric(rep(NA, nrow(data)))
data$SexLinkedAlphaZ <- as.numeric(rep(NA, nrow(data)))
data$SexLinkedAlphaW <- as.numeric(rep(NA, nrow(data)))
data$SexLinkedAlphaZW <- as.numeric(rep(NA, nrow(data)))
data$SexLinkedBetaA <- as.numeric(rep(NA, nrow(data)))
data$SexLinkedBetaZ <- as.numeric(rep(NA, nrow(data)))
data$SexLinkedBetaW <- as.numeric(rep(NA, nrow(data)))
data$SexLinkedBetaZW <- as.numeric(rep(NA, nrow(data)))
data$SexLinkedOmegaA <- as.numeric(rep(NA, nrow(data)))
data$SexLinkedOmegaZ <- as.numeric(rep(NA, nrow(data)))
data$SexLinkedOmegaW <- as.numeric(rep(NA, nrow(data)))
data$SexLinkedOmegaZW <- as.numeric(rep(NA, nrow(data)))
data$SexLinkedNIZ <- as.numeric(rep(NA, nrow(data)))
data$SexLinkedNIW <- as.numeric(rep(NA, nrow(data)))
data$SexLinkedDoSZ <- as.numeric(rep(NA, nrow(data)))
data$SexLinkedDoSW <- as.numeric(rep(NA, nrow(data)))

# Substitutions since species divergence
data$SpeciesAlphaA <- as.numeric(rep(NA, nrow(data)))
data$SpeciesAlphaZ <- as.numeric(rep(NA, nrow(data)))
data$SpeciesAlphaW <- as.numeric(rep(NA, nrow(data)))
data$SpeciesBetaA <- as.numeric(rep(NA, nrow(data)))
data$SpeciesBetaZ <- as.numeric(rep(NA, nrow(data)))
data$SpeciesBetaW <- as.numeric(rep(NA, nrow(data)))
data$SpeciesOmegaA <- as.numeric(rep(NA, nrow(data)))
data$SpeciesOmegaZ <- as.numeric(rep(NA, nrow(data)))
data$SpeciesOmegaW <- as.numeric(rep(NA, nrow(data)))
data$SpeciesNIA <- as.numeric(rep(NA, nrow(data)))
data$SpeciesNIZ <- as.numeric(rep(NA, nrow(data)))
data$SpeciesNIW <- as.numeric(rep(NA, nrow(data)))
data$SpeciesDoSA <- as.numeric(rep(NA, nrow(data)))
data$SpeciesDoSZ <- as.numeric(rep(NA, nrow(data)))
data$SpeciesDoSW <- as.numeric(rep(NA, nrow(data)))

# Substitutions between outgroup lineages
data$DSO <- as.numeric(rep(NA, nrow(data)))
data$AlphaO <- as.numeric(rep(NA, nrow(data)))
data$BetaO <- as.numeric(rep(NA, nrow(data)))
data$OmegaO <- as.numeric(rep(NA, nrow(data)))

data_other <- data[which(data[,14] != "Skylark_2021_Rasolark_2021"),]
data <- data[which(data[,14] == "Skylark_2021_Rasolark_2021"),]

### Add branch info for rate of sequence evoluion
for(i in 1:(nrow(data)/2)) {

  j=i+(nrow(data)/2)

  genePosInfo <-  data3[which(data3$Gene == data$geneID[i] & data3$Start != "missing"),]
  
  if(nrow(genePosInfo) > 0) {
    data$geneScaffold[i] <- genePosInfo$Scaffold[1]
    data$geneScaffold[j] <- data$geneScaffold[i]
    geneMax <- max(as.numeric(c(genePosInfo$Start, genePosInfo$End)))
    geneMin <- min(as.numeric(c(genePosInfo$Start, genePosInfo$End)))
    data$geneMidPos[i] <- geneMin + ((geneMax-geneMin) / 2)
    data$geneMidPos[j] <- data$geneMidPos[i]
  }

  if(nrow(data2[which(data2$geneID == data$geneID[i]),]) > 0) {

    geneinfo <- data2[which(data2$geneID == data$geneID[i]),]
    
    if(data$geneID[i]=="PTP4A1") {

      data$IngroupOmegaZ[i] <-  geneinfo$IngroupOmega[which(geneinfo$Branch == "BranchSkyZ")]
      data$IngroupOmegaW[i] <-  geneinfo$IngroupOmega[which(geneinfo$Branch == "BranchSkyW")]
      data$PNZ[i] <-  geneinfo$PN[which(geneinfo$Branch == "BranchSkyZ")]
      data$PNW[i] <-  geneinfo$PN[which(geneinfo$Branch == "BranchSkyW")]
      data$PSZ[i] <-  geneinfo$PS[which(geneinfo$Branch == "BranchSkyZ")]
      data$PSW[i] <-  geneinfo$PS[which(geneinfo$Branch == "BranchSkyW")]
      data$FullAlphaZ[i] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchA")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyZ")]
      data$FullBetaZ[i] <- geneinfo$Beta[which(geneinfo$Branch == "BranchA")] + geneinfo$Beta[which(geneinfo$Branch == "BranchSkyZ")]
      data$FullDSZ[i] <- geneinfo$DS[which(geneinfo$Branch == "BranchA")] + geneinfo$DS[which(geneinfo$Branch == "BranchSkyZ")]
      data$FullOmegaZ[i] <- (geneinfo$Beta[which(geneinfo$Branch == "BranchA")] + geneinfo$Beta[which(geneinfo$Branch == "BranchSkyZ")]) / (geneinfo$Alpha[which(geneinfo$Branch == "BranchA")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyZ")])
      data$FullNIZ[i] <-  data$IngroupOmegaZ[i] / data$FullOmegaZ[i]
      data$SexLinkedOmegaZ[i] <-  geneinfo$Beta[which(geneinfo$Branch == "BranchSkyZ")] / geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyZ")]
      data$SexLinkedNIZ[i] <- data$IngroupOmegaZ[i] / data$SexLinkedOmegaZ[i]
      data$FullAlphaW[i] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchA")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyW")]
      data$FullBetaW[i] <- geneinfo$Beta[which(geneinfo$Branch == "BranchA")] + geneinfo$Beta[which(geneinfo$Branch == "BranchSkyW")]
      data$FullDSW[i] <- geneinfo$DS[which(geneinfo$Branch == "BranchA")] + geneinfo$DS[which(geneinfo$Branch == "BranchSkyW")]
      data$FullOmegaW[i] <- (geneinfo$Beta[which(geneinfo$Branch == "BranchA")] + geneinfo$Beta[which(geneinfo$Branch == "BranchSkyW")]) / (geneinfo$Alpha[which(geneinfo$Branch == "BranchA")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyW")])
      data$FullNIW[i] <-  data$IngroupOmegaW[i] / data$FullOmegaW[i]
      data$DSO[i] <- geneinfo$DS[which(geneinfo$Branch == "BranchO")] + geneinfo$DS[which(geneinfo$Branch == "Node1")]
      data$AlphaO[i] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchO")] + geneinfo$Alpha[which(geneinfo$Branch == "Node1")]
      data$BetaO[i] <- geneinfo$Beta[which(geneinfo$Branch == "BranchO")] + geneinfo$Beta[which(geneinfo$Branch == "Node1")]
      data$OmegaO[i] <- data$BetaO[i] / data$AlphaO[i]
      data$SexLinkedOmegaW[i] <-  geneinfo$Beta[which(geneinfo$Branch == "BranchSkyW")] / geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyW")]
      data$SexLinkedNIW[i] <- data$IngroupOmegaW[i] / data$SexLinkedOmegaW[i]
      data$SexLinkedAlphaA[i] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchA")]
      data$SexLinkedBetaA[i] <- geneinfo$Beta[which(geneinfo$Branch == "BranchA")]
      data$SexLinkedOmegaA[i] <- geneinfo$Omega[which(geneinfo$Branch == "BranchA")]
      data$SpeciesAlphaZ[i] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyZ")]
      data$SpeciesAlphaW[i] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyW")]
      data$SpeciesBetaZ[i] <- geneinfo$Beta[which(geneinfo$Branch == "BranchSkyZ")]
      data$SpeciesBetaW[i] <- geneinfo$Beta[which(geneinfo$Branch == "BranchSkyW")]
      data$SpeciesOmegaZ[i] <- geneinfo$Omega[which(geneinfo$Branch == "BranchSkyZ")]
      data$SpeciesOmegaW[i] <- geneinfo$Omega[which(geneinfo$Branch == "BranchSkyW")]
      data$SpeciesNIZ[i] <- data$IngroupOmegaZ[i] / data$SpeciesOmegaZ[i]
      data$SpeciesNIW[i] <- data$IngroupOmegaW[i] / data$SpeciesOmegaW[i]
      data$FullDoSZ[i] <- data$FullBetaZ[i]/(data$FullBetaZ[i] + data$FullAlphaZ[i]) - data$PNZ[i]/(data$PNZ[i] + data$PSZ[i])
      data$FullDoSW[i] <- data$FullBetaW[i]/(data$FullBetaW[i] + data$FullAlphaW[i]) - data$PNW[i]/(data$PNW[i] + data$PSW[i])
      data$SexLinkedDoSZ[i] <- data$SexLinkedBetaZ[i]/(data$SexLinkedBetaZ[i] + data$SexLinkedAlphaZ[i]) - data$PNZ[i]/(data$PNZ[i] + data$PSZ[i])
      data$SexLinkedDoSW[i] <- data$SexLinkedBetaW[i]/(data$SexLinkedBetaW[i] + data$SexLinkedAlphaW[i]) - data$PNW[i]/(data$PNW[i] + data$PSW[i])
      data$SpeciesDoSZ[i] <- data$SpeciesBetaZ[i]/(data$SpeciesBetaZ[i] + data$SpeciesAlphaZ[i]) - data$PNZ[i]/(data$PNZ[i] + data$PSZ[i])
      data$SpeciesDoSW[i] <- data$SpeciesBetaW[i]/(data$SpeciesBetaW[i] + data$SpeciesAlphaW[i]) - data$PNW[i]/(data$PNW[i] + data$PSW[i])

      data$IngroupOmegaA[j] <-  geneinfo$IngroupOmega[which(geneinfo$Branch == "BranchRasoA")]
      data$PNA[j] <-  geneinfo$PN[which(geneinfo$Branch == "BranchRasoA")]
      data$PSA[j] <-  geneinfo$PS[which(geneinfo$Branch == "BranchRasoA")]
      data$FullAlphaA[j] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchA")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoA")]
      data$FullBetaA[j] <- geneinfo$Beta[which(geneinfo$Branch == "BranchA")] + geneinfo$Beta[which(geneinfo$Branch == "BranchRasoA")]
      data$FullDSA[j] <- geneinfo$DS[which(geneinfo$Branch == "BranchA")] + geneinfo$DS[which(geneinfo$Branch == "BranchRasoA")]
      data$FullOmegaA[j] <- (geneinfo$Beta[which(geneinfo$Branch == "BranchA")] + geneinfo$Beta[which(geneinfo$Branch == "BranchRasoA")]) / (geneinfo$Alpha[which(geneinfo$Branch == "BranchA")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoA")])
      data$FullNIA[j] <-  data$IngroupOmegaA[j] / data$FullOmegaA[j]
      data$DSO[j] <- geneinfo$DS[which(geneinfo$Branch == "BranchO")] + geneinfo$DS[which(geneinfo$Branch == "Node1")]
      data$AlphaO[j] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchO")] + geneinfo$Alpha[which(geneinfo$Branch == "Node1")]
      data$BetaO[j] <- geneinfo$Beta[which(geneinfo$Branch == "BranchO")] + geneinfo$Beta[which(geneinfo$Branch == "Node1")]
      data$OmegaO[j] <- data$BetaO[j] / data$AlphaO[j]
      data$SpeciesAlphaA[j] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoA")]
      data$SpeciesBetaA[j] <- geneinfo$Beta[which(geneinfo$Branch == "BranchRasoA")]
      data$SpeciesOmegaA[j] <- geneinfo$Omega[which(geneinfo$Branch == "BranchRasoA")]
      data$SpeciesNIA[j] <- data$IngroupOmegaA[j] / data$SpeciesOmegaA[j]
      data$FullDoSA[j] <- data$FullBetaA[j]/(data$FullBetaA[j] + data$FullAlphaA[j]) - data$PNA[j]/(data$PNA[j] + data$PSA[j])
      data$SpeciesDoSA[j] <- data$SpeciesBetaA[j]/(data$SpeciesBetaA[j] + data$SpeciesAlphaA[j]) - data$PNA[j]/(data$PNA[j] + data$PSA[j])

    } else {

      if(data$Region1[i] == "autosomal") {
        data$IngroupOmegaA[i] <-  geneinfo$IngroupOmega[which(geneinfo$Branch == "BranchSkyA")]
        data$PNA[i] <-  geneinfo$PN[which(geneinfo$Branch == "BranchSkyA")]
        data$PSA[i] <-  geneinfo$PS[which(geneinfo$Branch == "BranchSkyA")]
        data$FullAlphaA[i] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchA")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyA")]
        data$FullBetaA[i] <- geneinfo$Beta[which(geneinfo$Branch == "BranchA")] + geneinfo$Beta[which(geneinfo$Branch == "BranchSkyA")]
        data$FullDSA[i] <- geneinfo$DS[which(geneinfo$Branch == "BranchA")] + geneinfo$DS[which(geneinfo$Branch == "BranchSkyA")]
        data$FullOmegaA[i] <- (geneinfo$Beta[which(geneinfo$Branch == "BranchA")] + geneinfo$Beta[which(geneinfo$Branch == "BranchSkyA")]) / (geneinfo$Alpha[which(geneinfo$Branch == "BranchA")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyA")])
        data$FullNIA[i] <-  data$IngroupOmegaA[i] / data$FullOmegaA[i]
        data$DSO[i] <- geneinfo$DS[which(geneinfo$Branch == "BranchO")] + geneinfo$DS[which(geneinfo$Branch == "Node1")]
        data$AlphaO[i] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchO")] + geneinfo$Alpha[which(geneinfo$Branch == "Node1")]
        data$BetaO[i] <- geneinfo$Beta[which(geneinfo$Branch == "BranchO")] + geneinfo$Beta[which(geneinfo$Branch == "Node1")]
        data$OmegaO[i] <- data$BetaO[i] / data$AlphaO[i]
        data$SpeciesAlphaA[i] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyA")]
        data$SpeciesBetaA[i] <- geneinfo$Beta[which(geneinfo$Branch == "BranchSkyA")]
        data$SpeciesOmegaA[i] <- geneinfo$Omega[which(geneinfo$Branch == "BranchSkyA")]
        data$SpeciesNIA[i] <- data$IngroupOmegaA[i] / data$SpeciesOmegaA[i]
        data$FullDoSA[i] <- data$FullBetaA[i]/(data$FullBetaA[i] + data$FullAlphaA[i]) - data$PNA[i]/(data$PNA[i] + data$PSA[i])
        data$SpeciesDoSA[i] <- data$SpeciesBetaA[i]/(data$SpeciesBetaA[i] + data$SpeciesAlphaA[i]) - data$PNA[i]/(data$PNA[i] + data$PSA[i])

      } else if(data$Region1[i] == "sex_dropout") {
        data$IngroupOmegaZ[i] <-  geneinfo$IngroupOmega[which(geneinfo$Branch == "BranchSkyZ")]
        data$PNZ[i] <-  geneinfo$PN[which(geneinfo$Branch == "BranchSkyZ")]
        data$PSZ[i] <-  geneinfo$PS[which(geneinfo$Branch == "BranchSkyZ")]
        data$SpeciesAlphaZ[i] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyZ")]
        data$SpeciesBetaZ[i] <- geneinfo$Beta[which(geneinfo$Branch == "BranchSkyZ")]
        data$SpeciesOmegaZ[i] <- geneinfo$Omega[which(geneinfo$Branch == "BranchSkyZ")]
        data$SpeciesNIZ[i] <- data$IngroupOmegaZ[i] / data$SpeciesOmegaZ[i]
        data$SpeciesDoSZ[i] <- data$SpeciesBetaZ[i]/(data$SpeciesBetaZ[i] + data$SpeciesAlphaZ[i]) - data$PNZ[i]/(data$PNZ[i] + data$PSZ[i])

        if(data$Strata1[i] != "Z") {
        data$DSO[i] <- geneinfo$DS[which(geneinfo$Branch == "BranchO")] + geneinfo$DS[which(geneinfo$Branch == "Node1")]
        data$AlphaO[i] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchO")] + geneinfo$Alpha[which(geneinfo$Branch == "Node1")]
        data$BetaO[i] <- geneinfo$Beta[which(geneinfo$Branch == "BranchO")] + geneinfo$Beta[which(geneinfo$Branch == "Node1")]
        data$OmegaO[i] <- data$BetaO[i] / data$AlphaO[i]

          if(data$Region2[j] == "sex_phase" || data$Region2[j] == "partial_dropout") {
            data$FullAlphaZ[i] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchA")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyZ")]
            data$FullBetaZ[i] <- geneinfo$Beta[which(geneinfo$Branch == "BranchA")] + geneinfo$Beta[which(geneinfo$Branch == "BranchZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchSkyZ")]
            data$FullDSZ[i] <- geneinfo$DS[which(geneinfo$Branch == "BranchA")] + geneinfo$DS[which(geneinfo$Branch == "BranchZ")] + geneinfo$DS[which(geneinfo$Branch == "BranchSkyZ")]
            data$FullOmegaZ[i] <- (geneinfo$Beta[which(geneinfo$Branch == "BranchA")] + geneinfo$Beta[which(geneinfo$Branch == "BranchZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchSkyZ")]) / (geneinfo$Alpha[which(geneinfo$Branch == "BranchA")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyZ")])
            data$FullNIZ[i] <-  data$IngroupOmegaZ[i] / data$FullOmegaZ[i]
            data$SexLinkedOmegaZ[i] <-  (geneinfo$Beta[which(geneinfo$Branch == "BranchZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchSkyZ")]) / (geneinfo$Alpha[which(geneinfo$Branch == "BranchZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyZ")])
            data$SexLinkedNIZ[i] <- data$IngroupOmegaZ[i] / data$SexLinkedOmegaZ[i]
            data$SexLinkedAlphaA[i] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchA")]
            data$SexLinkedAlphaZ[i] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyZ")]
            data$SexLinkedBetaA[i] <- geneinfo$Beta[which(geneinfo$Branch == "BranchA")]
            data$SexLinkedBetaZ[i] <- geneinfo$Beta[which(geneinfo$Branch == "BranchZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchSkyZ")]
            data$SexLinkedOmegaA[i] <- geneinfo$Omega[which(geneinfo$Branch == "BranchA")]
          } else {
            data$FullAlphaZ[i] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchAZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyZ")]
            data$FullBetaZ[i] <- geneinfo$Beta[which(geneinfo$Branch == "BranchAZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchSkyZ")]
            data$FullDSZ[i] <- geneinfo$DS[which(geneinfo$Branch == "BranchAZ")] + geneinfo$DS[which(geneinfo$Branch == "BranchSkyZ")]
            data$FullOmegaZ[i] <- (geneinfo$Beta[which(geneinfo$Branch == "BranchAZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchSkyZ")]) / (geneinfo$Alpha[which(geneinfo$Branch == "BranchAZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyZ")])
            data$FullNIZ[i] <-  data$IngroupOmegaZ[i] / data$FullOmegaZ[i]
            data$FullDoSZ[i] <- data$FullBetaZ[i]/(data$FullBetaZ[i] + data$FullAlphaZ[i]) - data$PNZ[i]/(data$PNZ[i] + data$PSZ[i])
            data$SexLinkedDoSZ[i] <- data$SexLinkedBetaZ[i]/(data$SexLinkedBetaZ[i] + data$SexLinkedAlphaZ[i]) - data$PNZ[i]/(data$PNZ[i] + data$PSZ[i])
          }

        } else {
          if(data$Region2[j] == "sex_phase" || data$Region2[j] == "partial_dropout") {
            data$FullAlphaZ[i] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyZ")]
            data$FullBetaZ[i] <- geneinfo$Beta[which(geneinfo$Branch == "BranchZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchSkyZ")]
            data$FullDSZ[i] <- geneinfo$DS[which(geneinfo$Branch == "BranchZ")] + geneinfo$DS[which(geneinfo$Branch == "BranchSkyZ")]
            data$FullOmegaZ[i] <- (geneinfo$Beta[which(geneinfo$Branch == "BranchZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchSkyZ")]) / (geneinfo$Alpha[which(geneinfo$Branch == "BranchZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyZ")])
            data$DSO[i] <- geneinfo$DS[which(geneinfo$Branch == "BranchO")] + geneinfo$DS[which(geneinfo$Branch == "Node12")]
            data$AlphaO[i] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchO")] + geneinfo$Alpha[which(geneinfo$Branch == "Node12")]
            data$BetaO[i] <- geneinfo$Beta[which(geneinfo$Branch == "BranchO")] + geneinfo$Beta[which(geneinfo$Branch == "Node12")]
            data$OmegaO[i] <- data$BetaO[i] / data$AlphaO[i]
            data$FullNIZ[i] <-  data$IngroupOmegaZ[i] / data$FullOmegaZ[i]
            data$FullDoSZ[i] <- data$FullBetaZ[i]/(data$FullBetaZ[i] + data$FullAlphaZ[i]) - data$PNZ[i]/(data$PNZ[i] + data$PSZ[i])
            data$SexLinkedDoSZ[i] <- data$SexLinkedBetaZ[i]/(data$SexLinkedBetaZ[i] + data$SexLinkedAlphaZ[i]) - data$PNZ[i]/(data$PNZ[i] + data$PSZ[i])
          } else {
            data$FullAlphaZ[i] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchAZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyZ")]
            data$FullBetaZ[i] <- geneinfo$Beta[which(geneinfo$Branch == "BranchAZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchSkyZ")]
            data$FullDSZ[i] <- geneinfo$DS[which(geneinfo$Branch == "BranchAZ")] + geneinfo$DS[which(geneinfo$Branch == "BranchSkyZ")]
            data$FullOmegaZ[i] <- (geneinfo$Beta[which(geneinfo$Branch == "BranchAZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchSkyZ")]) / (geneinfo$Alpha[which(geneinfo$Branch == "BranchAZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyZ")])
            data$DSO[i] <- geneinfo$DS[which(geneinfo$Branch == "BranchO")] + geneinfo$DS[which(geneinfo$Branch == "Node1")]
            data$AlphaO[i] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchO")] + geneinfo$Alpha[which(geneinfo$Branch == "Node1")]
            data$BetaO[i] <- geneinfo$Beta[which(geneinfo$Branch == "BranchO")] + geneinfo$Beta[which(geneinfo$Branch == "Node1")]
            data$OmegaO[i] <- data$BetaO[i] / data$AlphaO[i]
            data$FullNIZ[i] <-  data$IngroupOmegaZ[i] / data$FullOmegaZ[i]
            data$FullDoSZ[i] <- data$FullBetaZ[i]/(data$FullBetaZ[i] + data$FullAlphaZ[i]) - data$PNZ[i]/(data$PNZ[i] + data$PSZ[i])
            data$SexLinkedDoSZ[i] <- data$SexLinkedBetaZ[i]/(data$SexLinkedBetaZ[i] + data$SexLinkedAlphaZ[i]) - data$PNZ[i]/(data$PNZ[i] + data$PSZ[i])
          }
        }

      } else if(data$Region1[i] == "sex_phase" || data$Region1[i] == "partial_dropout") {
        data$IngroupOmegaZ[i] <-  geneinfo$IngroupOmega[which(geneinfo$Branch == "BranchSkyZ")]
        data$IngroupOmegaW[i] <-  geneinfo$IngroupOmega[which(geneinfo$Branch == "BranchSkyW")]
        data$PNZ[i] <-  geneinfo$PN[which(geneinfo$Branch == "BranchSkyZ")]
        data$PNW[i] <-  geneinfo$PN[which(geneinfo$Branch == "BranchSkyW")]
        data$PSZ[i] <-  geneinfo$PS[which(geneinfo$Branch == "BranchSkyZ")]
        data$PSW[i] <-  geneinfo$PS[which(geneinfo$Branch == "BranchSkyW")]
        data$SpeciesAlphaZ[i] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyZ")]
        data$SpeciesAlphaW[i] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyW")]
        data$SpeciesBetaZ[i] <- geneinfo$Beta[which(geneinfo$Branch == "BranchSkyZ")]
        data$SpeciesBetaW[i] <- geneinfo$Beta[which(geneinfo$Branch == "BranchSkyW")]
        data$SpeciesOmegaZ[i] <- geneinfo$Omega[which(geneinfo$Branch == "BranchSkyZ")]
        data$SpeciesOmegaW[i] <- geneinfo$Omega[which(geneinfo$Branch == "BranchSkyW")]
        data$SpeciesNIZ[i] <- data$IngroupOmegaZ[i] / data$SpeciesOmegaZ[i]
        data$SpeciesNIW[i] <- data$IngroupOmegaW[i] / data$SpeciesOmegaW[i]
        data$SpeciesDoSZ[i] <- data$SpeciesBetaZ[i]/(data$SpeciesBetaZ[i] + data$SpeciesAlphaZ[i]) - data$PNZ[i]/(data$PNZ[i] + data$PSZ[i])
        data$SpeciesDoSW[i] <- data$SpeciesBetaW[i]/(data$SpeciesBetaW[i] + data$SpeciesAlphaW[i]) - data$PNW[i]/(data$PNW[i] + data$PSW[i])

        if(data$Strata1[i] != "Z") {
          data$FullAlphaZ[i] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchA")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyZ")]
          data$FullBetaZ[i] <- geneinfo$Beta[which(geneinfo$Branch == "BranchA")] + geneinfo$Beta[which(geneinfo$Branch == "BranchZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchSkyZ")]
          data$FullDSZ[i] <- geneinfo$DS[which(geneinfo$Branch == "BranchA")] + geneinfo$DS[which(geneinfo$Branch == "BranchZ")] + geneinfo$DS[which(geneinfo$Branch == "BranchSkyZ")]
          data$FullOmegaZ[i] <- (geneinfo$Beta[which(geneinfo$Branch == "BranchA")] + geneinfo$Beta[which(geneinfo$Branch == "BranchZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchSkyZ")]) / (geneinfo$Alpha[which(geneinfo$Branch == "BranchA")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyZ")])
          data$DSO[i] <- geneinfo$DS[which(geneinfo$Branch == "BranchO")] + geneinfo$DS[which(geneinfo$Branch == "Node1")]
          data$AlphaO[i] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchO")] + geneinfo$Alpha[which(geneinfo$Branch == "Node1")]
          data$BetaO[i] <- geneinfo$Beta[which(geneinfo$Branch == "BranchO")] + geneinfo$Beta[which(geneinfo$Branch == "Node1")]
          data$OmegaO[i] <- data$BetaO[i] / data$AlphaO[i]
          data$FullNIZ[i] <-  data$IngroupOmegaZ[i] / data$FullOmegaZ[i]
          data$SexLinkedOmegaZ[i] <-  (geneinfo$Beta[which(geneinfo$Branch == "BranchZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchSkyZ")]) / (geneinfo$Alpha[which(geneinfo$Branch == "BranchZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyZ")])
          data$SexLinkedNIZ[i] <- data$IngroupOmegaZ[i] / data$SexLinkedOmegaZ[i]
          data$SexLinkedAlphaA[i] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchA")]
          data$SexLinkedAlphaZ[i] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyZ")]
          data$SexLinkedBetaA[i] <- geneinfo$Beta[which(geneinfo$Branch == "BranchA")]
          data$SexLinkedBetaZ[i] <- geneinfo$Beta[which(geneinfo$Branch == "BranchZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchSkyZ")]
          data$SexLinkedOmegaA[i] <- geneinfo$Omega[which(geneinfo$Branch == "BranchA")]
          data$FullDoSZ[i] <- data$FullBetaZ[i]/(data$FullBetaZ[i] + data$FullAlphaZ[i]) - data$PNZ[i]/(data$PNZ[i] + data$PSZ[i])
          data$SexLinkedDoSZ[i] <- data$SexLinkedBetaZ[i]/(data$SexLinkedBetaZ[i] + data$SexLinkedAlphaZ[i]) - data$PNZ[i]/(data$PNZ[i] + data$PSZ[i])

          if(data$Region2[j] == "sex_phase" || data$Region2[j] == "partial_dropout") {
            data$FullAlphaW[i] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchA")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchW")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyW")]
            data$FullBetaW[i] <- geneinfo$Beta[which(geneinfo$Branch == "BranchA")] + geneinfo$Beta[which(geneinfo$Branch == "BranchW")] + geneinfo$Beta[which(geneinfo$Branch == "BranchSkyW")]
            data$FullDSW[i] <- geneinfo$DS[which(geneinfo$Branch == "BranchA")] + geneinfo$DS[which(geneinfo$Branch == "BranchW")] + geneinfo$DS[which(geneinfo$Branch == "BranchSkyW")]
            data$SexLinkedAlphaW[i] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchW")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyW")]
            data$SexLinkedAlphaZW[i] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchW")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyW")]
            data$SexLinkedBetaW[i] <- geneinfo$Beta[which(geneinfo$Branch == "BranchW")] + geneinfo$Beta[which(geneinfo$Branch == "BranchSkyW")]
            data$SexLinkedBetaZW[i] <- geneinfo$Beta[which(geneinfo$Branch == "BranchZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchSkyZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchW")] + geneinfo$Beta[which(geneinfo$Branch == "BranchSkyW")]
            data$FullOmegaW[i] <- (geneinfo$Beta[which(geneinfo$Branch == "BranchA")] + geneinfo$Beta[which(geneinfo$Branch == "BranchW")] + geneinfo$Beta[which(geneinfo$Branch == "BranchSkyW")]) / (geneinfo$Alpha[which(geneinfo$Branch == "BranchA")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchW")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyW")])
            data$FullNIW[i] <- data$IngroupOmegaW[i] / data$FullOmegaW[i]
            data$SexLinkedOmegaW[i] <- (geneinfo$Beta[which(geneinfo$Branch == "BranchW")] + geneinfo$Beta[which(geneinfo$Branch == "BranchSkyW")]) / (geneinfo$Alpha[which(geneinfo$Branch == "BranchW")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyW")])
            data$SexLinkedNIW[i] <- data$IngroupOmegaW[i] / data$SexLinkedOmegaW[i]
            data$SexLinkedOmegaZW[i] <- (geneinfo$Beta[which(geneinfo$Branch == "BranchZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchW")] + geneinfo$Beta[which(geneinfo$Branch == "BranchSkyZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchSkyW")]) / (geneinfo$Alpha[which(geneinfo$Branch == "BranchZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchW")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyW")])
            data$FullDoSW[i] <- data$FullBetaW[i]/(data$FullBetaW[i] + data$FullAlphaW[i]) - data$PNW[i]/(data$PNW[i] + data$PSW[i])
            data$SexLinkedDoSW[i] <- data$SexLinkedBetaW[i]/(data$SexLinkedBetaW[i] + data$SexLinkedAlphaW[i]) - data$PNW[i]/(data$PNW[i] + data$PSW[i])

          } else {
            data$FullAlphaW[i] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchA")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyW")]
            data$FullBetaW[i] <- geneinfo$Beta[which(geneinfo$Branch == "BranchA")] + geneinfo$Beta[which(geneinfo$Branch == "BranchSkyW")]
            data$FullDSW[i] <- geneinfo$DS[which(geneinfo$Branch == "BranchA")] + geneinfo$DS[which(geneinfo$Branch == "BranchSkyW")]
            data$SexLinkedAlphaW[i] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyW")]
            data$SexLinkedAlphaZW[i] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyW")]
            data$SexLinkedBetaW[i] <- geneinfo$Beta[which(geneinfo$Branch == "BranchSkyW")]
            data$SexLinkedBetaZW[i] <- geneinfo$Beta[which(geneinfo$Branch == "BranchZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchSkyZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchSkyW")]
            data$FullOmegaW[i] <- (geneinfo$Beta[which(geneinfo$Branch == "BranchA")] + geneinfo$Beta[which(geneinfo$Branch == "BranchSkyW")]) / (geneinfo$Alpha[which(geneinfo$Branch == "BranchA")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyW")])
            data$FullNIW[i] <-  data$IngroupOmegaW[i] / data$FullOmegaW[i]
            data$SexLinkedOmegaW[i] <- geneinfo$Beta[which(geneinfo$Branch == "BranchSkyW")] / geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyW")]
            data$SexLinkedNIW[i] <- data$IngroupOmegaW[i] / data$SexLinkedOmegaW[i]
            data$SexLinkedOmegaZW[i] <- (geneinfo$Beta[which(geneinfo$Branch == "BranchZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchSkyZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchSkyW")]) / (geneinfo$Alpha[which(geneinfo$Branch == "BranchZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyW")])
            data$FullDoSW[i] <- data$FullBetaW[i]/(data$FullBetaW[i] + data$FullAlphaW[i]) - data$PNW[i]/(data$PNW[i] + data$PSW[i])
            data$SexLinkedDoSW[i] <- data$SexLinkedBetaW[i]/(data$SexLinkedBetaW[i] + data$SexLinkedAlphaW[i]) - data$PNW[i]/(data$PNW[i] + data$PSW[i])
          }

        } else {
          data$FullAlphaZ[i] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyZ")]
          data$FullBetaZ[i] <- geneinfo$Beta[which(geneinfo$Branch == "BranchZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchSkyZ")]
          data$FullDSZ[i] <- geneinfo$DS[which(geneinfo$Branch == "BranchZ")] + geneinfo$DS[which(geneinfo$Branch == "BranchSkyZ")]
          data$FullOmegaZ[i] <- (geneinfo$Beta[which(geneinfo$Branch == "BranchZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchSkyZ")]) / (geneinfo$Alpha[which(geneinfo$Branch == "BranchZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyZ")])
          data$FullNIZ[i] <-  data$IngroupOmegaZ[i] / data$FullOmegaZ[i]
          data$FullDoSZ[i] <- data$FullBetaZ[i]/(data$FullBetaZ[i] + data$FullAlphaZ[i]) - data$PNZ[i]/(data$PNZ[i] + data$PSZ[i])
          data$SexLinkedDoSZ[i] <- data$SexLinkedBetaZ[i]/(data$SexLinkedBetaZ[i] + data$SexLinkedAlphaZ[i]) - data$PNZ[i]/(data$PNZ[i] + data$PSZ[i])

          if(data$Region2[j] == "sex_phase" || data$Region2[j] == "partial_dropout") {
            data$FullAlphaW[i] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchW")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyW")]
            data$FullBetaW[i] <- geneinfo$Beta[which(geneinfo$Branch == "BranchW")] + geneinfo$Beta[which(geneinfo$Branch == "BranchSkyW")]
            data$FullDSW[i] <- geneinfo$DS[which(geneinfo$Branch == "BranchW")] + geneinfo$DS[which(geneinfo$Branch == "BranchSkyW")]
            data$DSO[i] <- geneinfo$DS[which(geneinfo$Branch == "BranchO")] + geneinfo$DS[which(geneinfo$Branch == "Node24")]
            data$AlphaO[i] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchO")] + geneinfo$Alpha[which(geneinfo$Branch == "Node24")]
            data$BetaO[i] <- geneinfo$Beta[which(geneinfo$Branch == "BranchO")] + geneinfo$Beta[which(geneinfo$Branch == "Node24")]
            data$OmegaO[i] <- data$BetaO[i] / data$AlphaO[i]
            data$SexLinkedAlphaZW[i] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchW")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyW")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchO")] + geneinfo$Alpha[which(geneinfo$Branch == "Node23")]
            data$SexLinkedBetaZW[i] <- geneinfo$Beta[which(geneinfo$Branch == "BranchZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchSkyZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchW")] + geneinfo$Beta[which(geneinfo$Branch == "BranchSkyW")] + geneinfo$Beta[which(geneinfo$Branch == "BranchO")] + geneinfo$Beta[which(geneinfo$Branch == "Node23")]
            data$SexLinkedOmegaZW[i] <- data$SexLinkedBetaZW[i] / data$SexLinkedAlphaZW[i]
            data$FullOmegaW[i] <- (geneinfo$Beta[which(geneinfo$Branch == "BranchW")] + geneinfo$Beta[which(geneinfo$Branch == "BranchSkyW")]) / (geneinfo$Alpha[which(geneinfo$Branch == "BranchW")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyW")])
            data$FullNIW[i] <- data$IngroupOmegaW[i] / data$FullOmegaW[i]
            data$FullDoSW[i] <- data$FullBetaW[i]/(data$FullBetaW[i] + data$FullAlphaW[i]) - data$PNW[i]/(data$PNW[i] + data$PSW[i])
            data$SexLinkedDoSW[i] <- data$SexLinkedBetaW[i]/(data$SexLinkedBetaW[i] + data$SexLinkedAlphaW[i]) - data$PNW[i]/(data$PNW[i] + data$PSW[i])

          } else {
            data$FullAlphaW[i] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyW")]
            data$FullBetaW[i] <- geneinfo$Beta[which(geneinfo$Branch == "BranchSkyW")]
            data$FullDSW[i] <- geneinfo$DS[which(geneinfo$Branch == "BranchSkyW")]
            data$DSO[i] <- geneinfo$DS[which(geneinfo$Branch == "BranchO")] + geneinfo$DS[which(geneinfo$Branch == "Node13")]
            data$AlphaO[i] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchO")] + geneinfo$Alpha[which(geneinfo$Branch == "Node13")]
            data$BetaO[i] <- geneinfo$Beta[which(geneinfo$Branch == "BranchO")] + geneinfo$Beta[which(geneinfo$Branch == "Node13")]
            data$OmegaO[i] <- data$BetaO[i] / data$AlphaO[i]
            data$SexLinkedAlphaZW[i] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyW")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchO")] + geneinfo$Alpha[which(geneinfo$Branch == "Node12")]
            data$SexLinkedBetaZW[i] <- geneinfo$Beta[which(geneinfo$Branch == "BranchZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchSkyZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchSkyW")] + geneinfo$Beta[which(geneinfo$Branch == "BranchO")] + geneinfo$Beta[which(geneinfo$Branch == "Node12")]
            data$SexLinkedOmegaZW[i] <- data$SexLinkedBetaZW[i] / data$SexLinkedAlphaZW[i]
            data$FullOmegaW[i] <- geneinfo$Beta[which(geneinfo$Branch == "BranchSkyW")] / geneinfo$Alpha[which(geneinfo$Branch == "BranchSkyW")]
            data$FullNIW[i] <-  data$IngroupOmegaW[i] / data$FullOmegaW[i]
            data$FullDoSW[i] <- data$FullBetaW[i]/(data$FullBetaW[i] + data$FullAlphaW[i]) - data$PNW[i]/(data$PNW[i] + data$PSW[i])
            data$SexLinkedDoSW[i] <- data$SexLinkedBetaW[i]/(data$SexLinkedBetaW[i] + data$SexLinkedAlphaW[i]) - data$PNW[i]/(data$PNW[i] + data$PSW[i])
          }
        }
      }

      data$geneMidPos[j] <- geneMin + ((geneMax-geneMin) / 2)
      if(data$Region2[j] == "autosomal") {
        data$IngroupOmegaA[j] <-  geneinfo$IngroupOmega[which(geneinfo$Branch == "BranchRasoA")]
        data$PNA[j] <-  geneinfo$PN[which(geneinfo$Branch == "BranchRasoA")]
        data$PSA[j] <-  geneinfo$PS[which(geneinfo$Branch == "BranchRasoA")]
        data$FullAlphaA[j] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchA")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoA")]
        data$FullBetaA[j] <- geneinfo$Beta[which(geneinfo$Branch == "BranchA")] + geneinfo$Beta[which(geneinfo$Branch == "BranchRasoA")]
        data$FullDSA[j] <- geneinfo$DS[which(geneinfo$Branch == "BranchA")] + geneinfo$DS[which(geneinfo$Branch == "BranchRasoA")]
        data$FullOmegaA[j] <- (geneinfo$Beta[which(geneinfo$Branch == "BranchA")] + geneinfo$Beta[which(geneinfo$Branch == "BranchRasoA")]) / (geneinfo$Alpha[which(geneinfo$Branch == "BranchA")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoA")])
        data$FullNIA[j] <-  data$IngroupOmegaA[j] / data$FullOmegaA[j]
        data$DSO[j] <- geneinfo$DS[which(geneinfo$Branch == "BranchO")] + geneinfo$DS[which(geneinfo$Branch == "Node1")]
        data$AlphaO[j] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchO")] + geneinfo$Alpha[which(geneinfo$Branch == "Node1")]
        data$BetaO[j] <- geneinfo$Beta[which(geneinfo$Branch == "BranchO")] + geneinfo$Beta[which(geneinfo$Branch == "Node1")]
        data$OmegaO[j] <- data$BetaO[j] / data$AlphaO[j]
        data$SpeciesAlphaA[j] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoA")]
        data$SpeciesBetaA[j] <- geneinfo$Beta[which(geneinfo$Branch == "BranchRasoA")]
        data$SpeciesOmegaA[j] <- geneinfo$Omega[which(geneinfo$Branch == "BranchRasoA")]
        data$SpeciesNIA[j] <- data$IngroupOmegaA[j] / data$SpeciesOmegaA[j]
        data$FullDoSA[j] <- data$FullBetaA[j]/(data$FullBetaA[j] + data$FullAlphaA[j]) - data$PNA[j]/(data$PNA[j] + data$PSA[j])
        data$SpeciesDoSA[j] <- data$SpeciesBetaA[j]/(data$SpeciesBetaA[j] + data$SpeciesAlphaA[j]) - data$PNA[j]/(data$PNA[j] + data$PSA[j])

      } else if(data$Region2[j] == "sex_dropout") {
        data$IngroupOmegaZ[j] <-  geneinfo$IngroupOmega[which(geneinfo$Branch == "BranchRasoZ")]
        data$PNZ[j] <-  geneinfo$PN[which(geneinfo$Branch == "BranchRasoZ")]
        data$PSZ[j] <-  geneinfo$PS[which(geneinfo$Branch == "BranchRasoZ")]
        data$SpeciesAlphaZ[j] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoZ")]
        data$SpeciesBetaZ[j] <- geneinfo$Beta[which(geneinfo$Branch == "BranchRasoZ")]
        data$SpeciesOmegaZ[j] <- geneinfo$Omega[which(geneinfo$Branch == "BranchRasoZ")]
        data$SpeciesNIZ[j] <- data$IngroupOmegaZ[j] / data$SpeciesOmegaZ[j]
        data$SpeciesDoSZ[j] <- data$SpeciesBetaZ[j]/(data$SpeciesBetaZ[j] + data$SpeciesAlphaZ[j]) - data$PNZ[j]/(data$PNZ[j] + data$PSZ[j])

        if(data$Strata1[j] != "Z") {
          if(data$Region1[i] == "sex_phase" || data$Region1[i] == "partial_dropout") {
            data$FullAlphaZ[j] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchA")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoZ")]
            data$FullBetaZ[j] <- geneinfo$Beta[which(geneinfo$Branch == "BranchA")] + geneinfo$Beta[which(geneinfo$Branch == "BranchZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchRasoZ")]
            data$FullDSZ[j] <- geneinfo$DS[which(geneinfo$Branch == "BranchA")] + geneinfo$DS[which(geneinfo$Branch == "BranchZ")] + geneinfo$DS[which(geneinfo$Branch == "BranchRasoZ")]
            data$FullOmegaZ[j] <- (geneinfo$Beta[which(geneinfo$Branch == "BranchA")] + geneinfo$Beta[which(geneinfo$Branch == "BranchZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchRasoZ")]) / (geneinfo$Alpha[which(geneinfo$Branch == "BranchA")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoZ")])
            data$DSO[j] <- geneinfo$DS[which(geneinfo$Branch == "BranchO")] + geneinfo$DS[which(geneinfo$Branch == "Node1")]
            data$AlphaO[j] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchO")] + geneinfo$Alpha[which(geneinfo$Branch == "Node1")]
            data$BetaO[j] <- geneinfo$Beta[which(geneinfo$Branch == "BranchO")] + geneinfo$Beta[which(geneinfo$Branch == "Node1")]
            data$OmegaO[j] <- data$BetaO[j] / data$AlphaO[j]
            data$FullNIZ[j] <-  data$IngroupOmegaZ[j] / data$FullOmegaZ[j]
            data$SexLinkedOmegaZ[j] <-  (geneinfo$Beta[which(geneinfo$Branch == "BranchZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchRasoZ")]) / (geneinfo$Alpha[which(geneinfo$Branch == "BranchZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoZ")])
            data$SexLinkedNIZ[j] <- data$IngroupOmegaZ[j] / data$SexLinkedOmegaZ[j]
            data$SexLinkedAlphaA[j] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchA")]
            data$SexLinkedAlphaZ[j] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoZ")]
            data$SexLinkedBetaA[j] <- geneinfo$Beta[which(geneinfo$Branch == "BranchA")]
            data$SexLinkedBetaZ[j] <- geneinfo$Beta[which(geneinfo$Branch == "BranchZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchRasoZ")]
            data$SexLinkedOmegaA[j] <- geneinfo$Omega[which(geneinfo$Branch == "BranchA")]
          } else {
            data$FullAlphaZ[j] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchAZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoZ")]
            data$FullBetaZ[j] <- geneinfo$Beta[which(geneinfo$Branch == "BranchAZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchRasoZ")]
            data$FullDSZ[j] <- geneinfo$DS[which(geneinfo$Branch == "BranchAZ")] + geneinfo$DS[which(geneinfo$Branch == "BranchRasoZ")]
            data$FullOmegaZ[j] <- (geneinfo$Beta[which(geneinfo$Branch == "BranchAZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchRasoZ")]) / (geneinfo$Alpha[which(geneinfo$Branch == "BranchAZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoZ")])
            data$DSO[j] <- geneinfo$DS[which(geneinfo$Branch == "BranchO")] + geneinfo$DS[which(geneinfo$Branch == "Node1")]
            data$AlphaO[j] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchO")] + geneinfo$Alpha[which(geneinfo$Branch == "Node1")]
            data$BetaO[j] <- geneinfo$Beta[which(geneinfo$Branch == "BranchO")] + geneinfo$Beta[which(geneinfo$Branch == "Node1")]
            data$OmegaO[j] <- data$BetaO[j] / data$AlphaO[j]
            data$FullNIZ[j] <-  data$IngroupOmegaZ[j] / data$FullOmegaZ[j]
            data$FullDoSZ[j] <- data$FullBetaZ[j]/(data$FullBetaZ[j] + data$FullAlphaZ[j]) - data$PNZ[j]/(data$PNZ[j] + data$PSZ[j])
            data$SexLinkedDoSZ[j] <- data$SexLinkedBetaZ[j]/(data$SexLinkedBetaZ[j] + data$SexLinkedAlphaZ[j]) - data$PNZ[j]/(data$PNZ[j] + data$PSZ[j])
          }

        } else {
          if(data$Region1[i] == "sex_phase" || data$Region1[i] == "partial_dropout") {
            data$FullAlphaZ[j] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoZ")]
            data$FullBetaZ[j] <- geneinfo$Beta[which(geneinfo$Branch == "BranchZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchRasoZ")]
            data$FullDSZ[j] <- geneinfo$DS[which(geneinfo$Branch == "BranchZ")] + geneinfo$DS[which(geneinfo$Branch == "BranchRasoZ")]
            data$FullOmegaZ[j] <- (geneinfo$Beta[which(geneinfo$Branch == "BranchZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchRasoZ")]) / (geneinfo$Alpha[which(geneinfo$Branch == "BranchZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoZ")])
            data$DSO[j] <- geneinfo$DS[which(geneinfo$Branch == "BranchO")] + geneinfo$DS[which(geneinfo$Branch == "Node13")]
            data$AlphaO[j] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchO")] + geneinfo$Alpha[which(geneinfo$Branch == "Node13")]
            data$BetaO[j] <- geneinfo$Beta[which(geneinfo$Branch == "BranchO")] + geneinfo$Beta[which(geneinfo$Branch == "Node13")]
            data$OmegaO[j] <- data$BetaO[j] / data$AlphaO[j]
            data$FullNIZ[j] <-  data$IngroupOmegaZ[j] / data$FullOmegaZ[j]
            data$FullDoSZ[j] <- data$FullBetaZ[j]/(data$FullBetaZ[j] + data$FullAlphaZ[j]) - data$PNZ[j]/(data$PNZ[j] + data$PSZ[j])
            data$SexLinkedDoSZ[j] <- data$SexLinkedBetaZ[j]/(data$SexLinkedBetaZ[j] + data$SexLinkedAlphaZ[j]) - data$PNZ[j]/(data$PNZ[j] + data$PSZ[j])
          } else {
            data$FullAlphaZ[j] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchAZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoZ")]
            data$FullBetaZ[j] <- geneinfo$Beta[which(geneinfo$Branch == "BranchAZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchRasoZ")]
            data$FullDSZ[j] <- geneinfo$DS[which(geneinfo$Branch == "BranchAZ")] + geneinfo$DS[which(geneinfo$Branch == "BranchRasoZ")]
            data$FullOmegaZ[j] <- (geneinfo$Beta[which(geneinfo$Branch == "BranchAZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchRasoZ")]) / (geneinfo$Alpha[which(geneinfo$Branch == "BranchAZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoZ")])
            data$DSO[j] <- geneinfo$DS[which(geneinfo$Branch == "BranchO")] + geneinfo$DS[which(geneinfo$Branch == "Node1")]
            data$AlphaO[j] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchO")] + geneinfo$Alpha[which(geneinfo$Branch == "Node1")]
            data$BetaO[j] <- geneinfo$Beta[which(geneinfo$Branch == "BranchO")] + geneinfo$Beta[which(geneinfo$Branch == "Node1")]
            data$OmegaO[j] <- data$BetaO[j] / data$AlphaO[j]
            data$FullNIZ[j] <-  data$IngroupOmegaZ[j] / data$FullOmegaZ[j]
            data$FullDoSZ[j] <- data$FullBetaZ[j]/(data$FullBetaZ[j] + data$FullAlphaZ[j]) - data$PNZ[j]/(data$PNZ[j] + data$PSZ[j])
            data$SexLinkedDoSZ[j] <- data$SexLinkedBetaZ[j]/(data$SexLinkedBetaZ[j] + data$SexLinkedAlphaZ[j]) - data$PNZ[j]/(data$PNZ[j] + data$PSZ[j])
          }
        }

      } else if(data$Region2[j] == "sex_phase" || data$Region2[j] == "partial_dropout") {
        data$IngroupOmegaZ[j] <-  geneinfo$IngroupOmega[which(geneinfo$Branch == "BranchRasoZ")]
        data$IngroupOmegaW[j] <-  geneinfo$IngroupOmega[which(geneinfo$Branch == "BranchRasoW")]
        data$PNZ[j] <-  geneinfo$PN[which(geneinfo$Branch == "BranchRasoZ")]
        data$PNW[j] <-  geneinfo$PN[which(geneinfo$Branch == "BranchRasoW")]
        data$PSZ[j] <-  geneinfo$PS[which(geneinfo$Branch == "BranchRasoZ")]
        data$PSW[j] <-  geneinfo$PS[which(geneinfo$Branch == "BranchRasoW")]
        data$SpeciesAlphaZ[j] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoZ")]
        data$SpeciesAlphaW[j] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoW")]
        data$SpeciesBetaZ[j] <- geneinfo$Beta[which(geneinfo$Branch == "BranchRasoZ")]
        data$SpeciesBetaW[j] <- geneinfo$Beta[which(geneinfo$Branch == "BranchRasoW")]
        data$SpeciesOmegaZ[j] <- geneinfo$Omega[which(geneinfo$Branch == "BranchRasoZ")]
        data$SpeciesOmegaW[j] <- geneinfo$Omega[which(geneinfo$Branch == "BranchRasoW")]
        data$SpeciesNIZ[j] <- data$IngroupOmegaZ[j] / data$SpeciesOmegaZ[j]
        data$SpeciesNIW[j] <- data$IngroupOmegaW[j] / data$SpeciesOmegaW[j]
        data$SpeciesDoSZ[j] <- data$SpeciesBetaZ[j]/(data$SpeciesBetaZ[j] + data$SpeciesAlphaZ[j]) - data$PNZ[j]/(data$PNZ[j] + data$PSZ[j])
        data$SpeciesDoSW[j] <- data$SpeciesBetaW[j]/(data$SpeciesBetaW[j] + data$SpeciesAlphaW[j]) - data$PNW[j]/(data$PNW[j] + data$PSW[j])

        if(data$Strata1[j] != "Z") {
          data$FullAlphaZ[j] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchA")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoZ")]
          data$FullBetaZ[j] <- geneinfo$Beta[which(geneinfo$Branch == "BranchA")] + geneinfo$Beta[which(geneinfo$Branch == "BranchZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchRasoZ")]
          data$FullDSZ[j] <- geneinfo$DS[which(geneinfo$Branch == "BranchA")] + geneinfo$DS[which(geneinfo$Branch == "BranchZ")] + geneinfo$DS[which(geneinfo$Branch == "BranchRasoZ")]
          data$FullOmegaZ[j] <- (geneinfo$Beta[which(geneinfo$Branch == "BranchA")] + geneinfo$Beta[which(geneinfo$Branch == "BranchZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchRasoZ")]) / (geneinfo$Alpha[which(geneinfo$Branch == "BranchA")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoZ")])
          data$DSO[j] <- geneinfo$DS[which(geneinfo$Branch == "BranchO")] + geneinfo$DS[which(geneinfo$Branch == "Node1")]
          data$AlphaO[j] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchO")] + geneinfo$Alpha[which(geneinfo$Branch == "Node1")]
          data$BetaO[j] <- geneinfo$Beta[which(geneinfo$Branch == "BranchO")] + geneinfo$Beta[which(geneinfo$Branch == "Node1")]
          data$OmegaO[j] <- data$BetaO[j] / data$AlphaO[j]
          data$FullNIZ[j] <-  data$IngroupOmegaZ[j] / data$FullOmegaZ[j]
          data$SexLinkedOmegaZ[j] <-  (geneinfo$Beta[which(geneinfo$Branch == "BranchZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchRasoZ")]) / (geneinfo$Alpha[which(geneinfo$Branch == "BranchZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoZ")])
          data$SexLinkedNIZ[j] <- data$IngroupOmegaZ[j] / data$SexLinkedOmegaZ[j]
          data$SexLinkedAlphaA[j] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchA")]
          data$SexLinkedAlphaZ[j] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoZ")]
          data$SexLinkedBetaA[j] <- geneinfo$Beta[which(geneinfo$Branch == "BranchA")]
          data$SexLinkedBetaZ[j] <- geneinfo$Beta[which(geneinfo$Branch == "BranchZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchRasoZ")]
          data$SexLinkedOmegaA[j] <- geneinfo$Omega[which(geneinfo$Branch == "BranchA")]
          data$FullDoSZ[j] <- data$FullBetaZ[j]/(data$FullBetaZ[j] + data$FullAlphaZ[j]) - data$PNZ[j]/(data$PNZ[j] + data$PSZ[j])
          data$SexLinkedDoSZ[j] <- data$SexLinkedBetaZ[j]/(data$SexLinkedBetaZ[j] + data$SexLinkedAlphaZ[j]) - data$PNZ[j]/(data$PNZ[j] + data$PSZ[j])

          if(data$Region1[i] == "sex_phase" || data$Region1[i] == "partial_dropout") {
            data$FullAlphaW[j] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchA")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchW")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoW")]
            data$FullBetaW[j] <- geneinfo$Beta[which(geneinfo$Branch == "BranchA")] + geneinfo$Beta[which(geneinfo$Branch == "BranchW")] + geneinfo$Beta[which(geneinfo$Branch == "BranchRasoW")]
            data$FullDSW[j] <- geneinfo$DS[which(geneinfo$Branch == "BranchA")] + geneinfo$DS[which(geneinfo$Branch == "BranchW")] + geneinfo$DS[which(geneinfo$Branch == "BranchRasoW")]
            data$SexLinkedAlphaW[j] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchW")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoW")]
            data$SexLinkedAlphaZW[j] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchW")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoW")]
            data$SexLinkedBetaW[j] <- geneinfo$Beta[which(geneinfo$Branch == "BranchW")] + geneinfo$Beta[which(geneinfo$Branch == "BranchRasoW")]
            data$SexLinkedBetaZW[j] <- geneinfo$Beta[which(geneinfo$Branch == "BranchZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchRasoZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchW")] + geneinfo$Beta[which(geneinfo$Branch == "BranchRasoW")]
            data$FullOmegaW[j] <- (geneinfo$Beta[which(geneinfo$Branch == "BranchA")] + geneinfo$Beta[which(geneinfo$Branch == "BranchW")] + geneinfo$Beta[which(geneinfo$Branch == "BranchRasoW")]) / (geneinfo$Alpha[which(geneinfo$Branch == "BranchA")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchW")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoW")])
            data$FullNIW[j] <- data$IngroupOmegaW[j] / data$FullOmegaW[j]
            data$SexLinkedOmegaW[j] <- (geneinfo$Beta[which(geneinfo$Branch == "BranchW")] + geneinfo$Beta[which(geneinfo$Branch == "BranchRasoW")]) / (geneinfo$Alpha[which(geneinfo$Branch == "BranchW")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoW")])
            data$SexLinkedNIW[j] <- data$IngroupOmegaW[j] / data$SexLinkedOmegaW[j]
            data$SexLinkedOmegaZW[j] <- (geneinfo$Beta[which(geneinfo$Branch == "BranchZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchW")] + geneinfo$Beta[which(geneinfo$Branch == "BranchRasoZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchRasoW")]) / (geneinfo$Alpha[which(geneinfo$Branch == "BranchZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchW")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoW")])
            data$FullDoSW[j] <- data$FullBetaW[j]/(data$FullBetaW[j] + data$FullAlphaW[j]) - data$PNW[j]/(data$PNW[j] + data$PSW[j])
            data$SexLinkedDoSW[j] <- data$SexLinkedBetaW[j]/(data$SexLinkedBetaW[j] + data$SexLinkedAlphaW[j]) - data$PNW[j]/(data$PNW[j] + data$PSW[j])

          } else {
            data$FullAlphaW[j] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchA")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoW")]
            data$FullBetaW[j] <- geneinfo$Beta[which(geneinfo$Branch == "BranchA")] + geneinfo$Beta[which(geneinfo$Branch == "BranchRasoW")]
            data$FullDSW[j] <- geneinfo$DS[which(geneinfo$Branch == "BranchA")] + geneinfo$DS[which(geneinfo$Branch == "BranchRasoW")]
            data$SexLinkedAlphaW[j] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoW")]
            data$SexLinkedAlphaZW[j] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoW")]
            data$SexLinkedBetaW[j] <- geneinfo$Beta[which(geneinfo$Branch == "BranchRasoW")]
            data$SexLinkedBetaZW[j] <- geneinfo$Beta[which(geneinfo$Branch == "BranchZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchRasoZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchRasoW")]
            data$FullOmegaW[j] <- (geneinfo$Beta[which(geneinfo$Branch == "BranchA")] + geneinfo$Beta[which(geneinfo$Branch == "BranchRasoW")]) / (geneinfo$Alpha[which(geneinfo$Branch == "BranchA")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoW")])
            data$FullNIW[j] <-  data$IngroupOmegaW[j] / data$FullOmegaW[j]
            data$SexLinkedOmegaW[j] <- geneinfo$Beta[which(geneinfo$Branch == "BranchRasoW")] / geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoW")]
            data$SexLinkedNIW[j] <- data$IngroupOmegaW[j] / data$SexLinkedOmegaW[j]
            data$SexLinkedOmegaZW[j] <- (geneinfo$Beta[which(geneinfo$Branch == "BranchZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchRasoZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchRasoW")]) / (geneinfo$Alpha[which(geneinfo$Branch == "BranchZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoW")])
            data$FullDoSW[j] <- data$FullBetaW[j]/(data$FullBetaW[j] + data$FullAlphaW[j]) - data$PNW[j]/(data$PNW[j] + data$PSW[j])
            data$SexLinkedDoSW[j] <- data$SexLinkedBetaW[j]/(data$SexLinkedBetaW[j] + data$SexLinkedAlphaW[j]) - data$PNW[j]/(data$PNW[j] + data$PSW[j])
          }

        } else {
          data$FullAlphaZ[j] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoZ")]
          data$FullBetaZ[j] <- geneinfo$Beta[which(geneinfo$Branch == "BranchZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchRasoZ")]
          data$FullDSZ[j] <- geneinfo$DS[which(geneinfo$Branch == "BranchZ")] + geneinfo$DS[which(geneinfo$Branch == "BranchRasoZ")]
          data$FullOmegaZ[j] <- (geneinfo$Beta[which(geneinfo$Branch == "BranchZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchRasoZ")]) / (geneinfo$Alpha[which(geneinfo$Branch == "BranchZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoZ")])
          data$FullNIZ[j] <-  data$IngroupOmegaZ[j] / data$FullOmegaZ[j]
          data$FullDoSZ[j] <- data$FullBetaZ[j]/(data$FullBetaZ[j] + data$FullAlphaZ[j]) - data$PNZ[j]/(data$PNZ[j] + data$PSZ[j])
          data$SexLinkedDoSZ[j] <- data$SexLinkedBetaZ[j]/(data$SexLinkedBetaZ[j] + data$SexLinkedAlphaZ[j]) - data$PNZ[j]/(data$PNZ[j] + data$PSZ[j])

          if(data$Region1[i] == "sex_phase" || data$Region1[i] == "partial_dropout") {
            data$FullAlphaW[j] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchW")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoW")]
            data$FullBetaW[j] <- geneinfo$Beta[which(geneinfo$Branch == "BranchW")] + geneinfo$Beta[which(geneinfo$Branch == "BranchRasoW")]
            data$FullDSW[j] <- geneinfo$DS[which(geneinfo$Branch == "BranchW")] + geneinfo$DS[which(geneinfo$Branch == "BranchRasoW")]
            data$DSO[j] <- geneinfo$DS[which(geneinfo$Branch == "BranchO")] + geneinfo$DS[which(geneinfo$Branch == "Node24")]
            data$AlphaO[j] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchO")] + geneinfo$Alpha[which(geneinfo$Branch == "Node24")]
            data$BetaO[j] <- geneinfo$Beta[which(geneinfo$Branch == "BranchO")] + geneinfo$Beta[which(geneinfo$Branch == "Node24")]
            data$OmegaO[j] <- data$BetaO[j] / data$AlphaO[j]
            data$SexLinkedAlphaZW[j] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchW")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoW")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchO")] + geneinfo$Alpha[which(geneinfo$Branch == "Node24")]
            data$SexLinkedBetaZW[j] <- geneinfo$Beta[which(geneinfo$Branch == "BranchZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchRasoZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchW")] + geneinfo$Beta[which(geneinfo$Branch == "BranchRasoW")] + geneinfo$Beta[which(geneinfo$Branch == "BranchO")] + geneinfo$Beta[which(geneinfo$Branch == "Node24")]
            data$SexLinkedOmegaZW[j] <-   data$SexLinkedBetaZW[j] / data$SexLinkedAlphaZW[j]
            data$FullOmegaW[j] <- (geneinfo$Beta[which(geneinfo$Branch == "BranchW")] + geneinfo$Beta[which(geneinfo$Branch == "BranchRasoW")]) / (geneinfo$Alpha[which(geneinfo$Branch == "BranchW")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoW")])
            data$FullNIW[j] <- data$IngroupOmegaW[j] / data$FullOmegaW[j]
            data$FullDoSW[j] <- data$FullBetaW[j]/(data$FullBetaW[j] + data$FullAlphaW[j]) - data$PNW[j]/(data$PNW[j] + data$PSW[j])
            data$SexLinkedDoSW[j] <- data$SexLinkedBetaW[j]/(data$SexLinkedBetaW[j] + data$SexLinkedAlphaW[j]) - data$PNW[j]/(data$PNW[j] + data$PSW[j])

          } else {
            data$FullAlphaW[j] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoW")]
            data$FullBetaW[j] <- geneinfo$Beta[which(geneinfo$Branch == "BranchRasoW")]
            data$FullDSW[j] <- geneinfo$DS[which(geneinfo$Branch == "BranchRasoW")]
            data$DSO[j] <- geneinfo$DS[which(geneinfo$Branch == "BranchO")] + geneinfo$DS[which(geneinfo$Branch == "Node12")]
            data$AlphaO[j] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchO")] + geneinfo$Alpha[which(geneinfo$Branch == "Node12")]
            data$BetaO[j] <- geneinfo$Beta[which(geneinfo$Branch == "BranchO")] + geneinfo$Beta[which(geneinfo$Branch == "Node12")]
            data$OmegaO[j] <- data$BetaO[j] / data$AlphaO[j]
            data$SexLinkedAlphaZW[j] <- geneinfo$Alpha[which(geneinfo$Branch == "BranchZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoZ")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoW")] + geneinfo$Alpha[which(geneinfo$Branch == "BranchO")] + geneinfo$Alpha[which(geneinfo$Branch == "Node11")]
            data$SexLinkedBetaZW[j] <- geneinfo$Beta[which(geneinfo$Branch == "BranchZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchRasoZ")] + geneinfo$Beta[which(geneinfo$Branch == "BranchRasoW")] + geneinfo$Beta[which(geneinfo$Branch == "BranchO")] + geneinfo$Beta[which(geneinfo$Branch == "Node11")]
            data$SexLinkedOmegaZW[j] <-   data$SexLinkedBetaZW[j] / data$SexLinkedAlphaZW[j]
            data$FullOmegaW[j] <- geneinfo$Beta[which(geneinfo$Branch == "BranchRasoW")] / geneinfo$Alpha[which(geneinfo$Branch == "BranchRasoW")]
            data$FullNIW[j] <-  data$IngroupOmegaW[j] / data$FullOmegaW[j]
            data$FullDoSW[j] <- data$FullBetaW[j]/(data$FullBetaW[j] + data$FullAlphaW[j]) - data$PNW[j]/(data$PNW[j] + data$PSW[j])
            data$SexLinkedDoSW[j] <- data$SexLinkedBetaW[j]/(data$SexLinkedBetaW[j] + data$SexLinkedAlphaW[j]) - data$PNW[j]/(data$PNW[j] + data$PSW[j])
          }
        }
      }
    }
  }
}
data <- rbind(data, data_other)

### Add info for functional impacts on protein structure
data$FullImpactA <- as.numeric(rep(0, nrow(data)))
data$FullImpactZ <- as.numeric(rep(0, nrow(data)))
data$FullImpactW <- as.numeric(rep(0, nrow(data)))
data$SexLinkedImpactZ <- as.numeric(rep(0, nrow(data)))
data$SexLinkedImpactW <- as.numeric(rep(0, nrow(data)))
data$FullLOFA <- as.numeric(rep(0, nrow(data)))
data$FullLOFZ <- as.numeric(rep(0, nrow(data)))
data$FullLOFW <- as.numeric(rep(0, nrow(data)))
data$SexLinkedLOFZ <- as.numeric(rep(0, nrow(data)))
data$SexLinkedLOFW <- as.numeric(rep(0, nrow(data)))
data$snpEffWarnings <- as.numeric(rep(NA, nrow(data)))

# add high impacts
dataHigh <- read.delim(paste(WORKDIR1, "/D6_snpEff_substitutions/metadata/annotation_high_impacts.tsv", sep=""), head=F)

for(i in 1:(nrow(dataHigh))) {

  if(length(grep("Skylark", dataHigh[i,4])) > 0) {
    if(substr(dataHigh[i,4], nchar(dataHigh[i,4])-1, nchar(dataHigh[i,4])) == "_A") {
      data$FullImpactA[which(data$geneID == dataHigh[i,1] & data$Species == "Skylark")] <- dataHigh[i,3]
    } else if(substr(dataHigh[i,4], nchar(dataHigh[i,4])-1, nchar(dataHigh[i,4])) == "_Z")  {
      data$FullImpactZ[which(data$geneID == dataHigh[i,1] & data$Species == "Skylark")] <- dataHigh[i,3]
    } else if(substr(dataHigh[i,4], nchar(dataHigh[i,4])-1, nchar(dataHigh[i,4])) == "_W")  {
      data$FullImpactW[which(data$geneID == dataHigh[i,1] & data$Species == "Skylark")] <- dataHigh[i,3]
    } else if(substr(dataHigh[i,4], nchar(dataHigh[i,4])-11, nchar(dataHigh[i,4])) == "_Z_sexlinked")  {
      data$SexLinkedImpactZ[which(data$geneID == dataHigh[i,1] & data$Species == "Skylark")] <- dataHigh[i,3]
    } else if(substr(dataHigh[i,4], nchar(dataHigh[i,4])-11, nchar(dataHigh[i,4])) == "_W_sexlinked")  {
      data$SexLinkedImpactW[which(data$geneID == dataHigh[i,1] & data$Species == "Skylark")] <- dataHigh[i,3]
    }
  } else if(length(grep("Rasolark", dataHigh[i,4])) > 0) {
    if(substr(dataHigh[i,4], nchar(dataHigh[i,4])-1, nchar(dataHigh[i,4])) == "_A") {
      data$FullImpactA[which(data$geneID == dataHigh[i,1] & data$Species == "Rasolark")] <- dataHigh[i,3]
    } else if(substr(dataHigh[i,4], nchar(dataHigh[i,4])-1, nchar(dataHigh[i,4])) == "_Z")  {
      data$FullImpactZ[which(data$geneID == dataHigh[i,1] & data$Species == "Rasolark")] <- dataHigh[i,3]
    } else if(substr(dataHigh[i,4], nchar(dataHigh[i,4])-1, nchar(dataHigh[i,4])) == "_W")  {
      data$FullImpactW[which(data$geneID == dataHigh[i,1] & data$Species == "Rasolark")] <- dataHigh[i,3]
    } else if(substr(dataHigh[i,4], nchar(dataHigh[i,4])-11, nchar(dataHigh[i,4])) == "_Z_sexlinked")  {
      data$SexLinkedImpactZ[which(data$geneID == dataHigh[i,1] & data$Species == "Rasolark")] <- dataHigh[i,3]
    } else if(substr(dataHigh[i,4], nchar(dataHigh[i,4])-11, nchar(dataHigh[i,4])) == "_W_sexlinked")  {
      data$SexLinkedImpactW[which(data$geneID == dataHigh[i,1] & data$Species == "Rasolark")] <- dataHigh[i,3]
    }
  }
}


# add LOF impacts
dataLOF <- read.delim(paste(WORKDIR1, "/D6_snpEff_substitutions/metadata/annotation_LOF.tsv", sep=""), head=F)

for(i in 1:(nrow(dataLOF))) {

  if(length(grep("Skylark", dataLOF[i,2])) > 0) {
    if(substr(dataLOF[i,2], nchar(dataLOF[i,2])-1, nchar(dataLOF[i,2])) == "_A") {
      data$FullLOFA[which(data$geneID == dataLOF[i,1] & data$Species == "Skylark")] <- 1
    } else if(substr(dataLOF[i,2], nchar(dataLOF[i,2])-1, nchar(dataLOF[i,2])) == "_Z")  {
      data$FullLOFZ[which(data$geneID == dataLOF[i,1] & data$Species == "Skylark")] <- 1
    } else if(substr(dataLOF[i,2], nchar(dataLOF[i,2])-1, nchar(dataLOF[i,2])) == "_W")  {
      data$FullLOFW[which(data$geneID == dataLOF[i,1] & data$Species == "Skylark")] <- 1
    } else if(substr(dataLOF[i,2], nchar(dataLOF[i,2])-11, nchar(dataLOF[i,2])) == "_Z_sexlinked")  {
      data$SexLinkedLOFZ[which(data$geneID == dataLOF[i,1] & data$Species == "Skylark")] <- 1
    } else if(substr(dataLOF[i,2], nchar(dataLOF[i,2])-11, nchar(dataLOF[i,2])) == "_W_sexlinked")  {
      data$SexLinkedLOFW[which(data$geneID == dataLOF[i,1] & data$Species == "Skylark")] <- 1
    }
  } else if(length(grep("Rasolark", dataLOF[i,2])) > 0) {
    if(substr(dataLOF[i,2], nchar(dataLOF[i,2])-1, nchar(dataLOF[i,2])) == "_A") {
      data$FullLOFA[which(data$geneID == dataLOF[i,1] & data$Species == "Rasolark")] <- 1
    } else if(substr(dataLOF[i,2], nchar(dataLOF[i,2])-1, nchar(dataLOF[i,2])) == "_Z")  {
      data$FullLOFZ[which(data$geneID == dataLOF[i,1] & data$Species == "Rasolark")] <- 1
    } else if(substr(dataLOF[i,2], nchar(dataLOF[i,2])-1, nchar(dataLOF[i,2])) == "_W")  {
      data$FullLOFW[which(data$geneID == dataLOF[i,1] & data$Species == "Rasolark")] <- 1
    } else if(substr(dataLOF[i,2], nchar(dataLOF[i,2])-11, nchar(dataLOF[i,2])) == "_Z_sexlinked")  {
      data$SexLinkedLOFZ[which(data$geneID == dataLOF[i,1] & data$Species == "Rasolark")] <- 1
    } else if(substr(dataLOF[i,2], nchar(dataLOF[i,2])-11, nchar(dataLOF[i,2])) == "_W_sexlinked")  {
      data$SexLinkedLOFW[which(data$geneID == dataLOF[i,1] & data$Species == "Rasolark")] <- 1
    }
  }
}

# add snpEff warnings
datawarn <- read.delim(paste(WORKDIR1, "/D6_snpEff_substitutions/metadata/annotation_warnings.tsv", sep=""), head=F)

for(i in 1:(nrow(datawarn))) {
  if(length(which(!is.na(data$snpEffWarnings[which(data$geneID == datawarn[i,1])]) == T)) == 0) {
    data$snpEffWarnings[which(data$geneID == datawarn[i,1])] <- datawarn[i,2]
  } else {
    data$snpEffWarnings[which(data$geneID == datawarn[i,1])] <- paste(data$snpEffWarnings[which(data$geneID == datawarn[i,1])], datawarn[i,2], sep="&")
  }
}



data_other <- data[which(data[,14] != "Skylark_2021_Rasolark_2021"),]
data <- data[which(data[,14]== "Skylark_2021_Rasolark_2021"),]
data$Region1[which(data$Species == "Rasolark")] <- data$Region2[which(data$Species == "Rasolark")]
data$Strata1[which(data$Species == "Rasolark")] <- data$Strata2[which(data$Species == "Rasolark")]
data <- data[,-which(colnames(data) == "Region2" | colnames(data) == "Strata2")]
colnames(data)[which(colnames(data) == "Region1" | colnames(data) == "Strata1")] <- c("Region", "Strata")
data_other <- data_other[,-which(colnames(data) == "Region2" | colnames(data) == "Strata2")]
colnames(data_other)[which(colnames(data_other) == "Region1" | colnames(data_other) == "Strata1")] <- c("Region", "Strata")
data <- rbind(data, data_other)


data$Wdegeneration <- rep(NA, nrow(data))
data$Wdegeneration[which(data$FullLOFA > 0)] <- "W loss of function"
data$Wdegeneration[which(data$Region == "partial_dropout")] <- "W partially degenerated"
data$Wdegeneration[which(data$Region == "sex_dropout")] <- "W degenerated"
data$Wdegeneration[which(data$FullLOFW == 1)] <- "W loss of function"
data$Wdegeneration[which(is.na(data$Wdegeneration))] <- "W functional"

data$Zdegeneration <- rep(NA, nrow(data))
data$Zdegeneration[which(data$FullLOFZ == 1)] <- "Z loss of function"
data$Zdegeneration[which(is.na(data$Zdegeneration))] <- "Z functional"




# Add ancestral strata
data_ancestralZ <- as.data.frame(read.delim(paste(WORKDIR1, "/ancestralZ_strata/metadata/strata_match.tsv", sep=""), head=T))
data$Ancestral_Strata_Confidence_Rank <- rep(NA, nrow(data))

for(i in 1:(nrow(data_ancestralZ))) {
  if(!is.na(data_ancestralZ$Stratum[i])) {
    data$Strata[which(data$geneID == data_ancestralZ$Gene[i])] <-  data_ancestralZ$Stratum[i]
    data$Ancestral_Strata_Confidence_Rank[which(data$geneID == data_ancestralZ$Gene[i])] <- data_ancestralZ$Confidence_rank[i]
  }
}
data$Strata[which(data$Strata =="Z")] <- "Ancestral unknown"

# Add stratum age
data_Strata_Age <- as.data.frame(read.delim(paste(WORKDIR1, "/ancestral_generation_times/Strata_ages_generations.tsv", sep=""), head=T))
data$Strata_Age_Years <- as.numeric(rep(NA, nrow(data)))
data$Strata_Age_Generations <- as.numeric(rep(NA, nrow(data)))
for(i in 1:nrow(data_Strata_Age)) {
  data$Strata_Age_Years[which(data$Strata == data_Strata_Age$Strata[i])] <- data_Strata_Age$Age_years[i]
  data$Strata_Age_Generations[which(data$Strata == data_Strata_Age$Strata[i])] <- data_Strata_Age$Age_generations[i]
}

data$Strata_Age_Years[which(data$Strata == "PAR3")] <- 0
data$Strata_Age_Years[which(data$Strata == "PAR5")] <- 0
data$Strata_Age_Generations[which(data$Strata == "PAR3")] <- 0
data$Strata_Age_Generations[which(data$Strata == "PAR5")] <- 0

# Add haploinsufficiency scores
data_haploinsufficiency <- as.data.frame(read.delim(paste(WORKDIR1, "/haploinsufficiency_scores/haploinsufficiency_scores.tsv", sep=""), head=T))
colnames(data_haploinsufficiency) <- c("Gene", "Transcript", "Human_gene_name", "Human_gene_ID", "Method", "HI", "pHaplo")
data$HI <- rep(NA, nrow(data))
data$pHaplo <- rep(NA, nrow(data))
data$haploinsufficiency_Match_Method <- rep(NA, nrow(data))
data$Human_gene_ID <- rep(NA, nrow(data))

for(i in 1:(nrow(data_haploinsufficiency))) {
  if(!is.na(data_haploinsufficiency$Method[i])) {

    index <- which(data$geneID == data_haploinsufficiency$Gene[i])
    data$HI[index] <- data_haploinsufficiency$HI[i]
    data$pHaplo[index] <- data_haploinsufficiency$pHaplo[i]
    data$haploinsufficiency_Match_Method[index] <- data_haploinsufficiency$Method[i]
    data$Human_gene_ID[index] <- data_haploinsufficiency$Human_gene_ID[i]

  }

}

# Filter
# Filter1 Dont use for any rate of evolution estimates
data$Filter1 <- rep("OK", nrow(data))
data$Filter1[which(data$geneLengthTrimmed < 300)] <- "Too short"  # Trimmed length

# Filter1 Dont use for any rate of evolution estimates
data$Filter2 <- rep("OK", nrow(data))
data$Filter2[which(data$FullDSA < 0.01 | data$FullDSZ < 0.01 | data$FullDSW < 0.01 )] <- "Few synonymous substitutions"
data$Filter2[which(data$FullDSA == 0 | data$FullDSZ == 0 | data$FullDSW == 0)] <- "No synonymous substitutions"
data$Filter2[which(data$FullDSA >= 1 | data$FullDSZ >= 1 | data$FullDSW >= 1)] <- "Substitutions saturated"
data$Filter2[which((data$Region == "autosomal" & is.na(data$FullDSA)) | (data$Region != "autosomal" & (is.na(data$FullDSZ) | is.na(data$FullDSW))))] <- "No data"

# Filter 3 Dont use for anything
data$Filter3 <- rep("OK", nrow(data))
data$Filter3[grep("MULTIPLE_STOP_CODONS", data$snpEffWarnings)] <- "Multiple stop codons in DB"

# Filter 4 Dont use for anything
data$Filter4 <- rep("OK", nrow(data))
data$Filter4[which(data$Strata == "Ancestral unknown" | data$Strata == "PAR3unk" | data$Strata == "PAR5unk")] <- "Unknown Stratum"

# Filter 5 Dont use for anything
data$Filter5 <- rep("OK", nrow(data))
data$Filter5[which(data$Datasets != "Skylark_2021_Rasolark_2021")] <- "Gene not present in both datasets"


# Write organised table
write.table(data, paste(OUTDIR, "/", PROJECT1, "_", PROJECT2, "_organised_data1.tsv", sep=""), quote=F, sep='\t', row.names = F, col.names = T)
q()

