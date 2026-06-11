## Export variables and load libraries
rm(list=ls())
library(tidyverse)
library(patchwork)
library(ggridges)
sessionInfo()

options(scipen=999)
setwd("C:/Users/Simon JE/OneDrive - Lund University/Dokument/Simon/PhD/Projects/Rasolark_2021/Results")


### Rasolark
data_NRA <- read.delim("Neutral_popstats/Rasolark/autosomal_windows_10000_steps_10000_exon_dist_20000/Rasolark_2021_autosomal_windows_10000_steps_10000_exon_dist_20000_data.txt", sep="\t", head=T)
data_NRA$pi_abs[which(is.na(data_NRA$pi_abs))] <- 0
data_NRA <- data_NRA[which(data_NRA$N_callable_sites > 5000),]
data_NRA_mean <- mean(data_NRA$pi_abs)
data_NRZ <- read.delim("Neutral_popstats/Rasolark/homogametic_windows_10000_steps_10000_exon_dist_20000/Rasolark_2021_homogametic_windows_10000_steps_10000_exon_dist_20000_data.txt", sep="\t", head=T)
data_NRZ$pi_abs[which(is.na(data_NRZ$pi_abs))] <- 0
data_NRZ <- data_NRZ[which(data_NRZ$N_callable_sites > 5000),]
data_NRZ_mean <- mean(data_NRZ$pi_abs)
data_NRW <- read.delim("Neutral_popstats/Rasolark/heterogametic_windows_10000_steps_10000_exon_dist_20000/Rasolark_2021_heterogametic_windows_10000_steps_10000_exon_dist_20000_data.txt", sep="\t", head=T)
data_NRW$pi_abs[which(is.na(data_NRW$pi_abs))] <- 0
data_NRW <- data_NRW[which(data_NRW$N_callable_sites > 5000),]
data_NRW_mean <- mean(data_NRW$pi_abs)

data_ERA <- read.delim("0fold_popstats/Rasolark/autosomal_windows_10000_steps_10000_exon_dist_0fold/Rasolark_2021_autosomal_windows_10000_steps_10000_exon_dist_0fold_data.txt", sep="\t", head=T)
data_ERA$pi_abs[which(is.na(data_ERA$pi_abs))] <- 0
data_ERA <- data_ERA[which(data_ERA$N_callable_sites >= 5000),]
data_ERA_mean <- mean(data_ERA$pi_abs)
data_ERZ <- read.delim("0fold_popstats/Rasolark/homogametic_windows_10000_steps_10000_exon_dist_0fold/Rasolark_2021_homogametic_windows_10000_steps_10000_exon_dist_0fold_data.txt", sep="\t", head=T)
data_ERZ$pi_abs[which(is.na(data_ERZ$pi_abs))] <- 0
data_ERZ <- data_ERZ[which(data_ERZ$N_callable_sites > 5000),]
data_ERZ_mean <- mean(data_ERZ$pi_abs)
data_ERW <- read.delim("0fold_popstats/Rasolark/heterogametic_windows_10000_steps_10000_exon_dist_0fold/Rasolark_2021_heterogametic_windows_10000_steps_10000_exon_dist_0fold_data.txt", sep="\t", head=T)
data_ERW$pi_abs[which(is.na(data_ERW$pi_abs))] <- 0
data_ERW <- data_ERW[which(data_ERW$N_callable_sites > 5000),]
data_ERW_mean <- mean(data_ERW$pi_abs)

# Skylark
data_NSA <- read.delim("Neutral_popstats/Skylark/autosomal_windows_10000_steps_10000_exon_dist_20000/Skylark_2021_autosomal_windows_10000_steps_10000_exon_dist_20000_data.txt", sep="\t", head=T)
data_NSA$pi_abs[which(is.na(data_NSA$pi_abs))] <- 0
data_NSA <- data_NSA[which(data_NSA$N_callable_sites > 5000),]
data_NSA_mean <- mean(data_NSA$pi_abs)
data_NSZ <- read.delim("Neutral_popstats/Skylark/homogametic_windows_10000_steps_10000_exon_dist_20000/Skylark_2021_homogametic_windows_10000_steps_10000_exon_dist_20000_data.txt", sep="\t", head=T)
data_NSZ$pi_abs[which(is.na(data_NSZ$pi_abs))] <- 0
data_NSZ <- data_NSZ[which(data_NSZ$N_callable_sites > 5000),]
data_NSZ_mean <- mean(data_NSZ$pi_abs)
data_NSW <- read.delim("Neutral_popstats/Skylark/heterogametic_windows_10000_steps_10000_exon_dist_20000/Skylark_2021_heterogametic_windows_10000_steps_10000_exon_dist_20000_data.txt", sep="\t", head=T)
data_NSW$pi_abs[which(is.na(data_NSW$pi_abs))] <- 0
data_NSW <- data_NSW[which(data_NSW$N_callable_sites > 5000),]
data_NSW_mean <- mean(data_NSW$pi_abs)

data_ESA <- read.delim("0fold_popstats/Skylark/autosomal_windows_10000_steps_10000_exon_dist_0fold/Skylark_2021_autosomal_windows_10000_steps_10000_exon_dist_0fold_data.txt", sep="\t", head=T)
data_ESA$pi_abs[which(is.na(data_ESA$pi_abs))] <- 0
data_ESA <- data_ESA[which(data_ESA$N_callable_sites > 5000),]
data_ESA_mean <- mean(data_ESA$pi_abs)
data_ESZ <- read.delim("0fold_popstats/Skylark/homogametic_windows_10000_steps_10000_exon_dist_0fold/Skylark_2021_homogametic_windows_10000_steps_10000_exon_dist_0fold_data.txt", sep="\t", head=T)
data_ESZ$pi_abs[which(is.na(data_ESZ$pi_abs))] <- 0
data_ESZ <- data_ESZ[which(data_ESZ$N_callable_sites > 5000),]
data_ESZ_mean <- mean(data_ESZ$pi_abs)
data_ESW <- read.delim("0fold_popstats/Skylark/heterogametic_windows_10000_steps_10000_exon_dist_0fold/Skylark_2021_heterogametic_windows_10000_steps_10000_exon_dist_0fold_data.txt", sep="\t", head=T)
data_ESW$pi_abs[which(is.na(data_ESW$pi_abs))] <- 0
data_ESW <- data_ESW[which(data_ESW$N_callable_sites > 5000),]
data_ESW_mean <- mean(data_ESW$pi_abs)

data_gen_div <- rbind(data_NRA, data_NRZ, data_NRW, data_ERA, data_ERZ, data_ERW, data_NSA, data_NSZ, data_NSW, data_ESA, data_ESZ, data_ESW)
data_gen_div$Species <- factor(data_gen_div$Project, order=T, labels=c("Raso lark", "Skylark"), levels=c("Rasolark_2021", "Skylark_2021"))
data_gen_div$data_type <- factor(data_gen_div$data_type, order=T, labels=c("Autosomal", "Z", "W"), levels=c("autosomal", "homogametic", "heterogametic"))
data_gen_div$Data2 <- factor(data_gen_div$exon_distance, order=T, labels=c("Neutral", "Functional"), levels=c("20000", "0fold"))

gen_div_plot <- ggplot() +
  geom_violin(data=data_gen_div, aes(x=Data2, y=pi_abs, fill=data_type), color="black", width=1, position=position_dodge(1)) +
  geom_boxplot(data=data_gen_div, aes(x=Data2, y=pi_abs, fill=data_type), color="black", width=0.05, position=position_dodge(1)) +
  scale_fill_manual(name="Genomic region", values = c("Autosomal"="#E4EAF0", "Z"="#b30000", "W"="#fecc5c")) +
  scale_x_discrete(expand=c(0,0)) +
  facet_wrap(~Species, nrow=1, scales="free_y") +
  ylab("π") +
  theme_bw() +
  theme(legend.position = "right",
        legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=20, face="bold"), #change legend title font size
        legend.text = element_text(size=20), #change legend text font sizes
        strip.background = element_rect(color="black", fill="white", linewidth=1),
        strip.text = element_text(size=16, color="black"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=22),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=20, color="black"),
        axis.text.x = element_text(size=20, color="black"))

gen_div_plot

# Perform Bootstraps
Nboot <- 1000

gen_div <- as.data.frame(matrix(NA, Nboot*12, 4))
colnames(gen_div) <- c("pi_abs", "Data1", "Data2", "Species")

Ratios <- as.data.frame(matrix(NA, Nboot*38, 4))
colnames(Ratios) <- c("Ratio", "Data1", "Data2", "Species")

set.seed(123)
for(i in seq(1, Nboot*12, 12)) {
  
  # Record data
  gen_div$pi_abs[i] <- mean(sample(data_NRA$pi_abs, nrow(data_NRA), replace = TRUE))
  gen_div$Data1[i] <- "Autosomal"
  gen_div$Data2[i] <- "Neutral"
  gen_div$Species[i] <- "Raso lark"
  gen_div$pi_abs[i+1] <- mean(sample(data_NRZ$pi_abs, nrow(data_NRZ), replace = TRUE))
  gen_div$Data1[i+1] <- "Z"
  gen_div$Data2[i+1] <- "Neutral"
  gen_div$Species[i+1] <- "Raso lark"
  gen_div$pi_abs[i+2] <- mean(sample(data_NRW$pi_abs, nrow(data_NRW), replace = TRUE))
  gen_div$Data1[i+2] <- "W"
  gen_div$Data2[i+2] <- "Neutral"
  gen_div$Species[i+2] <- "Raso lark"
  gen_div$pi_abs[i+3] <- mean(sample(data_ERA$pi_abs, nrow(data_ERA), replace = TRUE))
  gen_div$Data1[i+3] <- "Autosomal"
  gen_div$Data2[i+3] <- "0-fold"
  gen_div$Species[i+3] <- "Raso lark"
  gen_div$pi_abs[i+4] <- mean(sample(data_ERZ$pi_abs, nrow(data_ERZ), replace = TRUE))
  gen_div$Data1[i+4] <- "Z"
  gen_div$Data2[i+4] <- "0-fold"
  gen_div$Species[i+4] <- "Raso lark"
  gen_div$pi_abs[i+5] <- mean(sample(data_ERW$pi_abs, nrow(data_ERW), replace = TRUE))
  gen_div$Data1[i+5] <- "W"
  gen_div$Data2[i+5] <- "0-fold"
  gen_div$Species[i+5] <- "Raso lark"
  gen_div$pi_abs[i+6] <- mean(sample(data_NSA$pi_abs, nrow(data_NSA), replace = TRUE))
  gen_div$Data1[i+6] <- "Autosomal"
  gen_div$Data2[i+6] <- "Neutral"
  gen_div$Species[i+6] <- "Skylark"
  gen_div$pi_abs[i+7] <- mean(sample(data_NSZ$pi_abs, nrow(data_NSZ), replace = TRUE))
  gen_div$Data1[i+7] <- "Z"
  gen_div$Data2[i+7] <- "Neutral"
  gen_div$Species[i+7] <- "Skylark"
  gen_div$pi_abs[i+8] <- mean(sample(data_NSW$pi_abs, nrow(data_NSW), replace = TRUE))
  gen_div$Data1[i+8] <- "W"
  gen_div$Data2[i+8] <- "Neutral"
  gen_div$Species[i+8] <- "Skylark"
  gen_div$pi_abs[i+9] <- mean(sample(data_ESA$pi_abs, nrow(data_ESA), replace = TRUE))
  gen_div$Data1[i+9] <- "Autosomal"
  gen_div$Data2[i+9] <- "0-fold"
  gen_div$Species[i+9] <- "Skylark"
  gen_div$pi_abs[i+10] <- mean(sample(data_ESZ$pi_abs, nrow(data_ESZ), replace = TRUE))
  gen_div$Data1[i+10] <- "Z"
  gen_div$Data2[i+10] <- "0-fold"
  gen_div$Species[i+10] <- "Skylark"
  gen_div$pi_abs[i+11] <- mean(sample(data_ESW$pi_abs, nrow(data_ESW), replace = TRUE))
  gen_div$Data1[i+11] <- "W"
  gen_div$Data2[i+11] <- "0-fold"
  gen_div$Species[i+11] <- "Skylark"
}

for(i in seq(1, Nboot*38, 38)) {
  
  # Resample
  resample_NRA <- sample(data_NRA$pi_abs, nrow(data_NRA), replace = TRUE)
  resample_NRZ <- sample(data_NRZ$pi_abs, nrow(data_NRZ), replace = TRUE)
  resample_NRW <- sample(data_NRW$pi_abs, nrow(data_NRW), replace = TRUE)
  resample_ERA <- sample(data_ERA$pi_abs, nrow(data_ERA), replace = TRUE)
  resample_ERZ <- sample(data_ERZ$pi_abs, nrow(data_ERZ), replace = TRUE)
  resample_ERW <- sample(data_ERW$pi_abs, nrow(data_ERW), replace = TRUE)
  resample_NSA <- sample(data_NSA$pi_abs, nrow(data_NSA), replace = TRUE)
  resample_NSZ <- sample(data_NSZ$pi_abs, nrow(data_NSZ), replace = TRUE)
  resample_NSW <- sample(data_NSW$pi_abs, nrow(data_NSW), replace = TRUE)
  resample_ESA <- sample(data_ESA$pi_abs, nrow(data_ESA), replace = TRUE)
  resample_ESZ <- sample(data_ESZ$pi_abs, nrow(data_ESZ), replace = TRUE)
  resample_ESW <- sample(data_ESW$pi_abs, nrow(data_ESW), replace = TRUE)

  # Record data
  Ratios$Ratio[i] <- mean(resample_NRZ)/mean(resample_NRA)
  Ratios$Data1[i] <- "Z:A"
  Ratios$Data2[i] <- "Neutral"
  Ratios$Species[i] <- "Raso lark"
  Ratios$Ratio[i+1] <- mean(resample_NRW)/mean(resample_NRA)
  Ratios$Data1[i+1] <- "W:A"
  Ratios$Data2[i+1] <- "Neutral"
  Ratios$Species[i+1] <- "Raso lark"
  Ratios$Ratio[i+2] <- mean(resample_NRW)/mean(resample_NRZ)
  Ratios$Data1[i+2] <- "W:Z"
  Ratios$Data2[i+2] <- "Neutral"
  Ratios$Species[i+2] <- "Raso lark"
  Ratios$Ratio[i+3] <- mean(resample_NSZ)/mean(resample_NSA)
  Ratios$Data1[i+3] <- "Z:A"
  Ratios$Data2[i+3] <- "Neutral"
  Ratios$Species[i+3] <- "Skylark"
  Ratios$Ratio[i+4] <- mean(resample_NSW)/mean(resample_NSA)
  Ratios$Data1[i+4] <- "W:A"
  Ratios$Data2[i+4] <- "Neutral"
  Ratios$Species[i+4] <- "Skylark"
  Ratios$Ratio[i+5] <- mean(resample_NSW)/mean(resample_NSZ)
  Ratios$Data1[i+5] <- "W:Z"
  Ratios$Data2[i+5] <- "Neutral"
  Ratios$Species[i+5] <- "Skylark"
  Ratios$Ratio[i+6] <- mean(resample_NRA)/mean(resample_NSA)
  Ratios$Data1[i+6] <- "RA:SA"
  Ratios$Data2[i+6] <- "Neutral"
  Ratios$Species[i+6] <- "Raso lark:Skylark"
  Ratios$Ratio[i+7] <- mean(resample_NRZ)/mean(resample_NSZ)
  Ratios$Data1[i+7] <- "RZ:SZ"
  Ratios$Data2[i+7] <- "Neutral"
  Ratios$Species[i+7] <- "Raso lark:Skylark"
  Ratios$Ratio[i+8] <- mean(resample_NRW)/mean(resample_NSW)
  Ratios$Data1[i+8] <- "RW:SW"
  Ratios$Data2[i+8] <- "Neutral"
  Ratios$Species[i+8] <- "Raso lark:Skylark"
  
  # Record data
  Ratios$Ratio[i+9] <- mean(resample_ERZ)/mean(resample_ERA)
  Ratios$Data1[i+9] <- "Z:A"
  Ratios$Data2[i+9] <- "0-fold"
  Ratios$Species[i+9] <- "Raso lark"
  Ratios$Ratio[i+10] <- mean(resample_ERW)/mean(resample_ERA)
  Ratios$Data1[i+10] <- "W:A"
  Ratios$Data2[i+10] <- "0-fold"
  Ratios$Species[i+10] <- "Raso lark"
  Ratios$Ratio[i+11] <- mean(resample_ERW)/mean(resample_ERZ)
  Ratios$Data1[i+11] <- "W:Z"
  Ratios$Data2[i+11] <- "0-fold"
  Ratios$Species[i+11] <- "Raso lark"
  Ratios$Ratio[i+12] <- mean(resample_ESZ)/mean(resample_ESA)
  Ratios$Data1[i+12] <- "Z:A"
  Ratios$Data2[i+12] <- "0-fold"
  Ratios$Species[i+12] <- "Skylark"
  Ratios$Ratio[i+13] <- mean(resample_ESW)/mean(resample_ESA)
  Ratios$Data1[i+13] <- "W:A"
  Ratios$Data2[i+13] <- "0-fold"
  Ratios$Species[i+13] <- "Skylark"
  Ratios$Ratio[i+14] <- mean(resample_ESW)/mean(resample_ESZ)
  Ratios$Data1[i+14] <- "W:Z"
  Ratios$Data2[i+14] <- "0-fold"
  Ratios$Species[i+14] <- "Skylark"
  Ratios$Ratio[i+15] <- mean(resample_ERA)/mean(resample_ESA)
  Ratios$Data1[i+15] <- "RA:SA"
  Ratios$Data2[i+15] <- "0-fold"
  Ratios$Species[i+15] <- "Raso lark:Skylark"
  Ratios$Ratio[i+16] <- mean(resample_ERZ)/mean(resample_ESZ)
  Ratios$Data1[i+16] <- "RZ:SZ"
  Ratios$Data2[i+16] <- "0-fold"
  Ratios$Species[i+16] <- "Raso lark:Skylark"
  Ratios$Ratio[i+17] <- mean(resample_ERW)/mean(resample_ESW)
  Ratios$Data1[i+17] <- "RW:SW"
  Ratios$Data2[i+17] <- "0-fold"
  Ratios$Species[i+17] <- "Raso lark:Skylark"
  
  
  # Record data
  Ratios$Ratio[i+18] <- mean(resample_ERA)/mean(resample_NRA)
  Ratios$Data1[i+18] <- "A:A"
  Ratios$Data2[i+18] <- "0-fold:Neutral"
  Ratios$Species[i+18] <- "Raso lark"
  Ratios$Ratio[i+19] <- mean(resample_ERZ)/mean(resample_NRZ)
  Ratios$Data1[i+19] <- "Z:Z"
  Ratios$Data2[i+19] <- "0-fold:Neutral"
  Ratios$Species[i+19] <- "Raso lark"
  Ratios$Ratio[i+20] <- mean(resample_ERW)/mean(resample_NRW)
  Ratios$Data1[i+20] <- "W:W"
  Ratios$Data2[i+20] <- "0-fold:Neutral"
  Ratios$Species[i+20] <- "Raso lark"
  Ratios$Ratio[i+21] <- mean(resample_ESA)/mean(resample_NSA)
  Ratios$Data1[i+21] <- "A:A"
  Ratios$Data2[i+21] <- "0-fold:Neutral"
  Ratios$Species[i+21] <- "Skylark"
  Ratios$Ratio[i+22] <- mean(resample_ESZ)/mean(resample_NSZ)
  Ratios$Data1[i+22] <- "Z:Z"
  Ratios$Data2[i+22] <- "0-fold:Neutral"
  Ratios$Species[i+22] <- "Skylark"
  Ratios$Ratio[i+23] <- mean(resample_ESW)/mean(resample_NSW)
  Ratios$Data1[i+23] <- "W:W"
  Ratios$Data2[i+23] <- "0-fold:Neutral"
  Ratios$Species[i+23] <- "Skylark"
  
  # Record data
  Ratios$Ratio[i+24] <- (mean(resample_ERA)/mean(resample_NRA))/(mean(resample_ESA)/mean(resample_NSA))
  Ratios$Data1[i+24] <- "RA:SA"
  Ratios$Data2[i+24] <- "0-fold:Neutral"
  Ratios$Species[i+24] <- "Raso lark:Skylark"
  Ratios$Ratio[i+25] <- (mean(resample_ERZ)/mean(resample_NRZ))/(mean(resample_ESZ)/mean(resample_NSZ))
  Ratios$Data1[i+25] <- "RZ:SZ"
  Ratios$Data2[i+25] <- "0-fold:Neutral"
  Ratios$Species[i+25] <- "Raso lark:Skylark"
  Ratios$Ratio[i+26] <- (mean(resample_ERW)/mean(resample_NRW))/(mean(resample_ESW)/mean(resample_NSW))
  Ratios$Data1[i+26] <- "RW:SW"
  Ratios$Data2[i+26] <- "0-fold:Neutral"
  Ratios$Species[i+26] <- "Raso lark:Skylark"
  
  # Record data
  Ratios$Ratio[i+27] <- (mean(resample_ERZ)/mean(resample_ERA))/(mean(resample_NRZ)/mean(resample_NRA))
  Ratios$Data1[i+27] <- "Z:A"
  Ratios$Data2[i+27] <- "0-fold:Neutral"
  Ratios$Species[i+27] <- "Raso lark"
  Ratios$Ratio[i+28] <- (mean(resample_ERW)/mean(resample_ERA))/(mean(resample_NRW)/mean(resample_NRA))
  Ratios$Data1[i+28] <- "W:A"
  Ratios$Data2[i+28] <- "0-fold:Neutral"
  Ratios$Species[i+28] <- "Raso lark"
  Ratios$Ratio[i+29] <- (mean(resample_ERW)/mean(resample_ERZ))/(mean(resample_NRW)/mean(resample_NRZ))
  Ratios$Data1[i+29] <- "W:Z"
  Ratios$Data2[i+29] <- "0-fold:Neutral"
  Ratios$Species[i+29] <- "Raso lark"
  Ratios$Ratio[i+30] <- (mean(resample_ESZ)/mean(resample_ESA))/(mean(resample_NRZ)/mean(resample_NSA))
  Ratios$Data1[i+30] <- "Z:A"
  Ratios$Data2[i+30] <- "0-fold:Neutral"
  Ratios$Species[i+30] <- "Skylark"
  Ratios$Ratio[i+31] <- (mean(resample_ESW)/mean(resample_ESA))/(mean(resample_NSW)/mean(resample_NSA))
  Ratios$Data1[i+31] <- "W:A"
  Ratios$Data2[i+31] <- "0-fold:Neutral"
  Ratios$Species[i+31] <- "Skylark"
  Ratios$Ratio[i+31] <- (mean(resample_ESW)/mean(resample_ESZ))/(mean(resample_NSW)/mean(resample_NSZ))
  Ratios$Data1[i+31] <- "W:Z"
  Ratios$Data2[i+31] <- "0-fold:Neutral"
  Ratios$Species[i+31] <- "Skylark"
  
  # Record data
  Ratios$Ratio[i+32] <- (mean(resample_ERZ)/mean(resample_ERA))/(mean(resample_ESZ)/mean(resample_ESA))
  Ratios$Data1[i+32] <- "Z:A"
  Ratios$Data2[i+32] <- "0-fold"
  Ratios$Species[i+32] <- "Raso lark:Skylark"
  Ratios$Ratio[i+33] <- (mean(resample_ERW)/mean(resample_ERA))/(mean(resample_ESW)/mean(resample_ESA))
  Ratios$Data1[i+33] <- "W:A"
  Ratios$Data2[i+33] <- "0-fold"
  Ratios$Species[i+33] <- "Raso lark:Skylark"
  Ratios$Ratio[i+34] <- (mean(resample_ERW)/mean(resample_ERZ))/(mean(resample_ESW)/mean(resample_ESZ))
  Ratios$Data1[i+34] <- "W:Z"
  Ratios$Data2[i+34] <- "0-fold"
  Ratios$Species[i+34] <- "Raso lark:Skylark"
  Ratios$Ratio[i+35] <- (mean(resample_NRZ)/mean(resample_NRA))/(mean(resample_NSZ)/mean(resample_NSA))
  Ratios$Data1[i+35] <- "Z:A"
  Ratios$Data2[i+35] <- "Neutral"
  Ratios$Species[i+35] <- "Raso lark:Skylark"
  Ratios$Ratio[i+36] <- (mean(resample_NRW)/mean(resample_NRA))/(mean(resample_NSW)/mean(resample_NSA))
  Ratios$Data1[i+36] <- "W:A"
  Ratios$Data2[i+36] <- "Neutral"
  Ratios$Species[i+36] <- "Raso lark:Skylark"
  Ratios$Ratio[i+37] <- (mean(resample_NRW)/mean(resample_NRZ))/(mean(resample_NSW)/mean(resample_NSZ))
  Ratios$Data1[i+37] <- "W:Z"
  Ratios$Data2[i+37] <- "Neutral"
  Ratios$Species[i+37] <- "Raso lark:Skylark"

}

gen_div$Species <- factor(gen_div$Species, order=T, levels=c("Raso lark", "Skylark"))
gen_div$Data1 <- factor(gen_div$Data1, order=T, levels=c("Autosomal", "Z", "W"))
gen_div$Data2 <- factor(gen_div$Data2, order=T, labels=c("Intergenic", "0-fold"), levels=c("Neutral", "0-fold"))

gen_div_plot2 <- ggplot() +
  geom_violin(data=gen_div, aes(x=Data2, y=pi_abs, fill=Data1), color="black", width=1, position=position_dodge(1)) +
  geom_boxplot(data=gen_div, aes(x=Data2, y=pi_abs, fill=Data1), color="black", width=0.05, position=position_dodge(1), show.legend = FALSE) +
  scale_fill_manual(name="Genomic region", values = c("Autosomal"="#E4EAF0", "Z"="#b30000", "W"="#fecc5c")) +
  scale_x_discrete(expand=c(0,0)) +
  facet_wrap(~Species, nrow=1, scales="free_y") +
  ylab("Bootstrapped π") +
  theme_bw() +
  theme(legend.position = "right",
        legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=20, face="bold"), #change legend title font size
        legend.text = element_text(size=20), #change legend text font sizes
        strip.background = element_rect(color="black", fill="white", linewidth=1),
        strip.text = element_text(size=16, color="black"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=22),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=15, color="black"),
        axis.text.x = element_text(size=20, color="black"))

gen_div_plot2

# Raso lark vs. Skylark differences
Ratios1 <- Ratios[which((Ratios$Species != "Raso lark" & Ratios$Species != "Skylark") & (Ratios$Data2 == "Neutral" | Ratios$Data2 == "0-fold") & (Ratios$Data1 == "RA:SA" | Ratios$Data1 == "RZ:SZ" | Ratios$Data1 == "RW:SW")),]
Ratios1$Data1 <- factor(Ratios1$Data1, order=T, levels=c("RA:SA", "RZ:SZ", "RW:SW"))
Ratios1$Data2 <- factor(Ratios1$Data2, order=T, labels=c("Intergenic", "0-fold"), levels=c("Neutral", "0-fold"))

data_NRA_mean/data_NSA_mean
data_NRZ_mean/data_NSZ_mean
data_NRW_mean/data_NSW_mean

data_ERA_mean/data_ESA_mean
data_ERZ_mean/data_ESZ_mean
data_ERW_mean/data_ESW_mean

# A
Ratios1A <- Ratios1[which(Ratios1$Data1 == "RA:SA"),]
Ratios1A_test <- rep(NA, Nboot)
for(i in 1:Nboot) {
  Ratios1A_test[i] <- Ratios1A$Ratio[i*2-1] - Ratios1A$Ratio[i*2]
}
diff <- data_NRA_mean/data_NSA_mean - data_ERA_mean/data_ESA_mean
diff
sd(Ratios1A_test)
quantile(Ratios1A_test, p=c(0.025, 0.975))
p_value_1A <- mean(abs(Ratios1A_test - mean(Ratios1A_test)) >= abs(diff))
p_value_1A

# Z
Ratios1B <- Ratios1[which(Ratios1$Data1 == "RZ:SZ"),]
Ratios1B_test <- rep(NA, Nboot)
for(i in 1:Nboot) {
  Ratios1B_test[i] <- Ratios1B$Ratio[i*2-1] - Ratios1B$Ratio[i*2]
}
diff <- data_NRZ_mean/data_NSZ_mean - data_ERZ_mean/data_ESZ_mean
diff
sd(Ratios1B_test)
quantile(Ratios1B_test, p=c(0.025, 0.975))
p_value_1B <- mean(abs(Ratios1B_test - mean(Ratios1B_test)) >= abs(diff))
p_value_1B

# W
Ratios1C <- Ratios1[which(Ratios1$Data1 == "RW:SW"),]
Ratios1C_test <- rep(NA, Nboot)
for(i in 1:Nboot) {
  Ratios1C_test[i] <- Ratios1C$Ratio[i*2-1] - Ratios1C$Ratio[i*2]
}
diff <- data_NRW_mean/data_NSW_mean - data_ERW_mean/data_ESW_mean
diff
sd(Ratios1C_test)
quantile(Ratios1C_test, p=c(0.025, 0.975))
p_value_1C <- mean(abs(Ratios1C_test - mean(Ratios1C_test)) >= abs(diff))
p_value_1C

ratio1 <- ggplot() +
  geom_violin(data=Ratios1, aes(x=Data2, y=Ratio, fill=Data1), color="black", width=1, position=position_dodge(1)) +
  geom_boxplot(data=Ratios1, aes(x=Data2, y=Ratio, fill=Data1), color="black", width=0.05, position=position_dodge(1)) +
  geom_hline(yintercept = 1, linewidth=1, linetype = 2, color = "#d7191c") +
  scale_fill_manual(name="Genomic region", values = c("RA:SA"="#E4EAF0", "RZ:SZ"="#b30000", "RW:SW"="#fecc5c"), labels=c("Autosomal", "Z", "W")) +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_continuous(limits=c(0, 0.5), expand=c(0,0)) +
  facet_wrap(~Data1, nrow=1, labeller = as_labeller(c("RA:SA"="Autosomal (Raso lark:Skylark)", "RZ:SZ"="Z (Raso lark:Skylark)", "RW:SW"="W (Raso lark:Skylark)"))) +
  xlab("Raso lark:Skylark") +
  ylab("Bootstrapped π ratio") +
  guides(fill="none") +
  theme_bw() +
  theme(legend.position = "right",
        legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=20, face="bold"), #change legend title font size
        legend.text = element_text(size=20), #change legend text font sizes
        strip.background = element_rect(color="black", fill="white", linewidth=1),
        strip.text = element_text(size=16, color="black"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=22),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=15, color="black"),
        axis.text.x = element_text(size=20, color="black"))


## Test deviation from expected ratios
Ratios2 <- Ratios[which((Ratios$Species == "Raso lark" | Ratios$Species == "Skylark") & Ratios$Data2 == "Neutral" & Ratios$Data1 != "W:Z"),]
intercepts <- data.frame(Data1 = c("Z:A", "W:A"), intercept = c(3/4, 1/4))
intercepts$Data1 <- factor(intercepts$Data1, order=T, levels=c("Z:A", "W:A"))
yscales <- data.frame(Data1 = c("Z:A", "W:A"), ymin = c(0, 0), ymax = c(1, 0.25/(0.75)))
Ratios2 <- Ratios2 |> left_join(yscales, by = "Data1")
Ratios2$Species <- factor(Ratios2$Species, order=T, levels=c("Raso lark", "Skylark"))
Ratios2$Data1 <- factor(Ratios2$Data1, order=T, levels=c("Z:A", "W:A"))

breaks_fun <- function(x) {
  if (max(x) > 0.5) {
    c(0.00, 0.25, 0.50, 0.75, 1.00)
  } else {
    c(0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30)
  }
}


ratio2 <- ggplot() +
  geom_violin(data=Ratios2, aes(x=Species, y=Ratio, fill=Species), color="black", width=1, position=position_dodge(1)) +
  geom_boxplot(data=Ratios2, aes(x=Species, y=Ratio, fill=Species), color="black", width=0.05, position=position_dodge(1)) +
  geom_hline(data=intercepts, aes(yintercept = intercept), linewidth=1, linetype = 2, color = "#d7191c") +
  scale_fill_manual(values = c("Raso lark" = "#92c5de", "Skylark" = "#f4a582")) +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_continuous(breaks = breaks_fun, limits = c(0, NA), expand = c(0,0)) + 
  geom_blank(data = Ratios2, aes(y = ymin)) +
  geom_blank(data = Ratios2, aes(y = ymax)) + 
  facet_wrap(~Data1, scales = "free_y", nrow=1, labeller = as_labeller(c("Z:A"="Intergenic (Z:Autosomal)", "W:A"="Intergenic (W:Autosomal)"))) +
  guides(fill="none", color="none") +
  ylab("Bootstrapped π ratio") +
  theme_bw() +
  theme(legend.position = "right",
        legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=20, face="bold"), #change legend title font size
        legend.text = element_text(size=20), #change legend text font sizestrip.background = element_rect(color="black", fill="white", linewidth=1),
        strip.background = element_rect(color="black", fill="white", linewidth=1),
        strip.text = element_text(size=16, color="black"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=22),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=15, color="black"),
        axis.text.x = element_text(size=20, color="black"))



## Test deviation from expected ratios
Ratios2B <- Ratios2[which(Ratios2$Data1 == "Z:A"),]
Ratios2BR <- Ratios2B[which(Ratios2B$Species == "Raso lark"),]
Ratios2BS <- Ratios2B[which(Ratios2B$Species == "Skylark"),]

# Raso Z:A
mean1 <- data_NRZ_mean / data_NRA_mean
mean1
mean2 <- mean1 - 0.75
mean2
quantile(Ratios2BR$Ratio, p=c(0.025, 0.975)) - 0.75
boot_null <- Ratios2BR$Ratio - mean(Ratios2BR$Ratio)
p_value_R2 <- mean(abs(boot_null) >= abs(mean2))
p_value_R2
p_value_R2*3

# Sky Z:A
mean1 <- data_NSZ_mean / data_NSA_mean
mean1
mean2 <- mean1 - 0.75
mean2
quantile(Ratios2BS$Ratio, p=c(0.025, 0.975)) - 0.75
boot_null <- Ratios2BS$Ratio - mean(Ratios2BS$Ratio)
p_value_S2 <- mean(abs(boot_null) >= abs(mean2))
p_value_S2
p_value_S2*3


# Raso Z:A vs. Sky Z:A
Ratios2BRS_test <- rep(NA, Nboot)
for(i in 1:Nboot) {
  Ratios2BRS_test[i] <- Ratios2B$Ratio[i*2-1] - Ratios2B$Ratio[i*2]
}
diff <- data_NRZ_mean/data_NRA_mean - data_NSZ_mean/data_NSA_mean
diff
sd(Ratios2BRS_test)
quantile(Ratios2BRS_test, p=c(0.025, 0.975))
p_value_RS2 <- mean(abs(Ratios2BRS_test - mean(Ratios2BRS_test)) >= abs(diff))
p_value_RS2
p_value_RS2 * 3

# Raso W:A and Sky W:A
data_NRW_mean/data_NRA_mean
data_NSW_mean/data_NSA_mean



### SFS
Rasolark_A_sfs <- read.table("DFE/SFS/Rasolark_A_sfs", skip=1)
Rasolark_A_sfs[1,] <- Rasolark_A_sfs[1,]/Rasolark_A_sfs[1,ncol(Rasolark_A_sfs)]
Rasolark_A_sfs[2,] <- Rasolark_A_sfs[2,]/Rasolark_A_sfs[2,ncol(Rasolark_A_sfs)]
Rasolark_A_sfs <- as.data.frame(t(Rasolark_A_sfs[,-ncol(Rasolark_A_sfs)]))
Rasolark_A_sfs$bin <- 1:nrow(Rasolark_A_sfs)
colnames(Rasolark_A_sfs) <- c("4-fold", "0-fold", "bin")
Rasolark_A_sfs <- Rasolark_A_sfs |> pivot_longer(
  cols=c("4-fold", "0-fold"), names_to = "Data", values_to = "value")
Rasolark_A_sfs$Species <- "Raso lark"
Rasolark_A_sfs$Region <- "Autosomal"

Rasolark_Z_sfs <- read.table("DFE/SFS/Rasolark_Z_sfs", skip=1)
Rasolark_Z_sfs[1,] <- Rasolark_Z_sfs[1,]/Rasolark_Z_sfs[1,ncol(Rasolark_Z_sfs)]
Rasolark_Z_sfs[2,] <- Rasolark_Z_sfs[2,]/Rasolark_Z_sfs[2,ncol(Rasolark_Z_sfs)]
Rasolark_Z_sfs <- as.data.frame(t(Rasolark_Z_sfs[,-ncol(Rasolark_Z_sfs)]))
Rasolark_Z_sfs$bin <- 1:nrow(Rasolark_Z_sfs)
colnames(Rasolark_Z_sfs) <- c("4-fold", "0-fold", "bin")
Rasolark_Z_sfs <- Rasolark_Z_sfs |> pivot_longer(
  cols=c("4-fold", "0-fold"), names_to = "Data", values_to = "value")
Rasolark_Z_sfs$Species <- "Raso lark"
Rasolark_Z_sfs$Region <- "Z"

Rasolark_W_sfs <- read.table("DFE/SFS/Rasolark_W_sfs", skip=1)
Rasolark_W_sfs[1,] <- Rasolark_W_sfs[1,]/Rasolark_W_sfs[1,ncol(Rasolark_W_sfs)]
Rasolark_W_sfs[2,] <- Rasolark_W_sfs[2,]/Rasolark_W_sfs[2,ncol(Rasolark_W_sfs)]
Rasolark_W_sfs <- as.data.frame(t(Rasolark_W_sfs[,-ncol(Rasolark_W_sfs)]))
Rasolark_W_sfs$bin <- 1:nrow(Rasolark_W_sfs)
colnames(Rasolark_W_sfs) <- c("4-fold", "0-fold", "bin")
Rasolark_W_sfs <- Rasolark_W_sfs |> pivot_longer(
  cols=c("4-fold", "0-fold"), names_to = "Data", values_to = "value")
Rasolark_W_sfs$Species <- "Raso lark"
Rasolark_W_sfs$Region <- "W"

Skylark_A_sfs <- read.table("DFE/SFS/Skylark_A_sfs", skip=1)
Skylark_A_sfs[1,] <- Skylark_A_sfs[1,]/Skylark_A_sfs[1,ncol(Skylark_A_sfs)]
Skylark_A_sfs[2,] <- Skylark_A_sfs[2,]/Skylark_A_sfs[2,ncol(Skylark_A_sfs)]
Skylark_A_sfs <- as.data.frame(t(Skylark_A_sfs[,-ncol(Skylark_A_sfs)]))
Skylark_A_sfs$bin <- 1:nrow(Skylark_A_sfs)
colnames(Skylark_A_sfs) <- c("4-fold", "0-fold", "bin")
Skylark_A_sfs <- Skylark_A_sfs |> pivot_longer(
  cols=c("4-fold", "0-fold"), names_to = "Data", values_to = "value")
Skylark_A_sfs$Species <- "Skylark"
Skylark_A_sfs$Region <- "Autosomal"

Skylark_Z_sfs <- read.table("DFE/SFS/Skylark_Z_sfs", skip=1)
Skylark_Z_sfs[1,] <- Skylark_Z_sfs[1,]/Skylark_Z_sfs[1,ncol(Skylark_Z_sfs)]
Skylark_Z_sfs[2,] <- Skylark_Z_sfs[2,]/Skylark_Z_sfs[2,ncol(Skylark_Z_sfs)]
Skylark_Z_sfs <- as.data.frame(t(Skylark_Z_sfs[,-ncol(Skylark_Z_sfs)]))
Skylark_Z_sfs$bin <- 1:nrow(Skylark_Z_sfs)
colnames(Skylark_Z_sfs) <- c("4-fold", "0-fold", "bin")
Skylark_Z_sfs <- Skylark_Z_sfs |> pivot_longer(
  cols=c("4-fold", "0-fold"), names_to = "Data", values_to = "value")
Skylark_Z_sfs$Species <- "Skylark"
Skylark_Z_sfs$Region <- "Z"

Skylark_W_sfs <- read.table("DFE/SFS/Skylark_W_sfs", skip=1)
Skylark_W_sfs[1,] <- Skylark_W_sfs[1,]/Skylark_W_sfs[1,ncol(Skylark_W_sfs)]
Skylark_W_sfs[2,] <- Skylark_W_sfs[2,]/Skylark_W_sfs[2,ncol(Skylark_W_sfs)]
Skylark_W_sfs <- as.data.frame(t(Skylark_W_sfs[,-ncol(Skylark_W_sfs)]))
Skylark_W_sfs$bin <- 1:nrow(Skylark_W_sfs)
colnames(Skylark_W_sfs) <- c("4-fold", "0-fold", "bin")
Skylark_W_sfs <- Skylark_W_sfs |> pivot_longer(
  cols=c("4-fold", "0-fold"), names_to = "Data", values_to = "value")
Skylark_W_sfs$Species <- "Skylark"
Skylark_W_sfs$Region <- "W"

data_SFS <- rbind(Rasolark_A_sfs, Rasolark_Z_sfs, Rasolark_W_sfs, Skylark_A_sfs,Skylark_Z_sfs, Skylark_W_sfs)
data_SFS$Data <- factor(data_SFS$Data, order=T, c("4-fold", "0-fold"))
data_SFS$Species <- factor(data_SFS$Species, order=T, c("Raso lark", "Skylark"))
data_SFS$Region <- factor(data_SFS$Region, order=T, c("Autosomal", "Z", "W"))


vline <- as.data.frame(c(36, 36, 27, 26, 9, 10))
colnames(vline) <- "max_haps"
vline$Species <- c("Raso lark", "Skylark", "Raso lark", "Skylark", "Raso lark", "Skylark")
vline$Region <- c("Autosomal", "Autosomal", "Z", "Z", "W", "W")
vline$Species <- factor(vline$Species, order=T, c("Raso lark", "Skylark"))
vline$Region <- factor(vline$Region, order=T, c("Autosomal", "Z", "W"))

sfs_plot <- ggplot() +
  geom_bar(data=data_SFS, aes(x=bin, y=value, fill=Data), color=scales::alpha("black", 0.5),  position_dodge2(preserve = "single"), stat="identity") +
  geom_vline(data=vline, aes(xintercept = max_haps), linetype="dashed", color="black", size=0.5) +
  scale_fill_manual(name="Site class", values = c("4-fold"="#440154FF", "0-fold"="#FDE725FF")) +
  facet_wrap(Region~Species, nrow=1, scales="free_y") +
  scale_x_continuous(limits=c(0,35), breaks=c(1, 18, 35)) +
  labs(x="Site-Frequency Spectrum", y="Frequency of segregating sites") +
  theme_bw() +
  theme(legend.position = "right",
        legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=20, face="bold"), #change legend title font size
        legend.text = element_text(size=20), #change legend text font size
        strip.background = element_rect(color="black", fill="white", linewidth=1),
        strip.text = element_text(size=16, color="black"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=22),
        axis.title.x = element_text(size=22),
       # axis.text.y = element_text(size=15, color="black"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=20, color="black"),
        axis.ticks.y = element_blank())
sfs_plot

### DFE

# Source functions from polyDFE postprocessing companion script
source("DFE/postprocessing.R")

# import data
polyDFEfiles <- dir("DFE/DFE_runs/Model_selection")
datasets <- unique(apply(str_split(polyDFEfiles, "_", simplify = TRUE)[, 1:2], 1, paste, collapse = "_"))
datasets
bins <- c(-10, -1, 0, 1, 10)

###
i="Rasolark_A"
est = list()
models <- polyDFEfiles[grep(i, polyDFEfiles)]
for(j in models) {
  est = c(est, parseOutput(paste("DFE/DFE_runs/Model_selection/", j, sep="")))
}
aic <- as.data.frame(getAICweights(est))
aic$analysis <-
  sapply(est, function(x) {
    x$values <- x$values[[1]]
    getModelName(x)
  })
aic$eps_ann <- 
  sapply(est, function(x) {
    x$values[[1]]["eps_an"]
  })
aic$alpha <- 
  sapply(est, function(x) {
    x$alpha[[1]]["alpha_dfe"]
  })
aic$grad <- 
  sapply(est, function(x) {
    x$criteria
  })
aic$analysis[which(aic$grad > 0.01)]
aic
best <- 10
aic$pval <- compareModels(est[rep(best, 12)], est[1:12], nested = TRUE)$LRT[,"p-value"]
alpha_dfe <- sum(sapply(1:length(est), function(i) aic[i, "weight"] * est[[i]]$alpha[[1]][["alpha_dfe"]]))
alpha_dfe
# 0.1158821513240845951787
select(aic, - analysis)
dfe_RA <- as.data.frame(t(t(sapply(est, getDiscretizedDFE, sRanges = bins))[best,]))
rownames(dfe_RA) <- aic$analysis[best]
colnames(dfe_RA) = toNames(bins)
dfe_RA$Species <- "Raso lark"
dfe_RA$Region <- "Autosomal"

###
i="Rasolark_Z"
est = list()
models <- polyDFEfiles[grep(i, polyDFEfiles)]
for(j in models) {
  est = c(est, parseOutput(paste("DFE/DFE_runs/Model_selection/", j, sep="")))
}
aic <- as.data.frame(getAICweights(est))
aic$analysis <-
  sapply(est, function(x) {
    x$values <- x$values[[1]]
    getModelName(x)
  })
aic$eps_ann <- 
  sapply(est, function(x) {
    x$values[[1]]["eps_an"]
  })
aic$alpha <- 
  sapply(est, function(x) {
    x$alpha[[1]]["alpha_dfe"]
  })
aic$grad <- 
  sapply(est, function(x) {
    x$criteria
  })
aic$analysis[which(aic$grad > 0.01)]
aic
best <- 10
aic$pval <- compareModels(est[rep(best, 12)], est[1:12], nested = TRUE)$LRT[,"p-value"]
alpha_dfe <- sum(sapply(1:length(est), function(i) aic[i, "weight"] * est[[i]]$alpha[[1]][["alpha_dfe"]]))
alpha_dfe
# 0.1306342458664643546573
select(aic, - analysis)
dfe_RZ <- as.data.frame(t(t(sapply(est, getDiscretizedDFE, sRanges = bins))[best,]))
rownames(dfe_RZ) <- aic$analysis[best]
colnames(dfe_RZ) = toNames(bins)
dfe_RZ$Species <- "Raso lark"
dfe_RZ$Region <- "Z"


###
i="Rasolark_W"
est = list()
models <- polyDFEfiles[grep(i, polyDFEfiles)]
for(j in models) {
  est = c(est, parseOutput(paste("DFE/DFE_runs/Model_selection/", j, sep="")))
}
aic <- as.data.frame(getAICweights(est))
aic$analysis <-
  sapply(est, function(x) {
    x$values <- x$values[[1]]
    getModelName(x)
  })
aic$eps_ann <- 
  sapply(est, function(x) {
    x$values[[1]]["eps_an"]
  })
aic$alpha <- 
  sapply(est, function(x) {
    x$alpha[[1]]["alpha_dfe"]
  })
aic$grad <- 
  sapply(est, function(x) {
    x$criteria
  })
aic$analysis[which(aic$grad > 0.01)]
aic
best <- 10
aic$pval <- compareModels(est[rep(best, 12)], est[1:12], nested = TRUE)$LRT[,"p-value"]
alpha_dfe <- sum(sapply(1:length(est), function(i) aic[i, "weight"] * est[[i]]$alpha[[1]][["alpha_dfe"]]))
alpha_dfe
# 0.1346945225147975488955
select(aic, - analysis)
dfe_RW <- as.data.frame(t(t(sapply(est, getDiscretizedDFE, sRanges = bins))[best,]))
rownames(dfe_RW) <- aic$analysis[best]
colnames(dfe_RW) = toNames(bins)
dfe_RW$Species <- "Raso lark"
dfe_RW$Region <- "W"


i="Skylark_A"
est = list()
models <- polyDFEfiles[grep(i, polyDFEfiles)]
for(j in models) {
  est = c(est, parseOutput(paste("DFE/DFE_runs/Model_selection/", j, sep="")))
}
aic <- as.data.frame(getAICweights(est))
aic$analysis <-
  sapply(est, function(x) {
    x$values <- x$values[[1]]
    getModelName(x)
  })
aic$eps_ann <- 
  sapply(est, function(x) {
    x$values[[1]]["eps_an"]
  })
aic$alpha <- 
  sapply(est, function(x) {
    x$alpha[[1]]["alpha_dfe"]
  })
aic$grad <- 
  sapply(est, function(x) {
    x$criteria
  })
aic$analysis[which(aic$grad > 0.01)]
aic
best <- 6
aic$pval <- compareModels(est[rep(best, 12)], est[1:12], nested = TRUE)$LRT[,"p-value"]
alpha_dfe <- sum(sapply(1:length(est), function(i) aic[i, "weight"] * est[[i]]$alpha[[1]][["alpha_dfe"]]))
alpha_dfe
#0.9924079988109463235091
select(aic, - analysis)
dfe_SA <- as.data.frame(t(t(sapply(est, getDiscretizedDFE, sRanges = bins))[best,]))
rownames(dfe_SA) <- aic$analysis[best]
colnames(dfe_SA) = toNames(bins)
dfe_SA$Species <- "Skylark"
dfe_SA$Region <- "Autosomal"


###
i="Skylark_Z"
est = list()
models <- polyDFEfiles[grep(i, polyDFEfiles)]
for(j in models) {
  est = c(est, parseOutput(paste("DFE/DFE_runs/Model_selection/", j, sep="")))
}
aic <- as.data.frame(getAICweights(est))
aic$analysis <-
  sapply(est, function(x) {
    x$values <- x$values[[1]]
    getModelName(x)
  })
aic$eps_ann <- 
  sapply(est, function(x) {
    x$values[[1]]["eps_an"]
  })
aic$alpha <- 
  sapply(est, function(x) {
    x$alpha[[1]]["alpha_dfe"]
  })
aic$grad <- 
  sapply(est, function(x) {
    x$criteria
  })
aic$analysis[which(aic$grad > 0.01)]
aic
best <- 5
aic$pval <- compareModels(est[rep(best, 12)], est[1:12], nested = TRUE)$LRT[,"p-value"]
alpha_dfe <- sum(sapply(1:length(est), function(i) aic[i, "weight"] * est[[i]]$alpha[[1]][["alpha_dfe"]]))
alpha_dfe
# 0.7143409575658318244606
select(aic, - analysis)
dfe_SZ <- as.data.frame(t(t(sapply(est, getDiscretizedDFE, sRanges = bins))[best,]))
rownames(dfe_SZ) <- aic$analysis[best]
colnames(dfe_SZ) = toNames(bins)
dfe_SZ$Species <- "Skylark"
dfe_SZ$Region <- "Z"

###
i="Skylark_W"
est = list()
models <- polyDFEfiles[grep(i, polyDFEfiles)]
for(j in models) {
  est = c(est, parseOutput(paste("DFE/DFE_runs/Model_selection/", j, sep="")))
}
aic <- as.data.frame(getAICweights(est))
aic$analysis <-
  sapply(est, function(x) {
    x$values <- x$values[[1]]
    getModelName(x)
  })
aic$eps_ann <- 
  sapply(est, function(x) {
    x$values[[1]]["eps_an"]
  })
aic$alpha <- 
  sapply(est, function(x) {
    x$alpha[[1]]["alpha_dfe"]
  })
aic$grad <- 
  sapply(est, function(x) {
    x$criteria
  })
aic$analysis[which(aic$grad > 0.01)]
aic
best <- 10
aic$pval <- compareModels(est[rep(best, 12)], est[1:12], nested = TRUE)$LRT[,"p-value"]
alpha_dfe <- sum(sapply(1:length(est), function(i) aic[i, "weight"] * est[[i]]$alpha[[1]][["alpha_dfe"]]))
alpha_dfe
# 0.4611299762934232271405
select(aic, - analysis)
dfe_SW <- as.data.frame(t(t(sapply(est, getDiscretizedDFE, sRanges = bins))[best,]))
rownames(dfe_SW) <- aic$analysis[best]
colnames(dfe_SW) = toNames(bins)
dfe_SW$Species <- "Skylark"
dfe_SW$Region <- "W"



########## Bootstraps from best models
i="Rasolark_A"
est = list()
polyDFEfiles <- dir(paste("DFE/DFE_runs/bootstraps/", i, "_DFE/", sep=""))
for(j in polyDFEfiles) {
  est = c(est, parseOutput(paste("DFE/DFE_runs/bootstraps/", i, "_DFE/", j, sep="")))
}
DFE <- t(sapply(est, getDiscretizedDFE, sRanges = bins))
colnames(DFE) = toNames(bins)
DFE_boot_RA <- data.frame(
  bin   = colnames(DFE),
  mean  = colMeans(DFE),
  low95 = apply(DFE, 2, quantile, probs = 0.025),
  upp95 = apply(DFE, 2, quantile, probs = 0.975))
DFE_boot_RA$Species <- "Raso lark"
DFE_boot_RA$Region <- "Autosomal"

i="Rasolark_Z"
est = list()
polyDFEfiles <- dir(paste("DFE/DFE_runs/bootstraps/", i, "_DFE/", sep=""))
for(j in polyDFEfiles) {
  est = c(est, parseOutput(paste("DFE/DFE_runs/bootstraps/", i, "_DFE/", j, sep="")))
}
DFE <- t(sapply(est, getDiscretizedDFE, sRanges = bins))
colnames(DFE) = toNames(bins)
DFE_boot_RZ <- data.frame(
  bin   = colnames(DFE),
  mean  = colMeans(DFE),
  low95 = apply(DFE, 2, quantile, probs = 0.025),
  upp95 = apply(DFE, 2, quantile, probs = 0.975))
DFE_boot_RZ$Species <- "Raso lark"
DFE_boot_RZ$Region <- "Z"

i="Rasolark_W"
est = list()
polyDFEfiles <- dir(paste("DFE/DFE_runs/bootstraps/", i, "_DFE/", sep=""))
for(j in polyDFEfiles) {
  est = c(est, parseOutput(paste("DFE/DFE_runs/bootstraps/", i, "_DFE/", j, sep="")))
}
DFE <- t(sapply(est, getDiscretizedDFE, sRanges = bins))
colnames(DFE) = toNames(bins)
DFE_boot_RW <- data.frame(
  bin   = colnames(DFE),
  mean  = colMeans(DFE),
  low95 = apply(DFE, 2, quantile, probs = 0.025),
  upp95 = apply(DFE, 2, quantile, probs = 0.975))
DFE_boot_RW$Species <- "Raso lark"
DFE_boot_RW$Region <- "W"

i="Skylark_A"
est = list()
polyDFEfiles <- dir(paste("DFE/DFE_runs/bootstraps/", i, "_DFE/", sep=""))
for(j in polyDFEfiles) {
  est = c(est, parseOutput(paste("DFE/DFE_runs/bootstraps/", i, "_DFE/", j, sep="")))
}
DFE <- t(sapply(est, getDiscretizedDFE, sRanges = bins))
colnames(DFE) = toNames(bins)
DFE_boot_SA <- data.frame(
  bin   = colnames(DFE),
  mean  = colMeans(DFE),
  low95 = apply(DFE, 2, quantile, probs = 0.025),
  upp95 = apply(DFE, 2, quantile, probs = 0.975))
DFE_boot_SA$Species <- "Skylark"
DFE_boot_SA$Region <- "Autosomal"

i="Skylark_Z"
est = list()
polyDFEfiles <- dir(paste("DFE/DFE_runs/bootstraps/", i, "_DFE/", sep=""))
for(j in polyDFEfiles) {
  est = c(est, parseOutput(paste("DFE/DFE_runs/bootstraps/", i, "_DFE/", j, sep="")))
}
DFE <- t(sapply(est, getDiscretizedDFE, sRanges = bins))
colnames(DFE) = toNames(bins)
DFE_boot_SZ <- data.frame(
  bin   = colnames(DFE),
  mean  = colMeans(DFE),
  low95 = apply(DFE, 2, quantile, probs = 0.025),
  upp95 = apply(DFE, 2, quantile, probs = 0.975))
DFE_boot_SZ$Species <- "Skylark"
DFE_boot_SZ$Region <- "Z"

i="Skylark_W"
est = list()
polyDFEfiles <- dir(paste("DFE/DFE_runs/bootstraps/", i, "_DFE/", sep=""))
for(j in polyDFEfiles) {
  est = c(est, parseOutput(paste("DFE/DFE_runs/bootstraps/", i, "_DFE/", j, sep="")))
}
DFE <- t(sapply(est, getDiscretizedDFE, sRanges = bins))
colnames(DFE) = toNames(bins)
DFE_boot_SW <- data.frame(
  bin   = colnames(DFE),
  mean  = colMeans(DFE),
  low95 = apply(DFE, 2, quantile, probs = 0.025),
  upp95 = apply(DFE, 2, quantile, probs = 0.975))
DFE_boot_SW$Species <- "Skylark"
DFE_boot_SW$Region <- "W"


# plot
DFE_data <- rbind(dfe_RA, dfe_RZ, dfe_RW, dfe_SA, dfe_SZ, dfe_SW)
DFE_data$Model <- rownames(DFE_data)
DFE_data <- DFE_data |> pivot_longer(
  cols=-c("Species", "Region", "Model"), names_to = "bin", values_to = "value")
DFE_data$Species <- factor(DFE_data$Species, order=T, levels=c("Raso lark", "Skylark"))
DFE_data$Region <- factor(DFE_data$Region, order=T, levels=c("Autosomal", "Z", "W"))
DFE_data$bin <- factor(DFE_data$bin, order=T, levels=unique(DFE_data$bin))
DFE_CI <- rbind(DFE_boot_RA, DFE_boot_RZ, DFE_boot_RW, DFE_boot_SA, DFE_boot_SZ, DFE_boot_SW)
DFE_CI$Species <- factor(DFE_CI$Species, order=T, levels=c("Raso lark", "Skylark"))
DFE_CI$Region <- factor(DFE_CI$Region, order=T, levels=c("Autosomal", "Z", "W"))
DFE_CI$bin <- factor(DFE_CI$bin, order=T, levels=unique(DFE_CI$bin))

DFE_plot <- ggplot() +
  geom_bar(data=DFE_data, aes(x=bin, y=value, fill=Species), color="black", position="dodge", width=1, stat="identity") +
  geom_point(data=DFE_CI, aes(x=bin, y=mean, group=Species), position=position_dodge(width=1), size=3) +
  geom_errorbar(data=DFE_CI, aes(x=bin, ymin=low95, ymax=upp95, group=Species), position="dodge", width=1, size=0.5) +
  scale_fill_manual(values = c("Raso lark" = "#92c5de", "Skylark" = "#f4a582")) +
  facet_wrap(~Region, nrow=1) +
  scale_y_continuous(limits=c(0,1), expand=c(0,0)) +
  labs(x=expression("Selection coefficient (" * S == 4 * N[e] * s * ")"), y="Distribution of fitness effects") +
  coord_cartesian(ylim = c(0, 1.02)) +
  theme_bw() +
  theme(legend.position = "right",
        legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=20, face="bold"), #change legend title font size
        legend.text = element_text(size=20), #change legend text font size
        strip.background = element_rect(color="black", fill="white", linewidth=1),
        strip.text = element_text(size=16, color="black"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=22),
        axis.title.x = element_text(size=22),
        axis.text.y = element_text(size=15, color="black"),
        axis.text.x = element_text(size=20, color="black"))

DFE_plot


### DOS ############################################
### Import evolutionary data
data <- read.delim("../../Skylark_2021/Results/Genes/Skylark_2021_Rasolark_2021_organised_data2.tsv", sep="\t", head=T)
options(scipen=999)

### Filter evolutionary data
data_evo <- data[which(data$Filter1=="OK" & data$Filter3=="OK" & data$Filter4=="OK" & data$Filter5=="OK"),]
data_evo$Strata2 <- data_evo$Strata
data_evo$Strata[which(data_evo$Strata == "PAR3" | data_evo$Strata == "PAR5")] <- "PAR"
data_evo$Strata <- factor(data_evo$Strata, order=T, labels=rev(c("S0", "S1", "S2", "S3", "4A", "3-a", "3-b", "5", "3-c", "PAR", "Autosomal")), levels=rev(c("S0", "S1", "S2", "S3", "4A", "3a", "3b", "5", "3c", "PAR", "autosomal")))
data_evo$Species <- factor(data_evo$Species, order=T, labels=rev(c("Raso lark", "Skylark")), levels=rev(c("Rasolark", "Skylark")))

data_evoA <- data_evo[which(data_evo$Region == "autosomal"),]
data_evoA <- data_evoA[which(data_evoA$SpeciesAlphaA > 0.001),]
data_evoA$SpeciesDNDS_DOS <- data_evoA$SpeciesBetaA / (data_evoA$SpeciesBetaA + data_evoA$SpeciesAlphaA)
data_evoA$SpeciesPNPS_DOS <- data_evoA$PNA / (data_evoA$PNA + data_evoA$PSA)
data_evoA$Data_type <- "Autosomal"
data_evoZ <- data_evo[which(data_evo$Region != "autosomal"),]
data_evoZ <- data_evoZ[which(data_evoZ$SpeciesAlphaZ > 0.001),]
data_evoZ$SpeciesDNDS_DOS <- data_evoZ$SpeciesBetaZ / (data_evoZ$SpeciesBetaZ + data_evoZ$SpeciesAlphaZ)
data_evoZ$SpeciesPNPS_DOS <- data_evoZ$PNZ / (data_evoZ$PNZ + data_evoZ$PSZ)
data_evoZ$Data_type <- "Z"
data_evoW <- data_evo[which(data_evo$Region != "autosomal"),]
data_evoW <- data_evoW[which(data_evoW$SpeciesAlphaW > 0.001),]
data_evoW$SpeciesDNDS_DOS <- data_evoW$SpeciesBetaW / (data_evoW$SpeciesBetaW + data_evoW$SpeciesAlphaW)
data_evoW$SpeciesPNPS_DOS <- data_evoW$PNW / (data_evoW$PNW + data_evoW$PSW)
data_evoW$Data_type <- "W"

data_evo2 <- rbind(data_evoA, data_evoZ, data_evoW)
data_evo2$Data_type <- factor(data_evo2$Data_type, order=T, levels=c("Autosomal", "Z", "W"))
data_evo2$SpeciesDOS <- data_evo2$SpeciesDNDS_DOS - data_evo2$SpeciesPNPS_DOS

DOS_RA <- data_evo2$SpeciesDOS[which(data_evo2$Species =="Raso lark" & data_evo2$Data_type == "Autosomal")]
DOS_RZ <- data_evo2$SpeciesDOS[which(data_evo2$Species =="Raso lark" & data_evo2$Data_type == "Z")]
DOS_RW <- data_evo2$SpeciesDOS[which(data_evo2$Species =="Raso lark" & data_evo2$Data_type == "W")]
DOS_SA <- data_evo2$SpeciesDOS[which(data_evo2$Species =="Skylark" & data_evo2$Data_type == "Autosomal")]
DOS_SZ <- data_evo2$SpeciesDOS[which(data_evo2$Species =="Skylark" & data_evo2$Data_type == "Z")]
DOS_SW <- data_evo2$SpeciesDOS[which(data_evo2$Species =="Skylark" & data_evo2$Data_type == "W")]

mean(DOS_RA)
mean(DOS_RZ)
mean(DOS_RW)

mean(DOS_SA)
mean(DOS_SZ)
mean(DOS_SW)


DOS_plot <- ggplot() +
  geom_histogram(data=data_evo2, aes(x=SpeciesDOS, group=interaction(Species, Data_type), fill=Species), color="black", bins=24, center = 0, position="identity", alpha=0.75) +
  scale_fill_manual(values = c("Raso lark" = "#92c5de", "Skylark" = "#f4a582")) +
  geom_vline(xintercept = 0, linewidth=1, linetype = 2, color = "black") +
  facet_wrap(~Data_type, nrow=1, scales = "free_y") +
  scale_x_continuous(limit = c(-1, 1)) +
  scale_y_continuous(expand = c(0.0, 0.0)) +
  labs(x="Direction of selection", y="Number of genes") +
  guides(fill="none", color="none") +
  theme_bw() +
  theme(legend.position = "right",
        legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=20, face="bold"), #change legend title font size
        legend.text = element_text(size=20), #change legend text font sizestrip.background = element_rect(color="black", fill="white", linewidth=1),
        strip.background = element_rect(color="black", fill="white", linewidth=1),
        strip.text = element_text(size=16, color="black"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=22),
        axis.title.x = element_text(size=22),
        axis.text.y = element_text(size=15, color="black"),
        axis.text.x = element_text(size=20, color="black"),
        axis.ticks.y = element_blank())

DOS_plot


# Assemble plots
row0 <- gen_div_plot2
row1 <- (ratio1 | plot_spacer() | ratio2) + plot_layout(widths = c(1, 0.025, 1))
row2 <- sfs_plot
row3 <- DFE_plot
row4 <- DOS_plot

fig3 <- row0 / row1 / row2 / row3 / row4 +
  plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 20, face="bold"), plot.tag.position = c(0.00, 1), legend.position = "bottom")


png("Figures/Figure3.png", width=7000, height=7000, res=300)
fig3
dev.off()
