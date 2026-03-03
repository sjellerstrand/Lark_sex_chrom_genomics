## Export variables and load libraries
rm(list=ls())
library(tidyverse)
library(patchwork)
options(scipen=999)
setwd("C:/Users/Simon JE/OneDrive - Lund University/Dokument/Simon/PhD/Projects/Skylark_2021/Results/Figures/Supplementary/PhaseWY evaluation/")


## Get sample data
Skylark_inds <- read.table("C:/Users/Simon JE/OneDrive - Lund University/Dokument/Simon/PhD/Projects/Skylark_2021/Results/Figures/Skylark_2021_sample_info.txt", sep="\t", head=T)
Rasolark_inds <- read.table("C:/Users/Simon JE/OneDrive - Lund University/Dokument/Simon/PhD/Projects/Skylark_2021/Results/Figures/Rasolark_2021_sample_info.txt", sep="\t", head=T)
samples <- rbind(Skylark_inds, Rasolark_inds)

### PCA
# Import data
dataS3pc <- read.table("Skylark_2021_3.eigenvec")
dataS3pc$Data <- "All variants"
dataS3eig <- read.table("Skylark_2021_3.eigenval")
dataR3pc <- read.table("Rasolark_2021_3.eigenvec")
dataR3pc$Data <- "All variants"
dataR3eig <- read.table("Rasolark_2021_3.eigenval")
dataSApc <- read.table("Skylark_2021_autosomal_CDS.eigenvec")
dataSApc$Data <- "Autosomal CDS variants"
dataSAeig <- read.table("Skylark_2021_autosomal_CDS.eigenval")
dataRApc <- read.table("Rasolark_2021_autosomal_CDS.eigenvec")
dataRApc$Data <- "Autosomal CDS variants"
dataRAeig <- read.table("Rasolark_2021_autosomal_CDS.eigenval")
dataSZpc <- read.table("Skylark_2021_homogametic_CDS.eigenvec")
dataSZpc$Data <- "Z CDS variants"
dataSZeig <- read.table("Skylark_2021_homogametic_CDS.eigenval")
dataRZpc <- read.table("Rasolark_2021_homogametic_CDS.eigenvec")
dataRZpc$Data <- "Z CDS variants"
dataRZeig <- read.table("Rasolark_2021_homogametic_CDS.eigenval")

eigenvals <- as.data.frame(rbind(t(dataS3eig), t(dataR3eig), t(dataSAeig), t(dataRAeig), t(dataSZeig), t(dataRZeig)))
eigenvals <- cbind(eigenvals, c("All variants", "All variants", "Autosomal CDS variants", "Autosomal CDS variants", "Z CDS variants", "Z CDS variants"))
eigenvals <- cbind(eigenvals, c("Skylark", "Raso lark", "Skylark", "Raso lark", "Skylark", "Raso lark"))
colnames(eigenvals) <- c(paste0("PC", 1:(ncol(eigenvals)-2)), "Data", "Species")

for(i in 1:nrow(eigenvals)) {
  eigenvals[i,1:(ncol(eigenvals)-2)] <- eigenvals[i,1:(ncol(eigenvals)-2)]/sum(eigenvals[i,1:(ncol(eigenvals)-2)])*100
  eigenvals[i,1:(ncol(eigenvals)-2)] <- eigenvals[i,1:(ncol(eigenvals)-2)]/sum(eigenvals[i,1:(ncol(eigenvals)-2)])*100
}

PC_data <- rbind(dataS3pc, dataR3pc, dataSApc, dataRApc, dataSZpc, dataRZpc)
PC_data <- PC_data[,-1]
colnames(PC_data) <- c("ID", paste0("PC", 1:(ncol(PC_data)-2)), "Data")
PC_data$Sex <- rep(NA, nrow(PC_data))
PC_data$Species <- rep(NA, nrow(PC_data))

for(i in 1:nrow(PC_data)) {
  PC_data$Sex[i] <- samples$Sex[which(samples$ID == PC_data$ID[i])]
  PC_data$Species[i] <- samples$Species[which(samples$ID == PC_data$ID[i])]
}

PC_data$Species[which(PC_data$Species == "Alauda arvensis")] <- "Skylark"
PC_data$Species[which(PC_data$Species == "Alauda razae")] <- "Raso lark"
PC_data$Species <- factor(PC_data$Species, order=T, levels=c("Skylark", "Raso lark"))

plot_list1 <- vector(mode="list", length=6)

for(i in 1:nrow(eigenvals)) {
  loop_data <- PC_data[which(PC_data$Data == eigenvals$Data[i] & PC_data$Species == eigenvals$Species[i]),]
  plot_list1[[i]] <- ggplot() +
    geom_point(data=loop_data, aes(x=PC1, y=PC2, color=Sex), position="identity", size=2) +
    scale_color_manual(values = c("#d7191c", "#2c7bb6"), limits=c("Female", "Male")) +
    facet_wrap(Species~Data) +
    labs(x=paste("PC1 (", sprintf("%.2f", eigenvals$PC1[i]), "%)", sep=""), y=paste("PC2 (", sprintf("%.2f", eigenvals$PC2[i]), "%)", sep="")) +
    scale_y_continuous(expand = c(0.05,0.05)) +
    scale_x_continuous(expand = c(0.05,0.05)) +
    theme_bw() +
    theme(legend.key.size = unit(1, 'line'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=20), #change legend title font size
          legend.text = element_text(size=20), #change legend text font size
          panel.spacing = unit(1, "lines"),
          strip.background = element_rect(color="black", fill="white", linewidth=1),
          strip.text.x = element_text(size=20, color="black"),
          strip.text.y = element_text(size=20, color="black"),
          axis.line = element_line(colour = "black"),
          axis.title.y = element_text(size=20),
          axis.title.x = element_text(size=20),
          axis.text.y = element_text(size=15, color="black"),
          axis.text.x =  element_text(size=15, color="black"))
}

png("PCAs.png", width=3000, height=4200, res=300)
(plot_list1[[1]] + plot_list1[[2]]) / (plot_list1[[3]] + plot_list1[[4]]) / (plot_list1[[5]] + plot_list1[[6]]) +
  plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 20, face="bold"), plot.tag.position=c(0.05,0.95))
dev.off()

### Folded SFS
dataS3 <- read.table("Skylark_2021_3.frq", head=T, row.names=NULL)
dataS3$Data <- "All variants"
dataS3$Species <- "Skylark"
dataR3 <- read.table("Rasolark_2021_3.frq", head=T, row.names=NULL)
dataR3$Data <- "All variants"
dataR3$Species <- "Raso lark"
dataSA <- read.table("Skylark_2021_autosomal_CDS.frq", head=T, row.names=NULL)
dataSA$Data <- "Autosomal CDS variants"
dataSA$Species <- "Skylark"
dataRA <- read.table("Rasolark_2021_autosomal_CDS.frq", head=T, row.names=NULL)
dataRA$Data <- "Autosomal CDS variants"
dataRA$Species <- "Raso lark"
dataSZ <- read.table("Skylark_2021_homogametic_CDS.frq", head=T, row.names=NULL)
dataSZ$Data <- "Z CDS variants"
dataSZ$Species <- "Skylark"
dataRZ <- read.table("Rasolark_2021_homogametic_CDS.frq", head=T, row.names=NULL)
dataRZ$Data <- "Z CDS variants"
dataRZ$Species <- "Raso lark"

SFS <- rbind(dataS3, dataR3, dataSA, dataRA, dataSZ, dataRZ)
colnames(SFS) <- c("CHROM", "POS", "N_ALLELES", "N_CHR", "FREQA", "FREQB", "Data", "Species")
SFS$UFSFS <- SFS |> dplyr::select(FREQA, FREQB) |> apply(1, function(z) min(z))

SFS$Species[which(SFS$Species == "Alauda arvensis")] <- "Skylark"
SFS$Species[which(SFS$Species == "Alauda razae")] <- "Raso lark"
SFS$Species <- factor(SFS$Species, order=T, levels=c("Skylark", "Raso lark"))

W_frequency <- as.data.frame(c(10/(18*2), 9/(18*2), 10/(18*2), 9/(18*2), 10/(18*1.5), 9/(18*1.5)))
colnames(W_frequency) <- "freq"


plot_list2 <- vector(mode="list", length=6)

i <- 1
loop_data <- SFS[which(SFS$Data == eigenvals$Data[i] & SFS$Species == eigenvals$Species[i]),]
plot_list2[[i]] <- ggplot() +
  geom_density(data=loop_data, aes(x=UFSFS, y=after_stat(count/sum(count))), adjust=1.5, linewidth=1.5) +
  geom_vline(data=loop_data, aes(xintercept=W_frequency$freq[1], color="#d7191c"), linewidth=1, linetype=2) +
  facet_wrap(Species~Data) +
  labs(x ="Minor allele frequency", y = "Proportion") +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.spacing = unit(1, "lines"),
        strip.background = element_rect(color="black", fill="white", linewidth=1),
        strip.text = element_text(size=20, color="black"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=15, color="black"),
        axis.text.x =  element_text(size=15, color="black"))

i <- 2
loop_data <- SFS[which(SFS$Data == eigenvals$Data[i] & SFS$Species == eigenvals$Species[i]),]
plot_list2[[i]] <- ggplot() +
  geom_density(data=loop_data, aes(x=UFSFS, y=after_stat(count/sum(count))), adjust=1.5, linewidth=1.5) +
  geom_vline(data=loop_data, aes(xintercept=W_frequency$freq[2], color="#d7191c"), linewidth=1, linetype=2) +
  facet_wrap(Species~Data) +
  labs(x ="Minor allele frequency", y = "Proportion") +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.spacing = unit(1, "lines"),
        strip.background = element_rect(color="black", fill="white", linewidth=1),
        strip.text = element_text(size=20, color="black"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=15, color="black"),
        axis.text.x =  element_text(size=15, color="black"))

i <- 3
loop_data <- SFS[which(SFS$Data == eigenvals$Data[i] & SFS$Species == eigenvals$Species[i]),]
plot_list2[[i]] <- ggplot() +
  geom_density(data=loop_data, aes(x=UFSFS, y=after_stat(count/sum(count))), adjust=1.5, linewidth=1.5) +
  geom_vline(data=loop_data, aes(xintercept=W_frequency$freq[3], color="#d7191c"), linewidth=1, linetype=2) +
  facet_wrap(Species~Data) +
  labs(x ="Minor allele frequency", y = "Proportion") +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.spacing = unit(1, "lines"),
        strip.background = element_rect(color="black", fill="white", linewidth=1),
        strip.text = element_text(size=20, color="black"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=15, color="black"),
        axis.text.x =  element_text(size=15, color="black"))

i <- 4
loop_data <- SFS[which(SFS$Data == eigenvals$Data[i] & SFS$Species == eigenvals$Species[i]),]
plot_list2[[i]] <- ggplot() +
  geom_density(data=loop_data, aes(x=UFSFS, y=after_stat(count/sum(count))), adjust=1.5, linewidth=1.5) +
  geom_vline(data=loop_data, aes(xintercept=W_frequency$freq[4], color="#d7191c"), linewidth=1, linetype=2) +
  facet_wrap(Species~Data) +
  labs(x ="Minor allele frequency", y = "Proportion") +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.spacing = unit(1, "lines"),
        strip.background = element_rect(color="black", fill="white", linewidth=1),
        strip.text = element_text(size=20, color="black"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=15, color="black"),
        axis.text.x =  element_text(size=15, color="black"))

i <- 5
loop_data <- SFS[which(SFS$Data == eigenvals$Data[i] & SFS$Species == eigenvals$Species[i]),]
plot_list2[[i]] <- ggplot() +
  geom_density(data=loop_data, aes(x=UFSFS, y=after_stat(count/sum(count))), adjust=1.5, linewidth=1.5) +
  geom_vline(data=loop_data, aes(xintercept=W_frequency$freq[5], color="#d7191c"), linewidth=1, linetype=2) +
  facet_wrap(Species~Data) +
  labs(x ="Minor allele frequency", y = "Proportion") +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.spacing = unit(1, "lines"),
        strip.background = element_rect(color="black", fill="white", linewidth=1),
        strip.text = element_text(size=20, color="black"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=15, color="black"),
        axis.text.x =  element_text(size=15, color="black"))

i <- 6
loop_data <- SFS[which(SFS$Data == eigenvals$Data[i] & SFS$Species == eigenvals$Species[i]),]
plot_list2[[i]] <- ggplot() +
  geom_density(data=loop_data, aes(x=UFSFS, y=after_stat(count/sum(count))), adjust=1.5, linewidth=1.5) +
  geom_vline(data=loop_data, aes(xintercept=W_frequency$freq[6], color="#d7191c"), linewidth=1, linetype=2) +
  facet_wrap(Species~Data) +
  labs(x ="Minor allele frequency", y = "Proportion") +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.spacing = unit(1, "lines"),
        strip.background = element_rect(color="black", fill="white", linewidth=1),
        strip.text = element_text(size=20, color="black"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=15, color="black"),
        axis.text.x =  element_text(size=15, color="black"))


png("SFS.png", width=3000, height=4500, res=300)
(plot_list2[[1]] + plot_list2[[2]]) / (plot_list2[[3]] + plot_list2[[4]]) / (plot_list2[[5]] + plot_list2[[6]]) +
  plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 20, face="bold"), plot.tag.position=c(0.05,0.95))
dev.off()


### Inbreeding coefficient
dataS3 <- read.table("Skylark_2021_3.het", head=T)
dataS3$Data <- "All variants"
dataR3 <- read.table("Rasolark_2021_3.het", head=T)
dataR3$Data <- "All variants"
dataSA <- read.table("Skylark_2021_autosomal_CDS.het", head=T)
dataSA$Data <- "Autosomal CDS variants"
dataRA <- read.table("Rasolark_2021_autosomal_CDS.het", head=T)
dataRA$Data <- "Autosomal CDS variants"

InbCoef <- rbind(dataS3, dataR3, dataSA, dataRA)
InbCoef$Sex <- rep(NA, nrow(InbCoef))
InbCoef$Species <- rep(NA, nrow(InbCoef))

for(i in 1:nrow(InbCoef)) {
  InbCoef$Sex[i] <- samples$Sex[which(samples$ID == InbCoef$INDV[i])]
  InbCoef$Species[i] <- samples$Species[which(samples$ID == InbCoef$INDV[i])]
}
InbCoef$Species[which(InbCoef$Species == "Alauda arvensis")] <- "Skylark"
InbCoef$Species[which(InbCoef$Species == "Alauda razae")] <- "Raso lark"
InbCoef$Species <- factor(InbCoef$Species, order=T, levels=c("Skylark", "Raso lark"))



plot_list3 <- vector(mode="list", length=6)
for(i in 1:nrow(eigenvals)) {
  loop_data <- InbCoef[which(InbCoef$Data == eigenvals$Data[i] & InbCoef$Species == eigenvals$Species[i]),]
  plot_list3[[i]] <- ggplot() +
    geom_histogram(data=loop_data, aes(x=F, fill=Sex), position="identity") +
    scale_fill_manual(values = c("#d7191c", "#2c7bb6"), limits=c("Female", "Male")) +
    facet_wrap(Species~Data) +
    labs(x=expression("Inbreeding Coefficient (F"[IS]*")"), y = "Count") +
    scale_y_continuous(expand = c(0,0), breaks=seq(0, 100, 1)) +
    theme_bw() +
    theme(legend.key.size = unit(1, 'line'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=20), #change legend title font size
          legend.text = element_text(size=20), #change legend text font size
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.spacing = unit(1, "lines"),
          strip.background = element_rect(color="black", fill="white", linewidth=1),
          strip.text = element_text(size=20, color="black"),
          axis.line = element_line(colour = "black"),
          axis.title.y = element_text(size=20),
          axis.title.x = element_text(size=20),
          axis.text.y = element_text(size=15, color="black"),
          axis.text.x =  element_text(size=15, color="black"))

}

png("Inbreeding_coefficient.png", width=3000, height=3000, res=300)
(plot_list3[[1]] + plot_list3[[2]]) / (plot_list3[[3]] + plot_list3[[4]]) +
  plot_layout(guides = "collect", axis_titles = "collect") + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 20, face="bold"), plot.tag.position=c(0.05,0.95))
dev.off()
