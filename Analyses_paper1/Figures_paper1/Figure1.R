## Load libraries
rm(list=ls())
library(tidyverse)
library(scales)
library(ggpattern)

options(scipen=999)
setwd("C:/Users/Simon JE/OneDrive - Lund University/Dokument/Simon/PhD/Projects/Skylark_2021/Results/")

#Figure 1 main
shift_right <- 50000000

### import data on synteny genome
genome <- read.delim("Coordinate lift over/Pmaj.fa.fai", sep="\t", head=F)
chrom_names <- read.delim("Coordinate lift over/Pmaj_synteny_chroms.txt", sep="\t", head=F)
colnames(chrom_names) <- c("scaffold", "name")
chrom_names <- chrom_names[-which(chrom_names$name == "MT"),]
chrom_names$size <- rep(NA, nrow(chrom_names))
chrom_names$start <- rep(NA, nrow(chrom_names))
chrom_names$end <- rep(NA, nrow(chrom_names))
chrom_names$CumMid <- rep(NA, nrow(chrom_names))
chrom_names <- chrom_names[c(1, 7, 2:6, 8:23, 31, 24:26, 28, 27, 29:30, 32),]

x <- shift_right
for(i in 1:nrow(chrom_names)) {
  chrom_names$size[i] <- genome[which(genome[,1] == chrom_names$scaffold[i]),2]
  chrom_names$start[i] <- x + 1
  chrom_names$end[i] <- chrom_names$size[i] + x 
  chrom_names$CumMid[i] <- (chrom_names$end[i] - chrom_names$start[i])/2 + chrom_names$start[i]
  x <- chrom_names$end[i]
}

# Skylark
# genome summary
data1 <- read.delim("Coordinate lift over/Skylark_2021_Pmaj_Genome_summary_final.txt", sep="\t", head=F)
colnames(data1) <- c("scaffold", "start" , "end", "data_type")
data1$CumStart <- rep(NA, nrow(data1))
data1$CumEnd <- rep(NA, nrow(data1))
for(i in 1:nrow(chrom_names)) {
  rows <- which(data1$scaffold == chrom_names$scaffold[i])
  data1$scaffold[rows] <- chrom_names$name[i]
  data1$CumStart[which(data1$scaffold == chrom_names$name[i])] <- data1$start[which(data1$scaffold == chrom_names$name[i])] + chrom_names$start[i] - 1
  data1$CumEnd[which(data1$scaffold == chrom_names$name[i])] <- data1$end[which(data1$scaffold == chrom_names$name[i])] + chrom_names$start[i] - 1
}
data1 <- data1[order(factor(data1$scaffold, levels=unique(chrom_names$name))),]
data1$scaffold <- factor(data1$scaffold, levels=chrom_names$name)
data1$data_type[which(data1$data_type == "both")] <- "Sex haplotype & coverage difference"
data1$data_type[which(data1$data_type == "hetgamdrop")] <- "Sex sequencing\ndepth difference"
data1$data_type[which(data1$data_type == "phase")] <- "Sex haplotype clustering"
data1$data_type[which(data1$data_type == "Missing data")] <- "No data"
data1$start <- as.numeric(data1$start)
data1$end <- as.numeric(data1$end)

# findzx
### Hetrozygosity
data2 <- read.delim("findzx/Skylark_Parus/tables/diffHeterozygosity.1000000bp.out", sep="\t", head=T)
data2 <- data2[order(factor(data2$chr, levels=unique(chrom_names$name))),]
data2$mid <- as.numeric(data2$start + (data2$end - data2$start + 1)/2)
data2$chr <- factor(data2$chr, levels=chrom_names$name)
data2$diff <- as.numeric(data2$diff)
data2$mid <- as.numeric(data2$mid)
data2$data <- "Sex heterozygosity\ndifference"
data2$CumPos <- rep(NA, nrow(data2))
for(i in 1:nrow(chrom_names)) {
  data2$CumPos[which(data2$chr == chrom_names$name[i])] <-  data2$mid[which(data2$chr == chrom_names$name[i])] + + chrom_names$start[i]
}

## Depth
data3 <- read.delim("findzx/Skylark_Parus/tables/diffGenomeCoverage.mismatch.unfiltered.1000000bp.out", sep="\t", head=T)
data3 <- data3[order(factor(data3$chr, levels=unique(chrom_names$name))),]
data3$mid <- as.numeric(data3$start + (data3$end - data3$start + 1)/2)
data3$chr <- factor(data3$chr, levels=chrom_names$name)
data3$diff <- as.numeric(data3$diff)
data3$mid <- as.numeric(data3$mid)
data3$data <- "Sex sequencing\ndepth difference"
data3$CumPos <- rep(NA, nrow(data3))
for(i in 1:nrow(chrom_names)) {
  data3$CumPos[which(data3$chr == chrom_names$name[i])] <-  data3$mid[which(data3$chr == chrom_names$name[i])] + chrom_names$start[i]
}

# Rasolark
# genome summary
data4 <- read.delim("Coordinate lift over/Rasolark_2021_Pmaj_Genome_summary_final.txt", sep="\t", head=F)
colnames(data4) <- c("scaffold", "start" , "end", "data_type")
data4$CumStart <- rep(NA, nrow(data4))
data4$CumEnd <- rep(NA, nrow(data4))
for(i in 1:nrow(chrom_names)) {
  rows <- which(data4$scaffold == chrom_names$scaffold[i])
  data4$scaffold[rows] <- chrom_names$name[i]
  data4$CumStart[which(data4$scaffold == chrom_names$name[i])] <- data4$start[which(data4$scaffold == chrom_names$name[i])] + chrom_names$start[i] - 1
  data4$CumEnd[which(data4$scaffold == chrom_names$name[i])] <- data4$end[which(data4$scaffold == chrom_names$name[i])] + chrom_names$start[i] - 1
}

data4 <- data4[order(factor(data4$scaffold, levels=unique(chrom_names$name))),]
data4$scaffold <- factor(data4$scaffold, levels=chrom_names$name)
data4$data_type[which(data4$data_type == "both")] <- "Sex haplotype & coverage difference"
data4$data_type[which(data4$data_type == "hetgamdrop")] <- "Sex sequencing\ndepth difference"
data4$data_type[which(data4$data_type == "phase")] <- "Sex haplotype clustering"
data4$data_type[which(data4$data_type == "Missing data")] <- "No data"
data4$start <- as.numeric(data4$start)
data4$end <- as.numeric(data4$end)

# findzx
### Hetrozygosity
data5 <- read.delim("findzx/Rasolark_Parus/tables/diffHeterozygosity.1000000bp.out", sep="\t", head=T)
data5 <- data5[order(factor(data5$chr, levels=unique(chrom_names$name))),]
data5$mid <- as.numeric(data5$start + (data5$end - data5$start + 1)/2)
data5$chr <- factor(data5$chr, levels=chrom_names$name)
data5$diff <- as.numeric(data5$diff)
data5$mid <- as.numeric(data5$mid)
data5$data <- "Sex heterozygosity\ndifference"
data5$CumPos <- rep(NA, nrow(data5))
for(i in 1:nrow(chrom_names)) {
  data5$CumPos[which(data5$chr == chrom_names$name[i])] <-  data5$mid[which(data5$chr == chrom_names$name[i])] + chrom_names$start[i]
}

## Depth
data6 <- read.delim("findzx/Rasolark_Parus/tables/diffGenomeCoverage.mismatch.unfiltered.1000000bp.out", sep="\t", head=T)
data6 <- data6[order(factor(data6$chr, levels=unique(chrom_names$name))),]
data6$mid <- as.numeric(data6$start + (data6$end - data6$start + 1)/2)
data6$chr <- factor(data6$chr, levels=chrom_names$name)
data6$diff <- as.numeric(data6$diff)
data6$mid <- as.numeric(data6$mid)
data6$data <- "Sex sequencing\ndepth difference"
data6$CumPos <- rep(NA, nrow(data6))
for(i in 1:nrow(chrom_names)) {
  data6$CumPos[which(data6$chr == chrom_names$name[i])] <-  data6$mid[which(data6$chr == chrom_names$name[i])] + chrom_names$start[i]
}

chrom_names$labels <- as.character(chrom_names$name)
chrom_names$labels[22:31] <- ""
chrom_names$labels[26] <- "21-28"

# Remove very small sections which makes the plot size smaller
data1 <- data1[which((data1$end-data1$start) > 100),]
data4 <- data4[which((data4$end-data4$start) > 100),]

data1reg1 <- data1[1:(min(which(data1$scaffold=="3"))-1),]
data1chr3 <- data1[which(data1$scaffold=="3"),]
data1chr4 <- data1[which(data1$scaffold=="4"),]
data1chr4A5 <- data1[which(data1$scaffold=="4A" | data1$scaffold=="5"),]
data1reg2 <- data1[(max(which(data1$scaffold=="5"))+1):(min(which(data1$scaffold=="Z"))-1),]
data1chrZ <- data1[which(data1$scaffold=="Z"),]

data2reg1 <- data2[1:(min(which(data2$chr=="3"))-1),]
data2chr3 <- data2[which(data2$chr=="3"),]
data2chr4 <- data2[which(data2$chr=="4"),]
data2chr4A5 <- data2[which(data2$chr=="4A" | data2$chr=="5"),]
data2reg2 <- data2[(max(which(data2$chr=="5"))+1):(min(which(data2$chr=="Z"))-1),]
data2chrZ <- data2[which(data2$chr=="Z"),]

data3reg1 <- data3[1:(min(which(data3$chr=="3"))-1),]
data3chr3 <- data3[which(data3$chr=="3"),]
data3chr4 <- data3[which(data3$chr=="4"),]
data3chr4A5 <- data3[which(data3$chr=="4A" | data3$chr=="5"),]
data3reg2 <- data3[(max(which(data3$chr=="5"))+1):(min(which(data3$chr=="Z"))-1),]
data3chrZ <- data3[which(data3$chr=="Z"),]

data4reg1 <- data4[1:(min(which(data4$scaffold=="3"))-1),]
data4chr3 <- data4[which(data4$scaffold=="3"),]
data4chr4 <- data4[which(data4$scaffold=="4"),]
data4chr4A5 <- data4[which(data4$scaffold=="4A" | data4$scaffold=="5"),]
data4reg2 <- data4[(max(which(data4$scaffold=="5"))+1):(min(which(data4$scaffold=="Z"))-1),]
data4chrZ <- data4[which(data4$scaffold=="Z"),]

data5reg1 <- data5[1:(min(which(data5$chr=="3"))-1),]
data5chr3 <- data5[which(data5$chr=="3"),]
data5chr4 <- data5[which(data5$chr=="4"),]
data5chr4A5 <- data5[which(data5$chr=="4A" | data5$chr=="5"),]
data5reg2 <- data5[(max(which(data5$chr=="5"))+1):(min(which(data5$chr=="Z"))-1),]
data5chrZ <- data5[which(data5$chr=="Z"),]

data6reg1 <- data6[1:(min(which(data6$chr=="3"))-1),]
data6chr3 <- data6[which(data6$chr=="3"),]
data6chr4 <- data6[which(data6$chr=="4"),]
data6chr4A5 <- data6[which(data6$chr=="4A" | data6$chr=="5"),]
data6reg2 <- data6[(max(which(data6$chr=="5"))+1):(min(which(data6$chr=="Z"))-1),]
data6chrZ <- data6[which(data6$chr=="Z"),]

chrom_namesreg1 <- chrom_names[1:(which(chrom_names$name=="3")-1),]
chrom_nameschr3 <- chrom_names[which(chrom_names$name=="3"),]
chrom_nameschr4 <- chrom_names[which(chrom_names$name=="4"),]
chrom_nameschr4A5 <- chrom_names[which(chrom_names$name=="4A" | chrom_names$name=="5"),]
chrom_namesreg2 <- chrom_names[(max(which(chrom_names$name=="5"))+1):(min(which(chrom_names$name=="Z"))-1),]
chrom_nameschrZ <- chrom_names[which(chrom_names$name=="Z"),]

strata3 <- chrom_nameschr3
strata3$size <- NA
strata3 <- rbind(strata3, strata3, strata3)
strata3$start <- c(11165000, 18425000, 22525000)+chrom_nameschr3$start[1]
strata3$end <- c(18425000, 22525000, 87275000)+chrom_nameschr3$start[1]
strata3$CumLabel <- (strata3$end-strata3$start)/2 + strata3$start
strata3$CumLabel[1] <- strata3$CumLabel[1] - 5500000
strata3$CumLabel[2] <- strata3$CumLabel[2] - 2500000
strata3$labels <- c("3-a", "3-b", "3-c")

yposadd <- 0

plot <- ggplot() +
  # PhaseWY results
  geom_rect(data=data1reg1, aes(xmin=CumStart, ymin=7, xmax=CumEnd, ymax=8, fill=data_type)) +
  geom_rect(data=data4reg1, aes(xmin=CumStart, ymin=4, xmax=CumEnd, ymax=5, fill=data_type)) +
  geom_rect(data=data1chr3, aes(xmin=CumStart, ymin=7+yposadd, xmax=CumEnd, ymax=8+yposadd, fill=data_type)) +
  geom_rect(data=data4chr3, aes(xmin=CumStart, ymin=4+yposadd, xmax=CumEnd, ymax=5+yposadd, fill=data_type)) +
  geom_rect(data=data1chr4, aes(xmin=CumStart, ymin=7, xmax=CumEnd, ymax=8, fill=data_type)) +
  geom_rect(data=data4chr4, aes(xmin=CumStart, ymin=4, xmax=CumEnd, ymax=5, fill=data_type)) +
  geom_rect(data=data1chr4A5, aes(xmin=CumStart, ymin=7+yposadd, xmax=CumEnd, ymax=8+yposadd, fill=data_type)) +
  geom_rect(data=data4chr4A5, aes(xmin=CumStart, ymin=4+yposadd, xmax=CumEnd, ymax=5+yposadd, fill=data_type)) +
  geom_rect(data=data1reg2, aes(xmin=CumStart, ymin=7, xmax=CumEnd, ymax=8, fill=data_type)) +
  geom_rect(data=data4reg2, aes(xmin=CumStart, ymin=4, xmax=CumEnd, ymax=5, fill=data_type)) +
  geom_rect(data=data1chrZ, aes(xmin=CumStart, ymin=7+yposadd, xmax=CumEnd, ymax=8+yposadd, fill=data_type)) +
  geom_rect(data=data4chrZ, aes(xmin=CumStart, ymin=4+yposadd, xmax=CumEnd, ymax=5+yposadd, fill=data_type)) +
  scale_fill_manual(name="PhaseWY", values = c("Autosomal"="#E4EAF0", "Sex haplotype clustering"="#fecc5c", "Sex sequencing\ndepth difference"="#b30000", "No data"="#404040", "#f03b20"),
                    limits = c("Autosomal", "Sex haplotype clustering", "Sex sequencing\ndepth difference", "No data")) +
  # Findzx results bars
  geom_rect(data=data3reg1, aes(xmin=CumPos-500000, xmax=CumPos+500000, ymin=6.5, ymax=((diff/((max(data3$diff, na.rm=T)-min(data3$diff, na.rm=T))*2))+6.5), color=data), lwd=1.5) +
  geom_rect(data=data2reg1, aes(xmin=CumPos-500000, xmax=CumPos+500000, ymin=5.5, ymax=((diff/((max(data2$diff, na.rm=T)-min(data2$diff, na.rm=T))*2))+5.5), color=data), lwd=1.5) +
  geom_rect(data=data6reg1, aes(xmin=CumPos-500000, xmax=CumPos+500000, ymin=3.5, ymax=((diff/((max(data6$diff, na.rm=T)-min(data6$diff, na.rm=T))*2))+3.5), color=data), lwd=1.5) +
  geom_rect(data=data5reg1, aes(xmin=CumPos-500000, xmax=CumPos+500000, ymin=2.5, ymax=((diff/((max(data5$diff, na.rm=T)-min(data5$diff, na.rm=T))*2))+2.5), color=data), lwd=1.5) +
  geom_rect(data=data3chr3, aes(xmin=CumPos-500000, xmax=CumPos+500000, ymin=6.5+yposadd, ymax=((diff/((max(data3$diff, na.rm=T)-min(data3$diff, na.rm=T))*2))+6.5+yposadd), color=data), lwd=1.5) +
  geom_rect(data=data2chr3, aes(xmin=CumPos-500000, xmax=CumPos+500000, ymin=5.5+yposadd, ymax=((diff/((max(data2$diff, na.rm=T)-min(data2$diff, na.rm=T))*2))+5.5+yposadd), color=data), lwd=1.5) +
  geom_rect(data=data6chr3, aes(xmin=CumPos-500000, xmax=CumPos+500000, ymin=3.5+yposadd, ymax=((diff/((max(data6$diff, na.rm=T)-min(data6$diff, na.rm=T))*2))+3.5+yposadd), color=data), lwd=1.5) +
  geom_rect(data=data5chr3, aes(xmin=CumPos-500000, xmax=CumPos+500000, ymin=2.5+yposadd, ymax=((diff/((max(data5$diff, na.rm=T)-min(data5$diff, na.rm=T))*2))+2.5+yposadd), color=data), lwd=1.5) +
  geom_rect(data=data3chr4, aes(xmin=CumPos-500000, xmax=CumPos+500000, ymin=6.5, ymax=((diff/((max(data3$diff, na.rm=T)-min(data3$diff, na.rm=T))*2))+6.5), color=data), lwd=1.5) +
  geom_rect(data=data2chr4, aes(xmin=CumPos-500000, xmax=CumPos+500000, ymin=5.5, ymax=((diff/((max(data2$diff, na.rm=T)-min(data2$diff, na.rm=T))*2))+5.5), color=data), lwd=1.5) +
  geom_rect(data=data6chr4, aes(xmin=CumPos-500000, xmax=CumPos+500000, ymin=3.5, ymax=((diff/((max(data6$diff, na.rm=T)-min(data6$diff, na.rm=T))*2))+3.5), color=data), lwd=1.5) +
  geom_rect(data=data5chr4, aes(xmin=CumPos-500000, xmax=CumPos+500000, ymin=2.5, ymax=((diff/((max(data5$diff, na.rm=T)-min(data5$diff, na.rm=T))*2))+2.5), color=data), lwd=1.5) +
  geom_rect(data=data3chr4A5, aes(xmin=CumPos-500000, xmax=CumPos+500000, ymin=6.5+yposadd, ymax=((diff/((max(data3$diff, na.rm=T)-min(data3$diff, na.rm=T))*2))+6.5+yposadd), color=data), lwd=1.5) +
  geom_rect(data=data2chr4A5, aes(xmin=CumPos-500000, xmax=CumPos+500000, ymin=5.5+yposadd, ymax=((diff/((max(data2$diff, na.rm=T)-min(data2$diff, na.rm=T))*2))+5.5+yposadd), color=data), lwd=1.5) +
  geom_rect(data=data6chr4A5, aes(xmin=CumPos-500000, xmax=CumPos+500000, ymin=3.5+yposadd, ymax=((diff/((max(data6$diff, na.rm=T)-min(data6$diff, na.rm=T))*2))+3.5+yposadd), color=data), lwd=1.5) +
  geom_rect(data=data5chr4A5, aes(xmin=CumPos-500000, xmax=CumPos+500000, ymin=2.5+yposadd, ymax=((diff/((max(data5$diff, na.rm=T)-min(data5$diff, na.rm=T))*2))+2.5+yposadd), color=data), lwd=1.5) +
  geom_rect(data=data3reg2, aes(xmin=CumPos-500000, xmax=CumPos+500000, ymin=6.5, ymax=((diff/((max(data3$diff, na.rm=T)-min(data3$diff, na.rm=T))*2))+6.5), color=data), lwd=1.5) +
  geom_rect(data=data2reg2, aes(xmin=CumPos-500000, xmax=CumPos+500000, ymin=5.5, ymax=((diff/((max(data2$diff, na.rm=T)-min(data2$diff, na.rm=T))*2))+5.5), color=data), lwd=1.5) +
  geom_rect(data=data6reg2, aes(xmin=CumPos-500000, xmax=CumPos+500000, ymin=3.5, ymax=((diff/((max(data6$diff, na.rm=T)-min(data6$diff, na.rm=T))*2))+3.5), color=data), lwd=1.5) +
  geom_rect(data=data5reg2, aes(xmin=CumPos-500000, xmax=CumPos+500000, ymin=2.5, ymax=((diff/((max(data5$diff, na.rm=T)-min(data5$diff, na.rm=T))*2))+2.5), color=data), lwd=1.5) +
  geom_rect(data=data3chrZ, aes(xmin=CumPos-500000, xmax=CumPos+500000, ymin=6.5+yposadd, ymax=((diff/((max(data3$diff, na.rm=T)-min(data3$diff, na.rm=T))*2))+6.5+yposadd), color=data), lwd=1.5) +
  geom_rect(data=data2chrZ, aes(xmin=CumPos-500000, xmax=CumPos+500000, ymin=5.5+yposadd, ymax=((diff/((max(data2$diff, na.rm=T)-min(data2$diff, na.rm=T))*2))+5.5+yposadd), color=data), lwd=1.5) +
  geom_rect(data=data6chrZ, aes(xmin=CumPos-500000, xmax=CumPos+500000, ymin=3.5+yposadd, ymax=((diff/((max(data6$diff, na.rm=T)-min(data6$diff, na.rm=T))*2))+3.5+yposadd), color=data), lwd=1.5) +
  geom_rect(data=data5chrZ, aes(xmin=CumPos-500000, xmax=CumPos+500000, ymin=2.5+yposadd, ymax=((diff/((max(data5$diff, na.rm=T)-min(data5$diff, na.rm=T))*2))+2.5+yposadd), color=data), lwd=1.5) +
  # Findzx results lines
  geom_line(data=data3reg1, aes(x=CumPos, y=((diff/((max(data3$diff, na.rm=T)-min(data3$diff, na.rm=T))*2))+6.5), color=data), key_glyph = "polygon", lwd=1.5) +
  geom_line(data=data2reg1, aes(x=CumPos, y=((diff/((max(data2$diff, na.rm=T)-min(data2$diff, na.rm=T))*2))+5.5), color=data), key_glyph = "polygon", lwd=1.5) +
  geom_line(data=data6reg1, aes(x=CumPos, y=((diff/((max(data6$diff, na.rm=T)-min(data6$diff, na.rm=T))*2))+3.5), color=data), key_glyph = "polygon", lwd=1.5) +
  geom_line(data=data5reg1, aes(x=CumPos, y=((diff/((max(data5$diff, na.rm=T)-min(data5$diff, na.rm=T))*2))+2.5), color=data), key_glyph = "polygon", lwd=1.5) +
  geom_line(data=data3chr3, aes(x=CumPos, y=((diff/((max(data3$diff, na.rm=T)-min(data3$diff, na.rm=T))*2))+6.5+yposadd), color=data), key_glyph = "polygon", lwd=1.5) +
  geom_line(data=data2chr3, aes(x=CumPos, y=((diff/((max(data2$diff, na.rm=T)-min(data2$diff, na.rm=T))*2))+5.5+yposadd), color=data), key_glyph = "polygon", lwd=1.5) +
  geom_line(data=data6chr3, aes(x=CumPos, y=((diff/((max(data6$diff, na.rm=T)-min(data6$diff, na.rm=T))*2))+3.5+yposadd), color=data), key_glyph = "polygon", lwd=1.5) +
  geom_line(data=data5chr3, aes(x=CumPos, y=((diff/((max(data5$diff, na.rm=T)-min(data5$diff, na.rm=T))*2))+2.5+yposadd), color=data), key_glyph = "polygon", lwd=1.5) +
  geom_line(data=data3chr4, aes(x=CumPos, y=((diff/((max(data3$diff, na.rm=T)-min(data3$diff, na.rm=T))*2))+6.5), color=data), key_glyph = "polygon", lwd=1.5) +
  geom_line(data=data2chr4, aes(x=CumPos, y=((diff/((max(data2$diff, na.rm=T)-min(data2$diff, na.rm=T))*2))+5.5), color=data), key_glyph = "polygon", lwd=1.5) +
  geom_line(data=data6chr4, aes(x=CumPos, y=((diff/((max(data6$diff, na.rm=T)-min(data6$diff, na.rm=T))*2))+3.5), color=data), key_glyph = "polygon", lwd=1.5) +
  geom_line(data=data5chr4, aes(x=CumPos, y=((diff/((max(data5$diff, na.rm=T)-min(data5$diff, na.rm=T))*2))+2.5), color=data), key_glyph = "polygon", lwd=1.5) +
  geom_line(data=data3chr4A5, aes(x=CumPos, y=((diff/((max(data3$diff, na.rm=T)-min(data3$diff, na.rm=T))*2))+6.5+yposadd), color=data), key_glyph = "polygon", lwd=1.5) +
  geom_line(data=data2chr4A5, aes(x=CumPos, y=((diff/((max(data2$diff, na.rm=T)-min(data2$diff, na.rm=T))*2))+5.5+yposadd), color=data), key_glyph = "polygon", lwd=1.5) +
  geom_line(data=data6chr4A5, aes(x=CumPos, y=((diff/((max(data6$diff, na.rm=T)-min(data6$diff, na.rm=T))*2))+3.5+yposadd), color=data), key_glyph = "polygon", lwd=1.5) +
  geom_line(data=data5chr4A5, aes(x=CumPos, y=((diff/((max(data5$diff, na.rm=T)-min(data5$diff, na.rm=T))*2))+2.5+yposadd), color=data), key_glyph = "polygon", lwd=1.5) +
  geom_line(data=data3reg2, aes(x=CumPos, y=((diff/((max(data3$diff, na.rm=T)-min(data3$diff, na.rm=T))*2))+6.5), color=data), key_glyph = "polygon", lwd=1.5) +
  geom_line(data=data2reg2, aes(x=CumPos, y=((diff/((max(data2$diff, na.rm=T)-min(data2$diff, na.rm=T))*2))+5.5), color=data), key_glyph = "polygon", lwd=1.5) +
  geom_line(data=data6reg2, aes(x=CumPos, y=((diff/((max(data6$diff, na.rm=T)-min(data6$diff, na.rm=T))*2))+3.5), color=data), key_glyph = "polygon", lwd=1.5) +
  geom_line(data=data5reg2, aes(x=CumPos, y=((diff/((max(data5$diff, na.rm=T)-min(data5$diff, na.rm=T))*2))+2.5), color=data), key_glyph = "polygon", lwd=1.5) +
  geom_line(data=data3chrZ, aes(x=CumPos, y=((diff/((max(data3$diff, na.rm=T)-min(data3$diff, na.rm=T))*2))+6.5+yposadd), color=data), key_glyph = "polygon", lwd=1.5) +
  geom_line(data=data2chrZ, aes(x=CumPos, y=((diff/((max(data2$diff, na.rm=T)-min(data2$diff, na.rm=T))*2))+5.5+yposadd), color=data), key_glyph = "polygon", lwd=1.5) +
  geom_line(data=data6chrZ, aes(x=CumPos, y=((diff/((max(data6$diff, na.rm=T)-min(data6$diff, na.rm=T))*2))+3.5+yposadd), color=data), key_glyph = "polygon", lwd=1.5) +
  geom_line(data=data5chrZ, aes(x=CumPos, y=((diff/((max(data5$diff, na.rm=T)-min(data5$diff, na.rm=T))*2))+2.5+yposadd), color=data), key_glyph = "polygon", lwd=1.5) +
  scale_color_manual(name="FindZX", values = c("Sex sequencing\ndepth difference"="#800026", "Sex heterozygosity\ndifference"="#f03b20")) +
  # Labels
  guides(fill = guide_legend(order = 1), color = guide_legend(order = 2)) +
  annotate(x=chrom_namesreg1$CumMid, label=as.factor(chrom_namesreg1$labels), y=8.75, size=12, geom="text") +
  #annotate(x=strata3$CumLabel, label=as.factor(strata3$labels), y=8.75+yposadd, size=12, geom="text") +
  annotate(x=chrom_nameschr3$CumMid, label=as.factor(chrom_nameschr3$labels), y=8.75, size=12, geom="text") +
  annotate(x=chrom_nameschr4$CumMid, label=as.factor(chrom_nameschr4$labels), y=8.75, size=12, geom="text") +
  annotate(x=chrom_nameschr4A5$CumMid, label=as.factor(chrom_nameschr4A5$labels), y=8.75+yposadd, size=12, geom="text") +
  annotate(x=chrom_namesreg2$CumMid, label=as.factor(chrom_namesreg2$labels), y=8.75, size=12, geom="text") +
  annotate(x=chrom_nameschrZ$CumMid, label=as.factor(chrom_nameschrZ$labels), y=8.75+yposadd, size=12, geom="text") +
  annotate(x=0, y=c(3.2, 4.7, 6.1, 7.6), size=seq(5.5, 6.25, 0.25), geom="text", 
           label=c("Raso lark\nFindZX", "Raso lark\nPhaseWY", "Skylark FindZX", "Skylark PhaseWY")) +
  # Vertical lines
  geom_segment(aes(x=c(shift_right, chrom_namesreg1$end[1:nrow(chrom_namesreg1)-1]), xend=c(shift_right, chrom_namesreg1$end[1:nrow(chrom_namesreg1)-1]), y=2.25, yend=8.25), linetype=3) +
  geom_segment(aes(x=c(chrom_nameschr3$start, chrom_nameschr3$end), xend=c(chrom_nameschr3$start, chrom_nameschr3$end), y=2.25, yend=8.25+yposadd), linetype=3) +
  geom_segment(aes(x=c(chrom_nameschr4A5$start[1], chrom_nameschr4A5$end[2]), xend=c(chrom_nameschr4A5$start[1], chrom_nameschr4A5$end[2]), y=2.25, yend=8.25+yposadd), linetype=3) +
  geom_segment(aes(x=chrom_nameschr4A5$start[2], xend=chrom_nameschr4A5$start[2], y=2.25+yposadd, yend=8.25+yposadd), linetype=3) +
  geom_segment(aes(x=c(chrom_namesreg2$end[1:nrow(chrom_namesreg2)-1]), xend=c(chrom_namesreg2$end[1:nrow(chrom_namesreg2)-1]), y=2.25, yend=8.25), linetype=3) +
  geom_segment(aes(x=c(chrom_nameschrZ$start, chrom_nameschrZ$end), xend=c(chrom_nameschrZ$start, chrom_nameschrZ$end), y=2.25, yend=8.25+yposadd), linetype=3) +
  # Vertical lines 3 substrata
  geom_segment(aes(x=c(strata3$start, strata3$end[3]), xend=c(strata3$start, strata3$end[3]), y=2.25+yposadd, yend=8.25+yposadd), linetype=3) +
  # Horisontal lines
  geom_segment(data=chrom_namesreg1, aes(x=min(start), xend=max(end), y=rep(2.5,nrow(chrom_namesreg1)), yend=rep(2.5,nrow(chrom_namesreg1)))) +
  geom_segment(data=chrom_namesreg1, aes(x=min(start), xend=max(end), y=rep(3.5,nrow(chrom_namesreg1)), yend=rep(3.5,nrow(chrom_namesreg1)))) +
  geom_segment(data=chrom_namesreg1, aes(x=min(start), xend=max(end), y=rep(4,nrow(chrom_namesreg1)), yend=rep(4,nrow(chrom_namesreg1)))) +
  geom_segment(data=chrom_namesreg1, aes(x=min(start), xend=max(end), y=rep(5,nrow(chrom_namesreg1)), yend=rep(5,nrow(chrom_namesreg1)))) +
  geom_segment(data=chrom_namesreg1, aes(x=min(start), xend=max(end), y=rep(5.5,nrow(chrom_namesreg1)), yend=rep(5.5,nrow(chrom_namesreg1)))) +
  geom_segment(data=chrom_namesreg1, aes(x=min(start), xend=max(end), y=rep(6.5,nrow(chrom_namesreg1)), yend=rep(6.5,nrow(chrom_namesreg1)))) +
  geom_segment(data=chrom_namesreg1, aes(x=min(start), xend=max(end), y=rep(7,nrow(chrom_namesreg1)), yend=rep(7,nrow(chrom_namesreg1)))) +
  geom_segment(data=chrom_namesreg1, aes(x=min(start), xend=max(end), y=rep(8,nrow(chrom_namesreg1)), yend=rep(8,nrow(chrom_namesreg1)))) +
  geom_segment(data=chrom_nameschr3, aes(x=min(start), xend=max(end), y=rep(2.5+yposadd,nrow(chrom_nameschr3)), yend=rep(2.5+yposadd,nrow(chrom_nameschr3)))) +
  geom_segment(data=chrom_nameschr3, aes(x=min(start), xend=max(end), y=rep(3.5+yposadd,nrow(chrom_nameschr3)), yend=rep(3.5+yposadd,nrow(chrom_nameschr3)))) +
  geom_segment(data=chrom_nameschr3, aes(x=min(start), xend=max(end), y=rep(4+yposadd,nrow(chrom_nameschr3)), yend=rep(4+yposadd,nrow(chrom_nameschr3)))) +
  geom_segment(data=chrom_nameschr3, aes(x=min(start), xend=max(end), y=rep(5+yposadd,nrow(chrom_nameschr3)), yend=rep(5+yposadd,nrow(chrom_nameschr3)))) +
  geom_segment(data=chrom_nameschr3, aes(x=min(start), xend=max(end), y=rep(5.5+yposadd,nrow(chrom_nameschr3)), yend=rep(5.5+yposadd,nrow(chrom_nameschr3)))) +
  geom_segment(data=chrom_nameschr3, aes(x=min(start), xend=max(end), y=rep(6.5+yposadd,nrow(chrom_nameschr3)), yend=rep(6.5+yposadd,nrow(chrom_nameschr3)))) +
  geom_segment(data=chrom_nameschr3, aes(x=min(start), xend=max(end), y=rep(7+yposadd,nrow(chrom_nameschr3)), yend=rep(7+yposadd,nrow(chrom_nameschr3)))) +
  geom_segment(data=chrom_nameschr3, aes(x=min(start), xend=max(end), y=rep(8+yposadd,nrow(chrom_nameschr3)), yend=rep(8+yposadd,nrow(chrom_nameschr3)))) +
  geom_segment(data=chrom_nameschr4, aes(x=min(start), xend=max(end), y=rep(2.5,nrow(chrom_nameschr4)), yend=rep(2.5,nrow(chrom_nameschr4)))) +
  geom_segment(data=chrom_nameschr4, aes(x=min(start), xend=max(end), y=rep(3.5,nrow(chrom_nameschr4)), yend=rep(3.5,nrow(chrom_nameschr4)))) +
  geom_segment(data=chrom_nameschr4, aes(x=min(start), xend=max(end), y=rep(4,nrow(chrom_nameschr4)), yend=rep(4,nrow(chrom_nameschr4)))) +
  geom_segment(data=chrom_nameschr4, aes(x=min(start), xend=max(end), y=rep(5,nrow(chrom_nameschr4)), yend=rep(5,nrow(chrom_nameschr4)))) +
  geom_segment(data=chrom_nameschr4, aes(x=min(start), xend=max(end), y=rep(5.5,nrow(chrom_nameschr4)), yend=rep(5.5,nrow(chrom_nameschr4)))) +
  geom_segment(data=chrom_nameschr4, aes(x=min(start), xend=max(end), y=rep(6.5,nrow(chrom_nameschr4)), yend=rep(6.5,nrow(chrom_nameschr4)))) +
  geom_segment(data=chrom_nameschr4, aes(x=min(start), xend=max(end), y=rep(7,nrow(chrom_nameschr4)), yend=rep(7,nrow(chrom_nameschr4)))) +
  geom_segment(data=chrom_nameschr4, aes(x=min(start), xend=max(end), y=rep(8,nrow(chrom_nameschr4)), yend=rep(8,nrow(chrom_nameschr4)))) +
  geom_segment(data=chrom_nameschr4A5, aes(x=min(start), xend=max(end), y=rep(2.5+yposadd,nrow(chrom_nameschr4A5)), yend=rep(2.5+yposadd,nrow(chrom_nameschr4A5)))) +
  geom_segment(data=chrom_nameschr4A5, aes(x=min(start), xend=max(end), y=rep(3.5+yposadd,nrow(chrom_nameschr4A5)), yend=rep(3.5+yposadd,nrow(chrom_nameschr4A5)))) +
  geom_segment(data=chrom_nameschr4A5, aes(x=min(start), xend=max(end), y=rep(4+yposadd,nrow(chrom_nameschr4A5)), yend=rep(4+yposadd,nrow(chrom_nameschr4A5)))) +
  geom_segment(data=chrom_nameschr4A5, aes(x=min(start), xend=max(end), y=rep(5+yposadd,nrow(chrom_nameschr4A5)), yend=rep(5+yposadd,nrow(chrom_nameschr4A5)))) +
  geom_segment(data=chrom_nameschr4A5, aes(x=min(start), xend=max(end), y=rep(5.5+yposadd,nrow(chrom_nameschr4A5)), yend=rep(5.5+yposadd,nrow(chrom_nameschr4A5)))) +
  geom_segment(data=chrom_nameschr4A5, aes(x=min(start), xend=max(end), y=rep(6.5+yposadd,nrow(chrom_nameschr4A5)), yend=rep(6.5+yposadd,nrow(chrom_nameschr4A5)))) +
  geom_segment(data=chrom_nameschr4A5, aes(x=min(start), xend=max(end), y=rep(7+yposadd,nrow(chrom_nameschr4A5)), yend=rep(7+yposadd,nrow(chrom_nameschr4A5)))) +
  geom_segment(data=chrom_nameschr4A5, aes(x=min(start), xend=max(end), y=rep(8+yposadd,nrow(chrom_nameschr4A5)), yend=rep(8+yposadd,nrow(chrom_nameschr4A5)))) +
  geom_segment(data=chrom_namesreg2, aes(x=min(start), xend=max(end), y=rep(2.5,nrow(chrom_namesreg2)), yend=rep(2.5,nrow(chrom_namesreg2)))) +
  geom_segment(data=chrom_namesreg2, aes(x=min(start), xend=max(end), y=rep(3.5,nrow(chrom_namesreg2)), yend=rep(3.5,nrow(chrom_namesreg2)))) +
  geom_segment(data=chrom_namesreg2, aes(x=min(start), xend=max(end), y=rep(4,nrow(chrom_namesreg2)), yend=rep(4,nrow(chrom_namesreg2)))) +
  geom_segment(data=chrom_namesreg2, aes(x=min(start), xend=max(end), y=rep(5,nrow(chrom_namesreg2)), yend=rep(5,nrow(chrom_namesreg2)))) +
  geom_segment(data=chrom_namesreg2, aes(x=min(start), xend=max(end), y=rep(5.5,nrow(chrom_namesreg2)), yend=rep(5.5,nrow(chrom_namesreg2)))) +
  geom_segment(data=chrom_namesreg2, aes(x=min(start), xend=max(end), y=rep(6.5,nrow(chrom_namesreg2)), yend=rep(6.5,nrow(chrom_namesreg2)))) +
  geom_segment(data=chrom_namesreg2, aes(x=min(start), xend=max(end), y=rep(7,nrow(chrom_namesreg2)), yend=rep(7,nrow(chrom_namesreg2)))) +
  geom_segment(data=chrom_namesreg2, aes(x=min(start), xend=max(end), y=rep(8,nrow(chrom_namesreg2)), yend=rep(8,nrow(chrom_namesreg2)))) +
  geom_segment(data=chrom_nameschrZ, aes(x=min(start), xend=max(end), y=rep(2.5+yposadd,nrow(chrom_nameschrZ)), yend=rep(2.5+yposadd,nrow(chrom_nameschrZ)))) +
  geom_segment(data=chrom_nameschrZ, aes(x=min(start), xend=max(end), y=rep(3.5+yposadd,nrow(chrom_nameschrZ)), yend=rep(3.5+yposadd,nrow(chrom_nameschrZ)))) +
  geom_segment(data=chrom_nameschrZ, aes(x=min(start), xend=max(end), y=rep(4+yposadd,nrow(chrom_nameschrZ)), yend=rep(4+yposadd,nrow(chrom_nameschrZ)))) +
  geom_segment(data=chrom_nameschrZ, aes(x=min(start), xend=max(end), y=rep(5+yposadd,nrow(chrom_nameschrZ)), yend=rep(5+yposadd,nrow(chrom_nameschrZ)))) +
  geom_segment(data=chrom_nameschrZ, aes(x=min(start), xend=max(end), y=rep(5.5+yposadd,nrow(chrom_nameschrZ)), yend=rep(5.5+yposadd,nrow(chrom_nameschrZ)))) +
  geom_segment(data=chrom_nameschrZ, aes(x=min(start), xend=max(end), y=rep(6.5+yposadd,nrow(chrom_nameschrZ)), yend=rep(6.5+yposadd,nrow(chrom_nameschrZ)))) +
  geom_segment(data=chrom_nameschrZ, aes(x=min(start), xend=max(end), y=rep(7+yposadd,nrow(chrom_nameschrZ)), yend=rep(7+yposadd,nrow(chrom_nameschrZ)))) +
  geom_segment(data=chrom_nameschrZ, aes(x=min(start), xend=max(end), y=rep(8+yposadd,nrow(chrom_nameschrZ)), yend=rep(8+yposadd,nrow(chrom_nameschrZ)))) +
  # Other
  coord_polar() +
  xlim(0, max(chrom_names$end) + shift_right) +
  ylim(0, 8.75+yposadd) +
  theme_void() +
  theme(
    legend.position.inside = c(0,0),
    legend.title = element_text(size=30),
    legend.text = element_text(size=20),
    legend.key.size = unit(1, units = "cm"),
    legend.spacing.y = unit(1, "cm"),
    axis.title.y = element_blank(),
    axis.text.x = element_blank())

jpeg("Figures/genome_summary.jpg", width=6000, height=6000, res=300)
plot
dev.off()


# Hypothetical sex chrom structure
struct <- as.data.frame(matrix(NA, 10, 5))
colnames(struct) <- c("category", "size", "label", "label_height", "label_pos")
tot_length <- 15.4 + 9.2 + 195.3 + 10.1 + 24.1
struct[1,] <- c("prop_unk_PAR5", 15.4, "PAR 5", 2, 15.4/2)
struct[2,] <- c("prop_PAR5", 9.2, "PAR 5", 2, 15.4 + 9.2/2)
struct[3,] <- c("prop_5", 36.3, "5", 2, 15.4 + 9.2 + 36.3/2)
struct[4,] <- c("prop_ancZ", 72.9, "Ancestral (S0, S1, S2, S3)", 2, 15.4 + 9.2 + 36.3 + 72.9/2)
struct[5,] <- c("prop_4A", 9.6, "4A", 2, 15.4 + 9.2 + 36.3 + 72.9 + 9.6/2)
struct[6,] <- c("prop_3a", 8.0, "3-a", 2, 15.4 + 9.2 + 36.3 + 72.9 + 9.6 + 8.0/2)
struct[7,] <- c("prop_3b", 3.6, "3-b", 2, 15.4 + 9.2 + 36.3 + 72.9 + 9.6 + 8.0 + 3.6/2)
struct[8,] <- c("prop_3c", 64.9, "3-c", 2, 15.4 + 9.2 + 36.3 + 72.9 + 9.6 + 8.0 + 3.6 + 64.9/2)
struct[9,] <- c("prop_PAR3", 10.1, "PAR 3", 2, 15.4 + 9.2 + 36.3 + 72.9 + 9.6 + 8.0 + 3.6 + 64.9 + 10.1/2)
struct[10,] <- c("prop_unk_PAR3", 24.1, "PAR 3", 2, 15.4 + 9.2 + 36.3 + 72.9 + 9.6 + 8.0 + 3.6 + 64.9 + 10.1 + 24.1/2)

struct$category <- factor(struct$category, order=T, levels = struct$category)  # Ensure factor levels
struct$size <- as.numeric(struct$size)
struct$label_height <- as.numeric(0.5)
struct$label_height[5] <- struct$label_height[6] - 1.1
struct$label_height[7] <- struct$label_height[6] - 1.1
struct$tick <- struct$label_height + 0.15
struct$label_pos <- as.numeric(struct$label_pos)/tot_length
struct$value <- struct$size/tot_length
struct$label2 <- paste(as.character(struct$size), "Mb", sep=" ")
struct$label2[6] <- "8.0 Mb"
struct$label2[1] <- "?"
struct$label2[10] <- "?"
struct$age <- c(NA, NA, "10.2 MY", "40.3 - 135.3 MY", "22 MY", "22 MY", "19.7 MY", "6.3 MY", NA, NA)
struct$age2 <- as.numeric( c(0, 0, 10.2, 87.8, 22, 22, 19.7, 6.3, 0, 0))

struct$value_norm <- rescale(struct$age2, to = c(0,1))
grad_pal <- gradient_n_pal(c("white", "gray30"))
struct$color_hex <- grad_pal(struct$value_norm)
#struct$color_hex[4)] <- "black"
color_vec <- setNames(struct$color_hex, struct$category)
pattern_vec <- setNames(rep("none", length(color_vec)), struct$category)
pattern_vec[c(1,10)] <- "stripe"

chrom_plot <- ggplot(struct, aes(x = 1, y = value, fill = category, pattern = category)) +
    geom_bar_pattern(stat = "identity", width = 0.5, color = "black", position = position_stack(reverse = TRUE),  pattern_fill="black", pattern_angle=45, pattern_density=0.05, pattern_spacing=0.05) +
    scale_fill_manual(values = color_vec) +
    scale_pattern_manual(values = pattern_vec) +
    geom_segment(aes(x=tick, xend=0.75, y=label_pos, yend=label_pos)) +
    annotate(y=struct$label_pos, label=struct$label, x=struct$label_height, size=6, geom="text") +
    annotate(y=struct$label_pos, label=struct$label2, x=struct$label_height-0.35, size=6, geom="text") +
    annotate(y=struct$label_pos, label=struct$age, x=struct$label_height-0.70, size=6, geom="text") +
    coord_flip(clip="off") +
    theme_void() +
    theme(legend.position = "none", 
          plot.margin = unit(c(1, 1, 1, 1), "lines"))

chrom_plot

jpeg("Figures/chrom_plot.jpeg", width=5000, height=750, res=300)
chrom_plot
dev.off()
