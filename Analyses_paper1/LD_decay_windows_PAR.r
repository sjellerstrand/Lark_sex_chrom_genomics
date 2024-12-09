#!/usr/bin/Rscript

## Export variables and load libraries
rm(list=ls())
options(scipen=999)
library(data.table)

args <- strsplit(commandArgs(trailingOnly=T), "=", fixed=T)
for(i in 1:length(args)) {
  assign(args[[i]][1], args[[i]][2])
}

SCAFFOLD_LENGTH <- as.numeric(SCAFFOLD_LENGTH)
WINDOW <- as.numeric(WINDOW)
STEP <- as.numeric(STEP)

data <- read.table(gzfile(paste(OUTDIR, "/", PROJECT, "_", SCAFFOLD, "_", DATA, ".LDdecay.LD.gz", sep="")), head=T)
colnames(data) <- c("scaffold", "snp1", "snp2","Dprime", "LOD","r2","CIlow","CIhi", "Dist")
snps <- unique(c(data$snp1, data$snp2))

# Create windows
startpos <- seq(1, SCAFFOLD_LENGTH - WINDOW, STEP)
endpos <- startpos + WINDOW - 1
if(endpos[length(endpos)] < SCAFFOLD_LENGTH) {
  startpos <- c(startpos, SCAFFOLD_LENGTH - WINDOW + 1)
  endpos <- c(endpos, SCAFFOLD_LENGTH)
}

# Store data
data2 <- as.data.frame(matrix(NA, length(startpos), 9))
colnames(data2) <- c("scaffold", "start", "end", "mid", "window_size", "window_N_snps", "r2", "CIlow","CIhi")
data2$scaffold <- SCAFFOLD

# Loop through windows
for (i in 1:nrow(data2)) {
  data2$start[i] <- startpos[i]
  data2$end[i] <- endpos[i]
  data2$mid[i] <- startpos[i] + (endpos[i] - startpos[i] + 1)/2
  data2$window_size[i] <- endpos[i] - startpos[i] + 1
  data2$window_N_snps[i] <-  length(which(between(snps, startpos[i], endpos[i])))
  data2$r2[i] <- mean(data$r2[which(between(data$snp1, startpos[i], endpos[i]) | between(data$snp2, startpos[i], endpos[i]))])
  data2$CIlow[i] <- mean(data$CIlow[which(between(data$snp1, startpos[i], endpos[i]) | between(data$snp2, startpos[i], endpos[i]))])
  data2$CIhi[i] <- mean(data$CIhi[which(between(data$snp1, startpos[i], endpos[i]) | between(data$snp2, startpos[i], endpos[i]))])
}

### Write file
write.table(data2, file=paste(OUTDIR, "/window_", WINDOW, "_step_", STEP, "/", PROJECT, "_", SCAFFOLD, "_", DATA, "_window_", WINDOW, "_step_", STEP, ".txt", sep=""), quote=FALSE, sep="\t", row.names = F, col.names = T)

quit()
