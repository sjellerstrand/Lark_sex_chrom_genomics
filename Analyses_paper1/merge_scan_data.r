#!/usr/bin/Rscript

## Export variables and load libraries
rm(list=ls())
options(scipen=999)

args <- strsplit(commandArgs(trailingOnly=T), "=", fixed=T)
for(i in 1:length(args)) {
  assign(args[[i]][1], args[[i]][2])
}

WINDOW <- as.numeric(WINDOW)
STEP <- as.numeric(STEP)


### Read in genome scan data
scan_data1 <- read.delim(paste(OUTDIR2, "/", PROJECT, "_", DATA, ".pi_tajimas_D", sep=""), sep=",", head=T)
if(DATA != "heterogametic") {
  scan_data2 <- read.delim(paste(OUTDIR2, "/", PROJECT, "_", DATA, ".pairwise_sex", sep=""), sep=",", head=T)
} else if(DATA == "heterogametic") {
  scan_data2 <- cbind(scan_data1[,1:5], matrix(NA, nrow(scan_data1), 4))
  colnames(scan_data2)[6:9] <- c("pi_Female", "pi_Male", "dxy_Female_Male", "Fst_Female_Male")
}
window_data <- read.delim(paste(OUTDIR1, "/windows/", DATA, "_windows_", WINDOW, "_steps_", STEP, ".txt", sep=""), sep="\t", head=T)

### Merge and adjust data
scan_data3 <- merge(scan_data1, scan_data2, by=c("scaffold", "start", "end", "mid", "sites"), sort=F)
scan_data <- merge(scan_data3, window_data, by=c("scaffold", "start", "end"), sort=F)
scan_data <- subset(scan_data, select=c(-l_all, -S_all))
scan_data <- cbind(scan_data, matrix(NA, nrow(scan_data), 4))
colnames(scan_data)[(ncol(scan_data)-3):ncol(scan_data)] <- c("pi_abs", "pi_Female_abs", "pi_Male_abs", "dxy_Female_Male_abs")
scan_data$pi_abs <- scan_data$thetaPi_all/scan_data$N_tot_sites
scan_data$pi_Female_abs <- (scan_data$pi_Female * scan_data$sites)/scan_data$N_tot_sites
scan_data$pi_Male_abs <- (scan_data$pi_Male * scan_data$sites)/scan_data$N_tot_sites
scan_data$dxy_Female_Male_abs <- (scan_data$dxy_Female_Male * scan_data$sites)/scan_data$N_tot_sites

### Set NAs to 0, correct mid coordinate, and add hetgam specific pi
no_snps <- which(scan_data$sites == 0 & (scan_data$N_callable_sites > 0 | !is.na(scan_data$N_callable_sites)))
if(DATA != "heterogametic") {
  scan_data[no_snps,c(6,7,8,9,10,11,12,13,17,18,19,20)] <- 0
} else if(DATA == "heterogametic") {
  scan_data[no_snps,c(6,7,8,9,17)] <- 0
  if(HETGAM == "Female") {
    scan_data$pi_Female <- scan_data$pi_all
    scan_data$pi_Female_abs <- scan_data$pi_abs
  } else if(HETGAM == "Male") {
    scan_data$pi_Male <- scan_data$pi_all
    scan_data$pi_Male_abs <- scan_data$pi_abs
  }
}
scan_data$mid <- round(scan_data$start + ((scan_data$end - scan_data$start)/2))

### Add data type
scan_data$Genome <- DATA

### Write data
write.table(scan_data, file=paste(OUTDIR2, "/", PROJECT, "_", DATA, "_windows_", WINDOW, "_steps_", STEP, "_data.txt", sep=""), quote=FALSE, sep='\t', row.names = F, col.names = T)

quit()

