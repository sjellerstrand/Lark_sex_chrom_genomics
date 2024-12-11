#!/usr/bin/Rscript

## Export variables and load libraries
rm(list=ls())

args <- strsplit(commandArgs(trailingOnly=T), "=", fixed=T)
for(i in 1:length(args)) {
  assign(args[[i]][1], args[[i]][2])
}

### Read in genome scan data
scan_data <- read.delim(paste(OUTDIR2, "/", PROJECT, "_", DATA, ".pi_tajimas_D", sep=""), sep=",", head=T)
window_data <- read.delim(paste(OUTDIR1, "/windows/", DATA, "_windows_", WINDOW, "_steps_", STEP, "_exon_dist_", EXON_DIST, ".txt", sep=""), sep="\t", head=T)

### Merge and adjust data
scan_data <- merge(scan_data, window_data, by=c("scaffold", "start", "end"), sort=F)
scan_data <- subset(scan_data, select=c(-l_all, -S_all))
scan_data <- cbind(scan_data, rep(NA, nrow(scan_data)), rep(NA, nrow(scan_data)), rep(NA, nrow(scan_data)), rep(NA, nrow(scan_data)))
colnames(scan_data)[(ncol(scan_data)-3):ncol(scan_data)] <- c("pi_abs", "data_type", "exon_distance", "Project")
scan_data$pi_abs <- scan_data$thetaPi_all/scan_data$N_callable_sites
scan_data$data_type <- DATA
scan_data$exon_distance <- EXON_DIST
scan_data$Project <- PROJECT
scan_data$mid <- round(scan_data$start + ((scan_data$end - scan_data$start)/2))

### Write data
write.table(scan_data, file=paste(OUTDIR2, "/", PROJECT, "_", DATA, "_windows_", WINDOW, "_steps_", STEP, "_exon_dist_", EXON_DIST, "_data.txt", sep=""), quote=FALSE, sep='\t', row.names = F, col.names = T)

quit()

