## Export variables and load libraries
rm(list=ls())
library(tidyverse)
library(patchwork)
library(scales)
library(viridis)
library(sf)
library(maps)
sessionInfo()

options(scipen=999)
setwd("C:/Users/Simon JE//OneDrive - Lund University/Dokument/Simon/PhD/Projects/Rasolark_2021/Results")

# Make map
world <- map_data("world")
layers <- st_read("Map/Alauda.shp")
species_list <- c("Alauda arvensis", "Alauda razae")
layers <- layers[layers$sci_name %in% species_list, ]

layers$sci_name <- factor(layers$sci_name , order=T, label=c("Skylark", "Raso lark"), level=c("Alauda arvensis", "Alauda razae"))

samples <- data.frame(
  longitude = c(-24.588, 5.609, 11.410),
  latitude  = c(16.618, 52.317, 43.255),
  label = c(
    "Raso~lark",
    "Skylark~\"'\"*italic(arvensis)*\"'\"",
    "Skylark~\"'\"*italic(cantarella)*\"'\""),
  species   = c("Raso lark", "Skylark", "Skylark")
)
samples_sf <- st_as_sf(samples, coords = c("longitude", "latitude"), crs = 4326)

plot_Map <- ggplot() +
  geom_polygon(data = world, aes(x=long, y = lat, group = group), fill = "grey95", color = "grey70", linewidth = 0.2) + 
  geom_sf(data = layers, aes(fill = sci_name), color = NA, alpha = 0.5, show.legend = FALSE) + # range
  geom_sf(data = samples_sf, aes(fill=species),  shape = 21, color = "black", size = 5, alpha=0.5, stroke = 1.5,) + # points
  scale_fill_manual(name = "Species", values=c("Raso lark" = "#92c5de", "Skylark" = "#f4a582")) +
  coord_sf(xlim = c(-25, 180), ylim = c(7, 85)) + # crop
  theme_void() + 
  theme(
    legend.position = "right",
    legend.key.size = unit(1, 'cm'),
    legend.title = element_text(size=15),
    legend.text = element_text(size=12))

png("Figures/Map.png", width=4000, height=2000, res=300)
plot_Map
dev.off()


mu <- 7.16*10^-9
gen_R <- 3.4143629
gen_S <- 2.833403383

min_rec <- 1960
max_rec <- 2020
min_int <- 1400
max_int <- 2020

sample_year <- 2017 - gen_R/2

# Sample info
sample_raso <- read.delim("Figures/Rasolark_2021_sample_info.txt", sep="\t", head=T)

# Census counts
census <- read.delim("demography/popsize_raso.txt", sep="\t", head=T)

census$Age <- 2018-census$Year

# Gone
gone_A <- read.delim("demography/gone/autosomal/Output_Ne_Rasolark_2021_gone_data", sep="\t", head=T, skip=1)
gone_A$Region <- "Autosomal"
gone_A$Age <- gone_A$Generation*gen_R
gone_A$Year <- 2018-gone_A$Age-gen_R/2

gone_Z <- read.delim("demography/gone/homogametic/Output_Ne_Rasolark_2021_gone_data", sep="\t", head=T, skip=1)
gone_Z$Region <- "Z"
gone_Z$Age <- gone_Z$Generation*gen_R/0.75
gone_Z$Year <- 2018-gone_Z$Age-(gen_R/0.75)/2


gone_A_rep <- cbind(rep(NA, nrow(gone_A)))
gone_Z_rep <- cbind(rep(NA, nrow(gone_Z)))
for(i in 1:1000) {
  gone_A_rep <- cbind(gone_A_rep, read.delim(paste("demography/gone/autosomal/TEMPORARY_FILES/outfileLD_TEMP/outfileLD_", i, "_GONE_Nebest", sep=""), sep="\t", head=F)[,2])
  gone_Z_rep <- cbind(gone_Z_rep, read.delim(paste("demography/gone/homogametic/TEMPORARY_FILES/outfileLD_TEMP/outfileLD_", i, "_GONE_Nebest", sep=""), sep="\t", head=F)[,2])
}
gone_A_rep <- gone_A_rep[,-1]
gone_Z_rep <- gone_Z_rep[,-1]

gone_A$Mean <- rep(NA, nrow(gone_A))
gone_A$Upper95 <- rep(NA, nrow(gone_A))
gone_A$Lower95 <- rep(NA, nrow(gone_A))
gone_Z$Mean <- rep(NA, nrow(gone_Z))
gone_Z$Upper95 <- rep(NA, nrow(gone_Z))
gone_Z$Lower95 <- rep(NA, nrow(gone_Z))

for(i in 1:nrow(gone_A)) {
  gone_A$Mean[i] <- mean(gone_A_rep[i,])
  gone_A$Upper95[i] <- quantile(gone_A_rep[i,], prob=0.975)
  gone_A$Lower95[i] <- quantile(gone_A_rep[i,], prob=0.025)
  gone_Z$Mean[i] <- mean(gone_Z_rep[i,])
  gone_Z$Upper95[i] <- quantile(gone_Z_rep[i,], prob=0.975)
  gone_Z$Lower95[i] <- quantile(gone_Z_rep[i,], prob=0.025)
}
gone <- rbind(gone_A, gone_Z)

# output for slim
raso <- gone[which(gone$Region=="Autosomal"),]
raso$Age_rev <- 1000000 - raso$Age
raso$gen <- ((raso$Age-1)/gen_R)
raso$gen_rev <- ((raso$Age_rev-1)/gen_R)

write.table(raso, "raso_demography_msmc.tsv", quote=F, sep="\t", row.names = F)


# msmc
msmc_data <- read.table("C:/Users/Simon JE/OneDrive - Lund University/Dokument/Simon/PhD/Projects/Rasolark_2021/Results/demography/msmc/Rasolark/autosomal/msmc_out.final.txt", head=T)
msmc_RA <- as.data.frame(matrix(NA, nrow(msmc_data)*2, 2))
colnames(msmc_RA) <- c("Age", "Ne")
for(i in 1:nrow(msmc_data)) {
  msmc_RA$Age[i*2-1] <- msmc_data$left_time_boundary[i]
  msmc_RA$Age[i*2] <- msmc_data$right_time_boundary[i]
  msmc_RA$Ne[(i*2-1):(i*2)] <- msmc_data$lambda[i]
}

msmc_RA <- msmc_RA[-nrow(msmc_RA),]
msmc_RA$Age <- (msmc_RA$Age/mu)*gen_R + 1
msmc_RA$Ne <- (1/(msmc_RA$Ne)/(2*mu))
msmc_RA$Region <- "Autosomal"
msmc_RA$Species <- "Raso lark"


msmc_data <- read.table("C:/Users/Simon JE/OneDrive - Lund University/Dokument/Simon/PhD/Projects/Rasolark_2021/Results/demography/msmc/Rasolark/homogametic/msmc_out.final.txt", head=T)
msmc_RZ <- as.data.frame(matrix(NA, nrow(msmc_data)*2, 2))
colnames(msmc_RZ) <- c("Age", "Ne")
for(i in 1:nrow(msmc_data)) {
  msmc_RZ$Age[i*2-1] <- msmc_data$left_time_boundary[i]
  msmc_RZ$Age[i*2] <- msmc_data$right_time_boundary[i]
  msmc_RZ$Ne[(i*2-1):(i*2)] <- msmc_data$lambda[i]
}

msmc_RZ <- msmc_RZ[-nrow(msmc_RZ),]
msmc_RZ$Age <- (msmc_RZ$Age/mu)*(gen_R/0.75) + 1
msmc_RZ$Ne <- (1/(msmc_RZ$Ne)/(2*mu))
msmc_RZ$Region <- "Z"
msmc_RZ$Species <- "Raso lark"

msmc_data <- read.table("C:/Users/Simon JE/OneDrive - Lund University/Dokument/Simon/PhD/Projects/Rasolark_2021/Results/demography/msmc/Skylark/autosomal/msmc_out.final.txt", head=T)
msmc_SA <- as.data.frame(matrix(NA, nrow(msmc_data)*2, 2))
colnames(msmc_SA) <- c("Age", "Ne")
for(i in 1:nrow(msmc_data)) {
  msmc_SA$Age[i*2-1] <- msmc_data$left_time_boundary[i]
  msmc_SA$Age[i*2] <- msmc_data$right_time_boundary[i]
  msmc_SA$Ne[(i*2-1):(i*2)] <- msmc_data$lambda[i]
}

msmc_SA <- msmc_SA[-nrow(msmc_SA),]
msmc_SA$Age <- (msmc_SA$Age/mu)*gen_S + 1
msmc_SA$Ne <- (1/(msmc_SA$Ne)/(2*mu))
msmc_SA$Region <- "Autosomal"
msmc_SA$Species <- "Skylark"


msmc_data <- read.table("C:/Users/Simon JE/OneDrive - Lund University/Dokument/Simon/PhD/Projects/Rasolark_2021/Results/demography/msmc/Skylark/homogametic/msmc_out.final.txt", head=T)
msmc_SZ <- as.data.frame(matrix(NA, nrow(msmc_data)*2, 2))
colnames(msmc_SZ) <- c("Age", "Ne")
for(i in 1:nrow(msmc_data)) {
  msmc_SZ$Age[i*2-1] <- msmc_data$left_time_boundary[i]
  msmc_SZ$Age[i*2] <- msmc_data$right_time_boundary[i]
  msmc_SZ$Ne[(i*2-1):(i*2)] <- msmc_data$lambda[i]
}

msmc_SZ <- msmc_SZ[-nrow(msmc_SZ),]
msmc_SZ$Age <- (msmc_SZ$Age/mu)*(gen_S/0.75) + 1
msmc_SZ$Ne <- (1/(msmc_SZ$Ne)/(2*mu))
msmc_SZ$Region <- "Z"
msmc_SZ$Species <- "Skylark"


# Import bootstraps
msmc_rep <- as.data.frame(rbind(rep(NA, 7)))
colnames(msmc_rep) <- c("Age", "Ne", "Region", "Species", "Mean", "Upper95", "Lower95")
for(i in c("Rasolark", "Skylark")) {
  if(i == "Rasolark") {
    gen <- gen_R
  } else if(i == "Skylark") {
    gen <- gen_S
  }
  for(j in c("autosomal", "homogametic")) {
    if(j == "autosomal") {
      if(i == "Rasolark") {
        msmc_data1 <- msmc_RA
      } else if(i == "Skylark") {
        msmc_data1 <- msmc_SA
      }
    } else if(j == "homogametic") {
      gen <- gen / 0.75
      if(i == "Rasolark") {
        msmc_data1 <- msmc_RZ
      } else if(i == "Skylark") {
        msmc_data1 <- msmc_SZ
      }
    }
    msmc_data_boot <- as.data.frame(matrix(NA, nrow(msmc_data1), 100))
    x <- 0
    for(k in 1:100) {
      x <- x+1
      msmc_data2 <- read.delim(paste("demography/msmc/", i, "/", j, "/Final_bootstraps/bootstraps_", k, ".final.txt", sep=""), sep="\t", head=T)
      boot <- rep(NA, nrow(msmc_data2)*2)
      for(l in 1:nrow(msmc_data2)) {
        boot[(l*2-1):(l*2)] <- msmc_data2$lambda[l]
      }
      boot <- boot[-length(boot)]
      msmc_data_boot[,k] <- (1/(boot)/(2*mu))
    }
    msmc_data_boot$Upper95 <- rep(NA, nrow(msmc_data_boot))
    msmc_data_boot$Lower95 <- rep(NA, nrow(msmc_data_boot))
    for(k in 1:nrow(msmc_data_boot)) {
      msmc_data_boot$Mean[k] <- mean(as.numeric(msmc_data_boot[k,1:100]))
      msmc_data_boot$Upper95[k] <- quantile(as.numeric(msmc_data_boot[k,1:100]), prob=0.975)
      msmc_data_boot$Lower95[k] <- quantile(as.numeric(msmc_data_boot[k,1:100]), prob=0.025)
    }
    msmc_rep <- rbind(msmc_rep, cbind(msmc_data1, "Mean" = msmc_data_boot$Mean, "Upper95" = msmc_data_boot$Upper95, "Lower95" = msmc_data_boot$Lower95))
  }
}

msmc <- msmc_rep[-1,]
msmc$Species <- factor(msmc$Species, order=T, levels=c("Raso lark", "Skylark"))
msmc$Region <- factor(msmc$Region, order=T, levels=c("Autosomal", "Z"))

# output for slim
raso <- msmc[which(msmc$Region=="Autosomal" & msmc$Species=="Raso lark"),]
raso$Age_rev <- 500000 - raso$Age
raso$gen <- ((raso$Age-1)/gen_R)
raso$gen_rev <- ((raso$Age_rev-1)/gen_R)

write.table(raso, "raso_demography_msmc.tsv", quote=F, sep="\t", row.names = F)

# linkage decay
### Import data
data_RA <- read.table(paste(getwd(), "/linkage_decay/Rasolark/Autosomal/Rasolark_2021_autosomal.LDdecay.stat.bins.tsv", sep=""), sep="\t", head=T)
data_RZ <- read.table(paste(getwd(), "/linkage_decay/Rasolark/Homogametic/Rasolark_2021_homogametic.LDdecay.stat.bins.tsv", sep=""), sep="\t", head=T)
data_SA <- read.table(paste(getwd(), "/linkage_decay/Skylark/Autosomal/Skylark_2021_autosomal.LDdecay.stat.bins.tsv", sep=""), sep="\t", head=T)
data_SZ <- read.table(paste(getwd(), "/linkage_decay/Skylark/Homogametic/Skylark_2021_homogametic.LDdecay.stat.bins.tsv", sep=""), sep="\t", head=T)

# Fit model
LB_RA <- median(data_RA$Mean_r2)
fit_RA <- nls(Mean_r2 ~ a * exp(-Mid / b) + LB_RA, start = list(a = max(data_RA$Mean_r2), b = 100000), data = data_RA)
data_RA$pred <- predict(fit_RA, data_RA)
LB_RZ <- median(data_RZ$Mean_r2)
fit_RZ <- nls(Mean_r2 ~ a * exp(-Mid / b) + LB_RZ, start = list(a = max(data_RZ$Mean_r2), b = 100000), data = data_RZ)
data_RZ$pred <- predict(fit_RZ, data_RZ)

LB_SA <- median(data_SA$Mean_r2)
fit_SA <- nls(Mean_r2 ~ a * exp(-Mid / b) + LB_SA, start = list(a = max(data_SA$Mean_r2), b = 100000), data = data_SA)
data_SA$pred <- predict(fit_SA, data_SA)
LB_SZ <- median(data_SZ$Mean_r2)
fit_SZ <- nls(Mean_r2 ~ a * exp(-Mid / b) + LB_SZ, start = list(a = max(data_SZ$Mean_r2), b = 100000), data = data_SZ)
data_SZ$pred <- predict(fit_SZ, data_SZ)

data_RA$Species <- "Raso lark"
data_RA$Region <- "Autosomal"
data_RZ$Species <- "Raso lark"
data_RZ$Region <- "Z"
data_SA$Species <- "Skylark"
data_SA$Region <- "Autosomal"
data_SZ$Species <- "Skylark"
data_SZ$Region <- "Z"

data_LD <- rbind(data_RA, data_RZ, data_SA, data_SZ)
data_LD$group <- paste(data_LD$Species, data_LD$Region, sep=" ")

### Runs of Homozygosity
ROH_A <- read.delim("ROH/autosomal/Rasolark_2021_autosomal_roh", sep="\t", head=T)
ROH_A$Region <- "Autosomal"
ROH_Z <- read.delim("ROH/homogametic//Rasolark_2021_homogametic_roh", sep="\t", head=T)
ROH_Z$Region <- "Z"
ROH <- rbind(ROH_A, ROH_Z)

ROH$Age <- rep(NA, nrow(ROH))
for(i in 1:nrow(sample_raso)) {
  ROH$Age[which(ROH$Sample == sample_raso$ID[i] & ROH$Region == "Autosomal")] <- 2018 + gen_R/2 - sample_raso$Year[i] + ROH$Age_generations[which(ROH$Sample == sample_raso$ID[i] & ROH$Region == "Autosomal")] * gen_R
  ROH$Age[which(ROH$Sample == sample_raso$ID[i] & ROH$Region == "Z")] <- 2018 + (gen_R/0.75)/2 - sample_raso$Year[i] + (ROH$Age_generations[which(ROH$Sample == sample_raso$ID[i] & ROH$Region == "Z")] * gen_R / 0.75)
}

ROH$Age[which(is.na(ROH$Age) & ROH$Region == "Z")] <- ROH$Age_generations[which(is.na(ROH$Age) & ROH$Region == "Z")] * gen_R / 0.75
ROH$Year<- 2018 - ROH$Age

ROH$Age_class <- rep(NA, nrow(ROH))
ROH$Year_class <- rep(NA, nrow(ROH))
ROH$Age_class[which(ROH$Age_generations < 2)] <- "< 2"
ROH$Year_class[which(ROH$Age_generations < 2 & ROH$Region == "Autosomal")] <- paste("<", signif(2*gen_R, 2), sep=" ")
ROH$Year_class[which(ROH$Age_generations < 2 & ROH$Region == "Z")] <- paste("<", signif(2*gen_R/0.75, 2), sep=" ")
ROH$Age_class[which(ROH$Age_generations >= 2 & ROH$Age_generations < 4)] <- "2 - 4"
ROH$Year_class[which(ROH$Age_generations >= 2 & ROH$Age_generations < 4 & ROH$Region == "Autosomal")] <- paste(signif(2*gen_R, 2), "-", signif(4*gen_R, 2), sep=" ")
ROH$Year_class[which(ROH$Age_generations >= 2 & ROH$Age_generations < 4 & ROH$Region == "Z")] <- paste(signif(2*gen_R/0.75, 2), "-", signif(4*gen_R/0.75, 2), sep=" ")
ROH$Age_class[which(ROH$Age_generations >= 4 & ROH$Age_generations < 8)] <- "4 - 8"
ROH$Year_class[which(ROH$Age_generations >= 4 & ROH$Age_generations < 8 & ROH$Region == "Autosomal")] <- paste(signif(4*gen_R, 2), "-", signif(8*gen_R, 2), sep=" ")
ROH$Year_class[which(ROH$Age_generations >= 4 & ROH$Age_generations < 8 & ROH$Region == "Z")] <- paste(signif(4*gen_R/0.75, 2), "-", signif(8*gen_R/0.75, 2), sep=" ")
ROH$Age_class[which(ROH$Age_generations >= 8 & ROH$Age_generations < 16)] <- "8 - 16"
ROH$Year_class[which(ROH$Age_generations >= 8 & ROH$Age_generations < 16 & ROH$Region == "Autosomal")] <- paste(signif(8*gen_R, 2), "-", signif(16*gen_R, 2), sep=" ")
ROH$Year_class[which(ROH$Age_generations >= 8 & ROH$Age_generations < 16 & ROH$Region == "Z")] <- paste(signif(8*gen_R/0.75, 2), "-", signif(16*gen_R/0.75, 2), sep=" ")
ROH$Age_class[which(ROH$Age_generations >= 16 & ROH$Age_generations < 32)] <- "16 - 32"
ROH$Year_class[which(ROH$Age_generations >= 16 & ROH$Age_generations < 32 & ROH$Region == "Autosomal")] <- paste(signif(16*gen_R, 2), "-", signif(32*gen_R, 2), sep=" ")
ROH$Year_class[which(ROH$Age_generations >= 16 & ROH$Age_generations < 32 & ROH$Region == "Z")] <- paste(signif(16*gen_R/0.75, 2), "-", signif(32*gen_R/0.75, 2), sep=" ")
ROH$Age_class[which(ROH$Age_generations >= 32 & ROH$Age_generations < 64)] <- "32 - 64"
ROH$Year_class[which(ROH$Age_generations >= 32 & ROH$Age_generations < 64 & ROH$Region == "Autosomal")] <- paste(signif(32*gen_R, 2), "-", signif(64*gen_R, 2), sep=" ")
ROH$Year_class[which(ROH$Age_generations >= 32 & ROH$Age_generations < 64 & ROH$Region == "Z")] <- paste(signif(32*gen_R/0.75, 2), "-", signif(64*gen_R/0.75, 2), sep=" ")
ROH$Age_class[which(ROH$Age_generations >= 64 & ROH$Age_generations < 128)] <- "64 - 128"
ROH$Year_class[which(ROH$Age_generations >= 64 & ROH$Age_generations < 128 & ROH$Region == "Autosomal")] <- paste(signif(64*gen_R, 2), "-", signif(128*gen_R, 2), sep=" ")
ROH$Year_class[which(ROH$Age_generations >= 64 & ROH$Age_generations < 128 & ROH$Region == "Z")] <- paste(signif(64*gen_R/0.75, 2), "-", signif(128*gen_R/0.75, 2), sep=" ")
ROH$Age_class[which(ROH$Age_generations >= 128 & ROH$Age_generations < 256)] <- "128 - 256"
ROH$Year_class[which(ROH$Age_generations >= 128 & ROH$Age_generations < 256 & ROH$Region == "Autosomal")] <- paste(signif(128*gen_R, 2), "-", signif(256*gen_R, 2), sep=" ")
ROH$Year_class[which(ROH$Age_generations >= 128 & ROH$Age_generations < 256 & ROH$Region == "Z")] <- paste(signif(128*gen_R/0.75, 2), "-", signif(256*gen_R/0.75, 2), sep=" ")

ROH <- ROH[-which(ROH$Age_generations >= 256),]
ROH$Age_class <- factor(ROH$Age_class, order=T, levels=c("< 2", "2 - 4", "4 - 8", "8 - 16", "16 - 32", "32 - 64", "64 - 128", "128 - 256"))
ROH$Year_class <- factor(ROH$Year_class, order=T, levels=c(paste("<", signif(2*gen_R, 2), sep=" "), paste("<", signif(2*gen_R/0.75, 2), sep=" "),
                                                           paste(signif(2*gen_R, 2), "-", signif(4*gen_R, 2), sep=" "), paste(signif(2*gen_R/0.75, 2), "-", signif(4*gen_R/0.75, 2), sep=" "),
                                                           paste(signif(4*gen_R, 2), "-", signif(8*gen_R, 2), sep=" "), paste(signif(4*gen_R/0.75, 2), "-", signif(8*gen_R/0.75, 2), sep=" "),
                                                           paste(signif(8*gen_R, 2), "-", signif(16*gen_R, 2), sep=" "), paste(signif(8*gen_R/0.75, 2), "-", signif(16*gen_R/0.75, 2), sep=" "),
                                                           paste(signif(16*gen_R, 2), "-", signif(32*gen_R, 2), sep=" "), paste(signif(16*gen_R/0.75, 2), "-", signif(32*gen_R/0.75, 2), sep=" "),
                                                           paste(signif(32*gen_R, 2), "-", signif(64*gen_R, 2), sep=" "), paste(signif(32*gen_R/0.75, 2), "-", signif(64*gen_R/0.75, 2), sep=" "),
                                                           paste(signif(64*gen_R, 2), "-", signif(128*gen_R, 2), sep=" "), paste(signif(64*gen_R/0.75, 2), "-", signif(128*gen_R/0.75, 2), sep=" "),
                                                           paste(signif(128*gen_R, 2), "-", signif(256*gen_R, 2), sep=" "), paste(signif(128*gen_R/0.75, 2), "-", signif(256*gen_R/0.75, 2), sep=" ")))


ROH_rect <- as.data.frame(matrix(NA, 2*1000/10, 4))
colnames(ROH_rect) <- c("Start", "End", "SumLength", "Region")
for(i in seq(1, nrow(ROH_rect), 2)) {
  ROH_rect$Start[i:(i+1)] <- 2017-i*5-5
  ROH_rect$End[i:(i+1)] <- 2017-i*5+5
  ROH_rect$SumLength[i] <- sum(ROH$Length[which(ROH$Year >= ROH_rect$Start[i] & ROH$Year < ROH_rect$End[i] & ROH$Region == "Autosomal")])
  ROH_rect$SumLength[i+1] <- sum(ROH$Length[which(ROH$Year >= ROH_rect$Start[i+1] & ROH$Year < ROH_rect$End[i+1] & ROH$Region == "Z")])
  ROH_rect$Region[i] <- "Autosomal"
  ROH_rect$Region[i+1] <- "Z"
}
ROH_rect <- rbind(ROH_rect, rep(NA,4))
ROH_rect[nrow(ROH_rect),1:3] <- c(1300, 1301, 1)
ROH_rect[nrow(ROH_rect),4] <- "W"
ROH_rect$Region <- factor(ROH_rect$Region, order=T, levels=c("Autosomal", "Z", "W"))


## ROH Ne
Tot_Size_A <- read.delim("../../Skylark_2021/Scripts/Rasolark_2021/metadata/findZX_regions/region_auto.bed", sep="\t", head=F)
Size_A <- sum(Tot_Size_A$V3 - Tot_Size_A$V2)
Tot_Size_Z <- read.delim("../../Skylark_2021/Scripts/Rasolark_2021/metadata/findZX_regions/region_Z.bed", sep="\t", head=F)
Tot_Size_Z <- rbind(Tot_Size_Z, read.delim("../../Skylark_2021/Scripts/Rasolark_2021/metadata/findZX_regions/region_4A.bed", sep="\t", head=F))
Tot_Size_Z <- rbind(Tot_Size_Z, read.delim("../../Skylark_2021/Scripts/Rasolark_2021/metadata/findZX_regions/region_3a.bed", sep="\t", head=F))
Tot_Size_Z <- rbind(Tot_Size_Z, read.delim("../../Skylark_2021/Scripts/Rasolark_2021/metadata/findZX_regions/region_3b.bed", sep="\t", head=F))
Tot_Size_Z <- rbind(Tot_Size_Z, read.delim("../../Skylark_2021/Scripts/Rasolark_2021/metadata/findZX_regions/region_5.bed", sep="\t", head=F))
Tot_Size_Z <- rbind(Tot_Size_Z, read.delim("../../Skylark_2021/Scripts/Rasolark_2021/metadata/findZX_regions/region_3c.bed", sep="\t", head=F))
Size_Z <- sum(Tot_Size_Z$V3 - Tot_Size_Z$V2)

# Calculate Ne
ROH_Ne <- as.data.frame(matrix(NA, length(levels(ROH$Age_class))*2, 7))
colnames(ROH_Ne) <- c("Age_class", "Region", "FrohMean", "Ne", "max_age", "Age_Start", "Age_End")
ROH_Ne$Age_class <- c(levels(ROH$Age_class), levels(ROH$Age_class))
ROH_Ne$Region <- c(rep("Autosomal", length(levels(ROH$Age_class))), rep("Z", length(levels(ROH$Age_class))))

for(i in 1:nrow(ROH_Ne)) {
  if(length(str_split(ROH_Ne$Age_class[i], " ", simplify=T)) == 2) {
    ROH_Ne$max_age[i] <- as.numeric(str_split(ROH_Ne$Age_class[i], " ", simplify=T)[1,2])
    if(ROH_Ne$Region[i] == "Autosomal") {
      ROH_Ne$Age_Start[i] <- as.numeric(0)
      ROH_Ne$Age_End[i] <- as.numeric(str_split(ROH_Ne$Age_class[i], " ", simplify=T)[1,2]) * gen_R
    } else {
      ROH_Ne$Age_Start[i] <- 0
      ROH_Ne$Age_End[i] <- as.numeric(str_split(ROH_Ne$Age_class[i], " ", simplify=T)[1,2]) * gen_R / 0.75
    }
  } else {
    ROH_Ne$max_age[i] <- as.numeric(str_split(ROH_Ne$Age_class[i], " ", simplify=T)[1,3])
    if(ROH_Ne$Region[i] == "Autosomal") {
      ROH_Ne$Age_Start[i] <- as.numeric(str_split(ROH_Ne$Age_class[i], " ", simplify=T)[1,1]) * gen_R
      ROH_Ne$Age_End[i] <- as.numeric(str_split(ROH_Ne$Age_class[i], " ", simplify=T)[1,3]) * gen_R
    } else {
      ROH_Ne$Age_Start[i] <- as.numeric(str_split(ROH_Ne$Age_class[i], " ", simplify=T)[1,1]) * gen_R / 0.75
      ROH_Ne$Age_End[i] <- as.numeric(str_split(ROH_Ne$Age_class[i], " ", simplify=T)[1,3]) * gen_R / 0.75
    }
  }
  
  if(ROH_Ne$Region[i] == "Autosomal") {
    size <- Size_A
  } else {
    size <- Size_Z
  }
  FrohSum <- 0
  for(k in unique(ROH$Sample[which(ROH$Region == ROH_Ne$Region[i])])) {
    FrohSum <- FrohSum + (sum(ROH$Length[which(ROH$Sample == k & ROH$Region == ROH_Ne$Region[i] & ROH$Age_class == ROH_Ne$Age_class[i])]) / size)
  }
  ROH_Ne$FrohMean[i] <- FrohSum/length(unique(ROH$Sample[which(ROH$Region == ROH_Ne$Region[i])]))
  ROH_Ne$Ne[i] <- as.numeric(-1/((((1-ROH_Ne$FrohMean[i])^(1/ROH_Ne$max_age[i]))-1)*2))
}

ROH_Ne2 <- ROH_Ne |> pivot_longer(cols=c(Age_Start,Age_End), names_to="Age_type", values_to="Age_gen")
ROH_Ne2$Age <- ROH_Ne2$Age_gen * gen_R
ROH_Ne2$Age[which(ROH_Ne2$Region == "Z")] <- ROH_Ne2$Age[which(ROH_Ne2$Region == "Z")] / 0.75
ROH_Ne2$Year <- 2017 - gen_R/2 - ROH_Ne2$Age
ROH_Ne2$Year[which(ROH_Ne2$Region == "Z")] <- 2017 - (gen_R/0.75)/2 - ROH_Ne2$Age[which(ROH_Ne2$Region == "Z")]


### Tajima's D
data_RA <- read.delim("Neutral_popstats/Rasolark/autosomal_windows_10000_steps_10000_exon_dist_20000/Rasolark_2021_autosomal_windows_10000_steps_10000_exon_dist_20000_data.txt", sep="\t", head=T)
data_RZ <- read.delim("Neutral_popstats/Rasolark/homogametic_windows_10000_steps_10000_exon_dist_20000/Rasolark_2021_homogametic_windows_10000_steps_10000_exon_dist_20000_data.txt", sep="\t", head=T)
data_RW <- read.delim("Neutral_popstats/Rasolark/heterogametic_windows_10000_steps_10000_exon_dist_20000/Rasolark_2021_heterogametic_windows_10000_steps_10000_exon_dist_20000_data.txt", sep="\t", head=T)
data_R <- rbind(data_RA, data_RZ, data_RW)
data_R$Species <- "Raso lark"

data_SA <- read.delim("Neutral_popstats/Skylark/autosomal_windows_10000_steps_10000_exon_dist_20000/Skylark_2021_autosomal_windows_10000_steps_10000_exon_dist_20000_data.txt", sep="\t", head=T)
data_SZ <- read.delim("Neutral_popstats/Skylark/homogametic_windows_10000_steps_10000_exon_dist_20000/Skylark_2021_homogametic_windows_10000_steps_10000_exon_dist_20000_data.txt", sep="\t", head=T)
data_SW <- read.delim("Neutral_popstats/Skylark/heterogametic_windows_10000_steps_10000_exon_dist_20000/Skylark_2021_heterogametic_windows_10000_steps_10000_exon_dist_20000_data.txt", sep="\t", head=T)
data_S <- rbind(data_SA, data_SZ, data_SW)
data_S$Species <- "Skylark"

tajima <- rbind(data_R, data_S)
tajima <- tajima[which(tajima$N_callable_sites > 5000 & tajima$N_tot_sites <= 20000),]
tajima$data_type <- factor(tajima$data_type, order=T, labels=c("Autosomal", "Z", "W"), levels=c("autosomal", "homogametic", "heterogametic"))
tajima$Species <- factor(tajima$Species, order=T, labels=rev(c("Skylark", "Raso lark")))


# Census
census_plot <- ggplot() +
  geom_line(data=census, aes(x=Year, y=est.size.individuals), color="black",linetype=1,  linewidth=2.5) +
  geom_line(data=census, aes(x=Year, y=est.size.individuals), color="#E4EAF0",linetype=1,  linewidth=1.5) +
  geom_point(data=census, aes(x=Year, y=est.size.individuals), color="black", shape=18, size=4) +
  geom_point(data=census, aes(x=Year, y=est.size.individuals), color="#E4EAF0", shape=18, size=3) +
  geom_vline(data=sample_raso, aes(xintercept=Year, color="#d7191c"), linewidth=1, linetype=2) +
  annotate(geom="text", x=c(2017-5.5, 2015-4, 2012-4) , y=c(1920, 1800, 1680), label=c("16 inds.", "1 ind.", "1 ind."), color="Black", size=5) +
  # Other
  guides(color="none", linetype="none") +
  ylab(expression("Census count")) +
  scale_x_reverse(expand = c(0, 0), breaks=c(2020,2010,2000,1990,1980,1970,1960), labels=c(2020,2010,2000,1990,1980,1970,1960)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim=c(max_rec, min_rec), ylim=c(0,2000)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.spacing = unit(1, "lines"),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=12, color="black"),
        axis.text.x = element_text(size=12, color="black"))
  

# Gone
gone_plot <- ggplot() +
  geom_line(data=gone, aes(x=Year, y=Geometric_mean, group=Region), color="black", linetype=1, linewidth=2.5) +
  geom_line(data=gone, aes(x=Year, y=Geometric_mean, color=Region, group=Region), linetype=1, linewidth=1.5) +
  scale_color_manual(name="Genomic region", values = c("Autosomal"="#E4EAF0", "Z"="#b30000")) +
  geom_vline(aes(xintercept=1462), color="#d7191c", linewidth=1, linetype=2) +
  annotate(geom="text", x=1540, y=870, label=expression("Human"), color="Black", size=5) +
  annotate(geom="text", x=1540, y=750, label=expression("settlement"), color="Black", size=5) +
  annotate(geom="text", x=1540, y=630, label=expression("1462"), color="Black", size=5) +
  # Other
  guides(color="none") +
  ylab(expression("N"[eLD]*"")) +
  scale_x_reverse(expand = c(0, 0), breaks=c(2000,1900,1800,1700,1600,1500,1400), labels=c(2000,1900,1800,1700,1600,1500,1400)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim=c(max_int, min_int), ylim=c(0,2000)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.spacing = unit(1, "lines"),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=12, color="black"),
        axis.text.x = element_text(size=12, color="black"))

# MSMC

msmc1 <- msmc[which(msmc$Region != "W"),]

msmc_plot <- ggplot() +
  geom_line(data=msmc1, aes(x=Age*1.02, y=Ne*1.04, group=interaction(Species, Region), linetype=Species), color="black", linewidth=1.5) +
  geom_line(data=msmc1, aes(x=Age*0.98, y=Ne*0.96, group=interaction(Species, Region), linetype=Species), color="black", linewidth=1.5) +
  geom_line(data=msmc1, aes(x=Age*1.02, y=Ne*0.96, group=interaction(Species, Region), linetype=Species), color="black", linewidth=1.5) +
  geom_line(data=msmc1, aes(x=Age*0.98, y=Ne*1.04, group=interaction(Species, Region), linetype=Species), color="black", linewidth=1.5) +
  geom_line(data=msmc1, aes(x=Age, y=Ne, group=interaction(Species, Region), color=Region, linetype=Species), linewidth=1.5) +
  scale_color_manual(name="Genomic region", values = c("Autosomal"="#E4EAF0", "Z"="#b30000")) +
  scale_linetype_manual(values = c("Raso lark"=1, "Skylark"=2)) +
  geom_vline(aes(xintercept=20000), color="#d7191c", linewidth=1, linetype=2) +
  annotate(geom="text", x=7000, y=55000000, label=expression("LGM"), color="Black", size=5) +
  annotate(geom="text", x=7000, y=25000000, label=expression("20 kya"), color="Black", size=5) +
  # Other
  guides(linetype=guide_legend(title="Species"), color="none") +
  xlab(expression("Years before present")) +
  ylab(expression("N"[eCOAL]*"")) +
  scale_x_continuous(expand = c(0,0), trans="log10", breaks=c(1000,10000,100000,1000000), labels=c("1 000","10 000","100 000","1 000 000")) +
  scale_y_continuous(expand = c(0,0), trans="log10", breaks=c(100,1000,10000,100000,1000000,10000000,100000000), labels=c(expression("10"^2*""), expression("10"^3*""), expression("10"^4*""),expression("10"^5*""),expression("10"^6*""),expression("10"^7*""),expression("10"^8*""))) +
  coord_cartesian(xlim=c(10^2.5, 1500000), ylim=c(100,100000000)) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(2, 'cm'), #change legend key width
        legend.title = element_text(size=15, face="bold"), #change legend title font size
        legend.text = element_text(size=15), #change legend text font size
        panel.spacing = unit(1, "lines"),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=12, color="black"),
        axis.text.x = element_text(size=12, color="black"))

LD_plot <- ggplot() +
  geom_line(data=data_LD, aes(x=Mid+2000, y=pred+0.002, group=interaction(Species, Region), linetype=Species), color="black", linewidth=1.5) +
  geom_line(data=data_LD, aes(x=Mid-2000, y=pred+0.002, group=interaction(Species, Region), linetype=Species), color="black", linewidth=1.5) +
  geom_line(data=data_LD, aes(x=Mid+2000, y=pred-0.002, group=interaction(Species, Region), linetype=Species), color="black", linewidth=1.5) +
  geom_line(data=data_LD, aes(x=Mid-2000, y=pred-0.002, group=interaction(Species, Region), linetype=Species), color="black", linewidth=1.5) +
  geom_line(data=data_LD, aes(x=Mid, y=pred, group=interaction(Species, Region), color=Region, linetype=Species), linewidth=1.5) +
  scale_color_manual(name="Genomic region", values = c("Autosomal"="#E4EAF0", "Z"="#b30000")) +
  scale_linetype_manual(values = c("Raso lark"=1, "Skylark"=2)) +
  # Other
  guides(color="none", linetype="none") +
  xlab(expression("Distance (Mb)")) +
  ylab(expression("LD (r"^2*")")) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0), labels=scales::label_number(scale = 1e-6)) +
  coord_cartesian(xlim=c(0, 2000000), ylim=c(0,0.5)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.spacing = unit(1, "lines"),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=12, color="black"),
        axis.text.x = element_text(size=12, color="black"))

# Runs of Homozygosity
ROH_bar <- ggplot() + 
  geom_rect(data=ROH_rect, aes(xmin=Start, xmax=End, ymin=0, ymax=SumLength, fill=Region), color="black") +
  scale_fill_manual(name="Genomic region", values = c("Autosomal"="#E4EAF0", "Z"="#b30000", "W"="#fecc5c")) +
  geom_vline(aes(xintercept=1462), color="#d7191c", linewidth=1, linetype=2) +
  annotate(geom="text", x=c(1462), y=c(800000000), label=expression(atop("Human", "settlement", "1462")), color="Black", size=4) +
  # Other
  xlab(expression("Year")) +
  ylab(expression("Sum of RoH (Mb)")) +
  scale_x_reverse(expand = c(0, 0), breaks=c(2000,1900,1800,1700,1600,1500,1400), labels=c(2000,1900,1800,1700,1600,1500,1400)) +
  scale_y_continuous(expand = c(0, 0), labels=scales::label_number(scale = 1e-6)) +
  coord_cartesian(xlim=c(max_int, min_int), ylim=c(0, 425000000)) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(0.7, 'cm'), #change legend key height
        legend.key.width = unit(0.7, 'cm'), #change legend key width
        legend.title = element_text(size=15, face="bold"), #change legend title font size
        legend.text = element_text(size=15), #change legend text font size,
        panel.spacing = unit(1, "lines"),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=12, color="black"),
        axis.text.x = element_text(size=12, color="black"))



# Tajima's D
tajima_plot <- ggplot() +
  geom_violin(data=tajima, aes(x=Species, y=TajD_all, fill=data_type), color="black", width=1, position=position_dodge(1)) +
  geom_boxplot(data=tajima, aes(x=Species, y=TajD_all, fill=data_type), color="black", width=0.05, position=position_dodge(1)) +
  scale_fill_manual(name="Genomic region", values = c("Autosomal"="#E4EAF0", "Z"="#b30000", "W"="#fecc5c")) +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_continuous(limits=c(-4,4), expand=c(0,0)) +
  ylab(expression("Tajima's D")) +
  guides(fill="none", color="none") +
  theme_bw() +
  theme(legend.position = "none",
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=15),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=12, color="black"),
        axis.text.x = element_text(size=15, color="black"))

fig1_plot <- census_plot + gone_plot + msmc_plot + LD_plot + ROH_bar + tajima_plot + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 20, face="bold"), plot.tag.position = c(0.01, 0.98), legend.position = "bottom", legend.box = "horizontal", plot.margin = margin(t = 5.5, r = 10, b = 5.5, l = 5.5))

png("Figures/Figure2.png", width=4500, height=3000, res=300)
fig1_plot
dev.off()


### Supplementary figures
Data_struct <- as.data.frame(matrix(c("Raso lark", "Raso lark", "Skylark", "Skylark", "Raso lark", "Skylark", "Autosomal", "Z", "Autosomal", "Z", "W", "W"), 6, 2))
colnames(Data_struct) <- c("Species", "Region")

# LD decay
plot_list1 <- vector(mode="list", length=4)
for(i in 1:(nrow(Data_struct)-2)) {
  loop_data <- data_LD[which(data_LD$Species == Data_struct$Species[i] & data_LD$Region == Data_struct$Region[i]),]
  
  plot_list1[[i]] <- ggplot()  +
    geom_point(data=loop_data, aes(x=Mid, y=Mean_r2, group=interaction(Species, Region)), color="black", size=0.5) +
    geom_line(data=loop_data, aes(x=Mid+2000, y=pred+0.002, group=interaction(Species, Region), linetype=Species), color="black", linewidth=1.5) +
    geom_line(data=loop_data, aes(x=Mid-2000, y=pred+0.002, group=interaction(Species, Region), linetype=Species), color="black", linewidth=1.5) +
    geom_line(data=loop_data, aes(x=Mid+2000, y=pred-0.002, group=interaction(Species, Region), linetype=Species), color="black", linewidth=1.5) +
    geom_line(data=loop_data, aes(x=Mid-2000, y=pred-0.002, group=interaction(Species, Region), linetype=Species), color="black", linewidth=1.5) +
    geom_line(data=loop_data, aes(x=Mid, y=pred, group=interaction(Species, Region), color=Region, linetype=Species), linewidth=1.5) +
    scale_color_manual(name="Genomic region", values = c("Autosomal"="#E4EAF0", "Z"="#b30000")) +
    scale_linetype_manual(values = c("Raso lark"=1, "Skylark"=2)) +
    facet_wrap(Species~Region) +
    # Other
    guides(color="none", linetype="none") +
    xlab(expression("Distance (Mb)")) +
    ylab(expression("LD (r"^2*")")) +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0), labels=scales::label_number(scale = 1e-6)) +
    coord_cartesian(xlim=c(0, 2000000), ylim=c(0,0.5)) +
    theme_bw() +
    theme(panel.spacing = unit(1, "lines"),
          strip.background = element_rect(color="black", fill="white", linewidth=1),
          strip.text = element_text(size=20, color="black"),
          axis.line = element_line(colour = "black"),
          axis.title.y = element_text(size=20),
          axis.title.x = element_text(size=20),
          axis.text.y = element_text(size=15, color="black"),
          axis.text.x = element_text(size=15, color="black"))
  
  
}

png("Figures/Supplementary/LD_decay.png", width=4000, height=4000, res=300)
(plot_list1[[1]] + plot_list1[[3]]) / (plot_list1[[2]] + plot_list1[[4]]) +
  plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 20, face="bold"), plot.tag.position=c(0.05,0.95))
dev.off()


# Gone
plot_list2 <- vector(mode="list", length=2)
gone$Species <- "Raso lark"
for(i in 1:(nrow(Data_struct)-2)) {
  loop_data <- gone[which(gone$Region == Data_struct$Region[i]),]
  plot_list2[[i]] <- ggplot() +
    geom_ribbon(data=loop_data, aes(x=Year, ymin=Lower95, ymax=Upper95, group=Region, fill=Region), alpha=0.5) +
    geom_line(data=loop_data, aes(x=Year, y=Geometric_mean, group=Region), color="black", linetype=1, linewidth=2.5) +
    geom_line(data=loop_data, aes(x=Year, y=Geometric_mean, color=Region, group=Region), linetype=1, linewidth=1.5) +
    geom_line(data=loop_data, aes(x=Year, y=Upper95, color=Region, group=Region), linetype=1, linewidth=1.5) +
    geom_line(data=loop_data, aes(x=Year, y=Lower95, color=Region, group=Region), linetype=1, linewidth=1.5) +
    scale_color_manual(name="Genomic region", values = c("Autosomal"="#E4EAF0", "Z"="#b30000")) +
    scale_fill_manual(name="Genomic region", values = c("Autosomal"="#E4EAF0", "Z"="#b30000")) +
    geom_vline(aes(xintercept=1462), color="#d7191c", linewidth=1, linetype=2) +
    annotate(geom="text", x=1540, y=max(loop_data$Upper95)/1.5, label=expression("Human"), color="Black", size=5) +
    annotate(geom="text", x=1540, y=max(loop_data$Upper95)/1.5 - max(loop_data$Upper95)/20, label=expression("settlement"), color="Black", size=5) +
    annotate(geom="text", x=1540, y=max(loop_data$Upper95)/1.5 - max(loop_data$Upper95)/10, label=expression("1462"), color="Black", size=5) +
    facet_wrap(Species~Region) +
    # Other
    guides(color="none", fill="none") +
    ylab(expression("N"[eLD]*"")) +
    scale_x_reverse(expand = c(0, 0), breaks=c(2000,1800,1600,1400,1200,1000), labels=c(2000,1800,1600,1400, 1200, 1000)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim=c(max_int, 1000)) +
    theme_bw() +
    theme(panel.spacing = unit(1, "lines"),
          strip.background = element_rect(color="black", fill="white", linewidth=1),
          strip.text = element_text(size=20, color="black"),
          axis.line = element_line(colour = "black"),
          axis.title.y = element_text(size=20),
          axis.title.x = element_text(size=20),
          axis.text.y = element_text(size=15, color="black"),
          axis.text.x = element_text(size=15, color="black"))
   
}
png("Figures/Supplementary/Gone.png", width=3000, height=4000, res=300)
(plot_list2[[1]] / plot_list2[[2]]) +
  plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 20, face="bold"), plot.margin = margin(c(5,10,5,5)), plot.tag.position=c(0.05,0.95))
dev.off()


# RoH Ne
RoH_Ne <- ggplot() +
  geom_line(data=ROH_Ne2, aes(x=Year, y=Ne, group=Region), color="black", linetype=1, linewidth=2.5) +
  geom_line(data=ROH_Ne2, aes(x=Year, y=Ne, color=Region, group=Region), linetype=1, linewidth=1.5) +
  scale_color_manual(name="Genomic region", values = c("Autosomal"="#E4EAF0", "Z"="#b30000")) +
  geom_vline(aes(xintercept=1462), color="#d7191c", linewidth=1, linetype=2) +
  annotate(geom="text", x=1540, y=470, label=expression("Human"), color="Black", size=5) +
  annotate(geom="text", x=1540, y=450, label=expression("settlement"), color="Black", size=5) +
  annotate(geom="text", x=1540, y=430, label=expression("1462"), color="Black", size=5) +
  # Other
  guides(fill="none") +
  ylab(expression("N"[eRoH]*"")) +
  scale_x_reverse(expand = c(0, 0), breaks=c(2000,1800,1600,1400), labels=c(2000,1800,1600,1400)) +
  scale_y_continuous(expand = c(0, 0),) +
  coord_cartesian(xlim=c(max_int, 1000), ylim=c(0,500)) +
  theme_bw() +
  theme(legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(0.7, 'cm'), #change legend key height
        legend.key.width = unit(0.7, 'cm'), #change legend key width
        legend.title = element_text(size=15, face="bold"), #change legend title font size
        legend.text = element_text(size=15), #change legend text font size,
        panel.spacing = unit(1, "lines"),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=12, color="black"),
        axis.text.x = element_text(size=12, color="black"))


png("Figures/Supplementary/RohNe.png", width=3000, height=2000, res=300)
RoH_Ne
dev.off()



# MSMC2
plot_list3 <- vector(mode="list", length=4)
for(i in 1:(nrow(Data_struct)-2)) {
  loop_data <- msmc[which(msmc$Species == Data_struct$Species[i] & msmc$Region == Data_struct$Region[i]),]
  plot_list3[[i]] <- ggplot() +
    geom_ribbon(data=loop_data, aes(x=Age, ymin=Lower95, ymax=Upper95, group=Region, fill=Region), alpha=0.5) +
    geom_line(data=loop_data, aes(x=Age, y=Upper95, group=interaction(Species, Region), color=Region, linetype=Species), linewidth=1.5) +
    geom_line(data=loop_data, aes(x=Age, y=Lower95, group=interaction(Species, Region), color=Region, linetype=Species), linewidth=1.5) +
    geom_line(data=loop_data, aes(x=Age*1.02, y=Ne*1.04, group=interaction(Species, Region), linetype=Species), color="black", linewidth=1.5) +
    geom_line(data=loop_data, aes(x=Age*0.98, y=Ne*0.96, group=interaction(Species, Region), linetype=Species), color="black", linewidth=1.5) +
    geom_line(data=loop_data, aes(x=Age*1.02, y=Ne*0.96, group=interaction(Species, Region), linetype=Species), color="black", linewidth=1.5) +
    geom_line(data=loop_data, aes(x=Age*0.98, y=Ne*1.04, group=interaction(Species, Region), linetype=Species), color="black", linewidth=1.5) +
    geom_line(data=loop_data, aes(x=Age, y=Ne, group=interaction(Species, Region), color=Region, linetype=Species), linewidth=1.5) +
    scale_color_manual(name="Genomic region", values = c("Autosomal"="#E4EAF0", "Z"="#b30000")) +
    scale_fill_manual(name="Genomic region", values = c("Autosomal"="#E4EAF0", "Z"="#b30000")) +
    scale_linetype_manual(values = c("Raso lark"=1, "Skylark"=2)) +
    geom_vline(aes(xintercept=20000), color="#d7191c", linewidth=1, linetype=2) +
    annotate(geom="text", x=7000, y=100000000/2, label=expression("LGM"), color="Black", size=5) +
    annotate(geom="text", x=7000, y=100000000/2  - 100000000/4, label=expression("20 kya"), color="Black", size=5) +
    facet_wrap(Species~Region) +
    # Other
    guides(linetype="none", fill="none", color="none") +
    xlab(expression("Years before present")) +
    ylab(expression("N"[eCOAL]*"")) +
    scale_x_continuous(expand = c(0,0), trans="log10", breaks=c(100,1000,10000,100000,1000000), labels=c("100","1 000","10 000","100 000","1 000 000")) +
    scale_y_continuous(expand = c(0,0), trans="log10", breaks=c(100,1000,10000,100000,1000000,10000000,100000000), labels=c(expression("10"^2*""), expression("10"^3*""), expression("10"^4*""),expression("10"^5*""),expression("10"^6*""),expression("10"^7*""),expression("10"^8*""))) +
    coord_cartesian(xlim=c(10, 10000000), ylim=c(100, 100000000)) +
    theme_bw() +
    theme_bw() +
    theme(panel.spacing = unit(1, "lines"),
          strip.background = element_rect(color="black", fill="white", linewidth=1),
          strip.text = element_text(size=20, color="black"),
          axis.line = element_line(colour = "black"),
          axis.title.y = element_text(size=20),
          axis.title.x = element_text(size=20),
          axis.text.y = element_text(size=15, color="black"),
          axis.text.x = element_text(size=15, color="black"))
}

png("Figures/Supplementary/msmc.png", width=4000, height=5000, res=300)
(plot_list3[[1]] + plot_list3[[3]]) / (plot_list3[[2]] + plot_list3[[4]]) +
  plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 20, face="bold"), plot.tag.position=c(0.05,0.95))
dev.off()


