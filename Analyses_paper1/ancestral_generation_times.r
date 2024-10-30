#!/usr/bin/Rscript

## Export variables and load libraries
rm(list=ls())
options(scipen=999)

library(tidyverse)
library(ape, lib.loc="/home/simonj/R_lib")
library(phytools, lib.loc="/home/simonj/R_lib")
library(geiger)
library(nlme)
library(evomap, lib.loc="/home/simonj/R_lib")
library(RRphylo, lib.loc="/home/simonj/R_lib")

args <- strsplit(commandArgs(trailingOnly=T), "=", fixed=T)
for(i in 1:length(args)) {
  assign(args[[i]][1], args[[i]][2])
}

ITERATIONS <- as.numeric(ITERATIONS)
BURNIN <- as.numeric(BURNIN)

### Import and prepare data
## Import generation times
data <- read.table(paste(RESOURCES, "/", "GenLength_Bird.tsv", sep=""), sep="\t", head=T, dec=",")
data <- as.data.frame(cbind(data$Scientific.name, data$Family, data$GenLength))
colnames(data) <- c("Species", "Family", "GenLength")
data$GenLength <- as.numeric(data$GenLength)

## Transform data
data$InvGenLength <- 1/(data$GenLength)

#Plot distribuions
jpeg(paste(OUTDIR,"/Generation_length_distribuion.jpeg", sep=""), width=1000, height=1000, quality=100)
hist(data$GenLength)
dev.off()
jpeg(paste(OUTDIR,"/Inverse_Generation_length_distribuion.jpeg", sep=""), width=1000, height=1000, quality=100)
hist(data$InvGenLength)
dev.off()

## Import tree
tree <- read.nexus(paste(RESOURCES, "/", "63k_dated_Stiller.tre", sep=""))
tree$tip.label[which(tree$tip.label == "Alaudala_cheleensis")] <- "Alauda_arvensis"
data2 <- data[0,]
data2 <- rbind(data2, matrix(NA, length(tree$tip.label), ncol(data)))

## Correct taxonomic and misspellings
data$Species <- gsub(" ", "_", data$Species)
data$Species[which(data$Species == "Buphagus_erythrorynchus")] <- "Buphagus_erythrorhynchus"
data$Species[which(data$Species == "Tychaedon_coryphoeus")] <- "Cercotrichas_coryphaeus"
data$Species[which(data$Species == "Larus_maculipennis")] <- "Chroicocephalus_maculipennis"
data$Species[which(data$Species == "Corvus_corone")] <- "Corvus_cornix"
data$Species[which(data$Species == "Eolophus_roseicapilla")] <- "Eolophus_roseicapillus"
data$Species[which(data$Species == "Horornis_flavolivaceus")] <- "Horornis_vulcanius"
data$Species[which(data$Species == "Urile_pelagicus")] <- "Phalacrocorax_pelagicus"
data$Species[which(data$Species == "Nannopterum_auritus")] <- "Phalacrocorax_auritus"
data$Species[which(data$Species == "Nannopterum_brasilianus")] <- "Phalacrocorax_brasilianus"
data$Species[which(data$Species == "Nannopterum_harrisi")] <- "Phalacrocorax_harrisi"
data$Species[which(data$Species == "Dryobates_pubescens")] <- "Picoides_pubescens"
data$Species[which(data$Species == "Phylloscopus_sibilatrix")] <- "Rhadina_sibilatrix"
data$Species[which(data$Species == "Smutsornis_africanus")] <- "Rhinoptilus_africanus"
data$Species[which(data$Species == "Saxicola_torquatus")] <- "Saxicola_maurus"

## Match species with generation time
colnames(data2) <- colnames(data)
for(i in 1:length(tree$tip.label)) {
  data2[i,] <- data[which(data$Species == tree$tip.label[i]),]
}

## Check names between files
rownames(data2) <- data2$Species
tree <- treedata(tree, data2, sort=T, warnings=T)$phy
data2 <- as.data.frame(treedata(tree, data2, sort=T, warnings=T)$data)
data2$GenLength <- as.numeric(data2$GenLength)
data2$InvGenLength <- as.numeric(data2$InvGenLength)

## Get mvBM rescaled tree
x <- as.numeric(data2$InvGenLength)
names(x) <- rownames(data2)
BMsigma2 <- ace(x=x, phy=tree, method="REML")$sigma2[1] # get BM sigma2 ('ace' requires the 'ape' package)
mvBMresults <- mvBM(x, tree, BMsigma2) # calculate rescaled branch lengths using mvBM
tree_mvBM <- tree
tree_mvBM$edge.length <- mvBMresults$rBL # create new tree with rescaled branch lengths
ML_mvBM_anc <- ace(x, tree_mvBM, method="REML") # get ancestral estimates using mvBM tree

## Map ML generation times on trees and plot
jpeg(paste(OUTDIR,"/Inverse_Generation_length_tree_ML.jpeg", sep=""), width=1000, height=1000, quality=100)
ML_tree_InvGenLength <- contMap(tree, x, anc.states=ML_mvBM_anc, plot=T, fsize=0.5, res=1000)
dev.off()
jpeg(paste(OUTDIR,"/Generation_length_tree_ML.jpeg", sep=""), width=1000, height=1000, quality=100)
ML_tree_GenLength <- contMap(tree, 1/x, anc.states=(1/ML_mvBM_anc$ace), plot=T, fsize=0.5, res=1000)
dev.off()

## Use mvBM rescaled tree in MCMC analysis
model_mvBM <- anc.Bayes(tree_mvBM, x, ngen=ITERATIONS) # 'anc.Bayes' requires the 'phytools' package; this may take a few minutes
MCMC_mvBM_sigma2 <- model_mvBM$mcmc[(BURNIN/100):nrow(model_mvBM$mcmc),2] # remove burnin
MCMC_mvBM_anc <- colMeans(model_mvBM$mcmc[(BURNIN/100):nrow(model_mvBM$mcmc),])[3:(3 + tree$Nnode - 1)] # These are the ancestral estimates
MCMC_mvBM_logLik_all <- model_mvBM$mcmc[,tree$Nnode + 3]
MCMC_mvBM_logLik_burnin <- model_mvBM$mcmc[(BURNIN/100):nrow(model_mvBM$mcmc),tree$Nnode + 3]

# Plot posterior log likelihood distribution
jpeg(paste(OUTDIR,"/MCMC_log_likelihood_disribution_all.jpeg", sep=""), width=1000, height=1000, quality=100)
hist(MCMC_mvBM_logLik_all)
dev.off()
jpeg(paste(OUTDIR,"/MCMC_log_likelihood_disribution_Burn_in.jpeg", sep=""), width=1000, height=1000, quality=100)
hist(MCMC_mvBM_logLik_burnin)
dev.off()

## Map MCMC generation times on trees and plot
jpeg(paste(OUTDIR,"/Inverse_Generation_length_tree_MCMC.jpeg", sep=""), width=1000, height=1000, quality=100)
MCMC_tree_InvGenLength <- contMap(tree, x, anc.states=MCMC_mvBM_anc, plot=T, fsize=0.5, res=1000)
dev.off()
jpeg(paste(OUTDIR,"/Generation_length_tree_MCMC.jpeg", sep=""), width=1000, height=1000, quality=100)
MCMC_tree_GenLength <- contMap(tree, 1/x, anc.states=(1/MCMC_mvBM_anc), plot=T, fsize=0.5, res=1000)
dev.off()

# Scale years to generations
root2node <- c(getMommy(tree, "Alauda_arvensis"))
GenLen <- as.data.frame(matrix(NA, (length(root2node)+1), 10))
colnames(GenLen) <- c("age_start", "age_end", "node_start", "node_end", "edge_length", "Anc_gen_len_ML", "Anc_gen_len_MCMC", "cumAge", "cumGenML", "cumGenMCMC")
GenLen[1,] <- c(rep(0, 2), which(tree$tip.label == "Alauda_arvensis"), root2node[1], rep(0, 6))
GenLen$edge_length[1] <- tree$edge.length[which(tree$edge[,1] == root2node[1] & tree$edge[,2] == which(tree$tip.label == "Alauda_arvensis"))]
GenLen$age_end[1] <- GenLen$edge_length[1]
GenLen$Anc_gen_len_ML[1] <- 1/ML_mvBM_anc$ace[which(names(ML_mvBM_anc$ace) == GenLen$node_end[1])]
GenLen$Anc_gen_len_MCMC[1] <- 1/MCMC_mvBM_anc[which(names(MCMC_mvBM_anc) == GenLen$node_end[1])]
GenLen$cumAge[1] <- GenLen$edge_length[1]
GenLen$cumGenML[1] <- GenLen$edge_length[1] / GenLen$Anc_gen_len_ML[1]
GenLen$cumGenMCMC[1] <- GenLen$edge_length[1] / GenLen$Anc_gen_len_MCMC[1]

for(i in 2:length(root2node)) {
  GenLen$age_start[i] <- GenLen$age_end[i-1]
  GenLen$node_start[i] <- root2node[i-1]
  GenLen$node_end[i] <- root2node[i]
  GenLen$edge_length[i] <- tree$edge.length[which(tree$edge[,1] == GenLen$node_end[i] & tree$edge[,2] == GenLen$node_start[i])]
  GenLen$age_end[i] <- GenLen$age_start[i] + GenLen$edge_length[i]
  GenLen$Anc_gen_len_ML[i] <- 1/ML_mvBM_anc$ace[which(names(ML_mvBM_anc$ace) == GenLen$node_end[i])]
  GenLen$Anc_gen_len_MCMC[i] <- 1/MCMC_mvBM_anc[which(names(MCMC_mvBM_anc) == GenLen$node_end[i])]
  GenLen$cumAge[i] <-   GenLen$age_end[i]
  GenLen$cumGenML[i] <- GenLen$cumGenML[i-1] + (GenLen$edge_length[i] / GenLen$Anc_gen_len_ML[i])
  GenLen$cumGenMCMC[i] <- GenLen$cumGenMCMC[i-1] + (GenLen$edge_length[i] / GenLen$Anc_gen_len_MCMC[i])
}
i <- i + 1
GenLen$age_start[i] <- GenLen$age_end[i-1]
GenLen$node_start[i] <- root2node[i-1]
GenLen$age_end[i] <- 135.250
GenLen$edge_length[i] <- GenLen$age_end[i] - GenLen$age_start[i]
GenLen$Anc_gen_len_ML[i] <- GenLen$Anc_gen_len_ML[i-1]
GenLen$Anc_gen_len_MCMC[i] <- GenLen$Anc_gen_len_MCMC[i-1]
GenLen$cumAge[i] <-   GenLen$age_end[i]
GenLen$cumGenML[i] <- GenLen$cumGenML[i-1] + (GenLen$edge_length[i] / GenLen$Anc_gen_len_ML[i])
GenLen$cumGenMCMC[i] <- GenLen$cumGenMCMC[i-1] + (GenLen$edge_length[i] / GenLen$Anc_gen_len_MCMC[i])

## Import strata ages
Strata_age <- read.table(paste(METADATA, "/Strata_ages.tsv", sep=""), sep="\t", head=T)
Strata_age$Age_years <- as.numeric(Strata_age$Age_years)
Strata_age$Age_years_Stiller <- as.numeric(Strata_age$Age_years_Stiller)
Strata_age$Age_generations <- as.numeric(rep(NA, nrow(Strata_age)))
Strata_age$Age_generations_Stiller <- as.numeric(rep(NA, nrow(Strata_age)))

## Fit a second degree polynomial model to generation time as a function of age, with fixed intercept at 0
fit <- lm(cumGenMCMC ~ -1 + cumAge + I(cumAge^2) + offset(rep(0, length(cumAge))), data=GenLen)
summary(fit)
sink(paste(OUTDIR, "/age_gen_fit.txt", sep=""))
print(summary(fit))
sink()
saveRDS(fit, paste(OUTDIR, "/age_gen_model.RDS", sep=""))

Age_years <- as.data.frame(Strata_age$Age_years)
colnames(Age_years) <- "cumAge"
Strata_age$Age_generations <- predict(fit, newdata=Age_years)

pred_seq <- as.data.frame(seq(0, max(Age_years$cumAge)))
colnames(pred_seq) <- "cumAge"
pred_seq$predGen <- predict(fit, newdata=pred_seq)

## Rescale strata ages from years to generations
for(i in 1:nrow(Strata_age)) {
  index <- which(GenLen$age_start <= Strata_age$Age_years[i] & GenLen$age_end >= Strata_age$Age_years[i])
  if(index == 1) {
    Strata_age$Age_generations_Stiller[i] <- Strata_age$Age_years_Stiller[i] * GenLen$Anc_gen_len_MCMC[index]
  } else {
      Strata_age$Age_generations_Stiller[i] <- GenLen$cumGenMCMC[index-1] + ((Strata_age$Age_years_Stiller[i] - GenLen$age_start[index]) * GenLen$Anc_gen_len_MCMC[index])
  }
}

# Write years to generations scale as table
write.table(GenLen, paste(OUTDIR, "/Years_to_generations_scale.tsv", sep=""), quote=F, sep='\t', row.names = F, col.names = T)

# Write rescaled strata ages as table
write.table(Strata_age, paste(OUTDIR, "/Strata_ages_generations.tsv", sep=""), quote=F, sep='\t', row.names = F, col.names = T)

## Save MCMC Analysis and R session in file
save.image(paste(OUTDIR, "/MCMC_analysis.RData", sep=""))

## Plot final graphs for supplementary material
jpeg(paste(OUTDIR,"/Generation_length_tree_ML_final.jpeg", sep=""), width=5000, height=5000, res=300)
ML_tree_GenLength <- contMap(tree, 1/x, anc.states=(1/ML_mvBM_anc), plot=T, fsize=0.5, res=1000)
dev.off()

jpeg(paste(OUTDIR,"/Generation_length_tree_MCMC_final.jpeg", sep=""), width=5000, height=5000, res=300)
MCMC_tree_GenLength <- contMap(tree, 1/x, anc.states=(1/MCMC_mvBM_anc), plot=T, fsize=0.5, res=1000)
dev.off()

plot_func <- ggplot() +
  geom_line(data=pred_seq, aes(x=cumAge, y=predGen, color="Fitted model")) +
  geom_point(data=GenLen, aes(x=cumAge, y=cumGenMCMC, color="Model data"), size=3) +
  geom_point(data=Strata_age, aes(x=Age_years, y=Age_generations, color="Predicted stratum ages"), shape=17, size=4) +
  scale_color_manual(values=c("Black", "Blue", "Red"), limits=c("Fitted model", "Model data", "Predicted stratum ages")) +
  labs(x ="Age (Millions of years)", y = "Age (Millions of genertions)") +
  labs(colour="Data") +
  theme_bw() +
  theme(legend.key.size = unit(1, 'line'), #change legend key size
       legend.key.height = unit(1, 'cm'), #change legend key height
       legend.key.width = unit(1, 'cm'), #change legend key width
       legend.title = element_text(size=20), #change legend title font size
       legend.text = element_text(size=20), #change legend text font size
       strip.text = element_text(size=20, color="black"),
       axis.line = element_line(colour = "black"),
       axis.title.y = element_text(size=20),
       axis.title.x = element_text(size=20),
       axis.text.y = element_text(size=20, color="black"),
       axis.text.x = element_text(size=20, color="black"))

jpeg(paste(OUTDIR,"/fitted_function.jpeg", sep=""), width=6000, height=3000, res=300)
plot_func
dev.off()
