## Export variables and load libraries
rm(list=ls())
library(tidyverse)
library(patchwork)
library(performance)
library(glmmTMB)
library(DHARMa)
library(multcomp)
library(MuMIn)
library(ggeffects)
library(car)
library(scales)
library(effectsize)
library(ggridges)
sessionInfo()

options(scipen=999)
setwd("C:/Users/Simon JE/OneDrive - Lund University/Dokument/Simon/PhD/Projects/Skylark_2021/Results/")
data <- read.delim("Genes/Skylark_2021_Rasolark_2021_organised_data1.tsv", sep="\t", head=T)


data_deg <- data[which(data$Filter3=="OK" & data$Filter4=="OK" & data$Filter5=="OK"),]
data_deg$Strata2 <- data_deg$Strata
data_deg$Strata[which(data_deg$Strata == "PAR3" | data_deg$Strata == "PAR5")] <- "PAR"
data_deg$Strata <- factor(data_deg$Strata, order=T, labels=rev(c("S0", "S1", "S2", "S3", "4A", "3-a", "3-b", "5", "3-c", "PAR", "Autosomal")), levels=rev(c("S0", "S1", "S2", "S3", "4A", "3a", "3b", "5", "3c", "PAR", "autosomal")))
data_deg$Species <- factor(data_deg$Species, order=T, labels=c("Skylark", "Raso lark"), levels=c("Skylark", "Rasolark"))


# Statistics
# Set up data

### Selection on W as a function of Z selection and other predictors

# Remove and missing data
data_deg5 <- data_deg[which(data_deg$Filter1=="OK"),]
data_deg5 <- data_deg5[which(!is.na(data_deg5$pHaplo)),]

data_deg5 <- data_deg5[which(data_deg5$Region != "autosomal" & data_deg5$Wdegeneration != "W degenerated" & data_deg5$Wdegeneration != "W partially degenerated"),]
data_deg5$SexLinkedDNDS_DOSW <- rep(NA, nrow(data_deg5))
data_deg5$SexLinkedDNDS_DOSZ <- rep(NA, nrow(data_deg5))
data_deg5$SexLinkedDNDS_DOSW <- data_deg5$SexLinkedBetaW/(data_deg5$SexLinkedBetaW + data_deg5$SexLinkedAlphaW)
data_deg5$SexLinkedDNDS_DOSZ <- data_deg5$SexLinkedBetaZ/(data_deg5$SexLinkedBetaZ + data_deg5$SexLinkedAlphaZ)
data_deg5 <- data_deg5[which(data_deg5$SexLinkedAlphaW < 1 & data_deg5$SexLinkedAlphaW > 0.001 & data_deg5$SexLinkedBetaW < 1 & data_deg5$SexLinkedAlphaZ < 1 & data_deg5$SexLinkedAlphaZ > 0.001 & data_deg5$SexLinkedBetaZ < 1),]

# Transform data
hist(data_deg5$SexLinkedDNDS_DOSZ, breaks=1000)
hist(data_deg5$SexLinkedDNDS_DOSW, breaks=1000)
data_deg5 <- data_deg5[which(data_deg5$SexLinkedDNDS_DOSW != 0 & data_deg5$SexLinkedDNDS_DOSZ != 0),]
data_deg5 <- data_deg5[which(!is.na(data_deg5$SexLinkedDNDS_DOSW) & !is.na(data_deg5$SexLinkedDNDS_DOSZ)),]
hist(data_deg5$SexLinkedDNDS_DOSZ, breaks=1000)
hist(data_deg5$SexLinkedDNDS_DOSW, breaks=1000)
genes <- names(which(table(data_deg5$geneID) == 2))
data_deg5 <- data_deg5[data_deg5$geneID %in% genes, ]
data_deg5$logSexLinkedDNDS_DOSW <- log10(data_deg5$SexLinkedDNDS_DOSW)

hist(data_deg5$logSexLinkedDNDS_DOSW)
qqnorm(data_deg5$logSexLinkedDNDS_DOSW)
qqline(data_deg5$logSexLinkedDNDS_DOSW)
hist(data_deg5$SexLinkedDNDS_DOSZ)
qqnorm(data_deg5$SexLinkedDNDS_DOSZ)
qqline(data_deg5$SexLinkedDNDS_DOSZ)

hist(data_deg5$geneLengthDataBase)
qqnorm(data_deg5$geneLengthDataBase)
qqline(data_deg5$geneLengthDataBase)
data_deg5$logGeneLen <- log10(data_deg5$geneLengthDataBase)
hist(data_deg5$logGeneLen)
qqnorm(data_deg5$logGeneLen)
qqline(data_deg5$logGeneLen)

# Is the functional variable evenly distributed?
table(data_deg5$Wdegeneration)
data_deg5$Wdegeneration <- factor(data_deg5$Wdegeneration, order=F, labels=c("W_functional", "W_loss_of_function"), levels=c("W functional", "W loss of function"))
table(data_deg5$Wdegeneration)
table(data_deg5$Strata)

### Test for colinearity of all factors
colin_model1 <- glmmTMB(logSexLinkedDNDS_DOSW ~
                          SexLinkedDNDS_DOSZ + pHaplo + logGeneLen + Species + Wdegeneration + Strata_Age_Generations,
                        family=gaussian, data = data_deg5, na.action = "na.fail", REML=F)
check_collinearity(colin_model1)

### Test model
globalmodel1 <- glmmTMB(logSexLinkedDNDS_DOSW ~
                          SexLinkedDNDS_DOSZ + pHaplo + logGeneLen + Wdegeneration + Species + Strata_Age_Generations +
                          SexLinkedDNDS_DOSZ:pHaplo + SexLinkedDNDS_DOSZ:logGeneLen + SexLinkedDNDS_DOSZ:Wdegeneration +
                          SexLinkedDNDS_DOSZ:Species + SexLinkedDNDS_DOSZ:Strata_Age_Generations +
                          pHaplo:logGeneLen + pHaplo:Wdegeneration + pHaplo:Species + pHaplo:Strata_Age_Generations +
                          logGeneLen:Wdegeneration + logGeneLen:Species + logGeneLen:Strata_Age_Generations +
                          Wdegeneration:Species + Wdegeneration:Strata_Age_Generations + Species:Strata_Age_Generations,
                        family=gaussian, data = data_deg5, na.action = "na.fail", REML=F)

options(na.action = "na.omit")

simulateResiduals(fittedModel = globalmodel1, plot = T) # Quantile deviations detected. Combined adjusted quantile test n.s.

combinations1 <- dredge(global.model=globalmodel1, rank="AIC")

# Use globalmodel1
print(combinations1)
# Top 11 models have delta < 2 (max delta 1.97, cum weight 0.039)
0.007+0.004*2+0.003*8


# Check top 4 models with Dharma
simulateResiduals(fittedModel = get.models(combinations1, subset = 1)[[1]], plot = T) # OK
simulateResiduals(fittedModel = get.models(combinations1, subset = 2)[[1]], plot = T) # OK
simulateResiduals(fittedModel = get.models(combinations1, subset = 3)[[1]], plot = T) # OK
simulateResiduals(fittedModel = get.models(combinations1, subset = 4)[[1]], plot = T) # OK
simulateResiduals(fittedModel = get.models(combinations1, subset = 5)[[1]], plot = T) # OK
simulateResiduals(fittedModel = get.models(combinations1, subset = 6)[[1]], plot = T) # OK
simulateResiduals(fittedModel = get.models(combinations1, subset = 7)[[1]], plot = T) # OK
simulateResiduals(fittedModel = get.models(combinations1, subset = 8)[[1]], plot = T) # OK
simulateResiduals(fittedModel = get.models(combinations1, subset = 9)[[1]], plot = T) # OK
simulateResiduals(fittedModel = get.models(combinations1, subset = 10)[[1]], plot = T) # OK
simulateResiduals(fittedModel = get.models(combinations1, subset = 11)[[1]], plot = T) # OK


# Get model average
FinalModel <- model.avg(get.models(combinations1, subset = delta <= 2))
summary(FinalModel)
#LogGeneLen*
#pHaplo*
#SexLinkedDNDS_DOSZ*
#logGeneLen:SexLinkedDNDS_DOSZ*
#pHaplo:SexLinkedDNDS_DOSZ**
#cond(Strata_Age_Generations:WdegenerationW_loss_of_function) **


terms <- sub("cond\\((.*)\\)", "\\1", colnames(FinalModel$coefficients))
terms[1] <- "(Intercept)"
ef_mat <- as.data.frame(matrix(0, ncol(FinalModel$coefficients), 11))
rownames(ef_mat) <- terms
eflow_mat <-as.data.frame(matrix(0, ncol(FinalModel$coefficients), 11))
rownames(eflow_mat) <- terms
efupp_mat <-as.data.frame(matrix(0, ncol(FinalModel$coefficients), 11))
rownames(efupp_mat) <- terms
weights <- FinalModel$msTable$weight


R2 <- rep(NA, 11)
R2_adj <- rep(NA, 11)
for(i in 1:11) {
  R2_0 <-  r2(get.models(combinations1, subset = i)[[1]])
  R2[i] <- R2_0$R2
  R2_adj[i] <- R2_0$R2_adjusted
  ef <- standardize_parameters(get.models(combinations1, subset = i)[[1]])
  for(j in 1:length(terms)) {
    if(length(which(ef$Parameter == terms[j])) > 0) {
      ef_mat[j,i] <- ef$Std_Coefficient[which(ef$Parameter == terms[j])]*weights[i]
      eflow_mat[j,i] <- ef$CI_low[which(ef$Parameter == terms[j])]*weights[i]
      efupp_mat[j,i] <- ef$CI_high[which(ef$Parameter == terms[j])]*weights[i]
    }
  }
}

min(R2) #0.2459441 = 0.246
max(R2) #[1] 0.2583375  = 0.258

min(R2_adj) #0.231443 = 0.231
max(R2_adj) #0.2375582 = 0.238



effects_df <- as.data.frame(matrix(NA, length(terms), 4))
colnames(effects_df) <- c("Parameter", "Std_Coefficien", "CI_low", "CI_high")
for(i in 1:length(terms)) {
  effects_df$Parameter[i] <- terms[i]
  effects_df$Std_Coefficien[i] <- sum(ef_mat[i,1:11])
  effects_df$CI_low[i] <- sum(eflow_mat[i,1:11])
  effects_df$CI_high[i] <- sum(efupp_mat[i,1:11])
}


#Plot as function of Z selection
pred_seq1a <- as.data.frame(seq(min(data_deg5$SexLinkedDNDS_DOSZ)*100, max(data_deg5$SexLinkedDNDS_DOSZ)*100)/100)
colnames(pred_seq1a) <- "SexLinkedDNDS_DOSZ"
pred_seq1a$pHaplo <- median(data_deg5$pHaplo)
pred_seq1a$Strata_Age_Generations <- median(data_deg5$Strata_Age_Generations)
pred_seq1a$logGeneLen <- median(data_deg5$logGeneLen)
pred_seq1a$Species <- "Skylark" 
pred_seq1a$Wdegeneration <- "W_functional"
pred_seq1b <- pred_seq1a
pred_seq1b$Wdegeneration <- "W_loss_of_function"
pred_seq1c <- pred_seq1a
pred_seq1c$Species <- "Raso lark"
pred_seq1d <- pred_seq1b
pred_seq1d$Species <- "Raso lark"
pred_seq1a[,c("Wsel", "SE")] <- predict(FinalModel, newdata=pred_seq1a, type="response", se.fit=T)
pred_seq1a[,c("Wsel", "SE")] <- predict(FinalModel, newdata=pred_seq1b, type="response", se.fit=T)
pred_seq1a[,c("Wsel", "SE")] <- predict(FinalModel, newdata=pred_seq1c, type="response", se.fit=T)
pred_seq1a[,c("Wsel", "SE")] <- predict(FinalModel, newdata=pred_seq1d, type="response", se.fit=T)
pred_seq_zsel <- pred_seq1a
pred_seq_zsel$Species <- NA
pred_seq_zsel$Wdegeneration <- NA
pred_seq_zsel$Wsel <- rowMeans(cbind(pred_seq1a$Wsel, pred_seq1c$Wsel, pred_seq1b$Wsel, pred_seq1d$Wsel))
pred_seq_zsel$SE <- rowMeans(cbind(pred_seq1a$SE, pred_seq1c$SE, pred_seq1b$SE, pred_seq1d$SE))


#Plot as function of age and W degeneration
pred_seq1a <- as.data.frame(seq(min(data_deg5$Strata_Age_Generations)*100, max(data_deg5$Strata_Age_Generations)*100)/100)
colnames(pred_seq1a) <- "Strata_Age_Generations"
pred_seq1a$pHaplo <- median(data_deg5$pHaplo)
pred_seq1a$logGeneLen <- median(data_deg5$logGeneLen)
pred_seq1a$SexLinkedDNDS_DOSZ <- median(data_deg5$SexLinkedDNDS_DOSZ)
pred_seq1a$Species <- "Skylark" 
pred_seq1a$Wdegeneration <- "W_functional"
pred_seq1b <- pred_seq1a
pred_seq1b$Wdegeneration <- "W_loss_of_function"
pred_seq1c <- pred_seq1a
pred_seq1c$Species <- "Raso lark"
pred_seq1d <- pred_seq1b
pred_seq1d$Species <- "Raso lark"
pred_seq1a$Wsel <- predict(FinalModel, newdata=pred_seq1a, type="response")
pred_seq1b$Wsel <- predict(FinalModel, newdata=pred_seq1b, type="response")
pred_seq1c$Wsel <- predict(FinalModel, newdata=pred_seq1c, type="response")
pred_seq1d$Wsel <- predict(FinalModel, newdata=pred_seq1d, type="response")
pred_seq1ac <- pred_seq1a
pred_seq1ac$Species <- NA
pred_seq1ac$Wsel <- rowMeans(cbind(pred_seq1a$Wsel, pred_seq1c$Wsel))
pred_seq1bd <- pred_seq1b
pred_seq1bd$Species <- NA
pred_seq1bd$Wsel <- rowMeans(cbind(pred_seq1b$Wsel, pred_seq1d$Wsel))
pred_seq_age_wdeg <- rbind(pred_seq1ac, pred_seq1bd)
pred_seq_age_wdeg$Wsel <- 10^pred_seq_age_wdeg$Wsel

#Plot as function of Z selection, varying pHaplo
pred_seq1a <- as.data.frame(seq(min(data_deg5$SexLinkedDNDS_DOSZ)*100, max(data_deg5$SexLinkedDNDS_DOSZ)*100)/100)
colnames(pred_seq1a) <- "SexLinkedDNDS_DOSZ"
pred_seq1a$logGeneLen <- median(data_deg5$logGeneLen)
pred_seq1a$Strata_Age_Generations <- median(data_deg5$Strata_Age_Generations)
pred_seq1a$Species <- "Skylark" 
pred_seq1a$Wdegeneration <- "W_functional"
pred_seq1b <- pred_seq1a
pred_seq1b$Wdegeneration <- "W_loss_of_function"
pred_seq1c <- pred_seq1a
pred_seq1c$Species <- "Raso lark"
pred_seq1d <- pred_seq1b
pred_seq1d$Species <- "Raso lark"

pHaplo_data <- matrix(NA, nrow(pred_seq1a), 6)
colnames(pHaplo_data) <- c("0.0","0.2","0.4","0.6","0.8","1.0")

for(i in 1:length(colnames(pHaplo_data))) {
  pred_seq1a$pHaplo <- as.numeric(colnames(pHaplo_data)[i])
  pred_seq1b$pHaplo <- as.numeric(colnames(pHaplo_data)[i])
  pred_seq1c$pHaplo <- as.numeric(colnames(pHaplo_data)[i])
  pred_seq1d$pHaplo <- as.numeric(colnames(pHaplo_data)[i])
  
  Wsel1 <- predict(FinalModel, newdata=pred_seq1a, type="response")
  Wsel2 <- predict(FinalModel, newdata=pred_seq1b, type="response")
  Wsel3 <- predict(FinalModel, newdata=pred_seq1c, type="response")
  Wsel4 <- predict(FinalModel, newdata=pred_seq1d, type="response")
  
  pHaplo_data[,i] <- rowMeans(cbind(Wsel1, Wsel2, Wsel3, Wsel4))
}
pred_seq1 <- as.data.frame(cbind(pred_seq1a$SexLinkedDNDS_DOSZ, pHaplo_data))
colnames(pred_seq1)[1] <- "Zsel"
pred_seq_phap_zsel <- pred_seq1 |> pivot_longer(cols=-"Zsel", names_to="pHaplo", values_to= "Wsel")
pred_seq_phap_zsel$Wsel <- 10^pred_seq_phap_zsel$Wsel


#Plot as function of Z selection, varying gene length
pred_seq1a <- as.data.frame(seq(min(data_deg5$SexLinkedDNDS_DOSZ)*100, max(data_deg5$SexLinkedDNDS_DOSZ)*100)/100)
colnames(pred_seq1a) <- "SexLinkedDNDS_DOSZ"
pred_seq1a$pHaplo <- median(data_deg5$pHaplo)
pred_seq1a$Strata_Age_Generations <- median(data_deg5$Strata_Age_Generations)
pred_seq1a$Species <- "Skylark" 
pred_seq1a$Wdegeneration <- "W_functional"
pred_seq1b <- pred_seq1a
pred_seq1b$Wdegeneration <- "W_loss_of_function"
pred_seq1c <- pred_seq1a
pred_seq1c$Species <- "Raso lark"
pred_seq1d <- pred_seq1b
pred_seq1d$Species <- "Raso lark"

GenLen_data <- matrix(NA, nrow(pred_seq1a), 6)
colnames(GenLen_data) <- c("log10(500)", "log10(1000)",  "log10(2000)", "log10(5000)", "log10(10000)", "log10(15000)")

for(i in 1:length(colnames(GenLen_data))) {
  pred_seq1a$logGeneLen <- eval(parse(text=colnames(GenLen_data)[i]))
  pred_seq1b$logGeneLen <- eval(parse(text=colnames(GenLen_data)[i]))
  pred_seq1c$logGeneLen <- eval(parse(text=colnames(GenLen_data)[i]))
  pred_seq1d$logGeneLen <- eval(parse(text=colnames(GenLen_data)[i]))
  
  Wsel1 <- predict(FinalModel, newdata=pred_seq1a, type="response")
  Wsel2 <- predict(FinalModel, newdata=pred_seq1b, type="response")
  Wsel3 <- predict(FinalModel, newdata=pred_seq1c, type="response")
  Wsel4 <- predict(FinalModel, newdata=pred_seq1d, type="response")
  
  GenLen_data[,i] <- rowMeans(cbind(Wsel1, Wsel2, Wsel3, Wsel4))
}
pred_seq1 <- as.data.frame(cbind(pred_seq1a$SexLinkedDNDS_DOSZ, GenLen_data))
colnames(pred_seq1)[1] <- "Zsel"
pred_seq_len_zsel <- pred_seq1 |> pivot_longer(cols=-"Zsel", names_to="logGeneLen", values_to= "Wsel")
pred_seq_len_zsel$Wsel <- 10^pred_seq_len_zsel$Wsel


#Plots

#Supplemenary
plot_len_zsel_int <- ggplot() +
  geom_line(data=pred_seq_len_zsel, aes(x=Zsel, y=Wsel, group=logGeneLen), linewidth=1, linetype=1) +
  #geom_vline(aes(xintercept=0.5), color="#d7191c", linewidth=1, linetype=2) +
  scale_x_continuous(expand = c(0.05,0.05), limits=c(0,1), breaks=c(0.0, 0.25, 0.50, 0.75, 1.0)) +
  scale_y_continuous(expand = c(0,0), limits=c(0,1), breaks=c(0.0, 0.25, 0.50, 0.75, 1.0)) +
  labs(x=expression(atop("dN/(dN+dS)", "Z-gametologs")), y=expression(atop("dN/(dN+dS)", "W-gametologs")), title = NULL) +
  annotate(geom="text", x=rev(c(0.95, 0.95, 0.95, 0.95, 0.95, 0.95)), y=rev(c(0.27, 0.33, 0.395, 0.48, 0.56, 0.63)), label=rev(c("GL = 500 bp", "GL = 1 000 bp", "GL = 2 000 bp", "GL = 5 000 bp", "GL = 10 000 bp", "GL = 15 000 bp")), color="Black", size=4) +
  labs(color='pHaplo\n') +
  theme_bw() +
  theme(legend.position= "none",
        panel.spacing = unit(0, "lines"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.text.y.left =element_text(size=15, color="black"),
        axis.text.y.right = element_text(size=15, color="black"),
        axis.text.x = element_text(size=15, color="black"),
        axis.ticks.y.left = element_blank())

jpeg("Figures/Supplementary/gene_len_Zsel_Wsel.jpg", width=2000, height=2000, res=300)
plot_len_zsel_int
dev.off()

plot_phap_zsel_int <- ggplot() +
  geom_line(data=pred_seq_phap_zsel, aes(x=Zsel, y=Wsel, group=pHaplo), linewidth=1, linetype=1) +
  #geom_vline(aes(xintercept=0.5), color="#d7191c", linewidth=1, linetype=2) +
  scale_x_continuous(expand = c(0.05,0.05), limits=c(0,1), breaks=c(0.0, 0.25, 0.50, 0.75, 1.0)) +
  scale_y_continuous(expand = c(0,0), limits=c(0,1), breaks=c(0.0, 0.25, 0.50, 0.75, 1.0)) +
  labs(x=expression(atop("dN/(dN+dS)", "Z-gametologs")), y=expression(atop("dN/(dN+dS)", "W-gametologs")), title = NULL) +
  annotate(geom="text", x=rev(c(0.95, 0.95, 0.95, 0.95, 0.95, 0.95)), y=rev(c(0.27, 0.32, 0.37, 0.42, 0.47, 0.52)), label=rev(c("pHaplo = 0.0", "pHaplo = 0.2", "pHaplo = 0.4", "pHaplo = 0.6", "pHaplo = 0.8", "pHaplo = 1.0")), color="Black", size=4) +
  labs(color='pHaplo\n') +
  theme_bw() +
  theme(legend.position= "none",
        panel.spacing = unit(0, "lines"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.text.y.left =element_text(size=15, color="black"),
        axis.text.y.right = element_text(size=15, color="black"),
        axis.text.x = element_text(size=15, color="black"),
        axis.ticks.y.left = element_blank())

jpeg("Figures/Supplementary/gene_phap_Zsel_Wsel.jpg", width=2000, height=2000, res=300)
plot_phap_zsel_int
dev.off()


## Main plots
#W degeneration and Z selection

data_deg5$Wdegeneration2 <- factor(data_deg5$Wdegeneration, order=T, labels=c("W functional", "W loss-of-function mutation"),  levels=c("W_functional", "W_loss_of_function"))
pred_seq_age_wdeg$Wdegeneration <- factor(pred_seq_age_wdeg$Wdegeneration, order=T, labels=c("W functional", "W loss-of-function mutation"),  levels=c("W_functional", "W_loss_of_function"))

plot_zsel <- ggplot() +
  geom_point(data=data_deg5, aes(x=SexLinkedDNDS_DOSZ, y=SexLinkedDNDS_DOSW), size=1.5, alpha=0.4) +
  geom_ribbon(data=pred_seq_zsel, aes(x=SexLinkedDNDS_DOSZ, ymin=10^(Wsel-SE*1.96), ymax=10^(Wsel+SE*1.96)), alpha=0.1) +
  geom_line(data=pred_seq_zsel, aes(x=SexLinkedDNDS_DOSZ, y=10^Wsel), linewidth=2, linetype=1) +
 # geom_vline(aes(xintercept=0.5), color="#d7191c", linewidth=1, linetype=2) +
  scale_x_continuous(expand = c(0.01,0.01), limits=c(0,1), breaks=c(0.0, 0.25, 0.50, 0.75, 1.0)) +
  scale_y_continuous(expand = c(0,0), limits=c(0,1), breaks=c(0.0, 0.25, 0.50, 0.75, 1.0)) +
  labs(x=expression(atop("dN/(dN+dS)", "Z-gametologs")), y=expression(atop("dN/(dN+dS)", "W-gametologs")), title = NULL) +
   labs(color='pHaplo\n') +
  theme_bw() +
  theme(legend.position= "none",
        panel.spacing = unit(0, "lines"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.text.y.left =element_text(size=15, color="black"),
        axis.text.y.right = element_text(size=15, color="black"),
        axis.text.x = element_text(size=15, color="black"),
        axis.ticks.y.left = element_blank())


plot_age_wdeg_int <- ggplot() +
  geom_jitter(data=data_deg5, aes(x=Strata_Age_Generations, y=SexLinkedDNDS_DOSW, group=Wdegeneration2, color=Wdegeneration2), size=2, height=0, width=0.05) +
  geom_line(data=pred_seq_age_wdeg, aes(x=Strata_Age_Generations, y=Wsel, group=Wdegeneration, color=Wdegeneration), linewidth=2, linetype=1) +
  scale_color_manual(name="Functionality", values = c("W functional"="#E4EAF0", "W loss-of-function mutation"="#fecc5c")) +
  geom_hline(aes(yintercept=median(data_deg5$SexLinkedDNDS_DOSZ)), color="#404040", linewidth=1, linetype=2) +
  scale_x_continuous(expand = c(0,0), limits=c(2,8), breaks=c(2, 4, 6, 8)) +
  scale_y_continuous(expand = c(0,0), limits=c(0,1), breaks=c(0.0, 0.25, 0.50, 0.75, 1.0)) +
  annotate(geom="text", x=rev(sort(unique(data_deg5$Strata_Age_Generations))), y=c(0.71, 0.71, 0.71, 0.71), label=rev(c("3-c", "5", "3-b", "3-a\n&\n4A")), color="Black", size=5) +
  labs(x=expression(atop("Age (million generations)")), y=expression(atop("dN/(dN+dS)", "W-gametologs")), title = NULL) +
  theme_bw() +
  theme(legend.position= "none",
        panel.spacing = unit(0, "lines"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.text.y.left =element_text(size=15, color="black"),
        axis.text.y.right = element_text(size=15, color="black"),
        axis.text.x = element_text(size=15, color="black"),
        axis.ticks.y.left = element_blank())







effects_df$Variable <- c("Intercept", "log10(Gene length) *", "pHaplo *", "Z selection *", "Strata age", "W functionality (loss-of-function mutation)",
                         "log10(Gene length) * Z selection *",  "pHaplo * Z selection **", "pHaplo * Strata age", "Strata age * W functionality (loss-of-function mutation) **",
                         "Strata age * Z selection", "Species (Raso lark)", "Species (Raso lark) * Z selection", "log10(Gene length) * Strata age",
                         "Species (Raso lark) * W functionality (loss-of-function mutation)", "pHaplo * W functionality (loss-of-function mutation)", 
                         "log10(Gene length) * W functionality (loss-of-function mutation)", "Z selection * W functionality (loss-of-function mutation)", "log10(Gene length) * pHaplo")

effects_df$Variable <- factor(effects_df$Variable, order=T,
                              levels=c("Intercept", "log10(Gene length) *", "pHaplo *", "Species (Raso lark)", "Strata age", "Z selection *",  "W functionality (loss-of-function mutation)",
                                       "log10(Gene length) * pHaplo", "log10(Gene length) * Strata age", "log10(Gene length) * Z selection *", "log10(Gene length) * W functionality (loss-of-function mutation)",
                                       "pHaplo * Strata age", "pHaplo * Z selection **", "pHaplo * W functionality (loss-of-function mutation)", "Species (Raso lark) * Z selection",
                                       "Species (Raso lark) * W functionality (loss-of-function mutation)", "Strata age * Z selection", "Strata age * W functionality (loss-of-function mutation) **",
                                        "Z selection * W functionality (loss-of-function mutation)"))

effects_df <- effects_df[which(effects_df$Variable != "Intercept"),]
effects_df <- effects_df[order(effects_df$Variable),]
effects_df$contVar <- rev(1:nrow(effects_df))

plot_Wsel_effects <- ggplot(effects_df, aes(x = contVar, y = Std_Coefficien)) +
  geom_point(size = 3, color = "Black") +  # Plot odds ratios as points
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2, color = "black") +  # Add CI bars
  geom_hline(yintercept = 0, linewidth=1, linetype = 2, color = "#d7191c") +  # Reference line at 0
  scale_x_continuous(labels=effects_df$Variable, breaks=effects_df$contVar, sec.axis = dup_axis(name=NULL)) +
  scale_y_continuous(expand = c(0.0,0.0), limits=c(-0.2, 0.6), breaks=c(-0.2, 0, 0.2, 0.4, 0.6)) +
  coord_flip() +  # Flip coordinates for a horizontal plot
  labs(x = expression(atop("","Predictors")), y = expression(atop("Standardized coefficients (β)", ""))) +
  theme_bw()+
  theme(legend.position= "none",
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.text.y.left = element_blank(),
        axis.text.y.right = element_text(size=15, color="black"),
        axis.text.x = element_text(size=15, color="black"),
        axis.ticks.y.left = element_blank())

### Distribution of selection per strata
data_hist5 <- data_deg
data_hist5$SexLinkedBeta_DOS <- rep(NA, nrow(data_hist5))
data_hist5$SexLinkedAlpha_DOS <- rep(NA, nrow(data_hist5))
data_hist5$SexLinkedDNDS_DOS <- rep(NA, nrow(data_hist5))
data_hist5$Wdegeneration <- factor(data_hist5$Wdegeneration, order=F, labels=rev(c("Z-gametologs", "W functional", "W loss-of-function mutation", "W degenerated", "W degenerated")), levels=rev(c("Z-gametologs", "W functional", "W loss of function", "W partially degenerated", "W degenerated")))
data_hist5$Strata <-  factor(data_hist5$Strata, order=T, labels=rev(c("4A", "3-a", "3-b", "5", "3-c")), levels=rev(c("4A", "3-a", "3-b", "5", "3-c")))
data_hist5 <- data_hist5[which(data_hist5$Filter1=="OK"),]
data_hist5 <- data_hist5[which(!is.na(data_hist5$pHaplo)),]


data_hist5 <- data_hist5[which(data_hist5$Region != "autosomal" & data_hist5$Wdegeneration != "W degenerated" & data_hist5$Wdegeneration != "W partially degenerated"),]
data_hist5$SexLinkedDNDS_DOSZ <- data_hist5$SexLinkedBetaZ/(data_hist5$SexLinkedBetaZ + data_hist5$SexLinkedAlphaZ)
data_hist5$SexLinkedDNDS_DOSW <- data_hist5$SexLinkedBetaW/(data_hist5$SexLinkedBetaW + data_hist5$SexLinkedAlphaW)
data_hist5 <- data_hist5[which(data_hist5$SexLinkedAlphaW < 1 & data_hist5$SexLinkedAlphaW > 0.001 & data_hist5$SexLinkedBetaW < 1 & data_hist5$SexLinkedAlphaZ < 1 & data_hist5$SexLinkedAlphaZ > 0.001 & data_hist5$SexLinkedBetaZ < 1),]


data_hist5 <- data_hist5[which(data_hist5$SexLinkedDNDS_DOSW != 0 & data_hist5$SexLinkedDNDS_DOSZ != 0),]
data_hist5 <- data_hist5[which(!is.na(data_hist5$SexLinkedDNDS_DOSW) & !is.na(data_hist5$SexLinkedDNDS_DOSZ)),]

data_hist5$SexLinkedBeta_DOSZ <- data_hist5$SexLinkedBetaZ
data_hist5$SexLinkedBeta_DOSW <- data_hist5$SexLinkedBetaW
data_hist5$SexLinkedAlpha_DOSZ <- data_hist5$SexLinkedAlphaZ
data_hist5$SexLinkedAlpha_DOSW <- data_hist5$SexLinkedAlphaW

data_hist5 <- rbind(data_hist5, data_hist5)
data_hist5$SexLinkedDNDS_DOS <- rep(NA, nrow(data_hist5))
data_hist5$SexLinkedBeta_DOS <- rep(NA, nrow(data_hist5))
data_hist5$SexLinkedAlpha_DOS <- rep(NA, nrow(data_hist5))
data_hist5$SexLinkedDNDS_DOS[1:(nrow(data_hist5)/2)] <- data_hist5$SexLinkedDNDS_DOSZ[1:(nrow(data_hist5)/2)]
data_hist5$SexLinkedBeta_DOS[1:(nrow(data_hist5)/2)] <- data_hist5$SexLinkedBeta_DOSZ[1:(nrow(data_hist5)/2)]
data_hist5$SexLinkedAlpha_DOS[1:(nrow(data_hist5)/2)] <- data_hist5$SexLinkedAlpha_DOSZ[1:(nrow(data_hist5)/2)]
data_hist5$Wdegeneration[1:(nrow(data_hist5)/2)] <- "Z-gametologs"

data_hist5$SexLinkedDNDS_DOS[(1+nrow(data_hist5)/2):nrow(data_hist5)] <- data_hist5$SexLinkedDNDS_DOSW[(1+nrow(data_hist5)/2):nrow(data_hist5)]
data_hist5$SexLinkedBeta_DOS[(1+nrow(data_hist5)/2):nrow(data_hist5)] <- data_hist5$SexLinkedBeta_DOSW[(1+nrow(data_hist5)/2):nrow(data_hist5)]
data_hist5$SexLinkedAlpha_DOS[(1+nrow(data_hist5)/2):nrow(data_hist5)] <- data_hist5$SexLinkedAlpha_DOSW[(1+nrow(data_hist5)/2):nrow(data_hist5)]


#Plot W sel
medians1 <- data_hist5 |> group_by(Strata, Wdegeneration) |> summarize(median_value = median(SexLinkedDNDS_DOS), .groups = 'drop')
medians_Z1 <- medians1[which(medians1$Wdegeneration == "Z-gametologs"),] |> mutate(next_median_value = lead(median_value), next_Strata = lead(Strata)) |> filter(!is.na(next_median_value)) 
medians_Wf1 <- medians1[which(medians1$Wdegeneration == "W functional"),] |> mutate(next_median_value = lead(median_value), next_Strata = lead(Strata)) |> filter(!is.na(next_median_value)) 
medians_Wl1 <- medians1[which(medians1$Wdegeneration == "W loss-of-function mutation"),] |> mutate(next_median_value = lead(median_value), next_Strata = lead(Strata)) |> filter(!is.na(next_median_value)) 


plot_Wsel_strata <- ggplot() +
  geom_density_ridges(data=data_hist5, aes(x=SexLinkedDNDS_DOS, y=Wdegeneration, fill=Wdegeneration), stat="binline", bins=20, alpha=1, scale=2, draw_baseline = FALSE) +
  geom_point(data=medians1, aes(x=median_value, y=Wdegeneration), color="black", size=5) +
  geom_point(data=medians1, aes(x=median_value, y=Wdegeneration, color=Wdegeneration), size=4) +
  scale_fill_manual(name="Functionality", values = c("Z-gametologs"= "#404040", "W functional"="#E4EAF0", "W loss-of-function mutation"="#fecc5c")) +
  scale_color_manual(name="Functionality", values = c("Z-gametologs"= "#404040", "W functional"="#E4EAF0", "W loss-of-function mutation"="#fecc5c")) +
  labs(x=expression(atop("dN/(dN+dS)")), y=NULL, title = NULL) +
  scale_x_continuous(expand = c(0,0), limits=c(0, 1), breaks=c(0.0, 0.25, 0.50, 0.75, 1.0)) +
  facet_wrap(~Strata, nrow=5) +
  theme_bw() +
  theme(legend.position= "none",
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=20),
        axis.text.y.left = element_blank(),
        axis.text.y.right = element_text(size=12, color="black"),
        axis.text.x = element_text(size=15, color="black"),
        axis.ticks.y.left = element_blank())

#Plot W DN
medians2 <- data_hist5 |> group_by(Strata, Wdegeneration) |> summarize(median_value = median(SexLinkedBeta_DOS), .groups = 'drop')
medians_Z2 <- medians2[which(medians2$Wdegeneration == "Z-gametologs"),] |> mutate(next_median_value = lead(median_value), next_Strata = lead(Strata)) |> filter(!is.na(next_median_value)) 
medians_Wf2 <- medians2[which(medians2$Wdegeneration == "W functional"),] |> mutate(next_median_value = lead(median_value), next_Strata = lead(Strata)) |> filter(!is.na(next_median_value)) 
medians_Wl2 <- medians2[which(medians2$Wdegeneration == "W loss-of-function mutation"),] |> mutate(next_median_value = lead(median_value), next_Strata = lead(Strata)) |> filter(!is.na(next_median_value)) 


plot_Wsel_beta <- ggplot() +
  geom_density_ridges(data=data_hist5, aes(x=log10(SexLinkedBeta_DOS), y=Wdegeneration, fill=Wdegeneration), stat="binline", bins=20, alpha=1, scale=2, draw_baseline = FALSE) +
  geom_point(data=medians2, aes(x=log10(median_value), y=Wdegeneration), color="black", size=5) +
  geom_point(data=medians2, aes(x=log10(median_value), y=Wdegeneration, color=Wdegeneration), size=4) +
  scale_fill_manual(name="Functionality", values = c("Z-gametologs"= "#404040", "W functional"="#E4EAF0", "W loss-of-function mutation"="#fecc5c")) +
  scale_color_manual(name="Functionality", values = c("Z-gametologs"= "#404040", "W functional"="#E4EAF0", "W loss-of-function mutation"="#fecc5c")) +
  labs(x=expression(atop("log"[10]*"(dN)", "Z- and W-gametologs")), y=NULL, title = NULL) +
  scale_x_continuous(expand = c(0,0), limits=c(-3, 0)) +
  facet_wrap(~Strata, nrow=5, strip.position = "left") +
  theme_bw() +
  theme(legend.position= "none",
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=20),
        axis.text.y.left = element_blank(),
        axis.text.y.right = element_text(size=12, color="black"),
        axis.text.x = element_text(size=15, color="black"),
        axis.ticks.y.left = element_blank())


#Plot DS
medians3 <- data_hist5 |> group_by(Strata, Wdegeneration) |> summarize(median_value = median(SexLinkedAlpha_DOS), .groups = 'drop')
medians_Z3 <- medians3[which(medians3$Wdegeneration == "Z-gametologs"),] |> mutate(next_median_value = lead(median_value), next_Strata = lead(Strata)) |> filter(!is.na(next_median_value)) 
medians_Wf3 <- medians3[which(medians3$Wdegeneration == "W functional"),] |> mutate(next_median_value = lead(median_value), next_Strata = lead(Strata)) |> filter(!is.na(next_median_value)) 
medians_Wl3 <- medians3[which(medians3$Wdegeneration == "W loss-of-function mutation"),] |> mutate(next_median_value = lead(median_value), next_Strata = lead(Strata)) |> filter(!is.na(next_median_value)) 


plot_Wsel_alpha <- ggplot() +
  geom_density_ridges(data=data_hist5, aes(x=log10(SexLinkedAlpha_DOS), y=Wdegeneration, fill=Wdegeneration), stat="binline", bins=20, alpha=1, scale=2, draw_baseline = FALSE) +
  geom_point(data=medians3, aes(x=log10(median_value), y=Wdegeneration), color="black", size=5) +
  geom_point(data=medians3, aes(x=log10(median_value), y=Wdegeneration, color=Wdegeneration), size=4) +
  scale_fill_manual(name="Functionality", values = c("Z-gametologs"= "#404040", "W functional"="#E4EAF0", "W loss-of-function mutation"="#fecc5c"), breaks=c("Z-gametologs", "W functional", "W loss-of-function mutation")) +
  scale_color_manual(name="Functionality", values = c("Z-gametologs"= "#404040", "W functional"="#E4EAF0", "W loss-of-function mutation"="#fecc5c")) +
  labs(x=expression(atop("log"[10]*"(dS)")), y=NULL, title = NULL) +
  scale_x_continuous(expand = c(0,0), limits=c(-3, 0)) +
  guides(color="none") +
  facet_wrap(~Strata, nrow=5, strip.position = "left") +
  theme_bw() +
  theme(legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=20, face="bold"), #change legend title font size
        legend.text = element_text(size=20), #change legend text font size,
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.y.left = element_text(size=20, color="black", angle=0),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=20),
        axis.text.y.left = element_blank(),
        axis.text.y.right = element_text(size=12, color="black"),
        axis.text.x = element_text(size=15, color="black"),
        axis.ticks.y.left = element_blank())

fig5_plot <- (plot_Wsel_alpha | plot_Wsel_beta | plot_Wsel_strata | (plot_zsel / plot_age_wdeg_int / plot_Wsel_effects)) + plot_layout(guides = "collect", axis_titles = "collect", widths = c(1, 1, 1, 2.3)) +
  plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 20, face="bold"), plot.tag.position=c(-0.025,1), plot.margin=margin(20,20,20,20))


jpeg("Figures/Figure5.jpg", width=8500, height=7000, res=300)
fig5_plot
dev.off()
