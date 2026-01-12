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
data_deg$bin_function <- rep(NA, nrow(data_deg))
data_deg$bin_function[which(data_deg$Wdegeneration == "W functional")] <- 0
data_deg$bin_function[which(data_deg$Wdegeneration != "W functional")] <- 1
data_deg$bin_function <- factor(data_deg$bin_function, levels=c(0,1))

### Selection on Z as a function of degeneration and other predictors

# Remove and missing data
data_deg4 <- data_deg[which(data_deg$Filter1=="OK"),]
data_deg4 <- data_deg4[which(!is.na(data_deg4$pHaplo)),]

data_deg4 <- data_deg4[which(data_deg4$Strata != "Autosomal"),]
data_deg4$FullDNDS_DOS <- rep(NA, nrow(data_deg4))
data_deg4$FullDNDS_DOS[which(data_deg4$Strata == "PAR")] <- data_deg4$FullBetaA[which(data_deg4$Strata == "PAR")]/(data_deg4$FullBetaA[which(data_deg4$Strata == "PAR")] + data_deg4$FullAlphaA[which(data_deg4$Strata == "PAR")])
data_deg4$FullDNDS_DOS[which(data_deg4$Strata != "PAR")] <- data_deg4$FullBetaZ[which(data_deg4$Strata != "PAR")]/(data_deg4$FullBetaZ[which(data_deg4$Strata != "PAR")] + data_deg4$FullAlphaZ[which(data_deg4$Strata != "PAR")])
data_deg4 <- data_deg4[which(data_deg4$Strata != "PAR" | (data_deg4$FullAlphaA < 1 & data_deg4$FullAlphaA > 0.001 & data_deg4$FullBetaA < 1)),]
data_deg4 <- data_deg4[which(data_deg4$Strata == "PAR" | (data_deg4$FullAlphaZ < 1 & data_deg4$FullAlphaZ > 0.001 & data_deg4$FullBetaZ < 1)),]

# Transform data
hist(data_deg4$FullDNDS_DOS, breaks=1000)
data_deg4 <- data_deg4[which(!is.na(data_deg4$FullDNDS_DOS) & data_deg4$FullDNDS_DOS != 0),]
data_deg4$logFullDNDS_DOS <- log10(data_deg4$FullDNDS_DOS)
hist(data_deg4$FullDNDS_DOS, breaks=1000)
genes <- names(which(table(data_deg4$geneID) == 2))
data_deg4 <- data_deg4[data_deg4$geneID %in% genes, ]

hist(data_deg4$logFullDNDS_DOS)
qqnorm(data_deg4$logFullDNDS_DOS)
qqline(data_deg4$logFullDNDS_DOS)

hist(data_deg4$geneLengthDataBase)
qqnorm(data_deg4$geneLengthDataBase)
qqline(data_deg4$geneLengthDataBase)
data_deg4$logGeneLen <- log10(data_deg4$geneLengthDataBase)
hist(data_deg4$logGeneLen)
qqnorm(data_deg4$logGeneLen)
qqline(data_deg4$logGeneLen)

# Is the functional variable evenly distributed?
table(data_deg4$bin_function)
table(data_deg4$Wdegeneration)
data_deg4$Wdegeneration[which(data_deg4$Wdegeneration == "W partially degenerated")] <- "W degenerated"
data_deg4$Wdegeneration[which(data_deg4$Strata == "PAR")] <- "PAR"
data_deg4$Wdegeneration <- factor(data_deg4$Wdegeneration, order=F, labels=c("PAR", "W_functional", "W_loss_of_function", "W_degenerated"), levels=c("PAR", "W functional", "W loss of function", "W degenerated"))
table(data_deg4$Wdegeneration)
table(data_deg4$Strata)

### Test for colinearity of all factors
colin_model1 <- glmmTMB(logFullDNDS_DOS ~
                          pHaplo + logGeneLen + Species + Wdegeneration + Strata_Age_Generations,
                        family=gaussian, data = data_deg4, na.action = "na.fail", REML=F)
check_collinearity(colin_model1)

colin_model1 <- glmmTMB(logFullDNDS_DOS ~
                          pHaplo + logGeneLen + Species + Wdegeneration + Strata,
                        family=gaussian, data = data_deg4, na.action = "na.fail", REML=F)
check_collinearity(colin_model1)

colin_model1 <- glmmTMB(logFullDNDS_DOS ~
                        pHaplo + logGeneLen + Species + Wdegeneration,
                        family=gaussian, data = data_deg4, na.action = "na.fail", REML=F)
check_collinearity(colin_model1)

### Test model
globalmodel1 <- glmmTMB(logFullDNDS_DOS ~
                          Wdegeneration + pHaplo + logGeneLen + Species +
                          Wdegeneration:pHaplo + Wdegeneration:logGeneLen +
                          Wdegeneration:Species +
                          pHaplo:logGeneLen + pHaplo:Species +
                          logGeneLen:Species,
                        family = gaussian, data = data_deg4, na.action = "na.fail", REML=F)

options(na.action = "na.omit")
simulateResiduals(fittedModel = globalmodel1, plot = T) # OK

combinations1 <- dredge(global.model=globalmodel1, rank="AIC")

# Use globalmodel1
print(combinations1)
# Top 3 models have delta < 2 (max delta 1.52, cum weight 0.655)
0.333+0.167+0.155

# Check top 2 models with Dharma
simulateResiduals(fittedModel = get.models(combinations1, subset = 1)[[1]], plot = T) 
simulateResiduals(fittedModel = get.models(combinations1, subset = 2)[[1]], plot = T)
simulateResiduals(fittedModel = get.models(combinations1, subset = 3)[[1]], plot = T)

testQuantiles(simulateResiduals(fittedModel = get.models(combinations1, subset = 1)[[1]]))
testQuantiles(simulateResiduals(fittedModel = get.models(combinations1, subset = 2)[[1]]))
testQuantiles(simulateResiduals(fittedModel = get.models(combinations1, subset = 3)[[1]]))

# Get model average
FinalModel <- model.avg(get.models(combinations1, subset = delta <= 2))
summary(FinalModel)
summary(get.models(combinations1, subset = 1)[[1]])
summary(get.models(combinations1, subset = 2)[[1]])
summary(get.models(combinations1, subset = 3)[[1]])

#logGeneLen*
#pHaplo**
#pHaplo:WdegenerationW loss of function**
#pHaplo:WdegenerationW partially/fully degenerated**
#logGeneLen:pHaplo***


# Calculate Effect sizes
standardize_parameters(FinalModel)
r2(FinalModel)
r2(get.models(combinations1, subset = 1)[[1]])
r2(get.models(combinations1, subset = 2)[[1]])
r2(get.models(combinations1, subset = 3)[[1]])


terms <- sub("cond\\((.*)\\)", "\\1", colnames(FinalModel$coefficients))
terms[1] <- "(Intercept)"
ef_mat <- as.data.frame(matrix(0, ncol(FinalModel$coefficients), 3))
rownames(ef_mat) <- terms
eflow_mat <-as.data.frame(matrix(0, ncol(FinalModel$coefficients), 3))
rownames(eflow_mat) <- terms
efupp_mat <-as.data.frame(matrix(0, ncol(FinalModel$coefficients), 3))
rownames(efupp_mat) <- terms
weights <- FinalModel$msTable$weight

for(i in 1:3) {
  ef <- standardize_parameters(get.models(combinations1, subset = i)[[1]])
  for(j in 1:length(terms)) {
    if(length(which(ef$Parameter == terms[j])) > 0) {
      ef_mat[j,i] <- ef$Std_Coefficient[which(ef$Parameter == terms[j])]*weights[i]
      eflow_mat[j,i] <- ef$CI_low[which(ef$Parameter == terms[j])]*weights[i]
      efupp_mat[j,i] <- ef$CI_high[which(ef$Parameter == terms[j])]*weights[i]
    }
  }
}



effects_df <- as.data.frame(matrix(NA, length(terms), 4))
colnames(effects_df) <- c("Parameter", "Std_Coefficien", "CI_low", "CI_high")
for(i in 1:length(terms)) {
  effects_df$Parameter[i] <- terms[i]
  effects_df$Std_Coefficien[i] <- sum(ef_mat[i,1:3])
  effects_df$CI_low[i] <- sum(eflow_mat[i,1:3])
  effects_df$CI_high[i] <- sum(efupp_mat[i,1:3])
}


#Plot as function of pHaplo and different categories
ggpredict(FinalModel, terms = c("pHaplo", "Wdegeneration"))
pred_seq2a <- as.data.frame(seq(0, 1*100)/100)
colnames(pred_seq2a) <- "pHaplo"
pred_seq2a$logGeneLen <- median(data_deg4$logGeneLen)
pred_seq2a$Wdegeneration <- "PAR"
pred_seq2a$Species <- "Skylark"
pred_seq2b <- pred_seq2a
pred_seq2b$Species <- "Raso lark"
pred_seq2a_pred <- predict(FinalModel, newdata=pred_seq2a, type="response")
pred_seq2b_pred <- predict(FinalModel, newdata=pred_seq2b, type="response")
pred_seq2_int1 <- pred_seq2a
pred_seq2_int1$logFullDNDS_DOS <- rowMeans(cbind(pred_seq2a_pred, pred_seq2b_pred))
pred_seq2a$Wdegeneration <- "W_functional"
pred_seq2b$Wdegeneration <- "W_functional"
pred_seq2a_pred <- predict(FinalModel, newdata=pred_seq2a, type="response")
pred_seq2b_pred <- predict(FinalModel, newdata=pred_seq2b, type="response")
pred_temp <- cbind(pred_seq2a, rowMeans(cbind(pred_seq2a_pred, pred_seq2b_pred)))
colnames(pred_temp)[5] <- "logFullDNDS_DOS"
pred_seq2_int1 <- rbind(pred_seq2_int1, pred_temp)
pred_seq2a$Wdegeneration <- "W_loss_of_function"
pred_seq2b$Wdegeneration <- "W_loss_of_function"
pred_seq2a_pred <- predict(FinalModel, newdata=pred_seq2a, type="response")
pred_seq2b_pred <- predict(FinalModel, newdata=pred_seq2b, type="response")
pred_temp <- cbind(pred_seq2a, rowMeans(cbind(pred_seq2a_pred, pred_seq2b_pred)))
colnames(pred_temp)[5] <- "logFullDNDS_DOS"
pred_seq2_int1 <- rbind(pred_seq2_int1, pred_temp)
pred_seq2a$Wdegeneration <- "W_degenerated"
pred_seq2b$Wdegeneration <- "W_degenerated"
pred_seq2a_pred <- predict(FinalModel, newdata=pred_seq2a, type="response")
pred_seq2b_pred <- predict(FinalModel, newdata=pred_seq2b, type="response")
pred_temp <- cbind(pred_seq2a, rowMeans(cbind(pred_seq2a_pred, pred_seq2b_pred)))
colnames(pred_temp)[5] <- "logFullDNDS_DOS"
pred_seq2_int1 <- rbind(pred_seq2_int1, pred_temp)
pred_seq2_int1$FullDNDS_DOS <- 10^pred_seq2_int1$logFullDNDS_DOS


#Plot as function of Z selection, varying pHaplo
pred_seq1a <- as.data.frame(seq(min(data_deg4$logGeneLen)*100, max(data_deg4$logGeneLen)*100)/100)
colnames(pred_seq1a) <- "logGeneLen"
pred_seq1a$Wdegeneration <- rep(NA, nrow(pred_seq1a))
pred_seq1a$Species <- "Skylark"
pred_seq1b <- pred_seq1a
pred_seq1b$Species <- "Raso lark"

pHaplo_data <- as.data.frame(matrix(NA, nrow(pred_seq1a), 6))
colnames(pHaplo_data) <- c("0.0","0.2","0.4","0.6","0.8","1.0")

for(i in 1:length(colnames(pHaplo_data))) {
  pred_seq1a$pHaplo <- as.numeric(colnames(pHaplo_data)[i])
  pred_seq1b$pHaplo <- as.numeric(colnames(pHaplo_data)[i])
  pred_seq1a$Wdegeneration <- "PAR"
  pred_seq1b$Wdegeneration <- "PAR"
  Zsel1 <- predict(FinalModel, newdata=pred_seq1a, type="response")
  Zsel2 <- predict(FinalModel, newdata=pred_seq1b, type="response")
  
  pred_seq1a$Wdegeneration <- "W_functional"
  pred_seq1b$Wdegeneration <- "W_functional"
  Zsel1 <- cbind(Zsel1, predict(FinalModel, newdata=pred_seq1a, type="response"))
  Zsel2 <- cbind(Zsel2, predict(FinalModel, newdata=pred_seq1b, type="response"))
  
  pred_seq1a$Wdegeneration <- "W_loss_of_function"
  pred_seq1b$Wdegeneration <- "W_loss_of_function"
  Zsel1 <- cbind(Zsel1, predict(FinalModel, newdata=pred_seq1a, type="response"))
  Zsel2 <- cbind(Zsel2, predict(FinalModel, newdata=pred_seq1b, type="response"))
  
  pred_seq1a$Wdegeneration <- "W_degenerated"
  pred_seq1b$Wdegeneration <- "W_degenerated"
  Zsel1 <- cbind(Zsel1, predict(FinalModel, newdata=pred_seq1a, type="response"))
  Zsel2 <- cbind(Zsel2, predict(FinalModel, newdata=pred_seq1b, type="response"))
  
  pHaplo_data[,i] <- rowMeans(cbind(Zsel1, Zsel2))
}
pred_seq1 <- as.data.frame(cbind(pred_seq1a$logGeneLen, pHaplo_data))
colnames(pred_seq1)[1] <- "logGeneLen"
pred_seq_phap_len <- pred_seq1 |> pivot_longer(cols=-"logGeneLen", names_to="pHaplo", values_to= "Zsel")
pred_seq_phap_len$FullDNDS_DOS <- 10^pred_seq_phap_len$Zsel


#Plots

### Supplementary
plot_curves_len_phap <- ggplot() +
  geom_line(data=pred_seq_phap_len, aes(x=logGeneLen, y=FullDNDS_DOS, group=pHaplo), color="black", linewidth=1) +
  labs(x =expression("log"[10]*"(Gene length)"), y=expression(atop("dN/(dN+dS)", "PAR and Z-gametologs"))) +
  annotate(geom="text", x=c(3.90, 4.0, 4.05, 4.1, 4.15, 4.2) , y=c(0.33, 0.27, 0.215, 0.164, 0.125, 0.09), label=c("pHaplo = 0.0", "pHaplo = 0.2", "pHaplo = 0.4", "pHaplo = 0.6", "pHaplo = 0.8", "pHaplo = 1.0"), color="Black", size=4) +
  scale_x_continuous(expand = c(0.05,0.05), breaks=c(2.5, 3.0, 3.5, 4.0, 4.5)) +
  scale_y_continuous(expand = c(0,0), limits=c(0.0, 0.4), breaks=c(0.0, 0.1, 0.2, 0.3, 0.4)) +
  guides(color=guide_legend(title="Functionality")) +
  theme_bw() +
  theme(legend.position= "none",
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.text.y.left =element_text(size=15, color="black"),
        axis.text.y.right = element_text(size=15, color="black"),
        axis.text.x = element_text(size=15, color="black"),
        axis.ticks.y.left = element_blank())

jpeg("Figures/Supplementary/gene_len_pHap_Zsel.jpg", width=2000, height=2000, res=300)
plot_curves_len_phap
dev.off()

## Main plots
#W degeneration and Z selection


data_deg4$Wdegeneration2 <- factor(data_deg4$Wdegeneration, order=T, labels=c("PAR", "W functional", "W loss-of-function mutation", "W exon loss"),  levels=c("PAR", "W_functional", "W_loss_of_function", "W_degenerated"))

plot_Zsel <- ggplot() +
  geom_violin(data=data_deg4, aes(x=Wdegeneration2, y=FullDNDS_DOS, group=Wdegeneration2, fill=Wdegeneration2), color="black", width=0.8, position=position_dodge(1)) +
  geom_boxplot(data=data_deg4, aes(x=Wdegeneration2, y=FullDNDS_DOS, group=Wdegeneration2, fill=Wdegeneration2), color="black", width=0.1, position=position_dodge(1)) +
  scale_fill_manual(name="W degeneration", values = c("PAR"="#404040", "W functional"="#E4EAF0", "W exon loss"="#b30000", "W loss-of-function mutation"="#fecc5c")) +
  scale_y_continuous(expand = c(0,0), limits=c(0,1), breaks=c(0.0, 0.25, 0.50, 0.75, 1.0)) +
  labs(y=expression(atop("dN/(dN+dS)", "PAR and Z-gametologs")), x="Functionality", title = NULL) +
  guides(color=guide_legend(title="Functionality")) +
  theme_bw() +
  theme(legend.position= "none",
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.text.y.left =element_text(size=15, color="black"),
        axis.text.y.right = element_text(size=15, color="black"),
        axis.text.x = element_text(size=15, color="black", angle = 20, vjust=0.5),
        axis.ticks.y.left = element_blank())

pred_seq2_int1$Wdegeneration <- factor(pred_seq2_int1$Wdegeneration, order=F, labels=c("Autosomal/PAR", "W functional", "W loss-of-function mutation", "W exon loss"), levels=c("Autosomal/PAR", "W_functional", "W_loss_of_function", "W_degenerated"))
pred_seq2_int1$Wdegeneration[which(is.na(pred_seq2_int1$Wdegeneration))] <- "Autosomal/PAR"


plot_Zsel_interaction <- ggplot() +
  geom_line(data=pred_seq2_int1, aes(x=pHaplo, y=FullDNDS_DOS, group=Wdegeneration, color=Wdegeneration), linewidth=2, key_glyph = draw_key_rect) +
  scale_color_manual(name="W degeneration", values = c("Autosomal/PAR"="#404040", "W functional"="#E4EAF0", "W exon loss"="#b30000", "W loss-of-function mutation"="#fecc5c")) +
  labs(y=expression(atop("dN/(dN+dS)", "PAR & Z-gametologs")), title = NULL) +
  guides(color=guide_legend(title="Functionality"), override.aes = list(shape = 1)) +
  #scale_y_continuous(expand = c(0,0), limits=c(0.05, 0.3), breaks=c(0.1, 0.2, 0.3)) +
  scale_x_continuous(expand = c(0,0), breaks=c(0.0, 0.25, 0.50, 0.75, 1.0)) +
  theme_bw() +
  theme(legend.position= "none",
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.text.y.left =element_text(size=15, color="black"),
        axis.text.y.right = element_text(size=15, color="black"),
        axis.text.x = element_text(size=15, color="black"),
        axis.ticks.y.left = element_blank())



effects_df$Variable <- c("Intercept ***", "log10(Gene length) *", "pHaplo *", 
                         "W functionality (functional)", "W functionality (loss-of-function mutation)", "W functionality (exon loss)", "log10(Gene length) * pHaplo ***",
                         "log10(Gene length) * W functionality (functional)", "log10(Gene length) * W functionality (loss-of-function mutation)", "log10(Gene length) * W functionality (exon loss)", 
                         "pHaplo * W functionality (functional)", "pHaplo * W functionality (loss-of-function mutation) **", "pHaplo * W functionality (exon loss) **",
                         "Species (Raso lark)")

effects_df$Variable <- factor(effects_df$Variable, order=T,
                              levels=c("Intercept ***", "log10(Gene length) *", "pHaplo *",  "Species (Raso lark)",
                                       "W functionality (functional)", "W functionality (loss-of-function mutation)", "W functionality (exon loss)", "log10(Gene length) * pHaplo ***",
                                       "log10(Gene length) * W functionality (functional)", "log10(Gene length) * W functionality (loss-of-function mutation)", "log10(Gene length) * W functionality (exon loss)", 
                                       "pHaplo * W functionality (functional)", "pHaplo * W functionality (loss-of-function mutation) **", "pHaplo * W functionality (exon loss) **"))

effects_df <- effects_df[which(effects_df$Variable != "Intercept ***"),]
effects_df <- effects_df[order(effects_df$Variable),]
effects_df$contVar <- rev(1:nrow(effects_df))

plot_Zsel_effects <- ggplot(effects_df, aes(x = contVar, y = Std_Coefficien)) +
  geom_point(size = 3, color = "Black") +  # Plot odds ratios as points
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2, color = "black") +  # Add CI bars
  geom_hline(yintercept = 0, linewidth=1, linetype = 2, color = "#d7191c") +  # Reference line at 0
  scale_x_continuous(labels=effects_df$Variable, breaks=effects_df$contVar, sec.axis = dup_axis(name=NULL)) +
  scale_y_continuous(expand = c(0,0), limits=c(-0.5, 1.0), breaks=c(-0.5, 0, 0.5, 1.0)) +
  coord_flip() +  # Flip coordinates for a horizontal plot
  labs(x = expression(atop("","Predictors")), y = "Standardized coefficients (β)") +
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
data_hist4 <- data_deg
data_hist4 <- data_hist4[which(data_hist4$Strata != "Autosomal"),]
data_hist4$FullBeta_DOS <- rep(NA, nrow(data_hist4))
data_hist4$FullAlpha_DOS <- rep(NA, nrow(data_hist4))
data_hist4$FullDNDS_DOS <- rep(NA, nrow(data_hist4))
data_hist4$Wdegeneration <- factor(data_hist4$Wdegeneration, order=F, labels=rev(c("PAR", "W functional", "W loss-of-function mutation", "W exon loss", "W exon loss")), levels=rev(c("Autosomal/PAR", "W functional", "W loss of function", "W partially degenerated", "W degenerated")))
data_hist4$Strata <- factor(data_hist4$Strata, order=T, labels=rev(c("S0", "S1", "S2", "S3", "4A", "3-a", "3-b", "5", "3-c", "PAR")), levels=rev(c("S0", "S1", "S2", "S3", "4A", "3-a", "3-b", "5", "3-c", "PAR")))
data_hist4 <- data_hist4[which(data_hist4$Filter1=="OK"),]
data_hist4 <- data_hist4[which(!is.na(data_hist4$pHaplo)),]

data_hist4$Wdegeneration[which(data_hist4$Strata == "PAR")] <- "PAR"
data_hist4$FullDNDS_DOS[which(data_hist4$Strata == "PAR")] <- data_hist4$FullBetaA[which(data_hist4$Strata == "PAR")]/(data_hist4$FullBetaA[which(data_hist4$Strata == "PAR")] + data_hist4$FullAlphaA[which(data_hist4$Strata == "PAR")])
data_hist4$FullDNDS_DOS[which(data_hist4$Strata != "PAR")] <- data_hist4$FullBetaZ[which(data_hist4$Strata != "PAR")]/(data_hist4$FullBetaZ[which(data_hist4$Strata != "PAR")] + data_hist4$FullAlphaZ[which(data_hist4$Strata != "PAR")])
data_hist4 <- data_hist4[which(data_hist4$Strata != "PAR" | (data_hist4$FullAlphaA < 1 & data_hist4$FullAlphaA > 0.001 & data_hist4$FullBetaA < 1)),]
data_hist4 <- data_hist4[which(data_hist4$Strata == "PAR" | (data_hist4$FullAlphaZ < 1 & data_hist4$FullAlphaZ > 0.001 & data_hist4$FullBetaZ < 1)),]

data_hist4 <- data_hist4[which(!is.na(data_hist4$FullDNDS_DOS) & data_hist4$FullDNDS_DOS != 0),]
data_hist4 <- data_hist4[which(!is.na(data_hist4$FullDNDS_DOS)),]

data_hist4$FullBeta_DOS[which(data_hist4$Strata == "PAR")] <- data_hist4$FullBetaA[which(data_hist4$Strata == "PAR")]
data_hist4$FullBeta_DOS[which(data_hist4$Strata != "PAR")] <- data_hist4$FullBetaZ[which(data_hist4$Strata != "PAR")]
data_hist4$FullAlpha_DOS[which(data_hist4$Strata == "PAR")] <- data_hist4$FullAlphaA[which(data_hist4$Strata == "PAR")]
data_hist4$FullAlpha_DOS[which(data_hist4$Strata != "PAR")] <- data_hist4$FullAlphaZ[which(data_hist4$Strata != "PAR")]

#Plot Z sel
medians1 <- data_hist4 |> group_by(Strata, Wdegeneration) |> summarize(median_value = median(FullDNDS_DOS), .groups = 'drop')

plot_Zsel_strata <- ggplot() +
  geom_density_ridges(data=data_hist4, aes(x=FullDNDS_DOS, y=Wdegeneration, fill=Wdegeneration), stat="binline", bins=20, alpha=1, scale=2, draw_baseline = FALSE) +
  geom_vline(data=medians1, aes(xintercept=0.118), color="#404040", linetype=2, linewidth=2) +
  geom_point(data=medians1, aes(x=median_value, y=Wdegeneration), color="black", size=5) +
  geom_point(data=medians1, aes(x=median_value, y=Wdegeneration, color=Wdegeneration), size=4) +
  scale_fill_manual(name="Functionality", values = c("PAR"= "#404040", "W functional"="#E4EAF0", "W exon loss"="#b30000", "W loss-of-function mutation"="#fecc5c")) +
  scale_color_manual(name="Functionality", values = c("PAR"= "#404040", "W functional"="#E4EAF0", "W exon loss"="#b30000", "W loss-of-function mutation"="#fecc5c")) +
  labs(x=expression(atop("dN/(dN+dS)")), y=NULL, title = NULL) +
  scale_x_continuous(expand = c(0,0), limits=c(0, 1), breaks=c(0.0, 0.25, 0.50, 0.75, 1.0)) +
  facet_wrap(~Strata, nrow=10) +
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
medians2 <- data_hist4 |> group_by(Strata, Wdegeneration) |> summarize(median_value = median(FullBeta_DOS), .groups = 'drop')

plot_Zsel_beta <- ggplot() +
  geom_density_ridges(data=data_hist4, aes(x=log10(FullBeta_DOS), y=Wdegeneration, fill=Wdegeneration), stat="binline", bins=20, alpha=1, scale=2, draw_baseline = FALSE) +
  geom_vline(data=medians2, aes(xintercept=log10(0.0290)), color="#404040", linetype=2, linewidth=2) +
  geom_point(data=medians2, aes(x=log10(median_value), y=Wdegeneration), color="black", size=5) +
  geom_point(data=medians2, aes(x=log10(median_value), y=Wdegeneration, color=Wdegeneration), size=4) +
  scale_fill_manual(name="Functionality", values = c("PAR"= "#404040", "W functional"="#E4EAF0", "W exon loss"="#b30000", "W loss-of-function mutation"="#fecc5c")) +
  scale_color_manual(name="Functionality", values = c("PAR"= "#404040", "W functional"="#E4EAF0", "W exon loss"="#b30000", "W loss-of-function mutation"="#fecc5c")) +
  labs(x=expression(atop("log"[10]*"(dN)", "PAR and Z-gametologs")), y=NULL, title = NULL) +
  scale_x_continuous(expand = c(0,0), limits=c(-3.5, 0)) +
  facet_wrap(~Strata, nrow=10, strip.position = "left") +
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
medians3 <- data_hist4 |> group_by(Strata, Wdegeneration) |> summarize(median_value = median(FullAlpha_DOS), .groups = 'drop')

Strata_labels <- rev(c("S0"="S0", "S1"="S1", "S2"="S2", "S3"="S3", "4A"="4A", "3-a"="3-a", "3-b"="3-b", "5"="5", "3-c"="3-c", "PAR"="PAR"))

plot_Zsel_alpha <- ggplot() +
  geom_density_ridges(data=data_hist4, aes(x=log10(FullAlpha_DOS), y=Wdegeneration, fill=Wdegeneration), stat="binline", bins=20, alpha=1, scale=2, draw_baseline = FALSE) +
  geom_vline(data=medians3, aes(xintercept=log10(0.201)), color="#404040", linetype=2, linewidth=2) +
  geom_point(data=medians3, aes(x=log10(median_value), y=Wdegeneration), color="black", size=5) +
  geom_point(data=medians3, aes(x=log10(median_value), y=Wdegeneration, color=Wdegeneration), size=4) +
  scale_fill_manual(name="Functionality", values = c("PAR"= "#404040", "W functional"="#E4EAF0", "W exon loss"="#b30000", "W loss-of-function mutation"="#fecc5c"), breaks=c("PAR", "W functional", "W loss-of-function mutation", "W exon loss")) +
  scale_color_manual(name="Functionality", values = c("PAR"= "#404040", "W functional"="#E4EAF0", "W exon loss"="#b30000", "W loss-of-function mutation"="#fecc5c")) +
  labs(x=expression(atop("log"[10]*"(dS)")), y=NULL, title = NULL) +
  scale_x_continuous(expand = c(0,0), limits=c(-3.5, 0)) +
  guides(color="none") +
  facet_wrap(~Strata, nrow=10, strip.position = "left", labeller=as_labeller(Strata_labels)) +
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

fig4_plot <- (plot_Zsel_alpha | plot_Zsel_beta | plot_Zsel_strata | (plot_Zsel/plot_Zsel_interaction/plot_Zsel_effects)) + plot_layout(guides = "collect", axis_titles = "collect", widths = c(1, 1, 1, 2.5)) +
  plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 20, face="bold"), plot.tag.position=c(-0.025,1), plot.margin=margin(20,20,20,20))


jpeg("Figures/Figure4.jpg", width=8500, height=7000, res=300)
fig4_plot
dev.off()

