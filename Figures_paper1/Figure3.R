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
sessionInfo()

options(scipen=999)
setwd("C:/Users/Simon JE/OneDrive - Lund University/Dokument/Simon/PhD/Projects/Skylark_2021/Results/")
data <- read.delim("Genes/Skylark_2021_Rasolark_2021_organised_data1.tsv", sep="\t", head=T)


data_deg <- data[which(data$Filter3=="OK" & data$Filter4=="OK" & data$Filter5=="OK"),]
data_deg$Strata2 <- data_deg$Strata
data_deg$Strata[which(data_deg$Strata == "PAR3" | data_deg$Strata == "PAR5")] <- "PAR"
data_deg$Strata <- factor(data_deg$Strata, order=T, labels=rev(c("S0", "S1", "S2", "S3", "4A", "3-a", "3-b", "5", "3-c", "PAR", "Autosomal")), levels=rev(c("S0", "S1", "S2", "S3", "4A", "3a", "3b", "5", "3c", "PAR", "autosomal")))
data_deg$Species <- factor(data_deg$Species, order=T, labels=c("Skylark", "Raso lark"), levels=c("Skylark", "Rasolark"))

deg_prop <- data_deg |> count(Species, Strata, Wdegeneration)
deg_prop$prop <- rep(NA, nrow(deg_prop))

for(j in unique(deg_prop$Strata)) {
  deg_prop$prop[which(deg_prop$Strata == j & deg_prop$Species == "Raso lark")] <- deg_prop$n[which(deg_prop$Strata == j & deg_prop$Species == "Raso lark")] / sum(deg_prop$n[which(deg_prop$Strata == j & deg_prop$Species == "Raso lark")])
  deg_prop$prop[which(deg_prop$Strata == j & deg_prop$Species == "Skylark")] <- deg_prop$n[which(deg_prop$Strata == j & deg_prop$Species == "Skylark")] / sum(deg_prop$n[which(deg_prop$Strata == j & deg_prop$Species == "Skylark")])
}

deg_prop$Wdegeneration <- factor(deg_prop$Wdegeneration, order=T, labels=c("Functional", "Loss-of-function mutation", "Partial exon loss", "Full exon loss"), levels=c("W functional", "W loss of function", "W partially degenerated", "W degenerated"))


# Statistics
# Set up data
data_deg$bin_function <- rep(NA, nrow(data_deg))
data_deg$bin_function[which(data_deg$Wdegeneration == "W functional")] <- 0
data_deg$bin_function[which(data_deg$Wdegeneration != "W functional")] <- 1
data_deg$bin_function <- factor(data_deg$bin_function, levels=c(0,1))



### Is there any difference in functional genes between genomic regions?
# Set up data
data_deg2 <- data_deg[which(data_deg$Region == "autosomal"),c("geneID", "Region", "Strata2", "Species", "FullLOFA")]
colnames(data_deg2)[5] <- "bin_function"
data_deg2$Strata3 <- data_deg2$Strata2
data_deg2$Wdegeneration <- rep(NA, nrow(data_deg2))
data_deg2$Wdegeneration[which(data_deg2$bin_function == 0)] <- "Functional"
data_deg2$Wdegeneration[which(data_deg2$bin_function == 1)] <- "Loss of function"
temp1 <- data_deg[which(data_deg$Region != "autosomal"),c("geneID", "Region", "Strata2", "Species", "Wdegeneration", "bin_function")]
temp1$Strata3 <- paste("W-", temp1$Strata2, sep="")
temp1$Wdegeneration[which(temp1$Wdegeneration == "W functional")] <- "Functional"
temp1$Wdegeneration[which(temp1$Wdegeneration == "W loss of function" )] <- "Loss of function"
temp1$Wdegeneration[which(temp1$Wdegeneration == "W functional")] <- "Functional"
temp1$Wdegeneration[which(temp1$Wdegeneration == "W partially degenerated" )] <- "Partially degenerated"
temp1$Wdegeneration[which(temp1$Wdegeneration == "W degenerated" )] <- "Degenerated"
temp2 <- data_deg[which(data_deg$Region != "autosomal"),c("geneID", "Region", "Strata2", "Species", "Zdegeneration")]
temp2$bin_function <- rep(NA, nrow(temp2))
temp2$bin_function[which(temp2$Zdegeneration == "Z functional")] <- 0
temp2$bin_function[which(temp2$Zdegeneration != "Z functional")] <- 1
temp2$Strata3 <- paste("Z-", temp2$Strata2, sep="")
colnames(temp2)[5] <- "Wdegeneration"
temp2$Wdegeneration[which(temp2$Wdegeneration == "Z functional")] <- "Functional"
temp2$Wdegeneration[which(temp2$Wdegeneration == "Z loss of function" )] <- "Loss of function"
data_deg2 <- rbind(data_deg2, temp1, temp2)
data_deg2$bin_function <- factor(data_deg2$bin_function, levels=c(0,1))
data_deg2$Strata2 <- factor(data_deg2$Strata2, levels=unique(data_deg2$Strata2))
data_deg2$Strata3 <- factor(data_deg2$Strata3, levels=unique(data_deg2$Strata3))


### Check data with plot
deg_prop2 <- data_deg2 |> count(Species, Strata3, Wdegeneration)
deg_prop2$prop <- rep(NA, nrow(deg_prop2))

for(j in unique(deg_prop2$Strata3)) {
  deg_prop2$prop[which(deg_prop2$Strata3 == j & deg_prop2$Species == "Raso lark")] <- deg_prop2$n[which(deg_prop2$Strata3 == j & deg_prop2$Species == "Raso lark")] / sum(deg_prop2$n[which(deg_prop2$Strata3 == j & deg_prop2$Species == "Raso lark")])
  deg_prop2$prop[which(deg_prop2$Strata3 == j & deg_prop2$Species == "Skylark")] <- deg_prop2$n[which(deg_prop2$Strata3 == j & deg_prop2$Species == "Skylark")] / sum(deg_prop2$n[which(deg_prop2$Strata3 == j & deg_prop2$Species == "Skylark")])
}

deg_prop2$Strata3 <- factor(deg_prop2$Strata3, order=T, labels=rev(c("W-S0", "Z-S0", "W-S1", "Z-S1", "W-S2", "Z-S2", "W-S3", "Z-S3", "W-4A", "Z-4A", "W-3-a", "Z-3-a", "W-3-b", "Z-3-b", "W-5", "Z-5", "W-3-c", "Z-3-c", "PAR 5", "PAR 3", "Autosomal")),
                            levels=rev(c("W-S0", "Z-S0", "W-S1", "Z-S1", "W-S2", "Z-S2", "W-S3", "Z-S3", "W-4A", "Z-4A", "W-3a", "Z-3a", "W-3b", "Z-3b", "W-5", "Z-5", "W-3c", "Z-3c", "PAR5", "PAR3", "autosomal")))

deg_prop2$Wdegeneration <- factor(deg_prop2$Wdegeneration, order=T, labels=c("Functional", "Loss-of-function mutation", "Partial exon loss", "Full exon loss"), levels=c("Functional", "Loss of function", "Partially degenerated", "Degenerated"))


### Test when patterns in strata significantly changes from autosomes and respective Z-W stratum pair 
# Make non-empty valuesfor ancestral strata to avoid empty cells
data_deg2 <- rbind(data_deg2, c("dummy", "dummy", "S0", "Skylark", 0, "W-S0", "Functional"))
data_deg2 <- rbind(data_deg2, c("dummy", "dummy", "S0", "Raso lark", 0, "W-S0", "Functional"))
data_deg2 <- rbind(data_deg2, c("dummy", "dummy", "S1", "Skylark", 0, "W-S1", "Functional"))
data_deg2 <- rbind(data_deg2, c("dummy", "dummy", "S1", "Raso lark", 0, "W-S1", "Functional"))
data_deg2 <- rbind(data_deg2, c("dummy", "dummy", "S2", "Skylark", 0, "W-S2", "Functional"))
data_deg2 <- rbind(data_deg2, c("dummy", "dummy", "S2", "Raso lark", 0, "W-S2", "Functional"))

# Set up data
data_deg2$Species <- factor(data_deg2$Species, order=T, labels=c("Skylark", "Rasolark"), levels=c("Skylark", "Raso lark"))
data_deg2$Strata3 <- factor(data_deg2$Strata3, order=T, labels=rev(c("WS0", "ZS0", "WS1", "ZS1", "WS2", "ZS2", "WS3", "ZS3", "W4A", "Z4A", "W3a", "Z3a", "W3b", "Z3b", "W5", "Z5", "W3c", "Z3c", "PAR5", "PAR3", "autosomal")),
                            levels=rev(c("W-S0", "Z-S0", "W-S1", "Z-S1", "W-S2", "Z-S2", "W-S3", "Z-S3", "W-4A", "Z-4A", "W-3a", "Z-3a", "W-3b", "Z-3b", "W-5", "Z-5", "W-3c", "Z-3c", "PAR5", "PAR3", "autosomal")))

# Run model
globalmodel1 <- glmmTMB(bin_function ~
                          Strata3 * Species,
                        family=binomial, data = data_deg2, na.action = "na.fail", REML=F)

options(na.action = "na.omit")

combinations1 <- dredge(global.model=globalmodel1, rank="AIC")
print(combinations1)
# Top 1 models have delta < 2 (max delta 0, cum weight 1)
1

# Check top 1 models with Dharma
simulateResiduals(fittedModel = get.models(combinations1, subset = 1)[[1]], plot = T) # OK

# Get model
FinalModel <- glmmTMB(bin_function ~
                        Strata3 + Species,
                      family=binomial, data = data_deg2, na.action = "na.fail", REML=F)
Anova(FinalModel)

# Post hoc test
post_hoc <- glht(FinalModel, linfct=mcp(Strata3=c(
  "autosomal - PAR3 = 0", "autosomal - PAR5 = 0", "PAR3 - PAR5 = 0", 
  "autosomal - Z3c = 0", "autosomal - W3c = 0", "autosomal - Z5 = 0", "autosomal - W5 = 0",
  "autosomal - Z3b = 0", "autosomal - W3b = 0", "autosomal - Z3a = 0", "autosomal - W3a = 0",
  "autosomal - Z4A = 0", "autosomal - W4A = 0", "autosomal - ZS3 = 0", "autosomal - WS3 = 0",
  "autosomal - ZS2 = 0", "autosomal - WS2 = 0", "autosomal - ZS1 = 0", "autosomal - WS1 = 0",
  "autosomal - ZS0 = 0", "autosomal - WS0 = 0", "Z3c - W3c = 0", "Z5 - W5 = 0",
  "Z3b - W3b = 0", "Z3a - W3a = 0", "Z4A - W4A = 0", "ZS3 - WS3 = 0",
  "ZS2 - WS2 = 0", "ZS1 - WS1 = 0", "ZS0 - WS0 = 0")))

summary(post_hoc, test=adjusted("bonferroni"))
format(0.0, scientific = T)
r2(FinalModel)

### Degeneration as a function of age and other predictors

# Remove Autosomes and missing data
data_deg3 <- data_deg[which(data_deg$Strata != "Autosomal"),]
data_deg3 <- data_deg3[which(!is.na(data_deg3$pHaplo) & !is.na(data_deg3$HI)),]

# Transform data
## gene length in zebra finch. A gene cant have a negative length, so transformation makes sense. (model wont converge unless I transform this). 
hist(data_deg3$geneLengthDataBase)
qqnorm(data_deg3$geneLengthDataBase)
qqline(data_deg3$geneLengthDataBase)
data_deg3$logGeneLen <- log10(data_deg3$geneLengthDataBase)
hist(data_deg3$logGeneLen)
qqnorm(data_deg3$logGeneLen)
qqline(data_deg3$logGeneLen)

# Is the response variable evenly distributed?
table(data_deg3$bin_function)
table(data_deg3$Strata)

### Test for colinearity of all factors
colin_model1 <- glmmTMB(bin_function ~
                          Strata_Age_Generations + pHaplo + HI + logGeneLen + Species,
                        family=binomial, data = data_deg3, na.action = "na.fail", REML=F)
check_collinearity(colin_model1)

### Test for correlation between haploinsufficiency scores
anova(lm(pHaplo ~ HI, data = data_deg3))

### Test two different global models with combinations of haploinsufficiency as HI or pHaplo
# pHaplo
globalmodel1 <- glmmTMB(bin_function ~
                          Strata_Age_Generations + pHaplo + logGeneLen + Species +
                          Strata_Age_Generations:pHaplo + Strata_Age_Generations:logGeneLen +
                          Strata_Age_Generations:Species +
                          pHaplo:logGeneLen + pHaplo:Species +
                          logGeneLen:Species,
                        family=binomial, data = data_deg3, na.action = "na.fail", REML=F)

# %HI
globalmodel2 <- glmmTMB(bin_function ~
                          Strata_Age_Generations + HI + logGeneLen + Species +
                          Strata_Age_Generations:HI + Strata_Age_Generations:logGeneLen +
                          Strata_Age_Generations:Species +
                          HI:logGeneLen + HI:Species +
                          logGeneLen:Species,
                        family=binomial, data = data_deg3, na.action = "na.fail", REML=F)


options(na.action = "na.omit")

combinations1 <- dredge(global.model=globalmodel1, rank="AIC")
combinations2 <- dredge(global.model=globalmodel2, rank="AIC")


print(combinations1) # Top5 AIC 1889.6 - 1891.5, top5 models AIC delta < 2
print(combinations2) # Top5 AIC 1903.0 - 1904.4, top10 models AIC delta < 2
# Use pHaplo

# Remove Autosomes and missing data
data_deg3 <- data_deg[which(data_deg$Strata != "Autosomal"),]
data_deg3 <- data_deg3[which(!is.na(data_deg3$pHaplo)),]
data_deg3$logGeneLen <- log10(data_deg3$geneLengthDataBase)

# Fit model with pHaplo
globalmodel1 <- glmmTMB(bin_function ~
                          Strata_Age_Generations + pHaplo + logGeneLen + Species +
                          Strata_Age_Generations:pHaplo + Strata_Age_Generations:logGeneLen +
                          Strata_Age_Generations:Species +
                          pHaplo:logGeneLen + pHaplo:Species +
                          logGeneLen:Species,
                        family=binomial, data = data_deg3, na.action = "na.fail", REML=F)

options(na.action = "na.omit")

combinations1 <- dredge(global.model=globalmodel1, rank="AIC")

# Use globalmodel1
print(combinations1)
# Top 5 models have delta < 2 (max delta 1.93, cum weight 0.557)
0.191+0.099+0.098+0.096+0.073

# Check top 5 models with Dharma
simulateResiduals(fittedModel = get.models(combinations1, subset = 1)[[1]], plot = T) # OK
simulateResiduals(fittedModel = get.models(combinations1, subset = 2)[[1]], plot = T) # OK
simulateResiduals(fittedModel = get.models(combinations1, subset = 3)[[1]], plot = T) # OK
simulateResiduals(fittedModel = get.models(combinations1, subset = 4)[[1]], plot = T) # OK
simulateResiduals(fittedModel = get.models(combinations1, subset = 5)[[1]], plot = T) # OK

# Get model average
FinalModel <- model.avg(get.models(combinations1, subset = delta <= 2))
summary(FinalModel)
#logGeneLen***
#pHaplo
#Strata_Age_Generations***
#logGeneLen:pHaplo*
#logGeneLen:Strata_Age_Generations**
#pHaplo:Strata_Age_Generations***


r2(get.models(combinations1, subset = 1)[[1]])
r2(get.models(combinations1, subset = 2)[[1]])
r2(get.models(combinations1, subset = 3)[[1]])
r2(get.models(combinations1, subset = 4)[[1]])
r2(get.models(combinations1, subset = 5)[[1]])

# Get odds_ratios
coef <- summary(FinalModel)$coefmat.full
odds_ratios <- exp(coef[, "Estimate"])
odds_ratios_df <- data.frame(
  Variable = rownames(coef),
  Odds_Ratio = odds_ratios,
  Lower_CI = exp(coef[, "Estimate"] - 1.96 * coef[, "Std. Error"]),
  Upper_CI = exp(coef[, "Estimate"] + 1.96 * coef[, "Std. Error"])
)
print(odds_ratios_df)

#Plot as function of age
ggpredict(FinalModel, terms = "Strata_Age_Generations")
pred_seq1a <- as.data.frame(seq(0, max(data_deg3$Strata_Age_Generations)*100)/100)
colnames(pred_seq1a) <- "Strata_Age_Generations"
pred_seq1a$pHaplo <- median(data_deg3$pHaplo)
pred_seq1a$logGeneLen <- median(data_deg3$logGeneLen)
pred_seq1a$Species <- "Skylark"
pred_seq1b <- pred_seq1a
pred_seq1b$Species <- "Raso lark"
pred_se <- predict(FinalModel, newdata=pred_seq1a, type="response", se.fit=T)
pred_seq1a$LOF <- pred_se$fit
pred_seq1a$SE <- pred_se$se.fit
pred_se <- predict(FinalModel, newdata=pred_seq1b, type="response", se.fit=T)
pred_seq1b$LOF <- pred_se$fit
pred_seq1a$SE <- pred_se$se.fit
pred_seq1_age <- pred_seq1a
pred_seq1_age$LOF <- rowMeans(cbind(pred_seq1a$LOF, pred_seq1b$LOF))
pred_seq1_age$SE <- rowMeans(cbind(pred_seq1a$SE, pred_seq1b$SE))
pred_seq1_age$Species <- "NA"


pred_seq2a <- as.data.frame(unique(data_deg3$Strata_Age_Generations))
colnames(pred_seq2a) <- c("Strata_Age_Generations")
pred_seq2a$pHaplo <- median(data_deg3$pHaplo)
pred_seq2a$logGeneLen <- median(data_deg3$logGeneLen)
pred_seq2a$Species <- "Skylark"
pred_seq2b <- pred_seq2a
pred_seq2b$Species <- "Raso lark"
pred_seq2a$LOF <- predict(FinalModel, newdata=pred_seq2a, type="response")
pred_seq2b$LOF <- predict(FinalModel, newdata=pred_seq2b, type="response")
pred_seq2_age <- pred_seq2a
pred_seq2_age$LOF <- rowMeans(cbind(pred_seq2a$LOF, pred_seq2b$LOF))
pred_seq2_age$Species <- "NA"

pred_seq2_age$Strata <- rep(NA, nrow(pred_seq2_age))
pred_seq2_age$prop <- rep(NA, nrow(pred_seq2_age))
pred_seq2_age$Strata[1] <- unique(as.character(data_deg3$Strata[which(data_deg3$Strata_Age_Generations == pred_seq2_age$Strata_Age_Generations[1])]))
pred_seq2_age$prop[1] <- sum(deg_prop$prop[which(deg_prop$Strata == pred_seq2_age$Strata[1] & deg_prop$Wdegeneration != "Functional")])/2
pred_seq2_age$Strata[2] <- unique(as.character(data_deg3$Strata[which(data_deg3$Strata_Age_Generations == pred_seq2_age$Strata_Age_Generations[2])]))
pred_seq2_age$prop[2] <- sum(deg_prop$prop[which(deg_prop$Strata == pred_seq2_age$Strata[2] & deg_prop$Wdegeneration != "Functional")])/2
pred_seq2_age$Strata[3] <- unique(as.character(data_deg3$Strata[which(data_deg3$Strata_Age_Generations == pred_seq2_age$Strata_Age_Generations[3])]))
pred_seq2_age$prop[3] <- sum(deg_prop$prop[which(deg_prop$Strata == pred_seq2_age$Strata[3] & deg_prop$Wdegeneration != "Functional")])/2
pred_seq2_age$Strata[4] <- "3-a"

pred_seq2_age$prop[4] <-sum(deg_prop$prop[which(deg_prop$Strata == pred_seq2_age$Strata[4] & deg_prop$Wdegeneration != "Functional")])/2
pred_seq2_age <- rbind(pred_seq2_age, pred_seq2_age[4,])
pred_seq2_age$Strata[5] <- unique(as.character(data_deg3$Strata[which(data_deg3$Strata_Age_Generations == pred_seq2_age$Strata_Age_Generations[5])]))
pred_seq2_age$prop[5] <- sum(deg_prop$prop[which(deg_prop$Strata == pred_seq2_age$Strata[5] & deg_prop$Wdegeneration != "Functional")])/2
pred_seq2_age$Strata[6] <- unique(as.character(data_deg3$Strata[which(data_deg3$Strata_Age_Generations == pred_seq2_age$Strata_Age_Generations[6])]))
pred_seq2_age$prop[6] <- sum(deg_prop$prop[which(deg_prop$Strata == pred_seq2_age$Strata[6] & deg_prop$Wdegeneration != "Functional")])/2
pred_seq2_age$Strata[7] <- unique(as.character(data_deg3$Strata[which(data_deg3$Strata_Age_Generations == pred_seq2_age$Strata_Age_Generations[7])]))
pred_seq2_age$prop[7] <- sum(deg_prop$prop[which(deg_prop$Strata == pred_seq2_age$Strata[7] & deg_prop$Wdegeneration != "Functional")])/2
pred_seq2_age$Strata[8] <- unique(as.character(data_deg3$Strata[which(data_deg3$Strata_Age_Generations == pred_seq2_age$Strata_Age_Generations[8])]))
pred_seq2_age$prop[8] <- sum(deg_prop$prop[which(deg_prop$Strata == pred_seq2_age$Strata[8] & deg_prop$Wdegeneration != "Functional")])/2
pred_seq2_age$Strata[9] <- unique(as.character(data_deg3$Strata[which(data_deg3$Strata_Age_Generations == pred_seq2_age$Strata_Age_Generations[9])]))
pred_seq2_age$prop[9] <- sum(deg_prop$prop[which(deg_prop$Strata == pred_seq2_age$Strata[9] & deg_prop$Wdegeneration != "Functional")])/2
pred_seq2_age$Strata[10] <- "4A"
pred_seq2_age$prop[10] <-sum(deg_prop$prop[which(deg_prop$Strata == pred_seq2_age$Strata[10] & deg_prop$Wdegeneration != "Functional")])/2

intercept95_age <- pred_seq1_age$Strata_Age_Generations[which(abs(pred_seq1_age$LOF-0.95) == min(abs(pred_seq1_age$LOF-0.95)))]

age_gen_model <- readRDS("age_gen_model.RDS")
age_seq <- as.data.frame(seq(0, max(135.25)*100)/100)
colnames(age_seq) <- "cumAge"
age_seq$Strata_Age_Generations <- predict(age_gen_model, newdata=age_seq)
agefit <- lm(cumAge ~ -1 + Strata_Age_Generations + I(Strata_Age_Generations^2) + offset(rep(0, length(Strata_Age_Generations))), data=age_seq)
ticks <- as.data.frame(c(0,10,20,30))

intercept_year <- (intercept95_age * coef(agefit)[1]) + (coef(agefit)[2] * (intercept95_age ^2))

#Plot as function of gene length
ggpredict(FinalModel, terms = "logGeneLen")
pred_seq1a <- as.data.frame(seq(min(data_deg3$logGeneLen)*100, max(data_deg3$logGeneLen)*100)/100)
colnames(pred_seq1a) <- "logGeneLen"
pred_seq1a$pHaplo <- median(data_deg3$pHaplo)
pred_seq1a$Species <- "Skylark"
pred_seq1b <- pred_seq1a
pred_seq1b$Species <- "Raso lark"

# Age 0
pred_seq1a$Strata_Age_Generations <- 0
pred_seq1b$Strata_Age_Generations <- 0
pred_seq1a$LOF <- predict(FinalModel, newdata=pred_seq1a, type="response")
pred_seq1b$LOF <- predict(FinalModel, newdata=pred_seq1b, type="response")
pred_seq1_length0 <- pred_seq1a
pred_seq1_length0$LOF <- rowMeans(cbind(pred_seq1a$LOF, pred_seq1b$LOF))
pred_seq1_length0$Species <- "NA"

# Age 5
pred_seq1a$Strata_Age_Generations <- 5
pred_seq1b$Strata_Age_Generations <- 5
pred_seq1a$LOF <- predict(FinalModel, newdata=pred_seq1a, type="response")
pred_seq1b$LOF <- predict(FinalModel, newdata=pred_seq1b, type="response")
pred_seq1_length5 <- pred_seq1a
pred_seq1_length5$LOF <- rowMeans(cbind(pred_seq1a$LOF, pred_seq1b$LOF))
pred_seq1_length5$Species <- "NA"

# Age 10
pred_seq1a$Strata_Age_Generations <- 10
pred_seq1b$Strata_Age_Generations <- 10
pred_seq1a$LOF <- predict(FinalModel, newdata=pred_seq1a, type="response")
pred_seq1b$LOF <- predict(FinalModel, newdata=pred_seq1b, type="response")
pred_seq1_length10 <- pred_seq1a
pred_seq1_length10$LOF <- rowMeans(cbind(pred_seq1a$LOF, pred_seq1b$LOF))
pred_seq1_length10$Species <- "NA"

# Age 15
pred_seq1a$Strata_Age_Generations <- 15
pred_seq1b$Strata_Age_Generations <- 15
pred_seq1a$LOF <- predict(FinalModel, newdata=pred_seq1a, type="response")
pred_seq1b$LOF <- predict(FinalModel, newdata=pred_seq1b, type="response")
pred_seq1_length15 <- pred_seq1a
pred_seq1_length15$LOF <- rowMeans(cbind(pred_seq1a$LOF, pred_seq1b$LOF))
pred_seq1_length15$Species <- "NA"

# Age 20
pred_seq1a$Strata_Age_Generations <- 20
pred_seq1b$Strata_Age_Generations <- 20
pred_seq1a$LOF <- predict(FinalModel, newdata=pred_seq1a, type="response")
pred_seq1b$LOF <- predict(FinalModel, newdata=pred_seq1b, type="response")
pred_seq1_length20 <- pred_seq1a
pred_seq1_length20$LOF <- rowMeans(cbind(pred_seq1a$LOF, pred_seq1b$LOF))
pred_seq1_length20$Species <- "NA"

#Plot as function of pHaplo
ggpredict(FinalModel, terms = "pHaplo")
pred_seq1a <- as.data.frame(seq(min(data_deg3$pHaplo)*100, max(data_deg3$pHaplo)*100)/100)
colnames(pred_seq1a) <- "pHaplo"
pred_seq1a$logGeneLen <- median(data_deg3$logGeneLen)
pred_seq1a$Species <- "Skylark"
pred_seq1b <- pred_seq1a
pred_seq1b$Species <- "Raso lark"

# Age 0
pred_seq1a$Strata_Age_Generations <- 0
pred_seq1b$Strata_Age_Generations <- 0
pred_seq1a$LOF <- predict(FinalModel, newdata=pred_seq1a, type="response")
pred_seq1b$LOF <- predict(FinalModel, newdata=pred_seq1b, type="response")
pred_seq1_pHaplo0 <- pred_seq1a
pred_seq1_pHaplo0$LOF <- rowMeans(cbind(pred_seq1a$LOF, pred_seq1b$LOF))
pred_seq1_pHaplo0$Species <- "NA"

# Age 5
pred_seq1a$Strata_Age_Generations <- 5
pred_seq1b$Strata_Age_Generations <- 5
pred_seq1a$LOF <- predict(FinalModel, newdata=pred_seq1a, type="response")
pred_seq1b$LOF <- predict(FinalModel, newdata=pred_seq1b, type="response")
pred_seq1_pHaplo5 <- pred_seq1a
pred_seq1_pHaplo5$LOF <- rowMeans(cbind(pred_seq1a$LOF, pred_seq1b$LOF))
pred_seq1_pHaplo5$Species <- "NA"

# Age 10
pred_seq1a$Strata_Age_Generations <- 10
pred_seq1b$Strata_Age_Generations <- 10
pred_seq1a$LOF <- predict(FinalModel, newdata=pred_seq1a, type="response")
pred_seq1b$LOF <- predict(FinalModel, newdata=pred_seq1b, type="response")
pred_seq1_pHaplo10 <- pred_seq1a
pred_seq1_pHaplo10$LOF <- rowMeans(cbind(pred_seq1a$LOF, pred_seq1b$LOF))
pred_seq1_pHaplo10$Species <- "NA"

# Age 15
pred_seq1a$Strata_Age_Generations <- 15
pred_seq1b$Strata_Age_Generations <- 15
pred_seq1a$LOF <- predict(FinalModel, newdata=pred_seq1a, type="response")
pred_seq1b$LOF <- predict(FinalModel, newdata=pred_seq1b, type="response")
pred_seq1_pHaplo15 <- pred_seq1a
pred_seq1_pHaplo15$LOF <- rowMeans(cbind(pred_seq1a$LOF, pred_seq1b$LOF))
pred_seq1_pHaplo15$Species <- "NA"

# Age 20
pred_seq1a$Strata_Age_Generations <- 20
pred_seq1b$Strata_Age_Generations <- 20
pred_seq1a$LOF <- predict(FinalModel, newdata=pred_seq1a, type="response")
pred_seq1b$LOF <- predict(FinalModel, newdata=pred_seq1b, type="response")
pred_seq1_pHaplo20 <- pred_seq1a
pred_seq1_pHaplo20$LOF <- rowMeans(cbind(pred_seq1a$LOF, pred_seq1b$LOF))
pred_seq1_pHaplo20$Species <- "NA"

# When is probability of loss 0.5?
pred_seq1_age$Strata_Age_Generations[which(abs(pred_seq1_age$LOF-0.5) == min(abs(pred_seq1_age$LOF-0.5)))]
# set age to 7.5


#Plot pHaplo as a function of gene length
ggpredict(FinalModel, terms = "logGeneLen")
pred_seq1a <- as.data.frame(seq(min(data_deg3$logGeneLen)*100, max(data_deg3$logGeneLen)*100)/100)
colnames(pred_seq1a) <- "logGeneLen"
pred_seq1a$Strata_Age_Generations <- 7.5
pred_seq1a$Species <- "Skylark"
pred_seq1b <- pred_seq1a
pred_seq1b$Species <- "Raso lark"

#Plot as function of gene length, varying pHaplo
ggpredict(FinalModel, terms = "logGeneLen")
pred_seq1a <- as.data.frame(seq(min(data_deg3$logGeneLen)*100, max(data_deg3$logGeneLen)*100)/100)
colnames(pred_seq1a) <- "logGeneLen"
pred_seq1a$Strata_Age_Generations <- 7.5
pred_seq1a$Species <- "Skylark"
pred_seq1b <- pred_seq1a
pred_seq1b$Species <- "Raso lark"

# pHaplo 0
pred_seq1a$pHaplo <- 0
pred_seq1b$pHaplo <- 0
pred_seq1a$LOF <- predict(FinalModel, newdata=pred_seq1a, type="response")
pred_seq1b$LOF <- predict(FinalModel, newdata=pred_seq1b, type="response")
pred_seq1_len_phap0 <- pred_seq1a
pred_seq1_len_phap0$LOF <- rowMeans(cbind(pred_seq1a$LOF, pred_seq1b$LOF))
pred_seq1_len_phap0$Species <- "NA"

# pHaplo 0.2
pred_seq1a$pHaplo <- 0.2
pred_seq1b$pHaplo <- 0.2
pred_seq1a$LOF <- predict(FinalModel, newdata=pred_seq1a, type="response")
pred_seq1b$LOF <- predict(FinalModel, newdata=pred_seq1b, type="response")
pred_seq1_len_phap02 <- pred_seq1a
pred_seq1_len_phap02$LOF <- rowMeans(cbind(pred_seq1a$LOF, pred_seq1b$LOF))
pred_seq1_len_phap02$Species <- "NA"

# pHaplo 0.4
pred_seq1a$pHaplo <- 0.4
pred_seq1b$pHaplo <- 0.4
pred_seq1a$LOF <- predict(FinalModel, newdata=pred_seq1a, type="response")
pred_seq1b$LOF <- predict(FinalModel, newdata=pred_seq1b, type="response")
pred_seq1_len_phap04 <- pred_seq1a
pred_seq1_len_phap04$LOF <- rowMeans(cbind(pred_seq1a$LOF, pred_seq1b$LOF))
pred_seq1_len_phap04$Species <- "NA"

# pHaplo 0.6
pred_seq1a$pHaplo <- 0.6
pred_seq1b$pHaplo <- 0.6
pred_seq1a$LOF <- predict(FinalModel, newdata=pred_seq1a, type="response")
pred_seq1b$LOF <- predict(FinalModel, newdata=pred_seq1b, type="response")
pred_seq1_len_phap06 <- pred_seq1a
pred_seq1_len_phap06$LOF <- rowMeans(cbind(pred_seq1a$LOF, pred_seq1b$LOF))
pred_seq1_len_phap06$Species <- "NA"

# pHaplo 0.8
pred_seq1a$pHaplo <- 0.8
pred_seq1b$pHaplo <- 0.8
pred_seq1a$LOF <- predict(FinalModel, newdata=pred_seq1a, type="response")
pred_seq1b$LOF <- predict(FinalModel, newdata=pred_seq1b, type="response")
pred_seq1_len_phap08 <- pred_seq1a
pred_seq1_len_phap08$LOF <- rowMeans(cbind(pred_seq1a$LOF, pred_seq1b$LOF))
pred_seq1_len_phap08$Species <- "NA"

# pHaplo 1
pred_seq1a$pHaplo <- 1
pred_seq1b$pHaplo <- 1
pred_seq1a$LOF <- predict(FinalModel, newdata=pred_seq1a, type="response")
pred_seq1b$LOF <- predict(FinalModel, newdata=pred_seq1b, type="response")
pred_seq1_len_phap1 <- pred_seq1a
pred_seq1_len_phap1$LOF <- rowMeans(cbind(pred_seq1a$LOF, pred_seq1b$LOF))
pred_seq1_len_phap1$Species <- "NA"




#Plots
#Supplementary regions degeneration
plot_regions1 <- ggplot() +
  geom_bar(data=deg_prop2, aes(x=Species, y=prop, fill=Wdegeneration), stat="identity") +
  guides(fill=guide_legend(title="Functionality")) +
  scale_fill_manual(values=c("#E4EAF0", "#f09b20", "#f03b20", "#b30000")) +
  labs(x ="Genomic regions (Left: Skylark, Right: Raso lark)", y = "Proportion of genes") +
  facet_wrap(~Strata3, nrow=1) +
  scale_y_continuous(expand = c(0,0), limits=c(0,1)) +
  theme_bw() +
  theme(legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=12), #change legend title font size
        legend.text = element_text(size=12), #change legend text font size
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.spacing = unit(0, "lines"),
        strip.background = element_rect(color="black", fill="white", linewidth=1),
        strip.text = element_text(size=8, color="black", face="bold"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=20, color="black"),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank())

jpeg("Figures/Supplementary/degeneration_all_regions.jpg", width=4500, height=3000, res=300)
plot_regions1
dev.off()



plot_curves_len_phap <- ggplot() +
  geom_line(data=pred_seq1_len_phap0, aes(x=logGeneLen, y=LOF, color="Probability of\nnon-functionality"), linewidth=1) +
  geom_line(data=pred_seq1_len_phap02, aes(x=logGeneLen, y=LOF, color="Probability of\nnon-functionality"), linewidth=1) +
  geom_line(data=pred_seq1_len_phap04, aes(x=logGeneLen, y=LOF, color="Probability of\nnon-functionality"), linewidth=1) +
  geom_line(data=pred_seq1_len_phap06, aes(x=logGeneLen, y=LOF, color="Probability of\nnon-functionality"), linewidth=1) +
  geom_line(data=pred_seq1_len_phap08, aes(x=logGeneLen, y=LOF, color="Probability of\nnon-functionality"), linewidth=1) +
  geom_line(data=pred_seq1_len_phap1, aes(x=logGeneLen, y=LOF, color="Probability of\nnon-functionality"), linewidth=1) +
  scale_color_manual(values = "black", limits="Probability of\nnon-functionality") +
  scale_x_continuous(breaks=c(2, 2.5, 3, 3.5, 4, 4.5), limits=c(2, 4.5)) +
  scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0, 1)) +
  labs(x =expression("log"[10]*"(Gene length)"), y = expression(atop("Probability of", "non-functionality")), title = "Age = 7.5 million generations") +
  annotate(geom="text", x=c(3.0, 3.3, 3.45, 3.6, 3.8, 4) , y=c(0.9, 0.8, 0.7, 0.55, 0.4,  0.25), label=c("pHaplo = 0.0", "pHaplo = 0.2", "pHaplo = 0.4", "pHaplo = 0.6", "pHaplo = 0.8", "pHaplo = 1.0"), color="Black", size=4) +
  guides(color=guide_legend(title="Data")) +
  theme_bw() +
  theme(legend.position= "none",
        legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=15), #change legend title font size
        legend.text = element_text(size=15), #change legend text font size
        panel.spacing = unit(0, "lines"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=15, color="black"),
        axis.text.x = element_text(size=15, color="black"))
  
jpeg("Figures/Supplementary/gene_len_pHap_interaction.jpg", width=2000, height=2000, res=300)
plot_curves_len_phap
dev.off()



## Main plots
#W degeneratiom

deg_prop$Strata <- factor(deg_prop$Strata, order=T, labels=rev(c("W-S0", "W-S1", "W-S2", "W-S3", "W-4A", "W-3-a", "W-3-b", "W-5", "W-3-c", "PAR", "Autosomal")),
                          levels=rev(c("S0", "S1", "S2", "S3", "4A", "3-a", "3-b", "5", "3-c", "PAR", "Autosomal")))

plot_regions2 <- ggplot() +
  geom_bar(data=deg_prop, aes(x=Species, y=prop, fill=Wdegeneration), stat="identity") +
  guides(fill=guide_legend(title="Functionality")) +
  scale_fill_manual(values=c("#E4EAF0", "#f09b20", "#f03b20", "#b30000")) +
  labs(x ="Sex-chromosome strata (Left: Skylark, Right: Raso lark)", y = expression(atop("Proportion", "of genes"))) +
  facet_wrap(~Strata, nrow=1) +
  scale_y_continuous(expand = c(0,0), limits=c(0,1)) +
  theme_bw() +
  theme(legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=15), #change legend title font size
        legend.text = element_text(size=15), #change legend text font size
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.spacing = unit(0, "lines"),
        strip.background = element_rect(color="black", fill="white", linewidth=1),
        strip.text = element_text(size=10, color="black"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=15, color="black"),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank())


plot_curve_age <- ggplot() +
  geom_ribbon(data=pred_seq1_age, aes(x=Strata_Age_Generations, ymin=LOF-SE*1.96, ymax=LOF+SE*1.96), alpha=0.1) +
  geom_line(data=pred_seq1_age, aes(x=Strata_Age_Generations, y=LOF-SE*1.96), alpha=0.1, linewidth=1) +
  geom_line(data=pred_seq1_age, aes(x=Strata_Age_Generations, y=LOF+SE*1.96), alpha=0.1, linewidth=1) +
  geom_line(data=pred_seq1_age, aes(x=Strata_Age_Generations, y=LOF, color="Probability of\nnon-functionality"), linewidth=1) +
  geom_point(data=pred_seq2_age, aes(x=Strata_Age_Generations, y=prop, color="Strata (Proportions)"), size=3.5) +
  geom_vline(aes(xintercept=intercept95_age, color="95% probability of\nnon-functionality"), linewidth=1, linetype=2) +
  scale_color_manual(values = c("black", "#2c7bb6", "#d7191c"), limits=c("Probability of\nnon-functionality", "Strata (Proportions)", "95% probability of\nnon-functionality")) +
  scale_x_continuous(sec.axis=sec_axis(LOF~ (. * coef(agefit)[1]) + (coef(agefit)[2] * (. ^2)), name="Age (million years in larks)", breaks=seq(0,140,20)), breaks=c(0,5,10,15,20,25,30)) +
  scale_y_continuous(limits=c(0,1)) +
  geom_text(data = pred_seq2_age, aes(x=Strata_Age_Generations, y=prop, label = Strata, hjust = c(-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5), vjust = c(1.5,1.5,1.5,0.8,1.5,1.5,-0.5,1.5,1.5,1.5)), color = "Black") +
  labs(x ="Age (million generations)", y = expression(atop("Probability of", "non-functionality")), title = NULL) +
  guides(color=guide_legend(title="Data")) +
  theme_bw() +
  theme(legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=15), #change legend title font size
        legend.text = element_text(size=15), #change legend text font size
        panel.spacing = unit(0, "lines"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=15, color="black"),
        axis.text.x = element_text(size=15, color="black"))

plot_curves_len <- ggplot() +
  geom_line(data=pred_seq1_length0, aes(x=logGeneLen, y=LOF, color="Probability of\nnon-functionality"), linewidth=1) +
  geom_line(data=pred_seq1_length5, aes(x=logGeneLen, y=LOF, color="Probability of\nnon-functionality"), linewidth=1) +
  geom_line(data=pred_seq1_length10, aes(x=logGeneLen, y=LOF, color="Probability of\nnon-functionality"), linewidth=1) +
  geom_line(data=pred_seq1_length15, aes(x=logGeneLen, y=LOF, color="Probability of\nnon-functionality"), linewidth=1) +
  geom_line(data=pred_seq1_length20, aes(x=logGeneLen, y=LOF, color="Probability of\nnon-functionality"), linewidth=1) +
  scale_color_manual(values = "black", limits="Probability of\nnon-functionality") +
  scale_x_continuous(breaks=c(2, 2.5, 3, 3.5, 4, 4.5), limits=c(2, 4.5)) +
  scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0, 1)) +
  labs(x =expression("log"[10]*"(Gene length)"), y = expression(atop("Probability of", "non-functionality")), title = NULL) +
  annotate(geom="text", x=c(4.25, 3.9, 3.25, 2.8, 2.25) , y=c(0.16, 0.36, 0.65, 0.90, 0.97), label=c("Age = 0 MG", "Age = 5 MG", "Age = 10 MG", "Age = 15 MG", "Age = 20 MG"), color="Black", size=4) +
  guides(color=guide_legend(title="Data")) +
  theme_bw() +
  theme(legend.position= "none",
        legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=15), #change legend title font size
        legend.text = element_text(size=15), #change legend text font size
        panel.spacing = unit(0, "lines"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=15, color="black"),
        axis.text.x = element_text(size=15, color="black"))

plot_curves_phap <- ggplot() +
  geom_line(data=pred_seq1_pHaplo0, aes(x=pHaplo, y=LOF, color="Probability of\nnon-functionality"), linewidth=1) +
  geom_line(data=pred_seq1_pHaplo5, aes(x=pHaplo, y=LOF, color="Probability of\nnon-functionality"), linewidth=1) +
  geom_line(data=pred_seq1_pHaplo10, aes(x=pHaplo, y=LOF, color="Probability of\nnon-functionality"), linewidth=1) +
  geom_line(data=pred_seq1_pHaplo15, aes(x=pHaplo, y=LOF, color="Probability of\nnon-functionality"), linewidth=1) +
  geom_line(data=pred_seq1_pHaplo20, aes(x=pHaplo, y=LOF, color="Probability of\nnon-functionality"), linewidth=1) +
  geom_vline(aes(xintercept=-1, color="95% Probability of\nnon-functionality"), linewidth=1, linetype=2) +
  scale_color_manual(values = "black", limits="Probability of\nnon-functionality") +
  scale_x_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0, 1)) +
  scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0, 1)) +
  labs(x = "pHaplo (probability)", y = expression(atop("Probability of", "non-functionality")), title = NULL) +
  annotate(geom="text", x=c(0.39, 0.52, 0.65, 0.78, 0.92) , y=c(0.01, 0.19, 0.55, 0.78, 0.90), label=c("Age = 0 MG", "Age = 5 MG", "Age = 10 MG", "Age = 15 MG", "Age = 20 MG"), color="Black", size=4) +
  guides(color=guide_legend(title="Data")) +
  theme_bw() +
  theme(legend.position= "none",
        legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=15), #change legend title font size
        legend.text = element_text(size=15), #change legend text font size
        panel.spacing = unit(0, "lines"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=15, color="black"),
        axis.text.x = element_text(size=15, color="black"))




odds_ratios_df$Variable <- c("Intercept ***", "log10(Gene length) ***", "pHaplo", "Species (Raso lark)", "Strata age ***",
                             "log10(Gene length) * pHaplo *", "log10(Gene length) * Strata age **", "pHaplo * Strata age ***", 
                             "Species (Raso lark) * Strata age", "log10(Gene length) * Species (Raso lark)", "pHaplo * Species (Raso lark)")



odds_ratios_df$Variable <- factor(odds_ratios_df$Variable, order=T,
                                  levels=c("Intercept ***", "log10(Gene length) ***", "pHaplo", "Species (Raso lark)", "Strata age ***",
                                           "log10(Gene length) * pHaplo *", "log10(Gene length) * Species (Raso lark)", 
                                           "log10(Gene length) * Strata age **", "pHaplo * Species (Raso lark)",
                                           "pHaplo * Strata age ***", "Species (Raso lark) * Strata age"))


odds_ratios_df$logOdds_Ratio <- log10(odds_ratios_df$Odds_Ratio)
odds_ratios_df$logLower_CI <- log10(odds_ratios_df$Lower_CI)
odds_ratios_df$logUpper_CI <- log10(odds_ratios_df$Upper_CI)
odds_ratios_df <- odds_ratios_df[which(odds_ratios_df$Variable != "Intercept ***"),]
odds_ratios_df <- odds_ratios_df[order(odds_ratios_df$Variable),]
odds_ratios_df$contVar <- rev(1:nrow(odds_ratios_df))



plot_OR <- ggplot(odds_ratios_df, aes(x = contVar, y = logOdds_Ratio)) +
  geom_point(size = 3, color = "Black") +  # Plot odds ratios as points
  geom_errorbar(aes(ymin = logLower_CI, ymax = logUpper_CI), width = 0.2, color = "black") +  # Add CI bars
  geom_hline(yintercept = log10(1), linewidth=1, linetype = 2, color = "#d7191c") +  # Reference line at OR = 0
  scale_x_continuous(labels=odds_ratios_df$Variable, breaks=odds_ratios_df$contVar, sec.axis = dup_axis(name=NULL)) +
  coord_flip() +  # Flip coordinates for a horizontal plot
  labs(x = expression(atop("","Predictors")), y = expression("log"[10]*"(Odds Ratio)")) +
  theme_bw()+
  theme(legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=15), #change legend title font size
        legend.text = element_text(size=15), #change legend text font size
        panel.spacing = unit(0, "lines"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.text.y.left = element_text(size=1, color="white"),
        axis.text.y.right = element_text(size=12, color="black"),
        axis.text.x = element_text(size=15, color="black"),
        axis.ticks.y.left = element_blank())


fig3_plot <- plot_regions2 / plot_curve_age / (plot_curves_len + plot_curves_phap + plot_OR) +
plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 20, face="bold"), plot.tag.position=c(0,0.95))
# Positions of subfigure tags + legends need to be adjusted manually

jpeg("Figures/Figure3.jpg", width=6000, height=3500, res=300)
fig3_plot
dev.off()
