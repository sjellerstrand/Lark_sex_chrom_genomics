## Export variables and load libraries
rm(list=ls())
library(tidyverse)
library(patchwork)
library(performance)
library(brms)
library(glmmTMB)
library(DHARMa)
library(ggeffects)
sessionInfo()

setwd("C:/Users/Simon JE/OneDrive - Lund University/Dokument/Simon/PhD/Projects/")

### Import evolutionary data
data_evo <- read.delim("Skylark_2021/Results/Genes/Skylark_2021_Rasolark_2021_organised_data2.tsv", sep="\t", head=T)
GERP <- read.delim("Rasolark_2021/Results/GERP/GERP_scores_genes.tsv", sep="\t", head=T)

### Import rescaling to generations
age_gen_model <- readRDS("Skylark_2021/Results/age_gen_model.RDS")
age_seq <- as.data.frame(seq(0, max(135.25)*100)/100)
colnames(age_seq) <- "cumAge"
age_seq$Strata_Age_Generations <- predict(age_gen_model, newdata=age_seq)
agefit <- lm(cumAge ~ -1 + Strata_Age_Generations + I(Strata_Age_Generations^2) + offset(rep(0, length(Strata_Age_Generations))), data=age_seq)

### Filter evolutionary data
data_evo <- data_evo[which(data_evo$Filter3=="OK" & data_evo$Filter4=="OK" & data_evo$Filter5=="OK"),]
data_evo$Strata2 <- data_evo$Strata
data_evo$Strata[which(data_evo$Strata == "PAR3" | data_evo$Strata == "PAR5")] <- "PAR"
data_evo$Strata <- factor(data_evo$Strata, order=T, labels=rev(c("S0", "S1", "S2", "S3", "4A", "3-a", "3-b", "5", "3-c", "PAR", "Autosomal")), levels=rev(c("S0", "S1", "S2", "S3", "4A", "3a", "3b", "5", "3c", "PAR", "autosomal")))
data_evo$Species <- factor(data_evo$Species, order=T, labels=c("Skylark", "Raso lark"), levels=c("Skylark", "Rasolark"))

### Import genotype data
data_gt <- read.delim("Rasolark_2021/Results/Genotypes/genotype_info.txt", sep="\t", head=T)
gene_strata <- unique(cbind(data_evo$geneID, as.character(data_evo$Strata2)))
gene_strata <- gene_strata[which(gene_strata[,2] != "Autosomal"), ]
for(i in 1:nrow(gene_strata)) {
  data_gt$Strata[which(data_gt$geneID == gene_strata[i,1])] <- gene_strata[i,2] 
}
data_gt$Heterozygosity_score_LOF[which(is.na(data_gt$Heterozygosity_score_LOF))] <- 0
data_gt$Heterozygosity_score_HIGH[which(is.na(data_gt$Heterozygosity_score_HIGH))] <- 0

### Filter genotype data
data_gt <- data_gt %>% filter(geneID %in% data_evo$geneID)
data_evo <- data_evo %>% filter(geneID %in% data_gt$geneID)
data_gt$Strata2 <- data_gt$Strata
data_gt$Strata[which(data_gt$Strata == "PAR3" | data_gt$Strata == "PAR5")] <- "PAR"
data_gt$Strata <- factor(data_gt$Strata, order=T, labels=rev(c("S0", "S1", "S2", "S3", "4A", "3-a", "3-b", "5", "3-c", "PAR", "Autosomal")), levels=rev(c("S0", "S1", "S2", "S3", "4A", "3a", "3b", "5", "3c", "PAR", "autosomal")))
data_gt$Species <- factor(str_split(data_gt$Project, "_", simplify=T)[,1], labels=c("Raso lark", "Skylark"), levels=c("Rasolark", "Skylark"))
data_gt$Sex <-factor(data_gt$Sex, levels=c("Female", "Male"))

# Add gene specific data
data_gt$pHaplo <- rep(NA, nrow(data_gt))
data_gt$geneLengthDataBase <- rep(NA, nrow(data_gt))
data_gt$Strata_Age_Generations <- rep(NA, nrow(data_gt))
data_gt$GERP_n <- rep(NA, nrow(data_gt))
data_gt$GERP_mean <- rep(NA, nrow(data_gt))
data_gt$GERP_sd <- rep(NA, nrow(data_gt))
data_gt$GERP_max <- rep(NA, nrow(data_gt))

unique_genes <- data_evo[which(data_evo$Species == "Skylark"),]

for(i in 1:nrow(unique_genes)) {
  indices <- which(data_gt$geneID == unique_genes$geneID[i])
  data_gt$Strata_Age_Generations[indices] <- unique_genes$Strata_Age_Generations[i]
  data_gt$geneLengthDataBase[indices] <- unique_genes$geneLengthDataBase[i]
  data_gt$pHaplo[indices] <- unique_genes$pHaplo[i]
  indices2 <- which(GERP$geneID == unique_genes$geneID[i])
  data_gt$GERP_n[indices] <- GERP$n[indices2]
  data_gt$GERP_mean[indices] <- GERP$mean[indices2]
  data_gt$GERP_sd[indices] <- GERP$sd[indices2]
  data_gt$GERP_max[indices] <- GERP$max[indices2]
}
rm(unique_genes)
hist(data_gt$GERP_n, breaks=1000, xlim=c(0,1000))
data_gt <- data_gt[which(data_gt$GERP_n > 100),]
hist(data_gt$GERP_n, breaks=1000, xlim=c(0,1000))
hist(data_gt$GERP_max, breaks=100)
hist(data_gt$GERP_mean, breaks=100)
plot(data_gt$GERP_max, data_gt$GERP_mean)
data_gt$logGeneLen <- log10(data_gt$geneLengthDataBase)
hist(data_gt$logGeneLen)
data_gt <- data_gt[which(!is.na(data_gt$pHaplo)),]

# Set up genotypes
data_gt$func_hom_score_LOF <- rep(0, nrow(data_gt))
data_gt$func_hom_score_LOF[which(data_gt$Heterozygosity_score_LOF == 0 & data_gt$LOF1 == 0 & data_gt$LOF2 == 0 & data_gt$W_Dropout == 0)] <- 1
data_gt$LOF_hom_score_LOF <- rep(0, nrow(data_gt))
data_gt$LOF_hom_score_LOF[which((data_gt$LOF1 == 1 & (data_gt$LOF2 == 1 & data_gt$W_Dropout == 0)) | (data_gt$LOF1 == 1 & data_gt$LOF2 == 0 & data_gt$W_Dropout == 1))] <- 1
data_gt$Mask_Het_score_LOF <- rep(0, nrow(data_gt))
data_gt$Mask_Het_score_LOF[which((data_gt$LOF1 == 0 & (data_gt$LOF2 == 1 | data_gt$W_Dropout == 1)) | (data_gt$LOF1 == 1 & data_gt$LOF2 == 0 & data_gt$W_Dropout == 0))] <- 1
data_gt$Func_Het_score_LOF <- data_gt$Heterozygosity_score_LOF
data_gt$Tot_Het_score_LOF <- rep(0, nrow(data_gt))
data_gt$Tot_Het_score_LOF[which(data_gt$Func_Het_score_LOF == 1 | data_gt$Mask_Het_score_LOF == 1)] <- 1





# Statistically test predictors for genotypes
data_gt2 <- data_gt[which(data_gt$Strata != "Autosomal"),]
data_gt2 <- data_gt2[,c("geneID", "Sex", "Species", "Individual", "Strata", "Strata_Age_Generations", "logGeneLen", "pHaplo", "GERP_max", "LOF_hom_score_LOF", "Mask_Het_score_LOF", "Func_Het_score_LOF", "func_hom_score_LOF", "Tot_Het_score_LOF")]
data_gt2long <- data_gt2 |> group_by(geneID, Sex, Species, Individual, Strata) |>
  pivot_longer(cols = c(func_hom_score_LOF, LOF_hom_score_LOF, Mask_Het_score_LOF, Func_Het_score_LOF, Tot_Het_score_LOF), names_to = "variable", values_to = "value")

data_gt2long <- data_gt2long[which(data_gt2long$value == 1 & data_gt2long$variable != "Tot_Het_score_LOF"),]

# Sanity check
any(duplicated(data_gt2long[, c("geneID", "Individual")]))
nrow(data_gt2long)
nrow(unique(data_gt2long[, c("geneID", "Individual")]))

data_gt2long$dummy <- rbinom(nrow(data_gt2long), 1, 0.5)

# Set reference categories
data_gt2long$variable <- factor(data_gt2long$variable, labels=c("func_hom", "LOF_hom", "Mask_Het", "Func_Het"), levels=c("func_hom_score_LOF", "LOF_hom_score_LOF", "Mask_Het_score_LOF", "Func_Het_score_LOF"))

# Check for collinearity
vif_model <- glm(dummy ~ Species + Sex + Strata_Age_Generations + GERP_max + pHaplo + logGeneLen,
                 family = binomial, data = data_gt2long)
check_collinearity(vif_model)

#### Decide reference category of responses (and other)
Model2 <- brm(variable ~
                Species + Sex + Strata_Age_Generations + GERP_max + pHaplo + logGeneLen +
                Species:Sex + Species:Strata_Age_Generations + Species:GERP_max + Species:pHaplo + Species:logGeneLen +
                Sex:Strata_Age_Generations + Sex:GERP_max + Sex:pHaplo + Sex:logGeneLen +
                Strata_Age_Generations:GERP_max + Strata_Age_Generations:pHaplo + Strata_Age_Generations:logGeneLen +
                GERP_max:pHaplo + GERP_max:logGeneLen +
                pHaplo:logGeneLen,
              family=categorical(), data = data_gt2long, cores = 4, file = "Rasolark_2021/Results/Model2_saved")
summary(Model2)
max(rhat(Model2))
neff_ratio(Model2)
min(summary(Model2)$fixed$Bulk_ESS)
min(summary(Model2)$fixed$Tail_ESS)
#plot(Model2)
pp_check(Model2, type = "bars")
posterior_summary(Model2)


#Posterior probability of direction
post <- posterior_samples(Model2)
colnames(post) <- sub("^b_", "", colnames(post))
terms <- colnames(post)
summary_df <- post |>
  summarise(across(everything(), list(mean = ~ mean(.x), lower = ~ quantile(.x, 0.025), upper = ~ quantile(.x, 0.975)))) |>
  pivot_longer(everything(), names_to = c("Variable", ".value"), names_pattern = "(.*)_(.*)")
ppd_df <- sapply(terms, function(term) {
  prop_pos <- mean(post[[term]] > 0)
  prop_neg <- mean(post[[term]] < 0)
  max(prop_pos, prop_neg)
})
ppd_df <- data.frame(Variable = names(ppd_df), PPD = as.numeric(ppd_df))
pdir_df <- sapply(terms, function(term) {
  prop_pos <- mean(post[[term]] > 0)
  2 * min(prop_pos, 1 - prop_pos)})
pdir_df <- data.frame(Variable = names(pdir_df), p_dir = as.numeric(pdir_df))

summary_df2 <- cbind(rownames(summary(Model2)$fixed), summary(Model2)$fixed[,c("Rhat", "Bulk_ESS", "Tail_ESS")])
colnames(summary_df2)[1] <- "Variable"

summary_df <- summary_df |>
  left_join(ppd_df, by = "Variable") |>
  left_join(pdir_df, by = "Variable") |>
  left_join(summary_df2, by = "Variable") |>
  arrange(Variable)
print(summary_df)

# Get odds_ratios
odds_ratios2 <- exp(posterior_summary(Model2)[, "Estimate"])
odds_ratios_df2 <- data.frame(
  Variable = rownames(posterior_summary(Model2)),
  Odds_Ratio = odds_ratios2,
  Lower_CI = exp(posterior_summary(Model2)[, "Q2.5"]),
  Upper_CI = exp(posterior_summary(Model2)[, "Q97.5"])
)

odds_ratios_df2$Variable <- sub("^b_", "", odds_ratios_df2$Variable)

print(odds_ratios_df2)
odds_ratios_df2 <- odds_ratios_df2 |>
  left_join(summary_df, by = "Variable") |>
  arrange(Variable)

# Significance
odds_ratios_df2$sig <- rep(NA, nrow(odds_ratios_df2))
odds_ratios_df2$sig[which(odds_ratios_df2$p_dir < 0.1)] <- "."
odds_ratios_df2$sig[which(odds_ratios_df2$p_dir < 0.05)] <- "*"
odds_ratios_df2$sig[which(odds_ratios_df2$p_dir < 0.01)] <- "**"
odds_ratios_df2$sig[which(odds_ratios_df2$p_dir < 0.001)] <- "***"

# Modify names and columns
odds_ratios_df2$data <- str_split(odds_ratios_df2$Variable, "_", simplify=T)[,1]
odds_ratios_df2 <- odds_ratios_df2[which(nchar(odds_ratios_df2$data) > 0),]
odds_ratios_df2$terms <- str_split(odds_ratios_df2$Variable, "_", simplify=T)[,2]
odds_ratios_df2 <- odds_ratios_df2[-which(odds_ratios_df2$Variable == "Intercept_muFuncHet" | odds_ratios_df2$Variable == "Intercept_muLOFhom" | odds_ratios_df2$Variable == "Intercept_muMaskHet" |
                                            odds_ratios_df2$Variable == "lp__" | odds_ratios_df2$Variable == "lprior"),]

odds_ratios_df2$data[which(odds_ratios_df2$data == "muLOFhom")] <- "Homozygous loss-of-function (LL)"
odds_ratios_df2$data[which(odds_ratios_df2$data == "muMaskHet")] <- "Heterozygous masking (FL)"
odds_ratios_df2$data[which(odds_ratios_df2$data == "muFuncHet")] <- "Heterozygosity functional (F'F'')"

odds_ratios_df2$logOdds_Ratio <- log10(odds_ratios_df2$Odds_Ratio)
odds_ratios_df2$logLower_CI <- log10(odds_ratios_df2$Lower_CI)
odds_ratios_df2$logUpper_CI <- log10(odds_ratios_df2$Upper_CI)
odds_ratios_df2 <- select(odds_ratios_df2, -c("Odds_Ratio", "Lower_CI", "Upper_CI", "mean", "lower", "upper"))

odds_ratios_df2$Variable2 <- rep(NA, nrow(odds_ratios_df2))

odds_ratios_df2$Variable2[which(odds_ratios_df2$data == "Homozygous loss-of-function (LL)")] <-
  c("Maximum GERP score", "Maximum GERP score × log10(Gene length)", "Maximum GERP score × pHaplo ***", "Intercept ***", "Sex (Male) ***",
    "Sex (Male) × Maximum GERP score ***",  "Sex (Male) × Strata age ***", "Sex (Male) × log10(Gene length) ***", "Sex (Male) × pHaplo ***", "Species (Skylark)",
    "Species (Skylark) × Maximum GERP score ***", "Species (Skylark) × Sex (Male) ***",  "Species (Skylark) × Strata age", "Species (Skylark) × log10(Gene length) ***",
    "Species (Skylark) × pHaplo **", "Strata age ***", "Strata age × Maximum GERP score ***", "Strata age × log10(Gene length) ***", "Strata age × pHaplo ***", "log10(Gene length) ***",
    "pHaplo ***", "pHaplo × log10(Gene length)")

odds_ratios_df2$Variable2[which(odds_ratios_df2$data == "Heterozygous masking (FL)")] <-
  c("Maximum GERP score ***", "Maximum GERP score × log10(Gene length) *", "Maximum GERP score × pHaplo ***", "Intercept **", "Sex (Male) ***",
    "Sex (Male) × Maximum GERP score",  "Sex (Male) × Strata age ***", "Sex (Male) × log10(Gene length)", "Sex (Male) × pHaplo **", "Species (Skylark)",
    "Species (Skylark) × Maximum GERP score **", "Species (Skylark) × Sex (Male) ***",  "Species (Skylark) × Strata age", "Species (Skylark) × log10(Gene length) ***",
    "Species (Skylark) × pHaplo *", "Strata age ***", "Strata age × Maximum GERP score ***", "Strata age × log10(Gene length) ***", "Strata age × pHaplo", "log10(Gene length) **",
    "pHaplo ***", "pHaplo × log10(Gene length) *")

odds_ratios_df2$Variable2[which(odds_ratios_df2$data == "Heterozygosity functional (F'F'')")] <-
  c("Maximum GERP score", "Maximum GERP score × log10(Gene length)", "Maximum GERP score × pHaplo ***", "Intercept ***", "Sex (Male) ***",
    "Sex (Male) × Maximum GERP score",  "Sex (Male) × Strata age ***", "Sex (Male) × log10(Gene length) ***", "Sex (Male) × pHaplo ***", "Species (Skylark) *",
    "Species (Skylark) × Maximum GERP score ***", "Species (Skylark) × Sex (Male) ***",  "Species (Skylark) × Strata age ***", "Species (Skylark) × log10(Gene length) ***",
    "Species (Skylark) × pHaplo **", "Strata age ***", "Strata age × Maximum GERP score *", "Strata age × log10(Gene length) **", "Strata age × pHaplo", "log10(Gene length) ***",
    "pHaplo *", "pHaplo × log10(Gene length) *")

write.table(odds_ratios_df2, "Rasolark_2021/Results/Model2_stats.tsv", sep="\t")

# Calculate posterior probability contrasts
age_seq <- seq(min(data_gt2long$Strata_Age_Generations), max(data_gt2long$Strata_Age_Generations), length.out = 7)[c(2,4,6)]
newdata <- expand.grid(
  Species = unique(data_gt2long$Species),
  Sex = unique(data_gt2long$Sex),
  Strata_Age_Generations = age_seq,
  GERP_max = median(data_gt2long$GERP_max, na.rm = TRUE),
  pHaplo = median(data_gt2long$pHaplo, na.rm = TRUE),
  logGeneLen = median(data_gt2long$logGeneLen, na.rm = TRUE))
newdata$id <- 1:nrow(newdata)
post_probs <- posterior_epred(Model2, newdata = newdata)
post_df <- as.data.frame.table(post_probs)
colnames(post_df) <- c("iteration", "row", "genotype", "prob")
post_df$row <- as.integer(post_df$row)
post_df <- post_df |> left_join(newdata, by = c("row" = "id"))

# Sex contrasts (Female - Male) within species and age
sex_in_species_diff <- post_df |>
  select(iteration, Species, Sex, Strata_Age_Generations, genotype, prob) |>
  pivot_wider(names_from = Sex, values_from = prob) |>
  mutate(female_minus_male = Female - Male)
sex_in_species_diff <- sex_in_species_diff |>
  group_by(Species, Strata_Age_Generations, genotype) |>
  summarise(
    mean_diff  = mean(female_minus_male),
    lower_diff = quantile(female_minus_male, 0.025),
    upper_diff = quantile(female_minus_male, 0.975),
    PPD        = max(mean(female_minus_male > 0), mean(female_minus_male < 0)),
    p_dir      = 2 * min(mean(female_minus_male > 0), mean(female_minus_male < 0)),
    .groups = "drop")

sex_in_species_diff <- sex_in_species_diff |> arrange(Species, genotype, Strata_Age_Generations)
print(sex_in_species_diff, n=100)

# Species contrasts (Raso lark - Skylark) within sex and age
species_in_sex_diff <- post_df |>
  select(iteration, Species, Sex, Strata_Age_Generations, genotype, prob) |>
  pivot_wider(names_from = Species, values_from = prob) |>
  mutate(raso_minus_skylark = `Raso lark` - Skylark)
species_in_sex_diff <- species_in_sex_diff |>
  group_by(Sex, Strata_Age_Generations, genotype) |>
  summarise(
    mean_diff  = mean(raso_minus_skylark),
    lower_diff = quantile(raso_minus_skylark, 0.025),
    upper_diff = quantile(raso_minus_skylark, 0.975),
    PPD        = max(mean(raso_minus_skylark > 0), mean(raso_minus_skylark < 0)),
    p_dir      = 2 * min(mean(raso_minus_skylark > 0), mean(raso_minus_skylark < 0)),
    .groups = "drop")
species_in_sex_diff <- species_in_sex_diff |> arrange(Sex, genotype, Strata_Age_Generations)
print(species_in_sex_diff, n=100)

# Sex contrasts (Female - Male) across species (Raso lark - Skylark) and age
sex_diff_wide <- post_df |>
  select(iteration, Species, Sex, Strata_Age_Generations, genotype, prob) |>
  pivot_wider(names_from = Sex, values_from = prob) |>
  mutate(female_minus_male = Female - Male)
sex_diff_interaction <- sex_diff_wide |>
  select(iteration, Species, Strata_Age_Generations, genotype, female_minus_male) |>
  pivot_wider(names_from = Species, values_from = female_minus_male) |>
  mutate(diff_of_diffs = `Raso lark` - Skylark)
sex_diff_interaction_summary <- sex_diff_interaction |>
  group_by(Strata_Age_Generations, genotype) |>
  summarise(
    mean_diff  = mean(diff_of_diffs),
    lower_diff = quantile(diff_of_diffs, 0.025),
    upper_diff = quantile(diff_of_diffs, 0.975),
    PPD        = max(mean(diff_of_diffs > 0), mean(diff_of_diffs < 0)),
    p_dir      = 2 * min(mean(diff_of_diffs > 0), mean(diff_of_diffs < 0)),
    .groups = "drop")
sex_diff_interaction_summary <- sex_diff_interaction_summary |> arrange(genotype, Strata_Age_Generations)
print(sex_diff_interaction_summary, n=100)

# Plot all probabilities as a function of age
age_seq_fine <- seq(min(data_gt2long$Strata_Age_Generations), max(data_gt2long$Strata_Age_Generations), length.out = 100)
newdata_fine <- expand.grid(
  Species = unique(data_gt2long$Species),
  Sex = unique(data_gt2long$Sex),
  Strata_Age_Generations = age_seq_fine,
  GERP_max = median(data_gt2long$GERP_max, na.rm = TRUE),
  pHaplo = median(data_gt2long$pHaplo, na.rm = TRUE),
  logGeneLen = median(data_gt2long$logGeneLen, na.rm = TRUE))
newdata_fine$id <- 1:nrow(newdata_fine)
post_probs <- posterior_epred(Model2, newdata = newdata_fine)
post_df <- as.data.frame.table(post_probs)
colnames(post_df) <- c("iteration", "row", "genotype", "prob")
post_df$row <- as.integer(post_df$row)
post_df <- post_df |> left_join(newdata_fine, by = c("row" = "id"))
prob_df <- post_df |>
  group_by(Species, Sex, Strata_Age_Generations, genotype) |>
  summarise(
    mean_prob  = mean(prob),
    lower_prob = quantile(prob, 0.025),
    upper_prob = quantile(prob, 0.975),
    .groups = "drop") |>
  mutate(analysis = paste(Species, "probabilities"), line_group = Sex)
sex_diff_df <- post_df |>
  pivot_wider(
    id_cols = c(iteration, Species, Strata_Age_Generations, genotype),
    names_from = Sex,
    values_from = prob) |>
  mutate(female_minus_male = Female - Male) |>
  group_by(Species, Strata_Age_Generations, genotype) |>
  summarise(
    mean_prob  = mean(female_minus_male),
    lower_prob = quantile(female_minus_male, 0.025),
    upper_prob = quantile(female_minus_male, 0.975),
    .groups = "drop") |>  mutate(analysis = "Sex diff within species", line_group = Species)
sex_diff_interaction_df <- post_df |>
  pivot_wider(
    id_cols = c(iteration, Species, Strata_Age_Generations, genotype),
    names_from = Sex,
    values_from = prob) |>
  mutate(female_minus_male = Female - Male) |>
  select(iteration, Species, Strata_Age_Generations, genotype, female_minus_male) |>
  pivot_wider(
    id_cols = c(iteration, Strata_Age_Generations, genotype),
    names_from = Species,
    values_from = female_minus_male) |>
  mutate(diff_of_diffs = `Raso lark` - Skylark) |>
  group_by(Strata_Age_Generations, genotype) |>
  summarise(
    mean_prob  = mean(diff_of_diffs),
    lower_prob = quantile(diff_of_diffs, 0.025),
    upper_prob = quantile(diff_of_diffs, 0.975),
    .groups = "drop") |>
  mutate(analysis = "Sex diff across species", line_group = genotype)
genotype_labels <- list('func_hom' = "Functional Homozygosity (FF)", 'LOF_hom'  = "Loss-of-function Homozygosity (LL)", 'Mask_Het' = "Masked Heterozygosity (FL)", 'Func_Het' = expression(paste("Functional Heterozygosit (", F^1, F^2, ")")))
genotype_labeller <- function(variable,value){
  return(genotype_labels[value])
}

Rasolark_genotypes <- ggplot() +
  geom_ribbon(data = prob_df |> filter(Species == "Raso lark"), aes(x = Strata_Age_Generations, y = mean_prob, ymin = lower_prob, ymax = upper_prob, group=Sex, fill = Sex), alpha = 0.5) +
  geom_line(data=prob_df |> filter(Species == "Raso lark"), aes(x=Strata_Age_Generations*1.01, y=mean_prob*1.02, group=Sex, linetype=Species), color="black", linewidth=1.5) +
  geom_line(data=prob_df |> filter(Species == "Raso lark"), aes(x=Strata_Age_Generations*0.99, y=mean_prob*0.98, group=Sex, linetype=Species), color="black", linewidth=1.5) +
  geom_line(data=prob_df |> filter(Species == "Raso lark"), aes(x=Strata_Age_Generations*1.01, y=mean_prob*0.98, group=Sex, linetype=Species), color="black", linewidth=1.5) +
  geom_line(data=prob_df |> filter(Species == "Raso lark"), aes(x=Strata_Age_Generations*0.99, y=mean_prob*1.02, group=Sex, linetype=Species), color="black", linewidth=1.5) +
  geom_line(data = prob_df |> filter(Species == "Raso lark"), aes(x = Strata_Age_Generations, y = mean_prob, group = Sex, color=Sex, linetype=Species), size = 1.5) +
  geom_vline(data = data.frame(xintercept = age_seq), aes(xintercept = xintercept), color="#d7191c", linewidth=1, linetype=2) +
  scale_linetype_manual(values = c("Raso lark"=1, "Skylark"=2)) +
  scale_color_manual(values = c("#d7191c", "#2c7bb6"), limits=c("Female", "Male")) +
  facet_wrap(~genotype, nrow = 4) +
  scale_x_continuous(sec.axis=sec_axis(LOF~ (. * coef(agefit)[1]) + (coef(agefit)[2] * (. ^2)), name=NULL, breaks=seq(0,140,20)), breaks=c(0,5,10,15,20,25,30), expand=c(0,0)) +
  scale_y_continuous(limits=c(0,1)) +
  labs(x = NULL, y = expression(atop("Probability")), color = "Sex", fill = "Sex") +
  guides(linetype = "none") +
  coord_cartesian(xlim=c(0, max(age_seq_fine))) +
  theme_bw() +
  theme(legend.position = "right",
        legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(2, 'cm'), #change legend key width
        legend.title = element_text(size=15, face="bold"), #change legend title font size
        legend.text = element_text(size=15), #change legend text font size
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.spacing = unit(1, "lines"),
        axis.title.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=15, color="black"),
        axis.text.x = element_text(size=16, color="black"))

Skylark_genotypes <- ggplot() +
  geom_ribbon(data = prob_df |> filter(Species == "Skylark"), aes(x = Strata_Age_Generations, y = mean_prob, ymin = lower_prob, ymax = upper_prob, group=Sex, fill = Sex), alpha = 0.5) +
  geom_line(data=prob_df |> filter(Species == "Skylark"), aes(x=Strata_Age_Generations*1.01, y=mean_prob*1.02, group=Sex, linetype=Species), color="black", linewidth=1.5) +
  geom_line(data=prob_df |> filter(Species == "Skylark"), aes(x=Strata_Age_Generations*0.99, y=mean_prob*0.98, group=Sex, linetype=Species), color="black", linewidth=1.5) +
  geom_line(data=prob_df |> filter(Species == "Skylark"), aes(x=Strata_Age_Generations*1.01, y=mean_prob*0.98, group=Sex, linetype=Species), color="black", linewidth=1.5) +
  geom_line(data=prob_df |> filter(Species == "Skylark"), aes(x=Strata_Age_Generations*0.99, y=mean_prob*1.02, group=Sex, linetype=Species), color="black", linewidth=1.5) +
  geom_line(data = prob_df |> filter(Species == "Skylark"), aes(x = Strata_Age_Generations, y = mean_prob, group = Sex, color=Sex, linetype=Species), size = 1.5) +
  geom_vline(data = data.frame(xintercept = age_seq), aes(xintercept = xintercept), color="#d7191c", linewidth=1, linetype=2) +
  scale_color_manual(values = c("#d7191c", "#2c7bb6"), limits=c("Female", "Male")) +
  scale_linetype_manual(values = c("Raso lark"=1, "Skylark"=2)) +
  facet_wrap(~genotype, nrow = 4, strip.position = "left") +
  scale_x_continuous(sec.axis=sec_axis(LOF~ (. * coef(agefit)[1]) + (coef(agefit)[2] * (. ^2)), name="Age (million years in larks)", breaks=seq(0,140,20)), breaks=c(0,5,10,15,20,25,30), expand=c(0,0)) +
  scale_y_continuous(limits=c(0,1)) +
  labs(x = "Age (million generations)", y = NULL, color = "Sex", fill = "Sex") +
  guides(fill="none", color="none", linetype = "none") +
  coord_cartesian(xlim=c(0, max(age_seq_fine))) +
  theme_bw() +
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        panel.spacing = unit(1, "lines"),
        axis.title.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=15, color="black"),
        axis.text.x = element_text(size=16, color="black"))

sex_diff_genotypes <- ggplot() +
  geom_ribbon(data = sex_diff_df, aes(x = Strata_Age_Generations, y = mean_prob, ymin = lower_prob, ymax = upper_prob, group=Species, fill=Species), alpha = 0.5) +
  geom_line(data=sex_diff_df, aes(x=Strata_Age_Generations*1.01, y=mean_prob*1.02, group=Species, linetype=Species), color="black", linewidth=1.5) +
  geom_line(data=sex_diff_df, aes(x=Strata_Age_Generations*0.99, y=mean_prob*0.98, group=Species, linetype=Species), color="black", linewidth=1.5) +
  geom_line(data=sex_diff_df, aes(x=Strata_Age_Generations*1.01, y=mean_prob*0.98, group=Species, linetype=Species), color="black", linewidth=1.5) +
  geom_line(data=sex_diff_df, aes(x=Strata_Age_Generations*0.99, y=mean_prob*1.02, group=Species, linetype=Species), color="black", linewidth=1.5) +
  geom_line(data = sex_diff_df, aes(x = Strata_Age_Generations, y = mean_prob, group=Species, linetype=Species, color=Species), size = 1.5) +
  geom_vline(data = data.frame(xintercept = age_seq), aes(xintercept = xintercept), color="#d7191c", linewidth=1, linetype=2) +
  geom_hline(yintercept = 0, color="black", linewidth=1, linetype=2) +
  scale_color_manual(values = c("Raso lark" = "#92c5de", "Skylark" = "#f4a582")) +
  scale_fill_manual(values = c("Raso lark" = "#92c5de", "Skylark" = "#f4a582")) +
  scale_linetype_manual(values = c("Raso lark"=1, "Skylark"=2)) +
  facet_wrap(~genotype, nrow = 4) +
  scale_x_continuous(sec.axis=sec_axis(LOF~ (. * coef(agefit)[1]) + (coef(agefit)[2] * (. ^2)), name=NULL, breaks=seq(0,140,20)), breaks=c(0,5,10,15,20,25,30), expand=c(0,0)) +
  scale_y_continuous(limits=c(-1,1)) +
  labs(x = NULL, y = NULL, color = "Species", fill = "Species") +
  coord_cartesian(xlim=c(0, max(age_seq_fine))) +
  theme_bw() +
  theme(legend.position = "right",
        legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(2, 'cm'), #change legend key width
        legend.title = element_text(size=15, face="bold"), #change legend title font size
        legend.text = element_text(size=15), #change legend text font size
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.spacing = unit(1, "lines"),
        axis.title.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=15, color="black"),
        axis.text.x = element_text(size=16, color="black"))


dist=0.4
yscales <- data.frame(genotype = c("func_hom", "LOF_hom", "Mask_Het", "Func_Het"),
                      ymin = c(median(sex_diff_interaction_df$mean_prob[which(sex_diff_interaction_df$genotype == "func_hom")])-dist,
                      median(sex_diff_interaction_df$mean_prob[which(sex_diff_interaction_df$genotype == "LOF_hom")])-dist,
                      median(sex_diff_interaction_df$mean_prob[which(sex_diff_interaction_df$genotype == "Mask_Het")])-dist,
                      median(sex_diff_interaction_df$mean_prob[which(sex_diff_interaction_df$genotype == "Func_Het")])-dist),
                      ymax = c(median(sex_diff_interaction_df$mean_prob[which(sex_diff_interaction_df$genotype == "func_hom")])+dist,
                      median(sex_diff_interaction_df$mean_prob[which(sex_diff_interaction_df$genotype == "LOF_hom")])+dist,
                      median(sex_diff_interaction_df$mean_prob[which(sex_diff_interaction_df$genotype == "Mask_Het")])+dist,
                      median(sex_diff_interaction_df$mean_prob[which(sex_diff_interaction_df$genotype == "Func_Het")])+dist))
sex_diff_interaction_df2 <- sex_diff_interaction_df |> left_join(yscales, by = "genotype")


sex_diff_species_genotypes <- ggplot() +
  geom_ribbon(data = sex_diff_interaction_df2, aes(x = Strata_Age_Generations, y = mean_prob, ymin = lower_prob, ymax = upper_prob), alpha = 0.5) +
  geom_line(data=sex_diff_interaction_df2, aes(x=Strata_Age_Generations*1.005, y=mean_prob*1.01), color="black", linewidth=1.3) +
  geom_line(data=sex_diff_interaction_df2, aes(x=Strata_Age_Generations*0.995, y=mean_prob*0.99), color="black", linewidth=1.3) +
  geom_line(data=sex_diff_interaction_df2, aes(x=Strata_Age_Generations*1.005, y=mean_prob*0.99), color="black", linewidth=1.3) +
  geom_line(data=sex_diff_interaction_df2, aes(x=Strata_Age_Generations*0.995, y=mean_prob*1.01), color="black", linewidth=1.3) +
  geom_line(data = sex_diff_interaction_df2, aes(x = Strata_Age_Generations, y = mean_prob), size = 1.5) +
  geom_vline(data = data.frame(xintercept = age_seq), aes(xintercept = xintercept), color="#d7191c", linewidth=1, linetype=2) +
  geom_hline(yintercept = 0, color="black", linewidth=1, linetype=2) +
  facet_wrap(~genotype, nrow = 4, strip.position = "right", scale="free_y") +
  scale_x_continuous(sec.axis=sec_axis(LOF~ (. * coef(agefit)[1]) + (coef(agefit)[2] * (. ^2)), name=NULL, breaks=seq(0,140,20)), breaks=c(0,5,10,15,20,25,30), expand=c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  geom_blank(data = sex_diff_interaction_df2, aes(y = ymin)) +
  geom_blank(data = sex_diff_interaction_df2, aes(y = ymax)) + 
  labs(x = NULL, y = "Raso lark (Female-Male) - Skylark (Female-Male)", color = "Sex", fill = "Sex") +
  guides(fill="none", color="none", linetype = "none") +
  coord_cartesian(xlim=c(0, max(age_seq_fine))) +
  theme_bw() +
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        panel.spacing = unit(1, "lines"),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        axis.text.y = element_text(size=15, color="black"),
        axis.text.x = element_text(size=15, color="black"))

### Assemble Figure 4 ##############################################################################

fig4 <- (Rasolark_genotypes | Skylark_genotypes | sex_diff_genotypes) +
  plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 20, face="bold"), plot.tag.position = c(0.05, 1.00), legend.position = "right")

png("Rasolark_2021/Results/Figures/Figure4.png", width=6000, height=4500, res=300)
fig4
dev.off()

### More plots ##############################################################

# plot logodds
labels2 <- as.data.frame(odds_ratios_df2$Variable[which(odds_ratios_df2$data == "Heterozygous masking (FL)")])
colnames(labels2) <- "term"
labels2$contVar <- rev(1:nrow(labels2))
odds_ratios_df2$contVar <- labels2$contVar
odds_ratios_df2 <- odds_ratios_df2[which(nchar(odds_ratios_df2$term) > 0 & odds_ratios_df2$term != "Intercept"),]

plot_OR2 <- ggplot(odds_ratios_df2, aes(x = contVar, y = logOdds_Ratio)) +
  geom_point(size = 3, color = "Black") +  # Plot odds ratios as points
  geom_errorbar(aes(ymin = logLower_CI, ymax = logUpper_CI), width = 0.2, color = "black") +  # Add CI bars
  geom_hline(yintercept = log10(1), linewidth=1, linetype = 2, color = "#d7191c") +  # Reference line at OR = 0
  scale_x_continuous(labels=labels2$term, breaks=labels2$contVar, sec.axis = dup_axis(name=NULL)) +
  facet_wrap(~data) +
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


# Plot pHaplo x GERP
pHaplo_vals <- seq(0, 1.0, 0.2)
GERP_vals <-  seq(min(data_gt2long$GERP_max), max(data_gt2long$GERP_max), length.out = 100)
newdata <- expand.grid(
  Species = unique(data_gt2long$Species),
  Sex = unique(data_gt2long$Sex),
  Strata_Age_Generations = median(data_gt2long$Strata_Age_Generations),
  GERP_max = as.numeric(GERP_vals),
  pHaplo = pHaplo_vals,
  logGeneLen = median(data_gt2long$logGeneLen))
newdata$id <- seq_len(nrow(newdata))
post_probs <- posterior_epred(Model2, newdata = newdata)
post_df <- as.data.frame.table(post_probs)
colnames(post_df) <- c("iteration", "row", "genotype", "prob")
post_df$row <- as.integer(post_df$row)
post_df <- post_df |> left_join(newdata, by = c("row" = "id"))
prob_summary_nosex <- post_df |>
  group_by(iteration, Species, genotype, GERP_max, pHaplo) |>
  summarise(prob = mean(prob), .groups = "drop") |>
  group_by(Species, genotype, GERP_max, pHaplo) |>
  summarise(
    mean_prob = mean(prob),
    lower_prob = quantile(prob, 0.025),
    upper_prob = quantile(prob, 0.975),
    .groups = "drop")
prob_func_hom <- prob_summary_nosex |> filter(genotype == "func_hom")

rect_df <- data.frame(
  xmin = min(data_gt2long$GERP_max),
  xmax = quantile(data_gt2long$GERP_max, probs = 0.05),
  ymin = 0,
  ymax = 0.75,
  Species = unique(prob_func_hom$Species))


funchom_pHaplo_GERP <- ggplot() +
  geom_ribbon(data = prob_func_hom, aes(x = GERP_max, y = mean_prob, ymin = lower_prob, ymax = upper_prob, group = pHaplo, fill=pHaplo), alpha = 0.5) +
  geom_line(data = prob_func_hom, aes(x = GERP_max, y = mean_prob, group = pHaplo, color=pHaplo), size = 1.5) +
  geom_rect(data = rect_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), inherit.aes = FALSE, alpha = 0.5) +
  geom_vline(data = data.frame(xintercept = median(data_gt2long$GERP_max)), aes(xintercept = xintercept), color="#d7191c", linewidth=1, linetype=2) +
  facet_wrap(~ Species) +
  scale_color_viridis_c(option = "plasma") +
  scale_fill_viridis_c(option = "plasma") +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(limits=c(0,0.75), expand=c(0,0), breaks=c(0.00, 0.25, 0.50, 0.75)) +
  labs(x = "Maximum GERP score", y = expression(atop("Probability", "Functional homozygosity")), color = "pHaplo") +
  theme_bw() +
  theme(legend.position = "right",
        legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(2, 'cm'), #change legend key width
        legend.title = element_text(size=15, face="bold"), #change legend title font size
        legend.text = element_text(size=15), #change legend text font size
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(size=14, face="bold"),
        panel.spacing = unit(3, "lines"),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        axis.text.y = element_text(size=15, color="black"),
        axis.text.x = element_text(size=15, color="black"))


png("Rasolark_2021/Results/Figures/Supplementary/funchom_pHaplo_GERP.png", width=3000, height=2000, res=300)
funchom_pHaplo_GERP
dev.off()


# Plot length as function of species
length_vals <-  seq(min(data_gt2long$logGeneLen), max(data_gt2long$logGeneLen), length.out = 100)
newdata <- expand.grid(
  Species = unique(data_gt2long$Species),
  Sex = unique(data_gt2long$Sex),
  Strata_Age_Generations = median(data_gt2long$Strata_Age_Generations),
  GERP_max = median(data_gt2long$GERP_max),
  pHaplo = median(data_gt2long$pHaplo),
  logGeneLen = as.numeric(length_vals))
newdata$id <- seq_len(nrow(newdata))
post_probs <- posterior_epred(Model2, newdata = newdata)
post_df <- as.data.frame.table(post_probs)
colnames(post_df) <- c("iteration", "row", "genotype", "prob")
post_df$row <- as.integer(post_df$row)
post_df <- post_df |> left_join(newdata, by = c("row" = "id"))
prob_summary_nosex <- post_df |>
  group_by(iteration, Species, genotype, logGeneLen) |>
  summarise(prob = mean(prob), .groups = "drop") |>
  group_by(Species, genotype, logGeneLen) |>
  summarise(
    mean_prob = mean(prob),
    lower_prob = quantile(prob, 0.025),
    upper_prob = quantile(prob, 0.975),
    .groups = "drop")
prob_func_hom_length <- prob_summary_nosex |> filter(genotype == "func_hom")

funchom_length <- ggplot() +
  geom_ribbon(data = prob_func_hom_length, aes(x = logGeneLen, y = mean_prob, ymin = lower_prob, ymax = upper_prob, group=Species, fill=Species), alpha = 0.5) +
  geom_line(data=prob_func_hom_length, aes(x=logGeneLen*1.002, y=mean_prob*1.002, group=Species, linetype=Species), color="black", linewidth=1.5) +
  geom_line(data=prob_func_hom_length, aes(x=logGeneLen*0.998, y=mean_prob*0.998, group=Species, linetype=Species), color="black", linewidth=1.5) +
  geom_line(data=prob_func_hom_length, aes(x=logGeneLen*1.002, y=mean_prob*0.998, group=Species, linetype=Species), color="black", linewidth=1.5) +
  geom_line(data=prob_func_hom_length, aes(x=logGeneLen*0.998, y=mean_prob*1.002, group=Species, linetype=Species), color="black", linewidth=1.5) +
  geom_line(data = prob_func_hom_length, aes(x = logGeneLen, y = mean_prob, group=Species, color=Species, linetype=Species), size = 1.5) +
  scale_color_manual(values = c("Raso lark" = "#92c5de", "Skylark" = "#f4a582")) +
  scale_fill_manual(values = c("Raso lark" = "#92c5de", "Skylark" = "#f4a582")) +
  scale_linetype_manual(values = c("Raso lark"=1, "Skylark"=2)) +
  scale_x_continuous(expand=c(0,0)) +
 scale_y_continuous(limits=c(0,1), expand=c(0,0), breaks=c(0.00, 0.25, 0.50, 0.75)) +
  labs(x = expression("log"[10]*"(Gene length)"), y = expression(atop("Probability", "Functional homozygosity")), color = "Species") +
  theme_bw() +
  theme(legend.position = "right",
        legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(2, 'cm'), #change legend key width
        legend.title = element_text(size=15, face="bold"), #change legend title font size
        legend.text = element_text(size=15), #change legend text font size
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(size=14, face="bold"),
        panel.spacing = unit(3, "lines"),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        axis.text.y = element_text(size=15, color="black"),
        axis.text.x = element_text(size=15, color="black"))

png("Rasolark_2021/Results/Figures/Supplementary/funchom_length.png", width=3000, height=2000, res=300)
funchom_length
dev.off()


### Statistically test predictors for proportion of LOF Z gametolog between sexes ###########################
data_gt3 <- data_gt[which(data_gt$Strata != "Autosomal" & data_gt$Strata != "PAR"),]

# Get proportion of W-LOF
WLOF <- data_gt3 |>
  filter(Sex == "Female") |>
  mutate(W_loss = (LOF2 == 1) | (W_Dropout == 1)) |>
  group_by(Species, geneID) |>
  summarise(
    W_LOF = mean(W_loss, na.rm = TRUE),
    n_females = sum(!is.na(LOF2) | !is.na(W_Dropout)),
    n_WLOF    = sum(W_loss, na.rm = TRUE),
    .groups = "drop")
data_gt3long <- data_gt3 |>
  mutate(LOF2 = ifelse(Sex == "Female", NA, LOF2)) |>
  pivot_longer(
    cols = c(LOF1, LOF2),
    names_to = "Allele_copy",
    values_to = "Z_LOF",
    values_drop_na = TRUE)
data_gt3long <- data_gt3long |>
  left_join(
    WLOF,
    by = c("Species", "geneID"))
table(data_gt3long$W_LOF, data_gt3long$Species)
data_gt3long$W_LOF[which(data_gt3long$W_LOF == 0.1)] <- 0
data_gt3long$W_LOF <- factor(data_gt3long$W_LOF, labels=c("Functional", "Loss-of-function"), levels=c(0, 1))


# Check for collinearity
vif_model <- glmmTMB(Z_LOF ~ Species + Sex + Strata_Age_Generations + W_LOF + GERP_max + pHaplo + logGeneLen,
                     family=binomial, data = data_gt3long, na.action = "na.fail", REML=F)
check_collinearity(vif_model)

m0 <- glmmTMB(Z_LOF ~
                Species + Sex + Strata_Age_Generations + W_LOF + GERP_max + pHaplo + logGeneLen +
                Species:Sex + Species:Strata_Age_Generations + Species:W_LOF + Species:GERP_max + Species:pHaplo + Species:logGeneLen +
                Sex:Strata_Age_Generations + Sex:W_LOF + Sex:GERP_max + Sex:pHaplo + Sex:logGeneLen +
                Strata_Age_Generations:GERP_max + Strata_Age_Generations:W_LOF + Strata_Age_Generations:pHaplo + Strata_Age_Generations:logGeneLen +
                W_LOF:GERP_max + W_LOF:pHaplo + W_LOF:logGeneLen +
                GERP_max:pHaplo + GERP_max:logGeneLen +
                pHaplo:logGeneLen,
              family=binomial, data = data_gt3long, na.action = "na.fail", REML=F)
simulateResiduals(fittedModel = m0, plot = T)
summary(m0)

m1 <- glmmTMB(Z_LOF ~
                Species + Sex + Strata_Age_Generations + W_LOF + GERP_max + pHaplo + logGeneLen +
                Species:Strata_Age_Generations + Species:W_LOF + Species:GERP_max + Species:pHaplo + Species:logGeneLen +
                Sex:Strata_Age_Generations + Sex:W_LOF + Sex:GERP_max + Sex:pHaplo + Sex:logGeneLen +
                Strata_Age_Generations:GERP_max + Strata_Age_Generations:W_LOF + Strata_Age_Generations:pHaplo + Strata_Age_Generations:logGeneLen +
                W_LOF:GERP_max + W_LOF:pHaplo + W_LOF:logGeneLen +
                GERP_max:pHaplo + GERP_max:logGeneLen +
                pHaplo:logGeneLen,
              family=binomial, data = data_gt3long, na.action = "na.fail", REML=F)
simulateResiduals(fittedModel = m1, plot = T)
anova(m0, m1)
summary(m1)

m2 <- glmmTMB(Z_LOF ~
                Species + Sex + Strata_Age_Generations + W_LOF + GERP_max + pHaplo + logGeneLen +
                Species:Strata_Age_Generations + Species:W_LOF + Species:GERP_max + Species:pHaplo + Species:logGeneLen +
                Sex:W_LOF + Sex:GERP_max + Sex:pHaplo + Sex:logGeneLen +
                Strata_Age_Generations:GERP_max + Strata_Age_Generations:W_LOF + Strata_Age_Generations:pHaplo + Strata_Age_Generations:logGeneLen +
                W_LOF:GERP_max + W_LOF:pHaplo + W_LOF:logGeneLen +
                GERP_max:pHaplo + GERP_max:logGeneLen +
                pHaplo:logGeneLen,
              family=binomial, data = data_gt3long, na.action = "na.fail", REML=F)
simulateResiduals(fittedModel = m2, plot = T)
anova(m1, m2)
summary(m2)

m3 <- glmmTMB(Z_LOF ~
                Species + Sex + Strata_Age_Generations + W_LOF + GERP_max + pHaplo + logGeneLen +
                Species:Strata_Age_Generations + Species:W_LOF + Species:GERP_max + Species:pHaplo + Species:logGeneLen +
                Sex:W_LOF + Sex:GERP_max + Sex:logGeneLen +
                Strata_Age_Generations:GERP_max + Strata_Age_Generations:W_LOF + Strata_Age_Generations:pHaplo + Strata_Age_Generations:logGeneLen +
                W_LOF:GERP_max + W_LOF:pHaplo + W_LOF:logGeneLen +
                GERP_max:pHaplo + GERP_max:logGeneLen +
                pHaplo:logGeneLen,
              family=binomial, data = data_gt3long, na.action = "na.fail", REML=F)
simulateResiduals(fittedModel = m3, plot = T)
anova(m2, m3)
summary(m3)

m4 <- glmmTMB(Z_LOF ~
                Species + Sex + Strata_Age_Generations + W_LOF + GERP_max + pHaplo + logGeneLen +
                Species:Strata_Age_Generations + Species:W_LOF + Species:GERP_max + Species:pHaplo + Species:logGeneLen +
                Sex:GERP_max + Sex:logGeneLen +
                Strata_Age_Generations:GERP_max + Strata_Age_Generations:W_LOF + Strata_Age_Generations:pHaplo + Strata_Age_Generations:logGeneLen +
                W_LOF:GERP_max + W_LOF:pHaplo + W_LOF:logGeneLen +
                GERP_max:pHaplo + GERP_max:logGeneLen +
                pHaplo:logGeneLen,
              family=binomial, data = data_gt3long, na.action = "na.fail", REML=F)
simulateResiduals(fittedModel = m4, plot = T)
anova(m3, m4)
summary(m4)

m5 <- glmmTMB(Z_LOF ~
                Species + Sex + Strata_Age_Generations + W_LOF + GERP_max + pHaplo + logGeneLen +
                Species:Strata_Age_Generations + Species:W_LOF + Species:GERP_max + Species:pHaplo + Species:logGeneLen +
                Sex:GERP_max +
                Strata_Age_Generations:GERP_max + Strata_Age_Generations:W_LOF + Strata_Age_Generations:pHaplo + Strata_Age_Generations:logGeneLen +
                W_LOF:GERP_max + W_LOF:pHaplo + W_LOF:logGeneLen +
                GERP_max:pHaplo + GERP_max:logGeneLen +
                pHaplo:logGeneLen,
              family=binomial, data = data_gt3long, na.action = "na.fail", REML=F)
simulateResiduals(fittedModel = m5, plot = T)
anova(m4, m5)
summary(m5)

m6 <- glmmTMB(Z_LOF ~
                Species + Sex + Strata_Age_Generations + W_LOF + GERP_max + pHaplo + logGeneLen +
                Species:Strata_Age_Generations + Species:W_LOF + Species:GERP_max + Species:logGeneLen +
                Strata_Age_Generations:GERP_max + Strata_Age_Generations:W_LOF + Strata_Age_Generations:pHaplo + Strata_Age_Generations:logGeneLen +
                W_LOF:GERP_max + W_LOF:pHaplo + W_LOF:logGeneLen +
                GERP_max:pHaplo + GERP_max:logGeneLen +
                pHaplo:logGeneLen,
              family=binomial, data = data_gt3long, na.action = "na.fail", REML=F)
simulateResiduals(fittedModel = m6, plot = T)
anova(m5, m6)
summary(m6)

m7 <- glmmTMB(Z_LOF ~
                Species + Strata_Age_Generations + W_LOF + GERP_max + pHaplo + logGeneLen +
                Species:Strata_Age_Generations + Species:W_LOF + Species:GERP_max + Species:logGeneLen +
                Strata_Age_Generations:GERP_max + Strata_Age_Generations:W_LOF + Strata_Age_Generations:pHaplo + Strata_Age_Generations:logGeneLen +
                W_LOF:GERP_max + W_LOF:pHaplo + W_LOF:logGeneLen +
                GERP_max:pHaplo + GERP_max:logGeneLen +
                pHaplo:logGeneLen,
              family=binomial, data = data_gt3long, na.action = "na.fail", REML=F)
simulateResiduals(fittedModel = m7, plot = T)
anova(m6, m7)
summary(m7)

m8 <- glmmTMB(Z_LOF ~
                Species + Strata_Age_Generations + W_LOF + GERP_max + pHaplo + logGeneLen +
                Species:Strata_Age_Generations + Species:W_LOF + Species:GERP_max + Species:logGeneLen +
                Strata_Age_Generations:GERP_max + Strata_Age_Generations:W_LOF + Strata_Age_Generations:pHaplo + Strata_Age_Generations:logGeneLen +
                W_LOF:GERP_max + W_LOF:pHaplo + W_LOF:logGeneLen +
                GERP_max:pHaplo + GERP_max:logGeneLen,
              family=binomial, data = data_gt3long, na.action = "na.fail", REML=F)
simulateResiduals(fittedModel = m8, plot = T)
anova(m7, m8)
summary(m8)

m9 <- glmmTMB(Z_LOF ~
                Species + Strata_Age_Generations + W_LOF + GERP_max + pHaplo + logGeneLen +
                Species:Strata_Age_Generations + Species:W_LOF + Species:GERP_max +
                Strata_Age_Generations:GERP_max + Strata_Age_Generations:W_LOF + Strata_Age_Generations:pHaplo + Strata_Age_Generations:logGeneLen +
                W_LOF:GERP_max + W_LOF:pHaplo + W_LOF:logGeneLen +
                GERP_max:pHaplo + GERP_max:logGeneLen,
              family=binomial, data = data_gt3long, na.action = "na.fail", REML=F)
simulateResiduals(fittedModel = m9, plot = T)
anova(m8, m9)
summary(m9)

m10 <- glmmTMB(Z_LOF ~
                 Species + Strata_Age_Generations + W_LOF + GERP_max + pHaplo + logGeneLen +
                 Species:W_LOF + Species:GERP_max +
                 Strata_Age_Generations:GERP_max + Strata_Age_Generations:W_LOF + Strata_Age_Generations:pHaplo + Strata_Age_Generations:logGeneLen +
                 W_LOF:GERP_max + W_LOF:pHaplo + W_LOF:logGeneLen +
                 GERP_max:pHaplo + GERP_max:logGeneLen,
               family=binomial, data = data_gt3long, na.action = "na.fail", REML=F)
simulateResiduals(fittedModel = m10, plot = T)
anova(m9, m10)
summary(m10)

Model3 <- m10
simulateResiduals(fittedModel = Model3, plot = T, n=1000)
summary(Model3)
r2(Model3)

# Get odds_ratios
coef3 <- summary(Model3)$coefficients$cond
odds_ratios3 <- exp(coef3[, "Estimate"])
odds_ratios_df3 <- data.frame(
  Variable = rownames(coef3),
  Odds_Ratio = odds_ratios3,
  Lower_CI = exp(coef3[, "Estimate"] - 1.96 * coef3[, "Std. Error"]),
  Upper_CI = exp(coef3[, "Estimate"] + 1.96 * coef3[, "Std. Error"])
)
print(odds_ratios_df3)


odds_ratios_df3$Variable <- c("Intercept ***", "Species (Skylark) ***", "Strata age ***", "W frequency (loss-of-function) ***", "Maximum GERP score ***", "pHaplo ***", "log10(Gene length)",
                              "Species (Skylark) × W frequency (loss-of-function) ***", "Species (Skylark) × Maximum GERP score *", 
                              "Strata age × Maximum GERP score ***", "Strata age × W frequency (loss-of-function) ***", "Strata age × pHaplo **", "Strata age × log10(Gene length) ***",
                              "W frequency (loss-of-function) × Maximum GERP score **", "W frequency (loss-of-function) × pHaplo ***", "W frequency (loss-of-function) × log10(Gene length) ***",
                              "Maximum GERP score × pHaplo ***",  "Maximum GERP score × log10(Gene length) *")

odds_ratios_df3$Variable <- factor(odds_ratios_df3$Variable, order=T, levels=c("Intercept ***", "Species (Skylark) ***", "Strata age ***", "W frequency (loss-of-function) ***", "Maximum GERP score ***", "pHaplo ***", "log10(Gene length)",
                              "Species (Skylark) × W frequency (loss-of-function) ***", "Species (Skylark) × Maximum GERP score *", 
                              "Strata age × W frequency (loss-of-function) ***", "Strata age × Maximum GERP score ***", "Strata age × pHaplo **", "Strata age × log10(Gene length) ***",
                              "W frequency (loss-of-function) × Maximum GERP score **", "W frequency (loss-of-function) × pHaplo ***", "W frequency (loss-of-function) × log10(Gene length) ***",
                              "Maximum GERP score × pHaplo ***",  "Maximum GERP score × log10(Gene length) *"))

odds_ratios_df3$logOdds_Ratio <- log10(odds_ratios_df3$Odds_Ratio)
odds_ratios_df3$logLower_CI <- log10(odds_ratios_df3$Lower_CI)
odds_ratios_df3$logUpper_CI <- log10(odds_ratios_df3$Upper_CI)
odds_ratios_df3 <- odds_ratios_df3[order(odds_ratios_df3$Variable),]
odds_ratios_df3 <- odds_ratios_df3[which(odds_ratios_df3$Variable != "Intercept ***"),]
odds_ratios_df3$contVar <- rev(1:nrow(odds_ratios_df3))

plot_OR3 <- ggplot(odds_ratios_df3, aes(x = contVar, y = logOdds_Ratio)) +
  geom_point(size = 3, color = "Black") +  # Plot odds ratios as points
  geom_errorbar(aes(ymin = logLower_CI, ymax = logUpper_CI), width = 0.2, color = "black") +  # Add CI bars
  geom_hline(yintercept = log10(1), linewidth=1, linetype = 2, color = "#d7191c") +  # Reference line at OR = 0
  scale_x_continuous(labels=odds_ratios_df3$Variable, breaks=odds_ratios_df3$contVar, sec.axis = dup_axis(name=NULL)) +
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


# Plot Z loss of function for W functional and W LOF
age_seq_fine <- seq(0, max(data_gt3long$Strata_Age_Generations), length.out = 100)
newdata3 <- expand.grid(
  Species = unique(data_gt3long$Species),
  W_LOF = unique(data_gt3long$W_LOF),
  Strata_Age_Generations = age_seq_fine,
  GERP_max = median(data_gt3long$GERP_max, na.rm = TRUE),
  pHaplo = median(data_gt3long$pHaplo, na.rm = TRUE),
  logGeneLen = median(data_gt3long$logGeneLen, na.rm = TRUE))
newdata3$id <- 1:nrow(newdata3)
pred <- predict(Model3, newdata = newdata3, type = "response", se.fit = TRUE)
newdata3$mean_prob  <- pred$fit
newdata3$lower_prob <- plogis(qlogis(pred$fit) - 1.96 * pred$se.fit)
newdata3$upper_prob <- plogis(qlogis(pred$fit) + 1.96 * pred$se.fit)
newdata3$W_LOF <- factor(newdata3$W_LOF, labels=c("W Functional", "W Loss-of-function"), levels=c("Functional", "Loss-of-function"))

# Get proportions
gt_prop3 <- data_gt3long |>
  group_by(Individual, Strata_Age_Generations, Species, W_LOF) |>
  summarise(n_total = n(),
    n_LOF   = sum(Z_LOF == 1),
    prop    = n_LOF / n_total,
    .groups = "drop") |>
  ungroup()

set.seed(123)
gt_prop3jit <- gt_prop3 |>
  mutate(x_jit = Strata_Age_Generations + runif(n(), -0.2, 0.2))
gt_prop3jit$W_LOF <- factor(gt_prop3jit$W_LOF, labels=c("W Functional", "W Loss-of-function"), levels=c("Functional", "Loss-of-function"))


ZLOF_species_WLOF <- ggplot() +
  geom_ribbon(data = newdata3, aes(x = Strata_Age_Generations, y = mean_prob, ymin = lower_prob, ymax = upper_prob, group=interaction(Species, W_LOF), fill=Species), alpha = 0.5) +
  geom_line(data=newdata3, aes(x=Strata_Age_Generations*1.01, y=mean_prob*1.01, group=interaction(Species, W_LOF), linetype=Species), color="black", linewidth=1.5) +
  geom_line(data=newdata3, aes(x=Strata_Age_Generations*0.99, y=mean_prob*0.99, group=interaction(Species, W_LOF), linetype=Species), color="black", linewidth=1.5) +
  geom_line(data=newdata3, aes(x=Strata_Age_Generations*1.01, y=mean_prob*0.99, group=interaction(Species, W_LOF), linetype=Species), color="black", linewidth=1.5) +
  geom_line(data=newdata3, aes(x=Strata_Age_Generations*0.99, y=mean_prob*1.01, group=interaction(Species, W_LOF), linetype=Species), color="black", linewidth=1.5) +
  geom_line(data = newdata3, aes(x = Strata_Age_Generations, y = mean_prob, group=interaction(Species, W_LOF), linetype=Species, color=Species), size = 1.5) +
  geom_point(data=gt_prop3jit, aes(x=x_jit, y=prop, group=interaction(Species, W_LOF)), color="black", size=4, alpha=0.25)+
  geom_point(data=gt_prop3jit, aes(x=x_jit, y=prop, group=interaction(Species, W_LOF), color=Species), size=3, alpha=0.5) +
  scale_color_manual(values = c("Raso lark" = "#92c5de", "Skylark" = "#f4a582")) +
  scale_fill_manual(values = c("Raso lark" = "#92c5de", "Skylark" = "#f4a582")) +
  scale_linetype_manual(values = c("Raso lark"=1, "Skylark"=2)) +
  facet_wrap(~ W_LOF, nrow=2) +
  scale_x_continuous(sec.axis=sec_axis(LOF~ (. * coef(agefit)[1]) + (coef(agefit)[2] * (. ^2)), name="Age (million years in larks)", breaks=seq(0,140,20)), breaks=c(0,5,10,15,20,25,30,35), expand=c(0,0), limits=c(0, 35)) +
  #scale_y_continuous(limits=c(0,0.3), expand=c(0,0), breaks=c(0, 0.1, 0.2, 0.3)) +
  labs(x = "Age (million generations)", y = expression(atop("Probability", "Z loss-of-function")), color="Species") +
  theme_bw() +
  theme(legend.position = "right",
        legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(2, 'cm'), #change legend key width
        legend.title = element_text(size=15, face="bold"), #change legend title font size
        legend.text = element_text(size=15), #change legend text font size
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(size=14, face="bold"),
        panel.spacing = unit(3, "lines"),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        axis.text.y = element_text(size=15, color="black"),
        axis.text.x = element_text(size=15, color="black"))


png("Rasolark_2021/Results/Figures/Supplementary/ZLOF_species_WLOF.png", width=3000, height=2500, res=300)
ZLOF_species_WLOF
dev.off()


# Plot pHaplo x GERP
pHaplo_vals <- seq(0, 1.0, 0.2)
GERP_vals <-  seq(min(data_gt2long$GERP_max), max(data_gt2long$GERP_max), length.out = 100)
newdata <- expand.grid(
  Species = "Raso lark",
  W_LOF = unique(data_gt3long$W_LOF),
  Strata_Age_Generations = median(data_gt3long$Strata_Age_Generations),
  GERP_max = as.numeric(GERP_vals),
  pHaplo = pHaplo_vals,
  logGeneLen = median(data_gt3long$logGeneLen))
newdata$id <- seq_len(nrow(newdata))
pred <- predict(Model3, newdata = newdata,type = "response", se.fit = TRUE)
newdata$mean_prob  <- pred$fit
newdata$lower_prob <- plogis(qlogis(pred$fit) - 1.96 * pred$se.fit)
newdata$upper_prob <- plogis(qlogis(pred$fit) + 1.96 * pred$se.fit)
newdata$W_LOF <- factor(newdata$W_LOF, labels=c("W Functional", "W Loss-of-function"), levels=c("Functional", "Loss-of-function"))

rect_df <- data.frame(
  xmin = min(data_gt2long$GERP_max),
  xmax = quantile(data_gt3long$GERP_max, probs = 0.05),
  ymin = 0,
  ymax = 1,
  Species = unique(prob_func_hom$Species))


ZLOF_funchom_pHaplo_GERP <- ggplot() +
  geom_ribbon(data = newdata, aes(x = GERP_max, y = mean_prob, ymin = lower_prob, ymax = upper_prob, group = pHaplo, fill=pHaplo), alpha = 0.5) +
  geom_line(data = newdata, aes(x = GERP_max, y = mean_prob, group = pHaplo, color=pHaplo), size = 1.5) +
  geom_vline(data = data.frame(xintercept = median(data_gt3$GERP_max)), aes(xintercept = xintercept), color="#d7191c", linewidth=1, linetype=2) +
  geom_rect(data = rect_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), inherit.aes = FALSE, alpha = 0.5) +
  facet_wrap(~W_LOF) +
  scale_color_viridis_c(option = "plasma") +
  scale_fill_viridis_c(option = "plasma") +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(limits=c(0,1), expand=c(0,0)) +
  labs(x = "Maximum GERP score", y = expression(atop("Probability", "Z loss-of-function")), color = "pHaplo") +
  theme_bw() +
  theme(legend.position = "right",
        legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(2, 'cm'), #change legend key width
        legend.title = element_text(size=15, face="bold"), #change legend title font size
        legend.text = element_text(size=15), #change legend text font size
        strip.background = element_blank(),
        strip.text = element_text(size=14, face="bold"),
        panel.spacing = unit(3, "lines"),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        axis.text.y = element_text(size=15, color="black"),
        axis.text.x = element_text(size=15, color="black"))


png("Rasolark_2021/Results/Figures/Supplementary/ZLOF_pHaplo_GERP.png", width=3000, height=2000, res=300)
ZLOF_funchom_pHaplo_GERP
dev.off()



# Statistically test predictors for guarding gametolog #####################################################
data_gt4 <- data_gt[which(data_gt$Strata != "Autosomal" & data_gt$Strata != "PAR" & data_gt$Sex == "Female" & !is.na(data_gt$Sheltering_gametolog_LOF)),]
data_gt4 <- data_gt4[,c("geneID", "Species", "Individual", "Strata", "Strata_Age_Generations", "logGeneLen", "pHaplo", "GERP_max", "Sheltering_gametolog_LOF")]

# Set reference categories
data_gt4$bin <- rep(0, nrow(data_gt4))
data_gt4$bin[which(data_gt4$Sheltering_gametolog_LOF == "W")] <- 1
table(data_gt4$Species, data_gt4$Strata)

# Check for collinearity
vif_model <- glmmTMB(bin ~ Species + Strata_Age_Generations + GERP_max + pHaplo + logGeneLen,
                     family=binomial, data = data_gt4, na.action = "na.fail", REML=F)
check_collinearity(vif_model)

m0 <- glmmTMB(bin ~
                Species + Strata_Age_Generations + GERP_max + pHaplo + logGeneLen +
                Species:Strata_Age_Generations + Species:GERP_max + Species:pHaplo + Species:logGeneLen +
                Strata_Age_Generations:GERP_max + Strata_Age_Generations:pHaplo + Strata_Age_Generations:logGeneLen +
                GERP_max:pHaplo + GERP_max:logGeneLen +
                pHaplo:logGeneLen,
              family=binomial, data = data_gt4, na.action = "na.fail", REML=F)
simulateResiduals(fittedModel = m0, plot = T)
summary(m0)

m1 <- glmmTMB(bin ~
                Species + Strata_Age_Generations + GERP_max + pHaplo + logGeneLen +
                Species:GERP_max + Species:pHaplo + Species:logGeneLen +
                Strata_Age_Generations:GERP_max + Strata_Age_Generations:pHaplo + Strata_Age_Generations:logGeneLen +
                GERP_max:pHaplo + GERP_max:logGeneLen +
                pHaplo:logGeneLen,
              family=binomial, data = data_gt4, na.action = "na.fail", REML=F)
simulateResiduals(fittedModel = m1, plot = T)
anova(m0, m1)
summary(m1)

m2 <- glmmTMB(bin ~
                Species + Strata_Age_Generations + GERP_max + pHaplo + logGeneLen +
                Species:GERP_max + Species:logGeneLen +
                Strata_Age_Generations:GERP_max + Strata_Age_Generations:pHaplo + Strata_Age_Generations:logGeneLen +
                GERP_max:pHaplo + GERP_max:logGeneLen +
                pHaplo:logGeneLen,
              family=binomial, data = data_gt4, na.action = "na.fail", REML=F)
simulateResiduals(fittedModel = m2, plot = T)
anova(m1, m2)
summary(m2)

Model4 <- m2
summary(Model4)
r2(Model4)

# Get odds_ratios
coef4 <- summary(Model4)$coefficients$cond
odds_ratios4 <- exp(coef4[, "Estimate"])
odds_ratios_df4 <- data.frame(
  Variable = rownames(coef4),
  Odds_Ratio = odds_ratios4,
  Lower_CI = exp(coef4[, "Estimate"] - 1.96 * coef4[, "Std. Error"]),
  Upper_CI = exp(coef4[, "Estimate"] + 1.96 * coef4[, "Std. Error"])
)
print(odds_ratios_df4)


odds_ratios_df4$Variable <- c("Intercept ***", "Species (Skylark) **", "Strata age .", "Maximum GERP score ***", "pHaplo **", "log10(Gene length) ***",
                              "Species (Skylark) × Maximum GERP score ***", "Species (Skylark) × log10(Gene length) ***", 
                              "Strata age × Maximum GERP score ***", "Strata age × pHaplo *", "Strata age × log10(Gene length) ***",
                              "Maximum GERP score × pHaplo ***", "Maximum GERP score × log10(Gene length) ***", "pHaplo × log10(Gene length) ***")


odds_ratios_df4$logOdds_Ratio <- log10(odds_ratios_df4$Odds_Ratio)
odds_ratios_df4$logLower_CI <- log10(odds_ratios_df4$Lower_CI)
odds_ratios_df4$logUpper_CI <- log10(odds_ratios_df4$Upper_CI)
odds_ratios_df4 <- odds_ratios_df4[which(odds_ratios_df4$Variable != "Intercept ***"),]
odds_ratios_df4$contVar <- rev(1:nrow(odds_ratios_df4))

plot_OR4 <- ggplot(odds_ratios_df4, aes(x = contVar, y = logOdds_Ratio)) +
  geom_point(size = 3, color = "Black") +  # Plot odds ratios as points
  geom_errorbar(aes(ymin = logLower_CI, ymax = logUpper_CI), width = 0.2, color = "black") +  # Add CI bars
  geom_hline(yintercept = log10(1), linewidth=1, linetype = 2, color = "#d7191c") +  # Reference line at OR = 0
  scale_x_continuous(labels=odds_ratios_df4$Variable, breaks=odds_ratios_df4$contVar, sec.axis = dup_axis(name=NULL)) +
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


### Plot frequency of genes masking LOF W
shetlering_frequency <- data_gt4 |>
  group_by(
    geneID,
    Species,
    Strata,
    Sheltering_gametolog_LOF) |>
  summarise(n = n(), .groups = "drop")
shetlering_frequency <- shetlering_frequency |>
  group_by(
    Species,
    Strata,
    Sheltering_gametolog_LOF,
    n) |>
  summarise(
    n_genes = n_distinct(geneID),
    .groups = "drop")
shetlering_frequency$Strata <- factor(shetlering_frequency$Strata, order=T, labels=rev(c("S0: 32.3 MG", "S1: 25.0 MG", "S2: 23.1 MG", "S3: 13.5 MG", "4A: 7.75 MG", "3-a: 7.75 MG", "3-b: 7.00 MG", "5: 3.71 MG", "3-c: 2.30 MG")), levels=rev(c("S0", "S1", "S2", "S3", "4A", "3-a", "3-b", "5", "3-c")))



sheltering_frequency_plot <- ggplot() +
  geom_bar(data=shetlering_frequency, aes(x=n, y=n_genes, fill=Sheltering_gametolog_LOF), color=scales::alpha("black", 1), position_dodge2(preserve = "single"), stat="identity") +
  scale_fill_manual(values = c("Z"="#b30000", "W"="#fecc5c")) +
  facet_wrap(Species~Strata, nrow=2, scales="free_y") +
  scale_x_continuous(limits=c(0,11), breaks=c(1, 5, 9)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5), labels = scales::label_number(accuracy = 1)) +
  labs(x="Guarding-Frequency Spectrum", y="Number of masking genes", fill="Guarding\ngametolog") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=18, face="bold"), #change legend title font size
        legend.text = element_text(size=16), #change legend text font size
        strip.background = element_rect(color="black", fill="white", linewidth=1),
        strip.text = element_text(size=16, color="black"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        axis.text.y = element_text(size=16, color="black"),
        axis.text.x = element_text(size=16, color="black"),
        axis.ticks.y = element_blank())


shetlering_frequency$Strata2 <- factor(shetlering_frequency$Strata, order=T, labels=rev(c("Ancestral: 13.5-32.3 MG", "Ancestral: 13.5-32.3 MG", "Ancestral: 13.5-32.3 MG", "Ancestral: 13.5-32.3 MG", "4A: 7.75 MG", "3-a: 7.75 MG", "3-b: 7.00 MG", "5: 3.71 MG", "3-c: 2.30 MG")), levels=rev(c("S0: 32.3 MG", "S1: 25.0 MG", "S2: 23.1 MG", "S3: 13.5 MG", "4A: 7.75 MG", "3-a: 7.75 MG", "3-b: 7.00 MG", "5: 3.71 MG", "3-c: 2.30 MG")))

sheltering_frequency_plot2 <- ggplot() +
  geom_bar(data=shetlering_frequency, aes(x=n, y=n_genes, fill=Sheltering_gametolog_LOF), color=scales::alpha("black", 1), position_dodge2(preserve = "single"), stat="identity") +
  scale_fill_manual(values = c("Z"="#b30000", "W"="#fecc5c")) +
  facet_wrap(Species~Strata2, nrow=2, scales="free_y") +
  scale_x_continuous(limits=c(0,11), breaks=c(1, 5, 9)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5), labels = scales::label_number(accuracy = 1)) +
  labs(x="Guarding-Frequency Spectrum", y="Number of masking genes", fill="Guarding\ngametolog") +
  theme_bw() +
  theme(legend.position = "right",
        legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=18, face="bold"), #change legend title font size
        legend.text = element_text(size=16), #change legend text font size
        strip.background = element_rect(color="black", fill="white", linewidth=1),
        strip.text = element_text(size=16, color="black"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        axis.text.y = element_text(size=16, color="black"),
        axis.text.x = element_text(size=16, color="black"),
        axis.ticks.y = element_blank())


sheltering_frequency_plot2

# Plot W guarding as heatmap between Raso lark and Skylark




######## Look into which genes are involved in guarding by W across species, and what their corresponding Z function is lost
data_gt4.2 <- data_gt[which(data_gt$Strata != "Autosomal" & data_gt$Strata != "PAR"),]

genes_with_W <- data_gt4.2 |> filter(Sheltering_gametolog_LOF == "W") |> pull(geneID) |> unique()

gene_counts <-  data_gt4.2 |>
  filter(geneID %in% genes_with_W) |>
  group_by(Species, geneID) |>
  summarise(n_LOF_W = sum(Hap2_region == "W" & LOF2 == 1, na.rm = TRUE), n_Hap_W = sum(Hap2_region == "W", na.rm = TRUE),
    n_LOF_W = sum(Hap2_region == "W" & LOF2 == 1 & Sheltering_gametolog_LOF == "W", na.rm = TRUE),
    n_Hap_W = sum(Hap2_region == "W" & Sheltering_gametolog_LOF == "W", na.rm = TRUE),
    n_LOF_Z = sum(Hap1_region == "Z" & LOF1 == 1, na.rm = TRUE) + sum(Hap2_region == "Z" & LOF2 == 1, na.rm = TRUE),
    n_Hap_Z = sum(Hap1_region == "Z", na.rm = TRUE) + sum(Hap2_region == "Z", na.rm = TRUE), .groups = "drop") |>
  mutate(prop_LOF_W = ifelse(n_Hap_W > 0, n_LOF_W / n_Hap_W, 0), prop_LOF_Z = ifelse(n_Hap_Z > 0, n_LOF_Z / n_Hap_Z, 0))


# List genes
matrix_n_LOF_Z_species <- gene_counts |>
  select(Species, geneID, prop_LOF_Z) |>
  pivot_wider(names_from = Species, values_from = prop_LOF_Z, values_fill = 0)

prop_W_LOF2_species <- data_gt4.2 |>
  group_by(Species, geneID) |>
  summarise(Strata = first(Strata),    prop_LOF2_W = ifelse(any(Hap2_region == "W"), mean(LOF2[Hap2_region == "W"], na.rm = TRUE), NA_real_), .groups = "drop")

prop_W_wide <- prop_W_LOF2_species |>
  select(Species, geneID, prop_LOF2_W) |>
  pivot_wider(names_from = Species,     values_from = prop_LOF2_W, names_prefix = "propW_")

strata_df <- prop_W_LOF2_species |> select(geneID, Strata) |> distinct()

merged_df <- matrix_n_LOF_Z_species |>
  left_join(prop_W_wide, by = "geneID") |>
  left_join(strata_df, by = "geneID")

# Plot heat map
gene_counts2 <- data_gt4.2 |>
  group_by(Species, geneID) |>
  filter(any(Sheltering_gametolog_LOF == "W")) |>
  summarise(n_LOF_W = sum(Hap2_region == "W" & LOF2 == 1, na.rm = TRUE), n_Hap_W = sum(Hap2_region == "W", na.rm = TRUE),
            n_LOF_Z = sum(Hap1_region == "Z" & LOF1 == 1, na.rm = TRUE) + sum(Hap2_region == "Z" & LOF2 == 1, na.rm = TRUE),
            n_Hap_Z = sum(Hap1_region == "Z", na.rm = TRUE) + sum(Hap2_region == "Z", na.rm = TRUE), .groups = "drop") |>
  mutate(prop_LOF_W = ifelse(n_Hap_W > 0, n_LOF_W / n_Hap_W, 0), prop_LOF_Z  = ifelse(n_Hap_Z > 0, n_LOF_Z / n_Hap_Z, 0))

gene_wide <- gene_counts2 |>
  select(Species, geneID, prop_LOF_Z) |>
  pivot_wider(names_from = Species, values_from = prop_LOF_Z, values_fill = 0)
colnames(gene_wide) <- c("geneID", "Raso lark", "Skylark")

pairwise_W_shelt <- ggplot(gene_wide, aes(x = `Raso lark`, y = Skylark)) +
  geom_bin2d(bins = 10) +   # counts of overlapping points
  scale_fill_viridis_c(option = "viridis") +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  labs(x = "Raso lark Z loss-of-function frequency", y = "Skylark Z loss-of-function frequency", fill = "Guarding by\nW count") +
  theme_bw() +
    theme(legend.position = "right",
          legend.key.size = unit(1, 'line'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=18, face="bold"), #change legend title font size
          legend.text = element_text(size=16), #change legend text font sizes
          strip.background = element_rect(color="black", fill="white", linewidth=1),
          strip.text = element_text(size=16, color="black"),
          axis.line = element_line(colour = "black"),
          axis.title.y = element_text(size=18),
          axis.title.x = element_text(size=18),
          axis.text.y = element_text(size=16, color="black"),
          axis.text.x = element_text(size=16, color="black"),
          axis.ticks.y = element_blank())
pairwise_W_shelt

# Plot as function of age
newdata4 <- expand.grid(
  Species = unique(data_gt4$Species),
  Strata_Age_Generations = seq(0, max(data_gt4$Strata_Age_Generations)*100)/100,
  GERP_max = median(data_gt4$GERP_max),
  pHaplo = median(data_gt4$pHaplo),
  logGeneLen = median(data_gt4$logGeneLen))
newdata4$id <- seq_len(nrow(newdata4))
pred <- predict(Model4, newdata = newdata4, type = "response", se.fit = TRUE)
newdata4$mean_prob  <- pred$fit
newdata4$lower_prob <- plogis(qlogis(pred$fit) - 1.96 * pred$se.fit)
newdata4$upper_prob <- plogis(qlogis(pred$fit) + 1.96 * pred$se.fit)

# Get proportions
gt_prop4 <- data_gt4 |>
  group_by(Individual, Strata_Age_Generations, Species, Sheltering_gametolog_LOF) |>
  summarise(n = n(), .groups = "drop") |>
  group_by(Individual, Strata_Age_Generations, Species) |>
  mutate(prop = n / sum(n)) |>
  ungroup() |>
  group_by(Individual, Species) |>
  complete(Strata_Age_Generations, Sheltering_gametolog_LOF, fill = list(n = 0, prop = 0)) |>
  ungroup()

gt_prop4W <- gt_prop4[which(gt_prop4$Sheltering_gametolog_LOF == "W"),]
gt_prop4W_strata <- gt_prop4W |>
  group_by(Strata_Age_Generations, Species) |>
  summarise(
    mean_prop = mean(prop, na.rm = TRUE),
    q025 = quantile(prop, 0.025, na.rm = TRUE),
    q975 = quantile(prop, 0.975, na.rm = TRUE),
    .groups = "drop")

set.seed(123)
gt_prop4W_jit <- gt_prop4W |>
  mutate(x_jit = Strata_Age_Generations + runif(n(), -0.2, 0.2))


sheltering_age <- ggplot() +
  geom_ribbon(data = newdata4, aes(x = Strata_Age_Generations, y = mean_prob, ymin = lower_prob, ymax = upper_prob, group=Species, fill=Species), alpha = 0.5) +
  geom_line(data=newdata4, aes(x=Strata_Age_Generations*1.02, y=mean_prob*1.001, group=Species, linetype=Species), color="black", linewidth=1.5) +
  geom_line(data=newdata4, aes(x=Strata_Age_Generations*0.98, y=mean_prob*0.995, group=Species, linetype=Species), color="black", linewidth=1.5) +
  geom_line(data=newdata4, aes(x=Strata_Age_Generations*1.02, y=mean_prob*0.995, group=Species, linetype=Species), color="black", linewidth=1.5) +
  geom_line(data=newdata4, aes(x=Strata_Age_Generations*0.98, y=mean_prob*1.001, group=Species, linetype=Species), color="black", linewidth=1.5) +
  geom_line(data = newdata4, aes(x = Strata_Age_Generations, y = mean_prob, group=Species, linetype=Species, color=Species), size = 1.5) +
  geom_point(data=gt_prop4W_jit, aes(x=x_jit, y=prop), color="black", size=4, alpha = 0.25)+
  geom_point(data=gt_prop4W_jit, aes(x=x_jit, y=prop, color=Species), size=3, alpha = 0.5) +
  scale_color_manual(values = c("Raso lark" = "#92c5de", "Skylark" = "#f4a582")) +
  scale_fill_manual(values = c("Raso lark" = "#92c5de", "Skylark" = "#f4a582")) +
  scale_linetype_manual(values = c("Raso lark"=1, "Skylark"=2)) +
  scale_x_continuous(sec.axis=sec_axis(LOF~ (. * coef(agefit)[1]) + (coef(agefit)[2] * (. ^2)), name="Age (million years in larks)", breaks=seq(0,140,20)), breaks=c(0,5,10,15,20,25,30,35), expand=c(0,0), limits=c(0, 35)) +
   scale_y_continuous(breaks=c(0, 0.1, 0.2, 0.3, 0.4, 0.5)) +
  coord_cartesian(ylim = c(0, 0.5), clip = "on") +
  labs(x = "Age (million generations)", y = expression(atop("Probability", "Functional W guarding Z loss-of-function")), color="Species") +
  theme_bw() +
  theme(legend.position = "right",
        legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(2, 'cm'), #change legend key width
        legend.title = element_text(size=18, face="bold"), #change legend title font size
        legend.text = element_text(size=16), #change legend text font size
        panel.spacing = unit(3, "lines"),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        axis.text.y = element_text(size=16, color="black"),
        axis.text.x = element_text(size=16, color="black"))

sheltering_age


fig5 <- (sheltering_frequency_plot2 / ((pairwise_W_shelt | sheltering_age) + plot_layout(widths = c(0.9, 2.1)))) +
  plot_layout(heights = c(2, 1), guides = "collect") + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 20, face = "bold"), plot.tag.position = c(0.01, 0.95), legend.position = "right")


png("Rasolark_2021/Results/Figures/Figure5.png", width=6000, height=4800, res=300)
fig5
dev.off()




# Plot as function of age and interaction with GERP
#pHaplo_vals <- seq(0, 1.0, 0.2)
GERP_vals   <- quantile(data_gt4$GERP_max, probs=c(0.05, 0.5, 0.95)) #seq(min(data_gt4$GERP_max), max(data_gt4$GERP_max), length.out = 5)
#GERP_vals   <-  quantile(data_gt4$GERP_max, probs=c(0.05, 0.5, 0.95))
GeneLen_vals <- log10(c(300, 1500, 5000, 10000))
newdata4 <- expand.grid(
  Species = "Raso lark",
  Strata_Age_Generations = seq(0, max(data_gt4$Strata_Age_Generations)*100)/100,
  GERP_max = GERP_vals,
  pHaplo = median(data_gt4$pHaplo),
  logGeneLen = GeneLen_vals)
newdata4$id <- seq_len(nrow(newdata4))
pred <- predict(Model4, newdata = newdata4, type = "response", se.fit = TRUE)
newdata4$mean_prob  <- pred$fit
newdata4$lower_prob <- plogis(qlogis(pred$fit) - 1.96 * pred$se.fit)
newdata4$upper_prob <- plogis(qlogis(pred$fit) + 1.96 * pred$se.fit)
#newdata4$GERP_max <- factor(newdata4$GERP_max, order=T, labels=c("Q0.05 (maximum GERP score) = 0.743", "Q0.50 (maximum GERP score) = 1.05", "Q0.975 (maximum GERP score) = 1.13"), levels=c(0.743, 1.050, 1.130))

newdata4$GeneLen <- 10^newdata4$logGeneLen
newdata4$GeneLen <- factor(newdata4$GeneLen, order=T, label=c("300 bp", "1,500 bp", "5,000 bp", "10,000 bp"), levels=c(300, 1500, 5000, 10000))

sheltering_length_GERP <- ggplot() +
  geom_ribbon(data = newdata4, aes(x = Strata_Age_Generations, y = mean_prob, ymin = lower_prob, ymax = upper_prob, group = GERP_max, fill=GERP_max), alpha = 0.5) +
  geom_line(data = newdata4, aes(x = Strata_Age_Generations, y = mean_prob, group = GERP_max, color=GERP_max), size = 1.5) +
  facet_wrap(~GeneLen) +
  scale_color_viridis_c(option = "plasma") +
  scale_fill_viridis_c(option = "plasma") +
  scale_linetype_manual(values = c("Raso lark"=1, "Skylark"=2)) +
  scale_x_continuous(sec.axis=sec_axis(LOF~ (. * coef(agefit)[1]) + (coef(agefit)[2] * (. ^2)), name="Age (million years in larks)", breaks=seq(0,140,20)), breaks=c(0,5,10,15,20,25,30,35), expand=c(0,0), limits=c(0, 35)) +
  scale_y_continuous(breaks=c(0, 0.2, 0.5, 0.75, 1.00)) +
  coord_cartesian(ylim = c(0, 1), clip = "on") +
  labs(x = "Age (million generations)", y = expression(atop("Probability", "Functional W guarding Z loss-of-function")), color=expression(atop("Maximum", "GERP score"))) +
  guides(fill="none") +
  theme_bw() +
  theme(legend.position = "right",
        legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(2, 'cm'), #change legend key width
        legend.title = element_text(size=15, face="bold"), #change legend title font size
        legend.text = element_text(size=15), #change legend text font size
        strip.placement = "inside",
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(size=14, face="bold"),
        panel.spacing = unit(3, "lines"),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        axis.text.y = element_text(size=15, color="black"),
        axis.text.x = element_text(size=15, color="black"))



png("Rasolark_2021/Results/Figures/Supplementary/sheltering_length_GERP.png", width=3000, height=2500, res=300)
sheltering_length_GERP
dev.off()




# Plot pHaplo x GERP x length
pHaplo_vals <- seq(0, 1.0, 0.2)
GERP_vals   <- seq(min(data_gt4$GERP_max), max(data_gt4$GERP_max), length.out = 100)
GeneLen_vals <- log10(c(300, 1500, 5000, 10000))
newdata <- expand.grid(
  Species = "Raso lark",
    GERP_max = GERP_vals,
    pHaplo = pHaplo_vals,
    logGeneLen = GeneLen_vals,
    Strata_Age_Generations = median(data_gt2$Strata_Age_Generations, na.rm = TRUE))
newdata$id <- seq_len(nrow(newdata))
pred <- predict(Model4, newdata = newdata,type = "response", se.fit = TRUE)
newdata$mean_prob  <- pred$fit
newdata$lower_prob <- plogis(qlogis(pred$fit) - 1.96 * pred$se.fit)
newdata$upper_prob <- plogis(qlogis(pred$fit) + 1.96 * pred$se.fit)

rect_df <- data.frame(
  xmin = min(data_gt4$GERP_max),
  xmax = quantile(data_gt4$GERP_max, probs = 0.05),
  ymin = 0,
  ymax = 1,
  Species = unique(prob_func_hom$Species))

newdata$GeneLen <- 10^newdata$logGeneLen
newdata$GeneLen <- factor(newdata$GeneLen, order=T, label=c("300 bp", "1,500 bp", "5,000 bp", "10,000 bp"), levels=c(300, 1500, 5000, 10000))

sheltering_pHaplo_GERP <- ggplot() +
  geom_ribbon(data = newdata, aes(x = GERP_max, y = mean_prob, ymin = lower_prob, ymax = upper_prob, group = pHaplo, fill=pHaplo), alpha = 0.5) +
  geom_line(data = newdata, aes(x = GERP_max, y = mean_prob, group = pHaplo, color=pHaplo), size = 1.5) +
  geom_vline(data = data.frame(xintercept = median(data_gt3$GERP_max)), aes(xintercept = xintercept), color="#d7191c", linewidth=1, linetype=2) +
  geom_rect(data = rect_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), inherit.aes = FALSE, alpha = 0.5) +
  facet_wrap(~GeneLen) +
  scale_color_viridis_c(option = "plasma") +
  scale_fill_viridis_c(option = "plasma") +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(limits=c(0,1), expand=c(0,0)) +
  labs(x = "Maximum GERP score", y = expression(atop("Probability", "Functional W guarding Z loss-of-function")), color = "pHaplo") +
  theme_bw() +
  theme(legend.position = "right",
        legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(2, 'cm'), #change legend key width
        legend.title = element_text(size=15, face="bold"), #change legend title font size
        legend.text = element_text(size=15), #change legend text font size
        strip.background = element_blank(),
        strip.text = element_text(size=14, face="bold"),
        panel.spacing = unit(3, "lines"),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        axis.text.y = element_text(size=15, color="black"),
        axis.text.x = element_text(size=15, color="black"))

sheltering_pHaplo_GERP

png("Rasolark_2021/Results/Figures/Supplementary/sheltering_pHaplo_GERP.png", width=3000, height=2500, res=300)
sheltering_pHaplo_GERP
dev.off()



