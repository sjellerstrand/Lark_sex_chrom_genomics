## Export variables and load libraries
rm(list=ls())
library(tidyverse)
library(patchwork)
options(scipen=999)
setwd("C:/Users/Simon JE/OneDrive - Lund University/Dokument/Simon/PhD/Projects/Rasolark_2021/Results/")
#setwd("C:/Users/si0630ja/OneDrive - Lund University/Dokument/Simon/PhD/Projects/Rasolark_2021/Results/")


## Get sample data
Skylark_inds <- read.table("C:/Users/Simon JE/OneDrive - Lund University/Dokument/Simon/PhD/Projects/Rasolark_2021/Results/Figures/Skylark_2021_sample_info.txt", sep="\t", head=T)
Skylark_inds$Population[which(Skylark_inds$Population == "IT")] <- "Italy"
Skylark_inds$Population[which(Skylark_inds$Population == "NL")] <- "Netherlands"

### PCA
# Import data
data_SApc <- read.table(paste(getwd(), "/Skylark_HWE/PCA/autosomal/Skylark_2021_autosomal.eigenvec", sep=""))
colnames(data_SApc) <- c("ID", "ID2", paste0("PC", 1:(ncol(data_SApc)-2)))
data_SApc$Sex <- rep(NA, nrow(data_SApc))
data_SApc$Population <- rep(NA, nrow(data_SApc))
data_SApc$Data <- "Autosomes"
data_SAeig <- t(read.table(paste(getwd(), "/Skylark_HWE/PCA/autosomal/Skylark_2021_autosomal.eigenval", sep="")))
colnames(data_SAeig) <- c(paste0("PC", 1:(ncol(data_SAeig))))
for(i in 1:nrow(data_SApc)) {
  data_SApc$Sex[i] <- Skylark_inds$Sex[which(Skylark_inds$ID == data_SApc$ID[i])]
  data_SApc$Population[i] <- Skylark_inds$Population[which(Skylark_inds$ID == data_SApc$ID[i])]
}

plot_SA <- ggplot() +
  geom_point(data=data_SApc, aes(x=PC1, y=PC2, color=Population, shape=Sex), position="identity", size=4) +
  scale_color_manual(values = c("#91cf60", "#fc8d59"), limits=c("Italy", "Netherlands")) +
  facet_wrap(~Data) +
  labs(x=paste("PC1 (", sprintf("%.2f", data_SAeig[1]), "%)", sep=""), y=paste("PC2 (", sprintf("%.2f", data_SAeig[2]), "%)", sep="")) +
  scale_y_continuous(expand = c(0.05,0.05)) +
  scale_x_continuous(expand = c(0.05,0.05)) +
  theme_bw() +
  theme(legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=20), #change legend title font size
        legend.text = element_text(size=20), #change legend text font size
        panel.spacing = unit(1, "lines"),
        strip.background = element_rect(color="black", fill="white", linewidth=1),
        strip.text.x = element_text(size=20, color="black"),
        strip.text.y = element_text(size=20, color="black"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=15, color="black"),
        axis.text.x =  element_text(size=15, color="black"))

png("Figures/Supplementary/Skylark_PCA.png", width=3000, height=2000, res=300)
plot_SA
dev.off()


### Admixture CV error

# Import data
data_SAcv <- read.table(paste(getwd(), "/Skylark_HWE/admixture/autosomal/Skylark_2021_autosomal_CVerror", sep=""))
colnames(data_SAcv) <- c("K", "CVerror")
data_SAcv$Data <- "Autosomal"


plot_SA_cv <- ggplot() +
  geom_point(data=data_SAcv, aes(x=K, y=CVerror), position="identity", size=4) +
  labs(x="K", y="CV error") +
  scale_y_continuous(expand = c(0.05,0.05)) +
  scale_x_continuous(expand = c(0.05,0.05)) +
  theme_bw() +
  theme(legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=20), #change legend title font size
        legend.text = element_text(size=20), #change legend text font size
        panel.spacing = unit(1, "lines"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=15, color="black"),
        axis.text.x =  element_text(size=15, color="black"))


### Admixture K=2

# Import data
data_SAK2 <- read.table(paste(getwd(), "/Skylark_HWE/admixture/autosomal/Skylark_2021_autosomal.2.Q", sep=""))
colnames(data_SAK2) <- c("Cluster 1", "Cluster 2")
ind_orderSA <- read.table(paste(getwd(), "/Skylark_HWE/admixture/autosomal/Skylark_2021_autosomal.fam", sep=""))
data_SAK2$ID <- ind_orderSA$V1
data_SAK2$Sex <- rep(NA, nrow(data_SAK2))
data_SAK2$Population <- rep(NA, nrow(data_SAK2))
for(i in 1:nrow(data_SAK2)) {
  data_SAK2$Sex[i] <- Skylark_inds$Sex[which(Skylark_inds$ID == data_SAK2$ID[i])]
  data_SAK2$Population[i] <- Skylark_inds$Population[which(Skylark_inds$ID == data_SAK2$ID[i])]
}
data_SAK2 <- data_SAK2 |> pivot_longer(cols=c("Cluster 1", "Cluster 2"), names_to="Cluster", values_to="Proportion")



plot_SA_ad <- ggplot() +
  geom_bar(data=data_SAK2, aes(x=ID, y=Proportion, fill=Cluster), stat="identity") +
  facet_grid(~ Population + Sex, scales = "free", space="free") +
  scale_fill_manual(values = c("#440154FF", "#FDE725FF")) +
  labs(x = "Sex & Location", title = "Autosomal SNPs, K=2", y = "Proportion ancestry") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme_bw() +
  theme(legend.position="none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), 
      panel.spacing = unit(0, "lines"),
      strip.background = element_rect(color="black", fill="white", linewidth=1),
      strip.text = element_text(size=20, color="black"),
      axis.line = element_line(colour = "black"),
      plot.title = element_text(size=20),
      axis.title.y = element_text(size=20),
      axis.title.x = element_text(size=20),
      axis.text.y = element_text(size=15, color="black"),
      axis.text.x = element_text(size=15, color="black", angle=90),
      axis.ticks.x=element_blank())


png("Figures/Supplementary/Skylark_Admixture.png", width=3000, height=6000, res=300)
plot_SA_ad / plot_SA_cv +
    plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 20, face="bold"), plot.tag.position=c(0,0.95))
dev.off()
  
