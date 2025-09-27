## Load libraries
rm(list=ls())
library(tidyverse)
library(patchwork)

options(scipen=999)
setwd("C:/Users/Simon JE/OneDrive - Lund University/Dokument/Simon/PhD/Projects/Skylark_2021/Results/Figures")

# Skylark genome summary
data7 <- read.delim("C:/Users/Simon JE/OneDrive - Lund University/Dokument/Simon/PhD/Projects/Skylark_2021/Results/Coordinate lift over/Skylark_2021_Genome_summary.bed", sep="\t", head=F)
colnames(data7) <- c("scaffold", "start", "end", "data_type")
data7 <- data7[which(data7$scaffold == "CADDXX010000137.1"),]
data7$start <- as.numeric(data7$start)
data7$end <- as.numeric(data7$end)
data7$data_type[which(data7$data_type == "Sex phase & depth difference")] <- "Sex haplotype & coverage difference"
data7$data_type[which(data7$data_type == "Sex depth difference")] <- "Sex sequencing\ndepth difference"
data7$data_type[which(data7$data_type == "Sex phase difference")] <- "Sex haplotype clustering"
data7$data_type[which(data7$data_type == "Missing data")] <- "No data"

data8 <- read.delim("C:/Users/Simon JE/OneDrive - Lund University/Dokument/Simon/PhD/Projects/Skylark_2021/Results/Coordinate lift over/Rasolark_2021_Genome_summary.bed", sep="\t", head=F)
colnames(data8) <- c("scaffold", "start", "end", "data_type")
data8 <- data8[which(data8$scaffold == "CADDXX010000137.1"),]
data8$start <- as.numeric(data8$start)
data8$end <- as.numeric(data8$end)
data8$data_type[which(data8$data_type == "Sex phase & depth difference")] <- "Sex haplotype & coverage difference"
data8$data_type[which(data8$data_type == "Sex depth difference")] <-"Sex sequencing\ndepth difference"
data8$data_type[which(data8$data_type == "Sex phase difference")] <- "Sex haplotype clustering"
data8$data_type[which(data8$data_type == "Missing data")] <- "No data"


# Dxy
data_dxy <- read.table("C:/Users/Simon JE/OneDrive - Lund University/Dokument/Simon/PhD/Projects/Skylark_2021/Results/PAR/Rasolark_2021_Skylark_2021_CADDXX010000137.1_dxy_all_data_window_2000_step_500.txt", sep="\t", head=T)
colnames(data_dxy) <- c("scaffold", "start", "end", "mid", "sites", "pi1", "pi2", "dxy", "Fst", "data_type")

data_dxy$data_type[which(data_dxy$data_type == "males")] <- "Skylark vs. Raso lark Z & PAR"
data_dxy <- data_dxy[-which(data_dxy$data_type == "females_Rasolark_2021"),]
data_dxy <- data_dxy[-which(data_dxy$data_type == "females_Skylark_2021"),]
data_dxy$data_type[which(data_dxy$data_type == "females_W")] <- "Skylark vs. Raso lark W"

data_dxy <- data_dxy[which(!is.na(data_dxy$dxy)),]
data_dxy$mid <- as.numeric(data_dxy$mid)
data_dxy$dxy <- as.numeric(data_dxy$dxy)
data_dxy$data_type <- factor(data_dxy$data_type, order=T, levels=c("Skylark vs. Raso lark W", "Skylark vs. Raso lark Z & PAR"))

# LD
data_LD <- read.delim("C:/Users/Simon JE/OneDrive - Lund University/Dokument/Simon/PhD/Projects/Skylark_2021/Results/PAR/Rasolark_2021_Skylark_2021_CADDXX010000137.1_LD_decay_all_data_window_5000_step_1250.txt", sep="\t", head=T)

data_LD <- data_LD[-which(data_LD$data_type == "Rasolark_2021_males"),]
data_LD$data_type[which(data_LD$data_type == "Rasolark_2021_females")] <- "Raso lark females"
data_LD <- data_LD[-which(data_LD$data_type == "Skylark_2021_males"),]
data_LD$data_type[which(data_LD$data_type == "Skylark_2021_females")] <- "Skylark females"
data_LD <- data_LD[-(which(data_LD$data_type == "data_type")),]

data_LD <- data_LD[which(!is.na(data_LD$r2)),]
data_LD$r2 <- as.numeric(data_LD$r2)
data_LD$mid <- as.numeric(data_LD$mid)
data_LD$data_type <- factor(data_LD$data_type, order=T, levels=c("Raso lark females", "Skylark females"))

# Gene features
data_genes <- read.delim("C:/Users/Simon JE/OneDrive - Lund University/Dokument/Simon/PhD/Projects/Skylark_2021/Results/PAR/CADDXX010000137.1.gtf", sep="\t", head=F)
colnames(data_genes) <- c("Scaffold", "Method", "Feature", "Start", "End", "vetej1", "Strandedness", "vetej2", "Attribtues")
data_genes$Gene <- str_split(str_split(data_genes$Attribtues, ";", simplify=T)[,1], " ", simplify=T)[,2]
data_genes <- data_genes[which(data_genes$Feature == "exon" | data_genes$Feature == "CDS" ),]
data_PTP4A1_exons <- data_genes[which(data_genes$Gene == "PTP4A1" & data_genes$Feature == "exon"),]
data_LOC_exons <- data_genes[which(data_genes$Gene == "LOC115494256" & data_genes$Feature == "exon"),]
data_PHF3_exons <- data_genes[which(data_genes$Gene == "PHF3" & data_genes$Feature == "exon"),]
data_PTP4A1_CDS <- data_genes[which(data_genes$Gene == "PTP4A1" & data_genes$Feature == "CDS"),]
data_LOC_CDS <- data_genes[which(data_genes$Gene == "LOC115494256" & data_genes$Feature == "CDS"),]
data_PHF3_CDS <- data_genes[which(data_genes$Gene == "PHF3" & data_genes$Feature == "CDS"),]

# Plot
minx <- 194300
maxx <- 275585


#maxx <- 279300

PAR_classS <- ggplot() +
  # PhaseWY results
  geom_rect(data=data7, aes(xmin=start, xmax=end,  ymin=0, ymax=1, fill=data_type)) +
  scale_fill_manual(values = c("#E4EAF0", "#fecc5c", "#b30000", "#404040"), limits = c("Autosomal", "Sex haplotype clustering", "Sex sequencing\ndepth difference", "No data")) +
  scale_color_manual(values = c("#E4EAF0", "#fecc5c", "#b30000", "#404040"), limits = c("Autosomal", "Sex haplotype clustering", "Sex sequencing\ndepth difference", "No data")) +
  # Horisontal lines
  geom_segment(aes(x=minx, xend=maxx, y=0, yend=0)) +
  geom_segment(aes(x=minx, xend=maxx, y=1, yend=1)) +
  # Legend
  guides(fill=guide_legend(title="PhaseWY classification")) +
  # Other
  ylab(expression(atop("","Skylark"))) +
  scale_x_continuous(limits = c(minx,maxx), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme_bw() +
  theme(
    legend.position.inside = c(1,20),
    legend.title = element_text(size=15),
    legend.text = element_text(size=15),
    legend.key.size = unit(1, units = "cm"),
    axis.title.y = element_text(size=15),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y=element_blank())

PAR_classR <- ggplot() +
  # PhaseWY results
  geom_rect(data=data8, aes(xmin=start, xmax=end,  ymin=0, ymax=1, fill=data_type)) +
  scale_fill_manual(values = c("#E4EAF0", "#fecc5c", "#b30000", "#404040"), limits = c("Autosomal", "Sex haplotype clustering", "Sex sequencing\ndepth difference", "No data")) +
  scale_color_manual(values = c("#E4EAF0", "#fecc5c", "#b30000", "#404040"), limits = c("Autosomal", "Sex haplotype clustering", "Sex sequencing\ndepth difference", "No data")) +
  # Horisontal lines
  geom_segment(aes(x=minx, xend=maxx, y=0, yend=0)) +
  geom_segment(aes(x=minx, xend=maxx, y=1, yend=1)) +
  # Legend
  guides(fill=guide_legend(title="PhaseWY classification")) +
  # Other
  ylab(expression(atop("","Raso lark"))) +
  scale_x_continuous(limits = c(minx,maxx), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme_bw() +
  theme(
    legend.position.inside = c(1,20),
    legend.title = element_text(size=15),
    legend.text = element_text(size=15),
    legend.key.size = unit(1, units = "cm"),
    axis.title.y = element_text(size=15),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y=element_blank())

PAR_genes <- ggplot() +
  # Horisontal lines
  geom_segment(aes(x=min(data_PTP4A1_exons$Start), xend=max(data_PTP4A1_exons$End), y=0.75, yend=0.75), color="#404040") +
  geom_segment(aes(x=min(data_LOC_exons$Start), xend=max(data_LOC_exons$End), y=0.625, yend=0.625), color="#404040") +
  geom_segment(aes(x=minx, xend=max(data_PHF3_exons$End), y=0.25, yend=0.25), color="#404040") +
  # Arrows
  ## PTP4A1
  geom_segment(aes(x=max(data_PTP4A1_exons$End), xend=239350, y=0.75, yend=0.75), color="#404040" ,  arrow = arrow(length = unit(0.3, "cm"))) +
  geom_segment(aes(x=max(data_PTP4A1_exons$End), xend=232400, y=0.75, yend=0.75), color="#404040" ,  arrow = arrow(length = unit(0.3, "cm"))) +
  ## LOC
  geom_segment(aes(x=min(data_LOC_exons$Start), xend=209000, y=0.625, yend=0.625), color="#404040" ,  arrow = arrow(length = unit(0.3, "cm"))) +
  ##PHF3
  geom_segment(aes(x=max(data_PHF3_exons$End), xend=206650, y=0.25, yend=0.25), color="#404040" ,  arrow = arrow(length = unit(0.3, "cm"))) +
  geom_segment(aes(x=max(data_PHF3_exons$End), xend=216900, y=0.25, yend=0.25), color="#404040" ,  arrow = arrow(length = unit(0.3, "cm"))) +
  # Exons
  geom_rect(data=data_PTP4A1_exons, aes(xmin=Start, xmax=End, ymin=0.625, ymax=0.875, fill="Non-coding"), color="#404040") +
  geom_rect(data=data_LOC_exons, aes(xmin=Start, xmax=End, ymin=0.50, ymax=0.75, fill="Non-coding"), color="#404040") +
  geom_rect(data=data_PHF3_exons, aes(xmin=Start, xmax=End, ymin=0.125, ymax=0.375, fill="Non-coding"), color="#404040") +
  # CDS
  geom_rect(data=data_PTP4A1_CDS, aes(xmin=Start, xmax=End, ymin=0.625, ymax=0.875, fill="Coding"), color="#404040",) +
  geom_rect(data=data_LOC_CDS, aes(xmin=Start, xmax=End, ymin=0.50, ymax=0.75, fill="Coding"), color="#404040") +
  geom_rect(data=data_PHF3_CDS, aes(xmin=Start, xmax=End, ymin=0.125, ymax=0.375, fill="Coding"), color="#404040") +
  scale_fill_manual(values = c("#7D7D7D", "#404040"), limits = c("Non-coding", "Coding")) +
  # Gene labels
  annotate(geom="text", x=max(data_PTP4A1_exons$End)+3000 , y=0.75, label="PTP4A1", color="Black") +
  annotate(geom="text", x=max(data_LOC_exons$End)+5200, y=0.625, label="LOC115494256", color="Black") +
  annotate(geom="text", x=max(data_PHF3_exons$End)+2300, y=0.25, label="PHF3", color="Black") +
  # Legend
  guides(fill=guide_legend(title="Exons")) +
  # Other
  scale_x_continuous(limits = c(minx,maxx), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  ylab(expression(atop("","Genes"))) + 
  theme_bw() +
  theme(legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=15), #change legend title font size
        legend.text = element_text(size=15), #change legend text font size
        panel.spacing = unit(1, "lines"),
        axis.title.y = element_text(size=15),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y=element_blank())

PAR_dxy <- ggplot() +
  geom_line(data=data_dxy, aes(x=mid, y=dxy, group=data_type, color=data_type, linetype=data_type), linewidth=1.2) +
  scale_color_manual(values = c("#ca0020", "#0571b0")) +
  scale_linetype_manual(values = c(6,6)) +
  # Legend
  guides(color=guide_legend(title="Absolute divergence"), linetype=guide_legend(title="Absolute divergence")) +
  # Other
  ylab(expression(atop("Absolute", "divergence (D"[XY]*")"))) +
  scale_x_continuous(limits = c(minx,maxx), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0.02,0.02)) +
  theme_bw() +
  theme(legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=15), #change legend title font size
        legend.text = element_text(size=15), #change legend text font size
        panel.spacing = unit(1, "lines"),
        axis.title.y = element_text(size=15),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=10, color="black"),
        axis.text.x = element_blank())

PAR_LD <- ggplot() + 
  geom_line(data=data_LD, aes(x=mid, y=r2, group=data_type, color=data_type, linetype=data_type), linewidth=1.2) +
  scale_color_manual(values = c("#92c5de", "#f4a582")) +
  scale_linetype_manual(values = c(1,2)) +
  # Legend
  guides(color=guide_legend(title="Linkage disequilibrium"), linetype=guide_legend(title="Linkage disequilibrium")) +
  # Other
  xlab("Scaffold CADDXX010000137.1\nposition (bp)") +
  ylab(expression(atop("Linkage","disequilibrium ( r"^2*")"))) +
  scale_x_continuous(limits = c(minx,maxx), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,0.7), expand = c(0.02,0.0), breaks= c(0.0, 0.2, 0.4, 0.6)) +
  theme_bw() +
  theme(legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=15), #change legend title font size
        legend.text = element_text(size=15), #change legend text font size
        panel.spacing = unit(1, "lines"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=10, color="black"),
        axis.text.x =  element_text(size=10, color="black"))


jpeg("Figure2.jpg", width=4500, height=2600, res=300)
((PAR_classS / PAR_classR ) + plot_layout(guides = "collect")) /  PAR_genes / PAR_dxy / PAR_LD + 
  plot_layout(heights = c(0.5, 0.5, 0.5, 1, 1)) + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 20, face="bold"), plot.tag.position=c(0.00,0.95))
dev.off()
