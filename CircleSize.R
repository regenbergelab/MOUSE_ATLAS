# circle size

rm(list = ls())
options(stringsAsFactors = FALSE)

# Load required libraries
library(reshape2)
library(ggplot2)
library(ggridges)
library(ggpubr)

# Import annotated bed file
rawdata <- read.csv("./inputfile/highMappingRate_clean_merged_145_medium_conf_labeled_circle.bed", sep = "\t", header = F)

# Filter out eccDNAs of non-chromosomal origin
data <- rawdata[rawdata$V1 %in% c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrY"),]
table(data$V1)

# Calculate circle length (for bed format, start is 0-based; end is 1-based)
data$length <- data$V3 - data$V2 -1
head(data)

# Plot circle length distribution - all
data$V14 <- factor(data$V14, levels = c("E","Y","A","O"))
p <- ggdensity(data, x="length", #color = "V14",
               alpha = 0, #palette = c("#666666", "#339900", "#FF9900", "#006699")
)+
  geom_vline(xintercept = c(2350), color="darkgrey", linetype="dashed", linewidth=0.5) +
  xlim(0, 20000) + 
  theme(legend.position="right") +
  labs(title="", x="Circle size (bp)", y="Density") +
  theme(axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12)) 
p
ggsave("11_circle_length_all.pdf", width = 6, height = 4)

# Plot circle length distribution by ages
data$V14 <- factor(data$V14, levels = c("E","Y","A","O"))
p <- ggdensity(data, x="length", color = "V14",
          alpha = 0, palette = c("#666666", "#339900", "#FF9900", "#006699")
          )+
  xlim(0, 20000) + 
  theme(legend.position="right") +
  labs(title="", x="Circle size (bp)", y="Density") +
  theme(axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12)) 
p
ggsave("11_circle_length_ages.pdf", width = 6, height = 4)

# Plot circle length distribution by tissues
data$V13 <- factor(data$V13, levels = c("BD","VF","LV","MS","SF","HP","PA","CT"))
p <- ggdensity(data, x="length", color = "V13",
               alpha = 0, palette = c("#85302c", "#cf864b", "#ddbc6e", "#a5b971", "#78b692", "#7db1d4", "#4980b0", "#595976")
)+
  xlim(0, 20000) + 
  theme(legend.position="right") +
  labs(title="", x="Circle size (bp)", y="Density") +
  theme(axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12)) 
p
ggsave("11_circle_length_tissues.pdf", width = 6, height = 4)
