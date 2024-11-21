rm(list = ls())
options(stringsAsFactors = FALSE)

library(reshape2)
library(ggplot2)
library(ggridges)
library(ggpubr)

rawdata <- read.csv("circle.bed", sep = "\t", header = F)
data <- rawdata[rawdata$V1 %in% c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrY"),]
data$length <- data$V3 - data$V2 - 1

# ======================= Young =======================
ggdensity(data, x="length", color = "V13", 
          alpha = 0, palette = c("#339900"))+ 
  xlim(0, 45000) + theme(legend.position="right") +
  labs(title="3 month old cortex", x="Circle size (bp)", y="Density") +
  theme(axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12))+
  theme(legend.position = 'none')

# ======================= Adult =======================
ggdensity(data, x="length", color = "V13",
          alpha = 0, palette = c("#FF9900"))+
  xlim(0, 45000) +  
  theme(legend.position="right") +
  labs(title="12 month old cortex", x="Circle size (bp)", y="Density") +
  theme(axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12))+
  theme(legend.position = 'none')

# ======================= Old =======================
ggdensity(data, x="length", color = "V13",
          alpha = 0, palette = c("#006699"))+
  xlim(0, 45000) + 
  theme(legend.position="right") +
  labs(title="22 month old cortex", x="Circle size (bp)", y="Density") +
  theme(axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12))+
  theme(legend.position = 'none')
