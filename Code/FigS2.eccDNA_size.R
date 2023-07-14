# Supplementary Fig2 - 11_circle_length

rm(list = ls())
options(stringsAsFactors = FALSE)

# Set working folder
mainDir <- "/Volumes/Denmark/Mouse_Atlas/Bioinformatics/Figures/Fig1"
subDir <- "11_circle_length"
# create a directory if it doesn't exist
ifelse(!dir.exists(file.path(mainDir, subDir)), dir.create(file.path(mainDir, subDir)), FALSE)
setwd(file.path(mainDir, subDir))
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

# ======================= Cortex =======================
# 3m ctx
ctx_3m <- data[c(data$V14=="Y" & data$V13=="CT"),]
# 12m ctx
ctx_12m <- data[c(data$V14=="A" & data$V13=="CT"),]
# 22m ctx
ctx_22m <- data[c(data$V14=="O" & data$V13=="CT"),]

# 3m plot
ggdensity(ctx_3m, x="length", color = "V13",
               alpha = 0, palette = c("#339900"))+
  xlim(0, 45000) + 
  theme(legend.position="right") +
  labs(title="3 month old cortex", x="Circle size (bp)", y="Density") +
  theme(axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12))+
  theme(legend.position = 'none')

ggsave("11_circle_length_3m_cortex.pdf", width = 5, height = 3)


# 12m plot
ggdensity(ctx_12m, x="length", color = "V13",
          alpha = 0, palette = c("#FF9900"))+
  xlim(0, 45000) +  
  theme(legend.position="right") +
  labs(title="12 month old cortex", x="Circle size (bp)", y="Density") +
  theme(axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12))+
  theme(legend.position = 'none')

ggsave("11_circle_length_12m_cortex.pdf", width = 5, height = 3)


# 22m plot
ggdensity(ctx_22m, x="length", color = "V13",
          alpha = 0, palette = c("#006699"))+
  xlim(0, 45000) + 
  theme(legend.position="right") +
  labs(title="22 month old cortex", x="Circle size (bp)", y="Density") +
  theme(axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12))+
  theme(legend.position = 'none')

ggsave("11_circle_length_22m_cortex.pdf", width = 5, height = 3)

# ======================= Hippocampus =======================
# 3m hipp
hipp_3m <- data[c(data$V14=="Y" & data$V13=="HP"),]
# 12m hipp
hipp_12m <- data[c(data$V14=="A" & data$V13=="HP"),]
# 22m hipp
hipp_22m <- data[c(data$V14=="O" & data$V13=="HP"),]

# 3m plot
ggdensity(hipp_3m, x="length", color = "V13",
          alpha = 0, palette = c("#339900"))+
  xlim(0, 45000) + 
  theme(legend.position="right") +
  labs(title="3 month old hippocampus", x="Circle size (bp)", y="Density") +
  theme(axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12))+
  theme(legend.position = 'none')

ggsave("11_circle_length_3m_hippocampus.pdf", width = 5, height = 3)


# 12m plot
ggdensity(hipp_12m, x="length", color = "V13",
          alpha = 0, palette = c("#FF9900"))+
  xlim(0, 45000) + 
  theme(legend.position="right") +
  labs(title="12 month old hippocampus", x="Circle size (bp)", y="Density") +
  theme(axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12))+
  theme(legend.position = 'none')

ggsave("11_circle_length_12m_hippocampus.pdf", width = 5, height = 3)


# 22m plot
ggdensity(hipp_22m, x="length", color = "V13",
          alpha = 0, palette = c("#006699"))+
  xlim(0, 45000) + 
  theme(legend.position="right") +
  labs(title="22 month old hippocampus", x="Circle size (bp)", y="Density") +
  theme(axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12))+
  theme(legend.position = 'none')

ggsave("11_circle_length_22m_hippocampus.pdf", width = 5, height = 3)


# ======================= Liver =======================
# 3m liver
liver_3m <- data[c(data$V14=="Y" & data$V13=="LV"),]
# 12m liver
liver_12m <- data[c(data$V14=="A" & data$V13=="LV"),]
# 22m liver
liver_22m <- data[c(data$V14=="O" & data$V13=="LV"),]

# 3m plot
ggdensity(liver_3m, x="length", color = "V13",
          alpha = 0, palette = c("#339900"))+
  xlim(0, 45000) + 
  theme(legend.position="right") +
  labs(title="3 month old liver", x="Circle size (bp)", y="Density") +
  theme(axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12))+
  theme(legend.position = 'none')

ggsave("11_circle_length_3m_liver.pdf", width = 5, height = 3)


# 12m plot
ggdensity(liver_12m, x="length", color = "V13",
          alpha = 0, palette = c("#FF9900"))+
  xlim(0, 45000) + 
  theme(legend.position="right") +
  labs(title="12 month old liver", x="Circle size (bp)", y="Density") +
  theme(axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12))+
  theme(legend.position = 'none')

ggsave("11_circle_length_12m_liver.pdf", width = 5, height = 3)


# 22m plot
ggdensity(liver_22m, x="length", color = "V13",
          alpha = 0, palette = c("#006699"))+
  xlim(0, 45000) + 
  theme(legend.position="right") +
  labs(title="22 month old liver", x="Circle size (bp)", y="Density") +
  theme(axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12))+
  theme(legend.position = 'none')

ggsave("11_circle_length_22m_liver.pdf", width = 5, height = 3)



# ======================= Pancreas =======================
# 3m pan
pan_3m <- data[c(data$V14=="Y" & data$V13=="PA"),]
# 12m pan
pan_12m <- data[c(data$V14=="A" & data$V13=="PA"),]
# 22m pan
pan_22m <- data[c(data$V14=="O" & data$V13=="PA"),]

# 3m plot
ggdensity(pan_3m, x="length", color = "V13",
          alpha = 0, palette = c("#339900"))+
  xlim(0, 45000) + 
  theme(legend.position="right") +
  labs(title="3 month old pancreas", x="Circle size (bp)", y="Density") +
  theme(axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12))+
  theme(legend.position = 'none')

ggsave("11_circle_length_3m_pancreas.pdf", width = 5, height = 3)


# 12m plot
ggdensity(pan_12m, x="length", color = "V13",
          alpha = 0, palette = c("#FF9900"))+
  xlim(0, 45000) + 
  theme(legend.position="right") +
  labs(title="12 month old pancreas", x="Circle size (bp)", y="Density") +
  theme(axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12))+
  theme(legend.position = 'none')

ggsave("11_circle_length_12m_pancreas.pdf", width = 5, height = 3)


# 22m plot
ggdensity(pan_22m, x="length", color = "V13",
          alpha = 0, palette = c("#006699"))+
  xlim(0, 45000) + 
  theme(legend.position="right") +
  labs(title="22 month old pancreas", x="Circle size (bp)", y="Density") +
  theme(axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12))+
  theme(legend.position = 'none')

ggsave("11_circle_length_22m_pancreas.pdf", width = 5, height = 3)



# ======================= Muscle =======================
# 3m muscle
muscle_3m <- data[c(data$V14=="Y" & data$V13=="MS"),]
# 12m muscle
muscle_12m <- data[c(data$V14=="A" & data$V13=="MS"),]
# 22m muscle
muscle_22m <- data[c(data$V14=="O" & data$V13=="MS"),]

# 3m plot
ggdensity(muscle_3m, x="length", color = "V13",
          alpha = 0, palette = c("#339900"))+
  xlim(0, 45000) + 
  theme(legend.position="right") +
  labs(title="3 month old muscle", x="Circle size (bp)", y="Density") +
  theme(axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12))+
  theme(legend.position = 'none')

ggsave("11_circle_length_3m_muscle.pdf", width = 5, height = 3)


# 12m plot
ggdensity(muscle_12m, x="length", color = "V13",
          alpha = 0, palette = c("#FF9900"))+
  xlim(0, 45000) + 
  theme(legend.position="right") +
  labs(title="12 month old muscle", x="Circle size (bp)", y="Density") +
  theme(axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12))+
  theme(legend.position = 'none')

ggsave("11_circle_length_12m_muscle.pdf", width = 5, height = 3)


# 22m plot
ggdensity(muscle_22m, x="length", color = "V13",
          alpha = 0, palette = c("#006699"))+
  xlim(0, 45000) + 
  theme(legend.position="right") +
  labs(title="22 month old muscle", x="Circle size (bp)", y="Density") +
  theme(axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12))+
  theme(legend.position = 'none')

ggsave("11_circle_length_22m_muscle.pdf", width = 5, height = 3)


# ======================= VAT(visceral adipose tissue) =======================
# 3m VAT
VAT_3m <- data[c(data$V14=="Y" & data$V13=="VF"),]
# 12m VAT
VAT_12m <- data[c(data$V14=="A" & data$V13=="VF"),]
# 22m VAT
VAT_22m <- data[c(data$V14=="O" & data$V13=="VF"),]

# 3m plot
ggdensity(VAT_3m, x="length", color = "V13",
          alpha = 0, palette = c("#339900"))+
  xlim(0, 45000) + 
  theme(legend.position="right") +
  labs(title="3 month old VAT", x="Circle size (bp)", y="Density") +
  theme(axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12))+
  theme(legend.position = 'none')

ggsave("11_circle_length_3m_VAT.pdf", width = 5, height = 3)


# 12m plot
ggdensity(VAT_12m, x="length", color = "V13",
          alpha = 0, palette = c("#FF9900"))+
  xlim(0, 45000) + 
  theme(legend.position="right") +
  labs(title="12 month old VAT", x="Circle size (bp)", y="Density") +
  theme(axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12))+
  theme(legend.position = 'none')

ggsave("11_circle_length_12m_VAT.pdf", width = 5, height = 3)


# 22m plot
ggdensity(VAT_22m, x="length", color = "V13",
          alpha = 0, palette = c("#006699"))+
  xlim(0, 45000) + 
  theme(legend.position="right") +
  labs(title="22 month old VAT", x="Circle size (bp)", y="Density") +
  theme(axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12))+
  theme(legend.position = 'none')

ggsave("11_circle_length_22m_VAT.pdf", width = 5, height = 3)



# ======================= SAT(subcutaneous adipose tissue) =======================
# 3m SAT
SAT_3m <- data[c(data$V14=="Y" & data$V13=="SF"),]
# 12m SAT
SAT_12m <- data[c(data$V14=="A" & data$V13=="SF"),]
# 22m SAT
SAT_22m <- data[c(data$V14=="O" & data$V13=="SF"),]

# 3m plot
ggdensity(SAT_3m, x="length", color = "V13",
          alpha = 0, palette = c("#339900"))+
  xlim(0, 45000) + 
  theme(legend.position="right") +
  labs(title="3 month old SAT", x="Circle size (bp)", y="Density") +
  theme(axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12))+
  theme(legend.position = 'none')

ggsave("11_circle_length_3m_SAT.pdf", width = 5, height = 3)


# 12m plot
ggdensity(SAT_12m, x="length", color = "V13",
          alpha = 0, palette = c("#FF9900"))+
  xlim(0, 45000) + 
  theme(legend.position="right") +
  labs(title="12 month old SAT", x="Circle size (bp)", y="Density") +
  theme(axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12))+
  theme(legend.position = 'none')

ggsave("11_circle_length_12m_SAT.pdf", width = 5, height = 3)


# 22m plot
ggdensity(SAT_22m, x="length", color = "V13",
          alpha = 0, palette = c("#006699"))+
  xlim(0, 45000) + 
  theme(legend.position="right") +
  labs(title="22 month old SAT", x="Circle size (bp)", y="Density") +
  theme(axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12))+
  theme(legend.position = 'none')

ggsave("11_circle_length_22m_SAT.pdf", width = 5, height = 3)


# ======================= Body (no head) =======================
# embryo body
embryo_body <- data[c(data$V14=="E" & data$V13=="BD"),]
ggdensity(embryo_body, x="length", color = "V13",
          alpha = 0, palette = c("#339900"))+
  xlim(0, 45000) + 
  theme(legend.position="right") +
  labs(title="Embryo body", x="Circle size (bp)", y="Density") +
  theme(axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12))+
  theme(legend.position = 'none')

ggsave("11_circle_length_embryo_body.pdf", width = 5, height = 3)
