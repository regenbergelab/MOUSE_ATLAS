rm(list = ls())
options(stringsAsFactors = FALSE)

library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(dplyr)

data <- read.csv("circle.bed", sep = "\t", header = F)

df <- data[data$V1 %in% c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrY"),]

df$ecc_Len <- df$V3 - df$V2

summary_df <- df %>%
  mutate(
    category = case_when(
      ecc_Len < 2000 ~ "< 2000",
      ecc_Len >= 2000 & ecc_Len <= 10000 ~ "2000 - 10000",
      ecc_Len > 10000 ~ "> 10000"
    )
  ) %>%
  group_by(category) %>%
  summarise(
    count = n(),
    proportion = n() / nrow(df)
  )

my_colors <- c("#fde8c6","#fbc979", "#e09148")

pie_chart <- ggplot(summary_df, aes(x = "", y = proportion, fill = category)) +
  geom_bar(stat = "identity", width = 1, color = "black") +  
  coord_polar("y", start = 0) +
  labs(title = "Ecc_Len Distribution", fill = "Category") +
  scale_fill_manual(values = my_colors) +
  theme_void()

df <- data.frame(ecc_Len = hh1$ecc_Len)
cuts <- c(0, 2000, 10000)
colors <- c("#fde8c6", "#fbc979")

histogram1 <- ggplot(df, aes(x = ecc_Len)) +
  geom_histogram(aes(y = ..density..), bins = 1000, fill = "skyblue", color = "#fbc979", alpha = 0.7) +  
  labs(title = "Histogram of ecc_Len", x = "ecc_Len", y = "Density") +
  theme_classic()+
  xlim(0, 10000)

histogram2 <- ggplot(df, aes(x = ecc_Len)) +
  geom_histogram(aes(y = ..density..), bins = 1000, fill = "skyblue", color = "#fde8c6", alpha = 0.7) +  
  labs(title = "Histogram of ecc_Len", x = "ecc_Len", y = "Density") +
  theme_classic(
    
  )+
  xlim(0, 10000)+
  geom_vline(xintercept = 2000, color = "black", linetype = "solid")
histogram3 <- ggplot(df, aes(x = ecc_Len)) +
  geom_histogram(aes(y = ..density..), bins = 1000, color = "#fde8c6", alpha = 0.7) +  
  labs(title = "Histogram of ecc_Len", x = "ecc_Len", y = "Density") +
  theme_classic()+
  xlim(100, 2000)+
  geom_vline(xintercept = c(200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800), color = "gray", linetype = "dashed")
