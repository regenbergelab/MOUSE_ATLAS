# Fig2e - piechart of circularized gene ratio with and without transcription

rm(list = ls())
options(stringsAsFactors = FALSE)

# Set working folder
mainDir <- "/Volumes/Denmark/Mouse_Atlas/Bioinformatics/Figures/Fig1"
subDir <- "08_piechart_circularized_gene_ratio"
# create a directory if it doesn't exist
ifelse(!dir.exists(file.path(mainDir, subDir)), dir.create(file.path(mainDir, subDir)), FALSE)
setwd(file.path(mainDir, subDir))
# Load required libraries
library(quantreg)
library(splines)
library(data.table)
library(dplyr)
library(ggplot2)
library(rGREAT)
library(ggpubr)
library(dplyr)

# --------------- Young --------------- 
df <- data.frame(group = c("Non-transcribed genes","Transcribed genes"), 
                 value = c(741, 5197))

# Compute the percentage of gene number
df <- df %>% arrange(desc(group)) %>% mutate(pct = round(value / sum(df$value), 4) *100)

# Pie chart
ggplot(df, aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start = 0)+
  theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold"))+
  scale_fill_manual(values=c("#CCCCCC", "#339900")) +
  theme(axis.text.x=element_blank()) + # Remove axis tick mark labels
  theme(legend.position="bottom")+
  geom_text(aes(x = 1.7, 
                label = paste0(pct, "%")), 
                position = position_stack(vjust = 0.5), 
                size = 3)+
  ggtitle("Young")
ggsave(filename = "08_piechart_circularized_gene_ratio_with@without_transcription_Young.pdf", width = 4, height = 4)


# --------------- Adult --------------- 
df <- data.frame(group = c("Non-transcribed genes","Transcribed genes"), 
                 value = c(604, 4353))

# Compute the percentage of gene number
df <- df %>% arrange(desc(group)) %>% mutate(pct = round(value / sum(df$value), 4) *100)

# Pie chart
ggplot(df, aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start = 0)+
  theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold"))+
  scale_fill_manual(values=c("#CCCCCC", "#FF9900")) +
  theme(axis.text.x=element_blank()) + # Remove axis tick mark labels
  theme(legend.position="bottom")+
  geom_text(aes(x = 1.7, 
                label = paste0(pct, "%")), 
            position = position_stack(vjust = 0.5), 
            size = 3)+
  ggtitle("Adult")
ggsave(filename = "08_piechart_circularized_gene_ratio_with@without_transcription_Adult.pdf", width = 4, height = 4)


# --------------- Old --------------- 
df <- data.frame(group = c("Non-transcribed genes","Transcribed genes"), 
                 value = c(308, 1343))

# Compute the percentage of gene number
df <- df %>% arrange(desc(group)) %>% mutate(pct = round(value / sum(df$value), 4) *100)


# Pie chart
ggplot(df, aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start = 0)+
  theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold"))+
  scale_fill_manual(values=c("#CCCCCC", "#006699")) +
  theme(axis.text.x=element_blank()) + # Remove axis tick mark labels
  theme(legend.position="bottom")+
  geom_text(aes(x = 1.7, 
                label = paste0(pct, "%")), 
            position = position_stack(vjust = 0.5), 
            size = 3)+
  ggtitle("Old")
ggsave(filename = "08_piechart_circularized_gene_ratio_with@without_transcription_Old.pdf", width = 4, height = 4)


