rm(list = ls())
options(stringsAsFactors = FALSE)

library(quantreg)
library(splines)
library(data.table)
library(dplyr)
library(ggplot2)
library(rGREAT)
library(ggpubr)
library(dplyr)

# --------------- Young --------------- 
df <- data.frame(group = c("Non-transcribed genes","Transcribed genes"), value = c(741, 5197))
df <- df %>% arrange(desc(group)) %>% mutate(pct = round(value / sum(df$value), 4) *100)
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


# --------------- Adult --------------- 
df <- data.frame(group = c("Non-transcribed genes","Transcribed genes"),  value = c(604, 4353))
df <- df %>% arrange(desc(group)) %>% mutate(pct = round(value / sum(df$value), 4) *100)
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

# --------------- Old --------------- 
df <- data.frame(group = c("Non-transcribed genes","Transcribed genes"), value = c(308, 1343))
df <- df %>% arrange(desc(group)) %>% mutate(pct = round(value / sum(df$value), 4) *100)
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
