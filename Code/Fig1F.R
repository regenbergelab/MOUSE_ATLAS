rm(list = ls())
options(stringsAsFactors = FALSE)

library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)

rawdata <- read.csv("circle.bed", sep = "\t", header = F)
chrLen <- read.csv("mm10_chrLen.csv")

# --------------- Young --------------- 
data <- rawdata[which(rawdata$V14 == "Y"),]
circle_num_per_chr <- as.data.frame(table(data$V1))
colnames(circle_num_per_chr) <- c("chrID", "eccNum")
input <- merge(circle_num_per_chr, chrLen, by = "chrID", all = FALSE)
input$chrLen2 <- input$chrLen/100000000
input$eccNum2 <- input$eccNum/1000
input$chrID2 <- c("1","10","11","12","13","14","15","16","17","18","19","2","3","4","5","6","7","8","9","X","Y")
p <- ggscatter(input,  x = "chrLen2", y = "eccNum2",
          add = "reg.line",  # Add Fit Curve
          conf.int = T,  # Add confidence interval
          add.params = list(color = "#339900", size = 1, fill = "#CCCCCC"), # Customize reg.line
          size = 1, color = "black", shape = 21, fill = "black")  + 
  stat_cor(method = "pearson", color = "black")+  # Adding correlation coefficients and p-values
  theme(legend.position="none") + 
  theme_classic()+
  labs(title="3 month", x = "Chromosome length (x 10^8 bp)", y = "eccDNA counts (x 1000)") 
p + geom_text_repel(label = input$chrID2)

# --------------- Adult --------------- 
data <- rawdata[which(rawdata$V14 == "A"),]
circle_num_per_chr <- as.data.frame(table(data$V1))
colnames(circle_num_per_chr) <- c("chrID", "eccNum")
input <- merge(circle_num_per_chr, chrLen, by = "chrID", all = FALSE)
input$chrLen2 <- input$chrLen/100000000
input$eccNum2 <- input$eccNum/1000
input$chrID2 <- c("1","10","11","12","13","14","15","16","17","18","19","2","3","4","5","6","7","8","9","X","Y")
p <- ggscatter(input,  x = "chrLen2", y = "eccNum2",
               add = "reg.line",  # Add Fit Curve
               conf.int = T,  # Add confidence interval
               add.params = list(color = "#FF9900", size = 1, fill = "#CCCCCC"), # Customize reg.line
               size = 1, color = "black", shape = 21, fill = "black")  + 
  stat_cor(method = "pearson", color = "black")+  # Adding correlation coefficients and p-values
  theme(legend.position="none") + 
  theme_classic()+
  labs(title="12 month", x = "Chromosome length (x 10^8 bp)", y = "eccDNA counts (x 1000)") 
p + geom_text_repel(label = input$chrID2)

# --------------- Old --------------- 
data <- rawdata[which(rawdata$V14 == "O"),]
circle_num_per_chr <- as.data.frame(table(data$V1))
colnames(circle_num_per_chr) <- c("chrID", "eccNum")
input <- merge(circle_num_per_chr, chrLen, by = "chrID", all = FALSE)
input$chrLen2 <- input$chrLen/100000000
input$eccNum2 <- input$eccNum/1000
input$chrID2 <- c("1","10","11","12","13","14","15","16","17","18","19","2","3","4","5","6","7","8","9","X","Y")
p <- ggscatter(input,  x = "chrLen2", y = "eccNum2",
               add = "reg.line",  # Add Fit Curve
               conf.int = T,  # Add confidence interval
               add.params = list(color = "#006699", size = 1, fill = "#CCCCCC"), # Customize reg.line
               size = 1, color = "black", shape = 21, fill = "black")  + 
  stat_cor(method = "pearson", color = "black")+  # Adding correlation coefficients and p-values
  theme(legend.position="none") + 
  theme_classic()+
  labs(title="22 month", x = "Chromosome length (x 10^8 bp)", y = "eccDNA counts (x 1000)") 
p + geom_text_repel(label = input$chrID2)
