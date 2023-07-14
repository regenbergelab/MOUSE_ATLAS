# Fig1f - Correlation between chromosome length and eccDNA counts

rm(list = ls())
options(stringsAsFactors = FALSE)

# Set working folder
mainDir <- "/Volumes/Denmark/Mouse_Atlas/Bioinformatics/Figures/Fig1"
subDir <- "21_chrLength_eccNumber"
# create a directory if it doesn't exist
ifelse(!dir.exists(file.path(mainDir, subDir)), dir.create(file.path(mainDir, subDir)), FALSE)
setwd(file.path(mainDir, subDir))

# Load required libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)

# Load required table
rawdata <- read.csv("../03_circle_number/inputfile/highMappingRate_clean_merged_145_medium_conf_labeled_circle.bed", sep = "\t", header = F)

# Import the length of each chromosome (mm10)
chrLen <- read.csv("./mm10_chrLen.csv")
chrLen

# --------------- Young --------------- 
data <- rawdata[which(rawdata$V14 == "Y"),]

# Calculate how many ecDNAs are dropped on each chromosome
circle_num_per_chr <- as.data.frame(table(data$V1))
colnames(circle_num_per_chr) <- c("chrID", "eccNum")
circle_num_per_chr

# merge two data frames
input <- merge(circle_num_per_chr, chrLen, by = "chrID", all = FALSE)
input$chrLen2 <- input$chrLen/100000000
input$eccNum2 <- input$eccNum/1000
input$chrID2 <- c("1","10","11","12","13","14","15","16","17","18","19","2","3","4","5","6","7","8","9","X","Y")
input

# Plot the correlation between chromosome length and circle number
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

ggsave(filename = "21_chrLen_eccNum_Young.pdf", width = 4, height = 4)


# --------------- Adult --------------- 
data <- rawdata[which(rawdata$V14 == "A"),]

# Calculate how many ecDNAs are dropped on each chromosome
circle_num_per_chr <- as.data.frame(table(data$V1))
colnames(circle_num_per_chr) <- c("chrID", "eccNum")
circle_num_per_chr

# merge two data frames
input <- merge(circle_num_per_chr, chrLen, by = "chrID", all = FALSE)
input$chrLen2 <- input$chrLen/100000000
input$eccNum2 <- input$eccNum/1000
input$chrID2 <- c("1","10","11","12","13","14","15","16","17","18","19","2","3","4","5","6","7","8","9","X","Y")
input

# Plot the correlation between chromosome length and circle number
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

ggsave(filename = "21_chrLen_eccNum_Adult.pdf", width = 4, height = 4)


# --------------- Old --------------- 
data <- rawdata[which(rawdata$V14 == "O"),]

# Calculate how many ecDNAs are dropped on each chromosome
circle_num_per_chr <- as.data.frame(table(data$V1))
colnames(circle_num_per_chr) <- c("chrID", "eccNum")
circle_num_per_chr

# merge two data frames
input <- merge(circle_num_per_chr, chrLen, by = "chrID", all = FALSE)
input$chrLen2 <- input$chrLen/100000000
input$eccNum2 <- input$eccNum/1000
input$chrID2 <- c("1","10","11","12","13","14","15","16","17","18","19","2","3","4","5","6","7","8","9","X","Y")
input

# Plot the correlation between chromosome length and circle number
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

ggsave(filename = "21_chrLen_eccNum_Old.pdf", width = 4, height = 4)
