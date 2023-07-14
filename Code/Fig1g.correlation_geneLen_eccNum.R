# Fig1g - Correlation between gene Length with eccDNA counts

rm(list = ls())
options(stringsAsFactors = FALSE)

# Set working folder
mainDir <- "/Volumes/Denmark/Mouse_Atlas/Bioinformatics/Figures/Fig1"
subDir <- "06_geneLength_eccNumber"
# create a directory if it doesn't exist
ifelse(!dir.exists(file.path(mainDir, subDir)), dir.create(file.path(mainDir, subDir)), FALSE)
setwd(file.path(mainDir, subDir))

# Load required libraries
library(splines)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)

# Load required Rdata
load("./inputfile/anno_eccDNA.RData")

# --------------- Young --------------- 
data <- rawdata[which(rawdata$V12 == "FS174"|rawdata$V12 == "FS177"|rawdata$V12 == "FS193"|rawdata$V12 == "FS185"|rawdata$V12 == "FS213"|rawdata$V12 == "FS201"),]

# Calculate the length of each protein-coding gene per MB
data$gene_length_MB <- (data$V17 - data$V16 + 1)/1000000
data$gene_length <- data$V17 - data$V16 + 1

# Calculate how many eccDNAs are dropped on each gene
circle_num_per_gene <- as.data.frame(table(data$V18))
colnames(circle_num_per_gene) <- c("gene_name", "ecc_number")
gene_len_MB <- data[,c("V18", "gene_length_MB")]
gene_len_MB_dup <- gene_len_MB[!duplicated(gene_len_MB),]
colnames(gene_len_MB_dup) <- c("gene_name", "gene_length_MB")
circle_num_per_gene_len <- merge(circle_num_per_gene, gene_len_MB_dup, by = "gene_name", all = FALSE)
circle_num_per_gene_len$ecc_number_per_MB <- circle_num_per_gene_len$ecc_number/circle_num_per_gene_len$gene_length_MB

# Import normalized count matrix
norm_counts <- read.csv("/Volumes/Denmark/Mouse_Atlas/Bioinformatics/RNAseq/Liver_01/normalized_counts.txt", sep = "\t", header = T)

single_norm_counts <- subset(norm_counts, select = c(X, s57,s61,s80,s69,s111,s89)) # Young age group
single_norm_counts$Young <- single_norm_counts$s57 + single_norm_counts$s61 + single_norm_counts$s80 + single_norm_counts$s69 + single_norm_counts$s111 + single_norm_counts$s89
single_norm_counts <- single_norm_counts[,c("X","Young")]
colnames(single_norm_counts) <- c("gene_name", "norm_count")
circle_num_per_gene_len_NormCount <- merge(circle_num_per_gene_len, single_norm_counts, by = "gene_name", all = FALSE)
input <- circle_num_per_gene_len_NormCount
input$norm_count_per_MB <- input$norm_count/input$gene_length_MB

# Plot the correlation between gene length and circle number
ggscatter(input,  x = "gene_length_MB", y = "ecc_number",
          add = "reg.line",  # Add Fit Curve
          conf.int = T,  # Add confidence interval
          add.params = list(color = "#CCCCCC", size = 0.5, fill = "#CCCCCC"), # Customize reg.line
          size = 1, color = "#339900", shape = 21, fill = "#339900")  + 
  stat_cor(method = "pearson", color = "black")+  # Adding correlation coefficients and p-values
  theme(legend.position="none") + 
  theme_classic()+
  labs(title="Young", x = "Gene length per MB", y = "EccDNA number") 
ggsave(filename = "06_geneLength_eccNumber_Young.pdf", width = 4, height = 4)

# --------------- Adult --------------- 
data <- rawdata[which(rawdata$V12 == "FS139"|rawdata$V12 == "FS141"|rawdata$V12 == "FS149"|rawdata$V12 == "FS160"|rawdata$V12 == "FS163"),] # Adult age group

# Calculate the length of each protein-coding gene per MB
data$gene_length_MB <- (data$V17 - data$V16 + 1)/1000000
data$gene_length <- data$V17 - data$V16 + 1

# Calculate how many eccDNAs are dropped on each gene
circle_num_per_gene <- as.data.frame(table(data$V18))
colnames(circle_num_per_gene) <- c("gene_name", "ecc_number")
gene_len_MB <- data[,c("V18", "gene_length_MB")]
gene_len_MB_dup <- gene_len_MB[!duplicated(gene_len_MB),]
colnames(gene_len_MB_dup) <- c("gene_name", "gene_length_MB")
circle_num_per_gene_len <- merge(circle_num_per_gene, gene_len_MB_dup, by = "gene_name", all = FALSE)
circle_num_per_gene_len$ecc_number_per_MB <- circle_num_per_gene_len$ecc_number/circle_num_per_gene_len$gene_length_MB

# Import normalized count matrix
norm_counts <- read.csv("./inputfile/normalized_counts.txt", sep = "\t", header = T)

single_norm_counts <- subset(norm_counts, select = c(X, s19,s21,s29,s40,s43)) # Adult age group
single_norm_counts$Adult <- single_norm_counts$s19 + single_norm_counts$s21 + single_norm_counts$s29 + single_norm_counts$s40 + single_norm_counts$s43
single_norm_counts <- single_norm_counts[,c("X","Adult")]
colnames(single_norm_counts) <- c("gene_name", "norm_count")
circle_num_per_gene_len_NormCount <- merge(circle_num_per_gene_len, single_norm_counts, by = "gene_name", all = FALSE)
input <- circle_num_per_gene_len_NormCount
input$norm_count_per_MB <- input$norm_count/input$gene_length_MB

# Plot the correlation between gene length and circle number
ggscatter(input,  x = "gene_length_MB", y = "ecc_number",
          add = "reg.line",  # Add Fit Curve
          conf.int = T,  # Add confidence interval
          add.params = list(color = "#CCCCCC", size = 0.5, fill = "#CCCCCC"), # Customize reg.line
          size = 1, color = "#FF9900", shape = 21, fill = "#FF9900")  + 
  stat_cor(method = "pearson", color = "black")+  # Adding correlation coefficients and p-values
  theme(legend.position="none") + 
  theme_classic()+
  labs(title="Adult", x = "Gene length per MB", y = "EccDNA number") 
ggsave(filename = "06_geneLength_eccNumber_Adult.pdf", width = 4, height = 4)

# --------------- Old --------------- 
data <- rawdata[which(rawdata$V12 == "FS254"|rawdata$V12 == "FS289"),] # Old age group

# Calculate the length of each protein-coding gene per MB
data$gene_length_MB <- (data$V17 - data$V16 + 1)/1000000
data$gene_length <- data$V17 - data$V16 + 1

# Calculate how many eccDNAs are dropped on each gene
circle_num_per_gene <- as.data.frame(table(data$V18))
colnames(circle_num_per_gene) <- c("gene_name", "ecc_number")
gene_len_MB <- data[,c("V18", "gene_length_MB")]
gene_len_MB_dup <- gene_len_MB[!duplicated(gene_len_MB),]
colnames(gene_len_MB_dup) <- c("gene_name", "gene_length_MB")
circle_num_per_gene_len <- merge(circle_num_per_gene, gene_len_MB_dup, by = "gene_name", all = FALSE)
circle_num_per_gene_len$ecc_number_per_MB <- circle_num_per_gene_len$ecc_number/circle_num_per_gene_len$gene_length_MB

# Import normalized count matrix
norm_counts <- read.csv("/Volumes/Denmark/Mouse_Atlas/Bioinformatics/RNAseq/Liver_01/normalized_counts.txt", sep = "\t", header = T)

single_norm_counts <- subset(norm_counts, select = c(X, s162,s200)) # Old age group
single_norm_counts$Old <- single_norm_counts$s162 + single_norm_counts$s200
single_norm_counts <- single_norm_counts[,c("X","Old")]
colnames(single_norm_counts) <- c("gene_name", "norm_count")
circle_num_per_gene_len_NormCount <- merge(circle_num_per_gene_len, single_norm_counts, by = "gene_name", all = FALSE)
input <- circle_num_per_gene_len_NormCount
input$norm_count_per_MB <- input$norm_count/input$gene_length_MB

# Plot the correlation between gene length and circle number
ggscatter(input,  x = "gene_length_MB", y = "ecc_number",
          add = "reg.line",  # Add Fit Curve
          conf.int = T,  # Add confidence interval, Default is 95% confidence interval. 
          add.params = list(color = "#CCCCCC", size = 0.5, fill = "#CCCCCC"), # Customize reg.line
          size = 1, color = "#006699", shape = 21, fill = "#006699")  + 
  stat_cor(method = "pearson", color = "black")+  # Adding correlation coefficients and p-values
  theme(legend.position="none") + 
  theme_classic()+
  labs(title="Old", x = "Gene length per MB", y = "EccDNA number") 
ggsave(filename = "06_geneLength_eccNumber_Old.pdf", width = 4, height = 4)
