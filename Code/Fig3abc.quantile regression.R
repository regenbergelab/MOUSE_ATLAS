# quantile regression

rm(list = ls())
options(stringsAsFactors = FALSE)

# Set working folder
mainDir <- "/Volumes/Denmark/Mouse_Atlas/Bioinformatics/Figures/Fig1"
subDir <- "07_spline_regression"
# create a directory if it doesn't exist
ifelse(!dir.exists(file.path(mainDir, subDir)), dir.create(file.path(mainDir, subDir)), FALSE)
setwd(file.path(mainDir, subDir))

# Load required libraries
library(quantreg)
library(splines)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)

# --------------- 3-month Liver(n=6)- --------------- 
# Load required Rdata
load("../06_geneLength_eccNumber/inputfile/anno_eccDNA.RData")

data <- rawdata[which(rawdata$V12 == "FS174"|rawdata$V12 == "FS177"|rawdata$V12 == "FS193"|rawdata$V12 == "FS185"|rawdata$V12 == "FS213"|rawdata$V12 == "FS201"),]

# Calculate the length of each protein-coding gene per MB
data$gene_length_MB <- (data$V17 - data$V16 + 1)/1000000
data$gene_length <- data$V17 - data$V16 + 1

# Calculate how many ecDNAs are dropped on each gene
circle_num_per_gene <- as.data.frame(table(data$V18))
colnames(circle_num_per_gene) <- c("gene_name", "ecc_number")
gene_len_MB <- data[,c("V18", "gene_length_MB")]
gene_len_MB_dup <- gene_len_MB[!duplicated(gene_len_MB),]
colnames(gene_len_MB_dup) <- c("gene_name", "gene_length_MB")
circle_num_per_gene_len <- merge(circle_num_per_gene, gene_len_MB_dup, by = "gene_name", all = FALSE)
circle_num_per_gene_len$ecc_number_per_MB <- circle_num_per_gene_len$ecc_number/circle_num_per_gene_len$gene_length_MB

# Import normalized count matrix
norm_counts <- read.csv("/Volumes/Denmark/Mouse_Atlas/Bioinformatics/RNAseq/Liver_01/normalized_counts.txt", sep = "\t", header = T)

single_norm_counts <- subset(norm_counts, select = c(X, s57,s61,s80,s69,s111,s89)) 
single_norm_counts$Young <- single_norm_counts$s57 + single_norm_counts$s61 + single_norm_counts$s80 + single_norm_counts$s69 + single_norm_counts$s111 + single_norm_counts$s89
single_norm_counts <- single_norm_counts[,c("X","Young")]

colnames(single_norm_counts) <- c("gene_name", "norm_count")
circle_num_per_gene_len_NormCount <- merge(circle_num_per_gene_len, single_norm_counts, by = "gene_name", all = FALSE)
input <- circle_num_per_gene_len_NormCount
input$norm_count_per_MB <- input$norm_count/input$gene_length_MB

# Take log() for both X-axis and Y-axis
input$log_norm_count_per_MB <-  log(input$norm_count_per_MB)
input$log_ecc_number_per_MB <-  log(input$ecc_number_per_MB)

# Remove -Inf value from the 'log_norm_count_per_MB' column
input_rm_Inf <- input[which(input$log_norm_count_per_MB != "-Inf"),]
df <- input_rm_Inf
df <- subset(df, select =c(gene_name, log_ecc_number_per_MB, log_norm_count_per_MB))
head(df)

# Conducted cubic spline quantile regression
tau_x <- c(0.1, 0.9) 
quantreg_fit <- rq(log_ecc_number_per_MB ~ bs(log_norm_count_per_MB, df=5), data = df, tau = tau_x)
predict_quantreg <- predict.rq(quantreg_fit, df) %>% data.frame
names(predict_quantreg) <- paste0("quantile_", tau_x)
df <- data.frame(df, predict_quantreg)
log_norm_count_per_MB <- subset(df, select =c(log_norm_count_per_MB))
variable_x<-data.frame(log_norm_count_per_MB)
predict_quantreg_plot <- predict.rq(quantreg_fit, variable_x) %>% data.frame
df_quantreg_plot <- data.frame(variable_x$log_norm_count_per_MB, predict_quantreg_plot)
names(df_quantreg_plot) <- c("log_norm_count_per_MB", paste0("quantile_",tau_x))

# Plot Young
ggplot(data = df, aes(x=log_norm_count_per_MB, y=log_ecc_number_per_MB)) +
  geom_point(alpha=0.8, color="gray", size = 0.5) +
  theme_classic() +
  stat_density_2d(alpha = 0.2, geom = "polygon", contour = TRUE, aes(fill = after_stat(level)), colour = "white", bins = 10) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  labs(title = "Young", 
       x = "Log(normalized counts per gene/MB)", 
       y = "Log(eccDNA number per gene/MB)") +
  geom_line(data=df_quantreg_plot, aes(x=log_norm_count_per_MB, y=quantile_0.1), lwd=0.8, linetype = "dashed", color="#339900")+
  geom_line(data=df_quantreg_plot, aes(x=log_norm_count_per_MB, y=quantile_0.9), lwd=0.8, linetype = "dashed", color="#339900")
ggsave("07_spline_regression_Young.pdf", width = 5, height = 5)

# Extract genes equal or less than 10% quantile(label="below")
df$label <- ifelse(df$log_ecc_number_per_MB <= df$quantile_0.1, "below",
                   ifelse(df$log_ecc_number_per_MB >= df$quantile_0.9, "above", "middle"))
write.csv(df, file = "10%_quantile_genelist_liver_young.csv", row.names = FALSE)


# ------------------------------- 12-month Liver (n=5) ------------------------------- 
rm(list=ls())

# Load required Rdata
load("../06_geneLength_eccNumber/inputfile/anno_eccDNA.RData")

data <- rawdata[which(rawdata$V12 == "FS139"|rawdata$V12 == "FS141"|rawdata$V12 == "FS149"|rawdata$V12 == "FS160"|rawdata$V12 == "FS163"),] # Adult age group

# Calculate the length of each protein-coding gene per MB
data$gene_length_MB <- (data$V17 - data$V16 + 1)/1000000
data$gene_length <- data$V17 - data$V16 + 1

# Calculate how many ecDNAs are dropped on each gene
circle_num_per_gene <- as.data.frame(table(data$V18))
colnames(circle_num_per_gene) <- c("gene_name", "ecc_number")
gene_len_MB <- data[,c("V18", "gene_length_MB")]
gene_len_MB_dup <- gene_len_MB[!duplicated(gene_len_MB),]
colnames(gene_len_MB_dup) <- c("gene_name", "gene_length_MB")
circle_num_per_gene_len <- merge(circle_num_per_gene, gene_len_MB_dup, by = "gene_name", all = FALSE)
circle_num_per_gene_len$ecc_number_per_MB <- circle_num_per_gene_len$ecc_number/circle_num_per_gene_len$gene_length_MB

# Import normalized count matrix
norm_counts <- read.csv("/Volumes/Denmark/Mouse_Atlas/Bioinformatics/RNAseq/Liver_01/normalized_counts.txt", sep = "\t", header = T)

single_norm_counts <- subset(norm_counts, select = c(X, s19,s21,s29,s40,s43)) # Adult age group
single_norm_counts$Adult <- single_norm_counts$s19 + single_norm_counts$s21 + single_norm_counts$s29 + single_norm_counts$s40 + single_norm_counts$s43
single_norm_counts <- single_norm_counts[,c("X","Adult")]

colnames(single_norm_counts) <- c("gene_name", "norm_count")
circle_num_per_gene_len_NormCount <- merge(circle_num_per_gene_len, single_norm_counts, by = "gene_name", all = FALSE)
input <- circle_num_per_gene_len_NormCount
input$norm_count_per_MB <- input$norm_count/input$gene_length_MB

# Take log() for both X-axis and Y-axis
input$log_norm_count_per_MB <-  log(input$norm_count_per_MB)
input$log_ecc_number_per_MB <-  log(input$ecc_number_per_MB)

# Remove -Inf value from the 'log_norm_count_per_MB' column
input_rm_Inf <- input[which(input$log_norm_count_per_MB != "-Inf"),]
df <- input_rm_Inf
df <- subset(df, select =c(gene_name, log_ecc_number_per_MB, log_norm_count_per_MB))
head(df)

# Conducted cubic spline quantile regression
tau_x <- c(0.1, 0.9) 
quantreg_fit <- rq(log_ecc_number_per_MB ~ bs(log_norm_count_per_MB, df=5), data = df, tau = tau_x)
predict_quantreg <- predict.rq(quantreg_fit, df) %>% data.frame
names(predict_quantreg) <- paste0("quantile_", tau_x)
df <- data.frame(df, predict_quantreg)
log_norm_count_per_MB <- subset(df, select =c(log_norm_count_per_MB))
variable_x<-data.frame(log_norm_count_per_MB)
predict_quantreg_plot <- predict.rq(quantreg_fit, variable_x) %>% data.frame
df_quantreg_plot <- data.frame(variable_x$log_norm_count_per_MB, predict_quantreg_plot)
names(df_quantreg_plot) <- c("log_norm_count_per_MB", paste0("quantile_",tau_x))

# Plot Adult
ggplot(data = df, aes(x=log_norm_count_per_MB, y=log_ecc_number_per_MB)) +
  geom_point(alpha=0.8, color="gray", size = 0.5) +
  theme_classic() +
  stat_density_2d(alpha = 0.2, geom = "polygon", contour = TRUE, aes(fill = after_stat(level)), colour = "white", bins = 10) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  labs(title = "Adult", 
       x = "Log(normalized counts per gene/MB)", 
       y = "Log(eccDNA number per gene/MB)") +
  geom_line(data=df_quantreg_plot, aes(x=log_norm_count_per_MB, y=quantile_0.1), lwd=0.8, linetype = "dashed", color="#FF9900")+
  geom_line(data=df_quantreg_plot, aes(x=log_norm_count_per_MB, y=quantile_0.9), lwd=0.8, linetype = "dashed", color="#FF9900")
ggsave("07_spline_regression_Adult.pdf", width = 5, height = 5)

# Extract genes equal or less than 10% quantile (label="below")
df$label <- ifelse(df$log_ecc_number_per_MB <= df$quantile_0.1, "below",
                   ifelse(df$log_ecc_number_per_MB >= df$quantile_0.9, "above", "middle"))
write.csv(df, file = "10%_quantile_genelist_liver_adult.csv", row.names = FALSE)


# ------------------------------- 22-month Liver (n=2) ------------------------------- 
rm(list=ls())

# Load required Rdata
load("../06_geneLength_eccNumber/inputfile/anno_eccDNA.RData")

data <- rawdata[which(rawdata$V12 == "FS254"|rawdata$V12 == "FS289"),] # Old age group

# Calculate the length of each protein-coding gene per MB
data$gene_length_MB <- (data$V17 - data$V16 + 1)/1000000
data$gene_length <- data$V17 - data$V16 + 1

# Calculate how many ecDNAs are dropped on each gene
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

# Take log() for both X-axis and Y-axis
input$log_norm_count_per_MB <-  log(input$norm_count_per_MB)
input$log_ecc_number_per_MB <-  log(input$ecc_number_per_MB)

# Remove -Inf value from the 'log_norm_count_per_MB' column
input_rm_Inf <- input[which(input$log_norm_count_per_MB != "-Inf"),]
df <- input_rm_Inf
df <- subset(df, select =c(gene_name, log_ecc_number_per_MB, log_norm_count_per_MB))
head(df)

# Conducted cubic spline quantile regression
tau_x <- c(0.1, 0.9) 
quantreg_fit <- rq(log_ecc_number_per_MB ~ bs(log_norm_count_per_MB, df=5), data = df, tau = tau_x)
predict_quantreg <- predict.rq(quantreg_fit, df) %>% data.frame
names(predict_quantreg) <- paste0("quantile_", tau_x)
df <- data.frame(df, predict_quantreg)
log_norm_count_per_MB <- subset(df, select =c(log_norm_count_per_MB))
variable_x<-data.frame(log_norm_count_per_MB)
predict_quantreg_plot <- predict.rq(quantreg_fit, variable_x) %>% data.frame
df_quantreg_plot <- data.frame(variable_x$log_norm_count_per_MB, predict_quantreg_plot)
names(df_quantreg_plot) <- c("log_norm_count_per_MB", paste0("quantile_",tau_x))

# Plot Old
ggplot(data = df, aes(x=log_norm_count_per_MB, y=log_ecc_number_per_MB)) +
  geom_point(alpha=0.8, color="gray", size = 0.5) +
  theme_classic() +
  stat_density_2d(alpha = 0.2, geom = "polygon", contour = TRUE, aes(fill = after_stat(level)), colour = "white", bins = 10) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  labs(title = "Old", 
       x = "Log(normalized counts per gene/MB)", 
       y = "Log(eccDNA number per gene/MB)") +
  geom_line(data=df_quantreg_plot, aes(x=log_norm_count_per_MB, y=quantile_0.1), lwd=0.8, linetype = "dashed", color="#006699")+
  geom_line(data=df_quantreg_plot, aes(x=log_norm_count_per_MB, y=quantile_0.9), lwd=0.8, linetype = "dashed", color="#006699")
ggsave("07_spline_regression_Old.pdf", width = 5, height = 5)

# Extract genes equal or less than 10% quantile (label="below")
df$label <- ifelse(df$log_ecc_number_per_MB <= df$quantile_0.1, "below",
                   ifelse(df$log_ecc_number_per_MB >= df$quantile_0.9, "above", "middle"))
write.csv(df, file = "10%_quantile_genelist_liver_old.csv", row.names = FALSE)


#------------------------------------------------3-month Ctx (n=5)------------------------------------------------

# ImporteccDNA bed file annotated by protein-coding genes
rawdata <- read.csv("/Volumes/Denmark/Mouse_Atlas/Bioinformatics/RNAseq/Liver_01/anno.bed", sep = "\t", header = F)
head(rawdata)

## FS206 (s103)
## FS232 (s133)
## FS245 (s152)
## FS144 (s154)
## FS129 (s9)

data <- rawdata[which(rawdata$V12 == "FS206"|rawdata$V12 == "FS232"|rawdata$V12 == "FS245"|rawdata$V12 == "FS144"|rawdata$V12 == "FS129"),]

#计算每个蛋白编码基因的长度（单位：MB)
data$gene_length_MB <- (data$V17 - data$V16 + 1)/1000000
data$gene_length <- data$V17 - data$V16 + 1

# 统计每个基因分别与多少个eccDNA有交集
circle_num_per_gene <- as.data.frame(table(data$V18))
colnames(circle_num_per_gene) <- c("gene_name", "ecc_number")
head(circle_num_per_gene)

# 将基因长度新增到circle_num_per_gene的最后一列
gene_len_MB <- data[,c("V18", "gene_length_MB")]
gene_len_MB_dup <- gene_len_MB[!duplicated(gene_len_MB),]
colnames(gene_len_MB_dup) <- c("gene_name", "gene_length_MB")
circle_num_per_gene_len <- merge(circle_num_per_gene, gene_len_MB_dup, by = "gene_name", all = FALSE)
circle_num_per_gene_len$ecc_number_per_MB <- circle_num_per_gene_len$ecc_number/circle_num_per_gene_len$gene_length_MB

# 如果希望将X轴改为log(gene length/MB), Y轴改为log(eccDNA number), 需要增加这一段代码
circle_num_per_gene_len$log_ecc_number <- log(circle_num_per_gene_len$ecc_number)
circle_num_per_gene_len$log_ecc_number_per_MB <- log(circle_num_per_gene_len$ecc_number_per_MB)
circle_num_per_gene_len$log_gene_length_MB <- log(circle_num_per_gene_len$gene_length_MB)
head(circle_num_per_gene_len)

# 读入Normalized后的count矩阵
norm_counts <- read.csv("/Volumes/Denmark/Mouse_Atlas/Bioinformatics/RNAseq/RNA_02/normalized_counts.txt", sep = "\t", header = T)
head(norm_counts)

## FS206 (s103)
## FS232 (s133)
## FS245 (s152)
## FS144 (s154)
## FS129 (s9)
single_norm_counts <- subset(norm_counts, select = c(X, s103,s133,s152,s154,s9)) 
single_norm_counts$Ctx_Y <- single_norm_counts$s103 + single_norm_counts$s133 + single_norm_counts$s152 + single_norm_counts$s154 + single_norm_counts$s9
single_norm_counts <- single_norm_counts[,c("X","Ctx_Y")]

colnames(single_norm_counts) <- c("gene_name", "norm_count")
head(single_norm_counts)

# 将每个基因的normalized_count值添加到circle_num_per_gene_len数据框的最后一列
circle_num_per_gene_len_NormCount <- merge(circle_num_per_gene_len, single_norm_counts, by = "gene_name", all = FALSE)
head(circle_num_per_gene_len_NormCount)

#简化变量名称
input <- circle_num_per_gene_len_NormCount

# 对横坐标进行gene_length_MB归一化
input$norm_count_per_MB <- input$norm_count/input$gene_length_MB

# 对横坐标取log()
input$log_norm_count_per_MB <-  log(input$norm_count_per_MB)

# 从input数据框中去掉所有log_norm_count_per_MB列是-Inf的行，以免这些点干扰相关系数R
input_rm_Inf <- input[which(input$log_norm_count_per_MB != "-Inf"),]

df <- input_rm_Inf
df <- subset(df, select =c(gene_name, log_ecc_number_per_MB, log_norm_count_per_MB))
head(df)

# tau是分位数回归的分位点
tau_x <- c( 0.1, 0.9) 
quantreg_fit <- rq(log_ecc_number_per_MB ~ bs(log_norm_count_per_MB, df=5), data = df, tau = tau_x)
predict_quantreg <- predict.rq(quantreg_fit, df) %>% data.frame
names(predict_quantreg) <- paste0("quantile_", tau_x)
df <- data.frame(df, predict_quantreg)
log_norm_count_per_MB <- subset(df, select =c(log_norm_count_per_MB))
variable_x <- data.frame(log_norm_count_per_MB)
predict_quantreg_plot <- predict.rq(quantreg_fit, variable_x) %>% data.frame
df_quantreg_plot <- data.frame(variable_x$log_norm_count_per_MB, predict_quantreg_plot)
names(df_quantreg_plot) <- c("log_norm_count_per_MB", paste0("quantile_",tau_x))

# Plot
ggplot(data = df, aes(x=log_norm_count_per_MB, y=log_ecc_number_per_MB)) +
  geom_point(alpha=0.8, color="gray", size = 0.5) +
  theme_classic() +
  stat_density_2d(alpha = 0.2, geom = "polygon", contour = TRUE, aes(fill = after_stat(level)), colour = "white", bins = 10) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  labs(title = "3-month cortex", 
       x = "Log(normalized counts per gene/MB)", 
       y = "Log(eccDNA number per gene/MB)") +
  geom_line(data=df_quantreg_plot, aes(x=log_norm_count_per_MB, y=quantile_0.1), lwd=0.8, linetype = "dashed", color="#339900")+
  geom_line(data=df_quantreg_plot, aes(x=log_norm_count_per_MB, y=quantile_0.9), lwd=0.8, linetype = "dashed", color="#339900")
ggsave("13_spline_regression_cortex_young.pdf", width = 5, height = 5)

# 提取小于或等于10%分位数的基因(label="below")
df$label <- ifelse(df$log_ecc_number_per_MB <= df$quantile_0.1, "below",
                   ifelse(df$log_ecc_number_per_MB >= df$quantile_0.9, "above", "middle"))
write.csv(df, file = "10%_quantile_genelist_cortex_young.csv", row.names = FALSE)


#------------------------------------------------12-month Ctx (n=3)------------------------------------------------

# ImporteccDNA bed file annotated by protein-coding genes
rawdata <- read.csv("/Volumes/Denmark/Mouse_Atlas/Bioinformatics/RNAseq/Liver_01/anno.bed", sep = "\t", header = F)
head(rawdata)

## FS217 (s115)
## FS138 (s18)
## FS169 (s52)
data <- rawdata[which(rawdata$V12 == "FS217" | rawdata$V12 == "FS138" | rawdata$V12 == "FS169"),]

#计算每个蛋白编码基因的长度（单位：MB)
data$gene_length_MB <- (data$V17 - data$V16 + 1)/1000000
data$gene_length <- data$V17 - data$V16 + 1

# 统计每个基因分别与多少个eccDNA有交集
circle_num_per_gene <- as.data.frame(table(data$V18))
colnames(circle_num_per_gene) <- c("gene_name", "ecc_number")
head(circle_num_per_gene)

# 将基因长度新增到circle_num_per_gene的最后一列
gene_len_MB <- data[,c("V18", "gene_length_MB")]
gene_len_MB_dup <- gene_len_MB[!duplicated(gene_len_MB),]
colnames(gene_len_MB_dup) <- c("gene_name", "gene_length_MB")
circle_num_per_gene_len <- merge(circle_num_per_gene, gene_len_MB_dup, by = "gene_name", all = FALSE)
circle_num_per_gene_len$ecc_number_per_MB <- circle_num_per_gene_len$ecc_number/circle_num_per_gene_len$gene_length_MB

# 如果希望将X轴改为log(gene length/MB), Y轴改为log(eccDNA number), 需要增加这一段代码
circle_num_per_gene_len$log_ecc_number <- log(circle_num_per_gene_len$ecc_number)
circle_num_per_gene_len$log_ecc_number_per_MB <- log(circle_num_per_gene_len$ecc_number_per_MB)
circle_num_per_gene_len$log_gene_length_MB <- log(circle_num_per_gene_len$gene_length_MB)
head(circle_num_per_gene_len)

# 读入Normalized后的count矩阵
norm_counts <- read.csv("/Volumes/Denmark/Mouse_Atlas/Bioinformatics/RNAseq/RNA_02/normalized_counts.txt", sep = "\t", header = T)
head(norm_counts)

## FS217 (s115)
## FS138 (s18)
## FS169 (s52)
single_norm_counts <- subset(norm_counts, select = c(X, s115,s18,s52)) 
single_norm_counts$Ctx_A <- single_norm_counts$s115 + single_norm_counts$s18 + single_norm_counts$s52
single_norm_counts <- single_norm_counts[,c("X","Ctx_A")]

colnames(single_norm_counts) <- c("gene_name", "norm_count")
head(single_norm_counts)

# 将每个基因的normalized_count值添加到circle_num_per_gene_len数据框的最后一列
circle_num_per_gene_len_NormCount <- merge(circle_num_per_gene_len, single_norm_counts, by = "gene_name", all = FALSE)
head(circle_num_per_gene_len_NormCount)

#简化变量名称
input <- circle_num_per_gene_len_NormCount

# 对横坐标进行gene_length_MB归一化
input$norm_count_per_MB <- input$norm_count/input$gene_length_MB

# 对横坐标取log()，以e为底
input$log_norm_count_per_MB <-  log(input$norm_count_per_MB)

# 从input数据框中去掉所有log_norm_count_per_MB列是-Inf的行，以免这些点干扰相关系数R
input_rm_Inf <- input[which(input$log_norm_count_per_MB != "-Inf"),]
df <- input_rm_Inf
df <- subset(df, select =c(gene_name, log_ecc_number_per_MB, log_norm_count_per_MB))
head(df)

# tau是分位数回归的分位点
tau_x <- c( 0.1, 0.9) 
quantreg_fit <- rq(log_ecc_number_per_MB ~ bs(log_norm_count_per_MB, df=5), data = df, tau = tau_x)
predict_quantreg <- predict.rq(quantreg_fit, df) %>% data.frame
names(predict_quantreg) <- paste0("quantile_", tau_x)
df <- data.frame(df, predict_quantreg)
log_norm_count_per_MB <- subset(df, select =c(log_norm_count_per_MB))
variable_x <- data.frame(log_norm_count_per_MB)
predict_quantreg_plot <- predict.rq(quantreg_fit, variable_x) %>% data.frame
df_quantreg_plot <- data.frame(variable_x$log_norm_count_per_MB, predict_quantreg_plot)
names(df_quantreg_plot) <- c("log_norm_count_per_MB", paste0("quantile_",tau_x))

# Plot
ggplot(data = df, aes(x=log_norm_count_per_MB, y=log_ecc_number_per_MB)) +
  geom_point(alpha=0.8, color="gray", size = 0.5) +
  theme_classic() +
  stat_density_2d(alpha = 0.2, geom = "polygon", contour = TRUE, aes(fill = after_stat(level)), colour = "white", bins = 10) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  labs(title = "12-month cortex", 
       x = "Log(normalized counts per gene/MB)", 
       y = "Log(eccDNA number per gene/MB)") +
  geom_line(data=df_quantreg_plot, aes(x=log_norm_count_per_MB, y=quantile_0.1), lwd=0.8, linetype = "dashed", color="#FF9900")+
  geom_line(data=df_quantreg_plot, aes(x=log_norm_count_per_MB, y=quantile_0.9), lwd=0.8, linetype = "dashed", color="#FF9900")
ggsave("13_spline_regression_cortex_adult.pdf", width = 5, height = 5)

# 提取小于或等于10%分位数的基因(label="below")
df$label <- ifelse(df$log_ecc_number_per_MB <= df$quantile_0.1, "below",
                   ifelse(df$log_ecc_number_per_MB >= df$quantile_0.9, "above", "middle"))
write.csv(df, file = "10%_quantile_genelist_cortex_adult.csv", row.names = FALSE)

# 读入注释完蛋白编码基因后的eccDNA bed文件
rawdata <- read.csv("/Volumes/Denmark/Mouse_Atlas/Bioinformatics/RNAseq/RNA_03/merged_anno.bed", sep = "\t", header = F)
head(rawdata)

#------------------------------------------------22-month Ctx (n=4)------------------------------------------------

# ImporteccDNA bed file annotated by protein-coding genes
rawdata <- read.csv("/Volumes/Denmark/Mouse_Atlas/Bioinformatics/RNAseq/Liver_01/anno.bed", sep = "\t", header = F)
head(rawdata)

## FS253 (s161)
## FS268 (s176)
## FS275 (s185)
## FS280 (s213)

data <- rawdata[which(rawdata$V12 == "FS253"|rawdata$V12 == "FS268"|rawdata$V12 == "FS275"|rawdata$V12 == "FS280"),]

#计算每个蛋白编码基因的长度（单位：MB)
data$gene_length_MB <- (data$V17 - data$V16 + 1)/1000000
data$gene_length <- data$V17 - data$V16 + 1

# 统计每个基因分别与多少个eccDNA有交集
circle_num_per_gene <- as.data.frame(table(data$V18))
colnames(circle_num_per_gene) <- c("gene_name", "ecc_number")
head(circle_num_per_gene)

# 将基因长度新增到circle_num_per_gene的最后一列
gene_len_MB <- data[,c("V18", "gene_length_MB")]
gene_len_MB_dup <- gene_len_MB[!duplicated(gene_len_MB),]
colnames(gene_len_MB_dup) <- c("gene_name", "gene_length_MB")
circle_num_per_gene_len <- merge(circle_num_per_gene, gene_len_MB_dup, by = "gene_name", all = FALSE)
circle_num_per_gene_len$ecc_number_per_MB <- circle_num_per_gene_len$ecc_number/circle_num_per_gene_len$gene_length_MB

# 如果希望将X轴改为log(gene length/MB), Y轴改为log(eccDNA number), 需要增加这一段代码
circle_num_per_gene_len$log_ecc_number <- log(circle_num_per_gene_len$ecc_number)
circle_num_per_gene_len$log_ecc_number_per_MB <- log(circle_num_per_gene_len$ecc_number_per_MB)
circle_num_per_gene_len$log_gene_length_MB <- log(circle_num_per_gene_len$gene_length_MB)
head(circle_num_per_gene_len)

# 读入Normalized后的count矩阵
norm_counts <- read.csv("/Volumes/Denmark/Mouse_Atlas/Bioinformatics/RNAseq/RNA_03/normalized_counts.txt", sep = "\t", header = T)
head(norm_counts)

## FS253 (s161)
## FS268 (s176)
## FS275 (s185)
## FS280 (s213)
single_norm_counts <- subset(norm_counts, select = c(X, s161,s176,s185,s213)) 
single_norm_counts$Ctx_O <- single_norm_counts$s161 + single_norm_counts$s176 + single_norm_counts$s185 + single_norm_counts$s213 
single_norm_counts <- single_norm_counts[,c("X","Ctx_O")]

colnames(single_norm_counts) <- c("gene_name", "norm_count")
head(single_norm_counts)

# 将每个基因的normalized_count值添加到circle_num_per_gene_len数据框的最后一列
circle_num_per_gene_len_NormCount <- merge(circle_num_per_gene_len, single_norm_counts, by = "gene_name", all = FALSE)
head(circle_num_per_gene_len_NormCount)

#简化变量名称
input <- circle_num_per_gene_len_NormCount

# 对横坐标进行gene_length_MB归一化
input$norm_count_per_MB <- input$norm_count/input$gene_length_MB

# 对横坐标取log()
input$log_norm_count_per_MB <-  log(input$norm_count_per_MB)

# 从input数据框中去掉所有log_norm_count_per_MB列是-Inf的行，以免这些点干扰相关系数R
input_rm_Inf <- input[which(input$log_norm_count_per_MB != "-Inf"),]

df <- input_rm_Inf
df <- subset(df, select =c(gene_name, log_ecc_number_per_MB, log_norm_count_per_MB))
head(df)

# tau是分位数回归的分位点
tau_x <- c( 0.1, 0.9) 
quantreg_fit <- rq(log_ecc_number_per_MB ~ bs(log_norm_count_per_MB, df=5), data = df, tau = tau_x)
predict_quantreg <- predict.rq(quantreg_fit, df) %>% data.frame
names(predict_quantreg) <- paste0("quantile_", tau_x)
df <- data.frame(df, predict_quantreg)
log_norm_count_per_MB <- subset(df, select =c(log_norm_count_per_MB))
variable_x <- data.frame(log_norm_count_per_MB)
predict_quantreg_plot <- predict.rq(quantreg_fit, variable_x) %>% data.frame
df_quantreg_plot <- data.frame(variable_x$log_norm_count_per_MB, predict_quantreg_plot)
names(df_quantreg_plot) <- c("log_norm_count_per_MB", paste0("quantile_",tau_x))

# Plot
ggplot(data = df, aes(x=log_norm_count_per_MB, y=log_ecc_number_per_MB)) +
  geom_point(alpha=0.8, color="gray", size = 0.5) +
  theme_classic() +
  stat_density_2d(alpha = 0.2, geom = "polygon", contour = TRUE, aes(fill = after_stat(level)), colour = "white", bins = 10) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  labs(title = "22-month cortex", 
       x = "Log(normalized counts per gene/MB)", 
       y = "Log(eccDNA number per gene/MB)") +
  geom_line(data=df_quantreg_plot, aes(x=log_norm_count_per_MB, y=quantile_0.1), lwd=0.8, linetype = "dashed", color="#006699")+
  geom_line(data=df_quantreg_plot, aes(x=log_norm_count_per_MB, y=quantile_0.9), lwd=0.8, linetype = "dashed", color="#006699")
ggsave("13_spline_regression_cortex_old.pdf", width = 5, height = 5)

# 提取小于或等于10%分位数的基因(label="below")
df$label <- ifelse(df$log_ecc_number_per_MB <= df$quantile_0.1, "below",
                   ifelse(df$log_ecc_number_per_MB >= df$quantile_0.9, "above", "middle"))
write.csv(df, file = "10%_quantile_genelist_cortex_old.csv", row.names = FALSE)

#------------------------------------------------3-month hippocampus (n=4)------------------------------------------------

# ImporteccDNA bed file annotated by protein-coding genes
rawdata <- read.csv("/Volumes/Denmark/Mouse_Atlas/Bioinformatics/RNAseq/Liver_01/anno.bed", sep = "\t", header = F)
head(rawdata)

## FS221 (s120)
## FS155 (s35)
## FS194 (s82)
## FS203 (s94)
data <- rawdata[which(rawdata$V12 == "FS221"|rawdata$V12 == "FS155"|rawdata$V12 == "FS194"|rawdata$V12 == "FS203"),]

#计算每个蛋白编码基因的长度（单位：MB)
data$gene_length_MB <- (data$V17 - data$V16 + 1)/1000000
data$gene_length <- data$V17 - data$V16 + 1

# 统计每个基因分别与多少个eccDNA有交集
circle_num_per_gene <- as.data.frame(table(data$V18))
colnames(circle_num_per_gene) <- c("gene_name", "ecc_number")
head(circle_num_per_gene)

# 将基因长度新增到circle_num_per_gene的最后一列
gene_len_MB <- data[,c("V18", "gene_length_MB")]
gene_len_MB_dup <- gene_len_MB[!duplicated(gene_len_MB),]
colnames(gene_len_MB_dup) <- c("gene_name", "gene_length_MB")
circle_num_per_gene_len <- merge(circle_num_per_gene, gene_len_MB_dup, by = "gene_name", all = FALSE)
circle_num_per_gene_len$ecc_number_per_MB <- circle_num_per_gene_len$ecc_number/circle_num_per_gene_len$gene_length_MB

# 如果希望将X轴改为log(gene length/MB), Y轴改为log(eccDNA number), 需要增加这一段代码
circle_num_per_gene_len$log_ecc_number <- log(circle_num_per_gene_len$ecc_number)
circle_num_per_gene_len$log_ecc_number_per_MB <- log(circle_num_per_gene_len$ecc_number_per_MB)
circle_num_per_gene_len$log_gene_length_MB <- log(circle_num_per_gene_len$gene_length_MB)
head(circle_num_per_gene_len)

# 读入Normalized后的count矩阵
norm_counts <- read.csv("/Volumes/Denmark/Mouse_Atlas/Bioinformatics/RNAseq/RNA_02/normalized_counts.txt", sep = "\t", header = T)
head(norm_counts)

## FS221 (s120)
## FS155 (s35)
## FS194 (s82)
## FS203 (s94)
single_norm_counts <- subset(norm_counts, select = c(X, s120,s35,s82,s94)) 
single_norm_counts$Hipp_Y <- single_norm_counts$s120 + single_norm_counts$s35 + single_norm_counts$s82 + single_norm_counts$s94
single_norm_counts <- single_norm_counts[,c("X","Hipp_Y")]

colnames(single_norm_counts) <- c("gene_name", "norm_count")
head(single_norm_counts)

# 将每个基因的normalized_count值添加到circle_num_per_gene_len数据框的最后一列
circle_num_per_gene_len_NormCount <- merge(circle_num_per_gene_len, single_norm_counts, by = "gene_name", all = FALSE)
head(circle_num_per_gene_len_NormCount)

#简化变量名称
input <- circle_num_per_gene_len_NormCount

# 对横坐标进行gene_length_MB归一化
input$norm_count_per_MB <- input$norm_count/input$gene_length_MB

# 对横坐标取log()
input$log_norm_count_per_MB <-  log(input$norm_count_per_MB)

# 从input数据框中去掉所有log_norm_count_per_MB列是-Inf的行，以免这些点干扰相关系数R
input_rm_Inf <- input[which(input$log_norm_count_per_MB != "-Inf"),]
df <- input_rm_Inf
df <- subset(df, select =c(gene_name, log_ecc_number_per_MB, log_norm_count_per_MB))
head(df)

# tau是分位数回归的分位点
tau_x <- c( 0.1, 0.9) 
quantreg_fit <- rq(log_ecc_number_per_MB ~ bs(log_norm_count_per_MB, df=5), data = df, tau = tau_x)
predict_quantreg <- predict.rq(quantreg_fit, df) %>% data.frame
names(predict_quantreg) <- paste0("quantile_", tau_x)
df <- data.frame(df, predict_quantreg)
log_norm_count_per_MB <- subset(df, select =c(log_norm_count_per_MB))
variable_x <- data.frame(log_norm_count_per_MB)
predict_quantreg_plot <- predict.rq(quantreg_fit, variable_x) %>% data.frame
df_quantreg_plot <- data.frame(variable_x$log_norm_count_per_MB, predict_quantreg_plot)
names(df_quantreg_plot) <- c("log_norm_count_per_MB", paste0("quantile_",tau_x))

# Plot
ggplot(data = df, aes(x=log_norm_count_per_MB, y=log_ecc_number_per_MB)) +
  geom_point(alpha=0.8, color="gray", size = 0.5) +
  theme_classic() +
  stat_density_2d(alpha = 0.2, geom = "polygon", contour = TRUE, aes(fill = after_stat(level)), colour = "white", bins = 10) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  labs(title = "3-month hippocampus", 
       x = "Log(normalized counts per gene/MB)", 
       y = "Log(eccDNA number per gene/MB)") +
  geom_line(data=df_quantreg_plot, aes(x=log_norm_count_per_MB, y=quantile_0.1), lwd=0.8, linetype = "dashed", color="#339900")+
  geom_line(data=df_quantreg_plot, aes(x=log_norm_count_per_MB, y=quantile_0.9), lwd=0.8, linetype = "dashed", color="#339900")
ggsave("14_spline_regression_hippocampus_young.pdf", width = 5, height = 5)

# 提取小于或等于10%分位数的基因(label="below")
df$label <- ifelse(df$log_ecc_number_per_MB <= df$quantile_0.1, "below",
                   ifelse(df$log_ecc_number_per_MB >= df$quantile_0.9, "above", "middle"))
write.csv(df, file = "10%_quantile_genelist_hippocampus_young.csv", row.names = FALSE)


#------------------------------------------------12 month Hipp的样本------------------------------------------------

# ImporteccDNA bed file annotated by protein-coding genes
rawdata <- read.csv("/Volumes/Denmark/Mouse_Atlas/Bioinformatics/RNAseq/Liver_01/anno.bed", sep = "\t", header = F)
head(rawdata)

## FS216 (s114)
## FS195 (s83)
## FS196 (s84)
## FS128 (s8)
data <- rawdata[which(rawdata$V12 == "FS216"|rawdata$V12 == "FS195"|rawdata$V12 == "FS196"|rawdata$V12 == "FS128"),]

#计算每个蛋白编码基因的长度（单位：MB)
data$gene_length_MB <- (data$V17 - data$V16 + 1)/1000000
data$gene_length <- data$V17 - data$V16 + 1

# 统计每个基因分别与多少个eccDNA有交集
circle_num_per_gene <- as.data.frame(table(data$V18))
colnames(circle_num_per_gene) <- c("gene_name", "ecc_number")
head(circle_num_per_gene)

# 将基因长度新增到circle_num_per_gene的最后一列
gene_len_MB <- data[,c("V18", "gene_length_MB")]
gene_len_MB_dup <- gene_len_MB[!duplicated(gene_len_MB),]
colnames(gene_len_MB_dup) <- c("gene_name", "gene_length_MB")
circle_num_per_gene_len <- merge(circle_num_per_gene, gene_len_MB_dup, by = "gene_name", all = FALSE)
circle_num_per_gene_len$ecc_number_per_MB <- circle_num_per_gene_len$ecc_number/circle_num_per_gene_len$gene_length_MB

# 如果希望将X轴改为log(gene length/MB), Y轴改为log(eccDNA number), 需要增加这一段代码
circle_num_per_gene_len$log_ecc_number <- log(circle_num_per_gene_len$ecc_number)
circle_num_per_gene_len$log_ecc_number_per_MB <- log(circle_num_per_gene_len$ecc_number_per_MB)
circle_num_per_gene_len$log_gene_length_MB <- log(circle_num_per_gene_len$gene_length_MB)
head(circle_num_per_gene_len)

# 读入Normalized后的count矩阵
norm_counts <- read.csv("/Volumes/Denmark/Mouse_Atlas/Bioinformatics/RNAseq/RNA_02/normalized_counts.txt", sep = "\t", header = T)
head(norm_counts)

## FS216 (s114)
## FS195 (s83)
## FS196 (s84)
## FS128 (s8)
single_norm_counts <- subset(norm_counts, select = c(X, s114,s83,s84,s8)) 
single_norm_counts$Hipp_A <- single_norm_counts$s114 + single_norm_counts$s83 + single_norm_counts$s84 + single_norm_counts$s8
single_norm_counts <- single_norm_counts[,c("X","Hipp_A")]

colnames(single_norm_counts) <- c("gene_name", "norm_count")
head(single_norm_counts)

# 将每个基因的normalized_count值添加到circle_num_per_gene_len数据框的最后一列
circle_num_per_gene_len_NormCount <- merge(circle_num_per_gene_len, single_norm_counts, by = "gene_name", all = FALSE)
head(circle_num_per_gene_len_NormCount)

#简化变量名称
input <- circle_num_per_gene_len_NormCount

# 对横、纵坐标分别进行gene_length_MB归一化
input$norm_count_per_MB <- input$norm_count/input$gene_length_MB

# 对横、纵坐标分别取log()
input$log_norm_count_per_MB <-  log(input$norm_count_per_MB)

# 从input数据框中去掉所有log_norm_count_per_MB列是-Inf的行，以免这些点干扰相关系数R
input_rm_Inf <- input[which(input$log_norm_count_per_MB != "-Inf"),]
df <- input_rm_Inf
df <- subset(df, select =c(gene_name, log_ecc_number_per_MB, log_norm_count_per_MB))
head(df)

# Conducted cubic spline quantile regression
tau_x <- c(0.1, 0.9) 
quantreg_fit <- rq(log_ecc_number_per_MB ~ bs(log_norm_count_per_MB, df=5), data = df, tau = tau_x)
predict_quantreg <- predict.rq(quantreg_fit, df) %>% data.frame
names(predict_quantreg) <- paste0("quantile_", tau_x)
df <- data.frame(df, predict_quantreg)
log_norm_count_per_MB <- subset(df, select =c(log_norm_count_per_MB))
variable_x<-data.frame(log_norm_count_per_MB)
predict_quantreg_plot <- predict.rq(quantreg_fit, variable_x) %>% data.frame
df_quantreg_plot <- data.frame(variable_x$log_norm_count_per_MB, predict_quantreg_plot)
names(df_quantreg_plot) <- c("log_norm_count_per_MB", paste0("quantile_",tau_x))

# Plot Adult
ggplot(data = df, aes(x=log_norm_count_per_MB, y=log_ecc_number_per_MB)) +
  geom_point(alpha=0.8, color="gray", size = 0.5) +
  theme_classic() +
  stat_density_2d(alpha = 0.2, geom = "polygon", contour = TRUE, aes(fill = after_stat(level)), colour = "white", bins = 10) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  labs(title = "12-month hippocampus", 
       x = "Log(normalized counts per gene/MB)", 
       y = "Log(eccDNA number per gene/MB)") +
  geom_line(data=df_quantreg_plot, aes(x=log_norm_count_per_MB, y=quantile_0.1), lwd=0.8, linetype = "dashed", color="#FF9900")+
  geom_line(data=df_quantreg_plot, aes(x=log_norm_count_per_MB, y=quantile_0.9), lwd=0.8, linetype = "dashed", color="#FF9900")
ggsave("14_spline_regression_hippocampus_adult.pdf", width = 5, height = 5)

# Extact the gene fall below 10% quantiles
head(df)

# 提取小于或等于10%分位数的基因(label="below")
df$label <- ifelse(df$log_ecc_number_per_MB <= df$quantile_0.1, "below",
                   ifelse(df$log_ecc_number_per_MB >= df$quantile_0.9, "above", "middle"))
write.csv(df, file = "10%_quantile_genelist_hippocampus_adult.csv", row.names = FALSE)



#------------------------------------------------22-month hippocampus (n=4)------------------------------------------------

# ImporteccDNA bed file annotated by protein-coding genes
rawdata <- read.csv("/Volumes/Denmark/Mouse_Atlas/Bioinformatics/RNAseq/Liver_01/anno.bed", sep = "\t", header = F)
head(rawdata)

## FS257 (s165)
## FS264 (s172)
## FS300 (s208)
## FS295 (s209)

data <- rawdata[which(rawdata$V12 == "FS257"|rawdata$V12 == "FS264"|rawdata$V12 == "FS300"|rawdata$V12 == "FS295"),]

#计算每个蛋白编码基因的长度（单位：MB)
data$gene_length_MB <- (data$V17 - data$V16 + 1)/1000000
data$gene_length <- data$V17 - data$V16 + 1

# 统计每个基因分别与多少个eccDNA有交集
circle_num_per_gene <- as.data.frame(table(data$V18))
colnames(circle_num_per_gene) <- c("gene_name", "ecc_number")
head(circle_num_per_gene)

# 将基因长度新增到circle_num_per_gene的最后一列
gene_len_MB <- data[,c("V18", "gene_length_MB")]
gene_len_MB_dup <- gene_len_MB[!duplicated(gene_len_MB),]
colnames(gene_len_MB_dup) <- c("gene_name", "gene_length_MB")
circle_num_per_gene_len <- merge(circle_num_per_gene, gene_len_MB_dup, by = "gene_name", all = FALSE)
circle_num_per_gene_len$ecc_number_per_MB <- circle_num_per_gene_len$ecc_number/circle_num_per_gene_len$gene_length_MB

# 如果希望将X轴改为log(gene length/MB), Y轴改为log(eccDNA number), 需要增加这一段代码
circle_num_per_gene_len$log_ecc_number <- log(circle_num_per_gene_len$ecc_number)
circle_num_per_gene_len$log_ecc_number_per_MB <- log(circle_num_per_gene_len$ecc_number_per_MB)
circle_num_per_gene_len$log_gene_length_MB <- log(circle_num_per_gene_len$gene_length_MB)
head(circle_num_per_gene_len)

# 读入Normalized后的count矩阵
norm_counts <- read.csv("/Volumes/Denmark/Mouse_Atlas/Bioinformatics/RNAseq/RNA_03/normalized_counts.txt", sep = "\t", header = T)
head(norm_counts)

## FS257 (s165)
## FS264 (s172)
## FS300 (s208)
## FS295 (s209)
single_norm_counts <- subset(norm_counts, select = c(X, s165,s172,s208,s209)) 
single_norm_counts$Hipp_O <- single_norm_counts$s165 + single_norm_counts$s172 + single_norm_counts$s208 + single_norm_counts$s209
single_norm_counts <- single_norm_counts[,c("X","Hipp_O")]

colnames(single_norm_counts) <- c("gene_name", "norm_count")
head(single_norm_counts)

# 将每个基因的normalized_count值添加到circle_num_per_gene_len数据框的最后一列
circle_num_per_gene_len_NormCount <- merge(circle_num_per_gene_len, single_norm_counts, by = "gene_name", all = FALSE)
head(circle_num_per_gene_len_NormCount)

#简化变量名称
input <- circle_num_per_gene_len_NormCount

# 对横坐标进行gene_length_MB归一化
input$norm_count_per_MB <- input$norm_count/input$gene_length_MB

# 对横坐标取log()
input$log_norm_count_per_MB <-  log(input$norm_count_per_MB)

# 从input数据框中去掉所有log_norm_count_per_MB列是-Inf的行，以免这些点干扰相关系数R
input_rm_Inf <- input[which(input$log_norm_count_per_MB != "-Inf"),]
df <- input_rm_Inf
df <- subset(df, select =c(gene_name, log_ecc_number_per_MB, log_norm_count_per_MB))
head(df)

# tau是分位数回归的分位点
tau_x <- c( 0.1, 0.9) 
quantreg_fit <- rq(log_ecc_number_per_MB ~ bs(log_norm_count_per_MB, df=5), data = df, tau = tau_x)
predict_quantreg <- predict.rq(quantreg_fit, df) %>% data.frame
names(predict_quantreg) <- paste0("quantile_", tau_x)
df <- data.frame(df, predict_quantreg)
log_norm_count_per_MB <- subset(df, select =c(log_norm_count_per_MB))
variable_x <- data.frame(log_norm_count_per_MB)
predict_quantreg_plot <- predict.rq(quantreg_fit, variable_x) %>% data.frame
df_quantreg_plot <- data.frame(variable_x$log_norm_count_per_MB, predict_quantreg_plot)
names(df_quantreg_plot) <- c("log_norm_count_per_MB", paste0("quantile_",tau_x))

# Plot
ggplot(data = df, aes(x=log_norm_count_per_MB, y=log_ecc_number_per_MB)) +
  geom_point(alpha=0.8, color="gray", size = 0.5) +
  theme_classic() +
  stat_density_2d(alpha = 0.2, geom = "polygon", contour = TRUE, aes(fill = after_stat(level)), colour = "white", bins = 10) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  labs(title = "22-month hippocampus", 
       x = "Log(normalized counts per gene/MB)", 
       y = "Log(eccDNA number per gene/MB)") +
  geom_line(data=df_quantreg_plot, aes(x=log_norm_count_per_MB, y=quantile_0.1), lwd=0.8, linetype = "dashed", color="#006699")+
  geom_line(data=df_quantreg_plot, aes(x=log_norm_count_per_MB, y=quantile_0.9), lwd=0.8, linetype = "dashed", color="#006699")
ggsave("14_spline_regression_hippocampus_old.pdf", width = 5, height = 5)

# 提取小于或等于10%分位数的基因(label="below")
df$label <- ifelse(df$log_ecc_number_per_MB <= df$quantile_0.1, "below",
                   ifelse(df$log_ecc_number_per_MB >= df$quantile_0.9, "above", "middle"))
write.csv(df, file = "10%_quantile_genelist_hippocampus_old.csv", row.names = FALSE)
