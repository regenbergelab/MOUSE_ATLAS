# Mouse Atlas - Fig1B - 25_circle_size(Pie Chart)

rm(list = ls())
options(stringsAsFactors = FALSE)

# Set working folder
mainDir <- "/Volumes/Denmark/Mouse_Atlas/Bioinformatics/Figures/Fig1"
subDir <- "25_circle_size"
# create a directory if it doesn't exist
ifelse(!dir.exists(file.path(mainDir, subDir)), dir.create(file.path(mainDir, subDir)), FALSE)
setwd(file.path(mainDir, subDir))

# Load required libraries

#-----test data

# Load necessary libraries
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

# Example data: proportions and labels
data <- data.frame(
  category = c("Category A", "Category B", "Category C", "Category D"),
  proportion = c(0.3, 0.25, 0.2, 0.25)
)

# Basic pie chart using ggplot2
pie_chart <- ggplot(data, aes(x = "", y = proportion, fill = category)) +
  geom_bar(stat = "identity", width = 1, color = "white") +  # Create pie chart as stacked bar chart
  coord_polar("y", start = 0) +  # Convert to polar coordinates
  labs(title = "Pie Chart Example", fill = "Category") +  # Add title and legend label
  theme_void()  # Remove background and grid lines

# Print the pie chart
print(pie_chart)

###----my data
# Load required Rdata
# 138 high mapping ratio samples(>=70%), 567963 eccDNA
hh <- read.csv("/Volumes/Denmark/Mouse_Atlas/Bioinformatics/Figures/Fig1/03_circle_number/inputfile/highMappingRate_clean_merged_145_medium_conf_labeled_circle.bed", sep = "\t", header = F)

# 只保留已知的21条常染色体上的环merged_75lowQC_medium_conf_labeled_circle.bed
hh1 <- hh[hh$V1 %in% c("chr1","chr2","chr3","chr4","chr5",
                       "chr6","chr7","chr8","chr9","chr10",
                       "chr11","chr12","chr13","chr14","chr15",
                       "chr16","chr17","chr18","chr19","chrX","chrY"),]

# Calculate the length of eccDNA
hh1$ecc_Len <- hh1$V3 - hh1$V2

# 导入必要的库
library(dplyr)
library(ggplot2)

df <- hh1

# 计算各个范围内的数量和比例
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

# 输出汇总结果
print(summary_df)

# 制作饼图展示比例
pie_chart <- ggplot(summary_df, aes(x = "", y = proportion, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  labs(title = "Ecc_Len Distribution", fill = "Category") +
  theme_void()

# 显示饼图
print(pie_chart)

# 定义颜色向量
my_colors <- c("#fde8c6","#fbc979", "#e09148") #"#c4c9d8" 


# 制作饼图展示比例，并设置颜色和边框
pie_chart <- ggplot(summary_df, aes(x = "", y = proportion, fill = category)) +
  geom_bar(stat = "identity", width = 1, color = "black") +  # 设置黑色边框
  coord_polar("y", start = 0) +
  labs(title = "Ecc_Len Distribution", fill = "Category") +
  scale_fill_manual(values = my_colors) +  # 设置颜色
  theme_void()

# 显示饼图
print(pie_chart)
ggsave(pie_chart, file = "pie_chart.pdf") # Fig1B-3



# 导入必要的库
library(ggplot2)

# 数据框
df <- data.frame(ecc_Len = hh1$ecc_Len)

# 设置分段和颜色映射
cuts <- c(0, 2000, 10000)
colors <- c("#fde8c6", "#fbc979")

# 使用ggplot绘制直方图，并添加描边
histogram1 <- ggplot(df, aes(x = ecc_Len)) +
  geom_histogram(aes(y = ..density..), bins = 1000, fill = "skyblue", color = "#fbc979", alpha = 0.7) +  
  # 设置直方图的填充颜色为"skyblue"，描边颜色为"black"，透明度为0.7
  labs(title = "Histogram of ecc_Len", x = "ecc_Len", y = "Density") +
  theme_classic()+
  xlim(0, 10000)

# 显示直方图
print(histogram1)

# 使用ggplot绘制直方图，并添加描边
histogram2 <- ggplot(df, aes(x = ecc_Len)) +
  geom_histogram(aes(y = ..density..), bins = 1000, fill = "skyblue", color = "#fde8c6", alpha = 0.7) +  
  # 设置直方图的填充颜色为"skyblue"，描边颜色为"black"，透明度为0.7
  labs(title = "Histogram of ecc_Len", x = "ecc_Len", y = "Density") +
  theme_classic()+
  xlim(0, 10000)+
  geom_vline(xintercept = 2000, color = "black", linetype = "solid")

# 显示直方图
print(histogram2)

ggsave(histogram1, file = "histogram1.pdf") # Fig1B-1
ggsave(histogram2, file = "histogram2.pdf") # Fig1B-2

# # Fig1C

# 使用ggplot绘制直方图，并添加描边
histogram3 <- ggplot(df, aes(x = ecc_Len)) +
  geom_histogram(aes(y = ..density..), bins = 1000, color = "#fde8c6", alpha = 0.7) +  
  # 设置直方图的填充颜色为"skyblue"，描边颜色为"black"，透明度为0.7
  labs(title = "Histogram of ecc_Len", x = "ecc_Len", y = "Density") +
  theme_classic()+
  xlim(100, 2000)+
  geom_vline(xintercept = c(200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800), color = "gray", linetype = "dashed")

# 显示直方图
print(histogram3)
ggsave(histogram3, file = "histogram3.pdf") # Fig1C
