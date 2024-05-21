# 1. Calculate the fraction of eccDNA below 2000 bp in length (16.8 %)
# 2. Calculate the average length of eccDNA (4786.5 bp)
# 3. Calculate the length range of eccDNA (varied from 8 bp to 258322 bp)

# Load required Rdata
# 138 high mapping ratio samples(>=70%)
hh <- read.csv("/Volumes/Denmark/Mouse_Atlas/Bioinformatics/Figures/Fig1/03_circle_number/inputfile/highMappingRate_clean_merged_145_medium_conf_labeled_circle.bed", sep = "\t", header = F)

# 75 low mapping ratio samples(<70%)
# hh <- read.csv("/Users/lx962456367/Desktop/test/merged_75lowQC_medium_conf_labeled_circle.bed", sep = "\t", header = F)

# 只保留已知的21条常染色体上的环merged_75lowQC_medium_conf_labeled_circle.bed
hh1 <- hh[hh$V1 %in% c("chr1","chr2","chr3","chr4","chr5",
                       "chr6","chr7","chr8","chr9","chr10",
                       "chr11","chr12","chr13","chr14","chr15",
                       "chr16","chr17","chr18","chr19","chrX","chrY"),]

#-------------- 1.

# 计算eccDNA的长度
hh1$ecc_len <- hh1$V3 - hh1$V2

# 计算大于等于2000的eccDNA数目
count_greater_2000 <- sum(hh1$ecc_len >= 2000)

# 计算小于2000的eccDNA数目
count_less_2000 <- sum(hh1$ecc_len < 2000)

# 获得总eccDNA数目
total_rows <- nrow(hh1)

# 计算占比
percent_greater_2000 <- (count_greater_2000 / total_rows) * 100
percent_less_2000 <- (count_less_2000 / total_rows) * 100

# 打印结果
cat("Number of rows with ecc_len >= 2000:", count_greater_2000, "\n")
# Number of rows with ecc_len >= 2000: 472424 

cat("Number of rows with ecc_len < 2000:", count_less_2000, "\n")
# Number of rows with ecc_len < 2000: 95539 

cat("Percentage of rows with ecc_len >= 2000:", percent_greater_2000, "%\n")
# Percentage of rows with ecc_len >= 2000: 83.17866 %

cat("Percentage of rows with ecc_len < 2000:", percent_less_2000, "%\n")
# Percentage of rows with ecc_len < 2000: 16.82134 %

#-------------- 2.

# 计算 ecc_len 列的最大值
max_ecc_len <- max(hh1$ecc_len)

# 计算 ecc_len 列的最小值
min_ecc_len <- min(hh1$ecc_len)

# 打印结果
cat("Maximum value in ecc_len:", max_ecc_len, "\n")
cat("Minimum value in ecc_len:", min_ecc_len, "\n")

#-------------- 3.

# 计算 ecc_len 列的平均值
mean_ecc_len <- mean(hh1$ecc_len, na.rm = TRUE)

# 计算 ecc_len 列的中位数
median_ecc_len <- median(hh1$ecc_len, na.rm = TRUE)

# 打印结果
cat("Average value of ecc_len:", mean_ecc_len, "\n")
cat("Median value of ecc_len:", median_ecc_len, "\n")

