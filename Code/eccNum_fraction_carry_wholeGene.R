# Calculate the number and fraction of the eccDNA that was large enough to carry whole genes, not only in VAT but also in other tissues
# (Liver: 260, 0.76%; Hippocampus: 255, 0.92%; Cortex: 412, 0.76%; Muscle: 162, 0.59%; SAT: 420, 1.28%; VAT: 196, 1.05%). 

# 读入注释完蛋白编码基因后的eccDNA bed文件
rawdata <- read.csv("/Volumes/Denmark/Mouse_Atlas/Bioinformatics/RNAseq/Liver_01/anno.bed", sep = "\t", header = F)
head(rawdata)

# 使用 ifelse() 函数和逻辑条件来创建新列
rawdata$carry_whole_gene <- ifelse(rawdata$V16 >= rawdata$V2 & rawdata$V17 <= rawdata$V3, "Yes", "No")

# 查看修改后的数据框
head(rawdata)

# 计算不同组织类型的 Yes 和 No 的数量
rawdata$V13 <- as.factor(rawdata$V13)
# 创建交叉表
cross_table <- table(rawdata$V13, rawdata$carry_whole_gene)
cross_table

# 计算每一行的总数
row_totals <- rowSums(cross_table)
row_totals

# 创建一个数据框来保存每一行的总数以及 Yes 和 No 的比例
summary_table <- data.frame(
  Organization = names(row_totals),   # 组织类型的名称
  Total = row_totals,                 # 每一行的总数
  Yes_Number = cross_table[,"Yes"],   # Yes 的数量
  No_Number = cross_table[,"No"],     # No 的数量
  Yes_Percent = (cross_table[, "Yes"] / row_totals) * 100,  # Yes 的比例
  No_Percent = (cross_table[, "No"] / row_totals) * 100     # No 的比例
)

# 打印总结数据框
print(summary_table)

summary_table
#   Organization Total Yes_Number No_Number Yes_Percent No_Percent
#BD           BD  2619         57      2562   2.1764032   97.82360
#CT           CT 53934        412     53522   0.7638966   99.23610
#HP           HP 27661        255     27406   0.9218756   99.07812
#LV           LV 34115        260     33855   0.7621281   99.23787
#MS           MS 27274        162     27112   0.5939723   99.40603
#PA           PA 64722        424     64298   0.6551095   99.34489
#SF           SF 32764        420     32344   1.2818948   98.71811
#VF           VF 18752        196     18556   1.0452218   98.95478

write.csv(summary_table, file = "summary_table_carry_whole_gene_for_each_tissue.csv")
