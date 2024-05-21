# Calculate the number and fraction of the eccDNA that was large enough to carry whole genes, not only in VAT but also in other tissues
# (Liver: 260, 0.76%; Hippocampus: 255, 0.92%; Cortex: 412, 0.76%; Muscle: 162, 0.59%; SAT: 420, 1.28%; VAT: 196, 1.05%). 

# 读入注释完蛋白编码基因后的eccDNA bed文件
rawdata <- read.csv("/Volumes/Denmark/Mouse_Atlas/Bioinformatics/RNAseq/Liver_01/anno.bed", sep = "\t", header = F)
head(rawdata)

# 使用 ifelse() 函数和逻辑条件来创建新列
rawdata$carry_whole_gene <- ifelse(rawdata$V16 >= rawdata$V2 & rawdata$V17 <= rawdata$V3, "Yes", "No")
rawdata

# 计算不同组织类型的 Yes 和 No 的数量
rawdata$V13 <- as.factor(rawdata$V13)

# 创建交叉表
cross_table <- table(rawdata$V13, rawdata$carry_whole_gene)
cross_table

# 计算每一行的总数
row_totals <- apply(cross_table, 1, sum)
row_totals

# 计算每一行的 Yes 和 No 的比例
percent_table <- prop.table(cross_table)
percent_table

# 你也可以创建一个数据框来保存每一行的总数和比例
summary_table <- data.frame(
  Tissue = names(row_totals),
  Total = row_totals,
  Yes_Number = cross_table[,"Yes"],
  No_Number = cross_table[,"Yes"],
  Yes_Percent = percent_table[,"Yes"] * 100,
  No_Percent = percent_table[,"No"] * 100
)

summary_table
# Tissue Total Yes_Number No_Number Yes_Percent No_Percent
# BD     BD  2619         57        57  0.02176894  0.9784564
# CT     CT 53934        412       412  0.15734740 20.4406491
# HP     HP 27661        255       255  0.09738735 10.4666572
# LV     LV 34115        260       260  0.09929690 12.9296023
# MS     MS 27274        162       162  0.06186961 10.3543754
# PA     PA 64722        424       424  0.16193033 24.5561238
# SF     SF 32764        420       420  0.16040269 12.3525346
# VF     VF 18752        196       196  0.07485459  7.0867435

write.csv(summary_table, file = "summary_table_carry_whole_gene_for_each_tissue.csv")

