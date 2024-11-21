rawdata <- read.csv("anno.bed", sep = "\t", header = F)
rawdata$carry_whole_gene <- ifelse(rawdata$V16 >= rawdata$V2 & rawdata$V17 <= rawdata$V3, "Yes", "No")
rawdata$V13 <- as.factor(rawdata$V13)
cross_table <- table(rawdata$V13, rawdata$carry_whole_gene)
row_totals <- rowSums(cross_table)
summary_table <- data.frame(
  Organization = names(row_totals),
  Total = row_totals,
  Yes_Number = cross_table[,"Yes"], 
  No_Number = cross_table[,"No"],
  Yes_Percent = (cross_table[, "Yes"] / row_totals) * 100,
  No_Percent = (cross_table[, "No"] / row_totals) * 100
)
print(summary_table)
