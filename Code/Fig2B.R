rm(list = ls())
options(stringsAsFactors = FALSE)

library(ggplot2)
library(ggpubr)

load("rld.RData")
load("meta.RData")

# --------------- Young --------------- 
group <- "young"
group_data <- meta[meta$sampletype == group,]
pcs <- plotPCA(rld, intgroup = "sampletype", returnData = TRUE)
group_pcs <- pcs[pcs$sampletype == group, c("PC1", "PC2")]
n <- nrow(group_pcs)
distances <- matrix(0, nrow = n, ncol = n)
for (i in 1:n) {
  for (j in 1:n) {
    distances[i,j] <- sqrt((group_pcs$PC1[i]-group_pcs$PC1[j])^2 + (group_pcs$PC2[i]-group_pcs$PC2[j])^2)
  }
}
young_dis <- as.numeric(distances)

# --------------- Adult --------------- 
group <- "adult"
group_data <- meta[meta$sampletype == group,]
pcs <- plotPCA(rld, intgroup = "sampletype", returnData = TRUE)
group_pcs <- pcs[pcs$sampletype == group, c("PC1", "PC2")]
n <- nrow(group_pcs)
distances <- matrix(0, nrow = n, ncol = n)
for (i in 1:n) {
  for (j in 1:n) {
    distances[i,j] <- sqrt((group_pcs$PC1[i]-group_pcs$PC1[j])^2 + (group_pcs$PC2[i]-group_pcs$PC2[j])^2)
  }
}
adult_dis <- as.numeric(distances)

# --------------- Old --------------- 
group <- "old"
group_data <- meta[meta$sampletype == group,]
pcs <- plotPCA(rld, intgroup = "sampletype", returnData = TRUE)
group_pcs <- pcs[pcs$sampletype == group, c("PC1", "PC2")]
n <- nrow(group_pcs)
distances <- matrix(0, nrow = n, ncol = n)
for (i in 1:n) {
  for (j in 1:n) {
    distances[i,j] <- sqrt((group_pcs$PC1[i]-group_pcs$PC1[j])^2 + (group_pcs$PC2[i]-group_pcs$PC2[j])^2)
  }
}
old_dis <- as.numeric(distances)
data <- data.frame(value = c(young_dis, adult_dis, old_dis), group = factor(rep(c("Young","Adult","Old"), c(36,25,36))))
data$group <- factor(data$group, levels=c("Young","Adult","Old"))
ggboxplot(data, x = "group", y = "value", 
          color = "group", add = "jitter",
          xlab = "Age", ylab = "Distance") + 
  stat_compare_means(comparisons = list(c("Young", "Adult"),
                                        c("Young", "Old"),
                                        c("Adult", "Old")),
                     method = "t.test",
                     label = "p.format", #"p.signif",
                     label.x = c(1.5, 2.5, 3.5), size = 3)+ 
  theme(legend.position="none")+
  scale_colour_manual(values = c("#339900", "#FF9900", "#006699"))+
  ylim(-10, 90)

stat.test <- data %>%
  t_test(value ~ group) %>%
  add_significance()
