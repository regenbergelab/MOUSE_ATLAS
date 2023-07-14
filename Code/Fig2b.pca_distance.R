# Mouse Atlas - Fig2b - pca_distance - Liver

rm(list = ls())
options(stringsAsFactors = FALSE)

# Set working folder
mainDir <- "/Volumes/Denmark/Mouse_Atlas/Bioinformatics/Figures/Fig1"
subDir <- "01_pca_distance"
# create a directory if it doesn't exist
ifelse(!dir.exists(file.path(mainDir, subDir)), dir.create(file.path(mainDir, subDir)), FALSE)
setwd(file.path(mainDir, subDir))

# Load required libraries
library(ggplot2)
library(ggpubr)

# Load required Rdata from DESeq2 pipeline
load("./inputfile/rld.RData")
load("./inputfile/meta.RData")

# --------------- Young --------------- 
group <- "young"
group_data <- meta[meta$sampletype == group,]
# Extract the coordinates of each point in PCA space
pcs <- plotPCA(rld, intgroup = "sampletype", returnData = TRUE)
group_pcs <- pcs[pcs$sampletype == group, c("PC1", "PC2")]
# Calculate the distance between any two points
n <- nrow(group_pcs)
distances <- matrix(0, nrow = n, ncol = n)
for (i in 1:n) {
  for (j in 1:n) {
    distances[i,j] <- sqrt((group_pcs$PC1[i]-group_pcs$PC1[j])^2 + (group_pcs$PC2[i]-group_pcs$PC2[j])^2)
  }
}
# Convert matrix to numeric
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

# Prepare the input data.frame for Box plot
data <- data.frame(value = c(young_dis, adult_dis, old_dis),
                   group = factor(rep(c("Young","Adult","Old"), c(36,25,36))))

# Box plot with p-values
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

ggsave("01_pca_distance_liver.pdf", width = 4, height = 4)

# Statistical test
stat.test <- data %>%
  t_test(value ~ group) %>%
  add_significance()
stat.test

#Returned value is a data frame with the following columns:
## .y.: the y variable used in the test.
## p: the p-value
## p.adj: the adjusted p-value. Default value for p.adjust.method = “holm”
## p.format: the formatted p-value
## p.signif: the significance level.
## method: the statistical test used to compare groups.
