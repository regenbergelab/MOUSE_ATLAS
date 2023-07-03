# Fig.2a - PCA plot (liver)

rm(list = ls())
options(stringsAsFactors = FALSE)

# Load required libraries
library(DESeq2)
library(ggplot2) 

# Load required Rdata from DESeq2 pipeline
load("./inputfile/rld.RData")

# Sample PCA plot for transformed data
PCAdata <- plotPCA(rld, 
                   intgroup = "sampletype", # Use age for grouping here
                   returnData=TRUE,  # Save the output of plotPCA() to a variable by specifying the returnData=TRUE argument
                   ntop = 500)  # By default the function uses the top 500 most variable genes

# PCA plot
ggplot(PCAdata, aes(x = PC1, y = PC2, color = sampletype)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("#FF9900", "#006699","#339900"))+
  geom_text_repel(aes(label = rownames(PCAdata))) +
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("PC1 (58% variance)") +
  ylab("PC2 (18% variance)") +
  theme_bw()

ggsave("02_pca.pdf", width = 6, height = 5)
