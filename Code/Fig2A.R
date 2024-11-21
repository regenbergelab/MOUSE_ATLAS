rm(list = ls())
options(stringsAsFactors = FALSE)

library(DESeq2)
library(ggplot2) 

load("rld.RData")

PCAdata <- plotPCA(rld, intgroup = "sampletype", returnData=TRUE, ntop = 500) 

ggplot(PCAdata, aes(x = PC1, y = PC2, color = sampletype)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("#FF9900", "#006699","#339900"))+
  geom_text_repel(aes(label = rownames(PCAdata))) +
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("PC1 (58% variance)") +
  ylab("PC2 (18% variance)") +
  theme_bw()
