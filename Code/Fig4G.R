library(clusterProfiler)
library("org.Mm.eg.db")
library(ggplot2)
library(ggnewscale)
library(enrichplot)

egosimp_up <- enrichGO(genelist_up, 
                       OrgDb=org.Mm.eg.db, 
                       ont='ALL',
                       pAdjustMethod='BH', 
                       pvalueCutoff=0.05, 
                       qvalueCutoff=0.2, 
                       keyType='SYMBOL')
edox <- setReadable(egosimp_up, 'org.Mm.eg.db', 'SYMBOL')
edox2 <- pairwise_termsim(edox)
treeplot(edox2)
