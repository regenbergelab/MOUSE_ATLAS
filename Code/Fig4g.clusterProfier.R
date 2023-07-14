# clusterProfiler

# Install packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Mm.eg.db")

# Import packages
library(clusterProfiler)
library("org.Mm.eg.db")
library(ggplot2)

# Genelist
# ---------- For the cortex, 37 common elements in "Young", "Adult" and "Old":
#genelist_up <- c("Ank2","App","Arnt2","Cacna1a","Camta1","Cdk14","Celf2","Cntn5","Cpq","Dach2","Dlgap1","Dst","Eif4g3","Ephb1","Fmn1","Frmpd4","Gpc5","Grm5","Iqcm","Iqsec1","L3mbtl4","Lrmda","Mycbp2","Ncald","Ntm","Pak3","Rp1","Rtn1","Scn8a","Sipa1l1","Srgap3","Sv2b","Tcf4","Tenm1","Trio","Tshz2","Vsnl1")

# -------- For the hippocampus, 37 common elements in "Young", "Adult" and "Old":
#genelist_up <- c("Ank2","App","Arnt2","Cacna1a","Camta1","Cdk14","Celf2","Cntn5","Cpq","Dach2","Dlgap1","Dst","Eif4g3","Ephb1","Fmn1","Frmpd4","Gpc5","Grm5","Iqcm","Iqsec1","L3mbtl4","Lrmda","Mycbp2","Ncald","Ntm","Pak3","Rp1","Rtn1","Scn8a","Sipa1l1","Srgap3","Sv2b","Tcf4","Tenm1","Trio","Tshz2","Vsnl")

# ---------- For the liver, all 910 protected genes
genelist_up_df <- read.csv("./liver_proteced_genellist.csv", sep = "\t", header = F)
genelist_up <- genelist_up_df[,1]

# ---------- For the cortex, all 962 protected genes
genelist_up_df <- read.csv("./cortex_proteced_genellist.csv", sep = "\t", header = F)
genelist_up <- genelist_up_df[,1]

# ---------- For the hippocampus, 734 XX protected genes
genelist_up_df <- read.csv("./hippocampus_proteced_genellist.csv", sep = "\t", header = F)
genelist_up <- genelist_up_df[,1]

# GO
egosimp_up <- enrichGO(genelist_up, 
                       OrgDb=org.Mm.eg.db, 
                       ont='ALL',
                       pAdjustMethod='BH', 
                       pvalueCutoff=0.05, 
                       qvalueCutoff=0.2, 
                       keyType='SYMBOL')

# Tree
library(ggnewscale)
library(enrichplot)
edox <- setReadable(egosimp_up, 'org.Mm.eg.db', 'SYMBOL')
edox2 <- pairwise_termsim(edox)
treeplot(edox2)
#treeplot(x, showCategory = 30, color = "p.adjust", label_format = 30, ...)

# Save PDF
#ggsave("Cortex_37_common_protected_gene.pdf", width = 12, height = 9)
#ggsave("Hippocampus_37_common_protected_gene.pdf", width = 12, height = 9)
ggsave("ClusterProfier4_910_protectedGenes_liver.pdf", width = 12, height = 9)
ggsave("ClusterProfier4_962_protectedGenes_cortex.pdf", width = 12, height = 9)
ggsave("ClusterProfier4_734_protected_gene_hippocampus.pdf", width = 12, height = 9)
