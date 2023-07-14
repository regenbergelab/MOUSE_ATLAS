library(data.table)
library(dplyr)
library(ggplot2)
library(quantreg)
library(splines)
library(writexl)
library(janitor)
library(readr)
library(rtracklayer)
library(base)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GenomicFeatures")
library(purrr)
library(tidyverse)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ensembldb)
library(IDPmisc)
library(GenomicFeatures)
library(GenomicRanges)
library(biomaRt)
library(BiocInstaller)
library(org.Mm.eg.db)
library(noisyr)
library(rstatix)
library(ggsignif)
library(readxl)
library(ggpubr)

## THIS ANALYSIS ONLY USES ECCDNA DERIVED FROM PROTEIN CODING GENES
#annotations from biomart download did not have nonconding genes and miRNA
#load gene annotations
#downloaded from http://ftp.ensembl.org/pub/release-100/gtf/mus_musculus/Mus_musculus.GRCm38.100.gtf.gz
setwd("C:/Users/vzl597/Downloads")
gtf <- rtracklayer::import(paste0("Mus_musculus.GRCm38.100.gtf.gz"))
gtf_df<-as.data.frame(gtf)
gtf_TEC <- gtf_df %>% 
  filter(gene_biotype == 'TEC')
gtf_df$gene_biotype %>% unique
#gene count in the    
gtf_df$gene_name %>% unique %>% length
gtf_df["gene_name_or_id"]<-gtf$gene_name
gtf_df$gene_name_or_id[is.na(gtf_df$gene_name_or_id)]<-gtf$gene_id[is.na(gtf_df["gene_name_or_id"])]
genes<-gtf_df[gtf_df$type=="gene",]
#check what is left that is duplicated
genes$gene_name_or_id[genes$gene_name_or_id %>% duplicated %>% which] %>% table
#It's safe to remove the duplicates
genes2<-genes[!duplicated(genes$gene_name_or_id),]
genes3<-genes2[!is.na(genes$gene_name),]
genes4 <- genes3[,-(4:27),drop=FALSE]
genes4
genes5 <- genes4[,-c(4,5,6,7)]
genes6 <- na.omit(genes5)
genes6
names(genes6)[names(genes6) == "seqnames"] <- "chr"
genes6$slen <- genes6$end - genes6$start

## FIND NUMBER OF EXONS PER GENE
#annotations from biomart download did not have nonconding genes and miRNA
#load gene annotations
#downloaded from http://ftp.ensembl.org/pub/release-100/gtf/mus_musculus/Mus_musculus.GRCm38.100.gtf.gz
setwd("C:/Users/vzl597/Downloads")
gtf <- rtracklayer::import(paste0("Mus_musculus.GRCm38.100.gtf.gz"))
txdb <- makeTxDbFromGFF("Mus_musculus.GRCm38.100.gtf.gz")
exons.list.per.gene <- exonsBy(txdb, by="gene")
transcripts.list.per.gene <- as.data.frame(transcriptsBy(txdb, by="gene"))
introns.list.per.gene <- intronsByTranscript(txdb)
introns <- unlist(introns.list.per.gene)
exonic.gene.sizes <- as.data.frame(width(GenomicRanges::reduce(exons.list.per.gene)))
count <- exonic.gene.sizes %>%
  group_by(group_name) %>%
  summarise(length(value))
GENOMEWIDEnexons <- as.data.frame(count)
count.transcripts <- transcripts.list.per.gene %>%
  group_by(group_name) %>%
  summarise(length(group_name))
GENOMEWIDEnexons <- as.data.frame(count)
GENOMEWIDEntranscripts <- as.data.frame(count.transcripts)

mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
genes <- GENOMEWIDEnexons$group_name
G_list <- getBM(filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "mgi_symbol"), values=genes, mart= mart)

GENOMEWIDEnexons1 <- left_join(G_list, GENOMEWIDEnexons, by=join_by(ensembl_gene_id ==  group_name))
exonnumb <- GENOMEWIDEnexons1[,-1]
names(exonnumb)[names(exonnumb) == "length(value)"] <- "number"
numberex <- as.numeric(as.character(exonnumb$number))
df3 <- cbind(exonnumb, numberex)
ExonAllGenes <- df3[,-2]
ExonAllGenes[order(ExonAllGenes$numberex, decreasing=T),]  

GENOMEWIDEntranscripts1 <- left_join(G_list, GENOMEWIDEntranscripts, by=join_by(ensembl_gene_id ==  group_name))
transcriptn <- GENOMEWIDEntranscripts1[,-1]
names(transcriptn)[names(transcriptn) == "length(group_name)"] <- "number"
transcriptn[order(transcriptn$number, decreasing=T),]  

genesinfo <- left_join(genes6, transcriptn, by=join_by(gene_name_or_id == mgi_symbol))
genesinfo2 <- left_join(genesinfo, ExonAllGenes, by=join_by(gene_name_or_id == mgi_symbol))

## Length of all intronic and exonic sequences per gene
exonic.gene.sizes <- as.data.frame(width(GenomicRanges::reduce(exons.list.per.gene)))
Exonlength <- exonic.gene.sizes %>%
  group_by(group_name) %>%
  summarise(sum(value))
Exonlength <- left_join(G_list, Exonlength, by=join_by(ensembl_gene_id ==  group_name))
TotalExonlength <- Exonlength[,-1]

GeneInfo <- left_join(genesinfo2, TotalExonlength, by=join_by(gene_name_or_id == mgi_symbol))
colnames(GeneInfo) <- c("chr", "start", "end", "Gene", "slen", "TranscriptNumber", "MaxExonNumber", "Exonlength")
GeneInfo$IntronicLength <- GeneInfo$slen - (GeneInfo$Exonlength - 1)
GeneInfo$RatioExTo <- (GeneInfo$Exonlength-1)*100/GeneInfo$slen
GeneInfo$RatioExIn <- (GeneInfo$Exonlength-1)*100/GeneInfo$IntronicLength

## FOR CORTEX
setwd("C:/Users/vzl597/Desktop/Xue Analysis")
cortex3 <- read_excel("Mb that were found in the 3 month and 12 month liver, hippocampus and cortex.xlsx", 4)
YoungDep <- cortex3[cortex3$label == "below",]
YoungDep <- left_join(YoungDep, GeneInfo, by = c("gene_name" = "Gene"))

YoungEnriched <- cortex3[cortex3$label == "above",]
YoungEnriched <- left_join(YoungEnriched, GeneInfo, by=c("gene_name" = "Gene"))

YoungBulk <- cortex3[cortex3$label == "middle",]
YoungBulk <- left_join(YoungBulk, GeneInfo, by=c("gene_name" = "Gene"))

setwd("C:/Users/vzl597/Desktop/Xue Analysis")
cortex12 <- read_excel("Mb that were found in the 3 month and 12 month liver, hippocampus and cortex.xlsx", 5)
AdultDep <- cortex12[cortex12$label == "below",]
AdultDep <- left_join(AdultDep, GeneInfo, by=c("gene_name" = "Gene"))

AdultEnriched <- cortex12[cortex12$label == "above",]
AdultEnriched <- left_join(AdultEnriched, GeneInfo, by=c("gene_name" = "Gene"))

AdultBulk <- cortex12[cortex12$label == "middle",]
AdultBulk <- left_join(AdultBulk, GeneInfo, by=c("gene_name" = "Gene"))

setwd("C:/Users/vzl597/Desktop/Xue Analysis")
cortex22 <- read_excel("Mb that were found in the 3 month and 12 month liver, hippocampus and cortex.xlsx", 6)
OldDep <- cortex22[cortex22$label == "below",]
OldDep <- left_join(OldDep, GeneInfo, by=c("gene_name" = "Gene"))

OldEnriched <- cortex22[cortex22$label == "above",]
OldEnriched <- left_join(OldEnriched, GeneInfo, by=c("gene_name" = "Gene"))

OldBulk <- cortex22[cortex22$label == "middle",]
OldBulk <- left_join(OldBulk, GeneInfo, by=c("gene_name" = "Gene"))

OldDep$Age <- "Old"
AdultDep$Age <- "Adult"
YoungDep$Age <- "Young"
DataFrameDepleted <- rbind(OldDep, AdultDep, YoungDep)

DataFrameDepleted <- DataFrameDepleted %>%
  mutate(Age = fct_relevel(Age,"Young", "Adult", "Old")) 


OldEnriched$Age <- "Old"
AdultEnriched$Age <- "Adult"
YoungEnriched$Age <- "Young"
DataFrameEnriched <- rbind(OldEnriched, AdultEnriched, YoungEnriched)

DataFrameEnriched <- DataFrameEnriched %>%
  mutate(Age = fct_relevel(Age,"Young", "Adult", "Old")) 

OldBulk$Age <- "Old"
AdultBulk$Age <- "Adult"
YoungBulk$Age <- "Young"
DataFrameBulk <- rbind(OldBulk, AdultBulk, YoungBulk)
DataFrameBulk <- DataFrameBulk %>%
  mutate(Age = fct_relevel(Age,"Young", "Adult", "Old")) 

## Manually figure out the reapeated genes. I did it using Excel
CommonGenes <- subset(GeneInfo, Gene %in% c("2210408I21Rik", "4930402H24Rik", "6330411D24Rik", "Astn1", "Cacna1a", "Cacnb2", "Caln1", "Cfap299", "Clasp2", "Eda", "Frmpd4", "Grid1", "Hpse2", "Il1rapl1", "Kidins220", "Kif1b", "Limch1", "Lrmda", "Nav1", "Nebl", "Ntrk2", "Pak3", "Prkce", "Rapgef5", "Rbms3", "Rgs6", "Rimbp2", "Shank2", "Slc1a2", "Spock1", "Synpr", "Tmem178b", "Tshz2", "Usp32", "Unc80", "Wwox", "Xpr1") )
AllDataFrame <- rbind(DataFrameBulk, DataFrameDepleted, DataFrameEnriched)

AllDataFrame$ExonsperMB <- AllDataFrame$MaxExonNumber*1000000/AllDataFrame$slen
AllDataFrame$Exonsdensity <- AllDataFrame$Exonlength*1000/AllDataFrame$IntronicLength
AllDataFrame$Introndensity <- AllDataFrame$IntronicLength/AllDataFrame$slen

Bulkset$Introndensity <- Bulkset$IntronicLength/Bulkset$slen
Depletedset$Introndensity <- Depletedset$IntronicLength/Depletedset$slen
Genestoplot$Introndensity <- Genestoplot$IntronicLength/Genestoplot$slen
AllDataFrame$Introndensity <- AllDataFrame$IntronicLength/AllDataFrame$slen
AllDataFrame$IsoformperExonNumberperMB <- AllDataFrame$TranscriptNumber/AllDataFrame$ExonsperMB  

gene_name <- c("2210408I21Rik", "4930402H24Rik", "6330411D24Rik", "Astn1", "Cacna1a", "Cacnb2", "Caln1", "Cfap299", "Clasp2", "Eda", "Frmpd4", "Grid1", "Hpse2", "Il1rapl1", "Kidins220", "Kif1b", "Limch1", "Lrmda", "Nav1", "Nebl", "Ntrk2", "Pak3", "Prkce", "Rapgef5", "Rbms3", "Rgs6", "Rimbp2", "Shank2", "Slc1a2", "Spock1", "Synpr", "Tmem178b", "Tshz2", "Usp32", "Unc80", "Wwox", "Xpr1")
RecurrentGenes <- data.frame(gene_name)
RecurrentGenesInfo <- left_join(RecurrentGenes, AllDataFrame, by = c("gene_name"= "gene_name"))
AllDataFrame <- data.frame(AllDataFrame)
AllDataFrametoPlot <- AllDataFrame[!AllDataFrame %in% RecurrentGenesInfo,]
RecurrentGenesInfo$eccDNAlevel <- "Recurrent"
ReadytoPlot <- rbind(AllDataFrametoPlot, RecurrentGenesInfo)

ReadytoPlot <- subset(ReadytoPlot, eccDNAlevel != "Enriched")
colnames(ReadytoPlot)
print(levels(ReadytoPlot$eccDNAlevel))

predefined_order <- c("Young", "Adult", "Old")

# Reorder the levels within the "Category" column based on the predefined order
ReadytoPlot$Age <- factor(ReadytoPlot$Age, levels = predefined_order)

ggplot(ReadytoPlot, aes(x = eccDNAlevel , y=Introndensity, fill=eccDNAlevel )) +
  geom_boxplot(outlier.shape = NA) + 
  scale_fill_manual(values = c("#999999", "#a7dcdd", "#3e7b7d")) +
  scale_y_continuous() +
  ylim(c(0.5,1.1))+
  labs(x = "Age and Condition", y = "Total Intron Length / Gene Length", title = "CORTEX - Intron density") +
  facet_grid(. ~ Age) 

# Perform Kruskal-Wallis test
stats <- compare_means(Introndensity ~ eccDNAlevel, data = ReadytoPlot, group.by = "Age", paired = FALSE, na.rm = TRUE, method="kruskal.test")

# Perform pairwise comparisons within each age group
pairwise <- lapply(unique(ReadytoPlot$Age), function(age) {
  subset <- ReadytoPlot[ReadytoPlot$Age == age, ]
  pairwise.wilcox.test(subset$Introndensity, subset$eccDNAlevel, p.adjust.method = "bonferroni")
})

ggplot(ReadytoPlot, aes(x = eccDNAlevel , y=log(TranscriptNumber,10), fill=eccDNAlevel )) +
  geom_boxplot(outlier.shape = NA) + 
  scale_fill_manual(values = c("#999999", "#a7dcdd", "#3e7b7d")) +
  scale_y_continuous() +
  ylim(c(0,1.5))+
  labs(x = "Age and Condition", y = "Number of Transcript Isoforms (log10)", title = "CORTEX - Number of Transcripts") +
  facet_grid(. ~ Age) 

# Perform Kruskal-Wallis test
stats <- compare_means(TranscriptNumber ~ eccDNAlevel, data = ReadytoPlot, group.by = "Age", paired = FALSE, na.rm = TRUE, method="kruskal.test")

# Perform pairwise comparisons within each age group
pairwise <- lapply(unique(ReadytoPlot$Age), function(age) {
  subset <- ReadytoPlot[ReadytoPlot$Age == age, ]
  pairwise.wilcox.test(subset$TranscriptNumber, subset$eccDNAlevel, p.adjust.method = "bonferroni")
})

ggplot(ReadytoPlot, aes(x = eccDNAlevel , y=IsoformperExonNumberperMB, fill=eccDNAlevel )) +
  geom_boxplot(outlier.shape = NA) + 
  scale_fill_manual(values = c("#999999", "#a7dcdd", "#3e7b7d")) +
  scale_y_continuous() +
  ylim(c(0,0.75))+
  labs(x = "Age and Condition", y = "Isoforms per Exon Number per Mb", title = "CORTEX - Number of Isoforms per Exon Number per Mb") +
  facet_grid(. ~ Age)

# Perform Kruskal-Wallis test
stats <- compare_means(IsoformperExonNumberperMB ~ eccDNAlevel, data = ReadytoPlot, group.by = "Age", paired = FALSE, na.rm = TRUE, method="kruskal.test")

# Perform pairwise comparisons within each age group
pairwise <- lapply(unique(ReadytoPlot$Age), function(age) {
  subset <- ReadytoPlot[ReadytoPlot$Age == age, ]
  pairwise.wilcox.test(subset$IsoformperExonNumberperMB, subset$eccDNAlevel, p.adjust.method = "bonferroni")
})

## FOR HIPPOCAMPUS
setwd("C:/Users/vzl597/Desktop/Xue Analysis")
hipp3 <- read_excel("Mb that were found in the 3 month and 12 month liver, hippocampus and cortex.xlsx", 7)
YoungDep1 <- hipp3[hipp3$label == "below",]
YoungDep1 <- left_join(YoungDep1, GeneInfo, by = c("gene_name" = "Gene"))

YoungEnriched1 <- hipp3[hipp3$label == "above",]
YoungEnriched1 <- left_join(YoungEnriched1, GeneInfo, by=c("gene_name" = "Gene"))

YoungBulk1 <- hipp3[hipp3$label == "middle",]
YoungBulk1 <- left_join(YoungBulk1, GeneInfo, by=c("gene_name" = "Gene"))

setwd("C:/Users/vzl597/Desktop/Xue Analysis")
hipp12 <- read_excel("Mb that were found in the 3 month and 12 month liver, hippocampus and cortex.xlsx", 8)
AdultDep <- hipp12[hipp12$label == "below",]
AdultDep <- left_join(AdultDep, GeneInfo, by=c("gene_name" = "Gene"))

AdultEnriched <- hipp12[hipp12$label == "above",]
AdultEnriched <- left_join(AdultEnriched, GeneInfo, by=c("gene_name" = "Gene"))

AdultBulk <- hipp12[hipp12$label == "middle",]
AdultBulk <- left_join(AdultBulk, GeneInfo, by=c("gene_name" = "Gene"))

setwd("C:/Users/vzl597/Desktop/Xue Analysis")
hipp22 <- read_excel("Mb that were found in the 3 month and 12 month liver, hippocampus and cortex.xlsx", 9)
OldDep <- hipp22[hipp22$label == "below",]
OldDep <- left_join(OldDep, GeneInfo, by=c("gene_name" = "Gene"))

OldEnriched <- hipp22[hipp22$label == "above",]
OldEnriched <- left_join(OldEnriched, GeneInfo, by=c("gene_name" = "Gene"))

OldBulk <- hipp22[hipp22$label == "middle",]
OldBulk <- left_join(OldBulk, GeneInfo, by=c("gene_name" = "Gene"))

OldDep$Age <- "Old"
AdultDep$Age <- "Adult"
YoungDep1$Age <- "Young"
DataFrameDepleted <- rbind(OldDep, AdultDep, YoungDep1)

DataFrameDepleted <- DataFrameDepleted %>%
  mutate(Age = fct_relevel(Age,"Young", "Adult", "Old"))

OldEnriched$Age <- "Old"
AdultEnriched$Age <- "Adult"
YoungEnriched1$Age <- "Young"
DataFrameEnriched <- rbind(OldEnriched, AdultEnriched, YoungEnriched1)

DataFrameEnriched <- DataFrameEnriched %>%
  mutate(Age = fct_relevel(Age,"Young", "Adult", "Old")) 

OldBulk$Age <- "Old"
AdultBulk$Age <- "Adult"
YoungBulk1$Age <- "Young"
DataFrameBulk <- rbind(OldBulk, AdultBulk, YoungBulk1)
DataFrameBulk <- DataFrameBulk %>%
  mutate(Age = fct_relevel(Age,"Young", "Adult", "Old")) 

## Manually figure out the reapeated genes. I did it using Excel
CommonGenes <- subset(GeneInfo, Gene %in% c("Ank2", "App", "Arnt2", "Cacna1a", "Camta1", "Cdk14", "Celf2", "Cntn5", "Cpq", "Dach2", "Dlgap1", "Dst", "Eif4g3", "Ephb1", "Fmn1", "Frmpd4", "Gpc5", "Grm5", "Iqcm", "Iqsec1", "L3mbtl4", "Lrmda", "Mycbp2", "Ncald", "Ntm", "Pak3", "Rp1", "Rtn1", "Scn8a", "Sipa1l1", "Srgap3", "Sv2b", "Tcf4", "Tenm1", "Trio", "Tshz2", "Vsnl1"))

DataFrameBulk$eccDNAlevel <- "Bulk"
DataFrameDepleted$eccDNAlevel <- "Depleted"
DataFrameEnriched$eccDNAlevel <- "Enriched"

AllDataFrame <- rbind(DataFrameBulk, DataFrameDepleted, DataFrameEnriched)

AllDataFrame$ExonsperMB <- AllDataFrame$MaxExonNumber*1000000/AllDataFrame$slen
AllDataFrame$Exonsdensity <- AllDataFrame$Exonlength*1000/AllDataFrame$IntronicLength
AllDataFrame$Introndensity <- AllDataFrame$IntronicLength/AllDataFrame$slen
AllDataFrame$Introndensity <- AllDataFrame$IntronicLength/AllDataFrame$slen
AllDataFrame$IsoformperExonNumberperMB <- AllDataFrame$TranscriptNumber/AllDataFrame$ExonsperMB  

gene_name <- c("Ank2", "App", "Arnt2", "Cacna1a", "Camta1", "Cdk14", "Celf2", "Cntn5", "Cpq", "Dach2", "Dlgap1", "Dst", "Eif4g3", "Ephb1", "Fmn1", "Frmpd4", "Gpc5", "Grm5", "Iqcm", "Iqsec1", "L3mbtl4", "Lrmda", "Mycbp2", "Ncald", "Ntm", "Pak3", "Rp1", "Rtn1", "Scn8a", "Sipa1l1", "Srgap3", "Sv2b", "Tcf4", "Tenm1", "Trio", "Tshz2", "Vsnl1")
RecurrentGenes <- data.frame(gene_name)
RecurrentGenesInfo <- left_join(RecurrentGenes, AllDataFrame, by = c("gene_name"= "gene_name"))
AllDataFrame <- data.frame(AllDataFrame)
AllDataFrametoPlot <- AllDataFrame[!AllDataFrame %in% RecurrentGenesInfo,]
RecurrentGenesInfo$eccDNAlevel <- "Recurrent"
ReadytoPlot <- rbind(AllDataFrametoPlot, RecurrentGenesInfo)

ReadytoPlot <- subset(ReadytoPlot, eccDNAlevel != "Enriched")
colnames(ReadytoPlot)
print(levels(ReadytoPlot$eccDNAlevel))

predefined_order <- c("Young", "Adult", "Old")

# Reorder the levels within the "Category" column based on the predefined order
ReadytoPlot$Age <- factor(ReadytoPlot$Age, levels = predefined_order)

ggplot(ReadytoPlot, aes(x = eccDNAlevel , y=Introndensity, fill=eccDNAlevel )) +
  geom_boxplot(outlier.shape = NA) + 
  scale_fill_manual(values = c("#999999", "#a7dcdd", "#3e7b7d")) +
  scale_y_continuous() +
  ylim(c(0.5,1.1))+
  labs(x = "Age and Condition", y = "Total Intron Length / Gene Length", title = "HIPPOCAMPUS - Intron density") +
  facet_grid(. ~ Age) 

# Perform Kruskal-Wallis test
stats <- compare_means(Introndensity ~ eccDNAlevel, data = ReadytoPlot, group.byhttp://127.0.0.1:16093/graphics/ab27c739-2d36-4983-930a-3350983a8add.png = "Age", paired = FALSE, na.rm = TRUE, method="kruskal.test")

# Perform pairwise comparisons within each age group
pairwise <- lapply(unique(ReadytoPlot$Age), function(age) {
  subset <- ReadytoPlot[ReadytoPlot$Age == age, ]
  pairwise.wilcox.test(subset$Introndensity, subset$eccDNAlevel, p.adjust.method = "bonferroni")
})

ggplot(ReadytoPlot, aes(x = eccDNAlevel , y=log(TranscriptNumber,10), fill=eccDNAlevel )) +
  geom_boxplot(outlier.shape = NA) + 
  scale_fill_manual(values = c("#999999", "#a7dcdd", "#3e7b7d")) +
  scale_y_continuous() +
  ylim(c(0,1.5))+
  labs(x = "Age and Condition", y = "Number of Transcript Isoforms (log10)", title = "HIPPOCAMPUS - Number of Transcripts") +
  facet_grid(. ~ Age) 

# Perform Kruskal-Wallis test
stats <- compare_means(TranscriptNumber ~ eccDNAlevel, data = ReadytoPlot, group.by = "Age", paired = FALSE, na.rm = TRUE, method="kruskal.test")

# Perform pairwise comparisons within each age group
pairwise <- lapply(unique(ReadytoPlot$Age), function(age) {
  subset <- ReadytoPlot[ReadytoPlot$Age == age, ]
  pairwise.wilcox.test(subset$TranscriptNumber, subset$eccDNAlevel, p.adjust.method = "bonferroni")
})

ggplot(ReadytoPlot, aes(x = eccDNAlevel , y=IsoformperExonNumberperMB, fill=eccDNAlevel )) +
  geom_boxplot(outlier.shape = NA) + 
  scale_fill_manual(values = c("#999999", "#a7dcdd", "#3e7b7d")) +
  scale_y_continuous() +
  ylim(c(0,0.75))+
  labs(x = "Age and Condition", y = "Isoforms per Exon Number per Mb", title = "HIPPOCAMPUS - Number of Isoforms per Exon Number per Mb") +
  facet_grid(. ~ Age)

# Perform Kruskal-Wallis test
stats <- compare_means(IsoformperExonNumberperMB ~ eccDNAlevel, data = ReadytoPlot, group.by = "Age", paired = FALSE, na.rm = TRUE, method="kruskal.test")

# Perform pairwise comparisons within each age group
pairwise <- lapply(unique(ReadytoPlot$Age), function(age) {
  subset <- ReadytoPlot[ReadytoPlot$Age == age, ]
  pairwise.wilcox.test(subset$IsoformperExonNumberperMB, subset$eccDNAlevel, p.adjust.method = "bonferroni")
})
