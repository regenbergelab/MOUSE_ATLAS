library(data.table)
library(dplyr)
library(ggplot2)
library(quantreg)
library(rGREAT)
library(splines)
library(writexl)
library(janitor)
library(readr)
library(rtracklayer)
library(base)
library(NLP)
library(purrr)
library(tidyverse)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ensembldb)
library(GenomicFeatures)
library(GenomicRanges)
library(biomaRt)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("rtracklayer")

## THIS ANALYSIS ONLY USES ECCDNA DERIVED FROM PROTEIN CODING GENES
#annotations from biomart download did not have nonconding genes and miRNA
#load gene annotations
#downloaded from http://ftp.ensembl.org/pub/release-100/gtf/mus_musculus/Mus_musculus.GRCm38.100.gtf.gz
setwd("C:/Users/vzl597/Downloads")
gtf <- rtracklayer::import(paste0("Mus_musculus.GRCm38.100.gtf.gz"))
gtf_df<-as.data.frame(gtf)
gtf_df
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

##FOR LIVER
## YOUNG DF FRAME
setwd("C:/Users/vzl597/Desktop/Xue Analysis")
liver3m <- read_excel("Mb that were found in the 3 month and 12 month liver, hippocampus and cortex.xlsx", 1)
lengthinfo <- left_join(liver3m, genes6, by = c("gene_name" = "gene_name_or_id"))
lengthinfo$slen <- lengthinfo$end - lengthinfo$start
lengthinfo$eccDNA <- round(exp(lengthinfo$log_ecc_number_per_MB ) * (lengthinfo$slen / 1000000))
lengthinfo$eccDNAperMB <- lengthinfo$eccDNA / (lengthinfo$slen/1000000)
lengthinfo$norm_count <- exp(lengthinfo$log_norm_count_per_MB)*(lengthinfo$slen / 1000000)
lengthinfo$ratio <- log(lengthinfo$eccDNAperMB / lengthinfo$norm_count,2)
median(lengthinfo$ratio)

liver3m <- lengthinfo
liver3m$age <- "Young"

liver12m <- read_excel("Mb that were found in the 3 month and 12 month liver, hippocampus and cortex.xlsx", 2)
lengthinfo <- left_join(liver12m, genes6, by = c("gene_name" = "gene_name_or_id"))
lengthinfo$slen <- lengthinfo$end - lengthinfo$start
lengthinfo$eccDNA <- round(exp(lengthinfo$log_ecc_number_per_MB ) * (lengthinfo$slen / 1000000))
lengthinfo$eccDNAperMB <- lengthinfo$eccDNA / (lengthinfo$slen/1000000)
lengthinfo$norm_count <- exp(lengthinfo$log_norm_count_per_MB)*(lengthinfo$slen / 1000000)
lengthinfo$ratio <- log(lengthinfo$eccDNAperMB / lengthinfo$norm_count,2)
median(lengthinfo$ratio)

liver12m <- lengthinfo
liver12m$age <- "Adult"

liver22m <- read_excel("Mb that were found in the 3 month and 12 month liver, hippocampus and cortex.xlsx", 3)
lengthinfo <- left_join(liver22m, genes6, by = c("gene_name" = "gene_name_or_id"))
lengthinfo$slen <- lengthinfo$end - lengthinfo$start
lengthinfo$eccDNA <- round(exp(lengthinfo$log_ecc_number_per_MB ) * (lengthinfo$slen / 1000000))
lengthinfo$eccDNAperMB <- lengthinfo$eccDNA / (lengthinfo$slen/1000000)
lengthinfo$norm_count <- exp(lengthinfo$log_norm_count_per_MB)*(lengthinfo$slen / 1000000)
lengthinfo$ratio <- log(lengthinfo$eccDNAperMB / lengthinfo$norm_count,2)
median(lengthinfo$ratio)

liver22m <- lengthinfo
liver22m$age <- "Old"

liver <- rbind(liver3m, liver12m, liver22m)
liver <- data.frame(liver)
liver$ecc_number_per_MB <- 2^(liver$log_ecc_number_per_MB)
liver$norm_count_per_MB  <- 2^(liver$log_norm_count_per_MB)

liver$check <- log(liver$ecc_number_per_MB,2)
liver$check2 <- log(liver$norm_count_per_MB,2)

liver$ratio <- (liver$ecc_number_per_MB/liver$norm_count_per_MB)

liver %>%
  group_by(age)
data.frame(liver)
ratiog <- liver[,c(14,15)]

ratiog$logratio <- log(ratiog$ratio,2)
data=ratiog

ratiog %>%
  group_by(age) %>%
  summarise(across(.cols=everything(), list(median=median,sd=sd)))


data$category <- factor(data$age, levels = c("Young", "Adult", "Old"))


my_comparisons <- list(c("Young", "Adult"), c("Adult", "Old"), c("Young", "Old"))
pairwise <- pairwise.wilcox.test(ratiog$logratio, ratiog$age, p.adjust.method = "bonferroni")

ggplot(data = ratiog) +
  aes(y = logratio, x = factor(age, level = c("Young", "Adult", "Old"))) +
  scale_fill_viridis_d(option = "D") +
  geom_point(aes(fill = age), shape = 21, size = 1.5, alpha = 0.2, position = "jitter", color = "black") +
  geom_violin(aes(fill = age), alpha = 0.75, size = 1, width = 1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label.y = 2) +
  ggtitle("Liver") +
  labs(y = "log2(eccDNA_count/norm_count_per_MB)", x = "Age") +
  stat_summary(fun = "median", geom = "crossbar", width = 0.5, colour = "white") +
  geom_signif(comparisons = pairwise$by, y_position = 3)

##FOR CORTEX
## YOUNG DF FRAME
setwd("C:/Users/vzl597/Desktop/Xue Analysis")
cortex3m <- read_excel("Mb that were found in the 3 month and 12 month liver, hippocampus and cortex.xlsx", 4)
lengthinfo <- left_join(cortex3m, genes6, by = c("gene_name" = "gene_name_or_id"))
lengthinfo$slen <- lengthinfo$end - lengthinfo$start
lengthinfo$eccDNA <- round(exp(lengthinfo$log_ecc_number_per_MB ) * (lengthinfo$slen / 1000000))
lengthinfo$eccDNAperMB <- lengthinfo$eccDNA / (lengthinfo$slen/1000000)
lengthinfo$norm_count <- exp(lengthinfo$log_norm_count_per_MB)*(lengthinfo$slen / 1000000)
lengthinfo$ratio <- log(lengthinfo$eccDNAperMB / lengthinfo$norm_count,2)
mean(lengthinfo$ratio)

cortex3m <- lengthinfo
cortex3m$age <- "Young"

cortex12 <- read_excel("Mb that were found in the 3 month and 12 month liver, hippocampus and cortex.xlsx", 5)
lengthinfo <- left_join(cortex12, genes6, by = c("gene_name" = "gene_name_or_id"))
lengthinfo$slen <- lengthinfo$end - lengthinfo$start
lengthinfo$eccDNA <- round(exp(lengthinfo$log_ecc_number_per_MB ) * (lengthinfo$slen / 1000000))
lengthinfo$eccDNAperMB <- lengthinfo$eccDNA / (lengthinfo$slen/1000000)
lengthinfo$norm_count <- exp(lengthinfo$log_norm_count_per_MB)*(lengthinfo$slen / 1000000)
lengthinfo$ratio <- log(lengthinfo$eccDNAperMB / lengthinfo$norm_count,2)
mean(lengthinfo$ratio)

cortex12m <- lengthinfo
cortex12m$age <- "Adult"
install.packages("writexl")

cortex22m <- read_excel("Mb that were found in the 3 month and 12 month liver, hippocampus and cortex.xlsx", 6)
lengthinfo <- left_join(cortex22m, genes6, by = c("gene_name" = "gene_name_or_id"))
lengthinfo$slen <- lengthinfo$end - lengthinfo$start
lengthinfo$eccDNA <- round(exp(lengthinfo$log_ecc_number_per_MB ) * (lengthinfo$slen / 1000000))
lengthinfo$eccDNAperMB <- lengthinfo$eccDNA / (lengthinfo$slen/1000000)
lengthinfo$norm_count <- exp(lengthinfo$log_norm_count_per_MB)*(lengthinfo$slen / 1000000)
lengthinfo$ratio <- log(lengthinfo$eccDNAperMB / lengthinfo$norm_count,2)
mean(lengthinfo$ratio)

cortex22m <- lengthinfo
cortex22m$age <- "Old"

cortex <- rbind(cortex3m, cortex12m, cortex22m)
cortex <- data.frame(cortex)
cortex$ecc_number_per_MB <- 2^(cortex$log_ecc_number_per_MB)
cortex$norm_count_per_MB  <- 2^(cortex$log_norm_count_per_MB)

cortex$check <- log(cortex$ecc_number_per_MB,2)
cortex$check2 <- log(cortex$norm_count_per_MB,2)

cortex$ratio <- (cortex$ecc_number_per_MB/cortex$norm_count_per_MB)

cortex %>%
  group_by(age)
data.frame(cortex)
ratiog <- cortex[,c(14,15)]

ratiog %>%
  group_by(age) %>%
  summarise(across(.cols=everything(), list(median=median,sd=sd)))

data=ratiog
data$category <- factor(data$age, levels = c("Young", "Adult", "Old"))

ratiog$logratio <- log(ratiog$ratio,2)

my_comparisons <- list(c("Young", "Adult"), c("Adult", "Old"), c("Young", "Old"))
pairwise <- pairwise.wilcox.test(ratiog$logratio, ratiog$age, p.adjust.method = "bonferroni")

ggplot(data = ratiog) +
  aes(y = logratio, x = factor(age, level = c("Young", "Adult", "Old"))) +
  scale_fill_viridis_d(option = "D") +
  geom_point(aes(fill = age), shape = 21, size = 1.5, alpha = 0.2, position = "jitter", color = "black") +
  geom_violin(aes(fill = age), alpha = 0.75, size = 1, width = 1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label.y = 2) +
  ggtitle("Cortex") +
  labs(y = "log2(eccDNA_count/norm_count_per_MB)", x = "Age") +
  stat_summary(fun = "median", geom = "crossbar", width = 0.5, colour = "white") +
  geom_signif(comparisons = pairwise$by, y_position = 3)



##FOR HIPPOCAMPUS
## YOUNG DF FRAME
setwd("C:/Users/vzl597/Desktop/Xue Analysis")
hip3m <- read_excel("Mb that were found in the 3 month and 12 month liver, hippocampus and cortex.xlsx", 7)
lengthinfo <- left_join(hip3m, genes6, by = c("gene_name" = "gene_name_or_id"))
lengthinfo$slen <- lengthinfo$end - lengthinfo$start
lengthinfo$eccDNA <- round(exp(lengthinfo$log_ecc_number_per_MB ) * (lengthinfo$slen / 1000000))
lengthinfo$eccDNAperMB <- lengthinfo$eccDNA / (lengthinfo$slen/1000000)
lengthinfo$norm_count <- exp(lengthinfo$log_norm_count_per_MB)*(lengthinfo$slen / 1000000)
lengthinfo$ratio <- log(lengthinfo$eccDNAperMB / lengthinfo$norm_count,2)
median(lengthinfo$ratio)

hip3m <- lengthinfo
hip3m$age <- "Young"

hip12 <- read_excel("Mb that were found in the 3 month and 12 month liver, hippocampus and cortex.xlsx", 8)
lengthinfo <- left_join(hip12, genes6, by = c("gene_name" = "gene_name_or_id"))
lengthinfo$slen <- lengthinfo$end - lengthinfo$start
lengthinfo$eccDNA <- round(exp(lengthinfo$log_ecc_number_per_MB ) * (lengthinfo$slen / 1000000))
lengthinfo$eccDNAperMB <- lengthinfo$eccDNA / (lengthinfo$slen/1000000)
lengthinfo$norm_count <- exp(lengthinfo$log_norm_count_per_MB)*(lengthinfo$slen / 1000000)
lengthinfo$ratio <- log(lengthinfo$eccDNAperMB / lengthinfo$norm_count,2)
median(lengthinfo$ratio)

hip12m <- lengthinfo
hip12m$age <- "Adult"

hip22 <- read_excel("Mb that were found in the 3 month and 12 month liver, hippocampus and cortex.xlsx", 9)
lengthinfo <- left_join(hip22, genes6, by = c("gene_name" = "gene_name_or_id"))
lengthinfo$slen <- lengthinfo$end - lengthinfo$start
lengthinfo$eccDNA <- round(exp(lengthinfo$log_ecc_number_per_MB ) * (lengthinfo$slen / 1000000))
lengthinfo$eccDNAperMB <- lengthinfo$eccDNA / (lengthinfo$slen/1000000)
lengthinfo$norm_count <- exp(lengthinfo$log_norm_count_per_MB)*(lengthinfo$slen / 1000000)
lengthinfo$ratio <- log(lengthinfo$eccDNAperMB / lengthinfo$norm_count,2)
median(lengthinfo$ratio)

hip22m <- lengthinfo
hip22m$age <- "Old"

hipp <- rbind(hip3m, hip12m, hip22m)
hipp <- data.frame(hipp)
hipp$ecc_number_per_MB <- 2^(hipp$log_ecc_number_per_MB)
hipp$norm_count_per_MB  <- 2^(hipp$log_norm_count_per_MB)

hipp$check <- log(hipp$ecc_number_per_MB,2)
hipp$check2 <- log(hipp$norm_count_per_MB,2)

hipp$ratio <- (hipp$ecc_number_per_MB/hipp$norm_count_per_MB)

hipp %>%
  group_by(age)
data.frame(hipp)
ratiog <- hipp[,c(14,15)]
ratiog$logratio <- log(ratiog$ratio,2)
ratiog %>%
  group_by(age) %>%
  summarise(across(.cols=everything(), list(median=median,sd=sd)))

data=ratiog
data$category <- factor(data$age, levels = c("Young", "Adult", "Old"))


my_comparisons <- list(c("Young", "Adult"), c("Adult", "Old"), c("Young", "Old"))
pairwise <- pairwise.wilcox.test(ratiog$logratio, ratiog$age, p.adjust.method = "bonferroni")

ggplot(data = ratiog) +
  aes(y = logratio, x = factor(age, level = c("Young", "Adult", "Old"))) +
  scale_fill_viridis_d(option = "D") +
  geom_point(aes(fill = age), shape = 21, size = 1.5, alpha = 0.2, position = "jitter", color = "black") +
  geom_violin(aes(fill = age), alpha = 0.75, size = 1, width = 1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label.y = 2) +
  ggtitle("Hippocampus") +
  labs(y = "log2(eccDNA_count/norm_count_per_MB)", x = "Age") +
  stat_summary(fun = "median", geom = "crossbar", width = 0.5, colour = "white") +
  geom_signif(comparisons = pairwise$by, y_position = 3)
