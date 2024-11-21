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

gtf <- rtracklayer::import(paste0("Mus_musculus.GRCm38.100.gtf.gz"))
gtf_df<-as.data.frame(gtf)
gtf_df$gene_biotype %>% unique
gtf_df$gene_name %>% unique %>% length
gtf_df["gene_name_or_id"]<-gtf$gene_name
gtf_df$gene_name_or_id[is.na(gtf_df$gene_name_or_id)]<-gtf$gene_id[is.na(gtf_df["gene_name_or_id"])]
genes<-gtf_df[gtf_df$type=="gene",]
genes$gene_name_or_id[genes$gene_name_or_id %>% duplicated %>% which] %>% table
genes2<-genes[!duplicated(genes$gene_name_or_id),]
genes3<-genes2[!is.na(genes$gene_name),]
genes4 <- genes3[,-(4:27),drop=FALSE]
genes5 <- genes4[,-c(4,5,6,7)]
genes6 <- na.omit(genes5)
names(genes6)[names(genes6) == "seqnames"] <- "chr"
lengthinfo <- left_join(liver3m, genes6, by = c("gene_name" = "gene_name_or_id"))
lengthinfo$slen <- lengthinfo$end - lengthinfo$start
lengthinfo$eccDNA <- round(exp(lengthinfo$log_ecc_number_per_MB ) * (lengthinfo$slen / 1000000))
lengthinfo$eccDNAperMB <- lengthinfo$eccDNA / (lengthinfo$slen/1000000)
lengthinfo$norm_count <- exp(lengthinfo$log_norm_count_per_MB)*(lengthinfo$slen / 1000000)
lengthinfo$ratio <- log(lengthinfo$eccDNAperMB / lengthinfo$norm_count,2)
median(lengthinfo$ratio)
liver3m <- lengthinfo
liver3m$age <- "Young"
lengthinfo <- left_join(liver12m, genes6, by = c("gene_name" = "gene_name_or_id"))
lengthinfo$slen <- lengthinfo$end - lengthinfo$start
lengthinfo$eccDNA <- round(exp(lengthinfo$log_ecc_number_per_MB ) * (lengthinfo$slen / 1000000))
lengthinfo$eccDNAperMB <- lengthinfo$eccDNA / (lengthinfo$slen/1000000)
lengthinfo$norm_count <- exp(lengthinfo$log_norm_count_per_MB)*(lengthinfo$slen / 1000000)
lengthinfo$ratio <- log(lengthinfo$eccDNAperMB / lengthinfo$norm_count,2)
median(lengthinfo$ratio)
liver12m <- lengthinfo
liver12m$age <- "Adult"
lengthinfo <- left_join(liver22m, genes6, by = c("gene_name" = "gene_name_or_id"))
lengthinfo$slen <- lengthinfo$end - lengthinfo$start
lengthinfo$eccDNA <- round(exp(lengthinfo$log_ecc_number_per_MB ) * (lengthinfo$slen / 1000000))
lengthinfo$eccDNAperMB <- lengthinfo$eccDNA / (lengthinfo$slen/1000000)
lengthinfo$norm_count <- exp(lengthinfo$log_norm_count_per_MB)*(lengthinfo$slen / 1000000)
lengthinfo$ratio <- log(lengthinfo$eccDNAperMB / lengthinfo$norm_count,2)
median(lengthinfo$ratio)
liver22m <- lengthinfo
liver22m$age <- "Old"
liver <- rbind(Young, Adult, Old)
liver <- data.frame(liver)
liver$ecc_number_per_MB <- 2^(liver$log_ecc_number_per_MB)
liver$norm_count_per_MB  <- 2^(liver$log_norm_count_per_MB)
liver$check <- log(liver$ecc_number_per_MB,2)
liver$check2 <- log(liver$norm_count_per_MB,2)
liver$ratio <- (liver$ecc_number_per_MB/liver$norm_count_per_MB)
liver %>% group_by(age)
ratiog <- liver[,c(14,15)]
ratiog$logratio <- log(ratiog$ratio,2)
ratiog %>% group_by(age) %>% summarise(across(.cols=everything(), list(median=median,sd=sd)))
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
