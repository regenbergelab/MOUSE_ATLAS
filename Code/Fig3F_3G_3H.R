# Load necessary libraries
library(dplyr)
library(ggplot2)
library(rtracklayer)
library(GenomicFeatures)
library(ggpubr)
library(rstatix)

# Load gene annotations from GTF file
gtf <- rtracklayer::import("Mus_musculus.GRCm38.100.gtf.gz")
gtf_df <- as.data.frame(gtf)

# Filter genes and remove duplicates
genes <- gtf_df[gtf_df$type == "gene", ]
genes <- genes[!duplicated(genes$gene_name), ]
genes <- genes[!is.na(genes$gene_name), ]
genes <- genes[, c("seqnames", "start", "end", "gene_name")]
colnames(genes) <- c("chr", "start", "end", "Gene")
genes$slen <- genes$end - genes$start

# Calculate exon and transcript information
txdb <- makeTxDbFromGFF("Mus_musculus.GRCm38.100.gtf.gz")
exons_per_gene <- exonsBy(txdb, by = "gene")
transcripts_per_gene <- as.data.frame(transcriptsBy(txdb, by = "gene"))

# Summarize exon and transcript counts
exon_counts <- as.data.frame(width(reduce(exons_per_gene))) %>%
  group_by(group_name) %>%
  summarise(exon_length = sum(value))

transcript_counts <- transcripts_per_gene %>%
  group_by(group_name) %>%
  summarise(transcript_number = n())

# Merge gene information with exon and transcript counts
gene_info <- left_join(genes, transcript_counts, by = c("Gene" = "group_name"))
gene_info <- left_join(gene_info, exon_counts, by = c("Gene" = "group_name"))
gene_info$intronic_length <- gene_info$slen - gene_info$exon_length
gene_info$exon_density <- gene_info$exon_length * 1000 / gene_info$intronic_length
gene_info$intron_density <- gene_info$intronic_length / gene_info$slen


# Load gene expression data for different age groups
load_gene_data <- function(age, tissue) {
  file_path <- paste0("./", age, "_", tissue, "/10%_quantile_genelist_", age, "_", tissue, ".csv")
  data <- read.csv(file_path)
  data <- left_join(data, gene_info, by = c("gene_name" = "Gene"))
  data$Age <- age
  return(data)
}

young_data <- load_gene_data("Young", "Liver")
adult_data <- load_gene_data("Adult", "Liver")
old_data <- load_gene_data("Old", "Liver")

# Combine data into a single dataframe
all_data <- bind_rows(
  young_data %>% filter(label == "below") %>% mutate(eccDNAlevel = "Depleted"),
  adult_data %>% filter(label == "below") %>% mutate(eccDNAlevel = "Depleted"),
  old_data %>% filter(label == "below") %>% mutate(eccDNAlevel = "Depleted"),
  young_data %>% filter(label == "middle") %>% mutate(eccDNAlevel = "Bulk"),
  adult_data %>% filter(label == "middle") %>% mutate(eccDNAlevel = "Bulk"),
  old_data %>% filter(label == "middle") %>% mutate(eccDNAlevel = "Bulk")
)

# Define enriched genes
common_genes <- c("Arhgef3", "Camta1", "Chcdh3", "Clybl", "Dmd", "Frmd4b", "Mast4", "Pcca", "Pcx", "Tbc1d5", "Zfand3")
common_genes_info <- subset(gene_info, Gene %in% common_genes)
common_genes_info$eccDNAlevel <- "Recurrent"
all_data <- bind_rows(all_data, common_genes_info)

# Calculate exon density and intron density
all_data$exons_per_mb <- all_data$MaxExonNumber * 1e6 / all_data$slen
all_data$exon_density <- all_data$exon_length * 1000 / all_data$intronic_length
all_data$intron_density <- all_data$intronic_length / all_data$slen

# Reorder age levels
all_data$Age <- factor(all_data$Age, levels = c("Young", "Adult", "Old"))

# ----------------------------- Fig3f: Intron Density -----------------------------
ggplot(all_data, aes(x = eccDNAlevel, y = intron_density, fill = eccDNAlevel)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("#999999", "#a7dcdd", "#3e7b7d")) +
  ylim(0.5, 1.1) +
  labs(x = "", y = "Total Intron Length per Gene Length") +
  facet_grid(. ~ Age) +
  theme_classic()
ggsave("Fig3f.pdf")

# ----------------------------- Fig3g: Transcript Isoforms -----------------------------
ggplot(all_data, aes(x = eccDNAlevel, y = log10(transcript_number), fill = eccDNAlevel)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("#999999", "#a7dcdd", "#3e7b7d")) +
  ylim(0, 1.5) +
  labs(y = "Number of Transcript Isoforms (log10)") +
  facet_grid(. ~ Age) +
  theme_classic()
ggsave("Fig3g.pdf")

# ----------------------------- Fig3h: Isoforms per Exon Number per Mb -----------------------------
all_data$isoforms_per_exon_per_mb <- all_data$transcript_number / (all_data$exon_length / 1e6)

ggplot(all_data, aes(x = eccDNAlevel, y = isoforms_per_exon_per_mb, fill = eccDNAlevel)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("#999999", "#a7dcdd", "#3e7b7d")) +
  ylim(0, 0.75) +
  labs(y = "Isoforms per Exon Number per Mb") +
  facet_grid(. ~ Age) +
  theme_classic()
ggsave("Fig3h.pdf")
