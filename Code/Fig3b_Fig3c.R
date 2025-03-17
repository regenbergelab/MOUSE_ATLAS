# Load necessary libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(eulerr)  # For Euler diagrams
library(readr)   # For reading CSV files

# Load gene lists for different age groups
young_genes <- read.csv("10%_quantile_genelist_3m.csv")
adult_genes <- read.csv("10%_quantile_genelist_12m.csv")
old_genes <- read.csv("10%_quantile_genelist_22m.csv")

# Add age labels to each dataset
young_genes$Age <- "Young"
adult_genes$Age <- "Adult"
old_genes$Age <- "Old"

# Combine all data into a single dataframe
all_genes <- rbind(young_genes, adult_genes, old_genes)

# Split genes into "above" and "below" groups based on label
above_genes <- all_genes[all_genes$label == "above", ]
below_genes <- all_genes[all_genes$label == "below", ]

# Create lists of genes for each age group
above_genes_list <- split(above_genes$gene_name, above_genes$Age)
below_genes_list <- split(below_genes$gene_name, below_genes$Age)

# Remove duplicates from each gene list
above_genes_list <- lapply(above_genes_list, unique)
below_genes_list <- lapply(below_genes_list, unique)

# ----------------------------- Fig3b: Euler Plot for Below Genes -----------------------------
# Calculate Euler diagram for "below" genes
euler_below <- euler(below_genes_list, shape = "ellipse")

pdf("Fig3b.pdf", height = 5, width = 5)
plot(euler_below, 
     fills = list(fill = c("#deebf7", "#deebf7", "#deebf7", "#6baed6", "#6baed6", "#6baed6", "#08519c"), alpha = 0.5),
     quantities = c(46, 182, 282, 29, 27, 79, 16))
dev.off()

# ----------------------------- Fig3c: Euler Plot for Above Genes -----------------------------
# Calculate Euler diagram for "above" genes
euler_above <- euler(above_genes_list, shape = "ellipse")

pdf("Fig3c.pdf", height = 5, width = 5)
plot(euler_above, 
     fills = list(fill = c("#feedde", "#feedde", "#feedde", "#fd8d3c", "#fd8d3c", "#fd8d3c", "#a63603"), alpha = 0.5),
     quantities = c(107, 285, 383, 4, 6, 14, 0))
dev.off()
