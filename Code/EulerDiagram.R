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
library(readxl)
library(purrr)
library(tidyverse)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ensembldb)
library(GenomicFeatures)
library(GenomicRanges)
library(biomaRt)
library(eulerr)
install.packages("eulerr")
## FOR CORTEX
setwd("C:/Users/vzl597/Desktop/Xue Analysis")
Young <- read_excel("Mb that were found in the 3 month and 12 month liver, hippocampus and cortex.xlsx", 4)
Adult <- read_excel("Mb that were found in the 3 month and 12 month liver, hippocampus and cortex.xlsx", 5)
Old <- read_excel("Mb that were found in the 3 month and 12 month liver, hippocampus and cortex.xlsx", 6)

Young$Age <- "Young"
Adult$Age <- "Adult"
Old$Age <- "Old"

AllData <- rbind(Young, Adult, Old)
above_genes <- AllData[AllData$label == "above",]
below_genes <- AllData[AllData$label == "below",]

above_genes_list <- split(above_genes$gene_name, above_genes$Age)
below_genes_list <- split(below_genes$gene_name, below_genes$Age)


# Make sure if there are duplicates and if so, remove them. 
above_genes_list$Young[duplicated(above_genes_list$Young)]
above_genes_list$Adult[duplicated(above_genes_list$Adult)]
above_genes_list$Old[duplicated(above_genes_list$Old)]

a1 <- above_genes_list$Young[!duplicated(above_genes_list$Young)]
above_genes_list[[3]] <- a1

a1 <- above_genes_list$Adult[!duplicated(above_genes_list$Adult)]
above_genes_list[[1]] <- a1

a1 <- above_genes_list$Old[!duplicated(above_genes_list$Old)]
above_genes_list[[2]] <- a1

# Make sure if there are duplicates and if so, remove them. 
below_genes_list$Young[duplicated(below_genes_list$Young)]
below_genes_list$Adult[duplicated(below_genes_list$Adult)]
below_genes_list$Old[duplicated(below_genes_list$Old)]

dev.new()

# Calcualte, print and plot euler plots
fit <- euler(above_genes_list, shape = "ellipse")
print(fit)
plot(fit)

# Calcualte, print and plot euler plots
fit <- euler(below_genes_list, shape = "ellipse")
print(fit)
plot(fit)

n_permutations <- 10000

# Example data frame of genes
genes_df <- AllData[["gene_name"]]

# Define group lengths
group_lengths <- c(138, 522, 438)

# Set number of random samples to generate
num_samples <- 1000

# Initialize empty vector to store overlap counts
overlap_counts <- rep(0, num_samples)

# Loop through each sample
for (i in 1:num_samples) {
  
  # Sample genes from data frame without replacement
  sampled_genes <- sample(genes_df, sum(group_lengths), replace = FALSE)
  
  # Split genes into three groups of specified lengths
  gene_groups <- split(sampled_genes, rep(1:3, group_lengths))
  
  # Calculate overlap between groups
  overlap_counts[i] <- length(Reduce(intersect, gene_groups))
}

# Calculate mean and standard deviation of overlap counts
mean_overlap <- mean(overlap_counts)
sd_overlap <- sd(overlap_counts)

# Print results
cat("Mean overlap:", mean_overlap, "\n")
cat("Standard deviation of overlap:", sd_overlap, "\n")



AllData

fit2 <- euler(below_genes_list, shape = "ellipse", quantities = TRUE)
print(fit2)
plot(fit2)

obs_overlap <- fit2$areas[7]

## FOR HIPPOCAMPUS
setwd("C:/Users/vzl597/Desktop/Xue Analysis")
Young <- read_excel("Mb that were found in the 3 month and 12 month liver, hippocampus and cortex.xlsx", 7)
Adult <- read_excel("Mb that were found in the 3 month and 12 month liver, hippocampus and cortex.xlsx", 8)
Old <- read_excel("Mb that were found in the 3 month and 12 month liver, hippocampus and cortex.xlsx", 9)

Young$Age <- "Young"
Adult$Age <- "Adult"
Old$Age <- "Old"

AllData <- rbind(Young, Adult, Old)
above_genes <- AllData[AllData$label == "above",]
below_genes <- AllData[AllData$label == "below",]

above_genes_list <- split(above_genes$gene_name, above_genes$Age)
below_genes_list <- split(below_genes$gene_name, below_genes$Age)


# Make sure if there are duplicates and if so, remove them. 
above_genes_list$Young[duplicated(above_genes_list$Young)]
above_genes_list$Adult[duplicated(above_genes_list$Adult)]
above_genes_list$Old[duplicated(above_genes_list$Old)]

a1 <- above_genes_list$Young[!duplicated(above_genes_list$Young)]
above_genes_list[[3]] <- a1

a1 <- above_genes_list$Adult[!duplicated(above_genes_list$Adult)]
above_genes_list[[1]] <- a1

# Make sure if there are duplicates and if so, remove them. 
below_genes_list$Young[duplicated(below_genes_list$Young)]
below_genes_list$Adult[duplicated(below_genes_list$Adult)]
below_genes_list$Old[duplicated(below_genes_list$Old)]

dev.new()

# Calcualte, print and plot euler plots
fit <- euler(above_genes_list, shape = "ellipse")
print(fit)
plot(fit)

# Calcualte, print and plot euler plots
fit <- euler(below_genes_list, shape = "ellipse")
print(fit)
plot(fit)

n_permutations <- 10000

# Example data frame of genes
genes_df <- AllData[["gene_name"]]

# Define group lengths
group_lengths <- c(138, 522, 438)

# Set number of random samples to generate
num_samples <- 1000

# Initialize empty vector to store overlap counts
overlap_counts <- rep(0, num_samples)

# Loop through each sample
for (i in 1:num_samples) {
  
  # Sample genes from data frame without replacement
  sampled_genes <- sample(genes_df, sum(group_lengths), replace = FALSE)
  
  # Split genes into three groups of specified lengths
  gene_groups <- split(sampled_genes, rep(1:3, group_lengths))
  
  # Calculate overlap between groups
  overlap_counts[i] <- length(Reduce(intersect, gene_groups))
}

# Calculate mean and standard deviation of overlap counts
mean_overlap <- mean(overlap_counts)
sd_overlap <- sd(overlap_counts)

# Print results
cat("Mean overlap:", mean_overlap, "\n")
cat("Standard deviation of overlap:", sd_overlap, "\n")



AllData

fit2 <- euler(below_genes_list, shape = "ellipse", quantities = TRUE)
print(fit2)
plot(fit2)

obs_overlap <- fit2$areas[7]

