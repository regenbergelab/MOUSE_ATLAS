library(GenomicRanges)

# Assume your data is stored in a list called "data" 
gr_list <- convert_to_GRanges(data)

# Reduce the GRanges objects, so that each range is unique. The final object is the eccDNA genomic coverage 
reduced_gr_list <- lapply(gr_list, GenomicRanges::reduce)

# Calculate the total length of each GRanges object
sum_lengths <- sapply(reduced_gr_list, function(x) sum(width(x)))

# Create a matrix of the sum of lengths between all GRanges objects 
sum_lengths_matrix <- matrix(nrow = length(sum_lengths), ncol = length(sum_lengths)) 
for (i in seq_along(sum_lengths)) {
for (j in seq_along(sum_lengths)) {
sum_lengths_matrix[i,j] <- sum_lengths[i] + sum_lengths[j] 
}
}

# Calculate the overlap matrix. Introduce min overlap between two samples is 1 bp. 
overlap_lengths <- matrix(0, nrow = length(reduced_gr_list), ncol = length(reduced_gr_list)) 
for (i in seq_along(reduced_gr_list)) {
for (j in seq_along(reduced_gr_list)) {
overlaps <- findOverlaps(reduced_gr_list[[i]], reduced_gr_list[[j]], type = "any")
overlap_lengths[i, j] <- sum(width(pintersect(reduced_gr_list[[i]][queryHits(overlaps)], reduced_gr_list[[j]][subjectHits(overlaps)]))) 
# replace 0s with 1s
overlap_lengths[i, j] <- ifelse(overlap_lengths[i, j] == 0, 1, overlap_lengths[i, j])
} 
}
overlap_lengths

# Normalize the overlap matrix
overlap_lengths_norm <- 2*overlap_lengths / sum_lengths_matrix

# Print the normalized overlap matrix 
overlap_lengths_norm
