## Reviewer requested summary of LDFP ref library / blast results of ASVs before control using 
## NTC samples and BLank DNA extract control samples to assess variaiton in ASVs across samples and controls

## firstly consolidating information in a workbook for supplementary material
precleanp1 <- p1data1.pr$p1pspool.lulu
precleanp2 <- p2data1.pr$p2pspool.lulu
str(precleanp1)
str(precleanp2)

sample_dataPNAS <- read.csv("Processed_data/sample_dataPNAS.csv")
str(sample_dataPNAS)

## need to remove samples that do not appear in this study
valid_names <- sample_dataPNAS$name
rownames_p1 <- rownames(precleanp1)
base_names_p1 <- sub("r[0-9]+$", "", rownames_p1)
keep_rows <- !grepl("^S_", rownames_p1) | base_names_p1 %in% valid_names
precleanp1_filtered <- precleanp1[keep_rows, ]

rownames_p2 <- rownames(precleanp2)
base_names_p2 <- sub("r[0-9]+$", "", rownames_p2)
keep_rows_p2 <- !grepl("^S_", rownames_p2) | base_names_p2 %in% valid_names
precleanp2_filtered <- precleanp2[keep_rows_p2, ]

precleanp1_filtered1 <- precleanp1_filtered[, colSums(precleanp1_filtered) > 0, drop = FALSE]
precleanp2_filtered1 <- precleanp2_filtered[, colSums(precleanp2_filtered) > 0, drop = FALSE]

p1_sequences <- colnames(precleanp1_filtered1)
p2_sequences <- colnames(precleanp2_filtered1)

colnames(precleanp1_filtered1) <- paste0("P1_ASV", seq_along(p1_sequences))
colnames(precleanp2_filtered1) <- paste0("P2_ASV", seq_along(p2_sequences))

write.csv(precleanp1_filtered1, "Processed_data/plateonePrefiltp.csv")
write.csv(precleanp2_filtered1, "Processed_data/platetwoPrefiltp.csv")

names(p1_sequences) <- colnames(precleanp1_filtered1)
names(p2_sequences) <- colnames(precleanp2_filtered1)

all_sequences <- union(p1_sequences, p2_sequences)

seq_mapping <- data.frame(
  P1_ASV_id = ifelse(all_sequences %in% p1_sequences, names(p1_sequences)[match(all_sequences, p1_sequences)], ""),
  P2_ASV_id = ifelse(all_sequences %in% p2_sequences, names(p2_sequences)[match(all_sequences, p2_sequences)], ""),
  Sequence  = all_sequences,
  stringsAsFactors = FALSE
)

df <- R1all.taxs.source.otu$dada.pspool.nc.lulu
source_with_seqs <- R1SppUnAssList$dada.pspool.nc.lulu
all_sources <- rbind(
  R1SppAssList$dada.pspool.nc.lulu,
  R1SppUnAssList$dada.pspool.nc.lulu
)

otu_to_seq <- setNames(rownames(all_sources), all_sources[, "OTU"])
df$sequence <- otu_to_seq[rownames(df)]

df$downstream_ASV_ID <- rownames(df)
joined_df <- merge(seq_mapping, df, by.x = "Sequence", by.y = "sequence", all.x = TRUE)
rownames(joined_df) <- joined_df$Original_OTU_ID
joined_df$OTU_num <- as.numeric(sub("OTU", "", joined_df$downstream_ASV_ID))
joined_df_sorted <- joined_df[order(joined_df$OTU_num), ]
joined_df_sorted$OTU_num <- NULL

library(Biostrings)
library(pwalign)

seqs <- DNAStringSet(joined_df_sorted$Sequence)
names(seqs) <- joined_df_sorted$downstream_ASV_ID

n <- length(seqs)
dna_dist <- matrix(NA, n, n, dimnames = list(names(seqs), names(seqs)))
submat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)

for (i in 1:n) {
  for (j in i:n) {
    aln <- pairwiseAlignment(seqs[[i]], seqs[[j]],
                             substitutionMatrix = submat,
                             gapOpening = 10,
                             gapExtension = 0.5)
    score_val <- -score(aln)
    dna_dist[i, j] <- score_val
    dna_dist[j, i] <- score_val
  }
}

hc <- hclust(as.dist(dna_dist), method = "average")
joined_df_sorted_similarity <- joined_df_sorted[hc$order, ]

seqs <- DNAStringSet(joined_df_sorted_similarity$Sequence)
names(seqs) <- joined_df_sorted_similarity$downstream_ASV_ID

n <- length(seqs)
seq_names <- names(seqs)
submat <- nucleotideSubstitutionMatrix(match = 0, mismatch = 1, baseOnly = TRUE)
gapOpen <- 10
gapExt <- 0.5

pairs <- combn(n, 2, simplify = FALSE)

library(parallel)
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
clusterExport(cl, c("seqs", "submat", "gapOpen", "gapExt"), envir = environment())
clusterEvalQ(cl, library(Biostrings))

dist_list <- parLapply(cl, pairs, function(idx) {
  i <- idx[1]
  j <- idx[2]
  aln <- Biostrings::pairwiseAlignment(seqs[[i]], seqs[[j]],
                                       substitutionMatrix = submat,
                                       gapOpening = gapOpen,
                                       gapExtension = gapExt)
  mismatches <- score(aln)
  list(i = i, j = j, dist = mismatches)
})

stopCluster(cl)

dna_dist <- matrix(0, n, n, dimnames = list(seq_names, seq_names))
for (res in dist_list) {
  i <- res$i
  j <- res$j
  d <- res$dist
  dna_dist[i, j] <- d
  dna_dist[j, i] <- d
}

otu_ids <- joined_df_sorted_similarity$downstream_ASV_ID
first_match <- second_match <- character(length(otu_ids))
first_diff <- second_diff <- numeric(length(otu_ids))

for (i in seq_along(otu_ids)) {
  this_otu <- otu_ids[i]
  distances <- dna_dist[this_otu, ]
  distances[this_otu] <- Inf
  closest <- order(abs(distances), decreasing = FALSE)[1:2]
  first_match[i] <- names(distances)[closest[1]]
  first_diff[i]  <- distances[closest[1]]
  second_match[i] <- names(distances)[closest[2]]
  second_diff[i]  <- distances[closest[2]]
}

joined_df_sorted_similarity$`1stClosestMatch` <- first_match
joined_df_sorted_similarity$`1stClosestDiff`  <- first_diff
joined_df_sorted_similarity$`2ndClosestMatch` <- second_match
joined_df_sorted_similarity$`2ndClosestDiff`  <- second_diff

count_controls <- function(otu_counts, sample_names) {
  present_samples <- sample_names[otu_counts > 0]
  pcr_N <- sum(grepl("^N_", present_samples))
  pcr_E <- sum(grepl("^E_", present_samples))
  base_names <- sub("r[0-9]+$", "", present_samples)
  sample_N <- length(unique(base_names[grepl("^N_", base_names)]))
  sample_E <- length(unique(base_names[grepl("^E_", base_names)]))
  return(c(pcr_N = pcr_N, sample_N = sample_N, pcr_E = pcr_E, sample_E = sample_E))
}

joined_df_sorted_similarity$pcr_N <- 0
joined_df_sorted_similarity$sample_N <- 0
joined_df_sorted_similarity$pcr_E <- 0
joined_df_sorted_similarity$sample_E <- 0

for (i in seq_len(nrow(joined_df_sorted_similarity))) {
  row <- joined_df_sorted_similarity[i, ]
  if (nzchar(row$P1_ASV_id) && row$P1_ASV_id %in% colnames(precleanp1_filtered1)) {
    counts <- count_controls(precleanp1_filtered1[, row$P1_ASV_id], rownames(precleanp1_filtered1))
    joined_df_sorted_similarity[i, c("pcr_N", "sample_N", "pcr_E", "sample_E")] <-
      joined_df_sorted_similarity[i, c("pcr_N", "sample_N", "pcr_E", "sample_E")] + counts
  }
  if (nzchar(row$P2_ASV_id) && row$P2_ASV_id %in% colnames(precleanp2_filtered1)) {
    counts <- count_controls(precleanp2_filtered1[, row$P2_ASV_id], rownames(precleanp2_filtered1))
    joined_df_sorted_similarity[i, c("pcr_N", "sample_N", "pcr_E", "sample_E")] <-
      joined_df_sorted_similarity[i, c("pcr_N", "sample_N", "pcr_E", "sample_E")] + counts
  }
}

joined_df_sorted_similarity$pcr_Nprop <- 0
joined_df_sorted_similarity$pcr_Eprop <- 0
joined_df_sorted_similarity$pcr_Nmaxcount <- 0
joined_df_sorted_similarity$pcr_Emaxcount <- 0

for (i in seq_len(nrow(joined_df_sorted_similarity))) {
  row <- joined_df_sorted_similarity[i, ]
  pcr_N_counts <- c()
  pcr_E_counts <- c()
  
  if (nzchar(row$P1_ASV_id) && row$P1_ASV_id %in% colnames(precleanp1_filtered1)) {
    counts <- count_controls(precleanp1_filtered1[, row$P1_ASV_id], rownames(precleanp1_filtered1))
    joined_df_sorted_similarity[i, c("pcr_N", "sample_N", "pcr_E", "sample_E")] <-
      joined_df_sorted_similarity[i, c("pcr_N", "sample_N", "pcr_E", "sample_E")] + counts
    pcr_N_counts <- c(pcr_N_counts, precleanp1_filtered1[grepl("^N_", rownames(precleanp1_filtered1)), row$P1_ASV_id])
    pcr_E_counts <- c(pcr_E_counts, precleanp1_filtered1[grepl("^E_", rownames(precleanp1_filtered1)), row$P1_ASV_id])
  }
  
  if (nzchar(row$P2_ASV_id) && row$P2_ASV_id %in% colnames(precleanp2_filtered1)) {
    counts <- count_controls(precleanp2_filtered1[, row$P2_ASV_id], rownames(precleanp2_filtered1))
    joined_df_sorted_similarity[i, c("pcr_N", "sample_N", "pcr_E", "sample_E")] <-
      joined_df_sorted_similarity[i, c("pcr_N", "sample_N", "pcr_E", "sample_E")] + counts
    pcr_N_counts <- c(pcr_N_counts, precleanp2_filtered1[grepl("^N_", rownames(precleanp2_filtered1)), row$P2_ASV_id])
    pcr_E_counts <- c(pcr_E_counts, precleanp2_filtered1[grepl("^E_", rownames(precleanp2_filtered1)), row$P2_ASV_id])
  }
  
  if (length(pcr_N_counts) > 0) {
    joined_df_sorted_similarity$pcr_Nprop[i] <- sum(pcr_N_counts > 0) / length(pcr_N_counts)
    joined_df_sorted_similarity$pcr_Nmaxcount[i] <- max(pcr_N_counts, na.rm = TRUE)
  }
  if (length(pcr_E_counts) > 0) {
    joined_df_sorted_similarity$pcr_Eprop[i] <- sum(pcr_E_counts > 0) / length(pcr_E_counts)
    joined_df_sorted_similarity$pcr_Emaxcount[i] <- max(pcr_E_counts, na.rm = TRUE)
  }
}

write.csv(joined_df_sorted_similarity, "Processed_data/summaryE&Nsamplesp.csv")
dna_dist[upper.tri(dna_dist)] <- NA
write.csv(dna_dist, file = "Processed_data/dna_distance_matrix.csv", quote = FALSE, na = "")

