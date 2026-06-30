# ============================================================
# DAPC - BIC optimization to determine the optimal number of clusters
# ============================================================

library(adegenet)
library(tidyverse)
library(readxl)
library(SNPRelate)

# ============================================================
# Input files
# ============================================================

# Metadata containing samples to exclude
metadata_flag <- read_excel("path/to/metadata_file.xlsx")

# Samples flagged for exclusion
flagged_samples <- metadata_flag %>%
  filter(`Flag-weirdposition` == "Y") %>%
  pull(sample_id_column)


# VCF input file
vcf_file <- "path/to/input_file.vcf"

# Temporary GDS file generated from VCF
gds_file <- "temporary_file.gds"


# Population metadata
population_metadata <- read_tsv("path/to/population_metadata.tsv")


# Extract sample identifiers
population_metadata$sample_short <- sapply(
  strsplit(sub("\\..*", "", population_metadata$sample.id), "_"),
  `[`,
  3
)

population_metadata$sample_short <- gsub("\\D", "", population_metadata$sample_short)


# Remove flagged samples
population_metadata <- population_metadata %>%
  filter(!(sample_short %in% flagged_samples))


# ============================================================
# Convert VCF to GDS format and extract genotype matrix
# ============================================================

snpgdsVCF2GDS(
  vcf_file,
  gds_file,
  method = "biallelic.only"
)

genofile <- snpgdsOpen(gds_file)


# Retrieve SNP and sample identifiers
snp_ids <- read.gdsn(index.gdsn(genofile, "snp.id"))
sample_ids <- read.gdsn(index.gdsn(genofile, "sample.id"))


# Keep only samples present in metadata
pattern <- paste0("^", population_metadata$sample.id, collapse = "|")

selected_samples <- sample_ids[
  grepl(pattern, sample_ids)
]


# Extract genotype matrix
genotype_matrix <- snpgdsGetGeno(
  genofile,
  sample.id = selected_samples
)


rownames(genotype_matrix) <- selected_samples
colnames(genotype_matrix) <- snp_ids


# ============================================================
# Convert genotype matrix to genind object
# ============================================================

genind_object <- df2genind(
  genotype_matrix,
  ploidy = 2,
  type = "codom",
  ind.names = selected_samples,
  loc.names = snp_ids,
  ncode = 1
)


# ============================================================
# DAPC clustering and BIC optimization
# ============================================================

pc_values <- seq(50, 250, by = 1)   # Number of PCs tested
n_replicates <- 20                  # Replicates per PC value
max_clusters <- 10                  # Maximum number of clusters tested


# Store results
bic_results <- data.frame(
  PC = integer(),
  Replicate = integer(),
  Best_K = integer()
)


set.seed(420)


for (pc in pc_values) {

  cat("Running analysis with", pc, "PCs\n")

  for (rep in 1:n_replicates) {

    set.seed(1000 + pc * 100 + rep)

    clustering <- find.clusters(
      genind_object,
      max.n.clust = max_clusters,
      n.pca = pc,
      choose.n.clust = FALSE
    )


    # Cluster number with minimum BIC
    optimal_k <- which.min(clustering$Kstat)


    bic_results <- rbind(
      bic_results,
      data.frame(
        PC = pc,
        Replicate = rep,
        Best_K = optimal_k
      )
    )
  }
}


# ============================================================
# Summary of optimal K values
# ============================================================

optimal_k_frequency <- table(bic_results$Best_K)


barplot(
  optimal_k_frequency,
  main = "Frequency of optimal K values (minimum BIC)",
  xlab = "Number of clusters (K)",
  ylab = "Number of occurrences"
)
