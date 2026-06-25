# =====================================================
# LIBRARIES
# =====================================================
library(tidyverse)
library(SNPRelate)
library(hierfstat)
library(adegenet)
library(ape)

# =====================================================
# PATHS (PROJECT-RELATIVE)
# =====================================================
base_dir <- "PROJECT_DIRECTORY"

data_dir <- file.path(base_dir, "data")
vcf_dir <- file.path(data_dir, "vcf")
meta_dir <- file.path(data_dir, "metadata")
cluster_dir <- file.path(data_dir, "clusters")

fig_dir <- file.path(base_dir, "results/figures")

# =====================================================
# LOAD METADATA
# =====================================================
metadata <- read_tsv(file.path(meta_dir,
                               "litterature_metadata_STRUCTURE_DAPC_results.tsv"))

metadataK4 <- read_tsv(file.path(meta_dir,
                                 "3D_metadata_STRUCTURE_DAPC_results.tsv"))

metadata_pop <- read_tsv(file.path(meta_dir,
                                   "PCA_litterature_withpopulationmetadata.tsv"))

# Update population labels using K4 clustering
metadata_pop$population[
  metadata_pop$sample.id %in% metadataK4$sample_id
] <- metadataK4$K4_cluster_majoritaire[
  match(metadata_pop$sample.id[metadata_pop$sample.id %in% metadataK4$sample_id],
        metadataK4$sample_id)
]

metadata <- left_join(metadata,
                      metadata_pop[, c(1, 4, 5)],
                      by = c("sample_id" = "sample.id"))

# =====================================================
# VCF → GDS
# =====================================================
vcf.fn <- file.path(vcf_dir, "dataset.vcf")
gds.fn <- file.path(vcf_dir, "temp.gds")

snpgdsVCF2GDS(vcf.fn, gds.fn, method = "biallelic.only")
genofile <- snpgdsOpen(gds.fn)

sample <- read.gdsn(index.gdsn(genofile, "sample.id"))
snps.ids <- read.gdsn(index.gdsn(genofile, "snp.id"))

# Standardize sample IDs
sample_new <- ifelse(
  grepl("^GBS", sample),
  gsub("\\D", "", sapply(strsplit(sub("\\..*", "", sample), "_"), `[`, 3)),
  sample
)

keep <- sample_new %in% metadata$sample_id
sample_filt <- sample[keep]

# =====================================================
# PCA
# =====================================================
pca <- snpgdsPCA(genofile,
                 autosome.only = FALSE,
                 eigen.cnt = 20,
                 sample.id = sample_filt)

pc.percent <- (pca$eigenval / sum(pca$eigenval)) * 100

pca_df <- data.frame(
  sample.id = pca$sample.id,
  PC1 = pca$eigenvect[, 1],
  PC2 = pca$eigenvect[, 2],
  PC3 = pca$eigenvect[, 3],
  PC4 = pca$eigenvect[, 4]
)

pca_df$sample_short <- ifelse(
  grepl("^GBS", pca_df$sample.id),
  gsub("\\D", "", sapply(strsplit(sub("\\..*", "", pca_df$sample.id), "_"), `[`, 3)),
  pca_df$sample.id
)

pca_df <- inner_join(pca_df, metadata, by = c("sample_short" = "sample_id"))

# =====================================================
# PCA PLOT (K2 EXAMPLE)
# =====================================================
p <- ggplot(pca_df,
            aes(PC1, PC2,
                color = as.factor(K2_cluster_majoritaire))) +
  geom_point(aes(alpha = K2_statut), size = 1) +
  stat_ellipse(aes(group = K2_cluster_majoritaire),
               linetype = "dashed",
               linewidth = 0.3) +
  labs(
    x = paste0("PC1 (", round(pc.percent[1], 1), "%)"),
    y = paste0("PC2 (", round(pc.percent[2], 1), "%)")
  ) +
  theme_classic()

print(p)

# =====================================================
# CLUSTERS
# =====================================================
load_cluster_files <- function(k, type) {
  lapply(1:k, function(i) {
    readLines(file.path(cluster_dir,
                        sprintf("litt_K%d_cluster%d_%s.txt",
                                k, i, type)))
  })
}

clusters <- list(
  K3_raw  = load_cluster_files(3, "majo"),
  K3_pure = load_cluster_files(3, "pur"),
  K4_raw  = load_cluster_files(4, "majo"),
  K4_pure = load_cluster_files(4, "pur"),
  K7_raw  = load_cluster_files(7, "majo"),
  K7_pure = load_cluster_files(7, "pur")
)

# =====================================================
# GENOTYPE MATRIX
# =====================================================
geno_mat <- snpgdsGetGeno(genofile, sample.id = sample)
rownames(geno_mat) <- sample_new
colnames(geno_mat) <- snps.ids

# =====================================================
# GENIND
# =====================================================
make_genind <- function(geno_mat, cluster_lists, labels) {
  
  inds <- intersect(unique(unlist(cluster_lists)), rownames(geno_mat))
  sub <- geno_mat[inds, , drop = FALSE]
  
  pop <- rep(NA, length(inds))
  names(pop) <- inds
  
  for (i in seq_along(cluster_lists)) {
    pop[cluster_lists[[i]]] <- labels[i]
  }
  
  df2genind(sub,
            ploidy = 2,
            type = "codom",
            pop = factor(pop[inds]),
            ncode = 1)
}

# =====================================================
# NEI DA + NJ + BOOTSTRAP
# =====================================================
bootstrap_da_nj <- function(hf, nboot = 100) {
  
  da <- genet.dist(hf, method = "Da")
  tree <- nj(as.dist(da))
  
  loci <- 2:ncol(hf)
  boots <- list()
  
  for (i in 1:nboot) {
    
    hf_b <- hf[, c(1, sample(loci, replace = TRUE))]
    
    res <- try({
      nj(as.dist(genet.dist(hf_b, method = "Da")))
    }, silent = TRUE)
    
    if (!inherits(res, "try-error")) {
      boots[[length(boots) + 1]] <- res
    }
  }
  
  list(
    tree = tree,
    bootstrap = prop.clades(tree, boots)
  )
}

# =====================================================
# PIPELINE
# =====================================================
run_pipeline <- function(geno_mat, cluster_list, name, nboot = 100) {
  
  cluster_list <- lapply(cluster_list,
                         function(x) x[x %in% rownames(geno_mat)])
  
  labels <- paste0(name, "_C", seq_along(cluster_list))
  
  gi <- make_genind(geno_mat, cluster_list, labels)
  hf <- genind2hierfstat(gi)
  
  if (length(cluster_list) < 3) {
    return(list(type = "da_only",
                da = genet.dist(hf, method = "Da")))
  }
  
  bootstrap_da_nj(hf, nboot)
}

# =====================================================
# RUN ALL
# =====================================================
results <- lapply(names(clusters), function(n) {
  run_pipeline(geno_mat, clusters[[n]], n)
})
names(results) <- names(clusters)

# =====================================================
# SAVE EXAMPLE TREE
# =====================================================
pdf(file.path(fig_dir, "K3_raw_tree.pdf"), 4, 4)
plot(results$K3_raw$tree)
dev.off()
