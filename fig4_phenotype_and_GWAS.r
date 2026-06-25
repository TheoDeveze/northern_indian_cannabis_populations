############################################
## GWAS + PHENOTYPES + STRUCTURE (SIMPLIFIÉ)
############################################

library(tidyverse)
library(readxl)
library(reshape2)
library(compositions)
library(vegan)
library(GAPIT)

# =========================================
# 1) DATA IMPORT
# =========================================

metadata_flag <- read_excel("data/metadata_samples.xlsx")

flagged_samples <- metadata_flag %>%
  filter(`Flag-weirdposition` == "Y") %>%
  pull(`sample_id`)

# Covariables (structure + environnement)
covar <- read_tsv("data/final_biogeo_matrix.tsv") %>%
  filter(!numero %in% c(78, 425)) %>%
  distinct(geometry, .keep_all = TRUE)

covar <- covar %>%
  mutate(across(1, as.character)) %>%
  select(1, 15:18)

# STRUCTURE Q -> ILR
Q <- as.matrix(covar[, 2:5])
Q <- Q / rowSums(Q)
Q[Q <= 0] <- 1e-6
Q_ilr <- ilr(Q)

covar <- cbind(covar[,1], as.data.frame(Q_ilr))

# Environnement transformé (déjà préparé)
env <- read_tsv("data/environment_scaled.tsv")
covar <- cbind(covar, env)

# PCA covariables (structure GWAS)
pca_cov <- read_tsv("data/pca_covar.tsv")

# =========================================
# 2) GENOTYPE + PHENOTYPE
# =========================================

geno <- read_tsv("data/genotype_filtered.dat")
geno$taxa <- gsub("\\D", "", geno$taxa)

geno <- geno %>%
  filter(taxa %in% covar$numero)

rownames(geno) <- geno$taxa
geno <- geno[covar$numero, ]

pheno <- read_excel("data/phenotypes.xlsx")

pheno <- pheno %>%
  filter(sample %in% covar$numero) %>%
  slice(match(covar$numero, sample))

# variables numériques
pheno_num <- pheno %>%
  select(sample, 11:13, 30:31) %>%
  mutate(across(-sample, as.numeric)) %>%
  mutate(
    inflo_density = nb_noeuds / long_inflo
  ) %>%
  select(-nb_noeuds, -long_inflo)

# =========================================
# 3) CLUSTERS STRUCTURE (K2 / K4)
# =========================================

k4_files <- list(
  c1 = "clusters/K4_c1.txt",
  c2 = "clusters/K4_c2.txt",
  c3 = "clusters/K4_c3.txt",
  c4 = "clusters/K4_c4.txt"
)

k4 <- lapply(k4_files, readLines)

# melt phenotype for plotting
df_long <- melt(pheno_num, id.vars = "sample")

df_long$cluster <- NA
df_long$cluster[df_long$sample %in% k4[[1]]] <- "C1"
df_long$cluster[df_long$sample %in% k4[[2]]] <- "C2"
df_long$cluster[df_long$sample %in% k4[[3]]] <- "C3"
df_long$cluster[df_long$sample %in% k4[[4]]] <- "C4"

# =========================================
# 4) BOXPLOTS (RAW)
# =========================================

df_long$Trait <- recode(df_long$variable,
                        "Internode_lenght" = "Internode length",
                        "Stem_circumference" = "Stem circumference",
                        "height" = "Height",
                        "inflo_density" = "Inflorescence density")

ggplot(df_long, aes(cluster, value, fill = cluster)) +
  geom_boxplot(linewidth = 0.2) +
  facet_wrap(~ Trait, scales = "free_y") +
  theme_bw()

# =========================================
# 5) TESTS NON PARAMETRIQUES
# =========================================

library(rstatix)
library(FSA)
library(multcompView)

kw_dunn <- function(df){
  
  kw <- df %>%
    group_by(Trait) %>%
    summarise(test = list(kruskal.test(value ~ cluster, data = cur_data())))
  
  dunn <- df %>%
    group_split(Trait) %>%
    map_df(~dunnTest(value ~ cluster, data = ., method = "bh")$res)
  
  list(kw = kw, dunn = dunn)
}

stats <- kw_dunn(df_long)

# =========================================
# 6) K2 ANALYSIS (SIMPLIFIÉ)
# =========================================

df_long$cluster_k2 <- NA
df_long$cluster_k2[df_long$sample %in% readLines("clusters/K2_c1.txt")] <- "C1"
df_long$cluster_k2[df_long$sample %in% readLines("clusters/K2_c2.txt")] <- "C2"

ggplot(df_long, aes(cluster_k2, value, fill = cluster_k2)) +
  geom_boxplot() +
  facet_wrap(~ Trait)

# =========================================
# 7) GWAS (GAPIT)
# =========================================

gwas <- GAPIT(
  Y = pheno_num,
  GD = geno,
  GM = read_tsv("data/map.nmap"),
  CV = pca_cov,
  KI = read_tsv("data/kinship.txt"),
  model = "Blink"
)

# =========================================
# 8) MANHATTAN PLOT (SIMPLIFIÉ)
# =========================================

gwas_res <- read.csv("data/gwas_results.csv")

gwas_res <- gwas_res %>%
  mutate(logp = -log10(P.value))

ggplot(gwas_res, aes(Pos, logp)) +
  geom_point(alpha = 0.5, size = 0.3) +
  theme_bw()
