
library(devtools)
library(tidyverse)
library(readxl)
# remotes::install_github("jiabowang/GAPIT", force = TRUE)
library(GAPIT)
library(reshape2)
library(compositions)
library(vegan)

# Set the project root directory before running the script
setwd("path/to/project")

# =====================================================
# 1. SAMPLE METADATA AND FLAGGED SAMPLES
# =====================================================

metadata_flag <- read_excel("1_tableau_bilan_echantillons.xlsx")

flagged_samples <- metadata_flag %>%
  filter(`Flag-weirdposition` == "Y") %>%
  pull(`numéro échantillon`)

# =====================================================
# 2. ENVIRONMENTAL AND POPULATION STRUCTURE COVARIATES
# =====================================================

# Load climatic and genomic covariates

map <- read_tsv(
  "F_VCF_finaux_pruned_unpruned_imputed/for_pheno_GWAS/3D_filtered_maf0_miss08_qual_noindel_biallelic_renamed_mind08_hz06_imputed_recode.nmap"
)
map <- as.data.frame(map)

# Load population structure and biogeographic covariates
# Samples flagged for environmental data issues are excluded

covar <- read_tsv("final3d_Q_biogeo_matrix.tsv")

covar <- covar %>%
  filter(!numero %in% c(78, 425)) %>% # Missing environmental values
  distinct(geometry, .keep_all = TRUE)

covar_filtered <- covar %>%
  mutate(across(1, as.character))

covar_filtered <- covar_filtered[, c(1, 15:18)]

# Apply an isometric log-ratio (ILR) transformation to Q proportions

Q <- as.matrix(covar_filtered[, 2:5])
Q <- Q / rowSums(Q)
Q[Q <= 0] <- 1e-6
Q_ilr <- ilr(Q)

covar_filtered <- cbind(
  covar_filtered[, 1],
  as.data.frame(Q_ilr)
)

# Load centered, scaled and transformed environmental variables

transformed_covar <- read_tsv(
  "GBS_3D_rerun/figures_last/tables/environmental_centered_scaled_transformed.tsv"
)

covar_filtered <- cbind(covar_filtered, transformed_covar)

# Load PCA scores used as GWAS covariates

pca_forcovar <- read_tsv(
  "GBS_3D_rerun/figures_last/tables/PCA_prunded_forcovarGWAS"
)

pca_forcovar_filtered <- pca_forcovar[
  pca_forcovar$sample_short %in% covar_filtered$numero,
]

pca_forcovar_filtered <- pca_forcovar_filtered[
  match(covar_filtered$numero, pca_forcovar_filtered$sample_short),
]

# Uncomment to exclude population structure covariates
# when using the three PCA axes
# covar_filtered <- covar_filtered[, 1:12]

# =====================================================
# 3. GENOTYPE DATA
# =====================================================

geno <- read_tsv(
  "F_VCF_finaux_pruned_unpruned_imputed/for_pheno_GWAS/3D_filtered_maf0_miss08_qual_noindel_biallelic_renamed_mind08_hz06_imputed_recode.dat"
)

geno$taxa <- sapply(
  strsplit(sub("\\..*", "", geno$taxa), "_"),
  `[`,
  3
)

geno$taxa <- gsub("\\D", "", geno$taxa)

geno_filtered <- geno %>%
  filter(taxa %in% covar_filtered$numero)

geno_filtered <- as.data.frame(geno_filtered)

rownames(geno_filtered) <- geno_filtered$taxa

# Match genotype rows to the order of the covariate table

geno_filtered <- geno_filtered[covar_filtered$numero, ]

# =====================================================
# 4. PHENOTYPE DATA
# =====================================================

pheno <- read_excel("GWAS/1_tableau-phéno.xlsx")

# Retain samples present in the covariate dataset

pheno_filtered <- pheno[
  pheno$sample %in% covar_filtered$numero,
]

# Reorder phenotype rows to match the covariate order

pheno_filtered <- pheno_filtered[
  match(pheno_filtered$sample, covar_filtered$numero),
]

# Select the phenotypic traits used in the analyses

pheno_filtered_num <- pheno_filtered[, c(1, 11, 12, 13, 30, 31)]

pheno_filtered_num <- lapply(
  pheno_filtered_num,
  function(x) as.numeric(as.character(x))
)

pheno_filtered_num <- as.data.frame(pheno_filtered_num)

# Calculate inflorescence density

pheno_filtered_num <- pheno_filtered_num %>%
  mutate(inflo_density = nb_noeuds / long_inflo) %>%
  select(-nb_noeuds, -long_inflo)

# Prepare the complete set of phenotypic traits for plotting

pheno_filtered_num2 <- pheno_filtered[
  ,
  c(1, 11, 12, 13, 22, 23, 24, 25, 26, 27, 30, 31)
]

pheno_filtered_num2 <- pheno_filtered_num2 %>%
  mutate(
    `Inflorescence density (flower/cm)` = nb_noeuds / long_inflo
  ) %>%
  rowwise() %>%
  mutate(
    `Leaflet length/width index` = mean(
      c(
        as.numeric(long_fol_centr_1) / as.numeric(larg_fol_centr_1),
        as.numeric(long_fol_centr_2) / as.numeric(larg_fol_centr_2),
        as.numeric(long_fol_centr_3) / as.numeric(larg_fol_centr_3)
      ),
      na.rm = TRUE
    )
  ) %>%
  ungroup() %>%
  select(
    -nb_noeuds,
    -long_inflo,
    -larg_fol_centr_1,
    -larg_fol_centr_2,
    -larg_fol_centr_3,
    -long_fol_centr_1,
    -long_fol_centr_2,
    -long_fol_centr_3
  )

# Convert phenotype data to long format

df_long <- melt(
  pheno_filtered_num2,
  id.vars = "sample",
  variable.name = "Trait",
  value.name = "Value"
)

# =====================================================
# 5. LOAD STRUCTURE CLUSTER ASSIGNMENTS
# =====================================================

# K2: individuals assigned to clusters 1 and 2, using maximum Q

k2_c1_raw <- readLines(
  "F_listes_clusters_STRUCTURE_3d_litt/3d_K2_cluster1_majo.txt"
)
k2_c2_raw <- readLines(
  "F_listes_clusters_STRUCTURE_3d_litt/3d_K2_cluster2_majo.txt"
)

# K2: individuals assigned to clusters 1 and 2, with Q > 0.7

k2_c1_pure <- readLines(
  "F_listes_clusters_STRUCTURE_3d_litt/3d_K2_cluster1_pur.txt"
)
k2_c2_pure <- readLines(
  "F_listes_clusters_STRUCTURE_3d_litt/3d_K2_cluster2_pur.txt"
)

# K4: individuals assigned to clusters 1, 2, 3 and 4,
# using maximum Q

k4_c1_raw <- readLines(
  "F_listes_clusters_STRUCTURE_3d_litt/3d_K4_cluster1_majo.txt"
)
k4_c2_raw <- readLines(
  "F_listes_clusters_STRUCTURE_3d_litt/3d_K4_cluster2_majo.txt"
)
k4_c3_raw <- readLines(
  "F_listes_clusters_STRUCTURE_3d_litt/3d_K4_cluster3_majo.txt"
)
k4_c4_raw <- readLines(
  "F_listes_clusters_STRUCTURE_3d_litt/3d_K4_cluster4_majo.txt"
)

# K4: individuals assigned to clusters 1, 2, 3 and 4,
# with Q > 0.7

k4_c1_pure <- readLines(
  "F_listes_clusters_STRUCTURE_3d_litt/3d_K4_cluster1_pur.txt"
)
k4_c2_pure <- readLines(
  "F_listes_clusters_STRUCTURE_3d_litt/3d_K4_cluster2_pur.txt"
)
k4_c3_pure <- readLines(
  "F_listes_clusters_STRUCTURE_3d_litt/3d_K4_cluster3_pur.txt"
)
k4_c4_pure <- readLines(
  "F_listes_clusters_STRUCTURE_3d_litt/3d_K4_cluster4_pur.txt"
)

# =====================================================
# 6. K4 CLUSTER ASSIGNMENTS
# =====================================================

# Assign cluster membership to each sample

df_long$cluster <- NA

df_long$cluster[df_long$sample %in% k4_c1_raw] <- "Cluster 1"
df_long$cluster[df_long$sample %in% k4_c2_raw] <- "Cluster 2"
df_long$cluster[df_long$sample %in% k4_c3_raw] <- "Cluster 3"
df_long$cluster[df_long$sample %in% k4_c4_raw] <- "Cluster 4"

df_long$cluster_pure <- NA

df_long$cluster_pure[df_long$sample %in% k4_c1_pure] <- "Cluster 1"
df_long$cluster_pure[df_long$sample %in% k4_c2_pure] <- "Cluster 2"
df_long$cluster_pure[df_long$sample %in% k4_c3_pure] <- "Cluster 3"
df_long$cluster_pure[df_long$sample %in% k4_c4_pure] <- "Cluster 4"

# =====================================================
# 7. K4 BOXPLOTS
# =====================================================

# Check cluster assignment completeness

table(df_long$cluster, useNA = "ifany")

df_long <- df_long %>%
  mutate(
    Trait = recode(
      Trait,
      "Internode_lenght" = "Internode length (cm)",
      "Stem_circumference" = "Stem circumference (cm)",
      "height" = "Height (cm)",
      "inflo_density" = "Inflorescence density (flower/cm)"
    )
  ) %>%
  rename(Cluster = cluster) %>%
  rename(`Cluster Q > 0.7` = cluster_pure)

# Dataset including all assigned individuals

df_cluster <- df_long %>%
  select(-`Cluster Q > 0.7`)

# Dataset restricted to individuals with Q > 0.7

df_core <- df_long %>%
  filter(
    !is.na(`Cluster Q > 0.7`),
    trimws(`Cluster Q > 0.7`) != ""
  ) %>%
  mutate(
    Cluster = paste0(`Cluster Q > 0.7`, " Q > 0.7")
  ) %>%
  select(-`Cluster Q > 0.7`)

# Combine the two datasets

df_plot <- bind_rows(df_cluster, df_core)

df_plot$Value <- as.numeric(df_plot$Value)
df_plot$Cluster <- factor(df_plot$Cluster)

df_cluster$Value <- as.numeric(df_cluster$Value)
df_cluster$Cluster <- factor(df_cluster$Cluster)

df_core$Value <- as.numeric(df_core$Value)
df_core$Cluster <- factor(df_core$Cluster)

# Define cluster colours

cluster_colors <- c(
  "Cluster 1" = "mediumpurple3",
  "Cluster 2" = "palegreen3",
  "Cluster 3" = "tomato3",
  "Cluster 4" = "steelblue3",
  "Cluster 1 Q > 0.7" = "purple",
  "Cluster 2 Q > 0.7" = "green",
  "Cluster 3 Q > 0.7" = "red",
  "Cluster 4 Q > 0.7" = "blue"
)

# Boxplots for all assigned individuals

boxplots_nostat_raw <- ggplot(
  df_cluster,
  aes(x = Cluster, y = Value, fill = Cluster)
) +
  geom_boxplot(
    outlier.alpha = 0.4,
    outlier.size = 0.5,
    linewidth = 0.2
  ) +
  facet_wrap(~ Trait, scales = "free_y") +
  scale_fill_manual(values = cluster_colors) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  ) +
  labs(x = "", y = "Value", fill = "")

boxplots_nostat_raw

# Boxplots for individuals with Q > 0.7

boxplots_nostat_pure <- ggplot(
  df_core,
  aes(x = Cluster, y = Value, fill = Cluster)
) +
  geom_boxplot(
    outlier.alpha = 0.4,
    outlier.size = 0.5,
    linewidth = 0.2
  ) +
  facet_wrap(~ Trait, scales = "free_y") +
  scale_fill_manual(values = cluster_colors) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  ) +
  labs(x = "", y = "Value", fill = "")

boxplots_nostat_pure

# =====================================================
# 8. K4 STATISTICAL TESTS
# =====================================================

library(rstatix)

# Count non-missing observations per trait and cluster

n_obs <- df_plot %>%
  group_by(Trait, Cluster) %>%
  summarise(
    n = sum(!is.na(Value)),
    .groups = "drop"
  )

n_obs

# Fit a linear model for each trait and extract residuals

residual_df <- df_plot %>%
  group_by(Trait) %>%
  group_modify(~ {
    model <- lm(
      Value ~ Cluster,
      data = .x,
      na.action = na.exclude
    )

    .x$residuals <- residuals(model)
    .x
  })

library(ggplot2)

# Inspect residual normality using Q-Q plots

ggplot(
  residual_df,
  aes(sample = residuals, color = Cluster)
) +
  stat_qq(alpha = 0.6) +
  stat_qq_line(color = "black") +
  facet_wrap(~ Trait, scales = "free") +
  theme_bw()

library(car)
library(dplyr)

# Test homogeneity of variances using Levene's test

levene_results <- df_plot %>%
  group_by(Trait) %>%
  group_modify(~ {
    test <- leveneTest(Value ~ Cluster, data = .x)

    data.frame(
      p_value = test$`Pr(>F)`[1]
    )
  })

levene_results

# Q-Q plots indicate deviations from normality.
# Variance homogeneity is observed only for the leaflet index.

library(FSA)
library(purrr)
library(multcompView)

# Run Kruskal-Wallis tests, Dunn post-hoc tests,
# and generate compact letter displays

run_kw_dunn_letters <- function(df) {

  # Kruskal-Wallis tests

  kw_results <- df %>%
    group_by(Trait) %>%
    summarise(
      kw = list(
        kruskal.test(
          Value ~ Cluster,
          data = pick(everything())
        )
      ),
      .groups = "drop"
    ) %>%
    mutate(
      statistic = sapply(kw, \(x) x$statistic),
      p.value = sapply(kw, \(x) x$p.value),
      df = sapply(kw, \(x) x$parameter)
    ) %>%
    select(-kw)

  print(kw_results)

  # Dunn post-hoc tests with Benjamini-Hochberg correction

  dunn_results <- df %>%
    group_split(Trait) %>%
    map_df(function(d) {

      res <- dunnTest(
        Value ~ Cluster,
        data = d,
        method = "bh"
      )$res

      res$Trait <- unique(d$Trait)
      res
    })

  # Generate compact letter displays

  letters_df <- dunn_results %>%
    group_split(Trait) %>%
    map_df(function(d) {

      trait_name <- unique(d$Trait)

      pvals <- d$P.adj
      names(pvals) <- gsub(" - ", "-", d$Comparison)

      letters <- multcompLetters(
        pvals,
        threshold = 0.05
      )$Letters

      data.frame(
        Trait = trait_name,
        Cluster = names(letters),
        Letters = letters
      )
    })

  # Preserve cluster factor order

  letters_df$Cluster <- factor(
    letters_df$Cluster,
    levels = levels(df$Cluster)
  )

  # Position statistical labels above the boxes

  letters_df <- letters_df %>%
    left_join(
      df %>%
        group_by(Trait) %>%
        summarise(
          y = max(Value, na.rm = TRUE) * 1.10,
          .groups = "drop"
        ),
      by = "Trait"
    )

  return(list(
    kw = kw_results,
    dunn = dunn_results,
    letters = letters_df
  ))
}

library(ggpubr)

stats_cluster <- run_kw_dunn_letters(df_cluster)
stats_core <- run_kw_dunn_letters(df_core)

# Add statistical letters to the K4 boxplots: all individuals

boxplots_stat_raw <- boxplots_nostat_raw +
  geom_text(
    data = stats_cluster$letters,
    aes(
      x = Cluster,
      y = y,
      label = Letters
    ),
    inherit.aes = FALSE,
    size = 1.5,
    color = "black"
  ) +
  labs(y = NULL) +
  theme_bw(base_size = 7) +
  theme(
    text = element_text(family = "sans", size = 8),
    strip.text = element_text(size = 6),
    axis.text = element_text(size = 4),
    legend.position = "none",
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 7),
    plot.title = element_text(size = 5, face = "bold")
  )

boxplots_stat_raw

ggsave(
  "Papiers/Article_structure/figures_vecto/figure4_phenotypeK4_raw_withDunntest_afterKruskal.pdf",
  boxplots_stat_raw,
  width = 8,
  height = 5,
  device = cairo_pdf,
  units = "cm",
  dpi = 600
)

# Add statistical letters to the K4 boxplots: Q > 0.7

boxplots_stat_pure <- boxplots_nostat_pure +
  geom_text(
    data = stats_core$letters,
    aes(
      x = Cluster,
      y = y,
      label = Letters
    ),
    inherit.aes = FALSE,
    size = 1.5,
    color = "black"
  ) +
  labs(y = NULL) +
  theme_bw(base_size = 7) +
  theme(
    text = element_text(family = "sans", size = 8),
    strip.text = element_text(size = 6),
    legend.position = "none",
    axis.text = element_text(size = 4),
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 7),
    plot.title = element_text(size = 5, face = "bold")
  )

boxplots_stat_pure

ggsave(
  "Papiers/Article_structure/figures_vecto/figure4_phenotypeK4_pure_withDunntest_afterKruskal.pdf",
  boxplots_stat_pure,
  width = 8,
  height = 5,
  device = cairo_pdf,
  units = "cm",
  dpi = 600
)

# =====================================================
# 9. K2 CLUSTER ASSIGNMENTS
# =====================================================

df_long$Cluster <- NA
df_long$Cluster[df_long$sample %in% k2_c1_raw] <- "Cluster 1"
df_long$Cluster[df_long$sample %in% k2_c2_raw] <- "Cluster 2"

df_long$`Cluster Q > 0.7` <- NA
df_long$`Cluster Q > 0.7`[
  df_long$sample %in% k2_c1_pure
] <- "Cluster 1"
df_long$`Cluster Q > 0.7`[
  df_long$sample %in% k2_c2_pure
] <- "Cluster 2"

# Split the data into all individuals and Q > 0.7 individuals

df_cluster <- df_long %>%
  select(-`Cluster Q > 0.7`)

df_core <- df_long %>%
  filter(
    !is.na(`Cluster Q > 0.7`),
    trimws(`Cluster Q > 0.7`) != ""
  ) %>%
  mutate(
    Cluster = paste0(`Cluster Q > 0.7`, " Q > 0.7")
  ) %>%
  select(-`Cluster Q > 0.7`)

df_cluster$Value <- as.numeric(df_cluster$Value)
df_core$Value <- as.numeric(df_core$Value)

df_cluster$Cluster <- factor(df_cluster$Cluster)
df_core$Cluster <- factor(df_core$Cluster)

# Define K2 cluster colours

cluster_colors <- c(
  "Cluster 1" = "steelblue3",
  "Cluster 2" = "tomato3",
  "Cluster 1 Q > 0.7" = "blue",
  "Cluster 2 Q > 0.7" = "red"
)

# =====================================================
# 10. K2 BOXPLOTS
# =====================================================

# Boxplots for all assigned individuals

boxplots_nostat_raw <- ggplot(
  df_cluster,
  aes(x = Cluster, y = Value, fill = Cluster)
) +
  geom_boxplot(
    outlier.alpha = 0.4,
    outlier.size = 0.5,
    linewidth = 0.2
  ) +
  facet_wrap(~ Trait, scales = "free_y") +
  scale_fill_manual(values = cluster_colors) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )

boxplots_nostat_raw

# Boxplots for individuals with Q > 0.7

boxplots_nostat_pure <- ggplot(
  df_core,
  aes(x = Cluster, y = Value, fill = Cluster)
) +
  geom_boxplot(
    outlier.alpha = 0.4,
    outlier.size = 0.5,
    linewidth = 0.2
  ) +
  facet_wrap(~ Trait, scales = "free_y") +
  scale_fill_manual(values = cluster_colors) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )

boxplots_nostat_pure

# =====================================================
# 11. K2 WILCOXON TESTS
# =====================================================

# Run Wilcoxon tests and generate compact letter displays

run_wilcox_letters <- function(df) {

  # Wilcoxon tests

  wilcox_results <- df %>%
    group_by(Trait) %>%
    summarise(
      test = list(
        wilcox.test(
          Value ~ Cluster,
          data = pick(everything())
        )
      ),
      .groups = "drop"
    ) %>%
    mutate(
      statistic = sapply(test, \(x) x$statistic),
      p.value = sapply(test, \(x) x$p.value)
    ) %>%
    select(-test)

  print(wilcox_results)

  # Generate compact letter displays for two clusters

  letters_df <- wilcox_results %>%
    rowwise() %>%
    mutate(
      letters = list({
        clusters <- unique(df$Cluster)

        if (p.value < 0.05) {
          setNames(letters[seq_along(clusters)], clusters)
        } else {
          setNames(rep("a", length(clusters)), clusters)
        }
      })
    ) %>%
    unnest_wider(letters) %>%
    pivot_longer(
      cols = -c(Trait, statistic, p.value),
      names_to = "Cluster",
      values_to = "Letters"
    )

  # Position statistical labels above the boxes

  letters_df <- letters_df %>%
    left_join(
      df %>%
        group_by(Trait) %>%
        summarise(
          y = max(Value, na.rm = TRUE) * 1.10,
          .groups = "drop"
        ),
      by = "Trait"
    )

  list(
    test = wilcox_results,
    letters = letters_df
  )
}

stats_k2_raw <- run_wilcox_letters(df_cluster)
stats_k2_pure <- run_wilcox_letters(df_core)

# Add statistical labels to the K2 boxplots: all individuals

boxplots_stat_k2_raw <- boxplots_nostat_raw +
  geom_text(
    data = stats_k2_raw$letters,
    aes(x = Cluster, y = y, label = Letters),
    inherit.aes = FALSE,
    size = 1.5,
    color = "black"
  ) +
  labs(y = NULL) +
  theme_bw(base_size = 7) +
  theme(
    text = element_text(family = "sans", size = 8),
    strip.text = element_text(size = 6),
    axis.text = element_text(size = 4),
    legend.position = "none"
  )

boxplots_stat_k2_raw

ggsave(
  "Papiers/Article_structure/figures_vecto/figure4_phenotypeK2_raw_withwilcoxtest.pdf",
  boxplots_stat_k2_raw,
  width = 8,
  height = 5,
  device = cairo_pdf,
  units = "cm",
  dpi = 600
)

# Add statistical labels to the K2 boxplots: Q > 0.7

boxplots_stat_k2_pure <- boxplots_nostat_pure +
  geom_text(
    data = stats_k2_pure$letters,
    aes(x = Cluster, y = y, label = Letters),
    inherit.aes = FALSE,
    size = 1.5,
    color = "black"
  ) +
  labs(y = NULL) +
  theme_bw(base_size = 7) +
  theme(
    text = element_text(family = "sans", size = 8),
    strip.text = element_text(size = 6),
    axis.text = element_text(size = 4),
    legend.position = "none"
  )

boxplots_stat_k2_pure

ggsave(
  "Papiers/Article_structure/figures_vecto/figure4_phenotypeK2_pure_withwilcoxtest.pdf",
  boxplots_stat_k2_pure,
  width = 8,
  height = 5,
  device = cairo_pdf,
  units = "cm",
  dpi = 600
)

# =====================================================
# 12. KINSHIP MATRIX
# =====================================================

# Load the centered IBS kinship matrix

# kinship_tuto_model <- read_tsv("GBS_3D_rerun/KSN.txt", col_names = FALSE)

kinship <- read_tsv(
  "GBS_3D_rerun/kinship_centeredIBS_8alleles.txt",
  col_names = FALSE
)

kinship <- as.data.frame(kinship)

kinship[, 1] <- sapply(
  strsplit(sub("\\..*", "", kinship[, 1]), "_"),
  `[`,
  3
)

kinship[, 1] <- gsub("\\D", "", kinship[, 1])

K_values <- kinship[, -1] %>%
  mutate_all(as.numeric)

# Assign row and column names

rownames(K_values) <- kinship[, 1]
colnames(K_values) <- kinship[, 1]

samples_to_keep <- covar_filtered$numero

# Retain the samples present in the covariate dataset

kinship_filtered <- K_values[
  samples_to_keep,
  samples_to_keep
]

kinship_filtered <- cbind(
  taxa = rownames(kinship_filtered),
  kinship_filtered
)

# =====================================================
# 13. VARIANCE PARTITIONING AND ENVIRONMENTAL RDA
# =====================================================

# Perform environmental and population structure analyses
# before the GWAS

keep <- complete.cases(
  pheno_filtered_num[, 2:5],
  covar_filtered[, 5:10],
  covar_filtered[, 2:4]
)

# Partition variance trait by trait

Q_ilr <- covar_filtered[, 2:4]
Y_scaled <- scale(pheno_filtered_num[, 2:5])
E_scaled <- covar_filtered[, 5:10]

Y <- Y_scaled[keep, ]
E <- E_scaled[keep, ]
G <- Q_ilr[keep, ]

# Retain PCA axes 1–3 for use as GWAS covariates

PC <- pca_forcovar_filtered[, 2:4]

PC_covar <- cbind(covar_filtered[, 1], PC)
PC_covar_biogeo <- cbind(PC_covar, covar_filtered[, 5:10])

# Variance partitioning between environmental and genetic effects

vp <- varpart(Y, E, G)
plot(vp)

E <- as.data.frame(E)
G <- as.data.frame(G)
Y <- as.data.frame(Y)

# Estimate the pure environmental effect, controlling for structure

rda_Epure <- rda(Y, E, G)
summary(rda_Epure)

# Estimate the pure population structure effect,
# controlling for environmental variables

rda_Gpure <- rda(Y, G, E)
summary(rda_Gpure)

# Test the significance of the constrained models

anova(rda_Epure, permutations = 999)
anova(rda_Gpure, permutations = 999)

# Equivalent model notation: Y ~ E + Condition(G)

# =====================================================
# 14. ENVIRONMENTAL RDA PLOT
# =====================================================

library(ggplot2)
library(vegan)
library(dplyr)
library(grid)

# Arrow scaling factor

arrow_scale <- 5

# Extract site scores

site_scores <- scores(
  rda_Epure,
  display = "sites",
  scaling = 2
) %>%
  as.data.frame()

site_scores$Sample <- rownames(site_scores)

# Extract environmental variable scores

var_scores <- scores(
  rda_Epure,
  display = "bp",
  scaling = 2
) %>%
  as.data.frame()

var_scores$Variable <- rownames(var_scores)

# Plot sites and environmental vectors

RDA <- ggplot() +
  # Sample scores
  geom_point(
    data = site_scores,
    aes(x = RDA1, y = RDA2),
    color = "blue1",
    size = 0.5,
    alpha = 0.5
  ) +

  # Environmental vectors, scaled for visualization
  geom_segment(
    data = var_scores,
    aes(
      x = 0,
      y = 0,
      xend = RDA1 * arrow_scale,
      yend = RDA2 * arrow_scale
    ),
    arrow = arrow(length = unit(0.25, "cm")),
    color = "black"
  ) +

  # Variable labels placed near the arrow tips
  geom_text(
    data = var_scores,
    aes(
      x = RDA1 * arrow_scale,
      y = RDA2 * arrow_scale,
      label = Variable
    ),
    color = "black",
    vjust = -0.5,
    hjust = 0.5
  ) +

  # Axis labels showing the percentage of constrained variance
  xlab(
    paste0(
      "RDA1 (",
      round(summary(rda_Epure)$cont$importance[2, 1] * 100, 1),
      "%)"
    )
  ) +
  ylab(
    paste0(
      "RDA2 (",
      round(summary(rda_Epure)$cont$importance[2, 2] * 100, 1),
      "%)"
    )
  ) +
  theme_bw()

RDA

ggsave(
  "Papiers/Article_structure/figures_vecto/figure4_phenotypeK2_RDA_env_sur_pheno.pdf",
  RDA,
  width = 8,
  height = 8,
  device = cairo_pdf,
  units = "cm",
  dpi = 600
)

# Preliminary interpretation: approximately 17% of the variance,
# including approximately 10% explained by the environment alone
# (approximately 4% for the pure environmental fraction).
# Additional covariate adjustment may be required.

# =====================================================
# 15. GWAS WITH GAPIT
# =====================================================

# Set the output directory for GWAS results

setwd(
  "GBS_3D_rerun/figures_last/tables/GWAS_output_wPCA&biogeo_kinship_corrected"
)

myGAPIT <- GAPIT(
  Y = pheno_filtered_num,
  GD = geno_filtered,
  GM = map,
  CV = PC_covar_biogeo, # Alternatives: PC_covar or covar_filtered (ILR Q)
  KI = kinship_filtered,
  model = c("MLMM"),
  SNP.MAF = 0.05
)

# Available models include:
# "Blink", "CMLM", "MLM", "SUPER", "MLMM", "FarmCPU"

# =====================================================
# 16. IMPORT AND COMBINE GWAS RESULTS
# =====================================================

folder <- "."

models <- c("BLINK", "MLM", "CMLM", "SUPER", "MLMM", "FarmCPU")

files_gwas <- unlist(
  lapply(models, function(mod) {

    list.files(
      path = folder,
      pattern = paste0(
        "GAPIT\\.Association\\.GWAS_Results\\.",
        mod,
        "\\..*\\(NYC\\)\\.csv$"
      ),
      full.names = TRUE
    )

  })
)

# Extract the GWAS model name from the filename

extract_model <- function(f) {

  fname <- basename(f)

  sub(
    "GAPIT\\.Association\\.GWAS_Results\\.([^.]+)\\..*",
    "\\1",
    fname
  )
}

# Extract the phenotype name from the filename

extract_trait <- function(f) {

  fname <- basename(f)

  sub(
    "GAPIT\\.Association\\.GWAS_Results\\.[^.]+\\.(.*?)\\(NYC\\)\\.csv",
    "\\1",
    fname
  )
}

# Import all GWAS results and identify the corresponding model and trait

gwas_list <- lapply(files_gwas, function(f) {

  df <- read.csv(f)

  df$Trait <- extract_trait(f)
  df$Model <- extract_model(f)

  # Specific correction for BLINK output column ordering

  if (df$Model[1] == "BLINK") {

    tmp <- df$MAF
    df$MAF <- df$nobs
    df$nobs <- tmp

  }

  df

})

gwas_all <- bind_rows(gwas_list)

# =====================================================
# 17. RENAME PHENOTYPIC TRAITS
# =====================================================

gwas_all <- gwas_all %>%
  mutate(
    Trait = dplyr::recode(
      Trait,
      "Internode_lenght" = "Internode length (cm)",
      "Stem_circumference" = "Stem circumference (cm)",
      "height" = "Height (cm)",
      "inflo_density" = "Inflorescence density (flower/cm)"
    )
  )

# =====================================================
# 18. PREPARE MANHATTAN PLOT COORDINATES
# =====================================================

gwas_all$logp <- -log10(gwas_all$P.value)

gwas_all$Chr_char <- as.character(gwas_all$Chr)

chr_order <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "X")

# Set chromosome factor levels in the desired order

gwas_all <- gwas_all %>%
  mutate(
    Chr_factor = factor(Chr_char, levels = chr_order)
  )

# Calculate cumulative chromosome coordinates

chr_info <- gwas_all %>%
  group_by(Chr_factor) %>%
  summarise(
    chr_len = max(Pos, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(factor(Chr_factor, levels = chr_order)) %>%
  mutate(
    chr_start = lag(cumsum(chr_len), default = 0)
  )

gwas_plot <- gwas_all %>%
  left_join(
    chr_info %>% select(Chr_factor, chr_start),
    by = "Chr_factor"
  ) %>%
  mutate(
    Pos_cum = Pos + chr_start,
    logp = -log10(P.value)
  )

axis_df <- chr_info %>%
  mutate(
    center = chr_start + chr_len / 2
  )

# =====================================================
# 19. BONFERRONI SIGNIFICANCE THRESHOLDS
# =====================================================

n_snp <- 8842 # Number of markers after MAF filtering at 0.05

bonf05 <- -log10(0.05 / n_snp)
bonf01 <- -log10(0.01 / n_snp)

# =====================================================
# 20. MANHATTAN PLOT
# =====================================================

# Define colours for the GWAS models

model_cols <- c(
  "BLINK" = "#D55E00",
  "FarmCPU" = "#0072B2",
  "MLM" = "#009E73",
  "CMLM" = "#CC79A7",
  "SUPER" = "#E69F00",
  "MLMM" = "#56B4E9"
)

# Recalculate chromosome coordinates for plotting

chr_info <- gwas_all %>%
  group_by(Chr_factor) %>%
  summarise(
    chr_len = max(Pos),
    .groups = "drop"
  ) %>%
  arrange(factor(Chr_factor, levels = chr_order)) %>%
  mutate(
    chr_start = lag(cumsum(chr_len), default = 0)
  )

gwas_plot <- gwas_all %>%
  left_join(
    chr_info %>% select(Chr_factor, chr_start),
    by = "Chr_factor"
  ) %>%
  mutate(
    Pos_cum = Pos + chr_start,
    logp = -log10(P.value)
  )

axis_df <- chr_info %>%
  mutate(
    center = chr_start + chr_len / 2
  )

# Calculate chromosome separation positions

vlines <- (
  axis_df$center[-1] +
    axis_df$center[-nrow(axis_df)]
) / 2

library(ggrastr)

# Plot GWAS associations across chromosomes and phenotypic traits

man <- ggplot(
  gwas_plot,
  aes(
    x = Pos_cum,
    y = logp,
    colour = Model
  )
) +
  geom_vline(
    xintercept = vlines,
    colour = "grey80",
    linewidth = 0.1
  ) +
  rasterise(
    geom_point(
      alpha = 0.6,
      size = 0.1
    ),
    dpi = 600
  ) +
  geom_hline(
    yintercept = bonf05,
    linetype = "solid",
    linewidth = 0.1
  ) +
  geom_hline(
    yintercept = bonf01,
    linetype = "dashed",
    linewidth = 0.1
  ) +
  facet_wrap(
    ~ Trait,
    ncol = 2,
    scales = "fixed"
  ) +
  scale_colour_manual(
    values = model_cols,
    name = "GWAS model"
  ) +
  scale_x_continuous(
    breaks = axis_df$center,
    labels = chr_order
  ) +
  labs(
    x = "Chromosome",
    y = expression(-log[10](P))
  ) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 4),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 6)
    # legend.position = "none"
  )

man

ggsave(
  "../../../../Papiers/Article_structure/figures_vecto/figure4_v2_manhattan_all_covarPCA_kinship_biogeo_legend.pdf",
  man,
  width = 8,
  height = 6,
  device = cairo_pdf,
  units = "cm",
  dpi = 600
)

# =====================================================
# 21. EXPORT SIGNIFICANT GWAS ASSOCIATIONS
# =====================================================

# Retain associations with P < 0.0001

gwas_sig <- gwas_all %>%
  filter(P.value < 0.0001)

library(openxlsx)

# Export significant associations to Excel

write.xlsx(
  gwas_sig,
  "../../../../Papiers/Article_structure/figures_vecto/gwas_significant_pvalue_0.0001.xlsx"
)

# =====================================================
# 22. STANDALONE LEGEND FOR THE MANHATTAN PLOT
# =====================================================

# Define trait labels and their corresponding colours

legend_df <- tibble(
  Trait = c(
    "Internode length (cm)",
    "Stem circumference (cm)",
    "Height (cm)",
    "Inflorescence density (flower/cm)"
  ),
  col_trait = c(
    "#1b9e77",
    "#d95f02",
    "#7570b3",
    "#e7298a"
  )
)

# Create the standalone legend

legend <- ggplot(
  legend_df,
  aes(x = 1, y = Trait)
) +
  geom_point(
    aes(color = Trait),
    size = 5
  ) +
  geom_text(
    aes(label = Trait),
    hjust = -0.1
  ) +
  scale_color_manual(
    values = setNames(legend_df$col_trait, legend_df$Trait)
  ) +
  theme_void() +
  theme(
    legend.position = "none"
  ) +
  xlim(1, 1.5)

ggsave(
  "../../../../Papiers/Article_structure/figures_vecto/figure4_manhattan_legend.pdf",
  legend,
  width = 6,
  height = 4,
  device = cairo_pdf,
  units = "cm",
  dpi = 600
)
