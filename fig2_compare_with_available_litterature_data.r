library("tidyverse")
library("ggplot2")
library("readxl")
library("dplyr")
library(SNPRelate)

# Project folders (run this script from the project root directory).
project_dir <- "."
figure_dir <- file.path(project_dir, "results", "figures", "population_structure")
table_dir <- file.path(project_dir, "results", "tables")
cluster_dir <- file.path(project_dir, "data", "cluster_lists", "STRUCTURE_3d_literature")
dapc_out_dir <- file.path(table_dir, "cluster_lists_DAPC_K7")
gds_dir <- file.path(project_dir, "data", "intermediate", "genotypes", "literature")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dapc_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gds_dir, recursive = TRUE, showWarnings = FALSE)

# Run this script from the project root; all file paths are relative.

# Contains all results for K = 2, 3, 4, and 7.
metadata <- read_tsv(file.path(project_dir, "data", "metadata", "literature_metadata_STRUCTURE_DAPC_results.tsv"))
metadataK4 <-  read_tsv(file.path(project_dir, "data", "metadata", "3D_metadata_STRUCTURE_DAPC_results.tsv"))
metadata_pop <- read_tsv(file.path(project_dir, "data", "metadata", "literature", "PCA_literature_with_population_metadata.tsv"))
metadata_pop$population[
  metadata_pop$sample.id %in% as.character(metadataK4$sample_id)
] <- metadataK4$K4_cluster_majoritaire[
  match(
    metadata_pop$sample.id[
      metadata_pop$sample.id %in% as.character(metadataK4$sample_id)
    ],
    as.character(metadataK4$sample_id)
  )
]

metadata <- left_join(metadata, metadata_pop[,c(1,4,5)], by = c("sample_id"="sample.id"))

out_dir <- dapc_out_dir

if(!dir.exists(out_dir)){
  dir.create(out_dir)
}


# Loop through clusters 1 to 7.
for(i in 1:7){
  
  # ---------------------------
  # Majority-assignment list
  # ---------------------------
  
  cluster_majo <- metadata %>%
    filter(DAPC_K7_cluster_majoritaire == i) %>%
    pull(sample_id)
  
  
  writeLines(
    cluster_majo,
    file.path(
      out_dir,
      paste0("litt_K7_cluster", i, "_majo_DAPC.txt")
    )
  )
  
  
  # ---------------------------
  # Pure-assignment list
  # ---------------------------
  
  cluster_pur <- metadata %>%
    filter(
      DAPC_K7_cluster_majoritaire == i,
      DAPC_K7_statut == "pure"
    ) %>%
    pull(sample_id)
  
  
  writeLines(
    cluster_pur,
    file.path(
      out_dir,
      paste0("litt_K7_cluster", i, "_pur_DAPC.txt")
    )
  )
  
}

vcf.fn <- file.path(project_dir, "data", "genotypes", "literature", "litterature_pruned_maf001_miss09_qual_noindel_biallelic_renamed_mind08_hz06.vcf")
gds.fn <- file.path(gds_dir, "1211.gds")  # Intermediate GDS file

# First, calculate genome-wide PCA.
snpgdsVCF2GDS(vcf.fn, gds.fn, method="biallelic.only")
genofile <- snpgdsOpen(gds.fn)
sample <- read.gdsn(index.gdsn(genofile, "sample.id"))
sample_new <- ifelse(
  grepl("^GBS", sample),
  gsub("\\D", "", sapply(strsplit(sub("\\..*", "", sample), "_"), `[`, 3)),
  sample
)
keep <- sample_new %in% metadata$sample_id
sample_filt <- sample[keep]
snps.ids <- read.gdsn(index.gdsn(genofile, "snp.id"))

pca <- snpgdsPCA(genofile, autosome.only = FALSE, eigen.cnt = 20, sample.id = sample_filt) # Calculate only the first 20 PCs.
pca$eigenval   <- pca$eigenval[1:20]
pca$varprop <- pca$eigenval / sum(pca$eigenval)

pc.percent <- pca$varprop * 100
pca_df <- data.frame(
  sample.id = pca$sample.id,
  PC1 = pca$eigenvect[, 1],  # First principal component
  PC2 = pca$eigenvect[, 2],
  PC3 = pca$eigenvect[, 3],  
  PC4 = pca$eigenvect[, 4] 
)
pca_df$sample_short <- ifelse(
  grepl("^GBS", pca_df$sample.id),
  gsub("\\D", "", sapply(strsplit(sub("\\..*", "", pca_df$sample.id), "_"), `[`, 3)),
  pca_df$sample.id)

#write_tsv(pca_df, "results/tables/PCA_pruned_for_covariate_GWAS.tsv")
pca_df <- inner_join(pca_df, metadata, by = c("sample_short"="sample_id"))

# Convert eigenvalues to the percentage of variance explained.

# Plot the PCA.
pve <- data.frame(pve = pca$varprop * 100)  # Percentage of variance explained.

theme_paper <- function() {
  theme(
    text = element_text(family = "sans", size = 8),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    strip.text = element_text(size = 8),
    plot.title = element_text(size = 9, face = "bold")
  )
}

pca_plot_K2 <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    color = as.factor(K2_cluster_majoritaire),
    shape = study
  )
) +
  geom_point(
    aes(alpha = K2_statut),
    size = 1.0
  ) +
  
  # Ellipses
  stat_ellipse(
    data = pca_df,
    aes(
      x = PC1,
      y = PC2,
      group = as.factor(K2_cluster_majoritaire),
      color = as.factor(K2_cluster_majoritaire)
    ),
    inherit.aes = FALSE,
    type = "norm",
    linetype = "dashed",
    linewidth = 0.3,
    level = 0.68
  ) +
  
  labs(
    x = paste0("PC1 (", round(pc.percent[1], 1), "%)"),
    y = paste0("PC2 (", round(pc.percent[2], 1), "%)"),
    color = "Cluster",
    shape = "Study"
  ) +
  
  scale_color_manual(values = c(
    "1" = "blue",
    "2" = "red"
  )) +
  scale_alpha_manual(values = c(
    "admixed" = 0.3,
    "pure" = 0.9
  )) +
  scale_shape_manual(
    values = c(
      "India" = 17,             # Triangle up
      "Sawler et al 2015" = 15, # Square
      "Soorni et al 2017" = 25, # Triangle down
      "Ren et al 2021" = 16,    # Circle
      "Chen et al 2022" = 9    # Diamond
    ),
    breaks = c(
      "India",
      "Sawler et al 2015",
      "Soorni et al 2017",
      "Ren et al 2021",
      "Chen et al 2022"
    ),
    labels = c(
      "India" = "India",
      "Sawler et al 2015" = "Sawler et al 2015",
      "Soorni et al 2017" = "Soorni et al 2017",
      "Ren et al 2021" = "Ren et al 2021",
      "Chen et al 2022" = "Chen et al 2022"
    )
  ) +
  coord_fixed() +
  
  theme_classic(base_size = 8) +
  theme_paper() +
  
  theme(
    legend.position = "none" # Set to "right" to display the legend or "none" to hide it.
  )

print(pca_plot_K2)
#ggsave("results/figures/population_structure/figure2_pca_K2_wo_leg.pdf", pca_plot_K2,
       #       width = 8,
       #       height = 8,
       #       device = cairo_pdf,
       #       units = "cm",
       #       dpi = 600)

pca_plot_K3 <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    color = as.factor(K3_cluster_majoritaire),
    shape = study
  )
) +
  geom_point(
    aes(alpha = K3_statut),
    size = 1.0
  ) +
  
  # Ellipses
  stat_ellipse(
    data = pca_df,
    aes(
      x = PC1,
      y = PC2,
      group = as.factor(K3_cluster_majoritaire),
      color = as.factor(K3_cluster_majoritaire)
    ),
    inherit.aes = FALSE,
    type = "norm",
    linetype = "dashed",
    linewidth = 0.3,
    level = 0.68
  ) +
  
  labs(
    x = paste0("PC1 (", round(pc.percent[1], 1), "%)"),
    y = paste0("PC2 (", round(pc.percent[2], 1), "%)"),
    color = "Cluster",
    shape = "Study"
  ) +
  
  scale_color_manual(values = c(
    "1" = "#1f77b4",  # blue
    "2" = "#d62728",  # red
    "3" = "#2ca02c"   # green
  )) +
  scale_alpha_manual(values = c(
    "admixed" = 0.3,
    "pure" = 0.9
  )) +
  scale_shape_manual(
    values = c(
      "India" = 17,             # Triangle up
      "Sawler et al 2015" = 15, # Square
      "Soorni et al 2017" = 25, # Triangle down
      "Ren et al 2021" = 16,    # Circle
      "Chen et al 2022" = 9    # Diamond
    ),
    breaks = c(
      "India",
      "Sawler et al 2015",
      "Soorni et al 2017",
      "Ren et al 2021",
      "Chen et al 2022"
    ),
    labels = c(
      "India" = "India",
      "Sawler et al 2015" = "Sawler et al 2015",
      "Soorni et al 2017" = "Soorni et al 2017",
      "Ren et al 2021" = "Ren et al 2021",
      "Chen et al 2022" = "Chen et al 2022"
    )
  ) +
  coord_fixed() +
  
  theme_classic(base_size = 8) +
  theme_paper() +
  
  theme(
    legend.position = "right" # Set to "right" to display the legend or "none" to hide it.
  )

print(pca_plot_K3)
#ggsave("results/figures/population_structure/figure2_pca_K3_wi_leg.pdf", pca_plot_K3,
       #      width = 10,
       #      height = 14,
       #      device = cairo_pdf,
       #       units = "cm",
       #      dpi = 600)

pca_plot_K4 <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    color = as.factor(K4_cluster_majoritaire),
    shape = study
  )
) +
  geom_point(
    aes(alpha = K4_statut),
    size = 1.0
  ) +
  
  # Ellipses
  stat_ellipse(
    data = pca_df,
    aes(
      x = PC1,
      y = PC2,
      group = as.factor(K4_cluster_majoritaire),
      color = as.factor(K4_cluster_majoritaire)
    ),
    inherit.aes = FALSE,
    type = "norm",
    linetype = "dashed",
    linewidth = 0.3,
    level = 0.68
  ) +
  
  labs(
    x = paste0("PC1 (", round(pc.percent[1], 1), "%)"),
    y = paste0("PC2 (", round(pc.percent[2], 1), "%)"),
    color = "Cluster",
    shape = "Study"
  ) +
  
  scale_color_manual(values = c(
    "1" = "#1f77b4",  # blue
    "2" = "#d62728",  # red
    "3" = "#2ca02c",  # green
    "4" = "#ff7f0e"   # orange
  )) +
  scale_alpha_manual(values = c(
    "admixed" = 0.3,
    "pure" = 0.9
  )) +
  scale_shape_manual(
    values = c(
      "India" = 17,             # Triangle up
      "Sawler et al 2015" = 15, # Square
      "Soorni et al 2017" = 25, # Triangle down
      "Ren et al 2021" = 16,    # Circle
      "Chen et al 2022" = 9    # Diamond
    ),
    breaks = c(
      "India",
      "Sawler et al 2015",
      "Soorni et al 2017",
      "Ren et al 2021",
      "Chen et al 2022"
    ),
    labels = c(
      "India" = "India",
      "Sawler et al 2015" = "Sawler et al 2015",
      "Soorni et al 2017" = "Soorni et al 2017",
      "Ren et al 2021" = "Ren et al 2021",
      "Chen et al 2022" = "Chen et al 2022"
    )
  ) +
  coord_fixed() +
  
  theme_classic(base_size = 8) +
  theme_paper() +
  
  theme(
    legend.position = "none" # Set to "right" to display the legend or "none" to hide it.
  )

print(pca_plot_K4)
#ggsave("results/figures/population_structure/figure2_pca_K4_wo_leg.pdf", pca_plot_K4,
#       width = 10,
#      height = 14,
#     device = cairo_pdf,
#    units = "cm",
#       dpi = 600)

pca_plot_K7 <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    color = as.factor(K7_cluster_majoritaire),
    shape = study
  )
) +
  geom_point(
    aes(alpha = K7_statut),
    size = 1.0
  ) +
  
  # Ellipses
  stat_ellipse(
    data = pca_df,
    aes(
      x = PC1,
      y = PC2,
      group = as.factor(K7_cluster_majoritaire),
      color = as.factor(K7_cluster_majoritaire)
    ),
    inherit.aes = FALSE,
    type = "norm",
    linetype = "dashed",
    linewidth = 0.3,
    level = 0.68
  ) +
  
  labs(
    x = paste0("PC1 (", round(pc.percent[1], 1), "%)"),
    y = paste0("PC2 (", round(pc.percent[2], 1), "%)"),
    color = "Cluster",
    shape = "Study"
  ) +
  
  scale_color_manual(values = c(
    "1" = "#1f77b4",  # blue
    "2" = "#d62728",  # red
    "3" = "#2ca02c",  # green
    "4" = "#ff7f0e",  # orange
    "5" = "#9467bd",  # purple
    "6" = "#8c564b",  # brown
    "7" = "#e377c2"   # pink
  )) +
  scale_alpha_manual(values = c(
    "admixed" = 0.3,
    "pure" = 0.9
  )) +
  scale_shape_manual(
    values = c(
      "India" = 17,             # Triangle up
      "Sawler et al 2015" = 15, # Square
      "Soorni et al 2017" = 25, # Triangle down
      "Ren et al 2021" = 16,    # Circle
      "Chen et al 2022" = 9    # Diamond
    ),
    breaks = c(
      "India",
      "Sawler et al 2015",
      "Soorni et al 2017",
      "Ren et al 2021",
      "Chen et al 2022"
    ),
    labels = c(
      "India" = "India",
      "Sawler et al 2015" = "Sawler et al 2015",
      "Soorni et al 2017" = "Soorni et al 2017",
      "Ren et al 2021" = "Ren et al 2021",
      "Chen et al 2022" = "Chen et al 2022"
    )
  ) +
  coord_fixed() +
  
  theme_classic(base_size = 8) +
  theme_paper() +
  
  theme(
    legend.position = "right" # Set to "right" to display the legend or "none" to hide it.
  )

print(pca_plot_K7)
#ggsave("results/figures/population_structure/figure2_pca_K7_wi_leg.pdf", pca_plot_K7,
#      width = 10,
#      height = 14,
#     device = cairo_pdf,
#       units = "cm",
       #     dpi = 600)

pca_plot_K7_DAPC <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    color = as.factor(DAPC_K7_cluster_majoritaire),
    shape = study
  )
) +
  geom_point(
    aes(alpha = DAPC_K7_statut),
    size = 1.0
  ) +
  
  # Ellipses
  stat_ellipse(
    data = pca_df,
    aes(
      x = PC1,
      y = PC2,
      group = as.factor(DAPC_K7_cluster_majoritaire),
      color = as.factor(DAPC_K7_cluster_majoritaire)
    ),
    inherit.aes = FALSE,
    type = "norm",
    linetype = "dashed",
    linewidth = 0.3,
    level = 0.68
  ) +
  
  labs(
    x = paste0("PC1 (", round(pc.percent[1], 1), "%)"),
    y = paste0("PC2 (", round(pc.percent[2], 1), "%)"),
    color = "Cluster",
    shape = "Study"
  ) +
  
  scale_color_manual(values = c(
    "1" = "#1f77b4",  # blue
    "2" = "#d62728",  # red
    "3" = "#2ca02c",  # green
    "4" = "#ff7f0e",  # orange
    "5" = "#9467bd",  # purple
    "6" = "#8c564b",  # brown
    "7" = "#e377c2"   # pink
  )) +
  scale_alpha_manual(values = c(
    "admixed" = 0.3,
    "pure" = 0.9
  )) +
  scale_shape_manual(
    values = c(
      "India" = 17,             # Triangle up
      "Sawler et al 2015" = 15, # Square
      "Soorni et al 2017" = 25, # Triangle down
      "Ren et al 2021" = 16,    # Circle
      "Chen et al 2022" = 9    # Diamond
    ),
    breaks = c(
      "India",
      "Sawler et al 2015",
      "Soorni et al 2017",
      "Ren et al 2021",
      "Chen et al 2022"
    ),
    labels = c(
      "India" = "India",
      "Sawler et al 2015" = "Sawler et al 2015",
      "Soorni et al 2017" = "Soorni et al 2017",
      "Ren et al 2021" = "Ren et al 2021",
      "Chen et al 2022" = "Chen et al 2022"
    )
  ) +
  coord_fixed() +
  
  theme_classic(base_size = 8) +
  theme_paper() +
  
  theme(
    legend.position = "right" # Set to "right" to display the legend or "none" to hide it.
  )

print(pca_plot_K7_DAPC)
ggsave("results/figures/population_structure/figure2_pca_K7_DAPC_wi_leg.pdf", pca_plot_K7_DAPC,
      width = 10,
      height = 14,
     device = cairo_pdf,
       units = "cm",
     dpi = 600)

pca_plot_metadata <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    color = as.factor(population),
    shape = study
  )
) +
  geom_point(
    aes(),
    size = 1.0,
    alpha= 0.9
  ) +
  
  labs(
    x = paste0("PC1 (", round(pc.percent[1], 1), "%)"),
    y = paste0("PC2 (", round(pc.percent[2], 1), "%)"),
    color = "Cluster",
    shape = "Study"
  ) +
  
  scale_color_manual(values = c(
    # Fixed clusters
    "1" = "purple", 
    "2" = "green",
    "3" = "red",  
    "4" = "blue",
    
    # Other groups (colors selected for contrast and a restrained palette)
    "Basal cannabis"   = "#332288",  # dark purple
    "Chinese hemp"     = "#88CCEE",  # light blue
    "Hemp-type"        = "#117733",  # dark green
    "Iran"             = "#AA4499",  # plum
    "Chinese wild"     = "#999933",  # olive
    "Drug-type"        = "#44AA99",  # teal
    "Drug-type feral"  = "#DDCC77"   # sand
  ),    breaks = c(
    "1", "2", "3", "4",
    "Basal cannabis",
    "Chinese hemp", "Hemp-type",
    "Iran",
    "Chinese wild",
    "Drug-type", "Drug-type feral"
  ),
  labels = c(
    "1" = "1 (Indian Himalayas)",
    "2" = "2 (Indian Himalayas)",
    "3" = "3 (Indian Himalayas)",
    "4" = "4 (Indian Himalayas)",
    "Hemp-type" = "Hemp-type",
    "Chinese hemp" = "Chinese hemp",
    "Chinese wild" = "Chinese wild",
    "Basal cannabis" = "Chinese basal",
    "Iran" = "Iranian",
    "Drug-type" = "Drug-type",
    "Drug-type feral" = "Drug-type feral"
  )
  ) +
  scale_shape_manual(
    values = c(
      "India" = 17,             # Triangle up
      "Sawler et al 2015" = 15, # Square
      "Soorni et al 2017" = 25, # Triangle down
      "Ren et al 2021" = 16,    # Circle
      "Chen et al 2022" = 9    # Diamond
    ),
    breaks = c(
      "India",
      "Sawler et al 2015",
      "Soorni et al 2017",
      "Ren et al 2021",
      "Chen et al 2022"
    ),
    labels = c(
      "India" = "India",
      "Sawler et al 2015" = "Sawler et al 2015",
      "Soorni et al 2017" = "Soorni et al 2017",
      "Ren et al 2021" = "Ren et al 2021",
      "Chen et al 2022" = "Chen et al 2022"
    )
  ) +
  coord_fixed() +
  
  theme_classic(base_size = 8) +
  theme_paper() +
  
  theme(
    legend.position = "right" # Set to "right" to display the legend or "none" to hide it.
  )

print(pca_plot_metadata)
#ggsave("results/figures/population_structure/figure2_pca_metadata_wi_leg.pdf", pca_plot_metadata,
#      width = 10,
#       height = 18,
#       device = cairo_pdf,
#      units = "cm",
#       dpi = 600)

## Nei's DA distances and STRUCTURE barplots

# Calculate Nei's DA distances and order the barplots according to the Nei distance tree.
# Filter the genind data and calculate Nei's DA distances.

# 1. Read cluster membership files
# ----------------------------

# Run this script from the project root; all file paths are relative.

# Individual lists for clusters 1 and 2, K2 majority assignment     
k2_c1_raw <- readLines("data/cluster_lists/STRUCTURE_3d_literature/litt_K2_cluster1_majo.txt")
k2_c2_raw <- readLines("data/cluster_lists/STRUCTURE_3d_literature/litt_K2_cluster2_majo.txt")
# Individual lists for clusters 1 and 2, K2 pure assignment 
k2_c1_pure <- readLines("data/cluster_lists/STRUCTURE_3d_literature/litt_K2_cluster1_pur.txt")
k2_c2_pure <- readLines("data/cluster_lists/STRUCTURE_3d_literature/litt_K2_cluster2_pur.txt")

# Individual lists for clusters 1, 2, and 3, K3 majority assignment 
k3_c1_raw <- readLines("data/cluster_lists/STRUCTURE_3d_literature/litt_K3_cluster1_majo.txt")
k3_c2_raw <- readLines("data/cluster_lists/STRUCTURE_3d_literature/litt_K3_cluster2_majo.txt")
k3_c3_raw <- readLines("data/cluster_lists/STRUCTURE_3d_literature/litt_K3_cluster3_majo.txt")

# Individual lists for clusters 1, 2, and 3, K3 pure assignment 
k3_c1_pure <- readLines("data/cluster_lists/STRUCTURE_3d_literature/litt_K3_cluster1_pur.txt")
k3_c2_pure <- readLines("data/cluster_lists/STRUCTURE_3d_literature/litt_K3_cluster2_pur.txt")
k3_c3_pure <- readLines("data/cluster_lists/STRUCTURE_3d_literature/litt_K3_cluster3_pur.txt")

# Individual lists for clusters 1, 2, 3, and 4, K4 majority assignment 
k4_c1_raw <- readLines("data/cluster_lists/STRUCTURE_3d_literature/litt_K4_cluster1_majo.txt")
k4_c2_raw <- readLines("data/cluster_lists/STRUCTURE_3d_literature/litt_K4_cluster2_majo.txt")
k4_c3_raw <- readLines("data/cluster_lists/STRUCTURE_3d_literature/litt_K4_cluster3_majo.txt")
k4_c4_raw <- readLines("data/cluster_lists/STRUCTURE_3d_literature/litt_K4_cluster4_majo.txt")
# Individual lists for clusters 1, 2, 3, and 4, K4 pure assignment
k4_c1_pure <- readLines("data/cluster_lists/STRUCTURE_3d_literature/litt_K4_cluster1_pur.txt")
k4_c2_pure <- readLines("data/cluster_lists/STRUCTURE_3d_literature/litt_K4_cluster2_pur.txt")
k4_c3_pure <- readLines("data/cluster_lists/STRUCTURE_3d_literature/litt_K4_cluster3_pur.txt")
k4_c4_pure <- readLines("data/cluster_lists/STRUCTURE_3d_literature/litt_K4_cluster4_pur.txt")


# Individual lists for clusters 1 to 7, K7 majority assignment 
k7_c1_raw <- readLines("data/cluster_lists/STRUCTURE_3d_literature/litt_K7_cluster1_majo.txt")
k7_c2_raw <- readLines("data/cluster_lists/STRUCTURE_3d_literature/litt_K7_cluster2_majo.txt")
k7_c3_raw <- readLines("data/cluster_lists/STRUCTURE_3d_literature/litt_K7_cluster3_majo.txt")
k7_c4_raw <- readLines("data/cluster_lists/STRUCTURE_3d_literature/litt_K7_cluster4_majo.txt")
k7_c5_raw <- readLines("data/cluster_lists/STRUCTURE_3d_literature/litt_K7_cluster5_majo.txt")
k7_c6_raw <- readLines("data/cluster_lists/STRUCTURE_3d_literature/litt_K7_cluster6_majo.txt")
k7_c7_raw <- readLines("data/cluster_lists/STRUCTURE_3d_literature/litt_K7_cluster7_majo.txt")
# Individual lists for clusters 1 to 7, K7 pure assignment 
k7_c1_pure <- readLines("data/cluster_lists/STRUCTURE_3d_literature/litt_K7_cluster1_pur.txt")
k7_c2_pure <- readLines("data/cluster_lists/STRUCTURE_3d_literature/litt_K7_cluster2_pur.txt")
k7_c3_pure <- readLines("data/cluster_lists/STRUCTURE_3d_literature/litt_K7_cluster3_pur.txt")
k7_c4_pure <- readLines("data/cluster_lists/STRUCTURE_3d_literature/litt_K7_cluster4_pur.txt")
k7_c5_pure <- readLines("data/cluster_lists/STRUCTURE_3d_literature/litt_K7_cluster5_pur.txt")
k7_c6_pure <- readLines("data/cluster_lists/STRUCTURE_3d_literature/litt_K7_cluster6_pur.txt")
k7_c7_pure <- readLines("data/cluster_lists/STRUCTURE_3d_literature/litt_K7_cluster7_pur.txt")

k7_c1_DAPC <- readLines("results/tables/cluster_lists_DAPC_K7/litt_K7_cluster1_majo_DAPC.txt")
k7_c2_DAPC <- readLines("results/tables/cluster_lists_DAPC_K7/litt_K7_cluster2_majo_DAPC.txt")
k7_c3_DAPC <- readLines("results/tables/cluster_lists_DAPC_K7/litt_K7_cluster3_majo_DAPC.txt")
k7_c4_DAPC <- readLines("results/tables/cluster_lists_DAPC_K7/litt_K7_cluster4_majo_DAPC.txt")
k7_c5_DAPC <- readLines("results/tables/cluster_lists_DAPC_K7/litt_K7_cluster5_majo_DAPC.txt")
k7_c6_DAPC <- readLines("results/tables/cluster_lists_DAPC_K7/litt_K7_cluster6_majo_DAPC.txt")
k7_c7_DAPC <- readLines("results/tables/cluster_lists_DAPC_K7/litt_K7_cluster7_majo_DAPC.txt")
library(hierfstat)

geno_mat <- snpgdsGetGeno(genofile, sample.id = sample)  # Individual-by-SNP matrix
rownames(geno_mat) <- sample_new
colnames(geno_mat) <- snps.ids

# =====================================================
# Libraries
# =====================================================
library(adegenet)
library(hierfstat)
library(ape)

# =====================================================
# Project directory
# =====================================================
# Run this script from the project root; all file paths are relative.

# =====================================================
# Genotype matrix (assumed to have already been created in the current R environment)
# =====================================================
# geno_mat <- snpgdsGetGeno(...)
# rownames(geno_mat) <- sample_new
# colnames(geno_mat) <- snps.ids


# =====================================================
# 1. LOAD ALL CLUSTER FILES INTO A STRUCTURED LIST
# =====================================================

load_cluster_files <- function(prefix, k, type) {
  lapply(1:k, function(i) {
    readLines(file.path(
      cluster_dir,
      sprintf("litt_K%d_cluster%d_%s.txt", k, i, type)
    ))
  })
}

clusters <- list(
  
  K2_raw  = load_cluster_files("litt", 2, "majo"),
  K2_pure = load_cluster_files("litt", 2, "pur"),
  
  K3_raw  = load_cluster_files("litt", 3, "majo"),
  K3_pure = load_cluster_files("litt", 3, "pur"),
  
  K4_raw  = load_cluster_files("litt", 4, "majo"),
  K4_pure = load_cluster_files("litt", 4, "pur"),
  
  K7_raw  = load_cluster_files("litt", 7, "majo"),
  K7_pure = load_cluster_files("litt", 7, "pur")
)


# =====================================================
# 2. LABEL GENERATOR
# =====================================================

make_labels <- function(name, k) {
  paste0(name, "_C", seq_len(k))
}


# =====================================================
# 3. GENIND BUILDER (SAFE VERSION)
# =====================================================

make_genind_from_clusters <- function(geno_mat, cluster_lists, cluster_labels) {
  
  inds <- unique(unlist(cluster_lists))
  inds <- intersect(inds, rownames(geno_mat))  # safety check
  
  sub_mat <- geno_mat[inds, , drop = FALSE]
  
  pop_vector <- rep(NA, length(inds))
  names(pop_vector) <- inds
  
  for (i in seq_along(cluster_lists)) {
    pop_vector[cluster_lists[[i]]] <- cluster_labels[i]
  }
  
  pop_vector <- pop_vector[inds]
  
  df2genind(
    sub_mat,
    ploidy = 2,
    type = "codom",
    ind.names = inds,
    pop = factor(pop_vector),
    ncode = 1
  )
}


# =====================================================
# 4. BOOTSTRAP NJ FUNCTION (UNCHANGED LOGIC)
# =====================================================

bootstrap_da_nj <- function(hf, nboot = 100) {
  
  da_obs <- genet.dist(hf, method = "Da")
  
  # Safety check
  if (nrow(as.matrix(da_obs)) < 3) {
    
    return(list(
      type = "da_only",
      da = da_obs
    ))
  }
  
  tree_obs <- nj(as.dist(da_obs))
  
  loci_cols <- 2:ncol(hf)
  trees_boot <- list()
  
  pb <- txtProgressBar(min = 0, max = nboot, style = 3)
  
  for (b in 1:nboot) {
    
    sampled_loci <- sample(loci_cols, replace = TRUE)
    
    hf_boot <- hf[, c(1, sampled_loci)]
    
    res <- try({
      da_boot <- genet.dist(hf_boot, method = "Da")
      nj(as.dist(da_boot))
    }, silent = TRUE)
    
    if (!inherits(res, "try-error")) {
      trees_boot[[length(trees_boot) + 1]] <- res
    }
    
    setTxtProgressBar(pb, b)
  }
  
  close(pb)
  
  bs <- prop.clades(tree_obs, trees_boot)
  
  list(
    tree = tree_obs,
    bootstrap = bs,
    boot_trees = trees_boot
  )
}


# =====================================================
# 5. FULL PIPELINE WRAPPER
# =====================================================

run_pipeline <- function(geno_mat, cluster_list, name, nboot = 100) {
  
  k <- length(cluster_list)
  
  # =========================
  # Clean cluster membership lists
  # =========================
  cluster_list <- lapply(cluster_list, function(x) {
    x[x %in% rownames(geno_mat)]
  })
  cluster_list <- cluster_list[sapply(cluster_list, length) > 0]
  
  labels <- paste0(name, "_C", seq_along(cluster_list))
  
  gi <- make_genind_from_clusters(
    geno_mat,
    cluster_list,
    labels
  )
  
  hf <- genind2hierfstat(gi)
  
  da <- genet.dist(hf, method = "Da")
  
  # =========================
  # K2 case: no tree is generated
  # =========================
  if (length(cluster_list) < 3) {
    
    message("K < 3 detected: returning Da matrix only for ", name)
    
    return(list(
      type = "da_only",
      da = da
    ))
  }
  
  # =========================
  # K >= 3 case: neighbor-joining tree and bootstrap
  # =========================
  bootstrap_da_nj(hf, nboot = nboot)
}

# =====================================================
# 6. RUN ALL CONDITIONS AUTOMATICALLY
# =====================================================

results <- list()

for (name in names(clusters)) {
  
  cat("Running:", name, "\n")
  
  results[[name]] <- run_pipeline(
    geno_mat = geno_mat,
    cluster_list = clusters[[name]],
    name = name,
    nboot = 100
  )
}


# =====================================================
# 7. PLOT FUNCTION
# =====================================================

plot_nj_bootstrap_pdf <- function(res,
                                  file = NULL,
                                  cex = 0.9,
                                  show_edge_length = TRUE,
                                  show_bootstrap = TRUE) {
  
  if (!is.null(file)) {
    pdf(file, width = 4, height = 4)
    on.exit(dev.off(), add = TRUE)
  }
  
  plot(
    res$tree,
    type = "phylogram",
    cex = cex,
    edge.width = 2,
    main = "NJ tree (Nei Da) + bootstrap"
  )
  
  if (show_edge_length) {
    edgelabels(round(res$tree$edge.length, 4),
               cex = 0.6, frame = "none", col = "darkblue")
  }
  
  if (show_bootstrap) {
    nodelabels(res$bootstrap,
               cex = 0.8, frame = "none")
  }
}


# =====================================================
# 8. EXPORT ALL PDFS 
pdf(
  "results/figures/population_structure/figure2_arbre_K3_raw_bootstraps_litt.pdf",
  width = 4,
  height = 4
)

plot_nj_bootstrap_pdf(results$K3_raw)

dev.off()

pdf(
  "results/figures/population_structure/figure2_arbre_K3_pure_bootstraps_litt.pdf",
  width = 4,
  height = 4
)

plot_nj_bootstrap_pdf(results$K3_pure)

dev.off()

pdf(
  "results/figures/population_structure/figure2_arbre_K4_raw_bootstraps_litt.pdf",
  width = 4,
  height = 4
)

plot_nj_bootstrap_pdf(results$K4_raw)

dev.off()

pdf(
  "results/figures/population_structure/figure2_arbre_K4_pure_bootstraps_litt.pdf",
  width = 4,
  height = 4
)

plot_nj_bootstrap_pdf(results$K4_pure)

dev.off()

pdf(
  "results/figures/population_structure/figure2_arbre_K7_raw_bootstraps_litt.pdf",
  width = 4,
  height = 4
)

plot_nj_bootstrap_pdf(results$K7_raw)

dev.off()

pdf(
  "results/figures/population_structure/figure2_arbre_K7_pure_bootstraps_litt.pdf",
  width = 4,
  height = 4
)

plot_nj_bootstrap_pdf(results$K7_pure)

dev.off()


# Pipeline for the Nei DA distance table
run_pipeline_table <- function(geno_mat, cluster_list, name) {
  
  # =========================
  # Clean cluster membership lists
  # =========================
  cluster_list <- lapply(cluster_list, function(x) {
    x[x %in% rownames(geno_mat)]
  })
  
  cluster_list <- cluster_list[sapply(cluster_list, length) > 0]
  
  labels <- paste0(name, "_C", seq_along(cluster_list))
  
  # =========================
  # Create the genind object
  # =========================
  gi <- make_genind_from_clusters(
    geno_mat,
    cluster_list,
    labels
  )
  
  hf <- genind2hierfstat(gi)
  
  # =========================
  # Raw Nei DA distances
  # =========================
  da <- genet.dist(
    hf,
    method = "Da"
  )
  
  # =========================
  # Convert the distance matrix to a table
  # =========================
  da_df <- as.data.frame(as.table(as.matrix(da)))
  
  colnames(da_df) <- c(
    "cluster1",
    "cluster2",
    "Da"
  )
  
  # =========================
  # Remove diagonal entries
  # =========================
  da_df <- da_df %>%
    dplyr::filter(cluster1 != cluster2)
  
  # =========================
  # Remove symmetric duplicates
  # =========================
  da_df <- da_df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      pair = paste(sort(c(cluster1, cluster2)), collapse = "_vs_")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(pair, .keep_all = TRUE) %>%
    dplyr::select(-pair)
  
  return(da_df)
}

# =====================================================
# RUN ALL CONDITIONS
# =====================================================

results_table <- list()

for (name in names(clusters)) {
  
  cat("Running:", name, "\n")
  
  results_table[[name]] <- run_pipeline_table(
    geno_mat = geno_mat,
    cluster_list = clusters[[name]],
    name = name
  )
}

# =====================================================
# COMBINE ALL RESULTS
# =====================================================

all_da <- dplyr::bind_rows(
  results_table,
  .id = "condition"
)

print(all_da)
## BARPLOTS

library(dplyr)
library(tidyr)
library(ggplot2)

# =========================
# Colors
# =========================

fill_values <- c(
  "1_pure" = "#1f77b4",  # blue
  "2_pure" = "#d62728",  # red
  "3_pure" = "#2ca02c",
  
  "1_admixed" = "steelblue3",
  "2_admixed" = "tomato3",
  "3_admixed" = "palegreen3"
)

# =========================
# Order of major clusters
# =========================

cluster_order <- c("1", "2", "3")

# =========================
# Sort individuals
# =========================

metadata_ordered <- metadata %>%
  
  mutate(
    K3_cluster_majoritaire = factor(
      K3_cluster_majoritaire,
      levels = cluster_order
    )
  ) %>%
  
  mutate(
    major_prop = case_when(
      K3_cluster_majoritaire == "1" ~ K3_cluster1,
      K3_cluster_majoritaire == "2" ~ K3_cluster2,
      K3_cluster_majoritaire == "3" ~ K3_cluster3
    )
  ) %>%
  
  arrange(
    K3_cluster_majoritaire,
    factor(K3_statut,
           levels = c("pure", "admixed")),
    desc(major_prop) 
  ) %>%
  mutate(
    ind = row_number()
  )

# =========================
# Convert to long format
# =========================

df_long <- metadata_ordered %>%
  
  dplyr::select(
    ind,
    K3_statut,
    K3_cluster1,
    K3_cluster2,
    K3_cluster3
  ) %>%
  
  pivot_longer(
    cols = starts_with("K3_cluster"),
    names_to = "cluster",
    values_to = "prop"
  ) %>%
  
  mutate(
    cluster = gsub("K3_cluster", "", cluster)
  )

# =========================
# Fixed segment order
# =========================

df_long$cluster <- factor(
  df_long$cluster,
  levels = c("1", "2", "3")
)

# =========================
# Independent color variable
# =========================

df_long <- df_long %>%
  
  mutate(
    fill_group = paste0(
      cluster,
      "_",
      ifelse(K3_statut == "admixed", "admixed", "pure")
    )
  )

# =========================
# Vertical STRUCTURE barplot
# =========================

p_k3 <- ggplot(
  df_long,
  aes(
    y = factor(ind, levels = unique(ind)),
    x = prop,
    fill = fill_group,
    order = cluster
  )
) +
  
  geom_bar(
    stat = "identity",
    width = 1,
    position = position_stack(reverse = TRUE),
    color = NA,
    linewidth = 0.05
  ) +
  
  scale_fill_manual(
    values = fill_values
  ) +
  
  scale_x_continuous(
    expand = c(0, 0)
  ) +
  scale_y_discrete(limits = rev) +
  
  coord_cartesian(
    xlim = c(0, 1)
  ) +
  
  theme_classic() +
  
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    panel.grid = element_blank(),
    axis.line = element_blank()
  )

p_k3

ggsave("results/figures/population_structure/figure2_barplot_K3.pdf", p_k3,
       width = 2,
       height = 8,
       device = cairo_pdf,
       units = "cm",
       dpi = 600)

# =========================
# Retrieve the final K4 order
# =========================

final_order <- metadata_ordered$ind

# To retrieve the original sample names:
# final_order_ids <- df2_ordered$sample_id

# =========================
# Apply exactly
# the same order to K2
# =========================

df2_k2 <- metadata_ordered %>%
  
  dplyr::select(
    ind,
    K2_statut,
    K2_cluster1,
    K2_cluster2
  )

# =========================
# Convert to long format K2
# =========================

df2_k2_long <- df2_k2 %>%
  
  pivot_longer(
    cols = starts_with("K2_cluster"),
    names_to = "cluster",
    values_to = "prop"
  ) %>%
  
  mutate(
    cluster = gsub("K2_cluster", "", cluster)
  )

# =========================
# Fixed segment order
# =========================

df2_k2_long$cluster <- factor(
  df2_k2_long$cluster,
  levels = c("1", "2")
)

# =========================
# Colors K2
# =========================

fill_values_k2 <- c(
  "1_pure" = "blue",
  "2_pure" = "red",
  
  "1_admixed" = "steelblue3",
  "2_admixed" = "tomato3"
)

# =========================
# Color variable
# =========================

df2_k2_long <- df2_k2_long %>%
  
  mutate(
    fill_group = paste0(
      cluster,
      "_",
      ifelse(K2_statut == "admixed", "admixed", "pure")
    )
  )

# =========================
# Plot K2
# =========================

p_k2 <- ggplot(
  df2_k2_long,
  aes(
    y = factor(ind, levels = final_order),
    x = prop,
    fill = fill_group,
    order = cluster
  )
) +
  
  geom_bar(
    stat = "identity",
    width = 1,
    position = position_stack(reverse = TRUE),
    color = NA,
    linewidth = 0.001
  ) +
  
  scale_fill_manual(
    values = fill_values_k2
  ) +
  
  scale_x_continuous(
    expand = c(0, 0)
  ) +
  scale_y_discrete(limits = rev) +
  coord_cartesian(
    xlim = c(0, 1)
  ) +
  
  theme_classic() +
  
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    panel.grid = element_blank(),
    axis.line = element_blank()
  )

p_k2

ggsave("results/figures/population_structure/figure2_barplot_K2.pdf", p_k2,
       width = 2,
       height = 8,
       device = cairo_pdf,
       units = "cm",
       dpi = 600)

df4_k4 <- metadata_ordered %>%
  
  dplyr::select(
    ind,
    K4_statut,
    K4_cluster1,
    K4_cluster2,
    K4_cluster3,
    K4_cluster4
  )

# =========================
# Convert to long format K4
# =========================

df4_k4_long <- df4_k4 %>%
  
  pivot_longer(
    cols = starts_with("K4_cluster"),
    names_to = "cluster",
    values_to = "prop"
  ) %>%
  
  mutate(
    cluster = gsub("K4_cluster", "", cluster)
  )

# =========================
# Fixed segment order
# =========================

df4_k4_long$cluster <- factor(
  df4_k4_long$cluster,
  levels = c("1", "2", "3", "4")
)

# =========================
# Colors K4
# =========================

fill_values_k4 <- c(
  "1_pure" = "blue",  # blue
  "2_pure" = "red",  # red
  "3_pure" = "#2ca02c",  # green
  "4_pure" = "orange",   # orange
  
  "1_admixed" = "steelblue3",
  "2_admixed" = "tomato3",
  "3_admixed" = "palegreen3",
  "4_admixed" = "peachpuff1"
)

# =========================
# Color variable
# =========================

df4_k4_long <- df4_k4_long %>%
  
  mutate(
    fill_group = paste0(
      cluster,
      "_",
      ifelse(K4_statut == "admixed", "admixed", "pure")
    )
  )

# =========================
# Plot K4
# =========================

p_k4 <- ggplot(
  df4_k4_long,
  aes(
    y = factor(ind, levels = final_order),
    x = prop,
    fill = fill_group,
    order = cluster
  )
) +
  
  geom_bar(
    stat = "identity",
    width = 1,
    position = position_stack(reverse = TRUE),
    color = NA,
    linewidth = 0.001
  ) +
  
  scale_fill_manual(
    values = fill_values_k4
  ) +
  
  scale_x_continuous(
    expand = c(0, 0)
  ) +
  scale_y_discrete(limits = rev) +
  coord_cartesian(
    xlim = c(0, 1)
  ) +
  
  theme_classic() +
  
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    panel.grid = element_blank(),
    axis.line = element_blank()
  )

p_k4

ggsave("results/figures/population_structure/figure2_barplot_K4.pdf", p_k4,
       width = 2,
       height = 8,
       device = cairo_pdf,
       units = "cm",
       dpi = 600)

df7_k7 <- metadata_ordered %>%
  
  dplyr::select(
    ind,
    K7_statut,
    K7_cluster1,
    K7_cluster2,
    K7_cluster3,
    K7_cluster4,
    K7_cluster5,
    K7_cluster6,
    K7_cluster7
  )

# =========================
# Convert to long format K7
# =========================

df7_k7_long <- df7_k7 %>%
  
  pivot_longer(
    cols = starts_with("K7_cluster"),
    names_to = "cluster",
    values_to = "prop"
  ) %>%
  
  mutate(
    cluster = gsub("K7_cluster", "", cluster)
  )

# =========================
# Fixed segment order
# =========================

df7_k7_long$cluster <- factor(
  df7_k7_long$cluster,
  levels = c("1", "2", "3", "4", "5", "6", "7")
)

# =========================
# Colors K7
# =========================

fill_values_k7 <- c(
  "1_pure" = "blue",  # blue
  "2_pure" = "red",  # red
  "3_pure" = "#2ca02c",  # green
  "4_pure" = "orange",   # orange
  "5_pure" = "#9467bd",  # purple
  "6_pure" = "#8c564b",  # brown
  "7_pure" = "#e377c2",   # pink
  
  "1_admixed" = "steelblue3",
  "2_admixed" = "tomato3",
  "3_admixed" = "palegreen3",
  "4_admixed" = "peachpuff1",
  "5_admixed" = "plum3",  # purple
  "6_admixed" = "burlywood3",  # brown
  "7_admixed" = "pink3"   # pink
)

# =========================
# Color variable
# =========================

df7_k7_long <- df7_k7_long %>%
  
  mutate(
    fill_group = paste0(
      cluster,
      "_",
      ifelse(K7_statut == "admixed", "admixed", "pure")
    )
  )

# =========================
# Plot K7
# =========================

p_k7 <- ggplot(
  df7_k7_long,
  aes(
    y = factor(ind, levels = final_order),
    x = prop,
    fill = fill_group,
    order = cluster
  )
) +
  
  geom_bar(
    stat = "identity",
    width = 1,
    position = position_stack(reverse = TRUE),
    color = NA,
    linewidth = 0.001
  ) +
  
  scale_fill_manual(
    values = fill_values_k7
  ) +
  
  scale_x_continuous(
    expand = c(0, 0)
  ) +
  scale_y_discrete(limits = rev) +
  coord_cartesian(
    xlim = c(0, 1)
  ) +
  
  theme_classic() +
  
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    panel.grid = element_blank(),
    axis.line = element_blank()
  )

p_k7

ggsave("results/figures/population_structure/figure2_barplot_K7.pdf", p_k7,
       width = 2,
       height = 8,
       device = cairo_pdf,
       units = "cm",
       dpi = 600)



### Nei DA analysis for K7 DAPC clusters 


# =====================================================
# 1. LOAD ALL CLUSTER FILES INTO A STRUCTURED LIST
# =====================================================

load_cluster_files <- function(prefix, k, type) {
  lapply(1:k, function(i) {
    readLines(file.path(
      dapc_out_dir,
      sprintf("litt_K%d_cluster%d_%s_DAPC.txt", k, i, type)
    ))
  })
}

clusters <- list(
  
  K7_raw  = load_cluster_files("litt", 7, "majo")
)


# =====================================================
# 2. LABEL GENERATOR
# =====================================================

make_labels <- function(name, k) {
  paste0(name, "_C", seq_len(k))
}


# =====================================================
# 3. GENIND BUILDER (SAFE VERSION)
# =====================================================

make_genind_from_clusters <- function(geno_mat, cluster_lists, cluster_labels) {
  
  inds <- unique(unlist(cluster_lists))
  inds <- intersect(inds, rownames(geno_mat))  # safety check
  
  sub_mat <- geno_mat[inds, , drop = FALSE]
  
  pop_vector <- rep(NA, length(inds))
  names(pop_vector) <- inds
  
  for (i in seq_along(cluster_lists)) {
    pop_vector[cluster_lists[[i]]] <- cluster_labels[i]
  }
  
  pop_vector <- pop_vector[inds]
  
  df2genind(
    sub_mat,
    ploidy = 2,
    type = "codom",
    ind.names = inds,
    pop = factor(pop_vector),
    ncode = 1
  )
}


# =====================================================
# 4. BOOTSTRAP NJ FUNCTION (UNCHANGED LOGIC)
# =====================================================

bootstrap_da_nj <- function(hf, nboot = 100) {
  
  da_obs <- genet.dist(hf, method = "Da")
  
  # Safety check
  if (nrow(as.matrix(da_obs)) < 3) {
    
    return(list(
      type = "da_only",
      da = da_obs
    ))
  }
  
  tree_obs <- nj(as.dist(da_obs))
  
  loci_cols <- 2:ncol(hf)
  trees_boot <- list()
  
  pb <- txtProgressBar(min = 0, max = nboot, style = 3)
  
  for (b in 1:nboot) {
    
    sampled_loci <- sample(loci_cols, replace = TRUE)
    
    hf_boot <- hf[, c(1, sampled_loci)]
    
    res <- try({
      da_boot <- genet.dist(hf_boot, method = "Da")
      nj(as.dist(da_boot))
    }, silent = TRUE)
    
    if (!inherits(res, "try-error")) {
      trees_boot[[length(trees_boot) + 1]] <- res
    }
    
    setTxtProgressBar(pb, b)
  }
  
  close(pb)
  
  bs <- prop.clades(tree_obs, trees_boot)
  
  list(
    tree = tree_obs,
    bootstrap = bs,
    boot_trees = trees_boot
  )
}


# =====================================================
# 5. FULL PIPELINE WRAPPER
# =====================================================

run_pipeline <- function(geno_mat, cluster_list, name, nboot = 100) {
  
  k <- length(cluster_list)
  
  # =========================
  # Clean cluster membership lists
  # =========================
  cluster_list <- lapply(cluster_list, function(x) {
    x[x %in% rownames(geno_mat)]
  })
  cluster_list <- cluster_list[sapply(cluster_list, length) > 0]
  
  labels <- paste0(name, "_C", seq_along(cluster_list))
  
  gi <- make_genind_from_clusters(
    geno_mat,
    cluster_list,
    labels
  )
  
  hf <- genind2hierfstat(gi)
  
  da <- genet.dist(hf, method = "Da")
  
  # =========================
  # K2 case: no tree is generated
  # =========================
  if (length(cluster_list) < 3) {
    
    message("K < 3 detected: returning Da matrix only for ", name)
    
    return(list(
      type = "da_only",
      da = da
    ))
  }
  
  # =========================
  # K >= 3 case: neighbor-joining tree and bootstrap
  # =========================
  bootstrap_da_nj(hf, nboot = nboot)
}

# =====================================================
# 6. RUN ALL CONDITIONS AUTOMATICALLY
# =====================================================

results <- list()

for (name in names(clusters)) {
  
  cat("Running:", name, "\n")
  
  results[[name]] <- run_pipeline(
    geno_mat = geno_mat,
    cluster_list = clusters[[name]],
    name = name,
    nboot = 100
  )
}


# =====================================================
# 7. PLOT FUNCTION
# =====================================================

plot_nj_bootstrap_pdf <- function(res,
                                  file = NULL,
                                  cex = 0.9,
                                  show_edge_length = TRUE,
                                  show_bootstrap = TRUE) {
  
  if (!is.null(file)) {
    pdf(file, width = 4, height = 4)
    on.exit(dev.off(), add = TRUE)
  }
  
  plot(
    res$tree,
    type = "phylogram",
    cex = cex,
    edge.width = 2,
    main = "NJ tree (Nei Da) + bootstrap"
  )
  
  if (show_edge_length) {
    edgelabels(round(res$tree$edge.length, 4),
               cex = 0.6, frame = "none", col = "darkblue")
  }
  
  if (show_bootstrap) {
    nodelabels(res$bootstrap,
               cex = 0.8, frame = "none")
  }
}


# =====================================================
# 8. EXPORT ALL PDFS 

pdf(
  "results/figures/population_structure/figure2_arbre_K7_DAPC_raw_bootstraps_litt.pdf",
  width = 4,
  height = 4
)

plot_nj_bootstrap_pdf(results$K7_raw)

dev.off()


# Pipeline for the Nei DA distance table
run_pipeline_table <- function(geno_mat, cluster_list, name) {
  
  # =========================
  # Clean cluster membership lists
  # =========================
  cluster_list <- lapply(cluster_list, function(x) {
    x[x %in% rownames(geno_mat)]
  })
  
  cluster_list <- cluster_list[sapply(cluster_list, length) > 0]
  
  labels <- paste0(name, "_C", seq_along(cluster_list))
  
  # =========================
  # Create the genind object
  # =========================
  gi <- make_genind_from_clusters(
    geno_mat,
    cluster_list,
    labels
  )
  
  hf <- genind2hierfstat(gi)
  
  # =========================
  # Raw Nei DA distances
  # =========================
  da <- genet.dist(
    hf,
    method = "Da"
  )
  
  # =========================
  # Convert the distance matrix to a table
  # =========================
  da_df <- as.data.frame(as.table(as.matrix(da)))
  
  colnames(da_df) <- c(
    "cluster1",
    "cluster2",
    "Da"
  )
  
  # =========================
  # Remove diagonal entries
  # =========================
  da_df <- da_df %>%
    dplyr::filter(cluster1 != cluster2)
  
  # =========================
  # Remove symmetric duplicates
  # =========================
  da_df <- da_df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      pair = paste(sort(c(cluster1, cluster2)), collapse = "_vs_")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(pair, .keep_all = TRUE) %>%
    dplyr::select(-pair)
  
  return(da_df)
}

# =====================================================
# RUN ALL CONDITIONS
# =====================================================

results_table <- list()

for (name in names(clusters)) {
  
  cat("Running:", name, "\n")
  
  results_table[[name]] <- run_pipeline_table(
    geno_mat = geno_mat,
    cluster_list = clusters[[name]],
    name = name
  )
}

# =====================================================
# COMBINE ALL RESULTS
# =====================================================

all_da <- dplyr::bind_rows(
  results_table,
  .id = "condition"
)

print(all_da)
