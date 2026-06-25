############################################################
# 0. LIBRARIES
############################################################

library(tidyverse)
library(readxl)
library(SNPRelate)
library(adegenet)
library(poppr)
library(sf)
library(terra)
library(elevatr)
library(ggrepel)
library(ggnewscale)
library(rnaturalearth)

############################################################
# 1. PROJECT PATHS (REPRODUCIBLE STRUCTURE)
############################################################

root <- "."
data_path <- file.path(root, "data")
fig_path <- file.path(root, "figures")

dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)

save_fig <- function(plot, name, w = 8, h = 8) {
  ggsave(file.path(fig_path, name),
         plot = plot,
         width = w,
         height = h,
         units = "cm",
         device = cairo_pdf,
         dpi = 600)
}

############################################################
# 2. METADATA AND SAMPLE FILTERING
############################################################

metadata <- read_excel(file.path(data_path, "metadata.xlsx"))

flagged_samples <- metadata %>%
  filter(`Flag-weirdposition` == "Y") %>%
  pull(`numéro échantillon`)

############################################################
# 3. GENETIC DATA IMPORT (VCF → GDS)
############################################################

vcf.fn <- file.path(data_path, "data.vcf")
gds.fn <- file.path(data_path, "data.gds")

snpgdsVCF2GDS(vcf.fn, gds.fn, method = "biallelic.only")
genofile <- snpgdsOpen(gds.fn)

popk4 <- read_tsv(file.path(data_path, "structure.tsv")) %>%
  mutate(sample_id = as.character(sample_id)) %>%
  filter(!sample_id %in% flagged_samples)

############################################################
# 4. PCA ANALYSIS
############################################################

run_pca <- function(genofile, popk4, n_pcs = 20) {
  
  sample <- read.gdsn(index.gdsn(genofile, "sample.id"))
  
  sample_short <- gsub("\\D", "",
                       sapply(strsplit(sub("\\..*", "", sample), "_"), `[`, 3))
  
  pca <- snpgdsPCA(genofile, sample.id = sample, eigen.cnt = n_pcs)
  
  pca_df <- tibble(
    sample.id = pca$sample.id,
    PC1 = pca$eigenvect[,1],
    PC2 = pca$eigenvect[,2],
    PC3 = pca$eigenvect[,3],
    PC4 = pca$eigenvect[,4]
  ) %>%
    mutate(
      sample_short = gsub("\\D", "",
                          sapply(strsplit(sub("\\..*", "", sample.id), "_"), `[`, 3))
    ) %>%
    inner_join(popk4, by = c("sample_short" = "sample_id"))
  
  list(df = pca_df, var = pca$varprop * 100)
}

pca_res <- run_pca(genofile, popk4)

# Add alpha transparency based on ancestry purity
pca_res$df <- pca_res$df %>%
  mutate(
    K2_alpha = ifelse(K2_statut == "pure", 1, 0.4),
    K4_alpha = ifelse(K4_statut == "pure", 1, 0.4)
  )

############################################################
# 5. PCA PLOTS
############################################################

theme_paper <- function() {
  theme_classic(base_size = 8) +
    theme(
      legend.position = "none",
      axis.title = element_blank()
    )
}

plot_pca <- function(df, var, cluster, alpha, colors) {
  
  ggplot(df, aes(PC1, PC2)) +
    geom_point(aes_string(color = cluster, alpha = alpha), size = 1) +
    stat_ellipse(aes_string(group = cluster),
                 linewidth = 0.3,
                 linetype = "dashed") +
    scale_color_manual(values = colors) +
    scale_alpha_identity() +
    labs(
      x = paste0("PC1 (", round(var[1],1), "%)"),
      y = paste0("PC2 (", round(var[2],1), "%)")
    ) +
    coord_equal() +
    theme_paper()
}

# K2 PCA
p_k2 <- plot_pca(
  pca_res$df, pca_res$var,
  "K2_cluster_majoritaire", "K2_alpha",
  c("1"="blue","2"="red")
)

# K4 PCA
p_k4 <- plot_pca(
  pca_res$df, pca_res$var,
  "K4_cluster_majoritaire", "K4_alpha",
  c("1"="purple","2"="green","3"="red","4"="blue")
)

save_fig(p_k2, "pca_k2.pdf")
save_fig(p_k4, "pca_k4.pdf")

############################################################
# 6. ENVIRONMENTAL DATA (CLIMATE EXTRACTION INPUT)
############################################################

# This dataset contains environmental variables used for biogeographic analyses
biogeo <- read_tsv(file.path(data_path, "enviroGN.txt"))

# Keep only samples present in genetic dataset
biogeo <- biogeo %>%
  filter(numero %in% pca_res$df$sample_short)

############################################################
# 7. SITE GEOMETRY AND ELEVATION EXTRACTION
############################################################

# Convert sampling sites to spatial object
sites_sf <- st_as_sf(
  data.frame(
    longitude = biogeo$longitude,
    latitude = biogeo$latitude
  ),
  coords = c("longitude", "latitude"),
  crs = 4326
)

# Extract elevation raster
elev <- get_elev_raster(sites_sf, z = 6, clip = "locations")

elev_df <- as.data.frame(elev[[1]], xy = TRUE)
colnames(elev_df) <- c("x", "y", "elevation")

############################################################
# 8. WHITTAKER BIOME PLOT
############################################################

# Base layout for Whittaker biome space
whittaker_base <- function() {
  ggplot() +
    labs(
      x = "Mean annual temperature (°C)",
      y = "Annual precipitation"
    )
}

plot_whittaker <- function(df, cluster, alpha, colors) {
  
  whittaker_base() +
    geom_point(
      data = df,
      aes(
        x = (tmaxM + tminM) / 2,
        y = pptM,
        fill = .data[[cluster]],
        alpha = .data[[alpha]]
      ),
      shape = 21,
      size = 1
    ) +
    scale_fill_manual(values = colors) +
    scale_alpha_identity() +
    theme_classic(base_size = 8) +
    theme(legend.position = "none")
}

# K2 Whittaker plot
p_whitt_k2 <- plot_whittaker(
  biogeo,
  "K2_cluster_majoritaire",
  "K2_statut",
  c("1"="blue","2"="red")
)

# K4 Whittaker plot
p_whitt_k4 <- plot_whittaker(
  biogeo,
  "K4_cluster_majoritaire",
  "K4_statut",
  c("1"="purple","2"="green","3"="red","4"="blue")
)

save_fig(p_whitt_k2, "whittaker_k2.pdf")
save_fig(p_whitt_k4, "whittaker_k4.pdf")

############################################################
# 9. MAPS (INDIA + SAMPLING SITES)
############################################################

# Load spatial boundaries
india <- ne_countries(country = "india", returnclass = "sf")
india2 <- ne_states(country = "india", returnclass = "sf")

plot_map <- function(df, cluster, colors) {
  
  ggplot() +
    geom_sf(data = india, fill = "white") +
    geom_sf(data = india2, fill = NA) +
    geom_point(
      data = df,
      aes(longitude, latitude,
          fill = .data[[cluster]]),
      shape = 21,
      size = 2
    ) +
    scale_fill_manual(values = colors) +
    theme_classic(base_size = 8) +
    theme(
      legend.position = "none",
      axis.title = element_blank()
    )
}

# K2 map
p_map_k2 <- plot_map(
  pca_res$df,
  "K2_cluster_majoritaire",
  c("1"="blue","2"="red")
)

# K4 map
p_map_k4 <- plot_map(
  pca_res$df,
  "K4_cluster_majoritaire",
  c("1"="purple","2"="green","3"="red","4"="blue")
)

save_fig(p_map_k2, "map_k2.pdf", 17, 13)
save_fig(p_map_k4, "map_k4.pdf", 17, 13)
