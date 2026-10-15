# =====================================================
# 1. LIBRARIES AND PROJECT PATHS
# =====================================================

library(tidyverse)
library(ggplot2)
library(readxl)
library(dplyr)
library(SNPRelate)
library(adegenet)
library(poppr)
library(readr)
library(Boruta)
library(Amelia)
library(randomForest)
library(sp)
library(pdp)
library(sf)
library(rnaturalearth)
library(ggspatial)
library(AOI)
library(climateR)
library(raster)
library(rasterVis)
library(RColorBrewer)
library(elevatr)
library(terra)
library(ggnewscale)
library(geodata)
library(hierfstat)
library(ape)
library(tidyr)
library(ggrepel)

# =====================================================
# PROJECT DIRECTORY STRUCTURE
# =====================================================

# Set the working directory to the project root before running this script
project_dir <- "."

# Input files
vcf.fn <- file.path(
  project_dir, "data", "genotypes",
  "pruned_forPCA_STRUCTURE_maf001_miss08_qual_noindel_biallelic_renamed_mind08_hz06.vcf"
)

gds.fn <- file.path(
  project_dir, "data", "genotypes",
  "pruned_forPCA_STRUCTURE_maf001_miss08_qual_noindel_biallelic_renamed_mind08_hz06.gds"
)

metadata.fn <- file.path(
  project_dir, "data", "metadata",
  "3D_metadata_STRUCTURE_DAPC_results.tsv"
)

environment.fn <- file.path(
  project_dir, "data", "environment",
  "1_enviroGN.txt"
)

cluster_dir <- file.path(
  project_dir, "data", "clusters",
  "F_listes_clusters_STRUCTURE_3d_litt"
)

# Output directories
figure_dir <- file.path(
  project_dir, "results", "figures", "vector"
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# =====================================================
# 2. LOAD METADATA
# =====================================================

# Replace with CLUMPAK results if necessary

popk4 <- read_tsv(metadata.fn)
popk4$sample_id <- as.character(popk4$sample_id)

# =====================================================
# 3. GENOME-WIDE PCA
# =====================================================

# Convert the VCF file to GDS format
# This creates the GDS file if it does not already exist.
# If the GDS file already exists, ensure it corresponds to the current VCF.

snpgdsVCF2GDS(
  vcf.fn,
  gds.fn,
  method = "biallelic.only"
)

genofile <- snpgdsOpen(gds.fn)

sample <- read.gdsn(
  index.gdsn(genofile, "sample.id")
)

# Extract sample identifiers
sample_new <- gsub(
  "\\D",
  "",
  sapply(
    strsplit(sub("\\..*", "", sample), "_"),
    `[`,
    3
  )
)

# Retain samples present in the metadata
keep <- sample_new %in% popk4$sample_id
sample_filt <- sample[keep]

snps.ids <- read.gdsn(
  index.gdsn(genofile, "snp.id")
)

# Compute the first 20 principal components
pca <- snpgdsPCA(
  genofile,
  autosome.only = FALSE,
  eigen.cnt = 20,
  sample.id = sample_filt
)

pca$eigenval <- pca$eigenval[1:20]
pca$varprop <- pca$eigenval / sum(pca$eigenval)

pc.percent <- pca$varprop * 100

# Extract the first four principal components
pca_df <- data.frame(
  sample.id = pca$sample.id,
  PC1 = pca$eigenvect[, 1],
  PC2 = pca$eigenvect[, 2],
  PC3 = pca$eigenvect[, 3],
  PC4 = pca$eigenvect[, 4]
)

# Match PCA sample identifiers to metadata
pca_df$sample_short <- sapply(
  strsplit(sub("\\..*", "", pca_df$sample.id), "_"),
  `[`,
  3
)

pca_df$sample_short <- gsub(
  "\\D",
  "",
  pca_df$sample_short
)

pca_df <- inner_join(
  pca_df,
  popk4,
  by = c("sample_short" = "sample_id")
)

pca_df$sample_short <- sub(
  "[^0-9].*$",
  "",
  pca_df$sample_short
)

pca_df$sample_short <- as.double(
  pca_df$sample_short
)

# Set cluster assignments as factors
pca_df$K2_cluster_majoritaire <- factor(
  pca_df$K2_cluster_majoritaire,
  levels = 1:2,
  labels = c("1", "2")
)

pca_df$K4_cluster_majoritaire <- factor(
  pca_df$K4_cluster_majoritaire,
  levels = 1:4,
  labels = c("1", "2", "3", "4")
)

# Percentage of variance explained
pve <- data.frame(
  pve = pca$varprop * 100
)

# Set point transparency according to assignment purity
pca_df$K2_alpha <- ifelse(
  pca_df$K2_statut == "pure",
  0.9,
  0.3
)

pca_df$K4_alpha <- ifelse(
  pca_df$K4_statut == "pure",
  0.9,
  0.3
)

# =====================================================
# 4. PUBLICATION THEME
# =====================================================

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

# =====================================================
# 5. PCA PLOT: K2
# =====================================================

# Flag individuals assigned to DAPC K5 cluster 3 as outliers
pca_df <- pca_df %>%
  mutate(
    K2_cluster_majoritaire = if_else(
      DAPC_K5_cluster_majoritaire == 3,
      "outlier",
      as.character(K2_cluster_majoritaire)
    )
  )

pca_K2 <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    color = K2_cluster_majoritaire
  )
) +
  geom_point(
    alpha = pca_df$K2_alpha,
    size = 1.0
  ) +

  # Add confidence ellipses
  stat_ellipse(
    aes(group = K2_cluster_majoritaire),
    type = "norm",
    linetype = "dashed",
    linewidth = 0.3,
    level = 0.68
  ) +

  labs(
    x = paste0("PC1 (", round(pc.percent[1], 1), "%)"),
    y = paste0("PC2 (", round(pc.percent[2], 1), "%)")
  ) +

  scale_color_manual(
    values = c(
      "1" = "blue",
      "2" = "red",
      "outlier" = "yellow"
    )
  ) +

  coord_fixed() +
  theme_classic(base_size = 8) +
  theme_paper() +
  theme(
    legend.position = "none"
  )

print(pca_K2)

# Export as a vector PDF
ggsave(
  file.path(figure_dir, "figure1_pca_K2.pdf"),
  pca_K2,
  width = 8,
  height = 8,
  device = cairo_pdf,
  units = "cm",
  dpi = 600
)

# =====================================================
# 6. PCA PLOT: K4
# =====================================================

pca_K4 <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    color = K4_cluster_majoritaire
  )
) +
  geom_point(
    alpha = pca_df$K4_alpha,
    size = 1.0
  ) +

  # Add confidence ellipses
  stat_ellipse(
    aes(group = K4_cluster_majoritaire),
    type = "norm",
    linetype = "dashed",
    linewidth = 0.3,
    level = 0.68
  ) +

  labs(
    x = paste0("PC1 (", round(pc.percent[1], 1), "%)"),
    y = paste0("PC2 (", round(pc.percent[2], 1), "%)")
  ) +

  scale_color_manual(
    values = c(
      "1" = "purple",
      "2" = "green",
      "3" = "red",
      "4" = "blue"
    )
  ) +

  coord_fixed() +
  theme_classic(base_size = 8) +
  theme_paper() +
  theme(
    legend.position = "none"
  )

print(pca_K4)

ggsave(
  file.path(figure_dir, "figure1_pca_K4.pdf"),
  pca_K4,
  width = 8,
  height = 8,
  device = cairo_pdf,
  units = "cm",
  dpi = 600
)

# =====================================================
# 7. ENVIRONMENTAL DATA AND BIOGEOGRAPHICAL METADATA
# =====================================================

data_clim <- read_tsv(environment.fn)

data_clim_filtered <- data_clim[
  data_clim$numero %in% pca_df$sample_short,
]

df2 <- merge(
  data_clim_filtered,
  popk4,
  by.x = "numero",
  by.y = "sample_id",
  all.x = TRUE
)

# Optional: export the combined environmental and population-genetic data
# write_tsv(
#   df2,
#   file.path(project_dir, "results", "final3d_Q_biogeo_matrix.tsv")
# )

# =====================================================
# 8. PACKAGE MANAGEMENT
# =====================================================

Set_pack <- function(packages) {

  # Check whether the package is installed
  if (!require(packages, character.only = TRUE)) {

    # Install the package if it is missing
    install.packages(
      packages,
      dependencies = TRUE
    )

    # Load the package
    library(packages, character.only = TRUE)

  } else {

    library(packages, character.only = TRUE)
  }
}

# Load packages required for mapping and environmental analyses
Set_pack("tidyverse")
Set_pack("readxl")
Set_pack("openxlsx")
Set_pack("sf")
Set_pack("sp")
Set_pack("ggplot2")
Set_pack("rnaturalearth")
Set_pack("ggspatial")
Set_pack("AOI")
Set_pack("climateR")
Set_pack("raster")
Set_pack("rasterVis")
Set_pack("dplyr")
Set_pack("plotbiomes")
Set_pack("RColorBrewer")
Set_pack("elevatr")
Set_pack("terra")

# Install development packages separately if required:
# remotes::install_github("mikejohnson51/AOI")
# remotes::install_github("mikejohnson51/climateR")
# remotes::install_github("valentinitnelav/plotbiomes")

# =====================================================
# 9. STUDY SITES AND ELEVATION DATA
# =====================================================

Tbl_Site <- df2

na_count <- sum(is.na(Tbl_Site$latitude))

Tbl_Site$latitude <- as.numeric(
  Tbl_Site$latitude
)

Tbl_Site$longitude <- as.numeric(
  Tbl_Site$longitude
)

# Load country and state boundaries
india <- ne_countries(
  country = "india",
  returnclass = "sf"
)

india2 <- ne_states(
  country = "india",
  returnclass = "sf"
)["geometry"]

# Convert sampling sites to spatial objects
sites_sp <- SpatialPointsDataFrame(
  coords = Tbl_Site[, c("longitude", "latitude")],
  data = Tbl_Site,
  proj4string = CRS("+proj=longlat +datum=WGS84")
)

sites_sf <- st_as_sf(sites_sp)

# Optional climate data extraction
# These operations may take some time.

# tmax_grid <- getTerraClimNormals(
#   sites_sf,
#   varname = "tmax",
#   scenario = "19812010",
#   month = 1:12
# )
#
# tmin_grid <- getTerraClimNormals(
#   sites_sf,
#   varname = "tmin",
#   scenario = "19812010",
#   month = 1:12
# )
#
# ppt_grid <- getTerraClimNormals(
#   sites_sf,
#   varname = "ppt",
#   scenario = "19812010",
#   month = 1:12
# )
#
# tmax <- terra::extract(tmax_grid[[1]], sites_sf) %>%
#   mutate(
#     tmaxM = rowMeans(
#       dplyr::select(., starts_with("tmax")),
#       na.rm = TRUE
#     )
#   )
#
# tmin <- terra::extract(tmin_grid[[1]], sites_sf) %>%
#   mutate(
#     tminM = rowMeans(
#       dplyr::select(., starts_with("tmin")),
#       na.rm = TRUE
#     )
#   )
#
# ppt <- terra::extract(ppt_grid[[1]], sites_sf) %>%
#   mutate(
#     pptM = rowMeans(
#       dplyr::select(., starts_with("ppt")),
#       na.rm = TRUE
#     )
#   )

# Extract elevation data for sampling sites
elev.grid <- get_elev_raster(
  sites_sf,
  z = 6
)

elev <- terra::extract(
  elev.grid[[1]],
  sites_sf
)

# Define the study area bounding box
aoi <- st_as_sf(
  st_as_sfc(
    st_bbox(
      c(
        xmin = 74,
        xmax = 82,
        ymin = 28,
        ymax = 34
      ),
      crs = st_crs(4326)
    )
  )
)

# Retrieve elevation data for the study area
# Zoom level controls spatial resolution
elev <- get_elev_raster(
  aoi,
  z = 6,
  clip = "bbox"
)

# Convert the elevation raster to a data frame for ggplot
elev_df <- as.data.frame(
  elev,
  xy = TRUE
)

colnames(elev_df) <- c(
  "x",
  "y",
  "elevation"
)

# =====================================================
# 10. REFERENCE CITIES
# =====================================================

cities_df <- data.frame(
  name = c(
    "Delhi",
    "Chandigarh",
    "Srinagar",
    "Leh",
    "Shimla",
    "Amritsar",
    "Dharamshala",
    "Roorkee"
  ),
  longitude = c(
    77.2090,
    76.7794,
    74.7973,
    77.5619,
    77.1734,
    74.8723,
    76.3213,
    77.8911
  ),
  latitude = c(
    28.6139,
    30.7333,
    34.0837,
    34.1526,
    31.1048,
    31.6340,
    32.2190,
    29.8543
  )
)

india_outline <- india2 %>%
  dplyr::summarise(
    geometry = sf::st_union(geometry)
  )

# Flag DAPC K5 cluster 3 as an outlier in K2 assignments
sites_sf <- sites_sf %>%
  mutate(
    K2_cluster_majoritaire = if_else(
      DAPC_K5_cluster_majoritaire == 3,
      "outlier",
      as.character(K2_cluster_majoritaire)
    )
  )

# =====================================================
# 11. BIOGEOGRAPHICAL MAP: K2
# =====================================================

Plot_Site_K2 <- ggplot() +

  geom_sf(
    data = india_outline,
    fill = "white",
    color = "black",
    linewidth = 1.2
  ) +

  geom_sf(
    data = india2,
    fill = NA,
    color = "black",
    linewidth = 0.5
  ) +

  # Add filled elevation contours
  geom_contour_filled(
    data = elev_df,
    aes(x = x, y = y, z = elevation),
    binwidth = 300,
    alpha = 0.5
  ) +

  scale_fill_grey(
    start = 0.9,
    end = 0.1,
    guide = "none"
  ) +

  # Create a new fill scale for population clusters
  ggnewscale::new_scale_fill() +

  # Add elevation contour lines
  geom_contour(
    data = elev_df,
    aes(x = x, y = y, z = elevation),
    binwidth = 300,
    color = "grey40",
    linewidth = 0.1
  ) +

  # Add reference cities
  geom_point(
    data = cities_df,
    aes(x = longitude, y = latitude),
    shape = 23,
    size = 2.5,
    fill = "black",
    color = "black",
    alpha = 0.8
  ) +

  ggrepel::geom_label_repel(
    data = cities_df,
    aes(
      x = longitude,
      y = latitude,
      label = name
    ),
    size = 3,
    fontface = "italic",
    color = "black",
    fill = "white",
    label.padding = unit(0.2, "lines"),
    nudge_x = 0.1,
    nudge_y = -0.1,
    segment.color = NA
  ) +

  # Add sampling sites coloured by K2 cluster assignment
  geom_point(
    data = sites_sf,
    aes(
      x = longitude,
      y = latitude,
      fill = as.character(K2_cluster_majoritaire)
    ),
    shape = 21,
    color = "black",
    size = 2
  ) +

  scale_fill_manual(
    values = c(
      "1" = "blue",
      "2" = "red",
      "outlier" = "yellow"
    )
  ) +

  annotation_scale(
    location = "bl",
    pad_x = unit(0, "cm")
  ) +

  coord_sf(
    xlim = c(74.5, 81.5),
    ylim = c(28.5, 33.5)
  ) +

  theme_classic(base_size = 8) +
  theme_paper() +

  theme(
    legend.position = "none",
    axis.title = element_blank()
  )

print(Plot_Site_K2)

# =====================================================
# 12. BIOGEOGRAPHICAL MAP: K4
# =====================================================

Plot_Site_K4 <- ggplot() +

  geom_sf(
    data = india_outline,
    fill = "white",
    color = "black",
    linewidth = 1.2
  ) +

  geom_sf(
    data = india2,
    fill = NA,
    color = "black",
    linewidth = 0.5
  ) +

  geom_contour_filled(
    data = elev_df,
    aes(x = x, y = y, z = elevation),
    binwidth = 300,
    alpha = 0.5
  ) +

  scale_fill_grey(
    start = 0.9,
    end = 0.1,
    guide = "none"
  ) +

  ggnewscale::new_scale_fill() +

  geom_contour(
    data = elev_df,
    aes(x = x, y = y, z = elevation),
    binwidth = 300,
    color = "grey40",
    linewidth = 0.1
  ) +

  geom_point(
    data = cities_df,
    aes(x = longitude, y = latitude),
    shape = 23,
    size = 2.5,
    fill = "black",
    color = "black",
    alpha = 0.8
  ) +

  ggrepel::geom_label_repel(
    data = cities_df,
    aes(
      x = longitude,
      y = latitude,
      label = name
    ),
    size = 3,
    fontface = "italic",
    color = "black",
    fill = "white",
    label.padding = unit(0.2, "lines"),
    nudge_x = 0.1,
    nudge_y = -0.1,
    segment.color = NA
  ) +

  geom_point(
    data = sites_sf,
    aes(
      x = longitude,
      y = latitude,
      fill = as.character(K4_cluster_majoritaire)
    ),
    shape = 21,
    color = "black",
    size = 2
  ) +

  scale_fill_manual(
    values = c(
      "1" = "purple",
      "2" = "green",
      "3" = "red",
      "4" = "blue"
    )
  ) +

  annotation_scale(
    location = "bl",
    pad_x = unit(0, "cm")
  ) +

  coord_sf(
    xlim = c(74.5, 81.5),
    ylim = c(28.5, 33.5)
  ) +

  theme_classic(base_size = 8) +
  theme_paper() +

  theme(
    legend.position = "none",
    axis.title = element_blank()
  )

print(Plot_Site_K4)

# =====================================================
# 13. BIOGEOGRAPHICAL MAP: K2 WITH ASSIGNMENT PURITY
# =====================================================

Plot_Site_K2_alpha <- ggplot() +

  geom_sf(
    data = india_outline,
    fill = "white",
    color = "black",
    linewidth = 1.2
  ) +

  geom_sf(
    data = india2,
    fill = NA,
    color = "black",
    linewidth = 0.5
  ) +

  geom_contour_filled(
    data = elev_df,
    aes(x = x, y = y, z = elevation),
    binwidth = 300,
    alpha = 0.5
  ) +

  scale_fill_grey(
    start = 0.9,
    end = 0.1,
    guide = "none"
  ) +

  ggnewscale::new_scale_fill() +

  geom_contour(
    data = elev_df,
    aes(x = x, y = y, z = elevation),
    binwidth = 300,
    color = "grey40",
    linewidth = 0.1
  ) +

  geom_point(
    data = cities_df,
    aes(x = longitude, y = latitude),
    shape = 23,
    size = 2.5,
    fill = "black",
    color = "black",
    alpha = 0.8
  ) +

  ggrepel::geom_label_repel(
    data = cities_df,
    aes(
      x = longitude,
      y = latitude,
      label = name
    ),
    size = 3,
    fontface = "italic",
    color = "black",
    fill = "white",
    label.padding = unit(0.2, "lines"),
    nudge_x = 0.1,
    nudge_y = -0.1,
    segment.color = NA
  ) +

  geom_point(
    data = sites_sf,
    aes(
      x = longitude,
      y = latitude,
      fill = as.character(K2_cluster_majoritaire),
      alpha = K2_statut
    ),
    shape = 21,
    color = "black",
    size = 2
  ) +

  scale_alpha_manual(
    values = c(
      "pure" = 1,
      "admixed" = 0.4
    )
  ) +

  scale_fill_manual(
    values = c(
      "1" = "blue",
      "2" = "red",
      "outlier" = "yellow"
    )
  ) +

  annotation_scale(
    location = "bl",
    pad_x = unit(0, "cm")
  ) +

  coord_sf(
    xlim = c(74.5, 81.5),
    ylim = c(28.5, 33.5)
  ) +

  theme_classic(base_size = 8) +
  theme_paper() +

  theme(
    legend.position = "none",
    axis.title = element_blank()
  )

print(Plot_Site_K2_alpha)

ggsave(
  file.path(figure_dir, "S2_map_K2_alpha_outliers.pdf"),
  Plot_Site_K2_alpha,
  width = 17,
  height = 13,
  device = cairo_pdf,
  units = "cm",
  dpi = 600
)

# =====================================================
# 14. BIOGEOGRAPHICAL MAP: K4 WITH ASSIGNMENT PURITY
# =====================================================

Plot_Site_K4_alpha <- ggplot() +

  geom_sf(
    data = india_outline,
    fill = "white",
    color = "black",
    linewidth = 1.2
  ) +

  geom_sf(
    data = india2,
    fill = NA,
    color = "black",
    linewidth = 0.5
  ) +

  geom_contour_filled(
    data = elev_df,
    aes(x = x, y = y, z = elevation),
    binwidth = 300,
    alpha = 0.5
  ) +

  scale_fill_grey(
    start = 0.9,
    end = 0.1,
    guide = "none"
  ) +

  ggnewscale::new_scale_fill() +

  geom_contour(
    data = elev_df,
    aes(x = x, y = y, z = elevation),
    binwidth = 300,
    color = "grey40",
    linewidth = 0.1
  ) +

  geom_point(
    data = cities_df,
    aes(x = longitude, y = latitude),
    shape = 23,
    size = 2.5,
    fill = "black",
    color = "black",
    alpha = 0.8
  ) +

  ggrepel::geom_label_repel(
    data = cities_df,
    aes(
      x = longitude,
      y = latitude,
      label = name
    ),
    size = 3,
    fontface = "italic",
    color = "black",
    fill = "white",
    label.padding = unit(0.2, "lines"),
    nudge_x = 0.1,
    nudge_y = -0.1,
    segment.color = NA
  ) +

  geom_point(
    data = sites_sf,
    aes(
      x = longitude,
      y = latitude,
      fill = as.character(K4_cluster_majoritaire),
      alpha = K4_statut
    ),
    shape = 21,
    color = "black",
    size = 2
  ) +

  scale_alpha_manual(
    values = c(
      "pure" = 1,
      "admixed" = 0.4
    )
  ) +

  scale_fill_manual(
    values = c(
      "1" = "purple",
      "2" = "green",
      "3" = "red",
      "4" = "blue"
    )
  ) +

  annotation_scale(
    location = "bl",
    pad_x = unit(0, "cm")
  ) +

  coord_sf(
    xlim = c(74.5, 81.5),
    ylim = c(28.5, 33.5)
  ) +

  theme_classic(base_size = 8) +
  theme_paper() +

  theme(
    legend.position = "none",
    axis.title = element_blank()
  )

print(Plot_Site_K4_alpha)

ggsave(
  file.path(figure_dir, "figure1_map_K4_alpha.pdf"),
  Plot_Site_K4_alpha,
  width = 17,
  height = 13,
  device = cairo_pdf,
  units = "cm",
  dpi = 600
)

# =====================================================
# 15. INDIA OVERVIEW MAP
# =====================================================

bbox <- st_bbox(
  c(
    xmin = 60,
    xmax = 95,
    ymin = 5,
    ymax = 40
  ),
  crs = st_crs(4326)
)

bbox_sf <- st_as_sfc(bbox)

# Download a global digital elevation model
# WorldClim resolution 10 corresponds to approximately 10 arc-minutes
elev <- geodata::worldclim_global(
  var = "elev",
  res = 10,
  path = tempdir()
)

# Crop the elevation raster to the study region
elev <- crop(
  elev,
  ext(60, 95, 5, 40)
)

# Convert the elevation raster to a data frame for ggplot
elev_df <- as.data.frame(
  elev,
  xy = TRUE,
  na.rm = TRUE
)

colnames(elev_df) <- c(
  "x",
  "y",
  "elevation"
)

# Download Natural Earth land boundaries
land <- ne_download(
  scale = 10,
  category = "physical",
  type = "land",
  returnclass = "sf"
)

# Crop the land polygons to the study region
land <- st_crop(
  land,
  bbox
)

# Create the overview map
Plot_India <- ggplot() +

  # Add filled elevation contours
  geom_contour_filled(
    data = elev_df,
    aes(x = x, y = y, z = elevation),
    binwidth = 1000,
    alpha = 0.45
  ) +

  scale_fill_grey(
    start = 0.95,
    end = 0.10,
    guide = "none"
  ) +

  ggnewscale::new_scale_fill() +

  # Add elevation contour lines
  geom_contour(
    data = elev_df,
    aes(x = x, y = y, z = elevation),
    binwidth = 1000,
    colour = "grey40",
    linewidth = 0.10
  ) +

  # Add land boundaries
  geom_sf(
    data = land,
    fill = NA,
    colour = "black",
    linewidth = 0.2
  ) +

  theme_classic(base_size = 8) +
  theme_paper() +

  theme(
    legend.position = "none",
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank()
  )

print(Plot_India)

ggsave(
  file.path(figure_dir, "figure1_map_India.pdf"),
  Plot_India,
  width = 3,
  height = 3,
  device = cairo_pdf,
  units = "cm",
  dpi = 600
)

# =====================================================
# 16. OPTIONAL ENVIRONMENTAL DATA EXTRACTION
# =====================================================

# Extract environmental variables such as elevation
# get_elev_raster() can be used for elevation data.

# elev.grid <- get_elev_raster(sites_sf, z = 6)
# elev <- terra::extract(elev.grid[[1]], sites_sf)
# getNASADEM(sites_sf, ID = "Site")

# Extract ISRIC soil data
# Soil data are available at approximately 250 m spatial resolution.

# bdod_grid <- getISRIC_soils(
#   sites_sf,
#   varname = "bdod_5-15cm_mean",
#   ID = "Site"
# )
#
# ph_grid <- getISRIC_soils(
#   sites_sf,
#   varname = "phh2o_5-15cm_mean",
#   ID = "Site"
# )
#
# sand_grid <- getISRIC_soils(
#   sites_sf,
#   varname = "sand_5-15cm_mean",
#   ID = "Site"
# )
#
# silt_grid <- getISRIC_soils(
#   sites_sf,
#   varname = "silt_5-15cm_mean",
#   ID = "Site"
# )
#
# clay_grid <- getISRIC_soils(
#   sites_sf,
#   varname = "clay_5-15cm_mean",
#   ID = "Site"
# )
#
# bdod <- terra::extract(
#   bdod_grid[[1]], sites_sf,
#   fun = mean, method = "bilinear"
# )
#
# ph <- terra::extract(
#   ph_grid[[1]], sites_sf,
#   fun = mean, method = "bilinear"
# )
#
# sand <- terra::extract(
#   sand_grid[[1]], sites_sf,
#   fun = mean, method = "bilinear"
# )
#
# silt <- terra::extract(
#   silt_grid[[1]], sites_sf,
#   fun = mean, method = "bilinear"
# )
#
# clay <- terra::extract(
#   clay_grid[[1]], sites_sf,
#   fun = mean, method = "bilinear"
# )

# =====================================================
# 17. WHITTAKER BIOME PLOTS
# =====================================================

biogeo <- df2

# The function whittaker_base_plot() must be defined or loaded
# from the relevant script/package before running this section.

# K2 Whittaker plot
Whitt_terra_K2 <- whittaker_base_plot() +

  ggnewscale::new_scale_fill() +

  geom_point(
    data = biogeo,
    aes(
      x = (tmaxM + tminM) / 2,
      y = pptM,
      fill = as.character(K2_cluster_majoritaire),
      alpha = K2_statut
    ),
    shape = 21,
    color = "black",
    size = 1
  ) +

  scale_alpha_manual(
    values = c(
      "pure" = 1,
      "admixed" = 0.4
    )
  ) +

  scale_fill_manual(
    values = c(
      "1" = "blue",
      "2" = "red"
    )
  ) +

  labs(
    x = "Mean annual temperature (°C)",
    y = "Annual precipitation (cm)"
  ) +

  # Set clean axis expansion
  scale_x_continuous(
    expand = expansion(mult = 0.02)
  ) +

  scale_y_continuous(
    expand = expansion(mult = 0.02)
  ) +

  theme_classic(base_size = 8) +
  theme_paper() +

  theme(
    legend.position = "none",
    panel.grid.major = element_line(
      color = "grey85",
      linewidth = 0.3
    ),
    panel.grid.minor = element_blank()
  )

print(Whitt_terra_K2)

ggsave(
  file.path(figure_dir, "figure3_whittaker_K2.pdf"),
  Whitt_terra_K2,
  width = 8,
  height = 8,
  device = cairo_pdf,
  units = "cm",
  dpi = 600
)

# K4 Whittaker plot
Whitt_terra_K4 <- whittaker_base_plot() +

  ggnewscale::new_scale_fill() +

  geom_point(
    data = biogeo,
    aes(
      x = (tmaxM + tminM) / 2,
      y = pptM,
      fill = as.character(K4_cluster_majoritaire),
      alpha = K4_statut
    ),
    shape = 21,
    color = "black",
    size = 1
  ) +

  scale_alpha_manual(
    values = c(
      "pure" = 1,
      "admixed" = 0.4
    )
  ) +

  scale_fill_manual(
    values = c(
      "1" = "purple",
      "2" = "green",
      "3" = "red",
      "4" = "blue"
    )
  ) +

  labs(
    x = "Mean annual temperature (°C)",
    y = "Annual precipitation (cm)"
  ) +

  scale_x_continuous(
    expand = expansion(mult = 0.02)
  ) +

  scale_y_continuous(
    expand = expansion(mult = 0.02)
  ) +

  theme_classic(base_size = 8) +
  theme_paper() +

  theme(
    legend.position = "none",
    panel.grid.major = element_line(
      color = "grey85",
      linewidth = 0.3
    ),
    panel.grid.minor = element_blank()
  )

print(Whitt_terra_K4)

ggsave(
  file.path(figure_dir, "figure3_whittaker_K4.pdf"),
  Whitt_terra_K4,
  width = 8,
  height = 8,
  device = cairo_pdf,
  units = "cm",
  dpi = 600
)

# Export the biome reference legend
Whitt_legend <- whittaker_base_plot()

print(Whitt_legend)

ggsave(
  file.path(figure_dir, "figure3_whittaker_legendbiomes.pdf"),
  Whitt_legend,
  width = 8,
  height = 8,
  device = cairo_pdf,
  units = "cm",
  dpi = 600
)

# =====================================================
# 18. PREPARE CLUSTER MEMBERSHIP MATRICES
# =====================================================

qmat_k2clean <- as.data.frame(
  popk4[, c(2:3)]
)

qmat_k4clean <- as.data.frame(
  popk4[, c(4:7)]
)

rownames(qmat_k2clean) <- popk4$sample_id
rownames(qmat_k4clean) <- popk4$sample_id

# =====================================================
# 19. LOAD CLUSTER MEMBERSHIP LISTS
# =====================================================

# Load individual membership lists for K2 and K4

# K2 majority-assignment groups
k2_c1_raw <- readLines(
  file.path(cluster_dir, "3d_K2_cluster1_majo.txt")
)

k2_c2_raw <- readLines(
  file.path(cluster_dir, "3d_K2_cluster2_majo.txt")
)

# K2 pure-assignment groups
k2_c1_pure <- readLines(
  file.path(cluster_dir, "3d_K2_cluster1_pur.txt")
)

k2_c2_pure <- readLines(
  file.path(cluster_dir, "3d_K2_cluster2_pur.txt")
)

# K4 majority-assignment groups
k4_c1_raw <- readLines(
  file.path(cluster_dir, "3d_K4_cluster1_majo.txt")
)

k4_c2_raw <- readLines(
  file.path(cluster_dir, "3d_K4_cluster2_majo.txt")
)

k4_c3_raw <- readLines(
  file.path(cluster_dir, "3d_K4_cluster3_majo.txt")
)

k4_c4_raw <- readLines(
  file.path(cluster_dir, "3d_K4_cluster4_majo.txt")
)

# K4 pure-assignment groups
k4_c1_pure <- readLines(
  file.path(cluster_dir, "3d_K4_cluster1_pur.txt")
)

k4_c2_pure <- readLines(
  file.path(cluster_dir, "3d_K4_cluster2_pur.txt")
)

k4_c3_pure <- readLines(
  file.path(cluster_dir, "3d_K4_cluster3_pur.txt")
)

k4_c4_pure <- readLines(
  file.path(cluster_dir, "3d_K4_cluster4_pur.txt")
)

# =====================================================
# 20. GENOTYPE MATRIX AND GENIND OBJECTS
# =====================================================

# Calculate Nei's Da genetic distances and order ancestry
# barplots according to the resulting tree.

geno_mat <- snpgdsGetGeno(
  genofile,
  sample.id = sample
)

rownames(geno_mat) <- sample_new
colnames(geno_mat) <- snps.ids

# Convert genotype matrices into genind objects
make_genind_from_clusters <- function(
    geno_mat,
    cluster_lists,
    cluster_labels
) {

  # Combine all individuals
  inds <- unlist(cluster_lists)

  # Keep unique individuals while preserving their order
  inds <- unique(inds)

  # Extract the genotype submatrix
  sub_mat <- geno_mat[
    inds,
    ,
    drop = FALSE
  ]

  # Create the population assignment vector
  pop_vector <- rep(
    NA,
    length(inds)
  )

  names(pop_vector) <- inds

  for (i in seq_along(cluster_lists)) {

    pop_vector[cluster_lists[[i]]] <- cluster_labels[i]
  }

  # Create the genind object
  gi <- df2genind(
    sub_mat,
    ploidy = 2,
    type = "codom",
    ind.names = inds,
    pop = factor(pop_vector),
    ncode = 1
  )

  return(gi)
}

# K2 majority-assignment groups
gi_k2_raw <- make_genind_from_clusters(
  geno_mat = geno_mat,
  cluster_lists = list(
    k2_c1_raw,
    k2_c2_raw
  ),
  cluster_labels = c(
    "K2_C1",
    "K2_C2"
  )
)

# K2 pure-assignment groups
gi_k2_pure <- make_genind_from_clusters(
  geno_mat = geno_mat,
  cluster_lists = list(
    k2_c1_pure,
    k2_c2_pure
  ),
  cluster_labels = c(
    "K2_C1",
    "K2_C2"
  )
)

# K4 majority-assignment groups
gi_k4_raw <- make_genind_from_clusters(
  geno_mat = geno_mat,
  cluster_lists = list(
    k4_c1_raw,
    k4_c2_raw,
    k4_c3_raw,
    k4_c4_raw
  ),
  cluster_labels = c(
    "K4_C1",
    "K4_C2",
    "K4_C3",
    "K4_C4"
  )
)

# K4 pure-assignment groups
gi_k4_pure <- make_genind_from_clusters(
  geno_mat = geno_mat,
  cluster_lists = list(
    k4_c1_pure,
    k4_c2_pure,
    k4_c3_pure,
    k4_c4_pure
  ),
  cluster_labels = c(
    "K4_C1",
    "K4_C2",
    "K4_C3",
    "K4_C4"
  )
)

# =====================================================
# 21. CONVERT TO HIERFSTAT AND CALCULATE NEI'S DA
# =====================================================

hf_K2_raw <- genind2hierfstat(gi_k2_raw)
hf_K2_pure <- genind2hierfstat(gi_k2_pure)
hf_K4_raw <- genind2hierfstat(gi_k4_raw)
hf_K4_pure <- genind2hierfstat(gi_k4_pure)

da_K2_raw <- genet.dist(
  hf_K2_raw,
  method = "Da"
)

da_K2_pure <- genet.dist(
  hf_K2_pure,
  method = "Da"
)

da_K4_raw <- genet.dist(
  hf_K4_raw,
  method = "Da"
)

da_K4_pure <- genet.dist(
  hf_K4_pure,
  method = "Da"
)

# =====================================================
# 22. NEI'S DA NEIGHBOUR-JOINING TREES WITH BOOTSTRAP
# =====================================================

bootstrap_da_nj <- function(
    hf,
    nboot = 100
) {

  # Build the observed tree
  da_obs <- genet.dist(
    hf,
    method = "Da"
  )

  tree_obs <- nj(
    as.dist(da_obs)
  )

  # Bootstrap loci
  loci_cols <- 2:ncol(hf)
  trees_boot <- list()

  pb <- txtProgressBar(
    min = 0,
    max = nboot,
    style = 3
  )

  for (b in 1:nboot) {

    sampled_loci <- sample(
      loci_cols,
      size = length(loci_cols),
      replace = TRUE
    )

    hf_boot <- hf[
      ,
      c(1, sampled_loci)
    ]

    res <- try({

      da_boot <- genet.dist(
        hf_boot,
        method = "Da"
      )

      nj(
        as.dist(da_boot)
      )

    }, silent = TRUE)

    if (!inherits(res, "try-error")) {

      trees_boot[[length(trees_boot) + 1]] <- res
    }

    setTxtProgressBar(
      pb,
      b
    )
  }

  close(pb)

  # Calculate bootstrap support
  bs <- prop.clades(
    tree_obs,
    trees_boot
  )

  return(
    list(
      tree = tree_obs,
      bootstrap = bs,
      boot_trees = trees_boot
    )
  )
}

# Run bootstrap analyses for K4
res_K4_raw <- bootstrap_da_nj(
  hf_K4_raw
)

res_K4_pure <- bootstrap_da_nj(
  hf_K4_pure
)

# =====================================================
# 23. PLOT NJ TREES WITH BOOTSTRAP SUPPORT
# =====================================================

plot_nj_bootstrap_pdf <- function(
    res,
    file = NULL,
    cex = 0.9,
    show_edge_length = TRUE,
    show_bootstrap = TRUE
) {

  tree <- res$tree
  bs <- res$bootstrap

  # Export as a vector PDF when a file is specified
  if (!is.null(file)) {

    pdf(
      file,
      width = 8,
      height = 6
    )

    on.exit(
      dev.off(),
      add = TRUE
    )
  }

  # Plot the base tree
  plot(
    tree,
    type = "phylogram",
    cex = cex,
    edge.width = 2,
    font = 1,
    main = "NJ tree (Nei Da) + bootstrap"
  )

  # Display branch lengths
  if (show_edge_length) {

    edgelabels(
      round(tree$edge.length, 4),
      cex = 0.7,
      frame = "none",
      col = "darkblue"
    )
  }

  # Display bootstrap values at internal nodes
  if (show_bootstrap) {

    nodelabels(
      bs,
      frame = "none",
      cex = 0.8,
      col = "black"
    )
  }

  return(
    invisible(tree)
  )
}

# Export K4 majority-assignment tree
pdf(
  file.path(
    figure_dir,
    "figure1_arbre_K4_raw_bootstraps_3d.pdf"
  ),
  width = 4,
  height = 4
)

plot_nj_bootstrap_pdf(
  res_K4_raw
)

dev.off()

# Export K4 pure-assignment tree
pdf(
  file.path(
    figure_dir,
    "figure1_arbre_K4_pure_bootstraps_3d.pdf"
  ),
  width = 4,
  height = 4
)

plot_nj_bootstrap_pdf(
  res_K4_pure
)

dev.off()

# =====================================================
# 24. STRUCTURE ANCESTRY BARPLOTS: K4
# =====================================================

# Order clusters from largest to smallest according to the
# selected cluster order, then apply the same individual
# ordering to K2.

library(dplyr)
library(tidyr)
library(ggplot2)

# Define colours
fill_values <- c(
  "1_pure" = "purple",
  "2_pure" = "green",
  "3_pure" = "red",
  "4_pure" = "blue",

  "1_admixed" = "mediumpurple3",
  "2_admixed" = "palegreen3",
  "3_admixed" = "tomato3",
  "4_admixed" = "steelblue3"
)

# Define the order of major-assignment clusters
cluster_order <- c(
  "4",
  "2",
  "1",
  "3"
