############################################
## ENVIRONMENT VS SPACE ON STRUCTURE Q (ILR)
############################################

# --- Libraries ---
library(corrplot)
library(vegan)
library(compositions)
library(adespatial)
library(eulerr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(knitr)
library(readr)

# --- Project-relative paths ---
# Run this script from the project root directory.
# Expected input: data/metadata/final3d_Q_biogeo_matrix.tsv
# Outputs: results/figures/population_structure/ and results/tables/
figure_dir <- file.path("results", "figures", "population_structure")
table_dir  <- file.path("results", "tables")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# --- 1) Load and clean data ---
data <- read_tsv(
  file.path("data", "metadata", "final3d_Q_biogeo_matrix.tsv")
)

# Filter samples with missing values and remove duplicate geometries
data <- data %>%
  filter(!numero %in% c(78, 425)) %>%
  distinct(geometry, .keep_all = TRUE)

# --- 2) Select variables ---
Q_cols       <- c("K4_cluster1", "K4_cluster2", "K4_cluster3", "K4_cluster4")
bioclim_cols <- c("elev", "tmaxM", "tminM", "pptM", "ph", "bdod")
soil_cols    <- c("clay", "sand", "silt")
coord_cols   <- c("longitude", "latitude")

data2 <- data %>%
  drop_na(all_of(c(Q_cols, bioclim_cols, soil_cols, coord_cols)))

# Initial environmental matrix
env <- cbind(
  data2[, bioclim_cols],
  data2[, soil_cols]
)

env_scaled <- as.data.frame(scale(env))  # Scaling is required for RDA

# --- 3) Correlation diagnostics among environmental variables ---
cor_env <- cor(env_scaled, method = "pearson")

corrplot(
  cor_env,
  method = "color",
  type = "upper",
  tl.cex = 0.8,
  addCoef.col = "black"
)

# --- 4) STRUCTURE Q: ILR transformation ---
Q <- as.matrix(data2[, Q_cols])
Q <- Q / rowSums(Q)
Q[Q <= 0] <- 1e-6
Q_ilr <- ilr(Q)

# --- 5) Soil variables: compositional (ILR) transformation ---
soil <- as.matrix(data2[, soil_cols])
soil <- soil / rowSums(soil)
soil[soil <= 0] <- 1e-6
soil_ilr <- ilr(soil)

# --- 5.1) Highly correlated variables: elevation and minimum/maximum temperature (PCA) ---
temp_elev <- env_scaled[, c("elev", "tminM", "tmaxM")]
pca_temp <- prcomp(temp_elev)
var_expl <- (pca_temp$sdev^2) / sum(pca_temp$sdev^2) * 100
var_expl

biplot(
  pca_temp,
  xlab = paste0("PC1 (", round(var_expl[1], 1), " %)"),
  ylab = paste0("PC2 (", round(var_expl[2], 1), " %)")
)

PC_temp <- pca_temp$x[, 1, drop = FALSE]
colnames(PC_temp) <- "gradient elev temp"

# Eigenvalues
eig <- pca_temp$sdev^2

# Percentage of variance explained
var_explained <- eig / sum(eig) * 100

# PC1 only
var_explained[1]

# Loadings on PC1
loadings_PC1 <- pca_temp$rotation[, 1]

# Display loadings
round(loadings_PC1, 3)

# --- 6) Environmental matrix (climate + soil) ---
# Select precipitation, soil bulk density, and pH, then add soil ILR coordinates and the elevation-temperature gradient.
# This selection assumes pptM, ph, and bdod are columns 8:10 in data2, as in the original script.
env <- cbind(
  data2[, c(8:10)],
  as.data.frame(soil_ilr),
  as.data.frame(PC_temp)
)
env <- env %>% rename(soil1 = V1, soil2 = V2)

env_scaled <- as.data.frame(scale(env))  # Scaling is required for RDA

# Optional: export the scaled environmental matrix
# write_tsv(
#   env_scaled,
#   file.path(table_dir, "environmental_centered_scaled_transformed.tsv")
# )

# --- 7) Correlation diagnostics among final environmental variables ---
cor_env <- cor(env_scaled, method = "pearson")

corrplot(
  cor_env,
  method = "color",
  type = "upper",
  tl.cex = 0.8,
  addCoef.col = "black"
)

soil_df <- cbind(data2[, Q_cols], data2[, soil_cols])
cor_soil <- cor(soil_df, method = "pearson")

corrplot(
  cor_soil,
  method = "color",
  type = "upper",
  tl.cex = 0.8,
  addCoef.col = "black"
)

# --- 8) RDA: environment only ---
rda_env_model <- rda(Q_ilr ~ ., data = env_scaled)
anova(rda_env_model, permutations = 999)

# --- 9) Spatial variables (dbMEMs) ---
coords <- as.matrix(data2[, coord_cols])
mem <- dbmem(coords)
mem_df <- as.data.frame(mem)

# --- 10) Forward selection of MEMs ---
rda_space_full <- rda(Q_ilr ~ ., data = mem_df)
R2adj_space <- RsquareAdj(rda_space_full)$adj.r.squared

fs_mem <- forward.sel(
  Y = Q_ilr,
  X = mem_df,
  adjR2thresh = R2adj_space,
  nperm = 999
)

mem_sel <- mem_df[, fs_mem$variables, drop = FALSE]

# --- 11) RDA: spatial component only ---
rda_space <- rda(Q_ilr ~ ., data = mem_sel)
anova(rda_space, permutations = 999)

# --- 12) Variation partitioning (environment vs space) ---
vp <- varpart(Q_ilr, env_scaled, mem_sel)
plot(vp, Xnames = c("Environment", "Space"))

# --- 13) Proportional Venn / Euler diagram ---
rda_both <- rda(Q_ilr ~ ., data = cbind(env_scaled, mem_sel))

R2_env   <- RsquareAdj(rda_env_model)$adj.r.squared
R2_space <- RsquareAdj(rda_space)$adj.r.squared
R2_both  <- RsquareAdj(rda_both)$adj.r.squared

shared <- (R2_env + R2_space) - R2_both
resid  <- 1 - R2_both

vals <- c(
  Environment = round(R2_env - shared, 2),
  Space = round(R2_space - shared, 2),
  "Environment&Space" = round(shared, 2)
  # Residual = round(1 - R2_both, 2)
)

fit <- euler(vals)

plot(
  fit,
  quantities = TRUE
)

# --- 14) RDA plot: environment ---
rda_env_plot_model <- rda(Q_ilr ~ ., data = env_scaled)

# Percentage of constrained variance explained by the first two axes
eig_env <- summary(rda_env_plot_model)$cont$importance[2, 1:2] * 100

# Site scores
sites_env <- as.data.frame(scores(rda_env_plot_model, display = "sites", scaling = 2))
sites_env$ID <- rownames(sites_env)

# Environmental variable scores
env_vectors <- as.data.frame(scores(rda_env_plot_model, display = "bp", scaling = 2))
env_vectors$var <- rownames(env_vectors)

# Assign each sample to its major STRUCTURE cluster using the untransformed Q matrix
cluster_major <- colnames(Q)[max.col(Q, ties.method = "first")]
alpha <- ifelse(apply(Q, 1, max) > 0.7, 0.9, 0.3)

sites_env$cluster <- cluster_major
sites_env$alpha <- alpha

env_vectors$var[env_vectors$var == "`gradient elev temp`"] <- "gradient elev temp"

rda_env_plot <- ggplot() +
  geom_point(
    data = sites_env,
    aes(x = RDA1, y = RDA2, color = cluster),
    alpha = sites_env$alpha,
    size = 0.8
  ) +
  scale_color_manual(values = c(
    "K4_cluster1" = "purple",
    "K4_cluster2" = "green",
    "K4_cluster3" = "red",
    "K4_cluster4" = "blue"
  )) +
  stat_ellipse(
    data = sites_env,
    aes(
      x = RDA1,
      y = RDA2,
      group = as.factor(cluster),
      color = cluster
    ),
    type = "norm",
    linetype = "dashed",
    linewidth = 0.3,
    level = 0.68
  ) +
  geom_segment(
    data = env_vectors,
    aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
    arrow = arrow(length = unit(0.25, "cm")),
    linewidth = 0.3,
    color = "grey55"
  ) +
  geom_text(
    data = env_vectors,
    aes(x = RDA1, y = RDA2, label = var),
    size = 2.5,
    vjust = -0.4
  ) +
  labs(
    x = paste0("RDA1 (", round(eig_env[1], 1), "%)"),
    y = paste0("RDA2 (", round(eig_env[2], 1), "%)")
  ) +
  theme_bw(base_size = 7) +
  theme(
    aspect.ratio = 1,
    legend.position = "right",
    text = element_text(family = "sans", size = 8),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    strip.text = element_text(size = 8),
    plot.title = element_text(size = 9, face = "bold")
  ) +
  guides(color = "none", alpha = "none")

rda_env_plot

ggsave(
  file.path(figure_dir, "figure3_rda_K4.pdf"),
  rda_env_plot,
  width = 8,
  height = 8,
  device = cairo_pdf,
  units = "cm",
  dpi = 600
)

# TODO: export vector figures, generate the K2 analysis, and interpret soil1 and soil2.

# --- 15) K2 analysis ---
Q2 <- data2$K2_cluster1

mod_env <- lm(
  Q2 ~ .,
  data = env_scaled
)

R2_env <- summary(mod_env)$adj.r.squared

mod_space_full <- lm(
  Q2 ~ .,
  data = mem_df
)

R2adj_space <- summary(mod_space_full)$adj.r.squared

fs_mem <- forward.sel(
  Y = Q2,
  X = mem_df,
  adjR2thresh = R2adj_space,
  nperm = 999
)

mem_sel <- mem_df[, fs_mem$variables, drop = FALSE]

mod_space <- lm(
  Q2 ~ .,
  data = mem_sel
)

R2_space <- summary(mod_space)$adj.r.squared

mod_both <- lm(
  Q2 ~ .,
  data = cbind(env_scaled, mem_sel)
)

R2_both <- summary(mod_both)$adj.r.squared

shared <- (R2_env + R2_space) - R2_both

vals <- c(
  Environment = round(R2_env - shared, 2),
  Space = round(R2_space - shared, 2),
  "Environment&Space" = round(shared, 2)
)

fit <- euler(vals)

plot(
  fit,
  quantities = TRUE
)

# --- 16) Relationship between K2 ancestry and environmental variables ---
df_plot <- cbind(Q2 = Q2, env_scaled)

df_long <- df_plot %>%
  pivot_longer(
    cols = colnames(env_scaled),
    names_to = "variable",
    values_to = "value"
  )

stats <- df_long %>%
  group_by(variable) %>%
  do({
    mod <- lm(Q2 ~ value, data = .)

    data.frame(
      R2 = summary(mod)$r.squared,
      p = summary(mod)$coefficients[2, 4]
    )
  }) %>%
  ungroup()

# Clean annotation labels
stats <- stats %>%
  mutate(
    label = paste0(
      "R² = ", round(R2, 2),
      "\n",
      "p = ", signif(p, 2)
    )
  )

# Position annotations within each facet
label_pos <- df_long %>%
  group_by(variable) %>%
  summarise(
    x = 0.98,
    y = min(value, na.rm = TRUE),
    .groups = "drop"
  )

stats <- left_join(
  stats,
  label_pos,
  by = "variable"
)

# Plot
p <- ggplot(
  df_long,
  aes(
    x = Q2,
    y = value
  )
) +
  geom_point(
    alpha = 0.7,
    size = 0.05
  ) +
  geom_smooth(
    method = "lm",
    linewidth = 0.2
  ) +
  geom_text(
    data = stats,
    aes(
      x = x,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 1,
    vjust = 0,
    size = 1
  ) +
  facet_wrap(
    ~ variable,
    scales = "free_y"
  ) +
  xlim(0, 1) +
  labs(
    y = NULL,
    x = "K2 cluster 1 ancestry"
  ) +
  scale_x_continuous(
    labels = function(x) format(round(x, 2), nsmall = 2)
  ) +
  theme_bw(base_size = 7) +
  theme(
    text = element_text(family = "sans", size = 8),
    strip.text = element_text(size = 6),
    axis.text = element_text(size = 4),
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 7),
    plot.title = element_text(size = 5, face = "bold")
  )

p

ggsave(
  file.path(figure_dir, "figure3_linear_reg_K2.pdf"),
  p,
  width = 8,
  height = 5,
  device = cairo_pdf,
  units = "cm",
  dpi = 600
)
