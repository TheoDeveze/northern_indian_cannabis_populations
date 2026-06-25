############################################
## ENVIRONMENT vs SPACE ON STRUCTURE Q (CLR)
############################################

# =========================
# Libraries
# =========================
library(corrplot)
library(vegan)
library(compositions)
library(adespatial)
library(eulerr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(readr)

# =========================
# 1. Load data
# =========================
# Remove missing samples and duplicates
data <- read_tsv("path/to/Q_biogeography_matrix.tsv") %>%
  filter(!numero %in% c(78, 425)) %>%
  distinct(geometry, .keep_all = TRUE)

# =========================
# 2. Define variable groups
# =========================
Q_cols       <- c("K4_cluster1","K4_cluster2","K4_cluster3","K4_cluster4")
bioclim_cols <- c("elev","tmaxM","tminM","pptM","ph","bdod")
soil_cols    <- c("clay","sand","silt")
coord_cols   <- c("longitude","latitude")

# Remove incomplete rows
data2 <- data %>%
  drop_na(all_of(c(Q_cols, bioclim_cols, soil_cols, coord_cols)))

# =========================
# 3. Environmental data scaling
# =========================
env <- cbind(data2[, bioclim_cols], data2[, soil_cols])
env_scaled <- scale(env)

# Correlation diagnostics (environment)
corrplot(cor(env_scaled), method = "color", type = "upper", tl.cex = 0.8)

# =========================
# 4. STRUCTURE Q → CLR transform
# =========================
Q <- as.matrix(data2[, Q_cols])
Q <- Q / rowSums(Q)
Q[Q <= 0] <- 1e-6
Q_ilr <- ilr(Q)

# =========================
# 5. Soil composition → ILR
# =========================
soil <- as.matrix(data2[, soil_cols])
soil <- soil / rowSums(soil)
soil[soil <= 0] <- 1e-6
soil_ilr <- ilr(soil)

# =========================
# 6. Temperature / elevation PCA (avoid collinearity)
# =========================
temp_elev <- env_scaled[, c("elev","tminM","tmaxM")]
pca_temp <- prcomp(temp_elev)

PC_temp <- pca_temp$x[,1, drop = FALSE]
colnames(PC_temp) <- "temp_elev_gradient"

# =========================
# 7. Final environmental matrix
# =========================
env_final <- cbind(
  data2[, c("pptM","bdod","ph")],
  soil_ilr,
  PC_temp
)

env_final <- scale(env_final)

# =========================
# 8. Environmental correlations
# =========================
corrplot(cor(env_final), method = "color", type = "upper", tl.cex = 0.8)

# =========================
# 9. RDA: environment effect
# =========================
rda_env <- rda(Q_ilr ~ ., data = as.data.frame(env_final))
anova(rda_env, permutations = 999)

# =========================
# 10. Spatial variables (dbMEM)
# =========================
coords <- as.matrix(data2[, coord_cols])
mem <- dbmem(coords)
mem_df <- as.data.frame(mem)

# =========================
# 11. Forward selection of spatial variables
# =========================
rda_space_full <- rda(Q_ilr ~ ., data = mem_df)
R2_space <- RsquareAdj(rda_space_full)$adj.r.squared

fs_mem <- forward.sel(
  Y = Q_ilr,
  X = mem_df,
  adjR2thresh = R2_space,
  nperm = 999
)

mem_sel <- mem_df[, fs_mem$variables, drop = FALSE]

# =========================
# 12. Spatial-only RDA
# =========================
rda_space <- rda(Q_ilr ~ ., data = mem_sel)
anova(rda_space, permutations = 999)

# =========================
# 13. Variation partitioning
# =========================
vp <- varpart(Q_ilr, env_final, mem_sel)
plot(vp)

# =========================
# 14. Euler diagram (shared effects)
# =========================
rda_both <- rda(Q_ilr ~ ., data = cbind(env_final, mem_sel))

R2_env   <- RsquareAdj(rda_env)$adj.r.squared
R2_space <- RsquareAdj(rda_space)$adj.r.squared
R2_both  <- RsquareAdj(rda_both)$adj.r.squared

shared <- (R2_env + R2_space) - R2_both

euler_vals <- c(
  Environment = R2_env - shared,
  Space = R2_space - shared,
  "Shared" = shared
)

plot(euler(euler_vals), quantities = TRUE)

# =========================
# 15. RDA plot (environment)
# =========================

sites <- as.data.frame(scores(rda_env, display = "sites", scaling = 2))
env_vectors <- as.data.frame(scores(rda_env, display = "bp", scaling = 2))

sites$cluster <- colnames(Q)[max.col(Q)]

eig <- summary(rda_env)$cont$importance[2,1:2] * 100

rda_plot <- ggplot() +
  geom_point(data = sites,
             aes(RDA1, RDA2, color = cluster),
             size = 0.8, alpha = 0.7) +
  stat_ellipse(data = sites,
               aes(RDA1, RDA2, group = cluster, color = cluster),
               type = "norm",
               linetype = "dashed",
               linewidth = 0.3,
               level = 0.68) +
  geom_segment(data = env_vectors,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "grey50",
               linewidth = 0.3) +
  geom_text(data = env_vectors,
            aes(RDA1, RDA2, label = rownames(env_vectors)),
            size = 2.5) +
  theme_bw(base_size = 8) +
  labs(
    x = paste0("RDA1 (", round(eig[1],1), "%)"),
    y = paste0("RDA2 (", round(eig[2],1), "%)")
  )

rda_plot

# =========================
# 16. K2: variation partitioning + regressions
# =========================

Q2 <- data2$K2_cluster1

mod_env <- lm(Q2 ~ ., data = as.data.frame(env_final))
R2_env <- summary(mod_env)$adj.r.squared

mod_space <- lm(Q2 ~ ., data = mem_sel)
R2_space <- summary(mod_space)$adj.r.squared

mod_both <- lm(Q2 ~ ., data = cbind(env_final, mem_sel))
R2_both <- summary(mod_both)$adj.r.squared

shared <- (R2_env + R2_space) - R2_both

euler_vals_k2 <- c(
  Environment = R2_env - shared,
  Space = R2_space - shared,
  Shared = shared
)

plot(euler(euler_vals_k2), quantities = TRUE)

# =========================
# 17. Environment vs K2 regression plots
# =========================

df_long <- cbind(Q2, env_final) %>%
  pivot_longer(-Q2, names_to = "variable", values_to = "value")

stats <- df_long %>%
  group_by(variable) %>%
  summarise(
    R2 = summary(lm(Q2 ~ value))$r.squared,
    p  = summary(lm(Q2 ~ value))$coefficients[2,4]
  ) %>%
  mutate(label = paste0("R²=", round(R2,2), "\n p=", signif(p,2)))

ggplot(df_long, aes(Q2, value)) +
  geom_point(alpha = 0.5, size = 0.3) +
  geom_smooth(method = "lm", linewidth = 0.3) +
  facet_wrap(~variable, scales = "free_y") +
  theme_bw(base_size = 8)
