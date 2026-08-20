############################################################
# run_all.R — Full NAPESHM SHMI pipeline
#
# This script executes the complete NAPESHM workflow:
#   • Load + clean raw NAPESHM data (soil, PLFA, treatments, sites)
#   • Fit Soil Health latent variable SEM (SH, SH.Bio, SH.CN, SH.Phys)
#   • (Optional) Bayesian optimization of SHMI weights
#   • Build final SHMI scores using optimized or expert weights
#   • Merge SHMI with soil health + yield data
#   • Generate manuscript figures (Fig. 2a–c)
#
# All functions are modularized in /NAPESHM/R/
############################################################

source("NAPESHM/R/setup.R")

# ---------------------------------------------------------------------------
# Step 1. Load and prepare data
#
# load_napes_data():
#   • Reads raw NAPESHM CSVs (sites, treatments, soil, PLFA)
#   • Constructs MGT_combo keys
#   • Removes outliers within climate zones
#   • Adjusts indicators for climate + texture (ACE, POXC, SOC, etc.)
#   • Returns finalDF containing all soil health indicators + metadata
# ---------------------------------------------------------------------------

finalDF <- load_napes_data()

# ---------------------------------------------------------------------------
# Step 2. Run Soil Health SEM
#
# run_lvSH_model():
#   • Fits latent variable SEM with three subcomponents:
#         SH.Bio   (biological)
#         SH.CN    (chemical/nutrient)
#         SH.Phys  (physical)
#   • Produces lv_SH (overall soil health)
#   • Appends latent scores to finalDF
# ---------------------------------------------------------------------------

lv_out <- run_lvSH_model(finalDF)
finalDF <- lv_out$finalDF

# ---------------------------------------------------------------------------
# Step 3. Bayesian optimization (optional)
#
# run_mbo():
#   • Uses mlrMBO to tune SHMI weights:
#         seasonal cover weights
#         diversity Hill number
#         disturbance transform
#         subindex weights (cover, diversity, invdist, animals)
#   • Objective: minimize MAE between SHMI and latent SH
#   • Produces mbo_result.rds + summary table
#
# extract_mbo_results():
#   • Returns top 100 models + consensus parameter set
# ---------------------------------------------------------------------------

if (FALSE) {
  sh_scores <- finalDF %>%
    dplyr::select(MGT_combo, land_use, SH, SH.Bio, SH.CN, SH.Phys)
  
  run_mbo(sh_scores)
  mbo_summary <- extract_mbo_results()
}

# current optimized SHMI settings used in manuscript
mySettings <- list(
  # cover
  w_winter = 0.130,
  w_spring = 0.129,
  w_summer = 0.513,
  w_fall   = 0.227,
  
  # diversity
  hill      = 2,
  max_div   = 8,
  
  # disturbance
  transform   = "raw",
  max_stir    = 205,
  
  # animals
  w_amend   = 1,
  w_animals = 1,
  
  # shmi weights
  w_cover    = 0.492,
  w_div      = 0.052,
  w_dist     = 0.131,
  w_ani      = 0.324
)

# ---------------------------------------------------------------------------
# Step 4. Build final SHMI scores using pre-optimized settings
#
# prepare_shmi_inputs():
#   • Reads Excel management sheets (Crop, Disturbance, Amendments, Animals)
#   • Harmonizes crop windows (plant → harvest → termination)
#   • Builds daily grid (crop presence + disturbance)
#   • Computes rotation bounds
#
# build_shmi():
#   • Computes SHMI subindices:
#         Cover, Diversity, InvDist, Animals
#   • Applies optimized weights
#   • Returns indicator_df with SHMI + subindices
# ---------------------------------------------------------------------------

shmi_inputs <- prepare_shmi_inputs("NAPESHM/data/NAPESHM_test.xlsx")

shmi_out <- build_shmi(shmi_inputs, settings = mySettings)

# Merge SHMI with SHMI scores with all meta-data
shmi_full <- shmi_out$indicator_df %>%
  left_join(finalDF, by = "MGT_combo")

saveRDS(shmi_full, "NAPESHM/data/shmi_full.rds")

# ---------------------------------------------------------------------------
# Step 5. Generate NAPESHM figures and analyses (Manuscript Fig. 2)
#
# Fig. 2a — Overall SHMI → lv_SH relationship across all sites
# Fig. 2b — Between-site median-centered SHMI → SH relationship
# Fig. 2c — Within-site deviation SHMI → SH relationship
#
# These correspond to the three manuscript regressions:
#   • Overall linear relationship
#   • Between-site (site medians)
#   • Within-site (site-centered deviations)
# ---------------------------------------------------------------------------

# --- Figure 2a. Overall relationship between SHMI and lv_SH across all sites
res <- lm(SH ~ SHMI, data = shmi_full)
summary(res)

fig2a <- ggplot(shmi_full, aes(x = SHMI, y = SH)) +
  geom_point(aes(fill = land_use),
             color = "grey50", shape = 21, size = 3, alpha = 0.5) +
  geom_smooth(formula = (y ~ x), method = "lm", se = TRUE) +
  stat_poly_eq(formula = (y ~ x), use_label("eq", "R2"), rr.digits = 3, size = 3) +
  labs(
    y = "Soil health (lv_SH)",
    x = "Soil health management index (SHMI)",
    fill = "Land use"
  ) +
  scale_x_continuous(limits = c(0, 100)) +
  scale_y_continuous(limits = c(0, 100)) +
  theme_bw(base_size = 12) +
  ggh4x::force_panelsizes(rows = unit(2, "in"),
                          cols = unit(3, "in"))

fig2a
ggsave(filename="NAPESHM/figs/fig2a.png", plot=fig2a)

# --- Figure 2b. Between-site median-centered relationship
site_means <- shmi_full %>% 
  group_by(sitecode) %>% 
  summarise(
    SH_mean   = median(SH),
    SHMI_mean = median(SHMI),
    land_use  = first(land_use),
    BD        = mean(BD)
  )

lm_between <- lm(SH_mean ~ SHMI_mean, data = site_means)
summary(lm_between)

fig2b <- ggplot(site_means, aes(x = SHMI_mean, y = SH_mean)) +
  geom_point(aes(fill = land_use),
             color = "grey50", shape = 21, size = 3, alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  stat_poly_eq(formula = (y ~ x), use_label("eq", "R2"),
               rr.digits = 3, size = 3) +
  scale_x_continuous(limits = c(0, 100)) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Soil health management index (SHMI)", 
       y = "Soil health (lv_SH)",
       fill = "Land use") +
  theme_bw(base_size = 12) +
  ggh4x::force_panelsizes(rows = unit(2, "in"),
                          cols = unit(3, "in"))

fig2b
ggsave(filename="NAPESHM/figs/fig2b.png", plot=fig2b)

# --- Figure 2c. Within-site median-centered relationship
df_within <- shmi_full %>%
  group_by(sitecode) %>%
  mutate(
    SH_center   = SH   - median(SH,   na.rm = TRUE),
    SHMI_center = SHMI - median(SHMI, na.rm = TRUE)
  ) %>%
  ungroup()

m_within_simple <- lm(SH_center ~ SHMI_center, data = df_within)
summary(m_within_simple)

fig2c <- ggplot(df_within, aes(x = SHMI_center, y = SH_center)) +
  geom_point(aes(fill = land_use),
             color = "grey50", shape = 21, size = 3, alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  stat_poly_eq(formula = (y ~ x), use_label("eq", "R2"),
               rr.digits = 3, size = 3) +
  labs(x = "Δ SHMI (within-site deviation)", 
       y = "Δ lv_SH (within-site deviation)",
       fill = "Land use") +
  theme_bw(base_size = 12) +
  ggh4x::force_panelsizes(rows = unit(2, "in"),
                          cols = unit(3, "in"))

fig2c
ggsave(filename="NAPESHM/figs/fig2c.png", plot=fig2c)

############################################################
# End run_all.R
############################################################