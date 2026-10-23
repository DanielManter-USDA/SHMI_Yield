################################################################################
# Compare_v110_v120.R
#
# SHMI v1.1.0 against the current v1.2.0 candidate (all-data calibration),
# overall and for each sub-index, on the same (rebuilt) NAPESHM workbook, so
# differences come only from definitions and weights.
#
# v1.1.0 is reproduced with SHMI 1.2.0 in legacy mode (verified earlier to give
# identical scores to the installed 1.1.0):
#   events window, no intensive-tillage end, animals in their start year only,
#   no clipping of crops to the window, gamma Diversity (max_div 10),
#   max_stir 342, v1.1.0 seasonal, organic-input and pillar weights.
#   NAPESHM records STIR only, so both versions use STIR (v1.1.0's official
#   method was EPA).
#
# v1.2.0 candidate: calendar-year window, intensive-tillage end, animals over
#   the whole grazing period, crops clipped to the window, max_stir 317.5,
#   gamma Diversity (cap 8), Organic Inputs by the "either" rule (a year counts
#   with an amendment or animals), and the weights of the all-data calibration
#   (Outputs/Calibration/all_floor010_results.rds). The Diversity parts
#   (within-year and rotation) are reported alongside.
#
# Output: Outputs/Summary/F13_v110_vs_v120.png, v110_vs_v120_scores.csv,
#         v110_vs_v120_summary.csv
################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(here)
library(SHMI)

OUT      <- function(...) here::here("Outputs", ...)
WORKBOOK <- here::here("Inputs", "NAPESHM_rebuilt.xlsx")   # adjust if named differently
RESULTS  <- OUT("Calibration", "all_floor010_results.rds")
dir.create(OUT("Summary"), showWarnings = FALSE)

# ---- v1.1.0 (legacy mode) ------------------------------------------------------
inp_v110 <- prepare_shmi_inputs(WORKBOOK, rotation_window = "events",
                                tillage_end = "none", verbose = FALSE)
v110 <- build_shmi(inp_v110, expert_mode = TRUE, settings = list(
  dist_meth = "STIR", max_stir = 342, ti_rep = "max",
  hill = 1, max_div = 10,
  w_winter = 0.1259, w_spring = 0.1260, w_summer = 0.3755, w_fall = 0.3726,
  orginput_rule = "weighted", w_amend = 0.6615, w_animal = 0.3385, animal_presence = "start",
  clip_to_rotation = FALSE,
  w_cover = 0.4481, w_diversity = 0.0904, w_invdist = 0.1431, w_orginput = 0.3184
))$indicator_df %>%
  select(MGT_combo, SHMI, Cover, Diversity, InvDist, OrgInput)

# ---- v1.2.0 candidate ------------------------------------------------------------
res <- readRDS(RESULTS)
st  <- res$final_settings                   # seasonal and pillar weights from the calibration

inp_v120 <- readRDS(OUT("shmi_inputs_v120.rds"))
comp <- compute_shmi_components(inp_v120, settings = list(
  dist_meth = "STIR", max_stir = 317.5, ti_rep = "max",
  hill = 1, max_div = 8, animal_presence = "span", orginput_rule = "either"))
dc <- compute_diversity_components(inp_v120$crop, inp_v120$rot_bounds, max_gamma = 8)

v120 <- comp %>%
  left_join(dc %>% select(MGT_combo, Gamma, Gamma_alpha_part, Gamma_beta_part), by = "MGT_combo") %>%
  mutate(across(c(Gamma, Gamma_alpha_part, Gamma_beta_part), ~ coalesce(.x, 0)),
         Cover     = st$w_winter * C_winter + st$w_spring * C_spring +
                     st$w_summer * C_summer + st$w_fall * C_fall,
         Diversity = Gamma,
         OrgInput  = O_any,
         SHMI      = st$w_cover * Cover + st$w_diversity * Diversity +
                     st$w_invdist * InvDist + st$w_orginput * OrgInput)
write.csv(v120 %>% select(MGT_combo, Diversity, Gamma_alpha_part, Gamma_beta_part),
          OUT("Summary", "v120_diversity_parts.csv"), row.names = FALSE)
v120 <- v120 %>% select(MGT_combo, SHMI, Cover, Diversity, InvDist, OrgInput)

# ---- Units, groups, and the comparison -------------------------------------------
units <- readRDS(OUT("sem_all_scores.rds")) %>%
  transmute(MGT_combo,
            group = case_when(system == "perennial" ~ "perennial",
                              grepl("^MX", sitecode) ~ "Mexico cropland",
                              TRUE ~ "US cropland"))

long <- inner_join(
  v110 %>% pivot_longer(-MGT_combo, names_to = "index", values_to = "v110"),
  v120 %>% pivot_longer(-MGT_combo, names_to = "index", values_to = "v120"),
  by = c("MGT_combo", "index")) %>%
  inner_join(units, by = "MGT_combo") %>%
  mutate(index = factor(index, c("SHMI", "Cover", "Diversity", "InvDist", "OrgInput")))

summary_tab <- long %>%
  group_by(index) %>%
  summarise(units = n(),
            pearson = cor(v110, v120), spearman = cor(v110, v120, method = "spearman"),
            mean_v110 = mean(v110), mean_v120 = mean(v120),
            mean_change = mean(v120 - v110),
            units_moved_gt10 = sum(abs(v120 - v110) > 10), .groups = "drop")
print(summary_tab, width = Inf)

labels <- summary_tab %>%
  mutate(label = sprintf("r = %.2f\nmean change %+.1f", pearson, mean_change))

p <- ggplot(long, aes(v110, v120)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_point(aes(colour = group), alpha = 0.6, size = 1.4) +
  geom_text(data = labels, aes(x = 0, y = 100, label = label),
            hjust = 0, vjust = 1, size = 3.2, inherit.aes = FALSE) +
  facet_wrap(~ index, nrow = 2) +
  coord_equal(xlim = c(0, 100), ylim = c(0, 100)) +
  labs(x = "SHMI v1.1.0 (legacy settings)", y = "SHMI v1.2.0 candidate (all-data calibration)",
       colour = NULL, title = "SHMI v1.1.0 vs the v1.2.0 candidate, overall and by sub-index",
       caption = "Dashed line = 1:1. Same rebuilt NAPESHM workbook for both versions; STIR disturbance.") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

dir.create(OUT("Summary"), showWarnings = FALSE)
ggsave(OUT("Summary", "F13_v110_vs_v120.png"), p, width = 10, height = 7.5, dpi = 300)
write.csv(long, OUT("Summary", "v110_vs_v120_scores.csv"), row.names = FALSE)
write.csv(summary_tab, OUT("Summary", "v110_vs_v120_summary.csv"), row.names = FALSE)
cat("\nSaved: Outputs/Summary/F13_v110_vs_v120.png, v110_vs_v120_scores.csv, v110_vs_v120_summary.csv\n")
