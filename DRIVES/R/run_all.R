############################################################
# run_all_DRIVES.R — Full DRIVES SHMI + Yield Stability Pipeline
#
# This script executes the complete DRIVES workflow:
#   • Load SHMI package + dependencies
#   • Build SHMI scores for each long-term DRIVES experiment
#   • Extract yield + N-rate data
#   • Compute rolling 10-year windows
#   • Fit nonlinear SHMI → yield models (per site + global)
#   • Fit mixed-effects models with AR(1) correlation
#   • Generate yield stability metrics
#   • Generate manuscript figures (global + site-level)
#
# All SHMI calculations use the SHMI package using default weights.
############################################################

#devtools::install_github ("DanielManter-USDA/SHMI")
library(SHMI)
library(dplyr)
library(readxl)
library(ggplot2)
library(purrr)
library(broom)
library(lmerTest)
library(nlme)
library(broom.mixed)
library(mgcv)
library(tidyr)
library(minpack.lm)
library(metafor)
library(scales)
library(mgcv)
library(tibble)
library(metafor)

############################################################
# Study 1 — NELITCSE (1986–2014)
#
# Notes:
#   • Era 1 ends in 2014
#   • Soil disturbance stops in 2014 (all NT afterwards)
#   • All plots managed at optimum N fertility
#   • 600-block excluded due to missing yield 1986–1995
############################################################

exclude <- c(
  "Drives_NELITCSE_601_R1T1P1", "Drives_NELITCSE_602_R3T1P2", "Drives_NELITCSE_603_R4T1P1",
  "Drives_NELITCSE_603_R2T1P1", "Drives_NELITCSE_604_R3T1P1", "Drives_NELITCSE_605_R1T1P1",
  "Drives_NELITCSE_605_R1T3P1", "Drives_NELITCSE_606_R3T1P1", "Drives_NELITCSE_606_R3T3P1",
  "Drives_NELITCSE_607_R4T1P1", "Drives_NELITCSE_607_R2T3P1", "Drives_NELITCSE_608_R3T1P2",
  "Drives_NELITCSE_608_R3T3P2", "Drives_NELITCSE_609_R1T1P1", "Drives_NELITCSE_609_R1T6P1",
  "Drives_NELITCSE_610_R4T1P1", "Drives_NELITCSE_610_R2T6P1", "Drives_NELITCSE_611_R3T1P2",
  "Drives_NELITCSE_611_R3T6P2", "Drives_NELITCSE_612_R3T1P1", "Drives_NELITCSE_612_R3T6P1",
  "Drives_NELITCSE_613_R3T4P1", "Drives_NELITCSE_613_R3T1P1", "Drives_NELITCSE_614_R3T4P2",
  "Drives_NELITCSE_614_R3T1P2", "Drives_NELITCSE_615_R2T4P1", "Drives_NELITCSE_615_R4T1P1",
  "Drives_NELITCSE_616_R1T4P1", "Drives_NELITCSE_616_R1T1P1", "Drives_NELITCSE_617_R3T1P1",
  "Drives_NELITCSE_617_R3T5P1", "Drives_NELITCSE_618_R3T1P2", "Drives_NELITCSE_618_R3T5P2",
  "Drives_NELITCSE_619_R4T1P1", "Drives_NELITCSE_619_R2T5P1", "Drives_NELITCSE_620_R1T1P1",
  "Drives_NELITCSE_620_R1T5P1", "Drives_NELITCSE_621_R1T1P1", "Drives_NELITCSE_621_R1T2P1",
  "Drives_NELITCSE_622_R3T1P2", "Drives_NELITCSE_622_R3T2P2", "Drives_NELITCSE_623_R4T1P1",
  "Drives_NELITCSE_623_R2T2P1", "Drives_NELITCSE_624_R3T1P1", "Drives_NELITCSE_624_R3T2P1"
)

inputs.nel <- prepare_shmi_inputs(
  "DRIVES/data/SHMI_template.NELITCSE.xlsx",
  start_date_override = "1986-01-01",
  end_date_override   = "2014-12-31",
  exclude             = exclude,
  calc_yield          = TRUE,
  calc_n_rate         = FALSE
)

shmi.nel <- build_shmi(inputs.nel)

yield.nel <- inputs.nel$yield %>%
  left_join(shmi.nel$indicator_df, by = c("MGT_combo", "MGT_study", "MGT_farm", "MGT_field", "MGT_trt")) %>% 
  tidyr::extract(
    MGT_trt,
    into = c("rotation", "tillage", "phase"),
    regex = "(R[0-9]+)(T[0-9]+)(P[0-9]+)",
    remove = FALSE
  ) %>% 
  mutate(SHMI = round(SHMI, 2)) %>% 
  mutate(
    MGT_trt_base = sub("P[0-9]+$", "", MGT_trt)
    
  )

window_starts <- 1986:2005
window_ends   <- window_starts + 9  # 10-year window
rolling_yield.nel <- purrr::map2_dfr(window_starts, window_ends, function(s, e) {
  yield.nel %>%
    filter(crop_year >= s, crop_year <= e) %>%
    mutate(window = paste0(s, "-", e))
}) %>% 
  mutate(window = factor(window))

############################################################
# Study 2 — NEMLTCRS (1983–2018)
#
# Notes:
#   • All rotations complete
#   • N-rate treatments explicit (N1, N2, N3)
############################################################

inputs.nem <- prepare_shmi_inputs(
  "DRIVES/data/SHMI_template.NEMLTCRS.xlsx",
  start_date_override = "1983-01-01",
  end_date_override   = "2018-12-31",
  calc_yield          = TRUE,
  calc_n_rate         = TRUE
)

shmi.nem <- build_shmi(inputs.nem)

yield.nem <- inputs.nem$yield %>%
  left_join(inputs.nem$n_rate, by = c("MGT_combo", "MGT_study", "MGT_farm", "MGT_field", "MGT_trt", "crop_year")) %>%
  left_join(shmi.nem$indicator_df, by = c("MGT_combo", "MGT_study", "MGT_farm", "MGT_field", "MGT_trt")) %>% 
  tidyr::extract(
    MGT_trt,
    into = c("rotation", "N_code", "phase"),
    regex = "(R[0-9]+)(N[0-9]+)(P[0-9]+)",
    remove = FALSE
  ) %>% 
  mutate(
    SHMI = round(SHMI, 2)
  ) %>% 
  mutate(
    MGT_trt_base = sub("P[0-9]+$", "", MGT_trt),
    MGT_trt_base = sub("N[0-9]+$", "", MGT_trt_base)
  )

window_starts <- 1983:2009
window_ends   <- window_starts + 9
rolling_yield.nem <- purrr::map2_dfr(window_starts, window_ends, function(s, e) {
  yield.nem %>%
    filter(crop_year >= s, crop_year <= e) %>%
    mutate(window = paste0(s, "-", e))
}) %>% 
  mutate(window = factor(window))

############################################################
# Study 3 — OHTVD (1964–2012)
############################################################

inputs.oh <- prepare_shmi_inputs(
  "DRIVES/data/SHMI_template.OHTVD_Ho.xlsx",
  start_date_override = "1964-01-01",
  end_date_override   = "2012-12-31",
  calc_yield          = TRUE,
  calc_n_rate         = FALSE
)

shmi.oh <- build_shmi(inputs.oh)

yield.oh <- inputs.oh$yield %>%
  left_join(shmi.oh$indicator_df, by = c("MGT_combo", "MGT_study", "MGT_farm", "MGT_field", "MGT_trt")) %>% 
  tidyr::extract(
    MGT_trt,
    into = c("rotation", "tillage", "phase"),
    regex = "(R[0-9]+)(T[0-9]+)(P[0-9]+)",
    remove = FALSE
  ) %>% 
  mutate(SHMI = round(SHMI, 2)) %>% 
  mutate(
    MGT_trt_base = sub("P[0-9]+$", "", MGT_trt)
  )

window_starts <- 1964:2003
window_ends   <- window_starts + 9
rolling_yield.oh <- purrr::map2_dfr(window_starts, window_ends, function(s, e) {
  yield.oh %>%
    filter(crop_year >= s, crop_year <= e) %>%
    mutate(window = paste0(s, "-", e))
}) %>% 
  mutate(window = factor(window))

############################################################
# Study 4 — WICST_ARS (1990–2020)
############################################################

inputs.wicst <- prepare_shmi_inputs(
  "DRIVES/data/SHMI_template.WICST_ARS.xlsx",
  start_date_override = "1990-01-01",
  end_date_override   = "2020-12-31",
  calc_yield          = TRUE,
  calc_n_rate         = FALSE,
  verbose             = TRUE
)

shmi.wicst <- build_shmi(inputs.wicst)

yield.wicst <- inputs.wicst$yield %>%
  left_join(shmi.wicst$indicator_df, by = c("MGT_combo", "MGT_study", "MGT_farm", "MGT_field", "MGT_trt")) %>% 
  tidyr::extract(
    MGT_trt,
    into = c("rotation", "phase"),
    regex = "(CS[0-9]+)(P[0-9]+)",
    remove = FALSE
  ) %>% 
  mutate(
    SHMI = round(SHMI, 2)
  ) %>% 
  mutate(
    MGT_trt_base = sub("P[0-9]+$", "", MGT_trt)
  )

window_starts <- 1990:2011
window_ends   <- window_starts + 9
rolling_yield.wicst <- purrr::map2_dfr(window_starts, window_ends, function(s, e) {
  yield.wicst %>%
    filter(crop_year >= s, crop_year <= e) %>%
    mutate(window = paste0(s, "-", e))
}) %>% 
  mutate(window = factor(window))

############################################################
# Yield Analysis
############################################################

rolling_yield.nel$MGT_field <- as.character(rolling_yield.nel$MGT_field)
rolling_yield.oh$MGT_field <- as.character(rolling_yield.oh$MGT_field)
rolling_yield.wicst$MGT_field <- as.character(rolling_yield.wicst$MGT_field)

rolling_yield_all <- dplyr::bind_rows(
  rolling_yield.nem,
  rolling_yield.nel,
  rolling_yield.oh,
  rolling_yield.wicst) %>%
  filter(CD_name == "corn") %>% 
  group_by(MGT_farm) %>% 
  tidyr::separate(window, into = c("start", "end"), sep = "-", remove = FALSE) %>%
  mutate(
    start = as.numeric(start),
    end   = as.numeric(end),
    window_midpoint = (start + end) / 2
  ) %>% 
  mutate(
    time = window_midpoint - min(window_midpoint, na.rm=TRUE) + 1,
  ) %>%
  ungroup() %>% 
  group_by(MGT_farm, MGT_trt_base, MGT_trt, N_code, MGT_field, window, time, window_midpoint) %>% 
  summarise(
    SHMI = mean(SHMI, na.rm = TRUE),
    yield_kg_ha = mean(yield_kg_ha, na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  mutate(
    yield_Mg_ha = yield_kg_ha / 1000,
    global_field_id = as.factor(paste(MGT_farm, MGT_field, sep = "_"))
  )

rolling_yield_all <- rolling_yield_all %>%
  mutate(
    N_class = case_when(
      
      # NEMLTCRS explicit N treatments
      MGT_farm == "NEMLTCRS" & N_code == "N1" ~ "low",
      MGT_farm == "NEMLTCRS" & N_code == "N2" ~ "medium",
      MGT_farm == "NEMLTCRS" & N_code == "N3" ~ "high",
      
      # NEL & OH always optimum N
      MGT_farm %in% c("NELITCSE", "OHTVD_Ho") ~ "high",
      
      # WICST rotation-based N fertility classes
      MGT_farm == "WICST_ARS" & MGT_trt_base == "CS1" ~ "high",
      MGT_farm == "WICST_ARS" & MGT_trt_base == "CS2" ~ "high",
      MGT_farm == "WICST_ARS" & MGT_trt_base == "CS3" ~ "low",
      MGT_farm == "WICST_ARS" & MGT_trt_base == "CS4" ~ "high",
      MGT_farm == "WICST_ARS" & MGT_trt_base == "CS5" ~ "low",
      
      # fallback
      TRUE ~ NA_character_
    )
  ) %>% 
  mutate(
    N_class = factor(N_class, levels = c("low", "medium", "high"))
  )

rolling_yield_all <- rolling_yield_all %>% 
  filter(time <= 30)

m_yield_nls <- nlsLM(
  yield_Mg_ha ~ Asym + (R0 - Asym) * exp(-k * SHMI),
  data = rolling_yield_all,
  start = list(
    Asym = min(rolling_yield_all$yield_Mg_ha),
    R0   = max(rolling_yield_all$yield_Mg_ha),
    k    = 0.05
  )
)
pars_yield <- coef(m_yield_nls)
pars_yield

rolling_yield_all$exp_SHMI <- pars_yield["Asym"] +
  (pars_yield["R0"] - pars_yield["Asym"]) * exp(-pars_yield["k"] * rolling_yield_all$SHMI)

m_yield_log <- nlme::lme(
  yield_Mg_ha ~ N_class + exp_SHMI * time,
  random = ~ 1 | global_field_id,
  correlation = corAR1(form = ~ time | global_field_id),
  data = rolling_yield_all
)

summary(m_yield_log)
anova(m_yield_log)
MuMIn::r.squaredGLMM(m_yield_log)

global_shmi_grid <- expand.grid(
  SHMI = seq(min(rolling_yield_all$SHMI), 
             max(rolling_yield_all$SHMI), 
             length.out = 300),
  time = mean(rolling_yield_all$time, na.rm = TRUE),
  window_midpoint = unique(rolling_yield_all$window_midpoint),
  N_class = unique(rolling_yield_all$N_class)
)

global_shmi_grid$exp_SHMI <- pars_yield["Asym"] +
  (pars_yield["R0"] - pars_yield["Asym"]) * exp(-pars_yield["k"] * global_shmi_grid$SHMI)

pred <- predict(m_yield_log, newdata = global_shmi_grid, level = 0, se.fit = TRUE)

global_shmi_grid$fit <- pred$fit
global_shmi_grid$se  <- pred$se.fit

global_shmi_plot <- global_shmi_grid %>%
  group_by(SHMI) %>%
  summarize(
    fit = mean(fit),
    se  = sqrt(mean(se^2)),
    .groups = "drop"
  )

agg_global <- rolling_yield_all %>%
  group_by(MGT_farm, MGT_trt_base, N_class) %>%
  summarize(
    SHMI_m = mean(SHMI, na.rm = TRUE),
    SHMI_se = sd(SHMI) / sqrt(n()),
    yield_m = mean(yield_Mg_ha, na.rm = TRUE),
    yield_se = sd(yield_Mg_ha) / sqrt(n()),
    .groups = "drop"
  ) %>% 
  mutate(
    N_class = factor(N_class, levels = c("low", "medium", "high", "low_organic", "high_organic"))
  )

fig3a <- ggplot() +
  geom_line(
    data = global_shmi_plot,
    aes(x = SHMI, y = fit),
    color = "black", linewidth = 1.2
  ) +
  geom_ribbon(
    data = global_shmi_plot,
    aes(x = SHMI, ymin = (fit - 1.96*se), ymax = (fit + 1.96*se)),
    fill = "blue", alpha = 0.15
  ) +
  geom_point(
    data = agg_global,
    aes(x = SHMI_m, y = yield_m, fill = N_class, shape = MGT_farm),
    size = 3, color = "black"
  ) +
  geom_errorbar(
    data = agg_global,
    aes(
      x = SHMI_m,
      ymin = yield_m - yield_se,
      ymax = yield_m + yield_se
    ),
    width = 0,
    orientation = "x"
  ) +
  geom_errorbar(
    data = agg_global,
    aes(
      y = yield_m,
      xmin = SHMI_m - SHMI_se,
      xmax = SHMI_m + SHMI_se
    ),
    width = 0,
    orientation = "y"
  ) +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  scale_fill_manual(values = c(
    "low"          = "#E69F00",
    "medium"       = "#CC79A7",
    "high"         = "#009E73",
    "low_organic"  = "#D55E00",
    "high_organic" = "#006400"
  )) +
  
  ylim(0, 15) +
  labs(
    x = "Soil Health Management Index (SHMI)",
    y = "Corn Yield (Mg ha\u207B\u00B9)",
    fill = "N class",
    shape = "Research Site"
  ) +
  
  theme_bw(base_size = 10) +
  theme(legend.position = "right") +
  ggh4x::force_panelsizes(rows = unit(2, "in"),
                          cols = unit(3, "in")) +
  
  guides(
    fill  = guide_legend(override.aes = list(shape = 21, size = 3)),
    shape = guide_legend(override.aes = list(fill = "grey80", size = 3))
  ) +
  theme(
    legend.box = "horizontal"
  )

fig3a
ggsave(filename="DRIVES/figs/fig3a.png", plot=fig3a)

# SHMI x time graph
plot_data_marginalized <- expand.grid(
  SHMI = seq(min(rolling_yield_all$SHMI, na.rm=TRUE), 
             max(rolling_yield_all$SHMI, na.rm=TRUE), 
             length.out = 200),
  time = c(1, 10, 30),
  window_midpoint = unique(rolling_yield_all$window_midpoint),
  N_class = unique(rolling_yield_all$N_class)
)

plot_data_marginalized$exp_SHMI <- pars_yield["Asym"] +
  (pars_yield["R0"] - pars_yield["Asym"]) * exp(-pars_yield["k"] * plot_data_marginalized$SHMI)

predictions_marginalized        <- predict(m_yield_log, newdata = plot_data_marginalized, level = 0, se.fit = TRUE)
plot_data_marginalized$pred_res <- predictions_marginalized$fit
plot_data_marginalized$se       <- predictions_marginalized$se.fit

plot_data_final <- plot_data_marginalized %>%
  group_by(SHMI, time) %>%
  summarize(
    pred_res = mean(pred_res, na.rm = TRUE),
    se       = mean(se, na.rm = TRUE),
    .groups = "drop"
  )

ribbon_ages <- c(1, 10, 30)

fig3b <- ggplot() +
  geom_line(data = plot_data_final, 
            aes(x = SHMI, y = pred_res, color = factor(time)), 
            linewidth = 1.2) +
  geom_ribbon(data = plot_data_final %>% filter(time %in% ribbon_ages), 
              aes(x = SHMI, 
                  ymin = (pred_res - 1.96 * se), 
                  ymax = (pred_res + 1.96 * se), 
                  fill = factor(time)), 
              alpha = 0.1, color = NA) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  ylim(0, 15) +
  labs(
    x = "Soil Health Management Index (SHMI)",
    y = "Corn Yield (Mg ha\u207B\u00B9)",
    # subtitle = "B) Interaction between SHMI and System Age",
    color = "System Age",
    fill = "System Age"
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "right") +
  ggh4x::force_panelsizes(rows = unit(2, "in"),
                          cols = unit(3, "in"))

fig3b
ggsave(filename="DRIVES/figs/fig3b.png", plot=fig3b)

# --- individual site analysis
# NELITCSE
dat <- filter(rolling_yield_all, MGT_farm == "NELITCSE")
m_yield_nls_nel <- nlsLM(
  yield_Mg_ha ~ Asym + (R0 - Asym) * exp(-k * SHMI),
  data = dat,
  start = list(
    Asym = min(dat$yield_Mg_ha),
    R0   = max(dat$yield_Mg_ha),
    k    = 0.05
  )
)
pars_nel <- coef(m_yield_nls_nel)
dat$exp_SHMI_NELITCSE <- pars_nel["Asym"] +
  (pars_nel["R0"] - pars_nel["Asym"]) * exp(-pars_nel["k"] * dat$SHMI)
m_yield_nel <- nlme::lme(
  yield_Mg_ha ~ exp_SHMI_NELITCSE * time,
  random = ~ 1 | global_field_id,
  correlation = corAR1(form = ~ time | global_field_id),
  data = dat
)
summary(m_yield_nel)
anova(m_yield_nel)
MuMIn::r.squaredGLMM(m_yield_nel)

# NEMLTCRS
dat <- filter(rolling_yield_all, MGT_farm == "NEMLTCRS") %>% 
  droplevels()
m_yield_nls_nem <- nlsLM(
  yield_Mg_ha ~ Asym + (R0 - Asym) * exp(-k * SHMI),
  data = dat,
  start = list(
    Asym = min(dat$yield_Mg_ha),
    R0   = max(dat$yield_Mg_ha),
    k    = 0.05
  )
)
pars_nem <- coef(m_yield_nls_nem)
dat$exp_SHMI_NEMLTCRS <- pars_nem["Asym"] +
  (pars_nem["R0"] - pars_nem["Asym"]) * exp(-pars_nem["k"] * dat$SHMI)
m_yield_nem <- nlme::lme(
  yield_Mg_ha ~ N_class + exp_SHMI_NEMLTCRS * time,
  random = ~ 1 | global_field_id,
  correlation = corAR1(form = ~ time | global_field_id),
  data = dat
)
summary(m_yield_nem)
anova(m_yield_nem)
MuMIn::r.squaredGLMM(m_yield_nem)

# OHTVD_Ho
dat <- filter(rolling_yield_all, MGT_farm == "OHTVD_Ho") %>% 
  droplevels()
m_yield_nls_oh <- nlsLM(
  yield_Mg_ha ~ Asym + (R0 - Asym) * exp(-k * SHMI),
  data = dat,
  start = list(
    Asym = min(dat$yield_Mg_ha),
    R0   = max(dat$yield_Mg_ha),
    k    = 0.05
  )
)
pars_oh <- coef(m_yield_nls_oh)
dat$exp_SHMI_OHTVD_Ho <- pars_oh["Asym"] +
  (pars_oh["R0"] - pars_oh["Asym"]) * exp(-pars_oh["k"] * dat$SHMI)
m_yield_oh <- nlme::lme(
  yield_Mg_ha ~ exp_SHMI_OHTVD_Ho * time,
  random = ~ 1 | global_field_id,
  correlation = corAR1(form = ~ time | global_field_id),
  data = dat
)
summary(m_yield_oh)
anova(m_yield_oh)
MuMIn::r.squaredGLMM(m_yield_oh)

# WICST_ARS
dat <- filter(rolling_yield_all, MGT_farm == "WICST_ARS") %>% 
  droplevels()
m_yield_nls_wi <- nlsLM(
  yield_Mg_ha ~ Asym + (R0 - Asym) * exp(-k * SHMI),
  data = dat,
  start = list(
    Asym = min(dat$yield_Mg_ha),
    R0   = max(dat$yield_Mg_ha),
    k    = 0.05
  )
)
pars_wi <- coef(m_yield_nls_wi)
dat$exp_SHMI_WICST_ARS <- pars_wi["Asym"] +
  (pars_wi["R0"] - pars_wi["Asym"]) * exp(-pars_wi["k"] * dat$SHMI)
m_yield_wi <- nlme::lme(
  yield_Mg_ha ~ N_class + exp_SHMI_WICST_ARS * time,
  random = ~ 1 | global_field_id,
  correlation = corAR1(form = ~ time | global_field_id),
  data = dat
)
summary(m_yield_wi)
anova(m_yield_wi)
MuMIn::r.squaredGLMM(m_yield_wi)

# Models and data objects
models <- list(
  NELITCSE = m_yield_nel,
  NEMLTCRS = m_yield_nem,
  OHTVD_Ho = m_yield_oh,
  WICST_ARS = m_yield_wi
)

get_shmi_effect_lme <- function(mod, data, pars, exp_name,
                                shmi_var = "SHMI", time_var = "time",
                                probs = c(0.25, 0.75), nclass_value = "high") {
  
  # 1. SHMI quantiles
  qvals <- quantile(data[[shmi_var]], probs = probs, na.rm = TRUE)
  
  # 2. Mean time
  tmean <- mean(data[[time_var]], na.rm = TRUE)
  
  # 3. Extract N_class levels from the model (NOT the data)
  mf <- model.frame(mod)
  has_nclass <- "N_class" %in% names(mf)
  
  if (has_nclass) {
    levels_N <- levels(mf$N_class)
    chosen <- if (nclass_value %in% levels_N) nclass_value else levels_N[1]
  }
  
  # 4. Build newdata with ONLY predictors used in the model
  newdat <- data.frame(
    time = rep(tmean, 2)
  )
  
  if (has_nclass) {
    newdat$N_class <- factor(chosen, levels = levels_N)
  }
  
  # 5. Compute exp_SHMI using correct variable name
  newdat[[exp_name]] <- pars["Asym"] +
    (pars["R0"] - pars["Asym"]) * exp(-pars["k"] * qvals)
  
  # 6. Predict fixed effects only
  pred <- predict(mod, newdata = newdat, level = 0)
  
  # 7. Contrast (75th - 25th)
  diff_est <- pred[2] - pred[1]
  
  # 8. SE using fixed-effects vcov
  X <- model.matrix(delete.response(terms(mod)), newdat)
  V <- summary(mod)$varFix
  L <- X[2, ] - X[1, ]
  diff_se <- sqrt(as.numeric(t(L) %*% V %*% L))
  
  # 9. Return tidy tibble
  tibble(
    q25   = qvals[1],
    q75   = qvals[2],
    effect = diff_est,
    se     = diff_se,
    lower  = diff_est - 1.96 * diff_se,
    upper  = diff_est + 1.96 * diff_se
  )
}

global_eff <- get_shmi_effect_lme(
  mod      = m_yield_log,
  data     = rolling_yield_all,
  pars     = pars_yield,
  exp_name = "exp_SHMI"
)

pars_list <- list(
  NELITCSE = pars_nel,
  NEMLTCRS = pars_nem,
  OHTVD_Ho = pars_oh,
  WICST_ARS = pars_wi
)

site_summaries <- imap_dfr(models, function(mod, site) {
  
  d_site <- filter(rolling_yield_all, MGT_farm == site)
  n_obs <- nrow(d_site)
  
  exp_name <- paste0("exp_SHMI_", site)
  
  eff <- tryCatch(
    get_shmi_effect_lme(
      mod      = mod,
      data     = d_site,
      pars     = pars_list[[site]],
      exp_name = exp_name
    ),
    error = function(e)
      tibble(q25=NA,q75=NA,effect=NA,se=NA,lower=NA,upper=NA)
  )
  
  tibble(
    site = site,
    n_obs = n_obs,
    q25 = eff$q25, q75 = eff$q75,
    effect = eff$effect, se = eff$se,
    lower = eff$lower, upper = eff$upper
  )
})

df_meta <- site_summaries %>% filter(!is.na(effect) & !is.na(se))
res_re <- rma(yi = effect, sei = se, data = df_meta, method = "REML")
pooled_est <- as.numeric(res_re$b)
pooled_se  <- sqrt(res_re$se^2)
pooled_ci  <- c(pooled_est - 1.96 * pooled_se, pooled_est + 1.96 * pooled_se)

pooled_row <- tibble(site = "Pooled (REML)", n_obs = NA_integer_,
                     q25 = NA_real_, q75 = NA_real_,
                     effect = pooled_est, se = pooled_se, lower = pooled_ci[1], upper = pooled_ci[2])
global_row <- tibble(site = "Global Model", n_obs = nrow(rolling_yield_all),
                     q25 = global_eff$q25, q75 = global_eff$q75,
                     effect = global_eff$effect, se = global_eff$se, lower = global_eff$lower, upper = global_eff$upper)

df_plot <- bind_rows(site_summaries, pooled_row, global_row)

new_levels <- c("Global Model", "Pooled (REML)", "NELITCSE", "NEMLTCRS", "OHTVD_Ho", "WICST_ARS")

df_plot <- df_plot %>%
  mutate(p_val = 2 * (1 - pnorm(abs(effect / se))),
         sig = case_when(
           is.na(p_val) ~ "",
           p_val < 0.001 ~ " *** ",
           p_val < 0.01  ~ " ** ",
           p_val < 0.05  ~ " * ",
           p_val < 0.1   ~ " + ",
           TRUE          ~ "")
  )

df_plot <- df_plot %>%
  mutate(site = factor(as.character(site), levels = new_levels)) %>%
  mutate(site = forcats::fct_rev(site)) %>% 
  mutate(effect = effect,
         se     = se,
         lower  = lower,
         upper  = upper)

fig3c <- ggplot(df_plot, aes(x = effect, y = site)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_errorbar(data = filter(df_plot, !(site %in% c("Global Model","Pooled (REML)"))),
                aes(y = site, xmin = lower, xmax = upper), width = 0.18, color = "grey40",
                orientation = "y") +
  geom_point(data = filter(df_plot, !(site %in% c("Global Model","Pooled (REML)"))),
             aes(size = pmax(2, sqrt(n_obs))), shape = 21, fill = "steelblue", color = "black", stroke = 0.3) +
  # annotate significance stars to the right of CI
  geom_text(data = filter(df_plot, site != "Global Model"),
            aes(x = upper + 0.03 * (max(df_plot$upper, na.rm = TRUE) - min(df_plot$lower, na.rm = TRUE)),
                label = sig),
            hjust = 0, size = 3, color = "black") +
  # Global row (emphasized)
  geom_errorbar(data = filter(df_plot, site == "Global Model"),
                aes(y = site, xmin = lower, xmax = upper), width = 0.12, 
                color = "darkgreen", linewidth = 0.8, orientation = "y") +
  geom_point(data = filter(df_plot, site == "Global Model"),
             aes(x = effect, y = site), shape = 23, size = 5, fill = "darkgreen", color = "black") +
  geom_text(data = filter(df_plot, site == "Global Model"),
            aes(x = upper + 0.03 * (max(upper, na.rm = TRUE) - min(lower, na.rm = TRUE)),
                label = sig),
            hjust = 0, size = 3.5, fontface = "bold", color = "darkgreen") +
  # Pooled row (distinct)
  geom_errorbar(data = filter(df_plot, site == "Pooled (REML)"),
                aes(y = site, xmin = lower, xmax = upper), 
                width = 0.10, color = "firebrick", linewidth = 0.8,
                orientation = "y") +
  geom_point(data = filter(df_plot, site == "Pooled (REML)"),
             aes(x = effect, y = site), shape = 23, size = 5.5, fill = "firebrick", color = "black") +
  scale_x_continuous(name = "Corn Yield (Mg ha\u207B\u00B9)", expand = c(0, 0.15)) +
  scale_size_continuous(range = c(2.5, 6), guide = "none") +
  theme_bw(base_size = 10) +
  labs(
    # subtitle = "C) Predicted yield differnce (75th-25th SHMI)"
  ) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 9, face = ifelse(levels(df_plot$site) == "Global Model", "bold", "plain")),
        panel.grid.major.y = element_blank()) +
  ggh4x::force_panelsizes(rows = unit(2, "in"),
                          cols = unit(3, "in"))  

fig3c
ggsave(filename="DRIVES/figs/fig3c.png", plot=fig3c)


############################################################
# Yield Stability Analysis
############################################################

yield.nel$MGT_field <- as.character(yield.nel$MGT_field)
yield.nem$MGT_field <- as.character(yield.nem$MGT_field)
yield.oh$MGT_field <- as.character(yield.oh$MGT_field)
yield.wicst$MGT_field <- as.character(yield.wicst$MGT_field)

yield_all <- dplyr::bind_rows(
  yield.nem,
  yield.nel,
  yield.oh,
  yield.wicst) %>% 
  filter(
    CD_name == "corn"
  ) %>% 
  group_by(MGT_farm) %>% 
  mutate(
    time = crop_year - min(crop_year, na.rm=TRUE) + 1,
  ) %>% 
  ungroup() %>% 
  mutate(
    yield_Mg_ha = yield_kg_ha / 1000,
    global_field_id = as.factor(paste(MGT_farm, MGT_field, sep = "_"))
  )

yield_all <- yield_all %>%
  mutate(
    N_class = case_when(
      
      # NEMLTCRS explicit N treatments
      MGT_farm == "NEMLTCRS" & N_code == "N1" ~ "low",
      MGT_farm == "NEMLTCRS" & N_code == "N2" ~ "medium",
      MGT_farm == "NEMLTCRS" & N_code == "N3" ~ "high",
      
      # NEL & OH always optimum N
      MGT_farm %in% c("NELITCSE", "OHTVD_Ho") ~ "high",
      
      # WICST rotation-based N fertility classes
      MGT_farm == "WICST_ARS" & MGT_trt_base == "CS1" ~ "high",
      MGT_farm == "WICST_ARS" & MGT_trt_base == "CS2" ~ "high",
      MGT_farm == "WICST_ARS" & MGT_trt_base == "CS3" ~ "low",
      MGT_farm == "WICST_ARS" & MGT_trt_base == "CS4" ~ "high",
      MGT_farm == "WICST_ARS" & MGT_trt_base == "CS5" ~ "low",
      
      # fallback
      TRUE ~ NA_character_
    )
  ) %>% 
  mutate(
    N_class = factor(N_class, levels = c("low", "low_organic", "medium", "high", "high_organic"))
  ) %>% 
  filter(
    !is.na(yield_kg_ha)
  )

# fit model to detrend yield
m_annual <- nlme::lme(
  yield_Mg_ha ~ N_class + time,
  random = ~ 1 | global_field_id,
  correlation = corAR1(form = ~ time | global_field_id),
  data = yield_all
)

# extract residuals and add to data-frame (note 17 samples were missing yield)
yield_all$resid_yield <- resid(m_annual)

yield_res <- yield_all %>%
  group_by(MGT_farm, MGT_trt, MGT_trt_base, N_class, N_code, MGT_field, global_field_id) %>%
  summarise(
    SHMI = mean(SHMI, na.rm = TRUE),
    resid_cv = sqrt(var(resid_yield)),
    .groups = "drop"
  ) %>% 
  filter(
    !is.na(resid_cv)
  )

m_global_res <- nlme::lme(
  resid_cv ~ N_class + SHMI,
  random = ~ 1 | global_field_id,
  data = yield_res
)
summary(m_global_res)
anova(m_global_res)
MuMIn::r.squaredGLMM(m_global_res)


# SHMI global effect - resist
global_shmi_grid <- expand.grid(
  SHMI = seq(min(yield_res$SHMI),
             max(yield_res$SHMI),
             length.out = 300),
  N_class = unique(yield_res$N_class),
  global_field_id = unique(yield_res$global_field_id)
)

pred <- predict(m_global_res, newdata = global_shmi_grid, level = 0, se.fit = TRUE)

global_shmi_grid$fit <- pred$fit
global_shmi_grid$se  <- pred$se.fit

global_shmi_plot <- global_shmi_grid %>%
  group_by(SHMI) %>%
  summarize(
    fit = mean(fit),
    se  = mean(se),
    .groups = "drop"
  )

agg_data <- yield_res %>%
  group_by(MGT_farm, MGT_trt_base, N_class) %>%
  summarize(
    SHMI_m = mean(SHMI),
    SHMI_se = sd(SHMI) / sqrt(n()),
    resid_cv_m = mean(resid_cv),
    resid_cv_se = sd(resid_cv) / sqrt(n()),
    .groups = "drop"
  )


fig4a <- ggplot() +
  geom_point(data = agg_data, aes(x=SHMI_m, y=resid_cv_m / 1000, fill=N_class, shape=MGT_farm),
             size = 3, color = "black") +
  geom_errorbar(
    data = agg_data,
    aes(
      x = SHMI_m,
      ymin = (resid_cv_m - resid_cv_se) / 1000,
      ymax = (resid_cv_m + resid_cv_se) / 1000
    ),
    width = 0,
    orientation = "x"
  ) +
  geom_errorbar(
    data = agg_data,
    aes(
      y = resid_cv_m / 1000,
      xmin = SHMI_m - SHMI_se,
      xmax = SHMI_m + SHMI_se
    ),
    width = 0,
    orientation = "y"
  ) +
  scale_shape_manual(values = c(21,22,23,24) ) +
  scale_fill_manual(values = c(
    "low"          = "#E69F00",
    "medium"       = "#CC79A7",
    "high"         = "#009E73"
  )) +
  geom_line(data = global_shmi_plot,
            aes(x = SHMI, y = fit / 1000),
            color = "black", linewidth = 1.2) +
  geom_ribbon(data = global_shmi_plot,
              aes(x = SHMI, ymin = (fit - 1.96*se) / 1000, ymax = (fit + 1.96*se) / 1000),
              fill = "blue", alpha = 0.15) +
  labs(
    x = "Soil Health Management Index (SHMI)",
    y = "Yield Variability (Mg ha\u207B\u00B9)",
    fill = "N class",
    shape = "Research Site"
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "right") +
  guides(
    fill  = guide_legend(override.aes = list(shape = 21, size = 3)),
    shape = guide_legend(override.aes = list(fill = "grey80", size = 3))
  ) +
  ggh4x::force_panelsizes(rows = unit(2, "in"),
                          cols = unit(3, "in")) +
  theme(
    legend.box = "horizontal"
  ) 

fig4a
ggsave(filename="DRIVES/figs/fig4a.png", plot=fig4a)


# Downside Risk Frequency
p10_by_group <- yield_all %>%
  group_by(MGT_farm, N_class) %>%
  summarise(
    p10 = quantile(yield_Mg_ha, 0.10, na.rm = TRUE),
    .groups = "drop"
  )

yield_all2 <- yield_all %>%
  filter(!is.na(yield_Mg_ha)) %>%
  left_join(p10_by_group, by = c("MGT_farm", "N_class")) %>%
  mutate(bad10 = yield_Mg_ha < p10)


failure_dat <- yield_all2 %>%
  group_by(MGT_farm, MGT_trt_base, N_class, global_field_id) %>%
  summarise(
    SHMI = median(SHMI, na.rm = TRUE),
    failure_rate10 = mean(bad10, na.rm = TRUE),
    n_years = n(),
    .groups = "drop"
  )

m_log_fail <- nlsLM(
  failure_rate10 ~ Asym + (R0 - Asym) * exp(-k * SHMI),
  data = failure_dat,
  start = list(
    Asym = min(failure_dat$failure_rate10, na.rm = TRUE),
    R0   = max(failure_dat$failure_rate10, na.rm = TRUE),
    k    = 0.5
  )
)

pars_fail <- coef(m_log_fail)
pars_fail

failure_dat$log_SHMI <- pars_fail["Asym"] + 
  (pars_fail["R0"] - pars_fail["Asym"]) * exp(-pars_fail["k"] * failure_dat$SHMI)

m_fail_log <- nlme::lme(
  failure_rate10 ~ N_class + log_SHMI,
  random = ~ 1 | global_field_id,
  data = failure_dat
)
summary(m_fail_log)
anova(m_fail_log)
MuMIn::r.squaredGLMM(m_fail_log)

global_shmi_grid <- expand.grid(
  SHMI = seq(min(failure_dat$SHMI), max(failure_dat$SHMI), length.out = 300),
  N_class = unique(failure_dat$N_class),
  global_field_id = unique(failure_dat$global_field_id)
)

global_shmi_grid$log_SHMI <- pars_fail["Asym"] + 
  (pars_fail["R0"] - pars_fail["Asym"]) * exp(-pars_fail["k"] * global_shmi_grid$SHMI)


pred <- predict(m_fail_log, newdata = global_shmi_grid, level = 0, se.fit = TRUE)

global_shmi_grid$fit <- pred$fit
global_shmi_grid$se  <- pred$se.fit

global_shmi_plot <- global_shmi_grid %>%
  group_by(SHMI) %>%
  summarize(
    fit = mean(fit),
    se  = mean(se),
    .groups = "drop"
  )

agg_data <- failure_dat %>%
  group_by(MGT_farm, MGT_trt_base, N_class) %>%
  summarize(
    SHMI_m = mean(SHMI),
    SHMI_se = sd(SHMI) / sqrt(n()),
    failure_rate10_m = mean(failure_rate10),
    failure_rate10_se = sd(failure_rate10) / sqrt(n()),
    .groups = "drop"
  )


fig4b <- ggplot() +
  geom_point(data = agg_data, aes(x=SHMI_m, y=failure_rate10_m, fill=N_class, shape=MGT_farm),
             size = 3, color = "black") +
  geom_errorbar(
    data = agg_data,
    aes(
      x = SHMI_m,
      ymin = failure_rate10_m - failure_rate10_se,
      ymax = failure_rate10_m + failure_rate10_se
    ),
    width = 0,
    orientation = "x"
  ) +
  geom_errorbar(
    data = agg_data,
    aes(
      y = failure_rate10_m,
      xmin = SHMI_m - SHMI_se,
      xmax = SHMI_m + SHMI_se
    ),
    width = 0,
    orientation = "y"
  ) +
  scale_shape_manual(values = c(21,22,23,24) ) +
  scale_fill_manual(values = c(
    "low"          = "#E69F00",
    "medium"       = "#CC79A7",
    "high"         = "#009E73"
  )) +
  geom_line(data = global_shmi_plot,
            aes(x = SHMI, y = fit),
            color = "black", linewidth = 1.2) +
  geom_ribbon(data = global_shmi_plot,
              aes(x = SHMI, ymin = fit - 1.96*se, ymax = fit + 1.96*se),
              fill = "blue", alpha = 0.15) +
  labs(
    x = "Soil Health Management Index (SHMI)",
    y = "Downside-Risk Frequency",
    fill = "N class", 
    shape = "Research Site"
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "right") +
  guides(
    fill  = guide_legend(override.aes = list(shape = 21, size = 3)),
    shape = guide_legend(override.aes = list(fill = "grey80", size = 3))
  ) +
  ggh4x::force_panelsizes(rows = unit(2, "in"),
                          cols = unit(3, "in")) +
  theme(
    legend.box = "horizontal"
  ) 

fig4b
ggsave(filename="DRIVES/figs/fig4b.png", plot=fig4b)

# --- stability
dat <- read.csv("DRIVES/data/TerraClimate_StressMetrics.csv")

stability <- yield_all %>% 
  left_join(
    dat %>% 
      dplyr::select(Site, Year, GS_PDSI, DroughtScore, CompoundStress),
    by = c("MGT_farm" = "Site", "crop_year" = "Year")
  ) %>% 
  group_by(MGT_farm) %>% 
  mutate(
    yr_class = case_when(
      DroughtScore > quantile(DroughtScore, 0.75, na.rm = TRUE) ~ "Bad",
      DroughtScore < quantile(DroughtScore, 0.25, na.rm = TRUE) ~ "Good",
      TRUE ~ "Average"
    )
  ) %>%
  ungroup()

stability <- stability %>%
  group_by(MGT_farm, MGT_trt, MGT_trt_base, N_class, global_field_id) %>%
  summarise(
    SHMI = mean(SHMI),
    median_yield = median(yield_Mg_ha, na.rm = TRUE),
    drought_yield = median(yield_Mg_ha[yr_class == "Bad"], na.rm = TRUE),
    good_yield    = median(yield_Mg_ha[yr_class == "Good"], na.rm = TRUE),
    retention     = drought_yield / median_yield,
    responsiveness = good_yield / median_yield,
    .groups = "drop"
  )

m_log <- nls(
  retention ~ L / (1 + exp(-k * (SHMI - x0))),
  data = stability,
  start = list(
    L  = 1.0,          # retention max is ~1
    k  = 0.05,         # slope guess
    x0 = mean(stability$SHMI)  # inflection near middle
  )
)

summary(m_log)

pars <- coef(m_log)
pars
stability$log_SHMI <- pars["L"] / (1 + exp(-pars["k"] * (stability$SHMI - pars["x0"])))

m_ret_log <- nlme::lme(
  retention ~ N_class + log_SHMI,
  random = ~ 1|global_field_id,
  data = stability
)
summary(m_ret_log)
anova(m_ret_log)
MuMIn::r.squaredGLMM(m_ret_log)

global_shmi_grid <- expand.grid(
  SHMI = seq(min(stability$SHMI), max(stability$SHMI), length.out = 300),
  N_class = unique(stability$N_class),
  global_field_id = unique(stability$global_field_id)
)

global_shmi_grid$log_SHMI <- pars["L"] / (1 + exp(-pars["k"] * (global_shmi_grid$SHMI - pars["x0"])))

pred <- predict(m_ret_log, newdata = global_shmi_grid, level = 0, se.fit = TRUE)

global_shmi_grid$fit <- pred$fit
global_shmi_grid$se  <- pred$se.fit

global_shmi_plot <- global_shmi_grid %>%
  group_by(SHMI) %>%
  summarize(
    fit = mean(fit),
    se  = mean(se),
    .groups = "drop"
  )

agg_data <- stability %>%
  group_by(MGT_farm, MGT_trt_base, N_class) %>%
  summarize(
    SHMI_m = mean(SHMI),
    SHMI_se = sd(SHMI) / sqrt(n()),
    retention_m = mean(retention),
    retention_se = sd(retention) / sqrt(n()),
    .groups = "drop"
  )

fig4c <- ggplot() +
  geom_point(data = agg_data,
             aes(x = SHMI_m, y = retention_m,
                 fill = N_class, shape = MGT_farm),
             size = 3, color = "black") +
  geom_errorbar(
    data = agg_data,
    aes(
      x = SHMI_m,
      ymin = retention_m - retention_se,
      ymax = retention_m + retention_se
    ),
    width = 0,
    orientation = "x"
  ) +
  geom_errorbar(
    data = agg_data,
    aes(
      y = retention_m,
      xmin = SHMI_m - SHMI_se,
      xmax = SHMI_m + SHMI_se
    ),
    width = 0,
    orientation = "y"
  ) +
  scale_shape_manual(values = c(21,22,23,24)) +
  scale_fill_manual(values = c(
    "low"          = "#E69F00",
    "medium"       = "#CC79A7",
    "high"         = "#009E73"
  )) +
  geom_line(data = global_shmi_plot,
            aes(x = SHMI, y = fit),
            color="black", linewidth = 1.2) +
  
  geom_ribbon(data = global_shmi_plot,
              aes(x = SHMI,
                  ymin = fit - 1.96*se,
                  ymax = fit + 1.96*se),
              fill = "blue", alpha = 0.15) +
  
  labs(
    x = "Soil Health Management Index (SHMI)",
    y = "Drought-Year Retention",
    fill = "N class", 
    shape = "Research Site"
  ) +
  
  theme_bw(base_size = 10) +
  theme(legend.position = "right") +
  
  guides(
    fill  = guide_legend(override.aes = list(shape = 21, size = 3)),
    shape = guide_legend(override.aes = list(fill = "grey80", size = 3))
  ) +
  
  ggh4x::force_panelsizes(rows = unit(2, "in"),
                          cols = unit(3, "in")) +
  
  theme(legend.box = "horizontal")

fig4c
ggsave(filename="DRIVES/figs/fig4c.png", plot=fig4c)

############################################################
# End run_all.R
############################################################