################################################################################
# C1e_labels_vs_SHMI.R
#
# If practices can be shown not to add up from practice labels and measured soil health
# (C1d), why score the record? Because the labels cannot predict soil health where it
# was not measured. Same systems as C1d (annual NAPESHM systems with all four labels;
# 277 at 66 sites), same leave-sites-out folds as the calibration (10 folds, 100 repeats,
# seed 2026 + repeat). In each fold, on the training sites only:
#   SHMI                          lv_SH ~ SHMI (published weights held fixed)
#   Practice labels               lv_SH ~ cover + organic + reduced + rotation
#   Practice labels + interactions  the same plus the six pairwise interactions
#   SHMI + labels                 lv_SH ~ SHMI + the four labels
# and the held-out sites are predicted. Measures as Fig. 2a: correlation within sites,
# share of sites (>= 3 systems) with a positive within-site rank correlation, correlation
# between site means, R2. Medians and 95% ranges over repeats; the share of repeats in
# which SHMI beats each alternative.
# SHMI's weights were calibrated on all sites, so its line here is slightly optimistic;
# in A4a refitting the weights in every fold made no material difference (Supplementary
# Data Table 4: R2 0.26 refitted vs 0.27 fixed; within-site r 0.60 both).
#
# Input   data/SHAPE_by_unit_with_SHMI.csv, Inputs/napeshm/treatments.csv
# Output  Outputs/Practice_additivity/labels_vs_SHMI.csv, labels_vs_SHMI_repeats.csv,
#         Fig_labels_vs_SHMI.png (Extended Data Fig. 6)
# Runs in under a minute.
################################################################################
suppressPackageStartupMessages({ library(dplyr); library(tidyr) })
OUT <- here::here("Outputs", "Practice_additivity")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
R_REP <- if (exists("LAB_REPEATS")) LAB_REPEATS else 100

d <- read.csv(here::here("data", "SHAPE_by_unit_with_SHMI.csv"))
trt <- read.csv(here::here("Inputs", "napeshm", "treatments.csv"), stringsAsFactors = FALSE) %>%
  mutate(MGT_combo = paste0("NAPESHM_", sitecode, "_", treatmentid),
         across(c(eu_tillage, eu_cover_crops, eu_crop_rot, eu_nutr_type), ~ tolower(trimws(.x)))) %>%
  select(MGT_combo, eu_tillage, eu_cover_crops, eu_crop_rot, eu_nutr_type)
a <- d %>% filter(system == "annual") %>% inner_join(trt, by = "MGT_combo") %>%
  transmute(sitecode, SH, SHMI,
            cover    = case_when(eu_cover_crops %in% c("grasses", "mix", "legume") ~ 1, eu_cover_crops == "none" ~ 0),
            organic  = case_when(eu_nutr_type %in% c("organic", "synthetic_and_organic") ~ 1, eu_nutr_type == "synthetic" ~ 0),
            reduced  = case_when(eu_tillage %in% c("none", "minimum") ~ 1, eu_tillage %in% c("intense", "moderate") ~ 0),
            rotation = case_when(eu_crop_rot %in% c("2", "3", "greater_than_3") ~ 1, eu_crop_rot == "1" ~ 0)) %>%
  filter(if_all(c(cover, organic, reduced, rotation), ~ !is.na(.x)))
P <- c("cover", "organic", "reduced", "rotation")
INT <- vapply(combn(P, 2, simplify = FALSE), function(pq) { a[[paste(pq, collapse = "_x_")]] <<- a[[pq[1]]] * a[[pq[2]]]; paste(pq, collapse = "_x_") }, "")
MODELS <- list("SHMI" = "SHMI", "Practice labels" = P, "Practice labels + interactions" = c(P, INT),
               "SHMI + labels" = c("SHMI", P))
cat(sprintf("%d annual systems at %d sites\n", nrow(a), n_distinct(a$sitecode)))

perf <- function(y, p, s) {
  sm <- tapply(y, s, mean); pm <- tapply(p, s, mean)
  yw <- y - ave(y, s); pw <- p - ave(p, s)
  ok <- vapply(split(seq_along(y), s), function(i) {
    if (length(i) < 3) return(NA_real_)
    if (sd(p[i]) == 0) return(0)
    as.numeric(cor(y[i], p[i], method = "spearman") > 0)
  }, numeric(1))
  c(R2 = 1 - sum((y - p)^2) / sum((y - mean(y))^2), between_r = cor(sm, pm), within_r = cor(yw, pw),
    sites_right = mean(ok, na.rm = TRUE), n_sites_ranked = sum(!is.na(ok)))
}
sites <- unique(a$sitecode)
reps <- bind_rows(lapply(seq_len(R_REP), function(r) {
  set.seed(2026 + r)
  fold <- setNames(sample(rep_len(1:10, length(sites))), sites)[a$sitecode]
  bind_rows(lapply(names(MODELS), function(m) {
    pr <- numeric(nrow(a))
    for (k in 1:10) {
      tr <- fold != k
      f <- lm(reformulate(MODELS[[m]], "SH"), data = a[tr, ])
      pr[!tr] <- predict(f, newdata = a[!tr, ])
    }
    tibble::as_tibble_row(c(repeat_id = r, perf(a$SH, pr, a$sitecode))) %>% mutate(model = m, .before = 1)
  }))
}))
write.csv(reps, file.path(OUT, "labels_vs_SHMI_repeats.csv"), row.names = FALSE)
q <- function(x, p) unname(quantile(x, p, na.rm = TRUE))
shmi <- reps %>% filter(model == "SHMI") %>% select(repeat_id, s_within = within_r, s_right = sites_right)
summ <- reps %>% left_join(shmi, by = "repeat_id") %>% group_by(model) %>%
  summarise(across(c(within_r, sites_right, between_r, R2), list(med = median, lo = ~ q(.x, 0.025), hi = ~ q(.x, 0.975))),
            SHMI_better_within = mean(s_within > within_r), SHMI_better_sites = mean(s_right > sites_right),
            n_sites_ranked = first(n_sites_ranked), .groups = "drop") %>%
  arrange(match(model, names(MODELS)))
write.csv(mutate(summ, across(where(is.double), ~ round(.x, 3))), file.path(OUT, "labels_vs_SHMI.csv"), row.names = FALSE)
cat("\nHeld-out sites (median and 95% range over repeats):\n")
print(as.data.frame(summ %>% transmute(model,
  within_r = sprintf("%.2f [%.2f, %.2f]", within_r_med, within_r_lo, within_r_hi),
  sites_right = sprintf("%.0f%% [%.0f, %.0f]", 100 * sites_right_med, 100 * sites_right_lo, 100 * sites_right_hi),
  between_r = sprintf("%.2f", between_r_med), R2 = sprintf("%.2f", R2_med),
  SHMI_better_within, SHMI_better_sites)), row.names = FALSE)
# ---- Extended Data Fig. 6: the checklist against the record, same measures as Fig. 2a ----
suppressPackageStartupMessages(library(ggplot2))
LAB <- c("SHMI" = "SHMI (management record)", "Practice labels" = "Four practice labels",
         "Practice labels + interactions" = "Labels + interactions", "SHMI + labels" = "SHMI + labels")
MET <- c(within_r = "r within sites", sites_right = "Sites ranked in the right direction (%)",
         between_r = "r between site means", R2 = "Overall R²")
fd <- bind_rows(lapply(names(MET), function(m) summ %>% transmute(model,
        metric = MET[[m]], med = .data[[paste0(m, "_med")]], lo = .data[[paste0(m, "_lo")]], hi = .data[[paste0(m, "_hi")]]))) %>%
  mutate(across(c(med, lo, hi), ~ ifelse(metric == MET[["sites_right"]], 100 * .x, .x)),
         label = factor(LAB[model], rev(LAB)), metric = factor(metric, MET),
         grp = ifelse(model == "SHMI", "SHMI", ifelse(grepl("^Practice", model), "Practice labels", "SHMI + labels")))
band <- filter(fd, model == "SHMI")
g <- ggplot(fd, aes(med, label, colour = grp)) +
  geom_rect(data = band, inherit.aes = FALSE, aes(xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf), fill = "#1b6ca8", alpha = 0.15) +
  geom_vline(data = band, aes(xintercept = med), colour = "#1b6ca8", linewidth = 0.3) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, orientation = "y") +
  geom_point(size = 2.4) +
  facet_wrap(~ metric, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = c(SHMI = "#1b6ca8", `Practice labels` = "#8a8a8a", `SHMI + labels` = "#e6a141"), guide = "none") +
  labs(x = "Median and 95% range across repeated fold assignments (sites withheld from fitting)", y = NULL) +
  theme_bw(base_size = 10) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"), panel.grid.minor = element_blank(),
        panel.spacing.x = unit(1.2, "lines"))
ggsave(file.path(OUT, "Fig_labels_vs_SHMI.png"), g, width = 11, height = 2.6, dpi = 300)
message("Labels vs SHMI done: ", OUT)
