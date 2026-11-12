################################################################################
# C1_Fig1_and_Table1.R
#
# Fig. 1 and Table 1 (Moore & Manter, in review).
#   Fig. 1a  path diagram: the formula weights of SHMI 2.0.0 (blue), the paths
#            estimated by least squares on the 354 US NAPESHM systems (black):
#            Cover/Roots ~ Diversity + Inverse Disturbance, and soil health (lv_SH)
#            ~ SHMI, with 95% intervals from a site-level bootstrap; the total effect
#            of each sub-index on SHMI (weight + path to Cover/Roots x its weight)
#            and their sum (the overstatement from adding practices one at a time).
#   Fig. 1b  practices vary: pairs of annual systems at the same NAPESHM site that differ
#            in one practice and match on every other practice label (tillage class,
#            cover crop, rotation length, nutrient source, residue, grazing, irrigation,
#            crop phase at sampling, crop category, other treatment). For each practice:
#            the change in the part of SHMI it targets, in SHMI and in measured soil
#            health (lv_SH), one dot per pair, medians and interquartile ranges.
#   Fig. 1c  practices do not add up in measured soil health (no SHMI): each practice's own
#            within-site effect on lv_SH, their sum, and all four together (from
#            C1d_practice_additivity.R, run first; the all-four bar is split among the
#            practices in proportion to their own effects, only its total is estimated)
#   Fig. 1d  how much of the data matched pairs can use (C1f, Fig_data_budget_panels.rds).
#            Layout since 7 Nov 2026: a / b / c | d (panel e dropped; FIG1E <- TRUE restores it)
#   Table 1  the example systems (scored with the SHMI package from example records below), their SHMI, and the expected corn yield change
#            relative to conventionally tilled corn-soybean at the average
#            within-site slope at 150 kg N/ha.
#
# Inputs  data/SHAPE_by_unit_with_SHMI.csv (sub-indices, SHMI and lv_SH by unit)
#         Inputs/napeshm/treatments.csv (NAPESHM practice labels, for Fig. 1b)
#         Outputs/Benchmark/benchmark_summary.csv (out-of-fold R2, optional)
#         Outputs/Practice_additivity/additivity_effects.csv (C1d, for Fig. 1c)
# Outputs Outputs/Manuscript_work/Fig1.png, Fig1.pdf, ED_Fig_path.png (the path diagram, Extended Data since 7 Nov 2026), Table1.csv, Table1_next_step.csv, Fig1_path_coefficients.csv,
#         Fig1_example_systems.csv, Fig1b_practice_pairs.csv (one row per pair),
#         Fig1b_practice_summary.csv
#
# The example records are constructed here: four years (2019-2022), corn planted
# 1 May and harvested 1 October, soybean planted 15 May and harvested 5 October,
# in a temperate climate whose growing season runs April-October (the climate of the
# SHMI overview vignette). Conventional tillage: a heavy tandem disk 10 days before
# planting; every system has a planter pass. Cereal rye: planted 15 October after
# harvest, terminated 25 April. Manure: one organic amendment every spring.
################################################################################
suppressPackageStartupMessages({
  library(SHMI); library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
})
stopifnot(packageVersion("SHMI") >= "2.0.0")
OUT <- here::here("Outputs", "Manuscript_work")   # working copies; C3 copies the paper's items to Outputs/Manuscript/numbered
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
SLOPE_10   <- 220      # kg/ha per 10 SHMI points at 150 kg N/ha, average age (B2; Fig. 2a); read from B2 below
RM <- here::here("Outputs", "yield_main_default_total_agcros", "results_main.csv")
if (file.exists(RM)) {
  rm_ <- read.csv(RM, stringsAsFactors = FALSE)
  s10 <- rm_$effect[rm_$figure == "Fig 1" & grepl("^Within sites, at 150", rm_$estimate)][1]
  if (length(s10) && !is.na(s10)) SLOPE_10 <- s10
}
message("Table 1 uses ", round(SLOPE_10, 1), " kg/ha per 10 SHMI points")
KG_PER_BU  <- 53.0     # kg/ha dry matter per bu/ac (25.4 kg/bu at 15.5% moisture)
PRICE_BU   <- 4.50     # $/bu
B          <- 2000     # site-bootstrap draws
set.seed(2026)

# ==============================================================================
# 1. Path coefficients (Fig. 1a)
# ==============================================================================
d <- read.csv(here::here("data", "SHAPE_by_unit_with_SHMI.csv"))
w <- shmi_weights()
paths <- function(D) {
  cr <- coef(lm(Cover ~ Diversity + InvDist, data = D))
  sh <- coef(lm(SH ~ SHMI, data = D))
  p <- c(div_to_cover = cr[["Diversity"]], inv_to_cover = cr[["InvDist"]], shmi_to_sh = sh[["SHMI"]])
  c(p, total_cover = w[["Cover"]],
    total_diversity = w[["Diversity"]] + p[["div_to_cover"]] * w[["Cover"]],
    total_invdist   = w[["InvDist"]]   + p[["inv_to_cover"]] * w[["Cover"]],
    sum_of_totals   = w[["Cover"]] + w[["Diversity"]] + p[["div_to_cover"]] * w[["Cover"]] +
                      w[["InvDist"]] + p[["inv_to_cover"]] * w[["Cover"]])
}
est <- paths(d)
by_site <- split(d, d$sitecode)
boot <- t(replicate(B, paths(do.call(rbind, by_site[sample(names(by_site), replace = TRUE)]))))
coef_tab <- tibble(term = names(est), estimate = unname(est),
                   lo95 = apply(boot, 2, quantile, 0.025), hi95 = apply(boot, 2, quantile, 0.975))
r2_in <- summary(lm(SH ~ SHMI, data = d))$r.squared
bf <- here::here("Outputs", "Benchmark", "benchmark_summary.csv")
r2_oof <- if (file.exists(bf)) { b <- read.csv(bf); b$R2[b$model == "SHMI"] } else NA_real_
coef_tab <- bind_rows(coef_tab, tibble(term = c("R2_in_sample", "R2_out_of_fold"),
                                       estimate = c(r2_in, r2_oof), lo95 = NA, hi95 = NA))
write.csv(coef_tab, file.path(OUT, "Fig1_path_coefficients.csv"), row.names = FALSE)
cat("Path coefficients (Fig. 1a):\n"); print(as.data.frame(coef_tab), digits = 3, row.names = FALSE)

# ==============================================================================
# 2. Example systems (Table 1)
# ==============================================================================
YEARS   <- 2019:2022
CLIMATE <- list(tavg = c(-5, -3, 3, 10, 16, 21, 24, 23, 18, 11, 4, -2),
                prec = c(20, 25, 50, 70, 110, 115, 95, 85, 70, 50, 35, 25))
SYSTEMS <- tribble(
  ~id,  ~label,                                         ~rotation, ~notill, ~rye,  ~manure,
  "CC", "Conventional-till continuous corn",            "CC",      FALSE,   FALSE, FALSE,
  "CS", "Conventional-till corn–soybean (reference)", "CS",   FALSE,   FALSE, FALSE,
  "NT", "+ no-till",                                    "CS",      TRUE,    FALSE, FALSE,
  "RY", "+ cereal rye cover crop",                      "CS",      FALSE,   TRUE,  FALSE,
  "MN", "+ manure",                                     "CS",      FALSE,   FALSE, TRUE,
  "AL", "+ no-till, rye and manure",                    "CS",      TRUE,    TRUE,  TRUE,
  # the rest of the corn-soybean factorial, for the "what to add next" table (Table1_next_step.csv)
  "NR", "+ no-till and rye",                            "CS",      TRUE,    TRUE,  FALSE,
  "NM", "+ no-till and manure",                         "CS",      TRUE,    FALSE, TRUE,
  "RM", "+ rye and manure",                             "CS",      FALSE,   TRUE,  TRUE)
TABLE1_IDS <- c("CC", "CS", "NT", "RY", "MN", "AL")
dt <- function(y, md) as.Date(sprintf("%d-%s", y, md))

records <- function(s) {
  cash <- if (s$rotation == "CC") rep("Corn", length(YEARS)) else rep(c("Corn", "Soybean"), length.out = length(YEARS))
  plant <- ifelse(cash == "Corn", "05-01", "05-15"); harv <- ifelse(cash == "Corn", "10-01", "10-05")
  crop <- tibble(MGT_combo = s$id, CD_name = cash, CD_cat = "cash",
                 crop_start = dt(YEARS, plant), crop_end = dt(YEARS, harv))
  if (s$rye) {
    crop <- bind_rows(crop,
      tibble(MGT_combo = s$id, CD_name = "Rye", CD_cat = "cover",
             crop_start = c(dt(YEARS[1], "01-01"), dt(YEARS, "10-15")),
             crop_end   = c(dt(YEARS, "04-25"), dt(max(YEARS), "12-31"))))
  }
  dist <- tibble(MGT_combo = s$id, SD_date = dt(YEARS, plant), SD_equip = "PLANTER, REGULAR")
  if (!s$notill)
    dist <- bind_rows(dist, tibble(MGT_combo = s$id, SD_date = dt(YEARS, plant) - 10,
                                   SD_equip = "HARROW, DISK, TANDEM, HEAVYDUTY"))
  amend <- if (s$manure) tibble(MGT_combo = s$id, SA_date = dt(YEARS, "04-01"), SA_cat = "Organic") else
    tibble(MGT_combo = character(), SA_date = as.Date(character()), SA_cat = character())
  harvests <- crop %>% filter(CD_cat == "cash") %>% transmute(MGT_combo, CD_name, CD_cat, harv_date = crop_end)
  list(crop = crop, dist = arrange(dist, SD_date), amend = amend, harvests = harvests)
}
recs <- lapply(split(SYSTEMS, SYSTEMS$id), records)
cat_rows <- function(k) bind_rows(lapply(recs, `[[`, k))
crop <- cat_rows("crop"); dist <- cat_rows("dist"); amend <- cat_rows("amend"); harvests <- cat_rows("harvests")
rot_bounds <- tibble(MGT_combo = SYSTEMS$id, rot_start = dt(min(YEARS), "01-01"), rot_end = dt(max(YEARS), "12-31"),
                     rot_start_yr = min(YEARS), rot_end_yr = max(YEARS))
climate <- tibble(MGT_combo = SYSTEMS$id)
climate[sprintf("tavg_%02d", 1:12)] <- as.list(CLIMATE$tavg)
climate[sprintf("prec_%02d", 1:12)] <- as.list(CLIMATE$prec)
animal <- tibble(MGT_combo = character(), AD_start_date = as.Date(character()), AD_end_date = as.Date(character()))

par <- shmi_parameters()
scores <- SYSTEMS %>% select(MGT_combo = id, label) %>%
  left_join(compute_cover(select(crop, MGT_combo, CD_name, crop_start, crop_end), rot_bounds, climate), by = "MGT_combo") %>%
  left_join(compute_diversity(select(crop, MGT_combo, CD_name, crop_start, crop_end), rot_bounds), by = "MGT_combo") %>%
  left_join(compute_disturbance(dist, rot_bounds, select(crop, MGT_combo, CD_name, crop_start, crop_end), harvests),
            by = "MGT_combo") %>%
  left_join(compute_orginput_components(rot_bounds, amend, animal), by = "MGT_combo") %>%
  mutate(SHMI_base = w[["Cover"]] * Cover + w[["Diversity"]] * Diversity + w[["InvDist"]] * InvDist,
         Amend_bonus = par$amend_gap * p_amend * (100 - SHMI_base),
         SHMI = SHMI_base + Amend_bonus)
write.csv(scores, file.path(OUT, "Fig1_example_systems.csv"), row.names = FALSE)
cat("\nExample systems (check against Table 1: 46, 48, 56, 76, 79, 93):\n")
print(as.data.frame(scores %>% select(label, Cover, Diversity, InvDist, SHMI_base, SHMI) %>%
                      mutate(across(where(is.numeric), ~ round(.x, 1)))), row.names = FALSE)

ref <- scores$SHMI[scores$MGT_combo == "CS"]
table1 <- scores %>% filter(MGT_combo %in% TABLE1_IDS) %>%
  # every column from the unrounded score (rounding SHMI first would feed the rounded value
  # into the yield change: +2 kg/ha for the reference itself)
  transmute(System = label, SHMI_r = round(SHMI), dSHMI = round(SHMI - ref),
            yield_kg_ha = ifelse(MGT_combo == "NT", NA, round(SLOPE_10 / 10 * (SHMI - ref))),
            bu_ac = round(yield_kg_ha / KG_PER_BU, 1),
            value_usd_ac = round(bu_ac * PRICE_BU)) %>%
  rename(SHMI = SHMI_r)
names(table1) <- c("System", "SHMI", "ΔSHMI", "Yield change, kg ha-1", "bu ac-1", "Value, $ ac-1")
write.csv(table1, file.path(OUT, "Table1.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# What to add next: the SHMI gain from adding no-till, a rye cover crop or yearly manure to
# each corn-soybean system (the full factorial of the three practices). Practices that act
# on different parts (rye: Cover/Roots and Diversity; no-till: Inverse Disturbance) add;
# the organic-input bonus closes a share of the gap to 100, so manure adds least, and
# makes the others add less, where the score is already high.
fac <- scores %>% filter(MGT_combo %in% c("CS", "NT", "RY", "MN", "NR", "NM", "RM", "AL")) %>%
  left_join(SYSTEMS %>% select(MGT_combo = id, notill, rye, manure), by = "MGT_combo")
key <- function(nt, ry, mn) fac$SHMI[fac$notill == nt & fac$rye == ry & fac$manure == mn]
nxt <- fac %>% rowwise() %>%
  mutate(`+ no-till` = if (notill) NA_real_ else key(TRUE, rye, manure) - SHMI,
         `+ rye cover crop` = if (rye) NA_real_ else key(notill, TRUE, manure) - SHMI,
         `+ manure` = if (manure) NA_real_ else key(notill, rye, TRUE) - SHMI) %>%
  ungroup() %>%
  mutate(System = c(CS = "Tilled corn–soybean", NT = "No-till", RY = "Rye", MN = "Manure", NR = "No-till + rye",
                    NM = "No-till + manure", RM = "Rye + manure", AL = "All three")[MGT_combo]) %>%
  arrange(SHMI) %>%
  transmute(System, SHMI = round(SHMI), across(c(`+ no-till`, `+ rye cover crop`, `+ manure`), ~ round(.x)))
write.csv(nxt, file.path(OUT, "Table1_next_step.csv"), row.names = FALSE, na = "", fileEncoding = "UTF-8")
cat("\nWhat to add next (SHMI gain):\n"); print(as.data.frame(nxt), row.names = FALSE)

# ==============================================================================
# 3. Figure
# ==============================================================================
BLUE <- "#1a5ea8"; BLACK <- "#2b2b2b"; ORANGE <- "#b8641a"
f2 <- function(x) sprintf("%.2f", x)
E <- setNames(coef_tab$estimate, coef_tab$term)
boxes <- tribble(
  ~key,                  ~x,   ~y,   ~w,   ~h,  ~fill,
  "Diversity",           13,   40.8, 19,   5.6, "#cfcbf5",
  "Cover/Roots",         38.5, 40.8, 15,   5.6, "#9fe1cb",
  "Inverse Disturbance", 12.75,23.1, 19.5, 5.6, "#e5e4df",
  "SHMI_base",           55,   23.1, 14,   5.6, "#b5d4f2",
  "SHMI",                71.5, 23.1, 11,   5.6, "#8fbde6",
  "Soil health",         91,   23.1, 14,   5.0, "#f1efe8",
  "Amendments",          64,   12.3, 14,   5.6, "#f9c677")
tot <- c(Diversity = E[["total_diversity"]], `Cover/Roots` = E[["total_cover"]], `Inverse Disturbance` = E[["total_invdist"]])
boxes <- boxes %>% mutate(total = unname(tot[key]))
arrows <- tribble(
  ~x,   ~y,   ~xend, ~yend, ~col,  ~lw, ~lab,                          ~lx,  ~ly,
  22.5, 40.8, 31,    40.8,  BLACK, 1.6, f2(E[["div_to_cover"]]),       26.7, 43.4,
  22.5, 23.5, 34.5,  38,    BLACK, 0.7, f2(E[["inv_to_cover"]]),       26.8, 29.4,
  46,   38.0, 51.5,  25.9,  BLUE,  1.8, f2(w[["Cover"]]),              50.6, 33.6,
  22.5, 39.2, 48,    24.4,  BLUE,  0.7, f2(w[["Diversity"]]),          40.0, 31.4,
  22.5, 23.1, 48,    23.1,  BLUE,  0.7, f2(w[["InvDist"]]),            38.5, 24.6,
  62,   23.1, 66,    23.1,  BLUE,  1.1, "",                            NA,   NA,
  64,   15.1, 64,    20.3,  ORANGE,1.1, "",                            NA,   NA,
  77,   23.1, 84,    23.1,  BLACK, 1.6, f2(E[["shmi_to_sh"]]),         80.5, 25.8)
pa <- ggplot() +
  geom_segment(data = arrows, aes(x = x, y = y, xend = xend, yend = yend, colour = col, linewidth = lw),
               arrow = arrow(length = unit(2.6, "mm"), type = "closed"), lineend = "round") +
  geom_label(data = filter(arrows, nzchar(lab)), aes(x = lx, y = ly, label = lab, colour = col),
             fontface = "bold", size = 3.6, label.size = 0, fill = "white", label.padding = unit(0.08, "lines")) +
  geom_tile(data = boxes, aes(x = x, y = y, width = w, height = h, fill = fill), colour = "#444444", linewidth = 0.5) +
  geom_text(data = boxes, aes(x = x, y = y + ifelse(is.na(total), 0, 1.0), label = key),
            fontface = "bold", size = 3.6) +
  geom_text(data = filter(boxes, !is.na(total)), aes(x = x, y = y - 1.4, label = paste("total", f2(total))),
            colour = BLUE, size = 3.2) +
  annotate("text", x = 64, y = 6.5, label = sprintf("SHMI = (1 − %.1fp) × SHMI_base + %dp", par$amend_gap, round(100 * par$amend_gap)),
           colour = ORANGE, fontface = "bold", size = 3.4) +
  annotate("text", x = 91, y = 18.6, hjust = 0.5, vjust = 1, size = 3, colour = "#444444", lineheight = 1.3,
           label = sprintf("R² (in-sample) = %.2f%s", r2_in,
                           if (is.na(r2_oof)) "" else sprintf("\nR² (out-of-fold) = %.2f", r2_oof))) +
  scale_colour_identity() + scale_fill_identity() + scale_linewidth_identity() +
  coord_fixed(xlim = c(0, 100), ylim = c(3, 46), expand = FALSE) +
  theme_void()

# ---- Fig. 1b: practices vary (matched within-site pairs) ------------------------------
trt <- read.csv(here::here("Inputs", "napeshm", "treatments.csv"), stringsAsFactors = FALSE) %>%
  mutate(MGT_combo = paste0("NAPESHM_", sitecode, "_", treatmentid),
         across(c(eu_tillage, eu_cover_crops, eu_crop_rot, eu_nutr_type, eu_residue, eu_grazing,
                  eu_irrigation, crop_phase, eu_crop_cat, eu_other_trt), ~ tolower(trimws(.x)))) %>%
  select(-sitecode)
an <- d %>% filter(system == "annual") %>% inner_join(trt, by = "MGT_combo")
ALL_LABELS <- c("eu_tillage", "eu_cover_crops", "eu_crop_rot", "eu_nutr_type", "eu_residue",
                "eu_grazing", "eu_irrigation", "crop_phase", "eu_crop_cat", "eu_other_trt")
PRACTICES <- tribble(
  ~practice,                          ~label,          ~target,       ~target_lab,           ~col,
  "Cover crop vs none",               "eu_cover_crops", "Cover",       "Cover/Roots",         "#4c9a2a",
  "Organic vs synthetic inputs",      "eu_nutr_type",   "Amend_bonus", "organic-input bonus", "#8e6c3a",
  "Reduced vs conventional tillage",  "eu_tillage",     "InvDist",     "Inverse Disturbance", "#b5651d",
  "Rotation vs monoculture",          "eu_crop_rot",    "Diversity",   "Diversity",           "#3a7dbf")
group_of <- list(   # "A" = with the practice, "B" = without, NA = neither
  eu_cover_crops = function(x) ifelse(x %in% c("grasses", "mix", "legume"), "A", ifelse(x == "none", "B", NA)),
  eu_nutr_type   = function(x) ifelse(x %in% c("organic", "synthetic_and_organic"), "A", ifelse(x == "synthetic", "B", NA)),
  eu_tillage     = function(x) ifelse(x %in% c("none", "minimum"), "A", ifelse(x %in% c("intense", "moderate"), "B", NA)),
  eu_crop_rot    = function(x) ifelse(x %in% c("2", "3", "greater_than_3"), "A", ifelse(x == "1", "B", NA)))
VARS <- c("Cover", "Diversity", "InvDist", "Amend_bonus", "SHMI", "SH")
pairs <- bind_rows(lapply(seq_len(nrow(PRACTICES)), function(i) {
  pr <- PRACTICES[i, ]; keys <- setdiff(ALL_LABELS, pr$label)
  g <- an %>% mutate(grp = group_of[[pr$label]](.data[[pr$label]])) %>% filter(!is.na(grp))
  A <- filter(g, grp == "A"); B <- filter(g, grp == "B")
  inner_join(A, B, by = c("sitecode", keys), suffix = c("_A", "_B"), relationship = "many-to-many") %>%
    transmute(practice = pr$practice, sitecode, unit_A = MGT_combo_A, unit_B = MGT_combo_B,
              # 10 Nov 2026: the practice's own label in each system, and the background the two
              # systems share (every other label is identical within a pair, by construction)
              with_practice = .data[[paste0(pr$label, "_A")]], without_practice = .data[[paste0(pr$label, "_B")]],
              !!!setNames(lapply(keys, function(k) rlang::sym(k)), paste0("bg_", sub("^eu_", "", keys))),
              !!!setNames(lapply(VARS, function(v) rlang::expr(!!rlang::sym(paste0(v, "_A")) - !!rlang::sym(paste0(v, "_B")))),
                          paste0("d_", VARS)))
}))
# Shared background used to mark Fig. 1a: tillage for cover-crop, organic-input and rotation
# pairs (no-till or minimum vs moderate or intense); cover crop for tillage pairs.
pairs <- pairs %>% mutate(
  background = case_when(
    practice != "Reduced vs conventional tillage" & bg_tillage %in% c("none", "minimum") ~ "No-till or minimum tillage",
    practice != "Reduced vs conventional tillage" & bg_tillage %in% c("moderate", "intense") ~ "Moderate or intense tillage",
    practice == "Reduced vs conventional tillage" & bg_cover_crops %in% c("grasses", "mix", "legume") ~ "With a cover crop",
    practice == "Reduced vs conventional tillage" & bg_cover_crops == "none" ~ "Without a cover crop",
    TRUE ~ "Other or not recorded"),
  bg_fill = background %in% c("No-till or minimum tillage", "With a cover crop"))
# Does the practice's soil-health change depend on the shared background? Per practice:
# pairs, sites and median change in lv_SH in each background, and a Wilcoxon rank-sum test
# between the two main backgrounds (descriptive: pairs share systems and sites).
bg_test <- pairs %>% filter(background != "Other or not recorded") %>% group_by(practice) %>%
  group_modify(function(g, k) {
    s <- g %>% group_by(background) %>% summarise(pairs = n(), sites = n_distinct(sitecode),
                                                    median_dSH = median(d_SH), .groups = "drop")
    p <- if (n_distinct(g$background) == 2 && min(table(g$background)) >= 3)
      suppressWarnings(wilcox.test(d_SH ~ background, data = g)$p.value) else NA_real_
    mutate(s, p_wilcoxon = p)
  }) %>% ungroup()
write.csv(mutate(bg_test, across(where(is.double), ~ round(.x, 3))), file.path(OUT, "Fig1a_background_effects.csv"), row.names = FALSE)
cat("\nFig. 1a, change in soil health by shared background:\n"); print(as.data.frame(bg_test), digits = 3, row.names = FALSE)
# every pair with its shared background (Supplementary Table 6); written once d_target exists
write_supp6 <- function() write.csv(pairs %>% transmute(Practice = practice, Site = sitecode, `System with practice` = unit_A,
                              `System without` = unit_B, `Practice label (with / without)` = paste(with_practice, "/", without_practice),
                              `Shared: tillage` = bg_tillage, `Shared: cover crop` = bg_cover_crops,
                              `Shared: rotation (crops)` = bg_crop_rot, `Shared: nutrient source` = bg_nutr_type,
                              `Shared: residue` = bg_residue, `Shared: grazing` = bg_grazing, `Shared: irrigation` = bg_irrigation,
                              `Shared: crop phase at sampling` = bg_crop_phase, `Shared: crop category` = bg_crop_cat,
                              `Shared: other treatment` = bg_other_trt,
                              `Change in targeted part` = round(d_target, 1), `Change in SHMI` = round(d_SHMI, 1),
                              `Change in soil health (lv_SH)` = round(d_SH, 1)),
          file.path(OUT, "SuppTable6_practice_pairs.csv"), row.names = FALSE, na = "")
q <- function(x, p) unname(quantile(x, p))
# change in the part each practice targets
pairs <- pairs %>% left_join(select(PRACTICES, practice, target), by = "practice") %>%
  mutate(d_target = case_when(target == "Cover" ~ d_Cover, target == "Diversity" ~ d_Diversity,
                              target == "InvDist" ~ d_InvDist, target == "Amend_bonus" ~ d_Amend_bonus)) %>%
  select(-target)
summ <- pairs %>% group_by(practice) %>%
  summarise(pairs = n(), sites = n_distinct(sitecode),
            target_med = median(d_target), target_min = min(d_target), target_max = max(d_target),
            shmi_med = median(d_SHMI), shmi_min = min(d_SHMI), shmi_max = max(d_SHMI),
            sh_med = median(d_SH), sh_share_pos = mean(d_SH > 0), .groups = "drop")
write.csv(pairs, file.path(OUT, "Fig1b_practice_pairs.csv"), row.names = FALSE)
write_supp6()
write.csv(mutate(summ, across(where(is.double), ~ round(.x, 2))), file.path(OUT, "Fig1b_practice_summary.csv"), row.names = FALSE)
cat("\nFig. 1b, matched practice pairs:\n"); print(as.data.frame(summ), digits = 3, row.names = FALSE)

COLS <- c("Change in the part it targets", "Change in SHMI", "Change in measured soil health")
# 10 Nov 2026: the middle column (change in SHMI) is no longer drawn; FIG1A_SHMI <- TRUE restores it
if (!exists("FIG1A_SHMI")) FIG1A_SHMI <- FALSE
COLS_SHOW <- if (FIG1A_SHMI) COLS else COLS[c(1, 3)]
long <- pairs %>% left_join(PRACTICES, by = "practice") %>%
  transmute(practice, col, target_lab, bg_fill, `Change in the part it targets` = d_target,
            `Change in SHMI` = d_SHMI, `Change in measured soil health` = d_SH) %>%
  pivot_longer(all_of(COLS), names_to = "panel", values_to = "value") %>%
  filter(panel %in% COLS_SHOW) %>%
  mutate(panel = factor(panel, COLS_SHOW),
         dot = ifelse(panel == COLS[1], col, ifelse(panel == COLS[2], "#222222", "#8e44ad")))
ylab <- summ %>% mutate(lab = sprintf("%s\n%d pairs, %d sites", practice, pairs, sites))
long <- long %>% mutate(practice = factor(practice, rev(PRACTICES$practice), rev(ylab$lab[match(PRACTICES$practice, ylab$practice)])))
med <- long %>% group_by(practice, panel) %>% summarise(m = median(value), lo = q(value, 0.25), hi = q(value, 0.75), .groups = "drop")
tlab <- long %>% filter(panel == COLS[1]) %>% distinct(practice, panel, target_lab, col)
pb <- ggplot(long, aes(value, practice)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3, colour = "grey50") +
  # filled: the pair shares no-till or minimum tillage (tillage pairs: shares a cover crop);
  # open: shares moderate or intense tillage (tillage pairs: no cover crop)
  geom_point(aes(colour = dot, fill = ifelse(bg_fill, dot, "white")), shape = 21, stroke = 0.6,
             position = position_jitter(height = 0.2, width = 0, seed = 4), size = 1.5, alpha = 0.8) +
  scale_fill_identity() +
  geom_segment(data = med, aes(x = lo, xend = hi, y = practice, yend = practice), linewidth = 1.6) +
  geom_point(data = med, aes(m, practice), shape = 21, fill = "white", size = 2.6, stroke = 0.8) +
  geom_text(data = tlab, aes(x = Inf, y = practice, label = target_lab, colour = col), hjust = 1.05, vjust = -1.1,
            size = 2.6, fontface = "italic") +
  facet_wrap(~ panel, nrow = 1, scales = "free_x") +
  scale_colour_identity() +
  labs(x = "Change between paired systems (points)", y = NULL) +
  theme_bw(base_size = 10) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())

# ---- Fig. 1c: practices do not add up in measured soil health (C1d) ------------------
ADD <- here::here("Outputs", "Practice_additivity", "additivity_effects.csv")
if (!file.exists(ADD)) stop("Fig. 1c needs ", ADD, ": run C1d_practice_additivity.R first")
ae <- read.csv(ADD, stringsAsFactors = FALSE); ge <- function(t) ae$estimate[ae$term == t]
P4  <- c("cover", "organic", "reduced", "rotation")
C4  <- setNames(PRACTICES$col, P4)                       # same colours as the practices in b
XL  <- c(cover = "Cover\ncrop", organic = "Organic\ninputs", reduced = "Reduced\ntillage", rotation = "Rotation",
         added = "Added\nup", together = "All four\ntogether")
m4 <- sapply(P4, ge); tog <- ge("together")
ia <- ae[ae$term == "interactions", ]
cd <- bind_rows(tibble(x = P4, practice = P4, value = m4),
                tibble(x = "added", practice = P4, value = m4),
                tibble(x = "together", practice = P4, value = m4 * tog / sum(m4))) %>%
  mutate(x = factor(XL[x], XL), practice = factor(practice, rev(P4)))
ctot <- tibble(x = factor(XL[c("added", "together")], XL), value = c(sum(m4), tog))
pc <- ggplot(cd, aes(x, value, fill = practice)) +
  geom_col(width = 0.68, colour = "white", linewidth = 0.3) +
  geom_text(data = ctot, aes(x, value, label = sprintf("+%.0f", value)), inherit.aes = FALSE,
            vjust = -0.4, size = 3.6, fontface = "bold") +
  annotate("text", x = 0.6, y = max(ctot$value) * 1.12, hjust = 0, vjust = 1, size = 3, colour = "grey30", lineheight = 1.1,
           label = sprintf("interactions %.0f\n(95%% CI %.0f to %.0f)", ia$estimate, ia$lower, ia$upper)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  scale_fill_manual(values = C4, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Change in measured soil health (lv_SH)") +
  theme_classic(base_size = 10)

# a how SHMI scores the system (path diagram) across the top, so SHMI's parts are introduced
# before b (practices vary: pairs) and c (they do not add up in measured soil health)
# 7 Nov 2026: Fig. 1e (C1f's SHMI on the systems pairs leave out) dropped; layout a / b / c | d.
# FIG1E <- TRUE before sourcing restores the earlier layout (a; b | c; d | e).
if (!exists("FIG1E")) FIG1E <- FALSE
DB <- here::here("Outputs", "Practice_additivity", "Fig_data_budget_panels.rds")
# 7 Nov 2026 (FIG1_MERGED, default TRUE): to fit Nature's page budget, Figs 1 and 2 merge into
# one soil-health figure: a practices vary (old 1b), b they do not add up (old 1c), c matched
# pairs use a minority (old 1d), d SHMI against its parts at held-out sites (old Fig. 2a, from
# A4b's Fig2a_panel.rds; run A4b first). The path diagram (old 1a) is saved alone as
# ED_Fig_path.png for Extended Data. FIG1_MERGED <- FALSE restores the layout below.
if (!exists("FIG1_MERGED")) FIG1_MERGED <- TRUE
F2A <- file.path(OUT, "Fig2a_panel.rds")
ggsave(file.path(OUT, "ED_Fig_path.png"), pa, width = 12.5, height = 4.6, dpi = 300, bg = "white")
# 10 Nov 2026 (FIG1_V3, default TRUE): a practices vary (target part and soil health only),
# b they do not add up, c SHMI against its parts at held-out sites (two within-site
# measures). The data-budget panel is dropped (its numbers are in the text).
if (!exists("FIG1_V3")) FIG1_V3 <- TRUE
# 10 Nov 2026 (FIG_V4, default TRUE): two soil-health figures.
#   Fig. 1  a how SHMI scores the system (path diagram); b SHMI against measured soil health:
#           all systems, between and within sites (C1h)
#   Fig. 2  a what a yes/no practice label and SHMI each explain within sites (C1g);
#           b how much of the data matched pairs can use (C1f); c SHMI against its parts at
#           held-out sites (A4b)
# "Practices do not add up" (old Fig. 1b) is saved alone for Extended Data (ED_Fig_additivity.png);
# the matched-pair panel (old Fig. 1a) is no longer drawn (its pairs are Supplementary Table 6).
if (!exists("FIG_V4")) FIG_V4 <- TRUE
PA2 <- here::here("Outputs", "Practice_additivity")
F1B <- file.path(PA2, "Fig_SHMI_vs_SH_panel.rds"); F2L <- file.path(PA2, "Fig_label_vs_SHMI_panel.rds")
ggsave(file.path(OUT, "ED_Fig_additivity.png"), pc, width = 6, height = 4.2, dpi = 300, bg = "white")
if (FIG_V4 && all(file.exists(c(F1B, F2L, DB, F2A)))) {
  f1b <- readRDS(F1B); f2l <- readRDS(F2L); db <- readRDS(DB); p2a <- readRDS(F2A)
  row_b <- ((db[[1]] + labs(title = NULL)) | wrap_elements(full = p2a)) + plot_layout(widths = c(1, 1.2))
  fig2 <- (wrap_elements(full = f2l) / row_b) + plot_layout(heights = c(1, 1)) +
    plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = "bold", size = 16))
  ggsave(file.path(OUT, "Fig2.png"), fig2, width = 12.5, height = 9, dpi = 300, bg = "white")
  fig <- (pa / wrap_elements(full = f1b)) + plot_layout(heights = c(1, 0.8)); FIG1_H <- 8.4
} else if (FIG_V4) {
  stop("Fig. 1/2 need ", paste(c(F1B, F2L, DB, F2A)[!file.exists(c(F1B, F2L, DB, F2A))], collapse = ", "),
       ": run C1h, C1g, C1f and A4b first (or set FIG_V4 <- FALSE for the old layout)")
} else if (FIG1_V3 && file.exists(F2A)) {
  p2a <- readRDS(F2A)
  fig <- pb / ((wrap_elements(full = pc) | wrap_elements(full = p2a)) + plot_layout(widths = c(1, 1.35))) +
    plot_layout(heights = c(1, 0.95)); FIG1_H <- 8.6
} else if (FIG1_MERGED && file.exists(DB) && file.exists(F2A)) {
  db <- readRDS(DB); p2a <- readRDS(F2A)
  row2 <- (wrap_elements(full = pc) | (db[[1]] + labs(title = NULL))) + plot_layout(widths = c(1, 1.4))
  fig <- (pb / row2 / wrap_elements(full = p2a)) + plot_layout(heights = c(0.9, 0.8, 0.85)); FIG1_H <- 11.5
} else if (file.exists(DB) && !FIG1E) {
  if (FIG1_MERGED) message("Fig. 1 drawn unmerged: ", F2A, " not found (run A4b_Fig2_soil_health.R first)")
  db <- readRDS(DB)
  # b (the pairs) across the full width; c (practices do not add up) beside d (data budget)
  # c wrapped so its axis is not aligned to b's long row labels; d without its working title
  row3 <- (wrap_elements(full = pc) | (db[[1]] + labs(title = NULL))) + plot_layout(widths = c(1, 1.4))
  fig <- (pa / pb / row3) + plot_layout(heights = c(0.95, 0.85, 0.75)); FIG1_H <- 12
} else if (file.exists(DB)) {
  db <- readRDS(DB)
  bottom <- (pb | pc) + plot_layout(widths = c(1.75, 1))
  row3 <- (db[[1]] | db[[2]]) + plot_layout(widths = c(1.15, 1))
  fig <- (pa / bottom / row3) + plot_layout(heights = c(0.95, 1, 0.7)); FIG1_H <- 11.6
} else {
  message("Fig. 1 drawn without d: ", DB, " not found (run C1f_paired_data_budget.R)")
  bottom <- (pb | pc) + plot_layout(widths = c(1.75, 1))
  fig <- (pa / bottom) + plot_layout(heights = c(0.95, 1)); FIG1_H <- 8.6
}
fig <- fig + plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = "bold", size = 16))
ggsave(file.path(OUT, "Fig1.png"), fig, width = 12.5, height = FIG1_H, dpi = 300, bg = "white")
tryCatch(ggsave(file.path(OUT, "Fig1.pdf"), fig, width = 12.5, height = FIG1_H, bg = "white", device = cairo_pdf),
         error = function(e) message("Fig1.pdf not written (cairo_pdf unavailable): ", conditionMessage(e)))
cat("\nSaved Fig1.png/.pdf, Table1.csv, Fig1_path_coefficients.csv and Fig1_example_systems.csv in", OUT, "\n")
