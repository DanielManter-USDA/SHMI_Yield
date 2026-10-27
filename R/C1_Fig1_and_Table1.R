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
#   Fig. 1b  SHMI of example systems, scored with the SHMI package from example
#            management records (below).
#   Table 1  the example systems, their SHMI, and the expected corn yield change
#            relative to conventionally tilled corn-soybean at the average
#            within-site slope at 150 kg N/ha.
#
# Inputs  data/SHAPE_by_unit_with_SHMI.csv (sub-indices, SHMI and lv_SH by unit)
#         Outputs/Benchmark/benchmark_summary.csv (out-of-fold R2, optional)
# Outputs Outputs/Manuscript/Fig1.png, Fig1.pdf, Table1.csv, Fig1_path_coefficients.csv,
#         Fig1_example_systems.csv
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
OUT <- here::here("Outputs", "Manuscript")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
SLOPE_10   <- 220      # kg/ha per 10 SHMI points at 150 kg N/ha, average age (B2; Fig. 3a)
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
# 2. Example systems (Fig. 1b, Table 1)
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
  "AL", "+ no-till, rye and manure",                    "CS",      TRUE,    TRUE,  TRUE)
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
table1 <- scores %>%
  transmute(System = label, SHMI = round(SHMI), dSHMI = round(SHMI - ref),
            yield_kg_ha = ifelse(MGT_combo == "NT", NA, round(SLOPE_10 / 10 * (SHMI - ref))),
            bu_ac = round(yield_kg_ha / KG_PER_BU, 1),
            value_usd_ac = round(bu_ac * PRICE_BU))
names(table1) <- c("System", "SHMI", "ΔSHMI", "Yield change, kg ha-1", "bu ac-1", "Value, $ ac-1")
write.csv(table1, file.path(OUT, "Table1.csv"), row.names = FALSE, fileEncoding = "UTF-8")

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

ladder <- scores %>% mutate(label = factor(label, rev(SYSTEMS$label)),
                            col = c(CC = "#9e9e9e", CS = "#9e9e9e", NT = "#8a8a80", RY = "#1d9e75",
                                    MN = "#b8731a", AL = BLUE)[MGT_combo])
pb <- ggplot(ladder, aes(SHMI, label)) +
  geom_col(aes(fill = col), width = 0.62, colour = "#444444", linewidth = 0.2) +
  geom_vline(xintercept = ref, linetype = 2, linewidth = 0.3) +
  geom_text(aes(label = round(SHMI)), hjust = -0.3, fontface = "bold", size = 3.4) +
  scale_fill_identity() +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = c(0, 0)) +
  labs(x = "SHMI (0–100)", y = NULL) +
  theme_classic(base_size = 10)

fig <- (pa / pb) + plot_layout(heights = c(1.15, 0.85)) +
  plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = "bold", size = 16))
ggsave(file.path(OUT, "Fig1.png"), fig, width = 10, height = 8.2, dpi = 300, bg = "white")
ggsave(file.path(OUT, "Fig1.pdf"), fig, width = 10, height = 8.2, bg = "white")
cat("\nSaved Fig1.png/.pdf, Table1.csv, Fig1_path_coefficients.csv and Fig1_example_systems.csv in", OUT, "\n")
