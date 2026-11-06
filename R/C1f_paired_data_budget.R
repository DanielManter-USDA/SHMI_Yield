################################################################################
# C1f_paired_data_budget.R
#
# How much of the soil-health data can practice-by-practice (paired) comparisons use,
# and how does SHMI do on the same systems and on the ones pairs leave out?
# NAPESHM US systems (354; lv_SH = SH). Pairing rule as Fig. 1b (C1): two systems at
# the same site that differ in one practice and match on every other recorded label
# (cover crop, tillage class, rotation length, nutrient source, residue, grazing,
# irrigation, crop phase at sampling, crop category, other treatment).
#   strict pairs  every other label matches (Fig. 1b)
#   loose pairs   same site only; other labels free
# Practices: cover crop (grasses, mix, legume vs none); organic inputs (organic, or
# organic + synthetic, vs synthetic); reduced tillage (none, minimum vs moderate,
# intense); rotation (2+ crops vs 1).
#
# 1 Data budget: pairs, systems, sites per practice; the union of all strict pairs;
#   SHMI: every system, every within-site comparison.
# 2 SHMI vs soil health within sites on each subset: all systems, annual systems,
#   systems in each practice's pairs, in any pair, and in no pair.
# 3 What explains the change in soil health across pairs: practice labels (one effect
#   per practice) vs the change in SHMI (one slope), on the strict pairs and on all
#   within-site pairs of annual systems split by how many of the four labels differ
#   (0 = identical checklist). In-sample fits (labels get four parameters, SHMI one);
#   the out-of-fold comparison is C1e (Extended Data Fig. 6). Site-bootstrap interval
#   for the R2 difference on the strict pairs.
#
# 4 Soil health per SHMI point from each part (all systems, within and between sites):
#   are points interchangeable? (data_budget_parts.csv)
#
# Input   data/SHAPE_by_unit_with_SHMI.csv, Inputs/napeshm/treatments.csv
# Output  Outputs/Practice_additivity/data_budget.csv, data_budget_subsets.csv,
#         data_budget_explain.csv, data_budget_pairs.csv, Fig_data_budget.png/.pdf
#         (Fig_data_budget_panels.rds: both panels, drawn as Fig. 1d,e by C1)
# Runs in a few minutes (bootstrap resamples).
################################################################################
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(patchwork) })
OUT <- here::here("Outputs", "Practice_additivity")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
B_BOOT <- if (exists("BUDGET_BOOT")) BUDGET_BOOT else 2000
set.seed(2026)

d <- read.csv(here::here("data", "SHAPE_by_unit_with_SHMI.csv"))
ALL_LABELS <- c("eu_tillage", "eu_cover_crops", "eu_crop_rot", "eu_nutr_type", "eu_residue",
                "eu_grazing", "eu_irrigation", "crop_phase", "eu_crop_cat", "eu_other_trt")
trt <- read.csv(here::here("Inputs", "napeshm", "treatments.csv"), stringsAsFactors = FALSE) %>%
  mutate(MGT_combo = paste0("NAPESHM_", sitecode, "_", treatmentid),
         across(all_of(ALL_LABELS), ~ tolower(trimws(as.character(.x))))) %>%
  select(MGT_combo, all_of(ALL_LABELS))
allsys <- d %>% left_join(trt, by = "MGT_combo")
an <- allsys %>% filter(system == "annual")
N_SYS <- nrow(d); N_SITES <- n_distinct(d$sitecode)

PRACTICES <- tibble::tribble(
  ~practice,         ~label,
  "Cover crop",      "eu_cover_crops",
  "Organic inputs",  "eu_nutr_type",
  "Reduced tillage", "eu_tillage",
  "Rotation",        "eu_crop_rot")
group_of <- list(   # "A" = with the practice, "B" = without, NA = neither
  eu_cover_crops = function(x) ifelse(x %in% c("grasses", "mix", "legume"), "A", ifelse(x == "none", "B", NA)),
  eu_nutr_type   = function(x) ifelse(x %in% c("organic", "synthetic_and_organic"), "A", ifelse(x == "synthetic", "B", NA)),
  eu_tillage     = function(x) ifelse(x %in% c("none", "minimum"), "A", ifelse(x %in% c("intense", "moderate"), "B", NA)),
  eu_crop_rot    = function(x) ifelse(x %in% c("2", "3", "greater_than_3"), "A", ifelse(x == "1", "B", NA)))

make_pairs <- function(df, label, keys) {
  g <- df %>% mutate(grp = group_of[[label]](.data[[label]])) %>% filter(!is.na(grp))
  inner_join(filter(g, grp == "A"), filter(g, grp == "B"), by = c("sitecode", keys),
             suffix = c("_A", "_B"), relationship = "many-to-many") %>%
    transmute(sitecode, unit_A = MGT_combo_A, unit_B = MGT_combo_B,
              d_SHMI = SHMI_A - SHMI_B, d_SH = SH_A - SH_B)
}
strict <- lapply(setNames(PRACTICES$label, PRACTICES$practice),
                 function(l) make_pairs(an, l, setdiff(ALL_LABELS, l)))
loose  <- lapply(setNames(PRACTICES$label, PRACTICES$practice), function(l) make_pairs(an, l, character()))
units_of <- function(p) unique(c(p$unit_A, p$unit_B))

# ---- 1 data budget --------------------------------------------------------------------------
budget_row <- function(p, practice, pairing) {
  u <- units_of(p)
  tibble(practice, pairing, comparisons = nrow(p), systems = length(u), pct_systems = 100 * length(u) / N_SYS,
         sites = n_distinct(p$sitecode),
         median_dSH = median(p$d_SH), share_dSH_pos = mean(p$d_SH > 0),
         dSHMI_min = min(p$d_SHMI), dSHMI_max = max(p$d_SHMI),
         r_dSHMI_dSH = if (nrow(p) > 2) cor(p$d_SHMI, p$d_SH) else NA_real_)
}
S_ALL <- bind_rows(strict, .id = "practice")
in_pair <- units_of(S_ALL)
n_by_site <- count(allsys, sitecode)$n
budget <- bind_rows(
  bind_rows(Map(budget_row, strict, names(strict), "strict")),
  bind_rows(Map(budget_row, loose, names(loose), "loose")),
  budget_row(S_ALL, "Any of the four", "strict"),
  tibble(practice = "SHMI", pairing = "all within-site comparisons", comparisons = sum(choose(n_by_site, 2)),
         systems = N_SYS, pct_systems = 100, sites = N_SITES))
write.csv(mutate(budget, across(where(is.double), ~ round(.x, 3))), file.path(OUT, "data_budget.csv"), row.names = FALSE)
write.csv(S_ALL, file.path(OUT, "data_budget_pairs.csv"), row.names = FALSE)
cat("\n== Data budget ==\n"); print(as.data.frame(budget %>% select(practice, pairing, comparisons, systems, pct_systems, sites)), digits = 3)

# ---- 2 SHMI vs soil health within sites, by subset -----------------------------------------
within_r <- function(df) {
  df <- df %>% group_by(sitecode) %>% filter(n() >= 2) %>%
    mutate(x = SHMI - mean(SHMI), y = SH - mean(SH)) %>% ungroup()
  tibble(systems = nrow(df), sites = n_distinct(df$sitecode), within_r = cor(df$x, df$y))
}
subs <- c(list("All systems" = allsys, "Annual systems" = an,
               "Systems in any strict pair" = filter(allsys, MGT_combo %in% in_pair),
               "Systems in no strict pair" = filter(allsys, !MGT_combo %in% in_pair)),
          setNames(lapply(strict, function(p) filter(allsys, MGT_combo %in% units_of(p))),
                   paste("Systems in", names(strict), "pairs")))
subsets <- bind_rows(lapply(subs, within_r), .id = "subset")
write.csv(mutate(subsets, within_r = round(within_r, 3)), file.path(OUT, "data_budget_subsets.csv"), row.names = FALSE)
cat("\n== SHMI vs soil health within sites ==\n"); print(as.data.frame(subsets), digits = 3)

# ---- 3 labels vs SHMI as explanations of the change in soil health ---------------------------
# uncentred R2: a pair difference of zero is the null (no intercept)
r2_0 <- function(X, y) { b <- qr.coef(qr(X), y); b[is.na(b)] <- 0; e <- y - X %*% b; 1 - sum(e^2) / sum(y^2) }
lab_X <- function(p) model.matrix(~ 0 + practice, p)
explain <- function(y, XL, xS) c(R2_labels = r2_0(XL, y), R2_SHMI = r2_0(cbind(xS), y),
                                 r_dSHMI_dSH = if (length(y) > 2 && sd(xS) > 0) cor(xS, y) else NA)
e_strict <- explain(S_ALL$d_SH, lab_X(S_ALL), S_ALL$d_SHMI)
sites <- unique(S_ALL$sitecode)
boot <- replicate(B_BOOT, {
  s <- sample(sites, replace = TRUE)
  p <- bind_rows(lapply(s, function(k) S_ALL[S_ALL$sitecode == k, ]))
  r2_0(cbind(p$d_SHMI), p$d_SH) - r2_0(lab_X(p), p$d_SH)
})

# all within-site pairs of annual systems with the four labels defined, by how many differ
P4 <- c("cover", "organic", "reduced", "rotation")
a4 <- an %>% transmute(sitecode, MGT_combo, SHMI, SH,
                       cover    = as.integer(group_of$eu_cover_crops(eu_cover_crops) == "A"),
                       organic  = as.integer(group_of$eu_nutr_type(eu_nutr_type) == "A"),
                       reduced  = as.integer(group_of$eu_tillage(eu_tillage) == "A"),
                       rotation = as.integer(group_of$eu_crop_rot(eu_crop_rot) == "A")) %>%
  filter(if_all(all_of(P4), ~ !is.na(.x)))
ap <- inner_join(a4, a4, by = "sitecode", suffix = c("_A", "_B"), relationship = "many-to-many") %>%
  filter(MGT_combo_A < MGT_combo_B)
DL <- sapply(P4, function(k) ap[[paste0(k, "_A")]] - ap[[paste0(k, "_B")]])
ap <- ap %>% mutate(d_SH = SH_A - SH_B, d_SHMI = SHMI_A - SHMI_B, n_diff = rowSums(DL != 0))
by_ndiff <- bind_rows(lapply(list("0 (identical checklist)" = ap$n_diff == 0, "1" = ap$n_diff == 1,
                                  "2 or more" = ap$n_diff >= 2, "all" = rep(TRUE, nrow(ap))), function(s)
  tibble::as_tibble_row(c(pairs = sum(s), explain(ap$d_SH[s], DL[s, , drop = FALSE], ap$d_SHMI[s])))),
  .id = "labels_that_differ")
explain_tab <- bind_rows(
  tibble(set = "Strict pairs", labels_that_differ = "1", pairs = nrow(S_ALL), R2_labels = e_strict[["R2_labels"]],
         R2_SHMI = e_strict[["R2_SHMI"]], r_dSHMI_dSH = e_strict[["r_dSHMI_dSH"]],
         R2_SHMI_minus_labels_lo = quantile(boot, 0.025), R2_SHMI_minus_labels_hi = quantile(boot, 0.975)),
  mutate(by_ndiff, set = "All within-site pairs, annual systems"))
write.csv(mutate(explain_tab, across(where(is.double), ~ round(.x, 3))), file.path(OUT, "data_budget_explain.csv"), row.names = FALSE)
cat("\n== What explains the change in soil health (in-sample) ==\n"); print(as.data.frame(explain_tab), digits = 3)

# ---- 4 is a SHMI point worth the same soil health whichever part it comes from? -----------
# SHMI = 0.687 Cover + 0.150 Diversity + 0.163 InvDist + amendment bonus (to within
# rounding), so lv_SH is regressed on the four contributions; equal coefficients would
# mean points are interchangeable. Across all systems near-equality is largely built in
# (the weights were fitted to this target); within sites (deviations from site means)
# and between sites (site means) are the informative tests. 95% intervals from 1,000
# site bootstrap resamples.
W <- c(Cover = 0.687, Diversity = 0.150, InvDist = 0.163)
pc4 <- d %>% transmute(sitecode, SH, SHMI, p_cov = W[["Cover"]] * Cover, p_div = W[["Diversity"]] * Diversity,
                       p_inv = W[["InvDist"]] * InvDist, p_am = Amend_bonus)
PARTS4 <- c(p_cov = "Cover/Roots", p_div = "Diversity", p_inv = "Inverse Disturbance", p_am = "Amendments")
part_fit <- function(df, scale) {
  if (scale == "within") {
    df <- df %>% group_by(sitecode) %>% mutate(across(c(SH, SHMI, all_of(names(PARTS4))), ~ .x - mean(.x))) %>% ungroup()
    X <- as.matrix(df[names(PARTS4)]); one <- sum(df$SHMI * df$SH) / sum(df$SHMI^2)
    b <- qr.coef(qr(X), df$SH)
  } else {
    if (scale == "between") df <- df %>% group_by(sitecode) %>% summarise(across(c(SH, SHMI, all_of(names(PARTS4))), mean))
    b <- coef(lm(reformulate(names(PARTS4), "SH"), df))[names(PARTS4)]; one <- coef(lm(SH ~ SHMI, df))[["SHMI"]]
  }
  c(b, one_slope = one)
}
site_ids <- unique(pc4$sitecode)
parts_tab <- bind_rows(lapply(c("all", "within", "between"), function(sc) {
  est <- part_fit(pc4, sc)
  bt <- replicate(1000, {
    s <- sample(site_ids, replace = TRUE)
    part_fit(bind_rows(lapply(seq_along(s), function(i) mutate(filter(pc4, sitecode == s[i]), sitecode = paste0(s[i], "_", i)))), sc)
  })
  tibble(scale = sc, term = c(PARTS4, one_slope = "SHMI (one slope)")[names(est)], estimate = unname(est),
         lower = apply(bt, 1, quantile, 0.025, na.rm = TRUE), upper = apply(bt, 1, quantile, 0.975, na.rm = TRUE))
}))
write.csv(mutate(parts_tab, across(where(is.double), ~ round(.x, 3))), file.path(OUT, "data_budget_parts.csv"), row.names = FALSE)
cat("\n== Soil health (lv_SH points) per SHMI point from each part ==\n"); print(as.data.frame(parts_tab), digits = 2)

# ---- figure -----------------------------------------------------------------------------------
PCOL <- c("Cover crop" = "#4c9a2a", "Organic inputs" = "#8e6c3a", "Reduced tillage" = "#b5651d",
          "Rotation" = "#3a7dbf", "Any of the four" = "grey45", "SHMI" = "black")
fb <- budget %>% filter(pairing %in% c("strict", "all within-site comparisons")) %>%
  mutate(practice = factor(practice, levels = names(PCOL)),
         lab = sprintf("%d systems, %d sites", systems, sites))
pa <- ggplot(fb, aes(pct_systems, practice, fill = practice)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = lab), hjust = -0.05, size = 3) +
  scale_fill_manual(values = PCOL, guide = "none") +
  scale_y_discrete(limits = rev) +
  scale_x_continuous(limits = c(0, 150), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  labs(x = "Share of the 354 systems usable (%)", y = NULL,
       title = "Matched pairs use a minority of the data") +
  theme_bw(base_size = 10) + theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
                                   plot.title = element_text(face = "bold", size = 10.5))

fs <- subsets %>% filter(subset %in% c("All systems", "Systems in any strict pair", "Systems in no strict pair")) %>%
  mutate(subset = recode(subset, "Systems in any strict pair" = "Systems in a matched pair",
                         "Systems in no strict pair" = "Systems with no one-practice partner")) %>%
  transmute(what = sprintf("%s\n(%d systems, %d sites)", subset, systems, sites), r = within_r) %>%
  bind_rows(tibble(what = sprintf("Pairs with identical checklists\n(%d pairs)", sum(ap$n_diff == 0)),
                   r = by_ndiff$r_dSHMI_dSH[by_ndiff$labels_that_differ == "0 (identical checklist)"]))
pb <- ggplot(fs, aes(r, what)) +
  geom_segment(aes(x = 0, xend = r, yend = what), colour = "grey60") +
  geom_point(size = 3) + geom_text(aes(label = sprintf("%.2f", r)), hjust = -0.4, size = 3) +
  scale_y_discrete(limits = rev) + scale_x_continuous(limits = c(0, 0.8), expand = c(0, 0)) +
  labs(x = "Within-site correlation of SHMI with soil health (r)", y = NULL,
       title = "SHMI uses the data pairs leave out") +
  theme_bw(base_size = 10) + theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
                                   plot.title = element_text(face = "bold", size = 10.5))
saveRDS(list(pa, pb), file.path(OUT, "Fig_data_budget_panels.rds"))   # Fig. 1d,e (C1)
g <- pa + pb + plot_layout(widths = c(1.15, 1)) + plot_annotation(tag_levels = "a")
ggsave(file.path(OUT, "Fig_data_budget.png"), g, width = 10.5, height = 3.6, dpi = 300, bg = "white")
tryCatch(ggsave(file.path(OUT, "Fig_data_budget.pdf"), g, width = 10.5, height = 3.6, device = cairo_pdf),
         error = function(e) message("PDF not written: ", conditionMessage(e)))
cat("Written to", OUT, "\n")
