################################################################################
# B2f_SHMI_yield_paired_budget.R
#
# The yield counterpart of C1f: in the ten long-term experiments, how much of the corn
# yield record can practice-by-practice (paired) comparisons use, and what does SHMI
# give on the same treatments and on the ones pairs leave out?
#
# Practice labels per treatment, from each experiment's SHMI workbook over the analysed
# period (first experiment year to the last corn year):
#   cover     any cover crop planted (CD_cat "Cover")
#   tillage   no-till (no pass other than the sprayer placeholder), conventional
#             (moldboard or disk plough) or reduced (any other tillage)
#   rotation  the set of non-cover crops grown (e.g. corn; corn+soybean; corn+oats+alfalfa)
#   manure    any manure or compost (biochar and fertilizer excluded)
# A treatment takes the most common label among its plots and phases.
# (labels in section 2 are signed differences: cover and manure 1/0, tillage rank with
# less tillage positive, number of crops in the rotation)
# Pairs: two treatments at the same experiment that differ in exactly one of the four.
#   strict   also the same "other" factor (residue removal, irrigation, organic, planter)
#            and recorded N within 25 kg/ha (fertilizer + manure N, legume credit not
#            counted, so rotation pairs at the same fertilizer rate still match)
#   one practice   any N and other factor
# Paired yield difference: within each site-year where both treatments have corn, the
# difference of their mean yields; averaged over those site-years.
#
# 1 Data budget: treatments, plot-years and experiments per practice and in the union of
#   strict pairs; SHMI: every treatment and plot-year (all within-site comparisons).
# 2 Across strict pairs and across all N- and other-matched pairs (by how many labels
#   differ, 0 = identical checklist): does the change in SHMI track the yield difference?
#   Practice labels (one effect per practice) vs the change in SHMI (one slope), in-sample.
# 3 The model of Fig. 3a refitted on the plot-years of treatments in any strict pair and on
#   those in no strict pair, against all plot-years.
#
# Input   Outputs/yield_main_default_total_agcros/B2_data.rds (B2), the SHMI workbooks
#         (Inputs/yield_workbooks and the AgCROS folder recorded in B2_data.rds)
# Output  Outputs/yield_paired_budget/yield_budget.csv, yield_pairs.csv, yield_explain.csv,
#         yield_subset_slopes.csv, treatment_labels.csv, Fig_yield_budget.png/.pdf
# Runs in several minutes (three model fits).
################################################################################
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(nlme); library(ggplot2); library(patchwork) })
setwd(here::here())
IN  <- "Outputs/yield_main_default_total_agcros/B2_data.rds"
OUT <- "Outputs/yield_paired_budget"
if (!file.exists(IN)) stop(IN, " not found: run B2_SHMI_yield_main.R first")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
ctrl  <- lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200)
N_TOL <- 25
D <- readRDS(IN); corn <- D$corn
corn <- corn %>% group_by(site) %>% mutate(SHMI_w = SHMI - mean(SHMI[!duplicated(MGT_combo)])) %>% ungroup()

# ---- practice labels from the workbooks ------------------------------------------------
wb_file <- function(site) {
  f <- c(file.path("Inputs", "yield_workbooks", paste0("SHMI_", site, ".xlsx")),
         if (!is.null(D$AGCROS_DIR)) file.path(D$AGCROS_DIR, paste0("SHMI_", site, ".xlsx")))
  f <- f[file.exists(f)]
  if (!length(f)) stop("No workbook found for ", site)
  f[1]
}
yr <- function(x) as.integer(format(as.Date(x), "%Y"))
labels_site <- function(site) {
  f <- wb_file(site)
  # an empty tab reads every column as logical; make the text columns text
  rd <- function(sh) readxl::read_excel(f, sheet = sh, skip = 3) %>%
    mutate(across(any_of(c("MGT_combo", "CD_cat", "CD_name", "SD_equip", "SA_cat", "SA_source")), as.character))
  combos <- unique(corn$MGT_combo[corn$site == site])
  y0 <- D$site_start[[site]]; y1 <- max(corn$year[corn$site == site])
  inside <- function(x) { y <- yr(x); !is.na(y) & y >= y0 & y <= y1 }
  cd <- rd("Crop_Diversity") %>% filter(MGT_combo %in% combos, inside(CD_plant_date))
  sd <- rd("Soil_Disturbance") %>% filter(MGT_combo %in% combos, inside(SD_date), !is.na(SD_equip),
                                          !grepl("^SPRAYER", SD_equip))
  sa <- rd("Soil_Amendments") %>% filter(MGT_combo %in% combos, inside(SA_date))
  manure_rx <- "manure|compost|dairy|cpm|poultry|swine|cattle|slurry"
  sa <- sa %>% filter(grepl(manure_rx, paste(SA_cat, SA_source), ignore.case = TRUE),
                      !grepl("biochar|charcoal|fertili", paste(SA_cat, SA_source), ignore.case = TRUE))
  tibble(MGT_combo = combos) %>%
    left_join(cd %>% group_by(MGT_combo) %>%
                summarise(cover = any(CD_cat %in% "Cover"),
                          rotation = paste(sort(unique(tolower(trimws(CD_name[!CD_cat %in% "Cover"])))), collapse = "+"),
                          .groups = "drop"), by = "MGT_combo") %>%
    # tillage class from the years of the record: no-till if fewer than 20% of years had
    # any tillage; otherwise conventional if a moldboard or disk plow was used in at least
    # half of the tilled years, else reduced. (Until 9 Nov 2026 one plow pass anywhere made
    # the whole record conventional: Hoytville's no-till plots, plowed once in 1983, were
    # labelled conventional, so pairs that differ in tillage counted as rotation-only pairs.)
    left_join(sd %>% mutate(year = yr(SD_date)) %>% group_by(MGT_combo, year) %>%
                summarise(conv = any(grepl("MOLDBOARD|PLOW, DISK", SD_equip)), .groups = "drop") %>%
                group_by(MGT_combo) %>%
                summarise(tillage = ifelse(n() < 0.2 * (y1 - y0 + 1), "no-till",
                                           ifelse(mean(conv) >= 0.5, "conventional", "reduced")),
                          .groups = "drop"), by = "MGT_combo") %>%
    left_join(sa %>% distinct(MGT_combo) %>% mutate(manure = TRUE), by = "MGT_combo") %>%
    mutate(site = site, cover = coalesce(cover, FALSE), tillage = coalesce(tillage, "no-till"),
           manure = coalesce(manure, FALSE), rotation = coalesce(na_if(rotation, ""), "none"))
}
lab_combo <- bind_rows(lapply(unique(corn$site), labels_site))
mode1 <- function(x) names(sort(table(x), decreasing = TRUE))[1]

# (corn already has columns such as rotation and tillage from the treatment codes;
# drop them here so the workbook labels are the ones used)
TRT <- corn %>% select(-any_of(c("cover", "tillage", "rotation", "manure"))) %>%
  left_join(select(lab_combo, -site), by = "MGT_combo") %>%
  group_by(site, trt) %>%
  summarise(SHMI = mean(SHMI[!duplicated(MGT_combo)]), plot_years = n(),
            N_rec = median(N_kg_rec, na.rm = TRUE), other = as.character(first(other)),
            cover = mode1(cover), tillage = mode1(tillage), rotation = mode1(rotation), manure = mode1(manure),
            .groups = "drop")
write.csv(TRT, file.path(OUT, "treatment_labels.csv"), row.names = FALSE)

# ---- pairs and paired yield differences ------------------------------------------------
PR <- c(cover = "Cover crop", tillage = "Tillage", rotation = "Rotation", manure = "Manure")
sy <- corn %>% group_by(site, trt, site_year) %>% summarise(y = mean(yield_kgha), .groups = "drop")
P <- inner_join(TRT, TRT, by = "site", suffix = c("_A", "_B"), relationship = "many-to-many") %>%
  filter(as.character(trt_A) < as.character(trt_B))
DIFF <- sapply(names(PR), function(k) P[[paste0(k, "_A")]] != P[[paste0(k, "_B")]])
P <- P %>% mutate(n_diff = rowSums(DIFF),
                  practice = ifelse(n_diff == 1, PR[names(PR)][max.col(DIFF * 1, ties.method = "first")], NA),
                  matched = other_A == other_B & abs(N_rec_A - N_rec_B) <= N_TOL,
                  strict = n_diff == 1 & matched,
                  d_SHMI = SHMI_A - SHMI_B)
dy <- inner_join(sy, sy, by = c("site", "site_year"), suffix = c("_A", "_B"), relationship = "many-to-many") %>%
  semi_join(P, by = c("site", "trt_A", "trt_B")) %>%
  group_by(site, trt_A, trt_B) %>% summarise(d_yield = mean(y_A - y_B), site_years = n(), .groups = "drop")
P <- P %>% left_join(dy, by = c("site", "trt_A", "trt_B"))
# signed label differences (A minus B): cover and manure 1/0, less tillage positive,
# more crops in the rotation positive
till_rank <- c("no-till" = 1, reduced = 2, conventional = 3)
n_crops <- function(r) ifelse(r == "none", 0, lengths(strsplit(r, "+", fixed = TRUE)))
SL <- with(P, cbind(cover = (cover_A == "TRUE") - (cover_B == "TRUE"),
                    tillage = unname(till_rank[tillage_B] - till_rank[tillage_A]),
                    rotation = n_crops(rotation_A) - n_crops(rotation_B),
                    manure = (manure_A == "TRUE") - (manure_B == "TRUE")))
# orient each single-practice pair so that A has the practice (cover, manure, less
# tillage, more crops); the labels flip with it
flip <- P$n_diff == 1 & rowSums(SL) < 0
P <- P %>% mutate(d_SHMI = ifelse(flip, -d_SHMI, d_SHMI), d_yield = ifelse(flip, -d_yield, d_yield))
SL[flip, ] <- -SL[flip, ]
P <- bind_cols(P, as_tibble(SL, .name_repair = ~ paste0("s_", .x)))
write.csv(select(P, site, trt_A, trt_B, n_diff, practice, matched, strict, d_SHMI, d_yield, site_years),
          file.path(OUT, "yield_pairs.csv"), row.names = FALSE)

# ---- 1 data budget ------------------------------------------------------------------------
PY_ALL <- nrow(corn); TR_ALL <- nrow(TRT); EXP_ALL <- n_distinct(TRT$site)
budget_row <- function(p, practice, pairing) {
  u <- TRT %>% filter(trt %in% c(as.character(p$trt_A), as.character(p$trt_B)))
  tibble(practice, pairing, pairs = nrow(p), treatments = nrow(u), plot_years = sum(u$plot_years),
         pct_plot_years = 100 * sum(u$plot_years) / PY_ALL, experiments = n_distinct(u$site),
         median_d_yield = if (nrow(p)) median(p$d_yield, na.rm = TRUE) else NA,
         share_d_yield_pos = if (nrow(p)) mean(p$d_yield > 0, na.rm = TRUE) else NA,
         r_dSHMI_dyield = if (sum(!is.na(p$d_yield)) > 2) cor(p$d_SHMI, p$d_yield, use = "complete.obs") else NA)
}
SP <- filter(P, strict); OP <- filter(P, n_diff == 1)
budget <- bind_rows(
  bind_rows(lapply(PR, function(k) budget_row(filter(SP, practice == k), k, "strict"))),
  bind_rows(lapply(PR, function(k) budget_row(filter(OP, practice == k), k, "one practice, any N"))),
  budget_row(SP, "Any of the four", "strict"),
  tibble(practice = "SHMI", pairing = "all within-site comparisons", pairs = nrow(P), treatments = TR_ALL,
         plot_years = PY_ALL, pct_plot_years = 100, experiments = EXP_ALL))
by_exp <- TRT %>% mutate(in_strict = trt %in% c(as.character(SP$trt_A), as.character(SP$trt_B))) %>%
  group_by(site) %>% summarise(treatments = n(), in_strict_pairs = sum(in_strict),
                               plot_years_in_pairs = sum(plot_years[in_strict]), plot_years = sum(plot_years),
                               .groups = "drop")
write.csv(mutate(budget, across(where(is.double), ~ round(.x, 3))), file.path(OUT, "yield_budget.csv"), row.names = FALSE)
write.csv(by_exp, file.path(OUT, "yield_budget_by_experiment.csv"), row.names = FALSE)
cat("\n== Yield data budget ==\n"); print(as.data.frame(select(budget, practice, pairing, pairs, treatments, plot_years, pct_plot_years, experiments)), digits = 3)
cat("\n== By experiment ==\n"); print(as.data.frame(by_exp))

# ---- 2 does the change in SHMI track the yield difference? ------------------------------------
r2_0 <- function(X, y) { b <- qr.coef(qr(X), y); b[is.na(b)] <- 0; e <- y - X %*% b; 1 - sum(e^2) / sum(y^2) }
explain <- function(p) {
  p <- filter(p, !is.na(d_yield))
  XL <- as.matrix(p[, paste0("s_", names(PR))])
  tibble(pairs = nrow(p), experiments = n_distinct(p$site),
         R2_labels = if (nrow(p) > 4) r2_0(XL, p$d_yield) else NA,
         R2_SHMI = if (nrow(p) > 2) r2_0(cbind(p$d_SHMI), p$d_yield) else NA,
         r_dSHMI_dyield = if (nrow(p) > 2 && sd(p$d_SHMI) > 0) cor(p$d_SHMI, p$d_yield) else NA,
         kg_per_10_SHMI = if (nrow(p) > 2) 10 * unname(coef(lm(d_yield ~ 0 + d_SHMI, p))) else NA)
}
MP <- filter(P, matched)
explain_tab <- bind_rows(
  mutate(explain(SP), set = "Strict single-practice pairs"),
  mutate(explain(filter(MP, n_diff == 0)), set = "Matched pairs, identical checklist"),
  mutate(explain(filter(MP, n_diff == 1)), set = "Matched pairs, 1 practice differs"),
  mutate(explain(filter(MP, n_diff >= 2)), set = "Matched pairs, 2+ practices differ"),
  mutate(explain(MP), set = "All matched pairs"))
write.csv(mutate(explain_tab, across(where(is.double), ~ round(.x, 3))), file.path(OUT, "yield_explain.csv"), row.names = FALSE)
cat("\n== What explains the paired yield difference (in-sample) ==\n"); print(as.data.frame(explain_tab), digits = 3)

# ---- 3 the Fig. 3a model on subsets ------------------------------------------------------------
varies  <- function(d, v) any(tapply(d[[v]], d$site_year, function(x) diff(range(x))) > 1e-8)
n_terms <- function(d) { v <- grep("^N[12]_", names(d), value = TRUE); v[vapply(v, function(k) varies(d, k), logical(1))] }
o_term  <- function(d) if (n_distinct(d$other) > 1) "other" else character()
# Subsets can make fixed effects aliased (e.g. an "other" level or an N term confounded
# with site-year once some treatments are gone), which lme cannot drop itself: the
# fixed-effect matrix is built here, aliased columns are dropped (QR pivoting), and
# lme is fitted on the rest.
fit_slope <- function(d, label) {
  d <- droplevels(d)
  f <- reformulate(c("site_year", n_terms(d), o_term(d), "SHMI_w",
                     if (varies(d, "N_c")) "SHMI_w:N_c", "SHMI_w:age_c"), "yield_kgha")
  X <- model.matrix(f, d)
  q <- qr(X); keep <- sort(q$pivot[seq_len(q$rank)])
  dropped <- setdiff(colnames(X), colnames(X)[keep])
  if (length(dropped)) message(label, ": dropped aliased terms: ", paste(head(dropped, 8), collapse = ", "))
  X <- X[, keep, drop = FALSE]
  if (!"SHMI_w" %in% colnames(X)) return(tibble(subset = label, note = "SHMI_w aliased"))
  d$X <- X
  m <- tryCatch(lme(yield_kgha ~ 0 + X, random = ~ 1 | trt/plot_id,
                    correlation = corCAR1(form = ~ year | trt/plot_id),
                    data = d, method = "REML", control = ctrl),
                error = function(e) { message(label, ": ", conditionMessage(e)); conditionMessage(e) })
  if (is.character(m)) return(tibble(subset = label, note = m))
  tt <- summary(m)$tTable["XSHMI_w", ]; qq <- qt(0.975, tt[["DF"]])
  tibble(subset = label, experiments = n_distinct(d$site), treatments = n_distinct(d$trt), plot_years = nrow(d),
         slope_per_10 = 10 * tt[["Value"]], lo = 10 * (tt[["Value"]] - qq * tt[["Std.Error"]]),
         hi = 10 * (tt[["Value"]] + qq * tt[["Std.Error"]]), p = tt[["p-value"]], note = "")
}
in_strict <- unique(c(as.character(SP$trt_A), as.character(SP$trt_B)))
slopes <- bind_rows(fit_slope(corn, "All plot-years (Fig. 3a)"),
                    fit_slope(filter(corn, as.character(trt) %in% in_strict), "Treatments in any strict pair"),
                    fit_slope(filter(corn, !as.character(trt) %in% in_strict), "Treatments in no strict pair"))
write.csv(mutate(slopes, across(where(is.double), ~ signif(.x, 4))), file.path(OUT, "yield_subset_slopes.csv"), row.names = FALSE)
# (the full-data row is refitted the same way; it should reproduce Fig. 3a, 220 kg/ha per 10 SHMI)
cat("\n== Fig. 3a model on subsets ==\n"); print(as.data.frame(slopes), digits = 3)

# ---- figure ------------------------------------------------------------------------------------
PCOL <- c("Cover crop" = "#4c9a2a", "Manure" = "#8e6c3a", "Tillage" = "#b5651d", "Rotation" = "#3a7dbf",
          "Any of the four" = "grey45", "SHMI" = "black")
fb <- budget %>% filter(pairing %in% c("strict", "all within-site comparisons")) %>%
  mutate(practice = factor(practice, levels = names(PCOL)),
         lab = sprintf("%d treatments, %s plot-years, %d exp.", treatments, format(plot_years, big.mark = ","), experiments))
pa <- ggplot(fb, aes(pct_plot_years, practice, fill = practice)) +
  geom_col(width = 0.7) + geom_text(aes(label = lab), hjust = -0.03, size = 2.8) +
  scale_fill_manual(values = PCOL, guide = "none") + scale_y_discrete(limits = rev) +
  scale_x_continuous(limits = c(0, 200), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  labs(x = "Share of corn plot-years usable (%)", y = NULL, title = "Share of the yield record matched pairs can use") +
  theme_bw(base_size = 10) + theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
                                   plot.title = element_text(face = "bold", size = 10.5))
pb <- ggplot(filter(slopes, !is.na(slope_per_10)), aes(slope_per_10, factor(subset, levels = rev(subset)))) +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.15, orientation = "y") + geom_point(size = 2.6) +
  geom_text(aes(label = sprintf("%s plot-years", format(plot_years, big.mark = ","))), vjust = -1.1, size = 2.8) +
  labs(x = expression("Within-site SHMI slope (kg ha"^-1*" per 10 SHMI)"), y = NULL,
       title = "The Fig. 3a model on each subset") +
  theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 10.5))
g <- pa + pb + plot_layout(widths = c(1.2, 1)) + plot_annotation(tag_levels = "a")
saveRDS(pa, file.path(OUT, "Fig_yield_budget_panel.rds"))
ggsave(file.path(OUT, "Fig_yield_budget.png"), g, width = 11, height = 3.6, dpi = 300, bg = "white")
tryCatch(ggsave(file.path(OUT, "Fig_yield_budget.pdf"), g, width = 11, height = 3.6, device = cairo_pdf),
         error = function(e) message("PDF not written: ", conditionMessage(e)))
cat("Written to", OUT, "\n")
