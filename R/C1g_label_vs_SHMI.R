################################################################################
# C1g_label_vs_SHMI.R   (10 Nov 2026; Fig. 2a and Extended Data Fig. 3)
#
# Is SHMI better than a yes/no practice label? For each practice of Fig. 1 (cover crop,
# organic inputs, reduced tillage, rotation), the annual NAPESHM systems at the sites
# that have systems both with and without the practice. Within sites (site means
# removed), measured soil health (lv_SH) is predicted by
#   label   the practice yes/no
#   SHMI    the score (published weights)
#   both    label + SHMI
# and the explained within-site variance is split into what the two share, what only
# SHMI explains and what only the label explains (shared = label + SHMI - both).
#
# Two versions:
#   in-sample   fitted on all sites; nested F/likelihood-ratio tests with site-clustered
#               standard errors (SHMI adds to the label? the label adds to SHMI?), and
#               the within-label test: SHMI centred on systems with the same label at the
#               same site (does SHMI sort systems the label calls the same?)
#   out-of-fold leave-sites-out cross-validation (the practice's sites in 10 folds, 100
#               repeats, seed 2026 + repeat): each model's slopes fitted on the training sites, the held-out
#               sites predicted, and the within-site R2 of the held-out predictions
#               (medians and 95% ranges over repeats). Fig. 2a uses these.
# REFIT_IN_FOLD = TRUE (default): in every fold SHMI's Cover/Roots and Inverse
# Disturbance weights are refitted on the NAPESHM systems at the training sites only
# (all 74 sites minus the held-out ones), so the held-out sites are truly unseen:
# SHMI_base = wC*Cover + 0.15*Diversity + (0.85 - wC)*InvDist, the amendment bonus as
# published (Diversity weight, growing-season share and amendment/livestock shares stay
# fixed, as in the calibration), wC chosen on a 0.40-0.85 grid to maximise the
# site-weighted R2 of lv_SH on SHMI. FALSE uses the published weights throughout.
#
# Input   data/SHAPE_by_unit_with_SHMI.csv, Inputs/napeshm/treatments.csv
# Output  Outputs/Practice_additivity/label_vs_SHMI.csv (both versions, the tests),
#         label_vs_SHMI_repeats.csv, label_vs_SHMI_systems.csv,
#         Fig_label_vs_SHMI.png and Fig_label_vs_SHMI_panel.rds (Fig. 2a),
#         Fig_label_within.png (Extended Data Fig. 3)
# Runs in under a minute.
################################################################################
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })
OUT <- here::here("Outputs", "Practice_additivity")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
R_REP <- if (exists("LAB_REPEATS")) LAB_REPEATS else 100
REFIT_IN_FOLD <- if (exists("LAB_REFIT_IN_FOLD")) LAB_REFIT_IN_FOLD else TRUE

d <- read.csv(here::here("data", "SHAPE_by_unit_with_SHMI.csv"))
trt <- read.csv(here::here("Inputs", "napeshm", "treatments.csv"), stringsAsFactors = FALSE) %>%
  mutate(MGT_combo = paste0("NAPESHM_", sitecode, "_", treatmentid),
         across(c(eu_tillage, eu_cover_crops, eu_crop_rot, eu_nutr_type), ~ tolower(trimws(.x)))) %>%
  select(MGT_combo, eu_tillage, eu_cover_crops, eu_crop_rot, eu_nutr_type)
an <- d %>% filter(system == "annual") %>% inner_join(trt, by = "MGT_combo")
PR <- tibble::tribble(
  ~practice,         ~col,      ~yes,                   ~no,
  "Cover crop",      "#4c9a2a", "cover crop",           "no cover crop",
  "Organic inputs",  "#8e6c3a", "organic inputs",       "synthetic only",
  "Reduced tillage", "#b5651d", "no-till or minimum",   "moderate or intense",
  "Rotation",        "#3a7dbf", "two or more crops",    "one crop")
lab_of <- list(
  "Cover crop"      = function(z) case_when(z$eu_cover_crops %in% c("grasses", "mix", "legume") ~ 1, z$eu_cover_crops == "none" ~ 0),
  "Organic inputs"  = function(z) case_when(z$eu_nutr_type %in% c("organic", "synthetic_and_organic") ~ 1, z$eu_nutr_type == "synthetic" ~ 0),
  "Reduced tillage" = function(z) case_when(z$eu_tillage %in% c("none", "minimum") ~ 1, z$eu_tillage %in% c("intense", "moderate") ~ 0),
  "Rotation"        = function(z) case_when(z$eu_crop_rot %in% c("2", "3", "greater_than_3") ~ 1, z$eu_crop_rot == "1" ~ 0))

# systems for one practice: sites with both levels
sys_for <- function(p) {
  an %>% mutate(lab = lab_of[[p]](an)) %>% filter(!is.na(lab)) %>%
    group_by(sitecode) %>% filter(n_distinct(lab) == 2) %>% ungroup() %>%
    transmute(practice = p, sitecode, MGT_combo, lab, SHMI, SH, Cover, Diversity, InvDist, p_amend)
}
S <- bind_rows(lapply(PR$practice, sys_for))
write.csv(S, file.path(OUT, "label_vs_SHMI_systems.csv"), row.names = FALSE)

dm <- function(x, g) x - ave(x, g)
# within-site least squares with site-clustered sandwich SEs
fit_w <- function(z, cols) {
  y <- dm(z$SH, z$sitecode)
  if (!length(cols)) return(list(rss = sum(y^2), k = 0))
  X <- sapply(cols, function(cn) dm(z[[cn]], z$sitecode)); X <- matrix(X, ncol = length(cols), dimnames = list(NULL, cols))
  b <- solve(crossprod(X), crossprod(X, y)); r <- as.numeric(y - X %*% b)
  XtXi <- solve(crossprod(X))
  meat <- Reduce(`+`, lapply(split(seq_along(y), z$sitecode), function(i) { u <- crossprod(X[i, , drop = FALSE], r[i]); u %*% t(u) }))
  V <- XtXi %*% meat %*% XtXi
  list(b = setNames(as.numeric(b), cols), se = setNames(sqrt(diag(V)), cols), rss = sum(r^2), k = length(cols), V = V)
}
in_sample <- function(z) {
  ns <- n_distinct(z$sitecode); n <- nrow(z)
  z$SHMI_inlab <- z$SHMI - ave(z$SHMI, interaction(z$sitecode, z$lab))
  f0 <- fit_w(z, character()); fl <- fit_w(z, "lab"); fs <- fit_w(z, "SHMI"); fb <- fit_w(z, c("lab", "SHMI"))
  fw <- fit_w(z, c("lab", "SHMI_inlab"))
  r2 <- function(f) 1 - f$rss / f0$rss
  pcl <- function(f, cn) 2 * pt(-abs(f$b[[cn]] / f$se[[cn]]), ns - 1)   # cluster-robust t, sites - 1 df
  lr <- function(a, b) pchisq(n * log(a$rss / b$rss), b$k - a$k, lower.tail = FALSE)
  gap <- z %>% group_by(sitecode, lab) %>% summarise(m = mean(SHMI), .groups = "drop") %>%
    pivot_wider(names_from = lab, values_from = m) %>% summarise(g = mean(`1` - `0`)) %>% pull(g)
  tibble(systems = n, sites = ns, with_practice = sum(z$lab == 1),
         label_effect = fl$b[["lab"]], label_P = pcl(fl, "lab"),
         SHMI_per10 = 10 * fs$b[["SHMI"]], SHMI_P = pcl(fs, "SHMI"),
         R2_label_in = r2(fl), R2_SHMI_in = r2(fs), R2_both_in = r2(fb),
         P_SHMI_adds_to_label = lr(fl, fb), P_label_adds_to_SHMI = lr(fs, fb),
         P_SHMI_adds_clustered = pcl(fb, "SHMI"), P_label_adds_clustered = pcl(fb, "lab"),
         SHMI_within_label_per10 = 10 * fw$b[["SHMI_inlab"]], SHMI_within_label_P = pcl(fw, "SHMI_inlab"),
         SHMI_gap = gap, SHMI_equivalent_effect = fs$b[["SHMI"]] * gap)
}
# out of fold: slopes from training sites, within-site R2 of the held-out predictions
oof <- function(z, rep) {
  sites <- unique(z$sitecode); set.seed(2026 + rep)
  fold <- setNames(sample(rep_len(1:10, length(sites))), sites)
  z$f <- fold[z$sitecode]
  pred <- matrix(NA_real_, nrow(z), 3, dimnames = list(NULL, c("label", "SHMI", "both")))
  for (k in unique(z$f)) {
    if (REFIT_IN_FOLD) {
      wC <- fit_wC(d[!d$sitecode %in% unique(z$sitecode[z$f == k]), ])
      WC_LOG <<- c(WC_LOG, wC); z$SHMI <- shmi_w(z, wC)
    }
    tr <- z[z$f != k, ]; te <- which(z$f == k)
    for (m in colnames(pred)) {
      cols <- switch(m, label = "lab", SHMI = "SHMI", both = c("lab", "SHMI"))
      b <- fit_w(tr, cols)$b
      X <- sapply(cols, function(cn) dm(z[[cn]][te], z$sitecode[te])); X <- matrix(X, ncol = length(cols))
      pred[te, m] <- as.numeric(X %*% b)
    }
  }
  y <- dm(z$SH, z$sitecode)
  tibble(repeat_id = rep, model = colnames(pred), R2 = apply(pred, 2, function(p) 1 - sum((y - p)^2) / sum(y^2)))
}
# SHMI from its parts with Cover/Roots weight wC (Inverse Disturbance 0.85 - wC)
shmi_w <- function(z, wC) { b <- wC * z$Cover + 0.15 * z$Diversity + (0.85 - wC) * z$InvDist; b + 0.60 * z$p_amend * (100 - b) }
chk <- cor(shmi_w(d, 0.687), d$SHMI, use = "complete.obs")
message(sprintf("SHMI rebuilt from its parts: r = %.5f with the published score", chk))
if (chk < 0.999) warning("SHMI rebuilt from its parts does not match the published score; check the columns")
WC_GRID <- seq(0.40, 0.85, by = 0.01)
fit_wC <- function(tr) {
  tr <- tr[stats::complete.cases(tr[, c("SH", "Cover", "Diversity", "InvDist", "p_amend")]), ]
  w <- 1 / ave(rep(1, nrow(tr)), tr$sitecode, FUN = length)
  r2 <- sapply(WC_GRID, function(g) { x <- shmi_w(tr, g); stats::cov.wt(cbind(x, tr$SH), wt = w, cor = TRUE)$cor[1, 2]^2 })
  WC_GRID[which.max(r2)]
}
WC_LOG <- c()
res <- list(); reps <- list()
for (p in PR$practice) {
  z <- filter(S, practice == p)
  rr <- bind_rows(lapply(seq_len(R_REP), function(r) oof(z, r))) %>% mutate(practice = p)
  reps[[p]] <- rr
  o <- rr %>% group_by(model) %>% summarise(med = median(R2), lo = quantile(R2, 0.025), hi = quantile(R2, 0.975), .groups = "drop")
  g <- function(m, v) o[[v]][o$model == m]
  res[[p]] <- in_sample(z) %>% mutate(practice = p, .before = 1,
    R2_label_oof = g("label", "med"), R2_label_oof_lo = g("label", "lo"), R2_label_oof_hi = g("label", "hi"),
    R2_SHMI_oof = g("SHMI", "med"), R2_SHMI_oof_lo = g("SHMI", "lo"), R2_SHMI_oof_hi = g("SHMI", "hi"),
    R2_both_oof = g("both", "med"),
    SHMI_beats_label_share = mean(rr$R2[rr$model == "SHMI"] > rr$R2[rr$model == "label"]))
}
res <- bind_rows(res) %>%
  mutate(shared_oof = R2_label_oof + R2_SHMI_oof - R2_both_oof,
         SHMI_only_oof = R2_both_oof - R2_label_oof, label_only_oof = R2_both_oof - R2_SHMI_oof)
if (REFIT_IN_FOLD) {
  message(sprintf("Cover/Roots weight refitted in each fold: median %.2f (range %.2f-%.2f; published 0.687)",
                  median(WC_LOG), min(WC_LOG), max(WC_LOG)))
  res$refit_in_fold <- TRUE; res$wC_fold_min <- min(WC_LOG); res$wC_fold_max <- max(WC_LOG)
}
write.csv(mutate(res, across(where(is.double), ~ signif(.x, 4))), file.path(OUT, "label_vs_SHMI.csv"), row.names = FALSE)
write.csv(bind_rows(reps), file.path(OUT, "label_vs_SHMI_repeats.csv"), row.names = FALSE)
cat("\nPractice label vs SHMI within sites (R2 out of fold, medians over repeats; tests in sample):\n")
print(as.data.frame(res %>% transmute(practice, systems, sites,
  label = round(R2_label_oof, 3), SHMI = round(R2_SHMI_oof, 3), shared = round(shared_oof, 3),
  SHMI_only = round(SHMI_only_oof, 3), label_only = round(label_only_oof, 3),
  P_SHMI_adds = signif(P_SHMI_adds_clustered, 2), P_label_adds = signif(P_label_adds_clustered, 2),
  SHMI_within_label_per10 = round(SHMI_within_label_per10, 1), P_within = signif(SHMI_within_label_P, 2))), row.names = FALSE)

# ---- Fig. 2a: what each explains (out of fold) -------------------------------------------
bars <- res %>% select(practice, shared_oof, SHMI_only_oof, label_only_oof) %>%
  pivot_longer(-practice, names_to = "part", values_to = "v") %>%
  mutate(v = 100 * pmax(v, 0),
         part = factor(recode(part, shared_oof = "Shared by label and SHMI", SHMI_only_oof = "SHMI only",
                              label_only_oof = "Label only"),
                       c("Label only", "SHMI only", "Shared by label and SHMI")),
         practice = factor(practice, rev(PR$practice)))
txt <- res %>% transmute(practice = factor(practice, rev(PR$practice)),
                         x = 100 * pmax(R2_both_oof, R2_label_oof, R2_SHMI_oof),
                         lab = sprintf("label %.0f%%   SHMI %.0f%%", 100 * R2_label_oof, 100 * R2_SHMI_oof))
ylab <- res %>% transmute(practice, l = sprintf("%s\n%d systems, %d sites", practice, systems, sites))
fig2a <- ggplot(bars, aes(v, practice, fill = part)) +
  geom_col(width = 0.68) +
  geom_text(data = txt, aes(x, practice, label = lab), inherit.aes = FALSE, hjust = -0.05, size = 3) +
  scale_fill_manual(values = c("Shared by label and SHMI" = "#bdbdbd", "SHMI only" = "#1f5f99", "Label only" = "#e69f00"),
                    breaks = c("Shared by label and SHMI", "SHMI only", "Label only"), name = NULL) +
  scale_y_discrete(labels = setNames(ylab$l, ylab$practice)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(x = "Within-site variance in measured soil health explained at held-out sites (%)", y = NULL) +
  theme_classic(base_size = 10) +
  theme(legend.position = "bottom", legend.key.size = unit(0.35, "cm"))
saveRDS(fig2a, file.path(OUT, "Fig_label_vs_SHMI_panel.rds"))
ggsave(file.path(OUT, "Fig_label_vs_SHMI.png"), fig2a, width = 6.5, height = 3.6, dpi = 300, bg = "white")

# ---- Extended Data Fig. 3: SHMI among systems the label calls the same ----------------------
w <- S %>% group_by(practice, sitecode, lab) %>%
  mutate(x = SHMI - mean(SHMI), y = SH - mean(SH)) %>% ungroup() %>%
  left_join(select(PR, practice, col), by = "practice") %>%
  left_join(res %>% transmute(practice, tag = sprintf("%+.1f per 10 SHMI points, P %s", SHMI_within_label_per10,
                                                       format.pval(SHMI_within_label_P, digits = 2, eps = 0.001))), by = "practice") %>%
  mutate(practice = factor(practice, PR$practice))
edw <- ggplot(w, aes(x, y)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) + geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_point(aes(colour = col, shape = factor(lab)), size = 1.4, alpha = 0.7) +
  geom_smooth(aes(colour = col), method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.9) +
  geom_text(data = distinct(w, practice, tag), aes(x = -Inf, y = Inf, label = tag), hjust = -0.05, vjust = 1.4, size = 2.8) +
  facet_wrap(~ practice, nrow = 1) +
  scale_colour_identity() + scale_shape_manual(values = c(`0` = 1, `1` = 16), labels = c("without", "with"), name = "Practice") +
  labs(x = "SHMI relative to systems with the same label at the same site",
       y = "Measured soil health (lv_SH),\nsame reference") +
  theme_bw(base_size = 10) + theme(strip.background = element_blank(), strip.text = element_text(face = "bold"),
                                   panel.grid.minor = element_blank(), legend.position = "bottom")
ggsave(file.path(OUT, "Fig_label_within.png"), edw, width = 11, height = 3.6, dpi = 300, bg = "white")
message("Wrote ", file.path(OUT, "label_vs_SHMI.csv"), " and the Fig. 2a / Extended Data Fig. 3 panels")
