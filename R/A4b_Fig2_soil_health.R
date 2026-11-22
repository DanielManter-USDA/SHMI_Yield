################################################################################
# A4b_Fig2_soil_health.R
#
# Manuscript Fig. 2: the whole system predicts soil health better than any part.
#   a  out-of-fold performance of SHMI and each part alone (unconstrained models
#      optional, FIG2_UNCONSTRAINED; from A4a_SHMI_benchmark.R): within sites first
#      (correlation, sites ranked in the right direction), then between site means
#      and overall R2; shaded band = SHMI's 95% range across repeated fold assignments
#   b  out-of-fold predicted vs measured soil health, all 354 systems
#   c  the same for the 74 site means
#   d  the same within sites (deviations from site means)
#
# Inputs
#   FIN/us_oof_predictions.csv          (A3a_SHMI_calibration_plus.R)
#   Outputs/Benchmark/benchmark_summary.csv (A4a_SHMI_benchmark.R)
# Output
#   Outputs/Manuscript_work/Fig2.png (and .pdf); Fig2a_panel.rds (panel a, drawn by C1 as
#   Fig. 1d since 7 Nov 2026) and ED_Fig_validation.png (b-d as an Extended Data figure)
#
# Run after A3a_SHMI_calibration_plus.R and A4a_SHMI_benchmark.R.
################################################################################
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(here)
})
if (!requireNamespace("patchwork", quietly = TRUE)) stop("Install patchwork: install.packages('patchwork')")

FIN   <- if (exists("FIN_DIR")) FIN_DIR else here::here("Outputs", "Final_plus")
BENCH <- here::here("Outputs", "Benchmark", "benchmark_summary.csv")
OUT   <- here::here("Outputs", "Manuscript_work")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# a-c: out-of-fold predictions
# ------------------------------------------------------------------------------
oof <- read.csv(file.path(FIN, "us_oof_predictions.csv"))
if (!"pred_oof" %in% names(oof)) stop("us_oof_predictions.csv has no pred_oof column")
d <- oof %>% transmute(sitecode, SH, pred = pred_oof)

tag_theme <- theme_classic(base_size = 11) +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0, margin = margin(b = 4)),
        plot.title.position = "plot")

panel <- function(x, y, tag, xlab, ylab) {
  r <- cor(x, y); slope <- unname(coef(lm(y ~ x))[2])
  ggplot(tibble(x = x, y = y), aes(x, y)) +
    geom_point(colour = "#1f6fab", alpha = 0.6, size = 1.4, stroke = 0) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "#d1495b", linewidth = 1) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.08, vjust = 1.15, size = 3.4, lineheight = 0.95,
             label = sprintf("r = %.2f\nR² = %.2f\nslope = %.2f", r, r^2, slope)) +
    labs(title = tag, x = xlab, y = ylab) + tag_theme
}
sm <- d %>% group_by(sitecode) %>% summarise(x = mean(pred), y = mean(SH), .groups = "drop")
wi <- d %>% group_by(sitecode) %>% mutate(x = pred - mean(pred), y = SH - mean(SH)) %>% ungroup()
pa <- panel(d$pred, d$SH, "b", "Predicted soil health", "Measured soil health")
pb <- panel(sm$x, sm$y, "c", "Site mean predicted soil health", "Site mean measured soil health")
pc <- panel(wi$x, wi$y, "d", "Predicted, minus site mean", "Measured, minus site mean")
cat(sprintf("b: n = %d systems | c: %d site means\n", nrow(d), nrow(sm)))

# ------------------------------------------------------------------------------
# a: SHMI against its parts and unconstrained models
# ------------------------------------------------------------------------------
bs <- read.csv(BENCH, check.names = FALSE)

ROWS <- tibble::tribble(                       # benchmark name -> label, group (top to bottom)
  ~model,                         ~label,                          ~group,
  "SHMI",                         "SHMI",                          "SHMI",
  "Cover/Roots",                  "Cover/Roots alone",             "Single part",
  "Diversity",                    "Diversity alone",               "Single part",
  "InvDist",                      "Inverse Disturbance alone",     "Single part")
# 6 Nov 2026: the unconstrained models ("OLS on all inputs", "Random forest on all inputs")
# are no longer drawn (a freely fitted model is not a transferable score); their results
# stay in benchmark_summary.csv and Supplementary Data Table 1. Set FIG2_UNCONSTRAINED <- TRUE
# before sourcing to draw them again.
if (exists("FIG2_UNCONSTRAINED") && isTRUE(FIG2_UNCONSTRAINED)) ROWS <- dplyr::bind_rows(ROWS, tibble::tribble(
  ~model,                         ~label,                          ~group,
  "OLS on all inputs",            "Linear model, all inputs",      "Unconstrained model",
  "Random forest on all inputs",  "Random forest, all inputs",     "Unconstrained model"))
miss <- setdiff(ROWS$model, bs$model)
if (length(miss)) stop("Not in benchmark_summary.csv: ", paste(miss, collapse = ", "))

METRICS <- c(within_r = "r within sites", sites_correct_of_all = "Sites ranked in the right direction (%)",
             between_r = "r between site means", R2 = "Overall R²")

dd <- bs %>% filter(model %in% ROWS$model) %>%
  select(model, all_of(unlist(lapply(names(METRICS), function(m) c(m, paste0(m, c("_lo", "_hi"))))))) %>%
  pivot_longer(-model, names_to = c("metric", "stat"), names_pattern = "^(.*?)(_lo|_hi)?$") %>%
  mutate(stat = case_when(is.na(stat) | stat == "" ~ "med", stat == "_lo" ~ "lo", TRUE ~ "hi")) %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  mutate(across(c(med, lo, hi), ~ ifelse(metric == "sites_correct_of_all", 100 * .x, .x))) %>%
  left_join(ROWS, by = "model") %>%
  mutate(label  = factor(label, rev(ROWS$label)),
         group  = factor(group, unique(ROWS$group)),
         metric = factor(METRICS[metric], METRICS))

# 10 Nov 2026: Fig. 1c shows the two within-site measures only (between-site r and R2 are
# given in the text and Supplementary Table 1); FIG1C_ALL_METRICS <- TRUE draws all four.
dd_all <- dd
if (!(exists("FIG1C_ALL_METRICS") && isTRUE(FIG1C_ALL_METRICS)))
  dd <- dd %>% filter(metric %in% METRICS[c("within_r", "sites_correct_of_all")]) %>%
    mutate(metric = factor(ifelse(metric == METRICS[["within_r"]], "r within sites", "Sites ranked in the\nright direction (%)"),
                           c("r within sites", "Sites ranked in the\nright direction (%)")))
X_LAB_1C <- if (exists("FIG1C_ALL_METRICS") && isTRUE(FIG1C_ALL_METRICS)) NULL else "Median and 95% range across repeated folds\n(sites withheld from fitting)"
band <- dd %>% filter(model == "SHMI")

grp_col <- c("SHMI" = "#1b6ca8", "Single part" = "#8a8a8a", "Unconstrained model" = "#e6a141")[levels(dd$group)]
pd <- ggplot(dd, aes(med, label, colour = group)) +
  geom_rect(data = band, inherit.aes = FALSE, aes(xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf),
            fill = "#1b6ca8", alpha = 0.15) +
  geom_vline(data = band, aes(xintercept = med), colour = "#1b6ca8", linewidth = 0.3) +
  geom_linerange(aes(xmin = lo, xmax = hi), linewidth = 0.6) +
  geom_point(size = 2.4) +
  facet_wrap(~ metric, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = grp_col, name = NULL) +
  labs(title = "a", x = if (exists("X_LAB_1C") && !is.null(X_LAB_1C)) X_LAB_1C else "Median and 95% range across repeated fold assignments (sites withheld from fitting)", y = NULL) +
  tag_theme +
  theme(legend.position = "bottom", strip.background = element_blank(),
        strip.text = element_text(size = 10, face = "bold"),
        panel.grid.major.y = element_line(colour = "grey92"),
        panel.border = element_rect(fill = NA, colour = "grey60"), axis.line = element_blank())

# Numbers for the text and caption
cat("\nPanel a (median [95% range])\n")
print(dd_all %>% arrange(metric, desc(label)) %>%
        transmute(metric, label, value = sprintf("%.2f [%.2f, %.2f]", med, lo, hi)) %>% as.data.frame(),
      row.names = FALSE)
if ("SHMI, weights refitted in each fold" %in% bs$model) {
  fx <- bs %>% filter(model %in% c("SHMI", "SHMI, weights refitted in each fold"))
  cat("\nPublished vs refitted weights (for the caption):\n")
  print(fx %>% select(model, R2, between_r, within_r, sites_correct_of_all) %>%
          mutate(across(where(is.numeric), ~ round(.x, 3))), row.names = FALSE)
}

# ------------------------------------------------------------------------------
# Assemble
# ------------------------------------------------------------------------------
# the comparison with the parts leads; the fit panels support it, smaller
fig <- patchwork::wrap_elements(full = pd) / (pa | pb | pc) + patchwork::plot_layout(heights = c(1, 0.75))
ggsave(file.path(OUT, "Fig2.png"), fig, width = 13, height = 8, dpi = 300)
ggsave(file.path(OUT, "Fig2.pdf"), fig, width = 13, height = 8)
cat("\nSaved Fig2.png/.pdf in", OUT, "\n")

# 7 Nov 2026 (main figures merged to fit Nature's page budget): panel a goes into Fig. 1
# (as 1d, drawn by C1 from this rds), and b-d become an Extended Data figure, relettered a-c.
saveRDS(pd + labs(title = NULL), file.path(OUT, "Fig2a_panel.rds"))
ed_val <- (pa + labs(title = "a")) | (pb + labs(title = "b")) | (pc + labs(title = "c"))
ggsave(file.path(OUT, "ED_Fig_validation.png"), ed_val, width = 13, height = 4.4, dpi = 300, bg = "white")
cat("Saved Fig2a_panel.rds (for C1) and ED_Fig_validation.png in", OUT, "\n")
