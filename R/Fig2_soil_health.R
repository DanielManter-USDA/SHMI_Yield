################################################################################
# Fig2_soil_health.R
#
# Manuscript Fig. 2: SHMI predicts soil health at sites withheld from
# calibration, and better than any of its parts.
#   a  out-of-fold predicted vs measured soil health, all 354 systems
#   b  the same for the 74 site means
#   c  the same within sites (deviations from site means)
#   d  out-of-fold performance of SHMI, each part alone, and unconstrained
#      models on the same inputs (from SHMI_benchmark.R), four metrics;
#      shaded band = SHMI's 95% range across repeated fold assignments
#
# Inputs
#   FIN/us_oof_predictions.csv          (SHMI_calibration_plus.R)
#   Outputs/Benchmark/benchmark_summary.csv (SHMI_benchmark.R)
# Output
#   Outputs/Manuscript/Fig2_soil_health.png (and .pdf)
#
# Run after SHMI_calibration_plus.R and SHMI_benchmark.R.
################################################################################
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(here)
})
if (!requireNamespace("patchwork", quietly = TRUE)) stop("Install patchwork: install.packages('patchwork')")

FIN   <- if (exists("FIN_DIR")) FIN_DIR else here::here("Outputs", "Final_plus")
BENCH <- here::here("Outputs", "Benchmark", "benchmark_summary.csv")
OUT   <- here::here("Outputs", "Manuscript")
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
pa <- panel(d$pred, d$SH, "a", "Predicted soil health", "Measured soil health")
pb <- panel(sm$x, sm$y, "b", "Site mean predicted soil health", "Site mean measured soil health")
pc <- panel(wi$x, wi$y, "c", "Predicted, minus site mean", "Measured, minus site mean")
cat(sprintf("a: n = %d systems | b: %d site means\n", nrow(d), nrow(sm)))

# ------------------------------------------------------------------------------
# d: SHMI against its parts and unconstrained models
# ------------------------------------------------------------------------------
bs <- read.csv(BENCH, check.names = FALSE)

ROWS <- tibble::tribble(                       # benchmark name -> label, group (top to bottom)
  ~model,                         ~label,                          ~group,
  "SHMI",                         "SHMI",                          "SHMI",
  "Cover/Roots",                  "Cover/Roots alone",             "Single part",
  "Diversity",                    "Diversity alone",               "Single part",
  "InvDist",                      "Inverse Disturbance alone",     "Single part",
  "OLS on all inputs",            "Linear model, all inputs",      "Unconstrained model",
  "Random forest on all inputs",  "Random forest, all inputs",     "Unconstrained model")
miss <- setdiff(ROWS$model, bs$model)
if (length(miss)) stop("Not in benchmark_summary.csv: ", paste(miss, collapse = ", "))

METRICS <- c(R2 = "Overall R²", between_r = "r between site means",
             within_r = "r within sites", sites_correct_of_all = "Sites ranked correctly (%)")

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
band <- dd %>% filter(model == "SHMI")

grp_col <- c("SHMI" = "#1b6ca8", "Single part" = "#8a8a8a", "Unconstrained model" = "#e6a141")
pd <- ggplot(dd, aes(med, label, colour = group)) +
  geom_rect(data = band, inherit.aes = FALSE, aes(xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf),
            fill = "#1b6ca8", alpha = 0.15) +
  geom_vline(data = band, aes(xintercept = med), colour = "#1b6ca8", linewidth = 0.3) +
  geom_linerange(aes(xmin = lo, xmax = hi), linewidth = 0.6) +
  geom_point(size = 2.4) +
  facet_wrap(~ metric, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = grp_col, name = NULL) +
  labs(title = "d", x = "Median and 95% range across repeated fold assignments", y = NULL) +
  tag_theme +
  theme(legend.position = "bottom", strip.background = element_blank(),
        strip.text = element_text(size = 10, face = "bold"),
        panel.grid.major.y = element_line(colour = "grey92"),
        panel.border = element_rect(fill = NA, colour = "grey60"), axis.line = element_blank())

# Numbers for the text and caption
cat("\nPanel d (median [95% range])\n")
print(dd %>% arrange(metric, desc(label)) %>%
        transmute(metric, label, value = sprintf("%.2f [%.2f, %.2f]", med, lo, hi)) %>% as.data.frame(),
      row.names = FALSE)
if ("SHMI, published weights (fixed)" %in% bs$model) {
  fx <- bs %>% filter(model %in% c("SHMI", "SHMI, published weights (fixed)"))
  cat("\nRefitted vs published weights (for the caption):\n")
  print(fx %>% select(model, R2, between_r, within_r, sites_correct_of_all) %>%
          mutate(across(where(is.numeric), ~ round(.x, 3))), row.names = FALSE)
}

# ------------------------------------------------------------------------------
# Assemble
# ------------------------------------------------------------------------------
fig <- (pa | pb | pc) / pd + patchwork::plot_layout(heights = c(1, 0.85))
ggsave(file.path(OUT, "Fig2_soil_health.png"), fig, width = 13, height = 8.6, dpi = 300)
ggsave(file.path(OUT, "Fig2_soil_health.pdf"), fig, width = 13, height = 8.6)
cat("\nSaved Fig2_soil_health.png/.pdf in", OUT, "\n")
