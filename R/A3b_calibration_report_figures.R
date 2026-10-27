################################################################################
# A3b_calibration_report_figures.R
#
# Report figures 5 and 6 and Table 10b for SHMI and SHMI_base (from
# A3a_SHMI_calibration_plus.R). Each figure has one row per score and three panels:
# (a) all treatments, (b) site means, (c) each treatment minus its site mean.
#   Figure 5  raw scores (0-100) against measured soil health (0-100); slopes are
#             soil-health points per score point
#   Figure 6  out-of-fold predicted soil health (each site predicted with weights
#             fitted without it, averaged over 100 repeats) against measured
# Table 10b agreement with measured soil health on the 1:1 line (bias, RMSE, r,
#           Lin's concordance, and the Kobayashi & Salam (2000) split into
#           offset, spread and scatter), US and Mexico, both scores
#
# Inputs (A3a_SHMI_calibration_plus.R): FIN/us_oof_predictions.csv,
#   FIN/shmi_scores.csv, FIN/final_results.xlsx (mexico_units). Outputs written
#   before the rename (shmi_plus_scores.csv, SHMI / SHMI_plus columns) are read
#   too and renamed: SHMI -> SHMI_base, SHMI_plus -> SHMI.
# Outputs: FIN/fig5_insample.png, fig6_oof.png, table10b_agreement.csv
# FIN defaults to Outputs/Final_plus; set FIN_DIR before sourcing to change it.
################################################################################
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(here)
})
if (!requireNamespace("patchwork", quietly = TRUE)) stop("Install patchwork: install.packages('patchwork')")
# Manuscript figures: each panel carries only its bold letter (a-f); titles go in the legend.
# options(shmi.ms_figures = FALSE) restores the descriptive panel titles and the figure title.
MS_FIG <- isTRUE(getOption("shmi.ms_figures", TRUE))
FIN <- if (exists("FIN_DIR")) FIN_DIR else here::here("Outputs", "Final_plus")

# Older outputs: SHMI (no bonus) -> SHMI_base, SHMI_plus -> SHMI, *_core -> *_base
new_names <- function(x) {
  x <- as.data.frame(x)
  if ("SHMI_plus" %in% names(x)) names(x)[names(x) == "SHMI"] <- "SHMI_base"
  names(x)[names(x) == "SHMI_plus"] <- "SHMI"
  names(x) <- sub("_core$", "_base", names(x))
  x
}
sc_file <- file.path(FIN, "shmi_scores.csv")
if (!file.exists(sc_file)) sc_file <- file.path(FIN, "shmi_plus_scores.csv")
oof <- new_names(read.csv(file.path(FIN, "us_oof_predictions.csv")))
sc  <- new_names(read.csv(sc_file))
d <- oof %>% inner_join(sc %>% select(MGT_combo, SHMI_base, SHMI), by = "MGT_combo")
stopifnot(nrow(d) == nrow(oof))

# One panel: points, fitted line, and r / R2 / slope in the corner
panel <- function(x, y, tag, title, xlab, ylab) {
  r <- cor(x, y); slope <- unname(coef(lm(y ~ x))[2])
  ggplot(tibble(x = x, y = y), aes(x, y)) +
    geom_point(colour = "#1f6fab", alpha = 0.6, size = 1.4, stroke = 0) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "#d1495b", linewidth = 1) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.08, vjust = 1.15, size = 3.4, lineheight = 0.95,
             label = sprintf("r = %.2f\nR² = %.2f\nslope = %.2f", r, r^2, slope)) +
    labs(title = if (MS_FIG) tag else sprintf("(%s) %s", tag, title), x = xlab, y = ylab) +
    theme_classic(base_size = 11) +
    theme(plot.title = if (MS_FIG) element_text(size = 14, face = "bold", hjust = 0, margin = margin(b = 4))
                       else element_text(size = 11),
          plot.title.position = "plot")
}

# Three panels for one x variable: all treatments, site means, within sites
three_panels <- function(d, xvar, xname, tags) {
  sm <- d %>% group_by(sitecode) %>% summarise(x = mean(.data[[xvar]]), y = mean(SH), .groups = "drop")
  wi <- d %>% group_by(sitecode) %>% mutate(x = .data[[xvar]] - mean(.data[[xvar]]), y = SH - mean(SH)) %>% ungroup()
  list(panel(d[[xvar]], d$SH, tags[1], sprintf("%s: all %d treatments", xname, nrow(d)), xname, "Soil health"),
       panel(sm$x, sm$y, tags[2], sprintf("%s: %d site means", xname, nrow(sm)),
             paste("Site mean", xname), "Site mean soil health"),
       panel(wi$x, wi$y, tags[3], sprintf("%s: within sites", xname),
             paste(xname, "minus site mean"), "Soil health minus site mean"))
}
plot_rows <- function(rows, title, file) {
  p <- patchwork::wrap_plots(unlist(rows, recursive = FALSE), ncol = 3)
  if (!MS_FIG) p <- p + patchwork::plot_annotation(title = title, theme = theme(plot.title = element_text(hjust = 0.5, size = 13)))
  ggsave(file, p, width = 13, height = 4.2 * length(rows), dpi = 300)
  invisible(p)
}

plot_rows(list(three_panels(d, "SHMI", "SHMI", c("a", "b", "c")),
               three_panels(d, "SHMI_base", "SHMI_base", c("d", "e", "f"))),
          "US sites: SHMI and SHMI_base against measured soil health (final weights, all US data)",
          file.path(FIN, "fig5_insample.png"))
d <- d %>% mutate(pred_SHMI = pred_oof, pred_SHMI_base = pred_oof_base)
plot_rows(list(three_panels(d, "pred_SHMI", "Predicted SH from SHMI", c("a", "b", "c")),
               three_panels(d, "pred_SHMI_base", "Predicted SH from SHMI_base", c("d", "e", "f"))),
          "US sites: each site predicted by weights fitted without it (out-of-fold)",
          file.path(FIN, "fig6_oof.png"))

# ==============================================================================
# Table 10b: agreement with measured soil health on the 1:1 line
# ==============================================================================
# Spreads use the population form (divide by n), as in Lin (1989) and
# Kobayashi & Salam (2000), so the three parts add exactly to the MSD.
agreement <- function(x, y, label) {
  n <- length(x); mx <- mean(x); my <- mean(y)
  sx <- sqrt(sum((x - mx)^2) / n); sy <- sqrt(sum((y - my)^2) / n)
  r <- cor(x, y); msd <- mean((x - y)^2)
  tibble(comparison = label, n = n, mean_x = mx, mean_SH = my, bias = mx - my,
         sd_x = sx, sd_SH = sy, RMSE = sqrt(msd), r = r,
         CCC = 2 * r * sx * sy / (sx^2 + sy^2 + (mx - my)^2),
         share_offset = (mx - my)^2 / msd, share_spread = (sx - sy)^2 / msd,
         share_scatter = 2 * sx * sy * (1 - r) / msd)
}
mx <- new_names(readxl::read_excel(file.path(FIN, "final_results.xlsx"), sheet = "mexico_units"))
agree <- bind_rows(
  agreement(d$SHMI,           d$SH,  "US treatments, raw SHMI"),
  agreement(d$SHMI_base,      d$SH,  "US treatments, raw SHMI_base"),
  agreement(d$pred_SHMI,      d$SH,  "US treatments, predicted SH from SHMI (out-of-fold)"),
  agreement(d$pred_SHMI_base, d$SH,  "US treatments, predicted SH from SHMI_base (out-of-fold)"),
  agreement(mx$SHMI,          mx$SH, "Mexico, raw SHMI"),
  agreement(mx$SHMI_base,     mx$SH, "Mexico, raw SHMI_base"),
  agreement(mx$pred,          mx$SH, "Mexico, predicted SH from SHMI"),
  agreement(mx$pred_base,     mx$SH, "Mexico, predicted SH from SHMI_base"))
stopifnot(all(abs(agree$share_offset + agree$share_spread + agree$share_scatter - 1) < 1e-8))
print(agree %>% mutate(across(where(is.numeric), ~ round(.x, 2))), width = Inf)
write.csv(agree, file.path(FIN, "table10b_agreement.csv"), row.names = FALSE)

cat("Saved fig5_insample.png, fig6_oof.png and table10b_agreement.csv in", FIN, "\n")
