################################################################################
# C1b_ED_Fig2.R
#
# Extended Data Fig. 2: the calibrated weights across the 1,000 calibration fits.
#   a  Cover/Roots, Diversity (fixed) and Inverse Disturbance (sum to 1)
#   b  the active-season share within Cover/Roots
# Violins and boxes show the fits; red points and bars the final value and its 95%
# range. Same panels as the working plot F1_weights.png (A3a_SHMI_calibration_plus.R),
# without titles.
#
# Drawn from the saved calibration, so the calibration need not be rerun:
#   Input   Outputs/Final_plus/final_results.rds   (A3a: $fits, $weights)
#   Output  Outputs/Manuscript/ED_Fig2.png, ED_Fig2.pdf
################################################################################
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })
if (!requireNamespace("patchwork", quietly = TRUE)) stop("Install patchwork: install.packages('patchwork')")
FIN <- here::here("Outputs", "Final_plus", "final_results.rds")
OUT <- here::here("Outputs", "Manuscript_work")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(FIN)) {
  message("Extended Data Fig. 2 not drawn: ", FIN, " not found (run Part A once).")
} else {
  fr  <- readRDS(FIN)
  RED <- "#d1495b"
  th  <- theme_bw(base_size = 11) + theme(plot.tag = element_text(face = "bold", size = 14))

  # a  the three weights (Diversity is fixed, so it has no spread across fits)
  lab <- c(cover = "Cover/Roots", diversity = "Diversity\n(fixed)", disturbance = "Inverse\nDisturbance")
  wl <- fr$fits %>% select(all_of(names(lab))) %>% pivot_longer(everything(), names_to = "weight") %>%
    mutate(weight = factor(lab[weight], lab))
  wf <- fr$weights %>% filter(weight %in% names(lab)) %>% mutate(weight = factor(lab[weight], lab))
  vary <- filter(wl, weight != lab[["diversity"]])
  pa <- ggplot(wl, aes(weight, value)) +
    geom_violin(data = vary, fill = "grey85", colour = NA, scale = "width") +
    geom_boxplot(data = vary, width = 0.15, outlier.size = 0.4) +
    geom_errorbar(data = wf, aes(y = final, ymin = lo95, ymax = hi95), width = 0.12, colour = RED) +
    geom_point(data = wf, aes(y = final), colour = RED, size = 3) +
    geom_text(data = wf, aes(y = final, label = sprintf("%.3f", final)), nudge_x = 0.32, colour = RED, size = 3.6) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = NULL, y = "Weight") + th

  # b  the active-season share within Cover/Roots
  gs <- fr$weights %>% filter(weight == "growing_share")
  fit <- data.frame(x = "Active-season share\nof Cover/Roots", growing_share = fr$fits$growing_share)
  pb <- ggplot(fit, aes(x, growing_share)) +
    geom_violin(fill = "#cfe8c8", colour = NA) +
    geom_boxplot(width = 0.15, outlier.size = 0.4) +
    geom_errorbar(data = gs, aes(x = fit$x[1], y = final, ymin = lo95, ymax = hi95), width = 0.12, colour = RED) +
    geom_point(data = gs, aes(x = fit$x[1], y = final), colour = RED, size = 3) +
    geom_text(data = gs, aes(x = fit$x[1], y = final, label = sprintf("%.3f", final)),
              nudge_x = 0.32, colour = RED, size = 3.6) +
    scale_y_continuous(limits = c(0.5, 1)) +
    labs(x = NULL, y = "Share") + th

  p <- patchwork::wrap_plots(pa, pb, widths = c(3, 1.2)) + patchwork::plot_annotation(tag_levels = "a")
  ggsave(file.path(OUT, "ED_Fig2.png"), p, width = 9, height = 4.2, dpi = 300)
  ggsave(file.path(OUT, "ED_Fig2.pdf"), p, width = 9, height = 4.2)
  message(sprintf("Extended Data Fig. 2 (%d fits): %s; active-season share %.3f (%.3f-%.3f)", nrow(fr$fits),
                  paste(sprintf("%s %.3f (%.3f-%.3f)", wf$weight, wf$final, wf$lo95, wf$hi95), collapse = "; "),
                  gs$final, gs$lo95, gs$hi95))
}
