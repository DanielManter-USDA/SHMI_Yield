################################################################################
# C1b_ED_Fig2.R
#
# Extended Data Fig. 2: the active-season share within Cover/Roots across the 1,000
# calibration fits (violin and box), with the final value and its 95% range in red.
# The three main weights (panel a of the working plot F1_weights.png made by
# A3a_SHMI_calibration_plus.R) are shown in Fig. 1a, so they are left out here.
#
# Drawn from the saved calibration, so the calibration need not be rerun:
#   Input   Outputs/Final_plus/final_results.rds   (A3a: $fits, $weights)
#   Output  Outputs/Manuscript/ED_Fig2.png, ED_Fig2.pdf
################################################################################
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })
FIN <- here::here("Outputs", "Final_plus", "final_results.rds")
OUT <- here::here("Outputs", "Manuscript")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(FIN)) {
  message("Extended Data Fig. 2 not drawn: ", FIN, " not found (run Part A once).")
} else {
  fr  <- readRDS(FIN)
  gs  <- fr$weights %>% filter(weight == "growing_share")
  fit <- data.frame(x = "", growing_share = fr$fits$growing_share)
  p <- ggplot(fit, aes(x, growing_share)) +
    geom_violin(fill = "#cfe8c8", colour = NA) +
    geom_boxplot(width = 0.15, outlier.size = 0.4) +
    geom_errorbar(data = gs, aes(x = "", y = final, ymin = lo95, ymax = hi95), width = 0.12, colour = "#d1495b") +
    geom_point(data = gs, aes(x = "", y = final), colour = "#d1495b", size = 3) +
    geom_text(data = gs, aes(x = "", y = final, label = sprintf("%.3f", final)),
              nudge_x = 0.3, colour = "#d1495b", size = 3.6) +
    scale_y_continuous(limits = c(0.5, 1)) +
    labs(x = NULL, y = "Active-season share of Cover/Roots") +
    theme_bw(base_size = 11) + theme(axis.ticks.x = element_blank())
  ggsave(file.path(OUT, "ED_Fig2.png"), p, width = 3.2, height = 4, dpi = 300)
  ggsave(file.path(OUT, "ED_Fig2.pdf"), p, width = 3.2, height = 4)
  message(sprintf("Extended Data Fig. 2: active-season share %.3f (95%% range %.3f-%.3f), %d fits",
                  gs$final, gs$lo95, gs$hi95, nrow(fr$fits)))
}
