################################################################################
# C1h_SHMI_vs_soil_health.R   (10 Nov 2026; Fig. 1b and Extended Data Fig. 4)
#
# SHMI against measured soil health (lv_SH) for the 354 NAPESHM US systems at 74 sites:
# all systems, site means (between sites) and deviations from site means (within sites),
# with r, R2 and the least-squares slope in each panel. The same three panels for SHAPE
# (the mean of the SHAPE indicator scores, 0-100) form Extended Data Fig. 4 (formerly
# panels a-c of the SHMI-and-SHAPE figure; panels d-f are now Fig. 1b).
#
# Input   data/SHAPE_by_unit_with_SHMI.csv
# Output  Outputs/Practice_additivity/Fig_SHMI_vs_SH_panel.rds (Fig. 1b), Fig_SHMI_vs_SH.png,
#         Fig_SHMI_vs_SHAPE.png (Extended Data Fig. 4), SHMI_vs_SH_SHAPE.csv
# Runs in seconds.
################################################################################
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(patchwork) })
OUT <- here::here("Outputs", "Practice_additivity")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
d <- read.csv(here::here("data", "SHAPE_by_unit_with_SHMI.csv")) %>%
  mutate(SHAPE = 100 * SHAPE_mean)

three <- function(yv, ylab) {
  z <- d %>% filter(!is.na(.data[[yv]]), !is.na(SHMI)) %>% transmute(sitecode, x = SHMI, y = .data[[yv]])
  sets <- list(
    all     = z,
    between = z %>% group_by(sitecode) %>% summarise(x = mean(x), y = mean(y), .groups = "drop"),
    within  = z %>% group_by(sitecode) %>% mutate(x = x - mean(x), y = y - mean(y)) %>% ungroup())
  XL <- c(all = "SHMI", between = "Site mean SHMI", within = "SHMI minus site mean")
  YL <- c(all = ylab, between = paste("Site mean", sub("^Soil health", "soil health", ylab)), within = paste(ylab, "minus site mean"))
  stats <- list(); pl <- list()
  for (k in names(sets)) {
    s <- sets[[k]]; r <- cor(s$x, s$y); b <- coef(lm(y ~ x, s))[["x"]]
    stats[[k]] <- tibble(outcome = yv, set = k, n = nrow(s), r = r, R2 = r^2, slope = b)
    pl[[k]] <- ggplot(s, aes(x, y)) +
      geom_point(colour = "steelblue", alpha = 0.6, size = 1) +
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "#c0392b", linewidth = 0.9) +
      annotate("text", x = -Inf, y = Inf, hjust = -0.08, vjust = 1.2, size = 3.1, lineheight = 1.05,
               label = sprintf("r = %.2f\nR² = %.2f\nslope = %.2f", r, r^2, b)) +
      labs(x = XL[[k]], y = YL[[k]]) + theme_classic(base_size = 10)
  }
  list(plots = pl, stats = bind_rows(stats))
}
sh <- three("SH", "Soil health (lv_SH)")
sp <- three("SHAPE", "SHAPE")
write.csv(mutate(bind_rows(sh$stats, sp$stats), across(where(is.double), ~ round(.x, 3))),
          file.path(OUT, "SHMI_vs_SH_SHAPE.csv"), row.names = FALSE)
print(as.data.frame(bind_rows(sh$stats, sp$stats)), digits = 3, row.names = FALSE)
fig1b <- sh$plots$all | sh$plots$between | sh$plots$within
saveRDS(fig1b, file.path(OUT, "Fig_SHMI_vs_SH_panel.rds"))
ggsave(file.path(OUT, "Fig_SHMI_vs_SH.png"), fig1b, width = 11, height = 3.6, dpi = 300, bg = "white")
ed <- (sp$plots$all | sp$plots$between | sp$plots$within) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 14))
ggsave(file.path(OUT, "Fig_SHMI_vs_SHAPE.png"), ed, width = 11, height = 3.6, dpi = 300, bg = "white")
