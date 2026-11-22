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
# 12 Nov 2026: panels for the figure scripts in R/figures (Fig. 1c and 3b)
PNL_DIR <- here::here("Outputs", "panels"); dir.create(PNL_DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(sh$plots, file.path(PNL_DIR, "panel_SHMI_vs_SH.rds"))      # list: all, between, within
saveRDS(sp$plots, file.path(PNL_DIR, "panel_SHMI_vs_SHAPE.rds"))

# ---- SHMI, its parts and individual indicators (16 Nov 2026) --------------------------------
# Within site (site means removed), across all systems: the correlation of SHMI and of each part
# alone (Cover/Roots, Diversity, Inverse Disturbance, amendment bonus) with lv_SH and with single
# measured indicators; and the best that indicator-specific weights could do, by non-negative
# least squares on the four parts, in sample and with each site left out in turn (weights refitted
# without the site). Asks whether separate principles, or weights tuned to each indicator, beat
# the one set of soil-health weights.
IND  <- c("lv_SH" = "SH", "Soil organic carbon" = "SOC_raw", "Wet aggregate stability" = "KSSL_WAS_raw",
          "Respiration" = "Resp_raw", "Active carbon (POXC)" = "AC_raw", "ACE protein" = "ACE_raw")
PRTS <- c("Cover/Roots" = "Cover", "Diversity" = "Diversity", "Inverse Disturbance" = "InvDist",
          "Amendments" = "Amend_bonus")
nnls_fit <- function(X, y) {                      # non-negative least squares, base R
  f <- function(b) sum((y - X %*% b)^2)
  optim(rep(0.1, ncol(X)), f, method = "L-BFGS-B", lower = 0)$par
}
ind_rows <- list()
for (nm in names(IND)) {
  yv <- IND[[nm]]
  if (!yv %in% names(d)) { message("  indicator not found: ", yv); next }
  z <- d %>% filter(!is.na(.data[[yv]]), !is.na(SHMI), if_all(all_of(unname(PRTS)), ~ !is.na(.x))) %>%
    group_by(sitecode) %>% filter(n() >= 2) %>% ungroup()
  wd <- function(v) z[[v]] - ave(z[[v]], z$sitecode)
  y <- wd(yv); X <- sapply(unname(PRTS), wd)
  b <- nnls_fit(X, y)
  pred <- numeric(length(y))
  for (s in unique(z$sitecode)) { m <- z$sitecode == s; pred[m] <- X[m, , drop = FALSE] %*% nnls_fit(X[!m, , drop = FALSE], y[!m]) }
  ind_rows[[nm]] <- tibble(indicator = nm, systems = nrow(z), sites = n_distinct(z$sitecode),
                           r_SHMI = cor(wd("SHMI"), y),
                           !!!setNames(lapply(unname(PRTS), function(p) cor(wd(p), y)), paste0("r_", names(PRTS))),
                           r_best_in = cor(X %*% b, y), r_best_site_out = cor(pred, y))
}
ind_tab <- bind_rows(ind_rows)
write.csv(mutate(ind_tab, across(where(is.double), ~ round(.x, 3))), file.path(OUT, "SHMI_parts_by_indicator.csv"), row.names = FALSE)
cat("\n== Within-site r with lv_SH and single indicators: SHMI, each part, indicator-specific weights ==\n")
print(as.data.frame(ind_tab), digits = 2)
