################################################################################
# C1_Fig1_path.R  (21 Nov 2026; replaces C1_Fig1_and_Table1.R, now in R/archive)
#
# Fig. 1b: how SHMI scores a system. The path diagram among its sub-indices:
#   the formula weights of SHMI 2.0.0 (blue) and the paths estimated by least squares on
#   the 354 US NAPESHM systems (black): Cover/Roots ~ Diversity + Inverse Disturbance and
#   soil health (lv_SH) ~ SHMI, with 95% intervals from a site-level bootstrap; the total
#   effect of each sub-index on SHMI (weight + path to Cover/Roots x its weight); in-sample
#   and out-of-fold R2 (weights refitted in each fold, A4a benchmark).
#
# Inputs  data/SHAPE_by_unit_with_SHMI.csv (sub-indices, SHMI and lv_SH by system)
#         Outputs/Benchmark/benchmark_summary.csv (out-of-fold R2; optional)
# Outputs Outputs/panels/panel_path.rds (drawn into Fig. 1 by R/figures/F1_soil_health_model.R)
#         Outputs/Manuscript_work/Fig1_path_coefficients.csv, Fig1b_path.png (working copy)
################################################################################
suppressPackageStartupMessages({
  library(SHMI); library(dplyr); library(tidyr); library(ggplot2)
})
stopifnot(packageVersion("SHMI") >= "2.0.0")
OUT <- here::here("Outputs", "Manuscript_work"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
PNL_DIR <- here::here("Outputs", "panels"); dir.create(PNL_DIR, recursive = TRUE, showWarnings = FALSE)
B <- 2000                       # site-bootstrap draws
set.seed(2026)
par <- shmi_parameters()        # amend_gap for the amendment formula in the diagram

# ==============================================================================
# 1. Path coefficients (Fig. 1b)
# ==============================================================================
d <- read.csv(here::here("data", "SHAPE_by_unit_with_SHMI.csv"))
w <- shmi_weights()
paths <- function(D) {
  cr <- coef(lm(Cover ~ Diversity + InvDist, data = D))
  sh <- coef(lm(SH ~ SHMI, data = D))
  p <- c(div_to_cover = cr[["Diversity"]], inv_to_cover = cr[["InvDist"]], shmi_to_sh = sh[["SHMI"]])
  c(p, total_cover = w[["Cover"]],
    total_diversity = w[["Diversity"]] + p[["div_to_cover"]] * w[["Cover"]],
    total_invdist   = w[["InvDist"]]   + p[["inv_to_cover"]] * w[["Cover"]],
    sum_of_totals   = w[["Cover"]] + w[["Diversity"]] + p[["div_to_cover"]] * w[["Cover"]] +
                      w[["InvDist"]] + p[["inv_to_cover"]] * w[["Cover"]])
}
est <- paths(d)
by_site <- split(d, d$sitecode)
boot <- t(replicate(B, paths(do.call(rbind, by_site[sample(names(by_site), replace = TRUE)]))))
coef_tab <- tibble(term = names(est), estimate = unname(est),
                   lo95 = apply(boot, 2, quantile, 0.025), hi95 = apply(boot, 2, quantile, 0.975))
r2_in <- summary(lm(SH ~ SHMI, data = d))$r.squared
bf <- here::here("Outputs", "Benchmark", "benchmark_summary.csv")
# Out-of-fold R2 with the weights refitted in each fold (0.26), as in the calibration;
# the "SHMI" row holds the published weights fixed, so its R2 (0.27) is not fully out of fold.
r2_oof <- if (file.exists(bf)) { b <- read.csv(bf); b$R2[b$model == "SHMI, weights refitted in each fold"] } else NA_real_
coef_tab <- bind_rows(coef_tab, tibble(term = c("R2_in_sample", "R2_out_of_fold"),
                                       estimate = c(r2_in, r2_oof), lo95 = NA, hi95 = NA))
write.csv(coef_tab, file.path(OUT, "Fig1_path_coefficients.csv"), row.names = FALSE)
cat("Path coefficients (Fig. 1b):\n"); print(as.data.frame(coef_tab), digits = 3, row.names = FALSE)

# ==============================================================================
# 2. Path diagram (Fig. 1b)
# ==============================================================================
BLUE <- "#1a5ea8"; BLACK <- "#2b2b2b"; ORANGE <- "#b8641a"
f2 <- function(x) sprintf("%.2f", x)
E <- setNames(coef_tab$estimate, coef_tab$term)
boxes <- tribble(
  ~key,                  ~x,   ~y,   ~w,   ~h,  ~fill,
  "Diversity",           13,   40.8, 19,   5.6, "#cfcbf5",
  "Cover/Roots",         38.5, 40.8, 15,   5.6, "#9fe1cb",
  "Inverse Disturbance", 12.75,23.1, 19.5, 5.6, "#e5e4df",
  "SHMI_base",           55,   23.1, 14,   5.6, "#b5d4f2",
  "SHMI",                71.5, 23.1, 11,   5.6, "#8fbde6",
  "Soil health",         91,   23.1, 14,   5.0, "#f1efe8",
  "Amendments",          64,   12.3, 14,   5.6, "#f9c677")
tot <- c(Diversity = E[["total_diversity"]], `Cover/Roots` = E[["total_cover"]], `Inverse Disturbance` = E[["total_invdist"]])
boxes <- boxes %>% mutate(total = unname(tot[key]))
arrows <- tribble(
  ~x,   ~y,   ~xend, ~yend, ~col,  ~lw, ~lab,                          ~lx,  ~ly,
  22.5, 40.8, 31,    40.8,  BLACK, 1.6, f2(E[["div_to_cover"]]),       26.7, 43.4,
  22.5, 23.5, 34.5,  38,    BLACK, 0.7, f2(E[["inv_to_cover"]]),       26.8, 29.4,
  46,   38.0, 51.5,  25.9,  BLUE,  1.8, f2(w[["Cover"]]),              50.6, 33.6,
  22.5, 39.2, 48,    24.4,  BLUE,  0.7, f2(w[["Diversity"]]),          40.0, 31.4,
  22.5, 23.1, 48,    23.1,  BLUE,  0.7, f2(w[["InvDist"]]),            38.5, 24.6,
  62,   23.1, 66,    23.1,  BLUE,  1.1, "",                            NA,   NA,
  64,   15.1, 64,    20.3,  ORANGE,1.1, "",                            NA,   NA,
  77,   23.1, 84,    23.1,  BLACK, 1.6, f2(E[["shmi_to_sh"]]),         80.5, 25.8)
pa <- ggplot() +
  geom_segment(data = arrows, aes(x = x, y = y, xend = xend, yend = yend, colour = col, linewidth = lw),
               arrow = arrow(length = unit(2.6, "mm"), type = "closed"), lineend = "round") +
  geom_label(data = filter(arrows, nzchar(lab)), aes(x = lx, y = ly, label = lab, colour = col),
             fontface = "bold", size = 3.6, label.size = 0, fill = "white", label.padding = unit(0.08, "lines")) +
  geom_tile(data = boxes, aes(x = x, y = y, width = w, height = h, fill = fill), colour = "#444444", linewidth = 0.5) +
  geom_text(data = filter(boxes, key != "SHMI_base"), aes(x = x, y = y + ifelse(is.na(total), 0, 1.0), label = key),
            fontface = "bold", size = 3.6) +
  # SHMI with "base" as a subscript (plotmath)
  geom_text(data = filter(boxes, key == "SHMI_base"), aes(x = x, y = y), label = "bold(SHMI[base])",
            parse = TRUE, size = 3.6) +
  geom_text(data = filter(boxes, !is.na(total)), aes(x = x, y = y - 1.4, label = paste("total", f2(total))),
            colour = BLUE, size = 3.2) +
  annotate("text", x = 64, y = 6.5, parse = TRUE,
           label = sprintf("bold(SHMI == (1 - %.1f*italic(p)) %%*%% SHMI[base] + %d*italic(p))",
                           par$amend_gap, round(100 * par$amend_gap)),
           colour = ORANGE, size = 3.4) +
  annotate("text", x = 91, y = 18.6, hjust = 0.5, vjust = 1, size = 3, colour = "#444444", lineheight = 1.3,
           label = sprintf("R² (in-sample) = %.2f%s", r2_in,
                           if (is.na(r2_oof)) "" else sprintf("\nR² (out-of-fold) = %.2f", r2_oof))) +
  scale_colour_identity() + scale_fill_identity() + scale_linewidth_identity() +
  coord_fixed(xlim = c(0, 100), ylim = c(3, 46), expand = FALSE) +
  theme_void()

saveRDS(pa, file.path(PNL_DIR, "panel_path.rds"))
ggsave(file.path(OUT, "Fig1b_path.png"), pa, width = 12.5, height = 4.6, dpi = 300, bg = "white")
message("Saved ", file.path("Outputs", "panels", "panel_path.rds"))
