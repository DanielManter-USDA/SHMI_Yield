################################################################################
# A5b_SuppFig2_weights_by_comparison.R
#
# Supplementary Data Fig. 2: SHMI weights fitted between sites, to all variation
# (the calibration) and within sites, at Diversity 0.10-0.25.
#   a-c  Inverse disturbance weight, Cover weight, growing-season share
#        (medians and 95% ranges across 1,000 fits); open red diamonds = 57:43
#        blend of the between- and within-site fits (shares of lv_SH variance)
#   d-f  Out-of-fold R2, between-site r, within-site r (medians of 100 repeats)
#
# Input  Outputs/Calib_scenarios/scenario_summary.csv (A5a_Run_within_calibration_scenarios.R)
# Output Outputs/Calib_scenarios/SuppFig2_weights_by_comparison.png / .pdf
################################################################################
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(here); library(patchwork) })

IN  <- here::here("Outputs", "Calib_scenarios", "scenario_summary.csv")
OUT <- dirname(IN)
BETWEEN_SHARE <- 0.57          # share of lv_SH variance between sites (US systems)

t <- read.csv(IN) %>%
  mutate(fit = as.character(fit),
         fit = factor(c(between = "Between\nsites", `0` = "All variation\n(calibrated)", `1` = "Within\nsites")[fit],
                      c("Between\nsites", "All variation\n(calibrated)", "Within\nsites")),
         div = factor(sprintf("%.2f", fixed_div))) %>%
  filter(!is.na(fit))
div_col <- setNames(c("#9ecae1", "#4292c6", "#2171b5", "#08306b")[seq_along(levels(t$div))], levels(t$div))
dodge <- position_dodge(width = 0.55)
theme_ed <- theme_classic(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 9), legend.position = "bottom",
        axis.text.x = element_text(size = 7.5))

# blend: 57% between-site fit + 43% within-site fit, at the all-variation position
blend <- t %>% filter(grepl("^Between|^Within", fit)) %>%
  select(div, fit, cover, disturbance) %>%
  mutate(w = ifelse(grepl("^Between", fit), BETWEEN_SHARE, 1 - BETWEEN_SHARE)) %>%
  group_by(div) %>% summarise(cover = sum(w * cover), disturbance = sum(w * disturbance), .groups = "drop") %>%
  mutate(fit = factor("All variation\n(calibrated)", levels(t$fit)))

wpanel <- function(v, title, ylim, show_blend = TRUE) {
  p <- ggplot(t, aes(fit, .data[[v]], colour = div)) +
    geom_linerange(aes(ymin = .data[[paste0(v, "_lo")]], ymax = .data[[paste0(v, "_hi")]]), position = dodge) +
    geom_point(position = dodge, size = 1.6)
  if (show_blend) p <- p + geom_point(data = blend, aes(fit, .data[[v]], group = div), position = dodge,
                                      shape = 23, size = 2.4, colour = "#d1495b", fill = NA, stroke = 0.7,
                                      inherit.aes = FALSE)
  p + scale_colour_manual(values = div_col, name = "Diversity weight") +
    coord_cartesian(ylim = ylim) + labs(x = NULL, y = NULL, title = title) + theme_ed
}
pa <- wpanel("disturbance", "a  Inverse disturbance weight", c(0, 0.42)) + labs(y = "Weight (median, 95% range)")
pb <- wpanel("cover", "b  Cover weight", c(0.3, 0.95))
pc <- wpanel("growing_share", "c  Growing-season share", c(0.3, 0.95), show_blend = FALSE)

short <- function(f) factor(sub("\n.*", "", f), c("Between", "All variation", "Within"))
ppanel <- function(v, title) {
  ggplot(t, aes(short(fit), .data[[v]], colour = div, group = div)) + geom_line() + geom_point(size = 1.4) +
    scale_colour_manual(values = div_col, name = "Diversity weight") +
    scale_x_discrete(labels = c("Between", "All", "Within")) +
    labs(x = NULL, y = NULL, title = title) + theme_ed
}
pd <- ppanel("R2", "d  Out-of-fold R²")
pe <- ppanel("between_r", "e  Between-site r")
pf <- ppanel("within_r", "f  Within-site r")

MS_FIG <- isTRUE(getOption("shmi.ms_figures", TRUE))   # manuscript figures: no titles/captions (legend in the manuscript)
fig <- (pa | pb | pc) / (pd | pe | pf) +
  plot_layout(guides = "collect") +
  plot_annotation(caption = if (MS_FIG) NULL else "Weights fitted to: between sites (site means), all variation, or within sites (deviations from site means). Red diamonds: 57:43 blend of the between- and within-site fits.",
                  theme = theme(plot.caption = element_text(size = 7.5, hjust = 0))) &
  theme(legend.position = "bottom")
ggsave(file.path(OUT, "SuppFig2_weights_by_comparison.png"), fig, width = 7.2, height = 5.2, dpi = 300)
ggsave(file.path(OUT, "SuppFig2_weights_by_comparison.pdf"), fig, width = 7.2, height = 5.2)
cat("Saved SuppFig2_weights_by_comparison.png/.pdf in", OUT, "\n")
