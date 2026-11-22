################################################################################
# C1c_ED_Fig5.R
#
# Extended Data Fig. 5: SHMI from the records a farm or a satellite can provide.
# Four record types only (the components removed one at a time are in
# Supplementary Data Table 2, written by C2_export_tables.R).
#   a  what each record keeps (one row per record type)
#   b  how SHMI from that record performs: against soil health (NAPESHM, out-of-fold,
#      published weights; same measures and order as Fig. 2a) and the within-site corn
#      yield slope (Fig. 3a model). Median and 95% range over repeats and seeds;
#      dashed line: full record.
#
# Drawn from the saved B9 results, so B9 (hours) need not be rerun:
#   Input   Outputs/Coarse/coarse_napeshm.csv, Outputs/Coarse/coarse_yield.csv
#   Output  Outputs/Manuscript_work/ED_Fig5.png, ED_Fig5.pdf
#           (also Outputs/Coarse/Fig_coarse_records.png, the working copy)
################################################################################
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })
if (!requireNamespace("patchwork", quietly = TRUE)) stop("Install patchwork: install.packages('patchwork')")
CO  <- here::here("Outputs", "Coarse")
OUT <- here::here("Outputs", "Manuscript_work")   # working copies; C3 copies the paper's items to Outputs/Manuscript/numbered
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
f_nap <- file.path(CO, "coarse_napeshm.csv"); f_yld <- file.path(CO, "coarse_yield.csv")

if (!file.exists(f_nap)) {
  message("Extended Data Fig. 5 not drawn: ", f_nap, " not found (run B9 once).")
} else {
  REC <- c(D0 = "Full record", FARM = "Typical farm record",
           RSM = "Remote sensing +\nreported manure, grazing", RS = "Remote sensing only")
  lv  <- rev(unname(REC))                                   # full record on top

  # ---- a: what each record keeps -------------------------------------------------
  # from the scenario definitions in B9 (SCEN: D0; FARM tillage,sjitter:14;
  # RSM cover,maincrop,tillage,jitter:14; RS the same with noamend)
  keep <- tibble::tribble(
    ~scenario, ~Crops,           ~`Cover crops`,  ~Dates,                 ~Tillage,         ~`Manure,\ngrazing`,
    "D0",      "all species",    "species",       "exact",                "each pass"     , "recorded",
    "FARM",    "all species",    "species",       "±14 d, farm",          "class",          "recorded",
    "RSM",     "main crop",      "presence",      "±14 d, field",         "class",          "reported",
    "RS",      "main crop",      "presence",      "±14 d, field",         "class",          "none")
  ka <- keep %>% pivot_longer(-scenario, names_to = "detail", values_to = "txt") %>%
    mutate(level = case_when(txt %in% c("all species", "species", "exact", "each pass", "recorded", "reported") ~ "As in the full record",
                             txt == "none" ~ "Missing", TRUE ~ "Simplified"),
           rec = factor(REC[scenario], lv), detail = factor(detail, names(keep)[-1]))
  pa <- ggplot(ka, aes(detail, rec, fill = level)) +
    geom_tile(colour = "white", linewidth = 1.2) +
    geom_text(aes(label = txt), size = 2.9, lineheight = 0.9) +
    scale_fill_manual(values = c("As in the full record" = "#c7dcef", "Simplified" = "#fde0b6", "Missing" = "#f4b6b6"),
                      breaks = c("As in the full record", "Simplified", "Missing"), name = NULL) +
    scale_x_discrete(position = "top", expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
    labs(title = "a", x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 14), plot.title.position = "plot",
          panel.grid = element_blank(), axis.text.x.top = element_text(face = "bold", size = 9),
          axis.text.y = element_text(size = 10, colour = "black"), legend.position = "bottom")

  # ---- b: performance -------------------------------------------------------------
  ns <- read.csv(f_nap, stringsAsFactors = FALSE) %>% filter(weights == "published", scenario %in% names(REC))
  M <- c("Within sites: r", "Sites ranked in the\nright direction (%)", "Between site means: r",
         "Overall R²", "Corn yield per 10 SHMI\n(kg ha⁻¹, within sites)")
  pl <- list(
    ns %>% transmute(scenario, measure = M[1], est = within_r_med, lo = within_r_lo, hi = within_r_hi),
    ns %>% transmute(scenario, measure = M[2], est = 100 * sites_correct_med, lo = 100 * sites_correct_lo, hi = 100 * sites_correct_hi),
    ns %>% transmute(scenario, measure = M[3], est = between_r_med, lo = between_r_lo, hi = between_r_hi),
    ns %>% transmute(scenario, measure = M[4], est = R2_med, lo = R2_lo, hi = R2_hi))
  if (file.exists(f_yld)) {
    pl[[5]] <- read.csv(f_yld, stringsAsFactors = FALSE) %>% filter(scenario %in% names(REC)) %>%
      transmute(scenario, measure = M[5], est = slope, lo = lower, hi = upper)
  } else message("  ", f_yld, " not found: yield column left out")
  fd <- bind_rows(pl) %>% mutate(rec = factor(REC[scenario], lv), measure = factor(measure, M))
  ref <- fd %>% filter(scenario == "D0") %>% select(measure, ref = est)
  pb <- ggplot(fd, aes(est, rec)) +
    geom_vline(data = ref, aes(xintercept = ref), linetype = "dashed", colour = "grey55") +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.25, orientation = "y", colour = "#1b6ca8") +
    geom_point(size = 2.4, colour = "#1b6ca8") +
    facet_wrap(~ measure, nrow = 1, scales = "free_x") +
    scale_x_continuous(expand = expansion(mult = 0.12), n.breaks = 4) +
    labs(title = "b", x = NULL, y = NULL) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 14), plot.title.position = "plot",
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          strip.background = element_blank(), strip.text = element_text(face = "bold", size = 9, hjust = 0),
          panel.spacing.x = unit(1.1, "lines"), panel.grid.minor = element_blank())

  # 12 Nov 2026: panels for R/figures/F6_records.R (Fig. 6), without their letter titles
  PNL_DIR <- here::here("Outputs", "panels"); dir.create(PNL_DIR, recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(a = pa + labs(title = NULL), b = pb + labs(title = NULL)), file.path(PNL_DIR, "panel_records.rds"))
  fig <- patchwork::wrap_plots(pa, pb, widths = c(1.35, 2.65))
  for (f in c(file.path(OUT, "ED_Fig5.png"), file.path(CO, "Fig_coarse_records.png")))
    ggsave(f, fig, width = 15, height = 3.6, dpi = 300)
  tryCatch(ggsave(file.path(OUT, "ED_Fig5.pdf"), fig, width = 15, height = 3.6, device = cairo_pdf),
           error = function(e) message("  ED_Fig5.pdf not written (cairo_pdf unavailable): ", conditionMessage(e)))

  cat("\nExtended Data Fig. 5 (median [95% range])\n")
  print(fd %>% arrange(measure, desc(rec)) %>%
          transmute(measure = gsub("\n", " ", measure), record = gsub("\n", " ", rec),
                    value = sprintf("%.3g [%.3g, %.3g]", est, lo, hi)) %>% as.data.frame(), row.names = FALSE)
  message("Saved ED_Fig5.png/.pdf in ", OUT)
}
