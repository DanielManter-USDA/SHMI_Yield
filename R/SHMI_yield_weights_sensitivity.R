################################################################################
# SHMI_yield_weights_sensitivity.R
#
# Does the research-site yield result depend on which SHMI weights are used?
# Runs SHMI_yield_main.R once per weight set, changing only the weights used to
# score the research plots (the Minnesota farms keep the official weights), and
# compares the within-site results:
#   - the within-site SHMI slope (kg/ha per 10 SHMI points at 150 kg N/ha, average
#     system age) and the SHMI x N and SHMI x age terms (the Fig 1-3 model)
#   - each site's own slope and the random-effects pool across sites
#   - the same slopes per within-site SD of SHMI, because a weight set that spreads
#     the plots further apart on the SHMI scale gives a smaller slope per point for
#     the same yield differences
#
# Weight sets (from the calibration scenarios, Diversity fixed at 0.15):
#   all      0.687 / 0.150 / 0.163, growing-season share 0.760 (official)
#   within   0.541 / 0.150 / 0.309, growing-season share 0.612 (fitted within sites)
#   between  0.792 / 0.150 / 0.058, growing-season share 0.780 (fitted to site means)
#
# RUN: put this next to SHMI_yield_main.R (the version with SHMI_WEIGHTS_OVERRIDE
# and OUT_OVERRIDE) and source it. Each weight set reruns the whole main script
# (downloads are cached), writing its usual outputs to Outputs/yield_main_<set>[_agcros]/.
# Results: Outputs/yield_weights[_agcros]/ weights_slopes.csv, weights_site_slopes.csv,
#          Fig_weights_site_slopes.png
################################################################################
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(purrr); library(ggplot2) })

MAIN_SCRIPT <- if (exists("MAIN_SCRIPT")) MAIN_SCRIPT else here::here("R", "SHMI_yield_main.R")
stopifnot(file.exists(MAIN_SCRIPT))
if (!any(grepl("SHMI_WEIGHTS_OVERRIDE", readLines(MAIN_SCRIPT))))
  stop("SHMI_yield_main.R needs the version with SHMI_WEIGHTS_OVERRIDE: ", MAIN_SCRIPT)
MAIN_SCRIPT <- normalizePath(MAIN_SCRIPT)

SETS <- list(
  all     = NULL,                                                            # official weights
  within  = list(weights = c(Cover = 0.541, Diversity = 0.150, InvDist = 0.309), growing_share = 0.612),
  between = list(weights = c(Cover = 0.792, Diversity = 0.150, InvDist = 0.058), growing_share = 0.780)
)
RUN_SETS <- if (exists("RUN_SETS")) RUN_SETS else c("all", "within", "between")

# Each set runs the main script in its own environment, whose parent is the global
# environment. The main script tests for *_OVERRIDE objects with exists(), which also
# looks in the global environment, so an override left there from an earlier run would
# apply to every set (a global SHMI_WEIGHTS_OVERRIDE would silently turn the "all" set
# into another weight set). Weights and output folder must therefore not be set globally;
# the scenario, AgCROS and refinement overrides apply to all three sets alike, so they
# are allowed but reported.
leak <- intersect(c("SHMI_WEIGHTS_OVERRIDE", "OUT_OVERRIDE"), ls(globalenv()))
if (length(leak)) stop("Remove these from the global environment first: ", paste(leak, collapse = ", "),
                       "  (rm(", paste(leak, collapse = ", "), "))")
shared <- intersect(c("SCENARIO_OVERRIDE", "INCLUDE_AGCROS_OVERRIDE", "AGCROS_REFINE_OVERRIDE"), ls(globalenv()))
if (length(shared)) message("Applied to every weight set: ", paste(shared, collapse = ", "))

# Results kept from each run (the full environments hold every site's inputs and would
# use several GB for three sets)
KEEP <- c("tt_full", "corn", "site_rows", "pool", "INCLUDE_AGCROS", "OUT")
run_set <- function(set) {
  cat(sprintf("\n==================== Weight set: %s ====================\n", set))
  env <- new.env(parent = globalenv())
  env$OUT_OVERRIDE <- file.path("Outputs", paste0("yield_main_", set))
  if (!is.null(SETS[[set]])) env$SHMI_WEIGHTS_OVERRIDE <- SETS[[set]]
  sys.source(MAIN_SCRIPT, envir = env)        # the main script sets the working directory itself
  # the site runs come from the disk cache when present; its key includes the weights,
  # so each set reuses only runs scored with its own weights
  out <- mget(KEEP, envir = env)
  rm(env); gc(verbose = FALSE)
  cat("Weight set", set, "->", out$OUT, "\n")
  out
}
envs <- setNames(lapply(RUN_SETS, run_set), RUN_SETS)
site_sets <- lapply(envs, function(e) sort(unique(as.character(e$corn$site))))
if (length(unique(site_sets)) > 1) stop("The weight sets used different sites: ",
                                        paste(names(site_sets), sapply(site_sets, paste, collapse = ","), collapse = " | "))
# relative to the main script's working directory; "_agcros" when the AgCROS sites were included
OUTW <- file.path("Outputs", paste0("yield_weights", if (isTRUE(envs[[1]]$INCLUDE_AGCROS)) "_agcros" else ""))
cat("Sites used:", paste(site_sets[[1]], collapse = ", "), "\n")
dir.create(OUTW, recursive = TRUE, showWarnings = FALSE)

# ---- overall within-site results --------------------------------------------------
overall <- imap_dfr(envs, function(e, set) {
  tt <- e$tt_full
  s  <- list(b = tt["SHMI_w", "Value"], se = tt["SHMI_w", "Std.Error"])   # slope at 150 kg N/ha, average age
  # within-site SD of SHMI: plots (one row per management unit) around their site mean
  sd_w <- e$corn %>% distinct(site, MGT_combo, SHMI_w) %>% pull(SHMI_w) %>% sd()
  tibble(weights = set,
         slope_per10 = 10 * s$b, lo = 10 * (s$b - qt(0.975, tt["SHMI_w", "DF"]) * s$se), hi = 10 * (s$b + qt(0.975, tt["SHMI_w", "DF"]) * s$se),
         p = tt["SHMI_w", "p-value"],
         SHMIxN_per100kgN = 10 * tt["SHMI_w:N_c", "Value"], p_SHMIxN = tt["SHMI_w:N_c", "p-value"],
         SHMIxage_per_yr = 10 * tt["SHMI_w:age_c", "Value"], p_SHMIxage = tt["SHMI_w:age_c", "p-value"],
         within_sd_SHMI = sd_w, slope_per_sd = s$b * sd_w,
         pooled_site_slope_per10 = as.numeric(e$pool$b), pooled_lo = e$pool$ci.lb, pooled_hi = e$pool$ci.ub,
         I2 = e$pool$I2)
})
write.csv(overall, file.path(OUTW, "weights_slopes.csv"), row.names = FALSE)
cat("\nWithin-site yield results by weight set (kg/ha; slope at 150 kg N/ha and average age)\n")
print(overall %>% mutate(across(where(is.numeric), ~ signif(.x, 3))) %>% as.data.frame(), row.names = FALSE)

# ---- each site's slope ---------------------------------------------------------------
site_sd <- function(e) e$corn %>% distinct(site, MGT_combo, SHMI_w) %>% group_by(site) %>%
  summarise(within_sd = sd(SHMI_w), SHMI_range = diff(range(SHMI_w)), .groups = "drop")
sites <- imap_dfr(envs, function(e, set) {
  e$site_rows %>% select(site = label, effect, se, p_val, n) %>%
    left_join(site_sd(e), by = "site") %>%
    mutate(weights = set, per_sd = effect / 10 * within_sd)
})
write.csv(sites, file.path(OUTW, "weights_site_slopes.csv"), row.names = FALSE)
cat("\nSite slopes (kg/ha per 10 SHMI) by weight set\n")
print(sites %>% select(site, weights, effect) %>% mutate(effect = round(effect)) %>%
        pivot_wider(names_from = weights, values_from = effect) %>% as.data.frame(), row.names = FALSE)
cat("\nWithin-site spread of SHMI (range among a site's treatments) by weight set\n")
print(sites %>% select(site, weights, SHMI_range) %>% mutate(SHMI_range = round(SHMI_range, 1)) %>%
        pivot_wider(names_from = weights, values_from = SHMI_range) %>% as.data.frame(), row.names = FALSE)

# ---- figure: site slopes and pooled slope by weight set ------------------------------
lab <- c(all = "All variation (official)", within = "Within sites", between = "Between sites")
cols <- c(all = "#7a7a7a", within = "#08519c", between = "#9ecae1")
pd <- bind_rows(
  sites %>% transmute(site, weights, est = effect, lo = effect - 1.96 * se, hi = effect + 1.96 * se),
  overall %>% transmute(site = "Pooled (Fig 1 model)", weights, est = slope_per10, lo, hi)) %>%
  mutate(site = factor(site, rev(c("Pooled (Fig 1 model)", sort(unique(sites$site))))),
         weights = factor(weights, names(lab)))
p <- ggplot(pd, aes(est, site, colour = weights)) +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.6), size = 0.3) +
  scale_colour_manual(values = cols, labels = lab, name = "SHMI weights") +
  labs(x = "Yield change per 10 SHMI points (kg/ha, 95% CI)", y = NULL,
       title = "Within-site yield slope under three SHMI weight sets") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(file.path(OUTW, "Fig_weights_site_slopes.png"), p, width = 7.5, height = 5, dpi = 300)
cat("\nSaved in", normalizePath(OUTW), "\n")
