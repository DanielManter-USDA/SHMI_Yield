################################################################################
# B2c_SHMI_record_length.R   (Fig. 3b)
#
# How many years of management record does the yield relationship need?
# B2 scores each plot once, over the experiment's whole period. Here each corn
# plot-year instead gets an SHMI computed from only the K crop years up to that
# harvest (crop year = 1 October to 30 September, so the window ends before the
# harvest and excludes the following fall's operations):
#   K = 1   that crop year's management only (the cover crop before the corn,
#           its tillage, the corn itself)
#   K = 2, 3, 4
#   full    B2's score (the whole period)
# Every version is fitted with B2's model (site x year, N terms, other, SHMI_w,
# SHMI_w x N, SHMI_w x age; random treatment/plot; AR(1)) on the SAME plot-years:
# those whose plot has at least max(K) crop years of record before the harvest.
# Compared by the slope per 10 points at 150 kg N/ha and by AIC (maximum likelihood,
# same data and fixed-effect structure, so the versions differ only in the score).
#
# Input   Outputs/yield_main_default_total_agcros/B2_data.rds (B2)
# Outputs Outputs/yield_record_length/
#   rolling_shmi.rds           SHMI by unit, harvest year and K (cached; slow to build)
#   record_length_results.csv  slope, CI, P, AIC and delta AIC by version
#   record_length_paired.csv   each window on its own largest sample vs the full score there
#   record_length_sites.csv    the same by experiment (each experiment's own model)
#   rolling_vs_full.csv        how far the K-year score is from the full score, by experiment
#   Fig_rolling_shmi.png       each treatment's K-year score by year (full score dashed)
#   Fig_record_length.png      Fig. 3b (assembled by B7): slope by record length (paired comparison)
################################################################################
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(nlme); library(SHMI)
})
setwd(here::here())
IN  <- "Outputs/yield_main_default_total_agcros/B2_data.rds"
OUT <- "Outputs/yield_record_length"
if (!file.exists(IN)) stop(IN, " not found: run B2_SHMI_yield_main.R first")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
K_SET   <- c(1, 2, 3, 4)
if (exists("RL_K_OVERRIDE")) K_SET <- RL_K_OVERRIDE                # e.g. c(1, 2, 4, 8)
if (exists("RL_K_OVERRIDE")) OUT <- paste0(OUT, "_K", paste(K_SET, collapse = "-"))   # keep each set's cache apart
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
N_REF   <- 150
REBUILD <- exists("RL_REBUILD") && isTRUE(RL_REBUILD)   # TRUE: recompute the rolling scores
ctrl <- lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200)

D <- readRDS(IN)
runs <- D$runs; corn <- D$corn

# ---- rolling SHMI ---------------------------------------------------------------------
yr <- function(d) as.integer(format(as.Date(d), "%Y"))
shmi_window <- function(inputs, climate, start, end) {       # as B2's shmi_for_period
  rb <- inputs$rot_bounds
  rb$rot_start <- pmax(as.Date(rb$rot_start), as.Date(start))
  rb$rot_end   <- pmin(as.Date(rb$rot_end),   as.Date(end))
  rb <- rb[rb$rot_start <= rb$rot_end, ]
  if (!nrow(rb)) return(NULL)
  rb$rot_start_yr <- yr(rb$rot_start); rb$rot_end_yr <- yr(rb$rot_end)
  inputs$rot_bounds <- rb
  inputs$mgt <- inputs$mgt[inputs$mgt$MGT_combo %in% rb$MGT_combo, ]
  tryCatch(build_shmi(inputs, climate)$indicator_df %>% select(MGT_combo, SHMI),
           error = function(e) { message("    ", start, " to ", end, ": ", conditionMessage(e)); NULL })
}
roll_file <- file.path(OUT, "rolling_shmi.rds")
if (file.exists(roll_file) && !REBUILD) {
  roll <- readRDS(roll_file)
} else {
  roll <- imap_dfr(runs, function(r, s) {
    yrs <- sort(unique(corn$year[corn$site == s]))
    if (!length(yrs)) return(NULL)
    message(s, ": ", length(yrs), " harvest years x ", length(K_SET), " windows")
    rec_start <- min(as.Date(r$inputs$rot_bounds$rot_start))
    tidyr::crossing(year = yrs, K = K_SET) %>%
      pmap_dfr(function(year, K) {
        st <- as.Date(sprintf("%d-10-01", year - K)); en <- as.Date(sprintf("%d-09-30", year))
        x <- shmi_window(r$inputs, r$climate, st, en)
        if (is.null(x)) return(NULL)
        x %>% mutate(site = s, year = year, K = K)
      })
  })
  saveRDS(roll, roll_file)
}

# a window is complete when it starts no earlier than the record (allowing the record to
# start up to three months into the first crop year, e.g. on 1 January): with K = 4, a
# plot's first three harvests have no complete window
rec_start <- vapply(runs, function(r) as.numeric(min(as.Date(r$inputs$rot_bounds$rot_start))), numeric(1))
roll <- roll %>% select(-any_of("complete")) %>%
  mutate(complete = as.numeric(as.Date(sprintf("%d-10-01", year - K))) >= rec_start[site] - 92)

# ---- how much does the score change with the window? ----------------------------------
# For each experiment and K: how far the K-year score is from the full-period score
# (RMS difference and mean difference, SHMI points), whether it ranks the units the
# same way in each harvest year (median Spearman correlation with the full score), and
# how much a unit's K-year score moves from year to year (median SD over years, which
# shows rotation phases and practice changes). Units and years with corn in B2.
full <- imap_dfr(runs, function(r, s) r$shmi %>% distinct(MGT_combo, SHMI) %>% transmute(MGT_combo, site = s, SHMI_full = SHMI))
cmp <- roll %>% semi_join(distinct(corn, MGT_combo, year), by = c("MGT_combo", "year")) %>%
  inner_join(full, by = c("MGT_combo", "site"))
rank_r <- function(a, b) if (length(a) >= 3 && sd(a) > 0 && sd(b) > 0) cor(a, b, method = "spearman") else NA_real_
chg <- cmp %>% group_by(site, K) %>%
  summarise(units = n_distinct(MGT_combo), plot_years = n(),
            rms_diff = sqrt(mean((SHMI - SHMI_full)^2)), mean_diff = mean(SHMI - SHMI_full),
            full_range = diff(range(SHMI_full)), .groups = "drop") %>%
  left_join(cmp %>% group_by(site, K, year) %>% summarise(rho = rank_r(SHMI, SHMI_full), .groups = "drop") %>%
              group_by(site, K) %>% summarise(rank_rho_median = median(rho, na.rm = TRUE), .groups = "drop"),
            by = c("site", "K")) %>%
  left_join(cmp %>% group_by(site, K, MGT_combo) %>% summarise(sd_years = sd(SHMI), .groups = "drop") %>%
              group_by(site, K) %>% summarise(sd_over_years_median = median(sd_years, na.rm = TRUE), .groups = "drop"),
            by = c("site", "K"))
write.csv(mutate(chg, across(where(is.double), ~ round(.x, 2))), file.path(OUT, "rolling_vs_full.csv"), row.names = FALSE)
cat("\nK-year score against the full-period score, by experiment:\n"); print(chg, n = Inf, width = Inf)

# figure: each treatment's K-year score by harvest year (mean over its plots), full score dashed
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  trt_of <- distinct(corn, MGT_combo, MGT_trt_base)
  fd <- cmp %>% inner_join(trt_of, by = "MGT_combo") %>% group_by(site, K, year, MGT_trt_base) %>%
    summarise(SHMI = mean(SHMI), SHMI_full = mean(SHMI_full), .groups = "drop") %>%
    mutate(K = factor(paste0("last ", K, ifelse(K == 1, " year", " years")), paste0("last ", K_SET, ifelse(K_SET == 1, " year", " years"))))
  g <- ggplot(fd, aes(year, SHMI, colour = MGT_trt_base)) +
    geom_line(aes(y = SHMI_full), linetype = "dashed", linewidth = 0.3) +
    geom_line(linewidth = 0.4) + facet_grid(site ~ K, scales = "free_x") +
    scale_y_continuous(limits = c(0, 100)) + labs(x = "Harvest year", y = "SHMI") +
    theme_bw(base_size = 8) + theme(legend.position = "none", strip.text.y = element_text(angle = 0))
  ggsave(file.path(OUT, "Fig_rolling_shmi.png"), g, width = 9, height = 1.2 * n_distinct(fd$site) + 0.8, dpi = 200)
}

# ---- data: B2's corn plot-years with each version of the score ------------------------
# common sample: plot-years whose largest window lies inside the record
ok <- roll %>% filter(K == max(K_SET), complete) %>% distinct(MGT_combo, year)
base <- corn %>% semi_join(ok, by = c("MGT_combo", "year"))
message(nrow(base), " of ", nrow(corn), " corn plot-years have ", max(K_SET), " crop years of record before harvest")
versions <- c(setNames(as.character(K_SET), paste0("K", K_SET)), full = "full")
data_for <- function(v) {
  d <- if (v == "full") base else
    base %>% select(-SHMI) %>%
      inner_join(roll %>% filter(K == as.numeric(v)) %>% select(MGT_combo, year, SHMI), by = c("MGT_combo", "year"))
  d %>% group_by(site) %>% mutate(SHMI_w = SHMI - mean(SHMI)) %>% ungroup() %>% droplevels()
}
dl <- lapply(versions, data_for)
stopifnot(length(unique(vapply(dl, nrow, 1L))) == 1)   # identical plot-years in every version

# ---- model (B2's, with aliased columns dropped as in B2b) -----------------------------
varies <- function(d, v) { d <- droplevels(d); any(tapply(as.numeric(d[[v]]), d$site_year, function(x) diff(range(x))) > 1e-8, na.rm = TRUE) }
n_terms <- function(d) { v <- grep("^N[12]_", names(d), value = TRUE); v[vapply(v, function(k) varies(d, k), logical(1))] }
fac_term <- function(d, v) if (n_distinct(d[[v]]) > 1 && varies(d, v)) v else character()
nm <- function(x) make.names(gsub("[^A-Za-z0-9_.]", "_", x), unique = TRUE)
fit <- function(d, rhs, method = "ML") {
  X <- model.matrix(reformulate(rhs), d)
  qx <- qr(X, tol = 1e-7); keep <- sort(qx$pivot[seq_len(qx$rank)])
  drop <- setdiff(colnames(X), colnames(X)[keep])
  if (any(grepl("^SHMI_w", drop))) stop("SHMI term aliased: ", paste(grep("^SHMI_w", drop, value = TRUE), collapse = ", "))
  X <- X[, keep, drop = FALSE]; colnames(X) <- nm(colnames(X))
  dd <- data.frame(yield_kgha = d$yield_kgha, trt = d$trt, plot_id = d$plot_id, year = d$year, X, check.names = FALSE)
  lme(reformulate(c("0", colnames(X)), "yield_kgha"), random = ~ 1 | trt/plot_id,
      correlation = corCAR1(form = ~ year | trt/plot_id), data = dd, method = method, control = ctrl)
}
rhs_for <- function(d) c("site_year", n_terms(d), fac_term(d, "other"), "SHMI_w",
                         if (varies(d, "N_c")) "SHMI_w:N_c", "SHMI_w:age_c")
one <- function(d, v) {
  m  <- fit(d, rhs_for(d), "ML")                      # AIC: same data, same terms
  mr <- fit(d, rhs_for(d), "REML")                    # slope as reported in B2
  tt <- summary(mr)$tTable["SHMI_w", ]
  b <- 10 * tt[["Value"]]; se <- 10 * tt[["Std.Error"]]; df <- tt[["DF"]]
  tibble(version = v, plot_years = nrow(d), sites = n_distinct(d$site), SHMI_sd_within = sd(d$SHMI_w),
         slope = b, lower = b - qt(0.975, df) * se, upper = b + qt(0.975, df) * se, p = tt[["p-value"]],
         slope_per_sd = b / 10 * sd(d$SHMI_w), AIC = AIC(m))
}
# each window on its own largest sample (plot-years with a complete K-year window),
# against the full score on the same plot-years: uses the early years that the common
# sample drops (with K = 1, every year)
paired <- map_dfr(K_SET, function(k) {
  okk <- roll %>% filter(K == k, complete) %>% distinct(MGT_combo, year)
  b <- corn %>% semi_join(okk, by = c("MGT_combo", "year"))
  dk <- b %>% select(-SHMI) %>% inner_join(roll %>% filter(K == k) %>% select(MGT_combo, year, SHMI), by = c("MGT_combo", "year"))
  b  <- b %>% semi_join(dk, by = c("MGT_combo", "year"))
  prep <- function(d) d %>% group_by(site) %>% mutate(SHMI_w = SHMI - mean(SHMI)) %>% ungroup() %>% droplevels()
  bind_rows(one(prep(dk), paste0("K", k)), one(prep(b), "full")) %>%
    mutate(window = k, dAIC_vs_full = AIC - AIC[version == "full"])
})
write.csv(mutate(paired, across(where(is.double), ~ signif(.x, 4))), file.path(OUT, "record_length_paired.csv"), row.names = FALSE)
cat("\nEach window on its own largest sample, against the full score on the same plot-years:\n"); print(paired, width = Inf)

res <- imap_dfr(dl, function(d, v) one(d, v)) %>%
  mutate(dAIC_vs_full = AIC - AIC[version == "full"])
write.csv(mutate(res, across(where(is.double), ~ signif(.x, 4))), file.path(OUT, "record_length_results.csv"), row.names = FALSE)
cat("\nCorn yield and SHMI from the last K crop years (same plot-years):\n"); print(res, width = Inf)

# each experiment's own model (slope at its own N; no interactions)
sites <- imap_dfr(dl, function(d, v) split(d, d$site) %>% imap_dfr(function(ds, s) {
  ds <- droplevels(ds)
  m <- tryCatch(fit(ds, c("site_year", n_terms(ds), fac_term(ds, "other"), "SHMI_w"), "REML"),
                error = function(e) { message("  ", s, " ", v, ": ", conditionMessage(e)); NULL })
  if (is.null(m)) return(NULL)
  tt <- summary(m)$tTable["SHMI_w", ]
  tibble(version = v, site = s, plot_years = nrow(ds), slope = 10 * tt[["Value"]], se = 10 * tt[["Std.Error"]], p = tt[["p-value"]])
}))
write.csv(mutate(sites, across(where(is.double), ~ signif(.x, 4))), file.path(OUT, "record_length_sites.csv"), row.names = FALSE)
# ---- Fig. 3b: the yield slope by length of record ------------------------
# each window on its own largest sample (filled) against the full-record score on the
# same plot-years (open); 95% CI
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  fd <- paired %>% mutate(score = ifelse(version == "full", "Full record, same plot-years", "Last K crop years only"),
                          score = factor(score, c("Last K crop years only", "Full record, same plot-years")),
                          x = factor(window, K_SET, ifelse(K_SET == 1, "1 season", paste(K_SET, "years"))))
  # plot-years under each window's name, not next to the points
  npy <- fd %>% filter(version != "full") %>% arrange(window)
  levels(fd$x) <- paste0(levels(fd$x), "\n(", format(npy$plot_years, big.mark = ","), " plot-years)")
  g <- ggplot(fd, aes(x, slope, colour = score, shape = score)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey55") +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15, position = position_dodge(width = 0.4)) +
    geom_point(size = 2.8, fill = "white", position = position_dodge(width = 0.4)) +
    scale_colour_manual(values = c("#1f6fab", "grey35"), name = NULL) +
    scale_shape_manual(values = c(16, 21), name = NULL) +
    labs(x = "Management record used to score SHMI before each harvest",
         y = expression("Corn yield change per 10 SHMI (kg ha"^-1*")")) +
    theme_bw(base_size = 10) + theme(legend.position = "bottom", panel.grid.minor = element_blank())
  ggsave(file.path(OUT, "Fig_record_length.png"), g, width = 5.5, height = 4, dpi = 300)
  ggsave(file.path(OUT, "Fig_record_length.pdf"), g, width = 5.5, height = 4)
}
message("Record-length test done: ", OUT)
