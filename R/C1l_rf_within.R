################################################################################
# C1l_rf_within.R  (21 Nov 2026)
#
# Stricter within-site benchmark for SHMI: does a flexible model of the same
# management records explain more of the WITHIN-site variation in lv_SH?
#
# The random forest in A4a was fitted to total variation, which is mostly between
# sites, so it may have spent its flexibility there. Here every model is fitted and
# scored on within-site variation only, with each site left out in turn (LOSO,
# the same scheme as C1j and C1k):
#
#   SHMI            site-centred lv_SH on site-centred SHMI (published weights),
#                   one common slope fitted without the held-out site
#   RF_centred      random forest: site-centred lv_SH on the seven site-centred inputs
#   RF_centred_ctx  as RF_centred, plus each input's site mean as context, so the
#                   forest can let a within-site difference depend on the site's
#                   overall management (e.g. adding cover where cover is already high)
#   RF_raw          the A4a forest (raw lv_SH on raw inputs), scored within sites,
#                   for comparison with A4a
#
# Inputs (all from the records): NonGrow, Grow, Diversity_plants, Richness, InvDist,
# p_amend, p_animal -- the same seven as A4a. Site weights 1/n as in A4a.
# Within-site R2 = 1 - SSE / SS of site-centred lv_SH, over all held-out predictions
# (each prediction is centred on its site, so between-site error is removed).
# Site bootstrap (B_BOOT) of the difference from SHMI, on the held-out predictions.
#
# Reads   Outputs/Final_plus/final_results.rds (A3a)
# Writes  Outputs/Benchmark_within/rf_within.csv, rf_within_predictions.csv
# Runs in a few minutes (4 x 74 forest fits).
################################################################################
suppressPackageStartupMessages({ library(dplyr); library(ranger) })
setwd(here::here())
OUT <- "Outputs/Benchmark_within"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
if (!exists("B_BOOT"))   B_BOOT   <- 2000
if (!exists("RF_TREES")) RF_TREES <- 500
set.seed(1)

res <- readRDS("Outputs/Final_plus/final_results.rds")
us  <- as_tibble(res$us)
VARS <- c("NonGrow", "Grow", "Diversity_plants", "Richness", "InvDist", "p_amend", "p_animal")
stopifnot(all(c("MGT_combo", "sitecode", "SH", VARS) %in% names(us)))
if (!"SHMI" %in% names(us)) {
  sc <- as_tibble(res$scores) %>% select(MGT_combo, SHMI)
  us <- us %>% left_join(sc, by = "MGT_combo")
}
us <- us %>% filter(!is.na(SH), !is.na(SHMI)) %>% mutate(sitecode = as.character(sitecode))
message(sprintf("C1l: %d systems at %d sites", nrow(us), n_distinct(us$sitecode)))

centre <- function(x, g) x - ave(x, g)
d <- us
d$y_c    <- centre(d$SH, d$sitecode)
d$SHMI_c <- centre(d$SHMI, d$sitecode)
for (v in VARS) {
  d[[paste0(v, "_c")]] <- centre(d[[v]], d$sitecode)       # within-site difference
  d[[paste0(v, "_m")]] <- ave(d[[v]], d$sitecode)          # site context (records only)
}
d$w <- 1 / as.numeric(table(d$sitecode)[d$sitecode])
V_C   <- paste0(VARS, "_c")
V_CTX <- c(V_C, paste0(VARS, "_m"))

rf_fit_pred <- function(tr, te, yvar, xvars) {
  f <- ranger(reformulate(xvars, yvar), data = as.data.frame(tr[, c(yvar, xvars)]),
              case.weights = tr$w, num.trees = RF_TREES, seed = 1, num.threads = 1)
  predict(f, data = as.data.frame(te[, xvars]))$predictions
}

pred <- d %>% select(MGT_combo, sitecode, y_c) %>%
  mutate(SHMI = NA_real_, RF_centred = NA_real_, RF_centred_ctx = NA_real_, RF_raw = NA_real_)
for (s in unique(d$sitecode)) {
  te <- d$sitecode == s; tr <- d[!te, ]; ts <- d[te, ]
  # SHMI: weighted slope through the origin on centred values
  b <- sum(tr$w * tr$SHMI_c * tr$y_c) / sum(tr$w * tr$SHMI_c^2)
  pred$SHMI[te]           <- b * ts$SHMI_c
  pred$RF_centred[te]     <- rf_fit_pred(tr, ts, "y_c", V_C)
  pred$RF_centred_ctx[te] <- rf_fit_pred(tr, ts, "y_c", V_CTX)
  p_raw <- rf_fit_pred(tr, ts, "SH", VARS)
  pred$RF_raw[te]         <- p_raw - mean(p_raw)                     # within-site part only
}
# centre every model's predictions on its site, so only within-site ranking is scored
MODELS <- c("SHMI", "RF_centred", "RF_centred_ctx", "RF_raw")
for (m in MODELS) pred[[m]] <- centre(pred[[m]], pred$sitecode)
write.csv(pred, file.path(OUT, "rf_within_predictions.csv"), row.names = FALSE)

r2w <- function(p, idx = seq_len(nrow(p)), m) {
  y <- p$y_c[idx]; 1 - sum((y - p[[m]][idx])^2) / sum(y^2)
}
est <- sapply(MODELS, function(m) r2w(pred, m = m))
sites <- unique(pred$sitecode)
boot <- t(replicate(B_BOOT, {
  pick <- sample(sites, length(sites), replace = TRUE)
  idx <- unlist(lapply(pick, function(s) which(pred$sitecode == s)))
  sapply(MODELS, function(m) r2w(pred, idx, m))
}))
out <- data.frame(
  model = MODELS, within_R2 = est,
  lo = apply(boot, 2, quantile, 0.025), hi = apply(boot, 2, quantile, 0.975),
  diff_from_SHMI = est - est[["SHMI"]],
  diff_lo = apply(boot - boot[, "SHMI"], 2, quantile, 0.025),
  diff_hi = apply(boot - boot[, "SHMI"], 2, quantile, 0.975),
  share_boot_SHMI_better = colMeans(boot[, "SHMI"] > boot),
  n_systems = nrow(d), n_sites = length(sites), row.names = NULL)
out$share_boot_SHMI_better[out$model == "SHMI"] <- NA
write.csv(out, file.path(OUT, "rf_within.csv"), row.names = FALSE)
print(out, digits = 3)
