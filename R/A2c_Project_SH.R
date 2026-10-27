################################################################################
# A2c_Project_SH.R
#
# Scores samples the SEM was NOT fitted on, using everything learned from the
# training units only (nothing is refitted):
#   1. outliers set to NA with the training quartile bounds (3 x IQR)
#   2. climate/texture part removed with the training adjustment models
#   3. adjusted indicators standardized with the training mean and SD
#   4. factor scores from the training SEM (regression scores, FIML)
#   5. 0-100 with the training minimum and range (values outside 0-100 are
#      kept: they mean the new sample lies beyond the training range)
# Default: the noMX (US) model applied to the Mexican samples.
#
# Input : Outputs/semDF_prepped.rds, Outputs/sem_noMX_projection.rds
# Output: Outputs/sem_noMX_projected_MX.rds (and .csv)
################################################################################

suppressPackageStartupMessages({ library(dplyr); library(lavaan); library(mgcv); library(here) })
OUT <- function(...) here::here("Outputs", ...)
source(here::here("R", "datasets.R"))

PROJ_FROM <- "noMX"                                   # model the scores come from
TARGET    <- function(d) is_mx(d)                     # samples to score

proj <- readRDS(OUT(paste0("sem_", PROJ_FROM, "_projection.rds")))
dat  <- readRDS(OUT("semDF_prepped.rds")) %>% mutate(sitecode = as.character(sitecode))
dat  <- dat[TARGET(dat), , drop = FALSE]
cat("Scoring", nrow(dat), "units at", length(unique(dat$sitecode)), "sites with the", PROJ_FROM, "model\n")

ind <- proj$indicators
d <- dat
# 1. outliers with the training bounds
n_out <- 0
for (v in ind) {
  b <- proj$bounds[[v]]; x <- d[[v]]
  bad <- !is.na(x) & (x < b[1] | x > b[2]); n_out <- n_out + sum(bad)
  x[bad] <- NA; d[[v]] <- x
}
cat("Values beyond the training outlier bounds (set to NA):", n_out, "\n")
# 2-3. adjustment with the training models, then training standardization
for (v in ind) {
  keep <- complete.cases(d[, c(v, proj$adj_vars[proj$adj_vars != "sitecode"])])
  out  <- rep(NA_real_, nrow(d))
  if (any(keep)) {
    nd <- as.data.frame(d[keep, , drop = FALSE]); nd$site_f <- factor(nd$sitecode)
    tt <- mgcv::predict.gam(proj$models[[v]], newdata = nd, type = "terms")   # gam models: use mgcv's method
    env <- rowSums(tt[, !grepl("site_f", colnames(tt)), drop = FALSE])
    out[keep] <- d[[v]][keep] - env
  }
  d[[v]] <- (out - proj$center[[v]]) / proj$sd[[v]]
}
# 4. factor scores from the training SEM
has_any <- rowSums(!is.na(d[, ind])) > 0
nd <- as.data.frame(d[has_any, c(ind, "sitecode")])
lv <- lavaan::lavPredict(proj$fit, newdata = nd, type = "lv")
# 5. training 0-100 scale
scores <- d %>% select(MGT_combo, system, sitecode) %>% mutate(n_missing = rowSums(is.na(d[, ind])))
for (nm in colnames(lv)) {
  raw <- rep(NA_real_, nrow(d)); raw[has_any] <- lv[, nm]
  scores[[paste0(nm, "_raw")]] <- raw
  scores[[nm]] <- 100 * (raw - proj$lv_min[[nm]]) / proj$lv_range[[nm]]
}
cat("\nProjected SH (training 0-100 scale):\n"); print(summary(scores$SH))
cat("Units outside 0-100:", sum(scores$SH < 0 | scores$SH > 100, na.rm = TRUE), "\n")

saveRDS(scores, OUT(paste0("sem_", PROJ_FROM, "_projected_MX.rds")))
write.csv(scores, OUT(paste0("sem_", PROJ_FROM, "_projected_MX.csv")), row.names = FALSE)
cat("Saved:", OUT(paste0("sem_", PROJ_FROM, "_projected_MX.rds")), "\n")
