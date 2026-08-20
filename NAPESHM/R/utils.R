############################################################
# utils.R — Shared helper functions
############################################################

# Outlier replacement within climate zones
replace_outliers <- function(x) {
  if (sum(!is.na(x)) < 20) return(x)
  Q1 <- quantile(x, 0.25, na.rm = TRUE)
  Q3 <- quantile(x, 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1
  x[(x < Q1 - 3 * IQR) | (x > Q3 + 3 * IQR)] <- NA
  x
}

# Min-max scaling
minmax_scale <- function(x, to = c(0, 100), na.rm = TRUE) {
  rng <- range(x, na.rm = na.rm)
  (x - rng[1]) / (rng[2] - rng[1]) * (to[2] - to[1]) + to[1]
}

# Indicator adjustment function
adjust_indicator <- function(df, var) {
  keep <- complete.cases(df[, c(var, "ClimateZ", "Clay_Silt",
                                "site_mean_temp", "site_precip", "hargreave_cmd")])
  f <- as.formula(paste0(var, " ~ site_mean_temp * Clay_Silt + site_precip * Clay_Silt"))
  m <- lm(f, data = df[keep,])
  pred <- predict(m)
  out <- rep(NA_real_, nrow(df))
  out[keep] <- resid(m) + mean(pred, na.rm = TRUE)
  out
}

############################################################
# End utils.R
############################################################
