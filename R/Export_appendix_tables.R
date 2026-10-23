################################################################################
# Export_appendix_tables.R
# Writes two CSV files from the US (noMX) soil-health model for the report
# appendix:
#   Outputs/adjustment_slopes_noMX.csv   climate/texture adjustment per indicator
#   Outputs/sem_noMX_standardized.csv    SEM standardized loadings and covariances
################################################################################
library(mgcv); library(lavaan)
proj <- readRDS(here::here("Outputs", "sem_noMX_projection.rds"))

# Climate/texture adjustment: slopes for MAT (per degree C) and MAP (per mm),
# evaluated at the average clay + silt, for each indicator
adj <- do.call(rbind, lapply(proj$indicators, function(v) {
  m <- proj$models[[v]]; d <- m$model; b <- coef(m); pt <- summary(m)$p.table
  int <- function(x) grep(paste0("^(", x, ":Clay_Silt|Clay_Silt:", x, ")$"), names(b), value = TRUE)
  cs <- mean(d$Clay_Silt)
  data.frame(indicator = v, n = nrow(d), r2 = summary(m)$r.sq,
             mean_indicator = mean(d[[v]]), sd_indicator = sd(d[[v]]), mean_clay_silt = cs,
             slope_MAT = b[["site_mean_temp"]] + b[[int("site_mean_temp")]] * cs,
             slope_MAP = b[["site_precip"]] + b[[int("site_precip")]] * cs,
             p_MAT = pt["site_mean_temp", 4], p_MAP = pt["site_precip", 4],
             p_MAT_x_texture = pt[int("site_mean_temp"), 4], p_MAP_x_texture = pt[int("site_precip"), 4])
}))
write.csv(adj, here::here("Outputs", "adjustment_slopes_noMX.csv"), row.names = FALSE)

# SEM: standardized loadings, second-order paths and residual covariances
write.csv(standardizedSolution(proj$fit), here::here("Outputs", "sem_noMX_standardized.csv"), row.names = FALSE)
cat("Saved adjustment_slopes_noMX.csv and sem_noMX_standardized.csv in Outputs\n")
