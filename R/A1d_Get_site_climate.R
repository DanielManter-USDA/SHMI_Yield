################################################################################
# A1d_Get_site_climate.R
#
# Monthly climate normals (WorldClim 2.1, 1970-2000) for each NAPESHM site,
# used to define climate-based Cover seasons: a month is ACTIVE at a site when
# its mean temperature is at or above a threshold (5 C by default; 10 C as a
# variant), DORMANT otherwise.
#
# Needs: geodata, terra (install.packages(c("geodata", "terra")))
# Downloads the global 10-arc-minute WorldClim layers once (~ tens of MB) to
# Inputs/worldclim/.
# Output: Outputs/site_climate_monthly.rds  (sitecode, tavg_01..12, prec_01..12)
################################################################################

library(dplyr)
library(tibble)
library(here)

RES     <- 10                        # arc-minutes: 10, 5, 2.5 (finer = larger download)
WC_DIR  <- here::here("Inputs", "worldclim")
dir.create(WC_DIR, showWarnings = FALSE, recursive = TRUE)

sites <- read.csv(here::here("Inputs", "napeshm", "sites.csv")) %>%
  select(sitecode, site_longitude, site_latitude, site_mean_temp)

tavg <- geodata::worldclim_global(var = "tavg", res = RES, path = WC_DIR)
prec <- geodata::worldclim_global(var = "prec", res = RES, path = WC_DIR)

pts <- cbind(sites$site_longitude, sites$site_latitude)
t_m <- terra::extract(tavg, pts)                 # one column per month (plus ID)
p_m <- terra::extract(prec, pts)
t_m <- as.matrix(t_m[, grep("tavg", names(t_m))])
p_m <- as.matrix(p_m[, grep("prec", names(p_m))])
stopifnot(ncol(t_m) == 12, ncol(p_m) == 12)

clim <- tibble(sitecode = sites$sitecode)
clim[sprintf("tavg_%02d", 1:12)] <- as.data.frame(round(t_m, 2))
clim[sprintf("prec_%02d", 1:12)] <- as.data.frame(round(p_m, 1))

# Sites on a grid cell without data (e.g. coastal): report; fix coordinates or
# fill from a neighbouring cell before use
miss <- clim$sitecode[rowSums(is.na(t_m)) > 0]
if (length(miss)) warning("No WorldClim data at: ", paste(miss, collapse = ", "))

# Check: WorldClim annual mean vs the site mean temperature in sites.csv
chk <- tibble(sitecode = sites$sitecode, site_mean_temp = sites$site_mean_temp,
              wc_annual_mean = round(rowMeans(t_m), 1))
cat("Correlation of annual mean temperature, sites.csv vs WorldClim:",
    round(cor(chk$site_mean_temp, chk$wc_annual_mean, use = "complete.obs"), 3), "\n")

# Active months per site at 5 and 10 C
act <- tibble(sitecode = sites$sitecode,
              active_5C  = rowSums(t_m >= 5,  na.rm = TRUE),
              active_10C = rowSums(t_m >= 10, na.rm = TRUE),
              # temperature and moisture: >= 5 C and not dry (P >= 2 x T, Bagnouls-Gaussen)
              active_5C_moist = rowSums(t_m >= 5 & p_m >= 2 * t_m, na.rm = TRUE),
              dry_months      = rowSums(p_m < 2 * t_m, na.rm = TRUE))
cat("\nNumber of active months per site:\n")
print(table(`5 C` = act$active_5C))
print(table(`10 C` = act$active_10C))
print(table(`5 C and not dry` = act$active_5C_moist))
cat("\nSites with at least one dry month (P < 2T):\n")
print(act %>% filter(dry_months > 0) %>% arrange(desc(dry_months)) %>% select(sitecode, dry_months, active_5C, active_5C_moist), n = Inf)
cat("\nSites with no dormant month at 5 C:",
    paste(act$sitecode[act$active_5C == 12], collapse = ", "), "\n")

saveRDS(clim, here::here("Outputs", "site_climate_monthly.rds"))
write.csv(left_join(clim, act, by = "sitecode"),
          here::here("Outputs", "site_climate_monthly.csv"), row.names = FALSE)
cat("\nSaved: Outputs/site_climate_monthly.rds and .csv\n")
