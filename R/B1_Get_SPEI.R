################################################################################
# B1_Get_SPEI.R
# Fill Outputs/yield/drought_years.csv with two growing-season SPEI values per
# site-year from SPEIbase (0.5 deg, monthly, 1901 onward), both from the 3-month
# file (spei03.nc), i.e. the standardized water balance (precipitation minus
# potential evapotranspiration) accumulated over the month and the two before it:
#   spei_mjj  SPEI-3 for July   = May-July   (to silking)
#   spei_jja  SPEI-3 for August = June-August (silking and early grain fill; the
#             window most strongly related to US corn yield in the literature)
#
# 1. Download spei03.nc from the latest SPEIbase release on Digital.CSIC, linked
#    from https://spei.csic.es/database.html
# 2. Put it in Inputs/climate/ (or change SPEI_FILE) and run this script.
#    drought_years.csv is created by B3_SHMI_yield_stability.R on its first run
#    (one row per site-year); this script fills both columns, overwriting
#    earlier values, and adds every year since 1950 for each site (so sites added
#    later, e.g. the AgCROS sites, are covered without rerunning the stability script).
# 3. INCLUDE_AGCROS = TRUE adds the AgCROS sites converted by agcros2shmi.py
#    (same switch as B2_SHMI_yield_main.R).
################################################################################
library(terra)
library(dplyr)
library(readxl)

setwd(here::here())                 # the project root (SHMI_Yield.Rproj)
SPEI_FILE <- "Inputs/climate/spei03.nc"
DFILE     <- "Outputs/yield/drought_years.csv"
WINDOWS   <- c(spei_mjj = "07", spei_jja = "08")     # column = SPEI-3 for that month

workbooks <- c(NELITCSE  = "Inputs/yield_workbooks/SHMI_NELITCSE.xlsx",
               NEMLTCRS  = "Inputs/yield_workbooks/SHMI_NEMLTCRS.xlsx",
               OHTVD_Ho  = "Inputs/yield_workbooks/SHMI_OHTVD_Ho.xlsx",
               OHTVD_Wo  = "Inputs/yield_workbooks/SHMI_OHTVD_Wo.xlsx",
               WICST_ARS = "Inputs/yield_workbooks/SHMI_WICST_ARS.xlsx",
               NEMEIRR   = "Inputs/yield_workbooks/SHMI_NEMEIRR.xlsx",
               NECCIRR   = "Inputs/yield_workbooks/SHMI_NECCIRR.xlsx")
INCLUDE_AGCROS <- TRUE
if (exists("INCLUDE_AGCROS_OVERRIDE")) INCLUDE_AGCROS <- INCLUDE_AGCROS_OVERRIDE
AGCROS_SITES <- c("INWLACRE", "IAAM7071", "PAUP")   # MNMOFS, INWLTPAC screened out (see B2_SHMI_yield_main.R)
if (exists("AGCROS_SITES_OVERRIDE")) AGCROS_SITES <- AGCROS_SITES_OVERRIDE
if (INCLUDE_AGCROS)
  workbooks <- c(workbooks, setNames(file.path("Inputs/yield_workbooks", paste0("SHMI_", AGCROS_SITES, ".xlsx")), AGCROS_SITES))

# site coordinates: mean of the plot latitudes/longitudes on each Mgt_Unit sheet
coords <- bind_rows(lapply(names(workbooks), function(s) {
  mu  <- read_excel(workbooks[[s]], sheet = "Mgt_Unit")
  lat <- grep("lat", names(mu), ignore.case = TRUE, value = TRUE)[1]
  lon <- grep("lon", names(mu), ignore.case = TRUE, value = TRUE)[1]
  tibble(site = s, lat = mean(as.numeric(mu[[lat]]), na.rm = TRUE),
         lon = mean(as.numeric(mu[[lon]]), na.rm = TRUE))
}))
print(coords)

r  <- rast(SPEI_FILE)
tt <- time(r)
if (is.null(tt) || all(is.na(tt)))                       # one layer per month from Jan 1901
  tt <- seq(as.Date("1901-01-15"), by = "month", length.out = nlyr(r))
pts <- vect(coords, geom = c("lon", "lat"), crs = "EPSG:4326")

extract_month <- function(mm, col) {
  k <- which(format(as.Date(tt), "%m") == mm)
  v <- terra::extract(r[[k]], pts, ID = FALSE)
  bind_rows(lapply(seq_len(nrow(coords)), function(i)
    tibble(site = coords$site[i], year = as.integer(format(as.Date(tt[k]), "%Y")),
           !!col := as.numeric(v[i, ]))))
}
spei <- Reduce(function(a, b) full_join(a, b, by = c("site", "year")),
               Map(extract_month, WINDOWS, names(WINDOWS)))

old <- if (file.exists(DFILE)) read.csv(DFILE) %>% select(site, year) else tibble(site = character(), year = integer())
d <- bind_rows(old, spei %>% filter(year >= 1950) %>% select(site, year)) %>% distinct() %>%
  left_join(spei, by = c("site", "year")) %>% arrange(site, year)
write.csv(d, DFILE, row.names = FALSE)

# drought years (SPEI <= -1, moderate or worse) per site and window, all years in the
# file (the stability script uses only the years with yields); the drought-year test
# needs at least 3 per site
d %>% group_by(site) %>%
  summarise(years = n(),
            missing_mjj = sum(is.na(spei_mjj)), drought_mjj = sum(spei_mjj <= -1, na.rm = TRUE),
            missing_jja = sum(is.na(spei_jja)), drought_jja = sum(spei_jja <= -1, na.rm = TRUE)) %>%
  print()
cat("May-July vs June-August SPEI, correlation across site-years:",
    round(cor(d$spei_mjj, d$spei_jja, use = "complete.obs"), 2), "\n")
