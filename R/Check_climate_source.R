################################################################################
# Check_climate_source.R
#
# Would SHMI's Cover seasons change if WorldClim were replaced by Daymet?
# A month is a growing month when its long-term mean temperature is >= 5 C and it
# is not dry (precipitation in mm >= 2 x temperature in C; Bagnouls-Gaussen), the
# rule in the SHMI package (irrigation, which switches off the dry rule for a
# management unit, is ignored here: this is a site-level check).
#
# Compares, for each US NAPESHM site, the 12 month labels from
#   WorldClim 2.1 normals, 1970-2000 (what SHMI uses)       Outputs/site_climate_monthly.rds
#   Daymet normals, 1980-2000 (same era, different source)  -> the source effect
#   Daymet normals, 1991-2020 (recent era)                  -> source + recent warming
# and reports how many site-months change label, which months, and which sites.
# Also reported for the temperature rule alone (>= 5 C), to show whether changes
# come from temperature or from the dry-month rule.
#
# Needs: daymetr. Daymet daily data are downloaded once per site and cached in
# Outputs/daymet_site_monthly.csv (Daymet starts in 1980).
# Outputs: Outputs/climate_source_check.csv (every site-month), printed summaries
################################################################################
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(here) })

GROW_T <- 5
CACHE_FILE <- here::here("Outputs", "daymet_site_monthly.csv")

sites <- read.csv(here::here("Inputs", "sites.csv")) %>%
  filter(grepl("^US", sitecode)) %>%
  select(sitecode, lon = site_longitude, lat = site_latitude)
wc <- readRDS(here::here("Outputs", "site_climate_monthly.rds")) %>% filter(sitecode %in% sites$sitecode)
cat("US sites:", nrow(sites), "\n")

# ---- Daymet monthly values per site and year (cached) -------------------------------
if (!file.exists(CACHE_FILE)) {
  if (!requireNamespace("daymetr", quietly = TRUE)) stop("install.packages('daymetr')")
  dm <- bind_rows(lapply(seq_len(nrow(sites)), function(i) {
    s <- sites[i, ]
    x <- tryCatch(daymetr::download_daymet(site = s$sitecode, lat = s$lat, lon = s$lon,
                                           start = 1980, end = 2020, internal = TRUE, silent = TRUE)$data,
                  error = function(e) { message("Daymet failed for ", s$sitecode, ": ", conditionMessage(e)); NULL })
    if (is.null(x)) return(NULL)
    x %>% mutate(date = as.Date(paste(year, yday), "%Y %j"), month = as.integer(format(date, "%m")),
                 tavg = (tmax..deg.c. + tmin..deg.c.) / 2) %>%
      group_by(year, month) %>% summarise(tavg = mean(tavg), prec = sum(prcp..mm.day.), .groups = "drop") %>%
      mutate(sitecode = s$sitecode)
  }))
  write.csv(dm, CACHE_FILE, row.names = FALSE)
}
dm <- read.csv(CACHE_FILE)

normals <- function(y0, y1, label) dm %>% filter(year >= y0, year <= y1) %>%
  group_by(sitecode, month) %>% summarise(tavg = mean(tavg), prec = mean(prec), .groups = "drop") %>%
  mutate(source = label)
wc_long <- wc %>% pivot_longer(-sitecode, names_to = c(".value", "month"), names_sep = "_") %>%
  mutate(month = as.integer(month), source = "WorldClim 1970-2000")
all <- bind_rows(wc_long, normals(1980, 2000, "Daymet 1980-2000"), normals(1991, 2020, "Daymet 1991-2020")) %>%
  mutate(grow_temp = tavg >= GROW_T,
         grow = grow_temp & prec >= 2 * tavg)

wide <- all %>% select(sitecode, month, source, tavg, prec, grow, grow_temp) %>%
  pivot_wider(names_from = source, values_from = c(tavg, prec, grow, grow_temp))
write.csv(wide, here::here("Outputs", "climate_source_check.csv"), row.names = FALSE)

# ---- summaries ----------------------------------------------------------------------
cmp <- function(ref, alt, rule = "grow") {
  a <- wide[[paste0(rule, "_", ref)]]; b <- wide[[paste0(rule, "_", alt)]]
  ok <- !is.na(a) & !is.na(b)
  tibble(sitecode = wide$sitecode[ok], month = wide$month[ok], ref = a[ok], alt = b[ok]) %>%
    mutate(change = case_when(ref == alt ~ "same", alt ~ "becomes growing", TRUE ~ "becomes non-growing"))
}
report <- function(alt) {
  for (rule in c("grow", "grow_temp")) {
    x <- cmp("WorldClim 1970-2000", alt, rule)
    cat(sprintf("\n==== WorldClim 1970-2000 vs %s (%s) ====\n", alt,
                if (rule == "grow") "temperature and dry-month rule" else "temperature rule only"))
    cat(sprintf("Site-months compared: %d; changed label: %d (%.1f%%); sites with any change: %d of %d\n",
                nrow(x), sum(x$change != "same"), 100 * mean(x$change != "same"),
                n_distinct(x$sitecode[x$change != "same"]), n_distinct(x$sitecode)))
    if (any(x$change != "same")) {
      cat("By month:\n")
      print(x %>% filter(change != "same") %>% count(month = month.abb[month], change) %>%
              pivot_wider(names_from = change, values_from = n, values_fill = 0) %>% as.data.frame(), row.names = FALSE)
      if (rule == "grow") {
        cat("Growing months per site (WorldClim -> alternative), sites that change:\n")
        print(x %>% group_by(sitecode) %>% summarise(WorldClim = sum(ref), alternative = sum(alt), .groups = "drop") %>%
                filter(WorldClim != alternative) %>% arrange(desc(abs(alternative - WorldClim))) %>% as.data.frame(), row.names = FALSE)
      }
    }
  }
}
report("Daymet 1980-2000")
report("Daymet 1991-2020")

# temperature agreement, for context
cat("\nMonthly mean temperature, Daymet minus WorldClim (C), across site-months:\n")
for (alt in c("Daymet 1980-2000", "Daymet 1991-2020")) {
  dlt <- wide[[paste0("tavg_", alt)]] - wide[["tavg_WorldClim 1970-2000"]]
  cat(sprintf("  %s: mean %+.2f, 95%% of differences within %+.2f to %+.2f\n", alt,
              mean(dlt, na.rm = TRUE), quantile(dlt, 0.025, na.rm = TRUE), quantile(dlt, 0.975, na.rm = TRUE)))
}
cat("\nSaved Outputs/climate_source_check.csv\n")
