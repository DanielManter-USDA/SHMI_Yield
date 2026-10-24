# ================================================================
# Get_site_ssurgo.R — soil series and taxonomy for each NAPESHM site
#   For each US site (sites.csv coordinates), finds the SSURGO map unit
#   at the point (or the map unit covering most of a small buffer around
#   it), takes the map unit's dominant component, and returns its series,
#   taxonomy, drainage class and surface-horizon properties.
#   Output: Outputs/site_soils_ssurgo.csv (+ .rds)
# ================================================================
suppressPackageStartupMessages({
  library(sf); library(dplyr); library(DBI); library(RSQLite); library(here)
})

SSURGO_GPKG <- Sys.getenv("GSSURGO_GPKG", here::here("Inputs", "gSSURGO_CONUS.gpkg"))     # <- path to the national gpkg
SITES_CSV   <- here::here("Inputs", "napeshm", "sites.csv")
OUT_CSV     <- here::here("Outputs", "site_soils_ssurgo.csv")
BUFFER_M    <- 100    # 0 = map unit at the exact point; > 0 = map unit covering most of a buffer

# ---- helpers ----------------------------------------------------------------
# SELECT ... WHERE key IN (...), chunked to stay under SQLite limits
sql_in <- function(con, table, key, values, chunk = 5000) {
  values <- unique(na.omit(as.character(values)))
  if (length(values) == 0) return(data.frame())
  chunks <- split(values, ceiling(seq_along(values) / chunk))
  out <- bind_rows(lapply(chunks, function(v) {
    DBI::dbGetQuery(con, sprintf("SELECT * FROM %s WHERE %s IN (%s)",
                                 table, key, paste0("'", v, "'", collapse = ",")))
  }))
  names(out) <- gsub("_(r|l|h)$", ".\\1", tolower(names(out)))
  out
}
col_or_na <- function(d, nm) if (nm %in% names(d)) d[[nm]] else NA

# SHAPE taxonomy groups (Nunes et al. 2021): USDA suborders grouped by inherent
# capacity to store soil organic carbon, from organic soils (S1) to dry,
# sandy or strongly weathered soils (S5)
SHAPE_GROUPS <- list(
  S1 = c("Fibrists", "Folists", "Hemists", "Histels", "Saprists", "Wassists"),
  S2 = c("Aquands", "Aquents", "Aquepts", "Aquods", "Aquox", "Cryods", "Humods", "Orthels", "Perox",
         "Torrands", "Tropepts", "Turbels", "Udands", "Udox", "Ustands"),
  S3 = c("Albolls", "Andepts", "Aquolls", "Aquults", "Cryands", "Cryepts", "Cryolls", "Gelepts", "Gelolls",
         "Humults", "Rendolls", "Umbrepts", "Ustox", "Vitrands", "Wassents", "Xerands"),
  S4 = c("Aqualfs", "Aquerts", "Boralfs", "Borolls", "Cryalfs", "Ochrepts", "Orthods", "Orthox", "Udalfs",
         "Udepts", "Uderts", "Udolls", "Usterts", "Ustolls", "Xeralfs", "Xerepts", "Xerolls", "Xerults"),
  S5 = c("Arents", "Argids", "Calcids", "Cambids", "Cryerts", "Cryids", "Durids", "Fluvents", "Gypsids",
         "Orthents", "Orthids", "Psamments", "Salids", "Torrerts", "Torrox", "Udults", "Ustalfs", "Ustepts",
         "Ustults", "Xererts")
)
shape_lookup <- setNames(rep(names(SHAPE_GROUPS), lengths(SHAPE_GROUPS)), tolower(unlist(SHAPE_GROUPS)))
# SSURGO spells oxisol suborders "Aquox", "Udox" etc.; accept the "-oxs" spelling too
shape_group_of <- function(suborder) {
  k <- tolower(trimws(suborder)); k <- sub("oxs$", "ox", k)
  unname(shape_lookup[k])
}

# ---- sites (US only; SSURGO covers the US) ----------------------------------
sites <- read.csv(SITES_CSV, stringsAsFactors = FALSE) %>%
  filter(!grepl("^MX", sitecode), !is.na(site_latitude), !is.na(site_longitude)) %>%
  transmute(sitecode, lat = site_latitude, lon = site_longitude)
pts <- st_as_sf(sites, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
message("US sites with coordinates: ", nrow(pts))

# ---- map unit at (or around) each site ---------------------------------------
lyrs <- suppressWarnings(st_layers(SSURGO_GPKG))
i    <- match("mupolygon", tolower(lyrs$name))
if (is.na(i)) stop("No MUPOLYGON layer in ", SSURGO_GPKG)
mu_layer <- lyrs$name[i]
pts_p    <- st_transform(pts, lyrs$crs[[i]])

mukey_for <- function(k) {
  pt   <- pts_p[k, ]
  area <- st_buffer(st_geometry(pt), max(BUFFER_M, 1))
  poly <- st_read(SSURGO_GPKG, layer = mu_layer, wkt_filter = st_as_text(area), quiet = TRUE)
  if (!nrow(poly)) return(tibble(sitecode = pt$sitecode, mukey = NA_character_, share = NA_real_))
  poly <- rename_with(poly, tolower, .cols = -all_of(attr(poly, "sf_column"))) %>%
    mutate(mukey = as.character(mukey)) %>% st_make_valid()
  if (BUFFER_M <= 0) {
    hit <- poly[lengths(st_intersects(poly, pt)) > 0, ]
    if (!nrow(hit)) hit <- poly[1, ]
    return(tibble(sitecode = pt$sitecode, mukey = hit$mukey[1], share = 1))
  }
  ov <- suppressWarnings(st_intersection(poly, area))
  ov$a <- as.numeric(st_area(ov))
  best <- ov %>% st_drop_geometry() %>% group_by(mukey) %>% summarise(a = sum(a), .groups = "drop") %>%
    mutate(share = a / sum(a)) %>% slice_max(a, n = 1, with_ties = FALSE)
  tibble(sitecode = pt$sitecode, mukey = best$mukey, share = best$share)
}
mu <- bind_rows(lapply(seq_len(nrow(pts_p)), mukey_for))
message("Sites matched to a map unit: ", sum(!is.na(mu$mukey)), " of ", nrow(mu))

# ---- tabular: map unit, dominant component, surface horizon -------------------
con <- DBI::dbConnect(RSQLite::SQLite(), SSURGO_GPKG)
on.exit(DBI::dbDisconnect(con), add = TRUE)
mapunit   <- sql_in(con, "mapunit",   "mukey", mu$mukey)
component <- sql_in(con, "component", "mukey", mu$mukey)
dom <- component %>% mutate(comppct.r = as.numeric(comppct.r)) %>%
  group_by(mukey) %>% slice_max(comppct.r, n = 1, with_ties = FALSE) %>% ungroup()
chorizon <- sql_in(con, "chorizon", "cokey", dom$cokey)
surf <- chorizon %>% mutate(hzdept.r = as.numeric(hzdept.r)) %>%
  group_by(cokey) %>% slice_min(hzdept.r, n = 1, with_ties = FALSE) %>% ungroup()

soils <- mu %>%
  left_join(mapunit %>% transmute(mukey = as.character(mukey), muname = col_or_na(mapunit, "muname")), by = "mukey") %>%
  left_join(dom %>% transmute(mukey = as.character(mukey), cokey,
                              series = col_or_na(dom, "compname"), comppct = comppct.r,
                              taxorder = col_or_na(dom, "taxorder"), taxsuborder = col_or_na(dom, "taxsuborder"),
                              taxgrtgroup = col_or_na(dom, "taxgrtgroup"), taxsubgrp = col_or_na(dom, "taxsubgrp"),
                              taxpartsize = col_or_na(dom, "taxpartsize"), drainagecl = col_or_na(dom, "drainagecl"),
                              hydgrp = col_or_na(dom, "hydgrp")), by = "mukey") %>%
  left_join(surf %>% transmute(cokey,
                               surf_om = as.numeric(col_or_na(surf, "om.r")),
                               surf_clay = as.numeric(col_or_na(surf, "claytotal.r")),
                               surf_sand = as.numeric(col_or_na(surf, "sandtotal.r")),
                               surf_ph = as.numeric(col_or_na(surf, "ph1to1h2o.r")),
                               surf_cec = as.numeric(col_or_na(surf, "cec7.r"))), by = "cokey") %>%
  mutate(shape_group = shape_group_of(taxsuborder),
         flag = case_when(is.na(mukey) ~ "no map unit found",
                          is.na(taxsuborder) | grepl("water|urban|pits", tolower(paste(muname, series))) ~ "check: no taxonomy (water, urban or miscellaneous area)",
                          share < 0.5 ~ "check: dominant map unit covers < 50% of the buffer",
                          is.na(shape_group) ~ "check: suborder not in the SHAPE groups",
                          TRUE ~ "")) %>%
  select(-cokey)

write.csv(soils, OUT_CSV, row.names = FALSE)
saveRDS(soils, sub("\\.csv$", ".rds", OUT_CSV))
message("Saved → ", OUT_CSV)
print(count(soils, shape_group, sort = TRUE))
print(count(soils, shape_group, taxsuborder, sort = TRUE), n = 40)
print(count(soils, drainagecl, sort = TRUE))
