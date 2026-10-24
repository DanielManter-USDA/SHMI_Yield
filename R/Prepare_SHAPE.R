################################################################################
# Prepare_SHAPE.R
#
# Writes the NAPESHM soil data in the SHAPE multi-indicator input format
# (Multi_Indicator_example.csv):
#
#   field_name, sample_num, MAPmm, MATC, ACE, KSSL_WAS, Resp, AC, SOC, Suborder, Texture
#
#   field_name  management unit (MGT_combo, e.g. NAPESHM_USAL01_137)
#   sample_num  1, 2, ... within the unit (one row per soil sample, or 1 when
#               LEVEL = "treatment")
#   MAPmm       mean annual precipitation (mm)      PRISM 30-year normals at the site
#   MATC        mean annual temperature (C)          (as SHAPE uses); see CLIMATE_SOURCE
#   ACE         a_ace     ACE protein (g/kg)
#   KSSL_WAS    a_stabi   wet aggregate stability, KSSL method (%)
#   Resp        a_resp    4-day respiration (mg CO2/g soil)
#   AC          a_act_c   active carbon, POXC (mg/kg)
#   SOC         b_soc     soil organic carbon (%)
#   Suborder    SHAPE taxonomy group S1-S5 of the site's dominant soil
#               (shape_group in Outputs/site_soils_ssurgo.csv, Get_site_ssurgo.R)
#   Texture     SHAPE texture group T1-T5 from the sample's measured sand,
#               silt and clay (b_sand, b_silt, b_clay) via the USDA texture class
#
# SHAPE groups (SHAPE app; Nunes et al. 2021)
#   T1 sand, loamy sand, sandy loam (< 8% clay)
#   T2 sandy loam (>= 8% clay), sandy clay loam, loam
#   T3 silt loam, silt
#   T4 sandy clay, clay loam, silty clay loam, silty clay, clay (< 60% clay)
#   T5 clay (>= 60% clay)
#   S1-S5: suborders as in Get_site_ssurgo.R (SHAPE_GROUPS)
#
# Inputs (R_NAPESHM project; paths relative to its root)
#   Outputs/semDF_prepped.rds          units kept after Prepare_SEM.R's exclusions
#   Inputs/napeshm/soildata.csv      soil data, one row per sample
#   Outputs/site_soils_ssurgo.csv      SHAPE taxonomy group per site
#   PRISM normals                      downloaded by the prism package (CLIMATE_SOURCE)
#   Outputs/site_climate_monthly.rds   monthly normals per site, used where PRISM has
#                                      no value; sites.csv if this file is missing
# Outputs
#   Outputs/SHAPE/SHAPE_multi_indicator_input.csv     the file for the SHAPE app
#   Outputs/SHAPE/SHAPE_multi_indicator_key.csv       field_name -> site, treatment,
#                                                     system, texture class, climate source, flags
#   Outputs/SHAPE/SHAPE_coordinates.csv               field_name, lat, long
################################################################################
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(here) })

IN  <- function(...) here::here("Inputs", ...)
OUT <- function(...) here::here("Outputs", ...)
dir.create(OUT("SHAPE"), showWarnings = FALSE, recursive = TRUE)

LEVEL      <- "sample"    # "sample": one row per soil sample; "treatment": unit means
SCOPE      <- "us"        # "us": US sites only (SHAPE is calibrated on US soils and
                          #   SSURGO has no Mexican soils); "all": keep Mexico (Suborder NA)
IN_SCOPE   <- FALSE       # TRUE: only the SHMI in-scope units (annual + fallow)
DROP_INCOMPLETE <- FALSE  # TRUE: drop rows missing any input (else kept, listed in the key)

# ---- units, climate, suborder -----------------------------------------------
semDF <- readRDS(OUT("semDF_prepped.rds"))
units <- semDF %>%
  transmute(MGT_combo, sitecode = as.character(sitecode), treatmentid = as.character(treatmentid),
            system, in_scope, site_precip, site_mean_temp)
if (IN_SCOPE) units <- filter(units, in_scope)
if (SCOPE == "us") units <- filter(units, !grepl("^MX", sitecode))

# Climate: SHAPE uses PRISM 30-year normals.
#   "prism"   download PRISM annual normals (prism package) and read them at each
#             site (PRISM covers the conterminous US only; other sites fall back
#             to the monthly normals below)
#   "monthly" the monthly normals used for SHMI (site_climate_monthly.rds)
# Either way the site coordinates are written to Outputs/SHAPE/SHAPE_coordinates.csv
# (field_name, lat, long) for looking PRISM values up by hand.
CLIMATE_SOURCE <- "prism"
PRISM_DIR <- OUT("SHAPE", "prism")          # download folder for the PRISM rasters
PRISM_RES <- "800m"                         # "800m" or "4km"

coords <- semDF %>%
  transmute(field_name = MGT_combo, sitecode = as.character(sitecode),
            lat = site_latitude, long = site_longitude) %>%
  filter(field_name %in% units$MGT_combo)
write.csv(select(coords, field_name, lat, long), OUT("SHAPE", "SHAPE_coordinates.csv"), row.names = FALSE)

prism_climate <- function(sites) {
  if (!requireNamespace("prism", quietly = TRUE) || !requireNamespace("terra", quietly = TRUE))
    stop("Install the prism and terra packages: install.packages(c('prism', 'terra'))")
  dir.create(PRISM_DIR, showWarnings = FALSE, recursive = TRUE)
  prism::prism_set_dl_dir(PRISM_DIR)
  # the annual-normal raster for a variable, found by file name: prism versions
  # name the downloads differently (".._30yr_normal_800mM4_annual_bil.bil" in
  # older ones; newer ones use a different pattern and .tif), so search the folder
  # rather than rely on prism's own lookup
  find_annual <- function(type) {
    f <- list.files(PRISM_DIR, pattern = "\\.(bil|tif)$", recursive = TRUE, full.names = TRUE)
    f <- f[grepl(paste0("(^|_)", type, "(_|$)"), tolower(basename(f)))]
    ann <- f[grepl("annual", tolower(basename(f)))]
    if (!length(ann))                      # no "annual" in the names: drop monthly files
      ann <- f[!grepl("(_|20[0-9]{2})(0[1-9]|1[0-2])(_|\\.)", basename(f))]
    if (length(ann) > 1) message("Several PRISM ", type, " files; using ", basename(ann[1]), " (of ",
                                 paste(basename(ann), collapse = ", "), ")")
    ann[1]
  }
  value_at <- function(type) {
    file <- find_annual(type)
    if (is.na(file)) {
      prism::get_prism_normals(type = type, resolution = PRISM_RES, annual = TRUE, keepZip = FALSE)
      file <- find_annual(type)
    }
    if (is.na(file)) stop("No PRISM ", type, " annual normals found in ", PRISM_DIR, ". Files there:\n",
                          paste(list.files(PRISM_DIR, recursive = TRUE), collapse = "\n"))
    message("PRISM ", type, ": ", file)
    r <- terra::rast(file)
    pts <- terra::vect(sites, geom = c("long", "lat"), crs = "EPSG:4326")
    terra::extract(r, terra::project(pts, terra::crs(r)))[[2]]
  }
  sites %>% mutate(MAPmm = value_at("ppt"), MATC = value_at("tmean")) %>% select(sitecode, MATC, MAPmm)
}

CLIMATE <- OUT("site_climate_monthly.rds")
clim <- if (file.exists(CLIMATE)) {
  readRDS(CLIMATE) %>%
    transmute(sitecode = as.character(sitecode),
              MATC  = rowMeans(across(sprintf("tavg_%02d", 1:12))),
              MAPmm = rowSums(across(sprintf("prec_%02d", 1:12))))
} else {
  message("No ", CLIMATE, ": using site_mean_temp and site_precip from sites.csv (whole degrees)")
  units %>% distinct(sitecode, MATC = site_mean_temp, MAPmm = site_precip)
}
clim$climate_source <- "monthly normals (SHMI)"
if (CLIMATE_SOURCE == "prism") {
  pr <- prism_climate(distinct(coords, sitecode, lat, long))
  ok <- !is.na(pr$MATC) & !is.na(pr$MAPmm)
  if (any(!ok)) message("No PRISM value (outside the conterminous US?) for: ",
                        paste(pr$sitecode[!ok], collapse = ", "), "; monthly normals kept there")
  clim <- clim %>%
    full_join(rename(pr[ok, ], MATC_p = MATC, MAPmm_p = MAPmm), by = "sitecode") %>%
    mutate(climate_source = ifelse(is.na(MATC_p), climate_source, paste("PRISM normals", PRISM_RES)),
           MATC = coalesce(MATC_p, MATC), MAPmm = coalesce(MAPmm_p, MAPmm)) %>%
    select(-MATC_p, -MAPmm_p)
}

SOILS <- OUT("site_soils_ssurgo.csv")
if (!file.exists(SOILS)) stop("No ", SOILS, ": run Get_site_ssurgo.R first (it assigns each site's SHAPE group)")
suborder <- read.csv(SOILS, stringsAsFactors = FALSE) %>%
  transmute(sitecode = as.character(sitecode), Suborder = shape_group, taxsuborder, ssurgo_flag = flag)

# ---- soil samples -----------------------------------------------------------
soil <- read.csv(IN("napeshm", "soildata.csv"), stringsAsFactors = FALSE) %>%
  mutate(MGT_combo = paste0("NAPESHM_", sitecode, "_", treatmentid)) %>%
  filter(MGT_combo %in% units$MGT_combo)
need <- c("a_ace", "a_stabi", "a_resp", "a_act_c", "b_soc", "b_sand", "b_silt", "b_clay")
miss <- setdiff(need, names(soil))
if (length(miss)) stop("soildata.csv lacks: ", paste(miss, collapse = ", "))

# order samples within a unit by whatever sample identifiers the file has
id_cols <- intersect(c("plotid", "plot_id", "plot", "rep", "replicate", "block", "sampleid", "sample_id",
                       "sample", "depth_top", "sampling_date", "date"), names(soil))
soil <- soil %>% arrange(MGT_combo, across(all_of(id_cols)))

if (LEVEL == "treatment") {
  soil <- soil %>% group_by(MGT_combo) %>%
    summarise(across(all_of(need), ~ mean(.x, na.rm = TRUE)), n_samples = n(), .groups = "drop") %>%
    mutate(across(all_of(need), ~ ifelse(is.nan(.x), NA_real_, .x)))
} else {
  soil <- soil %>% mutate(n_samples = 1L)
}

# ---- texture: USDA class from sand / silt / clay, then the SHAPE group --------
usda_texture <- function(sand, silt, clay) {
  tot <- sand + silt + clay                       # rescale to 100 when sums drift (lab rounding)
  sand <- 100 * sand / tot; silt <- 100 * silt / tot; clay <- 100 * clay / tot
  case_when(
    is.na(sand) | is.na(silt) | is.na(clay) ~ NA_character_,
    clay >= 40 & silt >= 40 ~ "silty clay",
    clay >= 40 & sand > 45  ~ "sandy clay",
    clay >= 40              ~ "clay",
    clay >= 35 & sand > 45  ~ "sandy clay",
    clay >= 27 & sand <= 20 ~ "silty clay loam",
    clay >= 27 & sand <= 45 ~ "clay loam",
    clay >= 27              ~ "sandy clay loam",
    clay >= 20 & sand > 45 & silt < 28 ~ "sandy clay loam",
    silt >= 80 & clay < 12  ~ "silt",
    silt >= 50              ~ "silt loam",
    clay >= 7 & silt >= 28 & sand <= 52 ~ "loam",
    silt + 1.5 * clay < 15  ~ "sand",
    silt + 2 * clay < 30    ~ "loamy sand",
    TRUE                    ~ "sandy loam")
}
shape_texture <- function(cls, clay) {
  case_when(
    is.na(cls) ~ NA_character_,
    cls %in% c("sand", "loamy sand") ~ "T1",
    cls == "sandy loam" & clay < 8   ~ "T1",
    cls %in% c("sandy loam", "sandy clay loam", "loam") ~ "T2",
    cls %in% c("silt loam", "silt") ~ "T3",
    cls == "clay" & clay >= 60 ~ "T5",
    TRUE ~ "T4")
}

# ---- assemble ---------------------------------------------------------------
shape <- soil %>%
  select(MGT_combo, all_of(need), n_samples) %>%
  left_join(units, by = "MGT_combo") %>%
  left_join(clim, by = "sitecode") %>%
  left_join(suborder, by = "sitecode") %>%
  mutate(texture_class = usda_texture(b_sand, b_silt, b_clay),
         Texture = shape_texture(texture_class, b_clay)) %>%
  group_by(MGT_combo) %>% mutate(sample_num = row_number()) %>% ungroup() %>%
  transmute(field_name = MGT_combo, sample_num, MAPmm, MATC,
            ACE = a_ace, KSSL_WAS = a_stabi, Resp = a_resp, AC = a_act_c, SOC = b_soc,
            Suborder, Texture,
            # key columns (not in the SHAPE file)
            sitecode, treatmentid, system, in_scope, n_samples, sand = b_sand, silt = b_silt, clay = b_clay,
            texture_class, taxsuborder, ssurgo_flag, climate_source)

SHAPE_COLS <- c("field_name", "sample_num", "MAPmm", "MATC", "ACE", "KSSL_WAS", "Resp", "AC", "SOC",
                "Suborder", "Texture")
shape <- shape %>%
  mutate(missing = apply(is.na(across(all_of(SHAPE_COLS))), 1, function(m) paste(SHAPE_COLS[m], collapse = "; ")),
         note = case_when(Suborder == "S1" ~ "S1 (organic soil): not covered by the SHAPE model",
                          nzchar(coalesce(ssurgo_flag, "")) ~ ssurgo_flag,
                          TRUE ~ ""))
if (DROP_INCOMPLETE) shape <- filter(shape, missing == "")

write.csv(shape[SHAPE_COLS], OUT("SHAPE", "SHAPE_multi_indicator_input.csv"), row.names = FALSE, na = "")
write.csv(shape, OUT("SHAPE", "SHAPE_multi_indicator_key.csv"), row.names = FALSE, na = "")

# ---- checks -----------------------------------------------------------------
cat("\nRows:", nrow(shape), "| units:", n_distinct(shape$field_name), "| sites:", n_distinct(shape$sitecode), "\n")
cat("\nSHAPE texture groups:\n");  print(table(shape$Texture, useNA = "ifany"))
cat("\nSHAPE taxonomy groups:\n"); print(table(shape$Suborder, useNA = "ifany"))
cat("\nRows missing an input:", sum(shape$missing != ""), "\n")
print(count(filter(shape, missing != ""), missing, sort = TRUE))
cat("\nIndicator ranges (check units against the template):\n")
print(summary(shape[c("MAPmm", "MATC", "ACE", "KSSL_WAS", "Resp", "AC", "SOC")]))
# texture groups that differ between samples of one site (expected for some sites;
# SHAPE scores each sample against its own texture group)
mixed <- shape %>% group_by(sitecode) %>% summarise(groups = paste(sort(unique(na.omit(Texture))), collapse = "/"),
                                                    .groups = "drop") %>% filter(grepl("/", groups))
if (nrow(mixed)) { cat("\nSites with samples in more than one texture group:\n"); print(mixed, n = Inf) }
message("Saved ", OUT("SHAPE", "SHAPE_multi_indicator_input.csv"))
