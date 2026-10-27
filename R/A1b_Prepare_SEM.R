################################################################################
# A1b_Prepare_SEM.R
#
# Purpose
#   Build the data set used to model soil health (SH) with a structural
#   equation model, and classify every management unit by cropping system so
#   the SHMI calibration population is defined once, reproducibly.
#
# Project layout (RStudio project R_NAPESHM; paths are relative to its root)
#   Inputs/   NAPESHM_test.xlsx   SHMI management workbook (corrected version)
#             napeshm/sites.csv   site climate and descriptors
#             napeshm/            treatments.csv, soildata.csv, plfa_ng.csv
#   Outputs/  everything this script writes
#   R/        scripts
#
# Steps
#   1. Read the SHMI workbook (prepare_shmi_inputs) and export data checks.
#   2. Classify each unit's cropping system from its management records:
#        annual            at least one annual crop in the evaluation window
#        fallow_reference  fallow only (bare-fallow reference treatments)
#        perennial         perennial crops only (pasture, rangeland, hay, ...)
#      In scope for SHMI 1.2.0: annual + fallow_reference.
#   3. Read site, treatment, soil and PLFA data; apply exclusions.
#   4. Build the CASH indicator variables used by the SEM.
#   5. Attach the system classification and save.
#
# Firewall
#   The SEM must not see management intensity. Only the system label and the
#   in-scope flag are attached here; they define WHICH units are modelled,
#   never HOW soil health is measured. Cover, Diversity, InvDist and OrgInput
#   are not added to the SEM data.
#
# Outputs (Outputs/)
#   shmi_inputs_v120.rds             SHMI inputs from the workbook
#   scope_v120.rds                   classification of every unit + metadata
#   in_scope_v120.rds                character vector of in-scope MGT_combo
#   data_checks_for_review.csv       records flagged by prepare_shmi_inputs()
#   semDF_prepped.rds / .csv         SEM data, ALL units, with `system` and
#                                    `in_scope`; the SEM script filters
################################################################################

library(dplyr)
library(tidyr)
library(here)
library(SHMI)

stopifnot(packageVersion("SHMI") >= "1.2.0")

# ==============================================================================
# 0. Settings
# ==============================================================================
IN  <- function(...) here::here("Inputs", ...)
OUT <- function(...) here::here("Outputs", ...)
dir.create(OUT(), showWarnings = FALSE)

WORKBOOK <- IN("NAPESHM_rebuilt.xlsx")   # rebuilt from the raw NAPESHM files

# Treatments removed from all analyses, with the reason for each group.
# Treatments without management records are reference soils (sod, remnant
# prairie, exclosures, benchmarks) and are kept: the rebuilt workbook enters
# them as perennial stands (conversion rule R2e).
EXCLUDE <- list(
  woody   = c("NAPESHM_USMI04_306", "NAPESHM_USNC03_368", "NAPESHM_USNC03_369",
              "NAPESHM_USNC03_370", "NAPESHM_USNC03_371", "NAPESHM_USNC03_372"),
  grazing = c("NAPESHM_USAZ01_167", "NAPESHM_USAZ01_168",
              "NAPESHM_USAZ03_169", "NAPESHM_USAZ03_170",
              "NAPESHM_USCO05_222", "NAPESHM_USCO05_223", "NAPESHM_USCO05_224",
              "NAPESHM_USSD03_465", "NAPESHM_USSD03_466", "NAPESHM_USSD03_467",
              "NAPESHM_USWY01_506", "NAPESHM_USWY01_507")
)
EXCLUDE_ALL <- unlist(EXCLUDE, use.names = FALSE)

EXCLUDE_SITES <- c("USWA03")        # in sites.csv but not in the treatment data

TRT_REF_YEAR <- 2019                # treatment age = REF_YEAR - start year + 1
TRT_AGE_CAP  <- 50                  # ages above this are set to the cap

# Build MGT_combo from site code and treatment id (same rule everywhere)
make_combo <- function(sitecode, treatmentid) paste0("NAPESHM_", sitecode, "_", treatmentid)


# ==============================================================================
# 1. SHMI management records
# ==============================================================================
shmi_inputs <- prepare_shmi_inputs(WORKBOOK, exclude = EXCLUDE_ALL)
stopifnot(identical(attr(shmi_inputs, "rotation_window"), "calendar"))
saveRDS(shmi_inputs, OUT("shmi_inputs_v120.rds"))


# ==============================================================================
# 2. Cropping-system classification (defines the SHMI calibration scope)
#    Uses the crop episodes clipped to each unit's evaluation window, i.e.
#    exactly the records the index is computed from.
# ==============================================================================
crop_eval <- clip_crop_to_rotation(shmi_inputs$crop, shmi_inputs$rot_bounds)

# Reference soils: unmanaged treatments (sod, remnant prairie, benchmarks,
# exclosures, native grassland), from the treatment descriptions
reference_ids <- read.csv(IN("napeshm", "treatments.csv")) %>%
  mutate(MGT_combo = make_combo(sitecode, treatmentid)) %>%
  filter(eu_crop_cat == "non_farmed" | grepl("not_grazed", eu_treat_desc)) %>%
  pull(MGT_combo)

systems <- crop_eval %>%
  group_by(MGT_combo) %>%
  summarise(
    has_annual = any(CD_cat == "Annual"),
    has_fallow = any(CD_cat == "Fallow"),
    has_perenn = any(CD_cat == "Perennial"),
    crop_cats  = paste(sort(unique(CD_cat)), collapse = " + "),
    .groups = "drop"
  ) %>%
  mutate(
    system = case_when(
      MGT_combo %in% reference_ids ~ "reference",
      has_annual                   ~ "annual",
      has_fallow & !has_perenn     ~ "fallow",
      TRUE                         ~ "perennial"
    ),
    in_scope = system %in% c("annual", "fallow")
  )

cat("\nCropping systems (from management records):\n")
print(count(systems, system, crop_cats))

in_scope_ids <- systems$MGT_combo[systems$in_scope]

saveRDS(in_scope_ids, OUT("in_scope_v120.rds"))
saveRDS(list(systems  = systems,
             in_scope = in_scope_ids,
             exclude  = EXCLUDE,
             rule     = paste("in scope = at least one annual crop in the evaluation",
                              "window, or a fallow-only reference treatment;",
                              "perennial-only systems excluded"),
             workbook = basename(WORKBOOK),
             shmi_version = as.character(packageVersion("SHMI")),
             created  = Sys.time()),
        OUT("scope_v120.rds"))

# Records flagged for review by the data owner, marked in or out of scope
checks <- shmi_inputs$assumptions %>%
  filter(level == "check") %>%
  mutate(in_scope = MGT_combo %in% in_scope_ids) %>%
  arrange(desc(in_scope), type, MGT_combo)
write.csv(checks, OUT("data_checks_for_review.csv"), row.names = FALSE)
cat("\nData checks for review:\n")
print(table(checks$type, in_scope = checks$in_scope))


# ==============================================================================
# 3. Site, treatment, soil and PLFA data
# ==============================================================================
SITE_COLS <- c("country", "state_prov", "siteid", "sitecode",
               "site_precip", "site_meanmax_temp", "site_meanmin_temp",
               "site_mean_temp", "hargreave_cmd", "mi", "site_elevation", "gdd10",
               "site_latitude", "site_longitude", "na_l1name",
               "site_animals", "site_corn", "site_cotton", "site_cover_crop",
               "site_crop_sys", "site_disturbance", "site_native",
               "site_nutr_rate", "site_nutr_type", "site_residue", "site_soy",
               "site_water", "site_wheat")

sites <- read.csv(IN("napeshm", "sites.csv"))
missing_cols <- setdiff(SITE_COLS, names(sites))
if (length(missing_cols)) stop("sites.csv lacks: ", paste(missing_cols, collapse = ", "))
sites <- sites %>%
  select(all_of(SITE_COLS)) %>%
  filter(!sitecode %in% EXCLUDE_SITES)

treatments <- read.csv(IN("napeshm", "treatments.csv")) %>%
  mutate(MGT_combo = make_combo(sitecode, treatmentid))

# Lab data: average replicate samples per treatment. Means of all-NA groups
# (NaN) are set to NA.
mean_by_unit <- function(df) {
  df %>%
    mutate(MGT_combo = make_combo(sitecode, treatmentid)) %>%
    group_by(MGT_combo) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    mutate(across(where(is.numeric), ~ ifelse(is.nan(.x), NA_real_, .x)))
}

soil <- read.csv(IN("napeshm", "soildata.csv")) %>%
  mean_by_unit() %>%
  select(MGT_combo, starts_with("a_"), starts_with("b_"), starts_with("d_"),
         shi_kfs_cmhr, slakes_gmean)

plfa <- read.csv(IN("napeshm", "plfa_ng.csv")) %>%
  mean_by_unit() %>%
  select(MGT_combo, mpng_total, starts_with("mt_"))

units <- sites %>%
  left_join(treatments, by = c("sitecode", "siteid")) %>%
  filter(!MGT_combo %in% EXCLUDE_ALL) %>%
  left_join(soil, by = "MGT_combo") %>%
  left_join(plfa, by = "MGT_combo") %>%
  mutate(treatmentid = as.character(treatmentid),
         trt_age     = pmin(TRT_REF_YEAR - treatment_start_year + 1, TRT_AGE_CAP))


# ==============================================================================
# 4. CASH indicators used by the SEM (plus texture and descriptors)
# ==============================================================================
semDF <- units %>%
  mutate(
    # biological activity
    Cmin   = d_co2_c,               # respiration (mineralizable C)
    Nmin   = b_pmn_am_mgkg,         # potentially mineralizable N
    ARYL   = d_aryl,                # arylsulfatase
    NAG    = d_beta_minidase,       # N-acetyl-beta-glucosaminidase
    # carbon and nitrogen pools
    SOC    = b_soc,
    TN     = b_total_n,
    POXC   = a_act_c,               # permanganate-oxidizable (active) C
    ACE    = a_ace,                 # autoclaved citrate-extractable protein
    # physical condition
    SSI    = b_ssi,                 # slake stability index
    WAS    = a_stabi,               # wet aggregate stability
    Ksat   = shi_kfs_cmhr,          # saturated hydraulic conductivity
    # other measurements (not in the SEM, kept for reference)
    TC     = b_total_c,
    PLFA   = mpng_total / 1000,
    BD     = a_bd_intact_av,
    pH     = b_ph,
    AWHC   = a_awhc_intact_vol,
    # texture (covariate for the climate/texture adjustment)
    Clay   = b_clay,
    Silt   = b_silt,
    Clay_Silt = b_clay + b_silt,
    TrtAge = trt_age
  )


# ==============================================================================
# 5. Attach the system classification, check coverage, save
# ==============================================================================
semDF <- semDF %>%
  left_join(systems %>% select(MGT_combo, system, in_scope, crop_cats), by = "MGT_combo") %>%
  mutate(system   = coalesce(system, "no_management_records"),
         in_scope = coalesce(in_scope, FALSE),
         sitecode = droplevels(as.factor(sitecode)))

# Coverage checks
no_mgmt <- semDF$MGT_combo[semDF$system == "no_management_records"]
no_soil <- setdiff(in_scope_ids, semDF$MGT_combo)
if (length(no_mgmt)) warning(length(no_mgmt), " units have soil data but no SHMI records: ",
                             paste(head(no_mgmt, 5), collapse = ", "))
if (length(no_soil)) warning(length(no_soil), " in-scope units have no site/soil record: ",
                             paste(head(no_soil, 5), collapse = ", "))

cat("\nUnits by cropping system (management records) and scope:\n")
print(table(semDF$system, in_scope = semDF$in_scope, useNA = "ifany"))
cat("\nSEM data:", nrow(semDF), "units at", nlevels(semDF$sitecode), "sites;",
    sum(semDF$in_scope), "in scope at",
    length(unique(semDF$sitecode[semDF$in_scope])), "sites\n")

# ---- SHAPE soil taxonomy groups (from A1a_Get_site_ssurgo.R) ------------------------
# S2 and S3 are merged (too few sites alone); S4 is the reference group.
# Mexican sites have no SSURGO data and stay NA (only the US soil-group run uses it).
SOILS_FILE <- OUT("site_soils_ssurgo.csv")
if (file.exists(SOILS_FILE)) {
  site_soils <- read.csv(SOILS_FILE, stringsAsFactors = FALSE) %>%
    transmute(sitecode = as.character(sitecode), shape_group)
  semDF <- semDF %>% mutate(sitecode = as.character(sitecode)) %>%
    left_join(site_soils, by = "sitecode") %>%
    mutate(shape_grp3 = factor(case_when(shape_group %in% c("S1", "S2", "S3") ~ "S2_3",
                                         shape_group == "S4" ~ "S4",
                                         shape_group == "S5" ~ "S5"),
                               levels = c("S4", "S2_3", "S5")))
  cat("\nSHAPE soil groups (units):\n"); print(table(semDF$shape_grp3, useNA = "ifany"))
  us_missing <- semDF %>% filter(!grepl("^MX", sitecode), is.na(shape_grp3)) %>% distinct(sitecode)
  if (nrow(us_missing)) warning("US sites without a SHAPE group: ", paste(us_missing$sitecode, collapse = ", "))
} else {
  message("No ", SOILS_FILE, ": SHAPE soil groups not added (run A1a_Get_site_ssurgo.R first)")
}

saveRDS(semDF, OUT("semDF_prepped.rds"))
write.csv(semDF, OUT("semDF_prepped.csv"), row.names = FALSE)
