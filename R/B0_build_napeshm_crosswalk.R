################################################################################
# B0_build_napeshm_crosswalk.R
#
# NAPESHM soil-health measurements for the treatments of the three yield sites that
# NAPESHM also sampled (Hoytville USOH02, Wooster USOH03, Arlington USWI01), for the
# mediation test in B5_SHMI_soil_mediation.R.
#
# Inputs
#   data/napeshm_yield_site_units.csv        NAPESHM unit -> yield site and treatment
#                                            (rotation, tillage, treatment code), by hand
#                                            from the NAPESHM and DRIVES treatment lists
#   Outputs/SHAPE/SHAPE_by_unit_with_SHMI.csv measured soil health (SH = lv_SH), SHAPE,
#                                            indicator values and SHMI per NAPESHM unit
#                                            (A6b_Merge_SHAPE_SHMI.R)
# Output
#   Outputs/yield/napeshm_mediators_crosswalk.csv
################################################################################
suppressPackageStartupMessages(library(dplyr))
setwd(here::here())

map   <- read.csv(here::here("data", "napeshm_yield_site_units.csv"), stringsAsFactors = FALSE)
shape <- read.csv(here::here("Outputs", "SHAPE", "SHAPE_by_unit_with_SHMI.csv"), stringsAsFactors = FALSE)

cw <- map %>%
  left_join(shape, by = c("napeshm_unit" = "MGT_combo")) %>%
  transmute(napeshm_unit, napeshm_site, site, rotation, tillage, unit, level,
            lv_SH        = round(SH, 2),           # latent soil health (SEM), 0-100
            SHAPE_mean   = round(SHAPE_mean, 3),   # mean SHAPE score of SOC, AC, Resp, WAS
            SOC_gkg      = round(10 * SOC_raw, 2), # % -> g/kg
            POXC_mgkg    = round(AC_raw, 1),       # active (permanganate-oxidizable) C
            Resp         = round(Resp_raw, 3),
            WAS          = round(KSSL_WAS_raw, 2), # water-stable aggregates
            ACE          = round(ACE_raw, 2),      # ACE protein
            SHMI_napeshm = round(SHMI, 1))
miss <- cw$napeshm_unit[is.na(cw$lv_SH)]
if (length(miss)) stop("NAPESHM units missing from SHAPE_by_unit_with_SHMI.csv: ", paste(miss, collapse = ", "))

dir.create(here::here("Outputs", "yield"), recursive = TRUE, showWarnings = FALSE)
write.csv(cw, here::here("Outputs", "yield", "napeshm_mediators_crosswalk.csv"), row.names = FALSE)
cat("Saved napeshm_mediators_crosswalk.csv:", nrow(cw), "treatments at", paste(unique(cw$site), collapse = ", "), "\n")
