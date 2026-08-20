############################################################
# data_prep.R — Load, clean, merge, and adjust NAPESHM data
############################################################

load_napes_data <- function() {
  
  # Load site metadata
  sitesDF <- read.csv("NAPESHM/data/sites.DKM.csv") %>%
    dplyr::select(country, state_prov, siteid, sitecode, ClimateZ,
                  site_precip, site_meanmax_temp, site_meanmin_temp,
                  site_mean_temp, hargreave_cmd,
                  site_animals, site_corn, site_cotton, site_cover_crop,
                  site_crop_sys, site_disturbance, site_native,
                  site_nutr_rate, site_nutr_type, site_residue,
                  site_soy, site_water, site_wheat) %>%
    filter(sitecode != "USWA03")
  
  # Treatment data
  trtDF <- read.csv("NAPESHM/data/treatments.DKM.csv") %>%
    mutate(MGT_combo = paste0("NAPESHM_", sitecode, "_", treatmentid))
  
  shmiDF <- sitesDF %>% left_join(trtDF, by = c("sitecode", "siteid"))
  
  # Exclusions
  EXCLUDE_ALL <- c(
    "NAPESHM_USMI04_306", "NAPESHM_USNC03_368", "NAPESHM_USNC03_369",
    "NAPESHM_USNC03_370", "NAPESHM_USNC03_371", "NAPESHM_USNC03_372",
    "NAPESHM_USNY01_395", "NAPESHM_USNY02_400", "NAPESHM_USNY03_403",
    "NAPESHM_USNY04_406", "NAPESHM_USOK01_429", "NAPESHM_USOK01_430",
    "NAPESHM_USOK01_431", "NAPESHM_USOK01_432"
  )
  
  shmiDF <- shmiDF %>% filter(!(MGT_combo %in% EXCLUDE_ALL))
  
  # Soil data
  soilDF <- read.csv("NAPESHM/data/soildata.csv") %>%
    mutate(MGT_combo = paste0("NAPESHM_", sitecode, "_", treatmentid)) %>%
    group_by(MGT_combo) %>%
    summarise_if(is.numeric, mean, na.rm = TRUE, .groups = "drop")
  
  shmiDF <- shmiDF %>% left_join(soilDF, by = c("MGT_combo", "siteid"))
  
  # PLFA data
  plfaDF <- read.csv("NAPESHM/data/plfa_ng.csv") %>%
    mutate(MGT_combo = paste0("NAPESHM_", sitecode, "_", treatmentid)) %>%
    group_by(MGT_combo) %>%
    summarise_if(is.numeric, mean, na.rm = TRUE, .groups = "drop")
  
  finalDF <- shmiDF %>% left_join(plfaDF, by = c("MGT_combo", "siteid", "euid"))
  
  # Treatment age
  finalDF$trt_age <- pmin((2019 - finalDF$treatment_start_year) + 1, 50)
  
  # CASH indicator creation
  finalDF <- finalDF %>%
    mutate(
      Clay = b_clay,
      Silt = b_silt,
      BD   = a_bd_intact_av,
      pH   = b_ph,
      ACE  = a_ace,
      POXC = a_act_c,
      AWHC = a_awhc_intact_vol,
      Nmin = b_pmn_am_mgkg,
      Ksat = shi_kfs_cmhr,
      Cmin = d_co2_c,
      WAS  = a_stabi,
      PLFA = mpng_total / 1000,
      ARYL = d_aryl,
      NAG  = d_beta_minidase,
      SOC  = b_soc,
      SSI  = b_ssi,
      TC   = b_total_c,
      TN   = b_total_n,
      Clay_Silt = b_clay + b_silt
    )
  
  # Outlier removal
  finalDF <- finalDF %>%
    group_by(ClimateZ) %>%
    mutate(across(c(ACE, POXC, Cmin, WAS, AWHC, SOC, PLFA, ARYL,
                    NAG, SSI, Nmin, Ksat, TC, TN), replace_outliers)) %>%
    ungroup()
  
  # Climate/texture adjustment
  for (v in c("ACE","POXC","Cmin","WAS","AWHC","SOC","PLFA","ARYL",
              "NAG","SSI","Nmin","Ksat","TC","TN")) {
    finalDF[[v]] <- adjust_indicator(finalDF, v)
  }
  
  return(finalDF)
}

############################################################
# End data_prep.R
############################################################
