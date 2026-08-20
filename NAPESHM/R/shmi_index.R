############################################################
# shmi_index.R — Build SHMI inputs + SHMI subindices
############################################################

library(dplyr)
library(tidyr)
library(lubridate)
library(purrr)
library(janitor)
library(scales)
library(stringr)

############################################################
# 1. prepare_shmi_inputs()
############################################################

prepare_shmi_inputs <- function(path,
                                exclude = NULL,
                                verbose = TRUE,
                                start_year_override = NULL,
                                end_year_override = NULL) {
  
  # ------------------------------------------------------------
  # 1. Exclusion lists
  # ------------------------------------------------------------
  EXCLUDE_WOODY      <- c("NAPESHM_USMI04_306", "NAPESHM_USNC03_368", "NAPESHM_USNC03_369",
                          "NAPESHM_USNC03_370", "NAPESHM_USNC03_371", "NAPESHM_USNC03_372")
  EXCLUDE_BENCHMARK  <- c("NAPESHM_USNY01_395", "NAPESHM_USNY02_400",
                          "NAPESHM_USNY03_403", "NAPESHM_USNY04_406")
  EXCLUDE_NOINF      <- c("NAPESHM_USOK01_429", "NAPESHM_USOK01_430",
                          "NAPESHM_USOK01_431", "NAPESHM_USOK01_432")
  default_exclude <- c(EXCLUDE_WOODY, EXCLUDE_BENCHMARK, EXCLUDE_NOINF)
  
  # If user supplies exclusions, combine them; otherwise use defaults
  exclude <- unique(c(default_exclude, exclude))
  
  if (verbose) {
    message("Excluding ", length(exclude), " management units:")
    message(paste("  -", exclude, collapse = "\n"))
  }
  
  # ------------------------------------------------------------
  # 2. Helper to read + filter sheets
  # ------------------------------------------------------------
  read_sheet <- function(sheet, skip=3, range=NULL, drop_cols = NULL) {
    df <- read_xlsx(path, sheet = sheet, skip = skip, range = range) %>%
      left_join(mgt, by = c("user_name", "MGT_combo")) %>%
      dplyr::select(-user_name, -MGT_id, all_of(drop_cols)) %>%
      janitor::remove_empty("rows") %>%
      filter(!(MGT_combo %in% exclude))
    
    if (verbose) {
      # Print included combos for this sheet
      included <- unique(df$MGT_combo)
      message("\nIncluded ", length(included), " management units for sheet '", sheet, "':")
      message(paste("  -", included, collapse = "\n"))
      
    }
    df
  }
  
  # ------------------------------------------------------------
  # 3. Load MGT first (needed for joins)
  # ------------------------------------------------------------
  mgt <- read_xlsx(path, sheet = "Mgt_Unit", skip = 3, range = NULL) %>%
    dplyr::select(user_name, MGT_combo) %>%
    janitor::remove_empty("rows") %>%
    filter(!(MGT_combo %in% exclude))
  
  if (verbose) {
    included <- unique(mgt$MGT_combo)
    message("\nIncluded ", length(included), " management units for sheet Mgt_Unit:")
    message(paste("  -", included, collapse = "\n"))
  }
  
  # ------------------------------------------------------------
  # 4. Load all sheets
  # ------------------------------------------------------------
  parse_date_safe <- function(x) {
    suppressWarnings({
      # Try ymd first
      d <- lubridate::ymd(x, quiet = TRUE)
      
      # If still NA, try mdy
      d[is.na(d)] <- lubridate::mdy(x[is.na(d)], quiet = TRUE)
      
      # If still NA, try dmy
      d[is.na(d)] <- lubridate::dmy(x[is.na(d)], quiet = TRUE)
      
      # If still NA, try Excel serial numbers
      suppressWarnings({
        nums <- suppressWarnings(as.numeric(x))
        d[is.na(d) & !is.na(nums)] <- as.Date(nums[is.na(d) & !is.na(nums)], origin = "1899-12-30")
      })
      
      # Final: return Date or NA
      d
    })
  }
  
  crop <- read_sheet("Crop_Diversity", skip = 3, drop_cols = c("CD_id")) %>%
    janitor::remove_empty("rows") %>%
    mutate(across(c(CD_plant_date, CD_harv_date, CD_term_date), parse_date_safe))
  
  if (!is.null(start_year_override)) {
    crop <- crop %>%
      filter(
        lubridate::year(CD_plant_date) >= start_year_override |
          lubridate::year(CD_harv_date)  >= start_year_override |
          lubridate::year(CD_term_date)  >= start_year_override
      )
  }
  
  if (!is.null(end_year_override)) {
    crop <- crop %>%
      filter(
        lubridate::year(CD_plant_date) <= end_year_override |
          lubridate::year(CD_harv_date)  <= end_year_override |
          lubridate::year(CD_term_date)  <= end_year_override
      )
  }
  
  dist <- read_sheet("Soil_Disturbance", skip = 3) %>%
    janitor::remove_empty("rows") %>%
    mutate(across(c(SD_date), parse_date_safe))
  
  if (!is.null(start_year_override)) {
    dist <- dist %>%
      filter(lubridate::year(SD_date) >= start_year_override)
  }
  
  if (!is.null(end_year_override)) {
    dist <- dist %>% filter(lubridate::year(SD_date) <= end_year_override)
  }
  
  amend <- read_sheet("Soil_Amendments", skip = 3) %>%
    janitor::remove_empty("rows") %>%
    mutate(across(c(SA_date), parse_date_safe))
  
  if (!is.null(start_year_override)) {
    amend <- amend %>%
      filter(lubridate::year(SA_date) >= start_year_override)
  }
  
  if (!is.null(end_year_override)) {
    amend <- amend %>% filter(lubridate::year(SA_date) <= end_year_override)
  }
  
  animal <- read_sheet("Animal_Diversity", skip = 3) %>%
    janitor::remove_empty("rows") %>%
    mutate(across(c(AD_start_date, AD_end_date), parse_date_safe))
  
  if (!is.null(start_year_override)) {
    animal <- animal %>%
      filter(
        lubridate::year(AD_start_date) >= start_year_override |
          lubridate::year(AD_end_date)   >= start_year_override
      )
  }
  
  if (!is.null(end_year_override)) {
    animal <- animal %>%
      filter(
        lubridate::year(AD_start_date) <= end_year_override |
          lubridate::year(AD_end_date)   <= end_year_override
      )
  }
  
  # ------------------------------------------------------------
  # 5. Harmonize crop windows
  # ------------------------------------------------------------
  harmonize_crop_windows <- function(crop) {
    
    crop %>%
      mutate(
        # ensure dates are Date
        CD_plant_date = as.Date(CD_plant_date),
        CD_harv_date  = as.Date(CD_harv_date),
        CD_term_date  = as.Date(CD_term_date),
      ) %>%
      group_by(MGT_combo) %>%
      # need a start and end year for each MGT_combo
      mutate(
        start_yr = suppressWarnings(
          min(
            lubridate::year(CD_plant_date),
            lubridate::year(CD_harv_date),
            lubridate::year(CD_term_date),
            na.rm = TRUE
          )
        ),
        end_yr = suppressWarnings(
          max(
            lubridate::year(CD_plant_date),
            lubridate::year(CD_harv_date),
            lubridate::year(CD_term_date),
            na.rm = TRUE
          )
        )
      ) %>%
      ungroup() %>%
      # now fill in missing planting/termination dates
      # for first and last crop in each MGT_combo
      group_by(MGT_combo, CD_seq_num) %>%
      mutate(
        # identify first/last crop in rotation
        is_first = CD_seq_num == min(CD_seq_num),
        is_last  = CD_seq_num == max(CD_seq_num),
        
        # raw crop windows across all species in the mixture
        # date may be missing (to be filled in later)
        start_raw = suppressWarnings(
          min(CD_plant_date, na.rm = TRUE)
        ),
        end_raw = suppressWarnings(
          max(CD_harv_date, CD_term_date, na.rm = TRUE)
        ),
        start_raw = dplyr::if_else(is.infinite(start_raw), NA_Date_, start_raw),
        end_raw   = dplyr::if_else(is.infinite(end_raw),   NA_Date_, end_raw)
      ) %>%
      # collapse to one row per event (mixture → event)
      summarize(
        CD_group = paste(unique(CD_group), collapse = " + "),
        CD_name  = paste(unique(CD_name), collapse = " + "),
        crop_start = dplyr::first(
          dplyr::if_else(
            is_first & is.na(start_raw),
            as.Date(paste0(start_yr, "-01-01")),   # fallback
            start_raw
          )
        ),
        crop_end = dplyr::first(
          dplyr::if_else(
            is_last & is.na(end_raw),
            as.Date(paste0(end_yr, "-12-31")),   # fallback
            end_raw
          )
        ),
        .groups = "drop"
      )
  }
  
  crop_harmonized <- harmonize_crop_windows(crop)
  
  if (!is.null(start_year_override)) {
    start_cut <- as.Date(paste0(start_year_override, "-01-01"))
    crop_harmonized <- crop_harmonized %>%
      mutate(
        crop_start = pmax(crop_start, start_cut)
      )
  }
  
  if (!is.null(end_year_override)) {
    end_cut <- as.Date(paste0(end_year_override, "-12-31"))
    crop_harmonized <- crop_harmonized %>%
      mutate(
        crop_end = pmin(crop_end, end_cut)
      )
  }
  
  # ------------------------------------------------------------
  # 6. Bounds helper
  # ------------------------------------------------------------
  compute_bounds <- function(crop_harmonized, dist, amend, animal) {
    df1 <- crop_harmonized %>%
      dplyr::select(MGT_combo, crop_start, crop_end) %>%
      group_by(MGT_combo) %>%
      summarize(start_year = min(lubridate::year(crop_start), na.rm = TRUE),
                end_year = max(lubridate::year(crop_end), na.rm = TRUE))
    df2 <- dist %>%
      dplyr::select(MGT_combo, SD_date) %>%
      group_by(MGT_combo) %>%
      summarize(start_year = min(lubridate::year(SD_date), na.rm = TRUE),
                end_year = max(lubridate::year(SD_date), na.rm = TRUE))
    
    df3 <- amend %>%
      dplyr::select(MGT_combo, SA_date) %>%
      group_by(MGT_combo) %>%
      summarize(start_year = min(lubridate::year(SA_date), na.rm = TRUE),
                end_year = max(lubridate::year(SA_date), na.rm = TRUE))
    
    df4 <- animal %>%
      dplyr::select(MGT_combo, AD_start_date, AD_end_date) %>%
      group_by(MGT_combo) %>%
      summarize(start_year = min(lubridate::year(AD_start_date), na.rm = TRUE),
                end_year = max(lubridate::year(AD_end_date), na.rm = TRUE))
    
    df5 <- rbind(df1, df2, df3, df4) %>%
      group_by(MGT_combo) %>%
      summarize(
        rot_start_yr = min(start_year, na.rm = TRUE),
        rot_end_yr = max(end_year, na.rm = TRUE)) %>%
      mutate(rot_start = as.Date(paste0(rot_start_yr, "-01-01")),
             rot_end = as.Date(paste0(rot_end_yr, "-12-31")))
    
    df5
    
  }
  
  rot_bounds <- compute_bounds(crop_harmonized, dist, amend, animal)
  
  if (!is.null(start_year_override)) {
    rot_bounds$rot_start_yr <- pmax(rot_bounds$rot_start_yr, start_year_override)
    rot_bounds$rot_start <- as.Date(paste0(rot_bounds$rot_start_yr, "-01-01"))
  }
  
  if (!is.null(end_year_override)) {
    rot_bounds$rot_end_yr <- pmin(rot_bounds$rot_end_yr, end_year_override)
    rot_bounds$rot_end <- as.Date(paste0(rot_bounds$rot_end_yr, "-12-31"))
  }
  
  # ------------------------------------------------------------
  # 10. Daily disturbance summary
  # ------------------------------------------------------------
  daily_dist <- dist %>%
    dplyr::select(MGT_combo, SD_date, SD_mixeff) %>%
    group_by(MGT_combo, SD_date) %>%
    summarize(
      SDmax = max(SD_mixeff, na.rm = TRUE),
      SDsum = sum(SD_mixeff, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(date = SD_date)
  
  
  # ------------------------------------------------------------
  # 11. Build daily grid
  # ------------------------------------------------------------
  
  build_daily_grid <- function(crop_harmonized, rot_bounds, daily_dist) {
    
    if (verbose) {
      message("Building the rotation grid...")
    }
    
    rot_grid <- rot_bounds %>%
      rowwise() %>%
      mutate(date = list(seq.Date(rot_start, rot_end, by = "day"))) %>%
      unnest(date) %>%
      ungroup()
    
    if (verbose) {
      message("Building the crop grid...")
    }
    
    crop_grid <- crop_harmonized %>%
      rowwise() %>%
      mutate(date = list(seq.Date(crop_start, crop_end, by = "day"))) %>%
      unnest(date) %>%
      ungroup() %>%
      mutate(crop_present = 1L)
    
    if (verbose) {
      message("Combining grids...")
    }
    
    daily_crop <- rot_grid %>%
      left_join(crop_grid, by = c("MGT_combo", "date"), relationship = "many-to-many") %>%
      mutate(crop_present = replace_na(crop_present, 0L))
    
    temp <- daily_crop %>%
      left_join(daily_dist, by = c("MGT_combo", "date")) %>%
      mutate(
        SDmax = replace_na(SDmax, 0),
        SDsum = replace_na(SDsum, 0)
      )
  }
  
  daily <- build_daily_grid(crop_harmonized, rot_bounds, daily_dist)
  
  # ------------------------------------------------------------
  # 12. Return everything in one clean list
  # ------------------------------------------------------------
  list(
    rot_bounds      = rot_bounds,
    crop_harmonized = crop_harmonized,
    daily           = daily,
    daily_dist      = daily_dist,
    mgt             = mgt,
    crop            = crop,
    dist            = dist,
    amend           = amend,
    animal          = animal
  )
}

############################################################
# 2. Cover/Roots subindex
############################################################

compute_w.cover <- function(daily, rot_bounds,
                            w_winter = 0.25,
                            w_spring = 0.25,
                            w_summer = 0.25,
                            w_fall   = 0.25) {
  
  crop_season <- daily %>%
    mutate(
      crop_present = if_else(
        !is.na(CD_name) & CD_name == "fallow",
        0L,
        coalesce(crop_present, 0L)
      ),
      month = lubridate::month(date),
      season = case_when(
        month %in% c(12, 1, 2)  ~ "winter",
        month %in% c(3, 4, 5)   ~ "spring",
        month %in% c(6, 7, 8)   ~ "summer",
        month %in% c(9, 10, 11) ~ "fall"
      )
    )
  
  crop_days <- crop_season %>%
    group_by(MGT_combo, season) %>%
    summarize(
      plant_days = sum(crop_present, na.rm = TRUE),
      days_possible = n(),
      .groups = "drop"
    )
  
  season_totals <- crop_days %>%
    tidyr::pivot_wider(
      names_from = season,
      values_from = c(plant_days, days_possible),
      values_fill = 0
    ) %>%
    mutate(
      prop_winter = if_else(days_possible_winter > 0,
                            plant_days_winter / days_possible_winter, 0),
      prop_spring = if_else(days_possible_spring > 0,
                            plant_days_spring / days_possible_spring, 0),
      prop_summer = if_else(days_possible_summer > 0,
                            plant_days_summer / days_possible_summer, 0),
      prop_fall   = if_else(days_possible_fall > 0,
                            plant_days_fall   / days_possible_fall, 0)
    )
  
  w_sum <- w_winter + w_spring + w_summer + w_fall
  w_winter <- w_winter / w_sum
  w_spring <- w_spring / w_sum
  w_summer <- w_summer / w_sum
  w_fall   <- w_fall   / w_sum
  
  season_totals %>%
    mutate(
      Cover = 100 * (
        w_winter * prop_winter +
          w_spring * prop_spring +
          w_summer * prop_summer +
          w_fall   * prop_fall
      )
    ) %>%
    dplyr::select(MGT_combo, Cover)
}


############################################################
# 3. Diversity subindex
############################################################

compute_rot_diversity <- function(crop_harmonized, daily,
                                  hill = 2, max_div = 8) {
  
  daily_sum <- daily %>%
    ungroup() %>%
    mutate(year = year(date)) %>%
    group_by(MGT_combo, CD_seq_num, CD_name, year) %>%
    summarize(days = sum(crop_present), .groups = "drop") %>%
    filter(!is.na(CD_seq_num))
  
  expand_mixtures <- function(df) {
    placeholder <- df %>%
      filter(str_detect(CD_name, "-species")) %>%
      mutate(
        mix_n = as.integer(str_extract(CD_name, "\\d+")),
        species = map(mix_n, ~ paste0("species_", seq_len(.x)))
      ) %>%
      unnest(species)
    
    realmix <- df %>%
      filter(!str_detect(CD_name, "-species")) %>%
      mutate(species = str_split(CD_name, "\\s*\\+\\s*")) %>%
      unnest(species) %>%
      mutate(species = str_trim(species))
    
    bind_rows(placeholder, realmix) %>%
      dplyr::select(MGT_combo, species, year, days)
  }
  
  species_expanded <- expand_mixtures(daily_sum)
  
  plant_days <- species_expanded %>%
    mutate(days = ifelse(species == "fallow", 0, days)) %>%
    group_by(MGT_combo, species, year) %>%
    summarize(days = sum(days), .groups = "drop")
  
  div_rot <- plant_days %>%
    group_by(MGT_combo, species) %>%
    summarize(days = sum(days), .groups = "drop") %>%
    group_by(MGT_combo) %>%
    mutate(p = days / sum(days)) %>%
    summarize(
      D = case_when(
        hill == 0 ~ sum(p > 0),
        hill == 1 ~ -sum(p * log(p), na.rm = TRUE),
        TRUE      ~ -log(sum(p^2, na.rm = TRUE))
      ),
      .groups = "drop"
    ) %>%
    mutate(D = ifelse(is.na(D), 0, D),
           D = pmin(D, max_div))
  
  div_rot %>%
    mutate(
      Diversity_raw = case_when(
        hill == 0 ~ D / max_div,
        TRUE      ~ D / log(max_div)
      ),
      Diversity = pmin(Diversity_raw, 1) * 100
    ) %>%
    dplyr::select(MGT_combo, Diversity)
}

############################################################
# 4. Disturbance subindex
############################################################

compute_avg_annual_disturbance <- function(daily_dist, rot_bounds,
                                           max_stir = 190,
                                           transform = c("sqrt", "log", "raw")) {
  
  transform <- match.arg(transform)
  
  full_years <- rot_bounds %>%
    mutate(year = map2(rot_start_yr, rot_end_yr, seq)) %>%
    unnest(year) %>%
    dplyr::select(MGT_combo, year)
  
  stir_years <- daily_dist %>%
    arrange(MGT_combo, date) %>%
    mutate(year = year(date)) %>%
    group_by(MGT_combo, year) %>%
    summarize(STIR = sum(SDsum), .groups = "drop")
  
  stir <- full_years %>%
    left_join(stir_years, by = c("MGT_combo", "year")) %>%
    mutate(STIR = replace_na(STIR, 0),
           STIR = pmin(STIR, max_stir))
  
  stir <- stir %>%
    mutate(STIR = case_when(
      transform == "raw"  ~ STIR,
      transform == "sqrt" ~ sqrt(STIR),
      transform == "log"  ~ log10(STIR + 1)
    )) %>%
    group_by(MGT_combo) %>%
    summarize(STIR = mean(STIR), .groups = "drop") %>%
    mutate(InvDist = 100 - rescale(STIR, to = c(0, 100)))
  
  stir %>% dplyr::select(MGT_combo, InvDist)
}

############################################################
# 5. Organic Inputs / Animals subindex
############################################################

compute_orginput <- function(rot_bounds,
                             amend,
                             animal,
                             w_amend = 1,
                             w_animal = 1) {
  
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(purrr)
  
  # ---------------------------------------------------------------------------
  # 1. Build year grid for each MGT_combo
  # ---------------------------------------------------------------------------
  rot_grid <- rot_bounds %>%
    mutate(year = map2(rot_start_yr, rot_end_yr, seq)) %>%
    unnest(year) %>%
    dplyr::select(MGT_combo, year)
  
  # ---------------------------------------------------------------------------
  # 2. Amendment events (binary per year)
  # ---------------------------------------------------------------------------
  amend_events <- amend %>%
    filter(SA_cat == "Organic") %>%
    mutate(year = year(SA_date)) %>%
    distinct(MGT_combo, year) %>%
    mutate(amend_present = 1)
  
  # ---------------------------------------------------------------------------
  # 3. Animal events (binary per year)
  # ---------------------------------------------------------------------------
  ani_events <- animal %>%
    mutate(year = year(AD_start_date)) %>%
    distinct(MGT_combo, year) %>%
    mutate(ani_present = 1)
  
  # ---------------------------------------------------------------------------
  # 4. Join and compute weighted organic-input presence
  # ---------------------------------------------------------------------------
  bio_events <- rot_grid %>%
    left_join(amend_events, by = c("MGT_combo", "year")) %>%
    left_join(ani_events,   by = c("MGT_combo", "year")) %>%
    mutate(
      amend_present = replace_na(amend_present, 0),
      ani_present   = replace_na(ani_present, 0),
      # weighted presence per year
      weighted_input = w_amend * amend_present +
        w_animal * ani_present
    )
  
  # ---------------------------------------------------------------------------
  # 5. Compute rotation length (years)
  # ---------------------------------------------------------------------------
  rot_meta <- rot_bounds %>%
    mutate(rot_years = rot_end_yr - rot_start_yr + 1) %>%
    dplyr::select(MGT_combo, rot_years)
  
  # ---------------------------------------------------------------------------
  # 6. Compute weighted events per year
  # ---------------------------------------------------------------------------
  animals_rate <- bio_events %>%
    group_by(MGT_combo) %>%
    summarise(
      total_weighted = sum(weighted_input),
      .groups = "drop"
    ) %>%
    left_join(rot_meta, by = "MGT_combo") %>%
    mutate(
      events_per_year = ifelse(rot_years > 0,
                               total_weighted / rot_years,
                               0)
    )
  
  # ---------------------------------------------------------------------------
  # 7. Scale to 0–100 and rename to Animals
  # ---------------------------------------------------------------------------
  animals_scaled <- animals_rate %>%
    mutate(
      Animals = scales::rescale(events_per_year, to = c(0, 100))
    ) %>%
    dplyr::select(MGT_combo, events_per_year, Animals)
  
  return(animals_scaled)
}
