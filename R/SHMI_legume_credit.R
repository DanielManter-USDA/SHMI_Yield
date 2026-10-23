################################################################################
# SHMI_legume_credit.R - legume N credits for the yield analyses
# Sourced by SHMI_yield_clean.R and SHMI_yield_onfarm.R.
#
# Corn following a legume needs less fertilizer N. The SHMI package reports only
# the N recorded as applied, so the yield analyses add a legume N credit to it:
# N used = recorded fertilizer and amendment N + the credit for the crop(s) that
# preceded the corn. Credits are in lb N/ac (converted to kg/ha), with low,
# default and high values for a sensitivity analysis.
#
# Values are in the range of Land-Grant guidelines (e.g. UW A2809); confirm them
# for your region. Soybean: 40 lb N/ac by default, the difference between the
# fertilizer N given to continuous corn (CS1) and to corn after soybean (CS2) at
# WICST in 2015-2020 (45 kg/ha); over 1990-2020 the difference was 16 lb/ac.
#
# Rules
#   - a credit comes from a legume (or legume-containing mix) that ended (harvest,
#     termination or plow-down) within LEGUME_LOOKBACK_DAYS before the corn was
#     planted, or on the planting day (the package ends a stand with no end date at
#     the next planting), including cover crops; with several, the largest credit is used
#     (credits are not added)
#   - second-year credit: corn after corn that itself followed alfalfa or a
#     perennial legume mix gets the "alfalfa_2nd" credit
#   - perennial_mix (e.g. the Wooster meadow, recorded as "perennial mix") is
#     credited at half the alfalfa values, assuming a legume-grass stand; check
################################################################################

LB_TO_KG <- 1.12085
LEGUME_LOOKBACK_DAYS <- 548
LEGUME_CREDIT_LB <- tibble::tribble(
  ~legume_class,   ~low, ~default, ~high,   ~description,
  "alfalfa",        100,    150,     190,   "alfalfa or alfalfa mix stand, first corn year",
  "alfalfa_2nd",      0,     50,      60,   "second corn year after alfalfa or a perennial legume mix",
  "perennial_mix",   50,     75,      95,   "perennial legume-grass mix (half the alfalfa credit)",
  "red_clover",      40,     80,     100,   "red clover or clover mix (plowed down, frost-seeded or cover)",
  "hairy_vetch",     40,     70,     100,   "hairy vetch cover crop",
  "other_legume",    20,     40,      60,   "other clovers, peas, legume-containing cover or forage mixes",
  "soybean",          0,     40,      50,   "soybean (or dry beans)")

# legume class of a crop name (NA if not a legume)
legume_class <- function(name) {
  x <- tolower(trimws(as.character(name)))
  dplyr::case_when(
    grepl("alfalfa", x)                                   ~ "alfalfa",
    grepl("perennial mix|meadow", x)                      ~ "perennial_mix",
    grepl("vetch", x)                                     ~ "hairy_vetch",
    grepl("red clover|clover, red|clover mix|sweet ?clover|yellow.?blossom", x) ~ "red_clover",
    grepl("soy|beans, dry|dry bean", x)                   ~ "soybean",
    grepl("clover|pea|legume|trefoil|lentil|lupin", x)    ~ "other_legume",
    TRUE                                                  ~ NA_character_)
}

# credit class for each corn crop (MGT_combo x year of planting), from the
# package's crop windows (inputs$crop: MGT_combo, CD_name, crop_start, crop_end)
legume_credit_class <- function(crop) {
  ep <- crop %>%
    dplyr::transmute(MGT_combo, name = tolower(trimws(as.character(CD_name))),
                     crop_start = as.Date(crop_start), crop_end = as.Date(crop_end)) %>%
    dplyr::mutate(cls = legume_class(name), is_corn = grepl("corn|maize", name))
  cn <- ep %>% dplyr::filter(is_corn) %>%
    dplyr::mutate(year = as.integer(format(crop_start, "%Y"))) %>%
    dplyr::group_by(MGT_combo, year) %>%
    dplyr::summarise(plant = min(crop_start), end = max(crop_end), .groups = "drop") %>%
    dplyr::arrange(MGT_combo, plant) %>%
    dplyr::group_by(MGT_combo) %>%
    dplyr::mutate(prev_corn_end = dplyr::lag(end)) %>%     # a legume before an earlier corn crop
    dplyr::ungroup()                                        # is that corn's, not this one's
  rank <- setNames(LEGUME_CREDIT_LB$default, LEGUME_CREDIT_LB$legume_class)
  first <- cn %>%
    dplyr::inner_join(dplyr::filter(ep, !is.na(cls)), by = "MGT_combo", relationship = "many-to-many") %>%
    dplyr::filter(crop_end <= plant, crop_end >= plant - LEGUME_LOOKBACK_DAYS,   # <=: a stand with no end
                                                                                  # date ends at the corn planting
                  is.na(prev_corn_end) | crop_end > prev_corn_end) %>%
    dplyr::group_by(MGT_combo, year) %>%
    dplyr::summarise(legume_class = cls[which.max(rank[cls])], .groups = "drop")
  out <- cn %>% dplyr::select(MGT_combo, year, plant) %>% dplyr::left_join(first, by = c("MGT_combo", "year")) %>%
    dplyr::arrange(MGT_combo, plant) %>%
    dplyr::group_by(MGT_combo) %>%
    dplyr::mutate(prev_corn_class = dplyr::lag(legume_class),
                  prev_corn_plant = dplyr::lag(plant)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(legume_class = dplyr::if_else(
      is.na(legume_class) & prev_corn_class %in% c("alfalfa", "perennial_mix") &
        !is.na(prev_corn_plant) & plant - prev_corn_plant <= 400,
      "alfalfa_2nd", legume_class))
  dplyr::select(out, MGT_combo, year, legume_class)
}

# credit in kg N/ha for a vector of classes and a scenario ("none", "low", "default", "high")
legume_credit_kg <- function(cls, scenario = "default") {
  if (scenario == "none") return(rep(0, length(cls)))
  v <- setNames(LEGUME_CREDIT_LB[[scenario]], LEGUME_CREDIT_LB$legume_class)
  out <- unname(v[cls]); out[is.na(out)] <- 0
  out * LB_TO_KG
}
