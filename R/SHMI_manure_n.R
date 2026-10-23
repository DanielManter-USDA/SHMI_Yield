################################################################################
# SHMI_manure_n.R - plant-available N from manure and other organic amendments
# Sourced by SHMI_yield_clean.R and SHMI_yield_onfarm.R.
#
# The input workbooks record organic amendments as TOTAL N (e.g. WICST liquid dairy
# manure: gal/ac x 24 lb total N/1000 gal). Only part of it is available to the crop
# it is applied for, so the yield analyses count organic N at its first-year
# availability (fertilizer N is counted in full):
#   N used = fertilizer N + organic total N x availability + legume credit
# Values are first-year availability as a share of total N, in the range of
# Land-Grant guidance (e.g. UW A2809, Nutrient application guidelines for field,
# vegetable and fruit crops in Wisconsin); confirm them for your region. Low and
# high are for the sensitivity analysis; "none" counts no manure N (fertilizer only);
# "total" counts organic N in full in the year applied, which for systems manured
# every year approximates steady-state availability (this year's share plus carry-over
# from earlier applications). Carry-over is not counted separately in any scenario.
################################################################################

MANURE_AVAIL <- tibble::tribble(
  ~manure_class,  ~none, ~low, ~default, ~high, ~description,
  "liquid_dairy",  0, 0.30,  0.40,    0.50,  "dairy slurry or liquid dairy manure",
  "solid_dairy",   0, 0.20,  0.30,    0.40,  "solid or semi-solid dairy or beef manure",
  "swine_liquid",  0, 0.40,  0.50,    0.60,  "liquid swine manure",
  "poultry",       0, 0.40,  0.50,    0.60,  "poultry manure or litter (not composted)",
  "compost",       0, 0.10,  0.20,    0.30,  "composted manure",
  "other_organic", 0, 0.20,  0.30,    0.40,  "other or unknown organic amendment")

# manure class of an amendment source name. Site-specific names first: at WICST,
# "Dairy" is liquid dairy (gal/ac), "manure" semi-solid dairy (ton/ac) and "CPM"
# composted poultry manure (see the SA_notes in SHMI_WICST_ARS.xlsx)
manure_class <- function(source) {
  x <- tolower(trimws(as.character(source)))
  dplyr::case_when(
    x == "cpm" | grepl("compost", x)                        ~ "compost",
    x == "dairy" | grepl("slurry|liquid", x) & !grepl("swine|hog|pig", x) ~ "liquid_dairy",
    grepl("swine|hog|pig", x)                               ~ "swine_liquid",
    grepl("poultry|chicken|turkey|litter", x)               ~ "poultry",
    x == "manure" | grepl("dairy|beef|cattle|solid", x)     ~ "solid_dairy",
    TRUE                                                    ~ "other_organic")
}

# availability (share of total N) for a vector of sources and a scenario
# ("none", "low", "default", "high" or "total")
manure_avail <- function(source, scenario = "default") {
  if (scenario == "total") return(rep(1, length(source)))
  v <- setNames(MANURE_AVAIL[[scenario]], MANURE_AVAIL$manure_class)
  unname(v[manure_class(source)])
}

# N for each corn crop (MGT_combo x crop year), from the package's applications
# (inputs$n_applications: each application assigned to the first harvest after it).
# N_rec = fertilizer N + available organic N for `scenario`; N_rec_<scenario> for
# every scenario (sensitivity) and N_organic_total (organic N as total N). A crop
# with no convertible N record is absent (NA), as before.
MANURE_SCENARIOS <- c("none", "low", "default", "high", "total")   # every scenario must be listed here
corn_n_available <- function(inputs, scenario = "default") {
  ap <- inputs$n_applications
  if (is.null(ap) || !"crop_year" %in% names(ap))
    stop("N applications by crop (n_applications$crop_year) not found: update the SHMI package and use calc_n_rate = TRUE")
  if (!"SA_cat" %in% names(ap)) ap$SA_cat <- NA_character_
  if (!"SA_source" %in% names(ap)) ap$SA_source <- NA_character_
  ap <- ap %>%
    dplyr::filter(!is.na(crop_year), grepl("corn|maize", crop, ignore.case = TRUE)) %>%
    dplyr::mutate(organic = grepl("organic|manure", SA_cat, ignore.case = TRUE))
  sum_n <- function(x) if (all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE)
  out <- ap %>% dplyr::group_by(MGT_combo, year = as.integer(crop_year)) %>%
    dplyr::summarise(N_organic_total = sum(N_kg_ha[organic], na.rm = TRUE), .groups = "drop")
  for (sc in MANURE_SCENARIOS) {
    s <- ap %>%
      dplyr::mutate(N_av = N_kg_ha * ifelse(organic, manure_avail(SA_source, sc), 1)) %>%
      dplyr::group_by(MGT_combo, year = as.integer(crop_year)) %>%
      dplyr::summarise(!!paste0("N_rec_", sc) := sum_n(N_av), .groups = "drop")
    out <- dplyr::left_join(out, s, by = c("MGT_combo", "year"))
  }
  if (!paste0("N_rec_", scenario) %in% names(out))
    stop("Manure scenario '", scenario, "' is not in MANURE_SCENARIOS")
  out$N_rec <- out[[paste0("N_rec_", scenario)]]
  out
}
