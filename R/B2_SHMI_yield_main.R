################################################################################
# B2_SHMI_yield_main.R
# Does corn yield rise with the Soil Health Management Index (SHMI)?
#
# Within-site analysis of ten long-term experiments (Moore & Manter, in review):
#   Fig. 3a   within sites, treatments with higher SHMI yield more than other
#             treatments at the same site in the same year (panel Fig1_within_sites.png)
#   Fig. 3b   SHMI x N: the within-site slope against N, from 50 kg N/ha (Fig2_SHMI_x_N.png)
#   Fig. 3c   SHMI x system age: the slope at 1, 10 and 30 years (Fig3_SHMI_x_age.png)
#   Extended Data Fig. 4 and Supplementary Data Table 1 (section 7b): each
#             experiment's SHMI range and own slope, the random-slope average,
#             contrast groups, heterogeneity and SHMI-scale checks
#   Supplementary Data Fig. 3: effect sizes of SHMI, N and the other treatments,
#             each over its typical range (Fig_effect_sizes)
# Panels are assembled by B7_Assemble_manuscript_figures.R. Year conditions,
# drought and worst years are in B3_SHMI_yield_stability.R.
#
# WHAT SHMI IS (SHMI R package, version 2.0.0)
#   SHMI_base = 0.687 Cover/Roots + 0.150 Diversity + 0.163 Inverse Disturbance
#   SHMI      = SHMI_base + 0.60 x p_amend x (100 - SHMI_base)
#   (p_amend = share of years with an organic amendment). 0-100, scored from
#   management records only, so it is independent of the yields tested here.
#
# DATA
#   DRIVES (built by hand): OHTVD_Ho 1964-2012, OHTVD_Wo 1962-2021, WICST_ARS 1990-2020.
#   AgCROS (converted by python/agcros2shmi/agcros2shmi.py): NELITCSE 1986-2014,
#     NEMLTCRS 1983-2006 (before its 2007 change to no-till), NEMEIRR 2001-2015,
#     NECCIRR 2010-2018, and (INCLUDE_AGCROS = TRUE) INWLACRE 2004-2009, IAAM7071
#     2008-2012 and PAUP 2008-2012, whose treatments are grouped into systems from
#     the AgCROS descriptors (SHMI_<site>_units.csv). MNMOFS and INWLTPAC were
#     screened out (see AGCROS_SITES); MNSPReap is left out (its three locations
#     share one set of coordinates in AgCROS).
#   Corn years only. County yields: USDA-NASS county corn grain yield, same county
#   and year (irrigated county yield at the two irrigated sites; agricultural-district
#   yield where the county is not published). Yields: dry matter, kg/ha
#   (bu/ac x 25.40117 x 2.471054 x 0.845).
#
# NITROGEN (N, kg/ha, for each corn crop)
#   N = fertilizer N + plant-available manure N + legume credit
#   - each application counts toward the first harvest after it, so fall
#     applications go to the next crop (SHMI package, n_applications)
#   - manure and other organic N (recorded as total N) at its first-year
#     availability: liquid dairy 40%, solid dairy 30%, swine 50%, compost 20%
#     (SHMI_manure_n.R)
#   - legume credit for the crop before the corn: soybean 40 lb N/ac, alfalfa 150
#     in the first corn year, etc. (SHMI_legume_credit.R)
#   - no N record: 0 at NEMLTCRS, WICST and NECCIRR (zero-N plots, organic
#     systems); the Ohio sites have no N records and were managed at the optimum,
#     entered as 180 lb N/ac; NELITCSE and NEMEIRR fill a gap with that
#     site-year's median
#
# MODEL (all corn plot-years)
#   yield ~ site x year + N per site + other treatment factors
#           + SHMI_w + SHMI_w:N + SHMI_w:age + (1 | treatment/plot), AR(1) over years
#   SHMI_w = SHMI minus its site mean, so the slope compares treatments at the same
#   site; site x year removes the weather, county yield, soil, the site's yield level
#   and long-term trends. The slope is reported at 150 kg N/ha and the average system
#   age. Other treatment factors (not scored by SHMI): residue removal (NEMEIRR,
#   IAAM7071, PAUP), irrigation x residue removal (NECCIRR), the modified no-till
#   planter (INWLACRE) and organic management (WICST CS3, CS5), so the organic
#   systems' yield penalty is not charged to SHMI. MEI_SLOPE_BY_RESIDUE = TRUE lets
#   NEMEIRR's slope differ by residue level (off in the global model).
#   Plot-years with a scaled residual beyond +/- 4 in their own model are set aside
#   as single-plot failures (excluded_plot_years.csv).
#
# ROBUSTNESS (section 6b): SHMI vs a full treatment factor (AIC; SHMI_vs_treatment.csv)
#   and the slope refitted leaving out each experiment in turn (results_main.csv).
#
# FIGURES: y = corn yield minus the county average (kg/ha). Figures are saved
#   without titles or captions; options(shmi.ms_figures = FALSE) keeps them.
#
# RUN: from the project (here::here()); SHMI_legume_credit.R and SHMI_manure_n.R
#   are sourced from R/. NASS yields are read from data/cache (a NASS key is needed
#   only to download them again: NASS_API_KEY=... in ~/.Renviron).
#
# Note: an earlier version also compared Minnesota farms (Discovery Farms MN);
#   that comparison is not in the paper. Its settings (CLIMATE_IN_MAIN, FARM_N)
#   and NCCPI helpers remain in the code but are not used for the paper's results.
################################################################################
library(SHMI)        # SHMI scores, crop windows, N by crop
library(dplyr)
library(tidyr)
library(purrr)
library(nlme)        # mixed models with AR(1) (research sites)
library(lme4)        # mixed models (farms)
library(lmerTest)    # p-values for lme4 models (Satterthwaite df)
library(metafor)     # pooling the site slopes (Fig 4)
library(ggplot2)
library(ggh4x)       # fixed panel sizes
library(readxl)
library(sf)
library(soilDB)
library(daymetr)

setwd(here::here())                 # the project root (SHMI_Yield.Rproj)
OUT   <- "Outputs/yield_main_default_total"     # this script's figures and tables
# To rerun with other SHMI weights at the research sites (R/B5b_SHMI_yield_weights_sensitivity.R),
# set before sourcing: SHMI_WEIGHTS_OVERRIDE <- list(weights = c(Cover =, Diversity =, InvDist =),
# growing_share = ) and OUT_OVERRIDE <- "<folder>".
if (exists("OUT_OVERRIDE")) OUT <- OUT_OVERRIDE
# ---- AgCROS sites converted with agcros2shmi.py -------------------------------------
INCLUDE_AGCROS <- TRUE          # TRUE: add the AgCROS sites below to the research sites
if (exists("INCLUDE_AGCROS_OVERRIDE")) INCLUDE_AGCROS <- INCLUDE_AGCROS_OVERRIDE
AGCROS_SITES <- c("INWLACRE", "IAAM7071", "PAUP")      # the 10-site analysis
if (exists("AGCROS_SITES_OVERRIDE")) AGCROS_SITES <- AGCROS_SITES_OVERRIDE
# Screened out after conversion (sensitivity runs only, add back here to include):
#   MNMOFS   - fertilizer N missing in 2008-09 and only starter N recorded on many plots;
#              organic plots have no manure N in AgCROS.
#   INWLTPAC - split-N plots become reduced-N plots in 2008, N records incomplete (split N
#              records one of two applications, modified planter about half the design N).
#              INWLACRE shares the 2008 redesign (split N -> reduced N, modified planter ->
#              intermittent tillage); AGCROS_REFINE below handles both sites alike.
# AGCROS_REFINE (FALSE = the 10-site analysis as run): TRUE ends the West Lafayette period
# in 2007 (before the redesign), gives every corn plot there the site-year median N (all
# treatments get the same design N), adds biochar (IAAM7071) and twin rows (PAUP) as
# "other" factors, and matches the split-N "other" factor on the N descriptor (with FALSE
# it is matched on the treatment descriptor only, where "split" does not appear).
AGCROS_REFINE <- FALSE
if (exists("AGCROS_REFINE_OVERRIDE")) AGCROS_REFINE <- AGCROS_REFINE_OVERRIDE
AGCROS_END       <- if (AGCROS_REFINE) list(INWLACRE = 2007, INWLTPAC = 2007) else list()
AGCROS_UNIFORM_N <- if (AGCROS_REFINE) c("INWLACRE", "INWLTPAC") else character()
# Sites whose organic systems are left out of the yield analysis: at MNMOFS half the corn
# plots are organic, their manure has no N in AgCROS (manure N unknown), and organic
# management (weeds, N timing) is confounded with the amendment bonus that raises their
# SHMI. Conventional plots keep the 2- vs 4-year rotation x tillage x N contrasts.
AGCROS_CONVENTIONAL_ONLY <- c("MNMOFS")
# Treatment factors at AgCROS sites that SHMI does not score, read from the AgCROS
# treatment descriptors (pattern -> "other" level), in addition to residue removal and
# organic management. West Lafayette: split N application (NS) and the no-till planter
# fitted with row cleaners (NM); N rate itself is in the N curve. Ames: biochar
# (Charcoal1, Charcoal2). University Park: twin-row planting.
AGCROS_OTHER <- list(
  INWLACRE = c("split application" = "splitN", "modified planter" = "planter"),
  INWLTPAC = c("split application" = "splitN", "modified planter" = "planter"))
if (AGCROS_REFINE) AGCROS_OTHER <- c(AGCROS_OTHER, list(IAAM7071 = c("charcoal" = "biochar"),
                                                        PAUP     = c("twin row" = "twinrow")))
AGCROS_DIR   <- "Inputs/yield_workbooks"   # SHMI_<site>.xlsx and SHMI_<site>_units.csv
if (INCLUDE_AGCROS && !grepl("_agcros$", OUT)) OUT <- paste0(OUT, "_agcros")   # keep the 7-site results (also for OUT_OVERRIDE)
CACHE <- "Outputs/yield"          # cached downloads (shared with the other scripts)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)

source(here::here("R", "SHMI_legume_credit.R"))    # legume N credits
source(here::here("R", "SHMI_manure_n.R"))         # plant-available manure N
LEGUME_SCENARIO <- "default"      # legume credits: none / low / default / high
MANURE_SCENARIO <- "total"        # manure availability: none / low / default / high / total (total = manuscript)
# To run under other scenarios from another script, set before sourcing, e.g.
#   SCENARIO_OVERRIDE <- c(legume = "high", manure = "total")
# (and rm(SCENARIO_OVERRIDE) afterwards). Editing the two lines above also works.
if (exists("SCENARIO_OVERRIDE")) {
  LEGUME_SCENARIO <- SCENARIO_OVERRIDE[["legume"]]; MANURE_SCENARIO <- SCENARIO_OVERRIDE[["manure"]]
}
stopifnot(LEGUME_SCENARIO %in% c("none", "low", "default", "high"),
          MANURE_SCENARIO %in% c("none", "low", "default", "high", "total"))
message("N accounting: legume credits = ", LEGUME_SCENARIO, ", manure availability = ", MANURE_SCENARIO)

# ---- settings -------------------------------------------------------------------
N_REF        <- 150     # kg N/ha at which the within-site SHMI slope is reported
AGES         <- c(1, 10, 30)    # system ages drawn in Fig 3 (years)
N_LINES      <- c(50, 100, 150, 200) # N rates marked in Fig 2 (kg N/ha)
N_FIG_MIN    <- 50      # Fig 2 starts here: below it only NEMLTCRS compared SHMI levels
N_SUB        <- c(50, 150)      # N rates of the forest-plot sub-rows (kg N/ha)
CLIMATE_IN_MAIN <- TRUE         # MAT and PPT in the Minnesota farm model (Fig 5); FALSE = sensitivity only
FARM_N          <- TRUE   # N terms in the farm model (Fig 5); FALSE: farms are at optimum N
MEI_SLOPE_BY_RESIDUE <- FALSE   # TRUE: NEMEIRR's SHMI slope may differ by residue removal
CLIM     <- if (CLIMATE_IN_MAIN) c("MAT_c", "PPT100") else character()
ADJ_LAB  <- paste(c("county yield", if (FARM_N) "N", "NCCPI", if (CLIMATE_IN_MAIN) c("MAT", "PPT")), collapse = ", ")
ALT_CLIM <- if (CLIMATE_IN_MAIN) "without MAT and PPT" else "+ MAT and PPT"   # the sensitivity row
RESID_CUT    <- 4       # scaled residual beyond which a plot-year is set aside
MIN_LEVELS   <- 2       # a research site needs at least 2 SHMI levels
OPT_N_KG     <- 180 * 1.12085                        # 180 lb N/ac in kg/ha (Ohio)
BU_TO_KG_DM  <- 25.40117 * 2.471054 * 0.845          # bu/ac (56 lb, 15.5% moisture) -> kg/ha dry matter
ctrl <- lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200)

site_cols <- c(NELITCSE = "#0072B2", NEMLTCRS = "#D55E00", OHTVD_Ho = "#56B4E9", OHTVD_Wo = "#CC79A7",
               WICST_ARS = "#009E73", NEMEIRR = "#E69F00", NECCIRR = "#8C510A")   # Okabe-Ito based; one green
if (INCLUDE_AGCROS)
  site_cols <- c(site_cols, MNMOFS = "#999999", INWLACRE = "#000000", INWLTPAC = "#882255",
                 IAAM7071 = "#F0E442", PAUP = "#6A3D9A")
FARM_COL <- "darkorange"
yr <- function(d) as.integer(format(as.Date(d), "%Y"))


################################################################################
# 1. RESEARCH SITES: SHMI, corn yields and N
################################################################################

# N for each corn crop: fertilizer in full, manure at first-year availability
corn_n <- function(inputs) corn_n_available(inputs, MANURE_SCENARIO)

# SHMI over one period of a site's record (the plots' evaluation window narrowed)
shmi_for_period <- function(inputs, climate, start, end) {
  rb <- inputs$rot_bounds
  rb$rot_start <- pmax(as.Date(rb$rot_start), as.Date(start))
  rb$rot_end   <- pmin(as.Date(rb$rot_end),   as.Date(end))
  rb <- rb[rb$rot_start <= rb$rot_end, ]
  rb$rot_start_yr <- yr(rb$rot_start); rb$rot_end_yr <- yr(rb$rot_end)
  inputs$rot_bounds <- rb
  inputs$mgt <- inputs$mgt[inputs$mgt$MGT_combo %in% rb$MGT_combo, ]
  res <- if (exists("SHMI_WEIGHTS_OVERRIDE")) {
    build_shmi(inputs, climate, weights = SHMI_WEIGHTS_OVERRIDE$weights,
               growing_share = SHMI_WEIGHTS_OVERRIDE$growing_share)
  } else build_shmi(inputs, climate)
  res$indicator_df %>%
    mutate(period_start = yr(start), period_end = yr(end))
}

# annual corn yields joined to SHMI (one row per plot and year); the treatment code
# is split into its parts (rotation, tillage, ...) with trt_regex
corn_yields <- function(inputs, shmi, trt_regex, trt_parts) {
  inputs$yield %>%
    filter(tolower(trimws(CD_name)) == "corn", !is.na(yield_kg_ha)) %>%
    group_by(MGT_combo, CD_name, year) %>%
    summarise(yield_kgha = mean(yield_kg_ha), .groups = "drop") %>%
    inner_join(shmi, by = "MGT_combo", relationship = "many-to-many") %>%
    filter(year >= period_start, year <= period_end) %>%
    tidyr::extract(MGT_trt, into = trt_parts, regex = trt_regex, remove = FALSE) %>%
    mutate(MGT_field = as.character(MGT_field),
           MGT_trt_base = sub("P[0-9]+$", "", MGT_trt))      # treatment without rotation phase
}

# one site: read and check the workbook, climate normals, SHMI, yields, N
run_site_fresh <- function(file, start, end, exclude = NULL, trt_regex, trt_parts) {
  inputs  <- prepare_shmi_inputs(file, start_date_override = start, end_date_override = end,
                                 rotation_window = "events", exclude = exclude,
                                 calc_yield = TRUE, calc_n_rate = TRUE)
  climate <- get_shmi_climate(inputs)
  shmi    <- shmi_for_period(inputs, climate, start, end)
  list(inputs = inputs, climate = climate, shmi = shmi,
       yield = corn_yields(inputs, shmi, trt_regex, trt_parts) %>%
         mutate(year = as.integer(year)) %>%
         left_join(corn_n(inputs), by = c("MGT_combo", "year")) %>%
         left_join(legume_credit_class(inputs$crop), by = c("MGT_combo", "year")))
}

# ---- cached site runs ---------------------------------------------------------------
# A site is scored again only when something that changes its result has changed: the
# workbook (path and modification time), the period, excluded units, the treatment
# regex, the manure and legume scenarios, the SHMI weights override or the SHMI package
# version. Results are kept in memory (SITE_CACHE) and on disk (Outputs/yield/site_runs),
# so a restart or a crash part-way through reuses the sites already done.
# Force a fresh run with SITE_CACHE_USE <- FALSE, or delete Outputs/yield/site_runs.
if (!exists("SITE_CACHE_USE")) SITE_CACHE_USE <- TRUE   # FALSE (set before sourcing) rescores every site
SITE_CACHE_DIR <- file.path(CACHE, "site_runs")
if (!exists("SITE_CACHE", inherits = FALSE) || !is.environment(SITE_CACHE)) SITE_CACHE <- new.env()
# also part of the key: when the SHMI package was installed, and the code that scores a
# site (these functions plus the legume and manure helper scripts), so a package fix or a
# change to that code reruns the sites even without a new version number
SHMI_INSTALLED <- as.character(file.mtime(system.file("Meta", "package.rds", package = "SHMI")))
# (removeSource: the same key whether the script is sourced, run line by line or sys.source'd)
code_text <- function(f) deparse(utils::removeSource(f))
SITE_CODE_HASH <- rlang::hash(list(code_text(run_site_fresh), code_text(corn_yields), code_text(shmi_for_period),
                                   code_text(corn_n),
                                   unname(tools::md5sum(here::here("R", c("SHMI_legume_credit.R", "SHMI_manure_n.R"))))))
run_site <- function(file, start, end, exclude = NULL, trt_regex, trt_parts,
                     coarse = exists("SCORE_WB_DIR_OVERRIDE"), alt = exists("SHMI_FN_OVERRIDE")) {
  if (alt) return(run_site_alt(file, start, end, exclude, trt_regex, trt_parts, coarse = coarse))
  if (coarse) return(run_site_coarse(file, start, end, exclude, trt_regex, trt_parts))
  w <- if (exists("SHMI_WEIGHTS_OVERRIDE")) paste(unlist(SHMI_WEIGHTS_OVERRIDE), collapse = ",") else "official"
  key <- rlang::hash(list(normalizePath(file, mustWork = TRUE), as.numeric(file.mtime(file)), start, end,
                          sort(exclude), trt_regex, trt_parts, MANURE_SCENARIO, LEGUME_SCENARIO, w,
                          as.character(utils::packageVersion("SHMI")),
                          SHMI_INSTALLED,          # a reinstall of the package (same version number)
                          SITE_CODE_HASH))         # the scoring code in this script and its helpers
  rds <- file.path(SITE_CACHE_DIR, paste0(tools::file_path_sans_ext(basename(file)), "_", key, ".rds"))
  if (SITE_CACHE_USE) {
    if (!is.null(SITE_CACHE[[key]])) { message("Using cached run (memory): ", basename(file)); return(SITE_CACHE[[key]]) }
    if (file.exists(rds)) {
      message("Using cached run (disk): ", basename(file))
      SITE_CACHE[[key]] <- readRDS(rds); return(SITE_CACHE[[key]])
    }
  }
  r <- run_site_fresh(file, start, end, exclude, trt_regex, trt_parts)
  SITE_CACHE[[key]] <- r
  dir.create(SITE_CACHE_DIR, recursive = TRUE, showWarnings = FALSE)
  saveRDS(r, rds)
  r
}

# Coarse-record test (R/B9_SHMI_coarse_records.R): with SCORE_WB_DIR_OVERRIDE set, SHMI is
# scored from the degraded copy of each workbook in that folder (same file name), while
# yields, N, climate and everything else come from the full record (the cached full run).
# Not part of SITE_CODE_HASH, so adding it leaves the full-record cache valid.
run_site_coarse <- function(file, start, end, exclude, trt_regex, trt_parts) {
  f2 <- file.path(SCORE_WB_DIR_OVERRIDE, basename(file))
  if (!file.exists(f2)) stop("No degraded workbook for ", basename(file), " in ", SCORE_WB_DIR_OVERRIDE, call. = FALSE)
  full <- run_site(file, start, end, exclude, trt_regex, trt_parts, coarse = FALSE, alt = FALSE)
  w <- if (exists("SHMI_WEIGHTS_OVERRIDE")) paste(unlist(SHMI_WEIGHTS_OVERRIDE), collapse = ",") else "official"
  key <- rlang::hash(list("coarse", normalizePath(f2, mustWork = TRUE), as.numeric(file.mtime(f2)),
                          normalizePath(file, mustWork = TRUE), as.numeric(file.mtime(file)), start, end,
                          sort(exclude), trt_regex, trt_parts, MANURE_SCENARIO, LEGUME_SCENARIO, w,
                          as.character(utils::packageVersion("SHMI")), SHMI_INSTALLED, SITE_CODE_HASH))
  if (!is.null(SITE_CACHE[[key]])) return(SITE_CACHE[[key]])
  message("Scoring the coarse record: ", f2)
  inp2 <- prepare_shmi_inputs(f2, start_date_override = start, end_date_override = end,
                              rotation_window = "events", exclude = exclude,
                              calc_yield = TRUE, calc_n_rate = TRUE)
  shmi <- shmi_for_period(inp2, full$climate, start, end)
  miss <- setdiff(unique(full$shmi$MGT_combo), shmi$MGT_combo)
  if (length(miss)) stop("Coarse record lost units: ", paste(head(miss, 5), collapse = ", "), call. = FALSE)
  r <- list(inputs = full$inputs, climate = full$climate, shmi = shmi,
            yield = corn_yields(full$inputs, shmi, trt_regex, trt_parts) %>%
              mutate(year = as.integer(year)) %>%
              left_join(corn_n(full$inputs), by = c("MGT_combo", "year")) %>%
              left_join(legume_credit_class(full$inputs$crop), by = c("MGT_combo", "year")))
  SITE_CACHE[[key]] <- r
  r
}

# Another score in place of SHMI (used by an earlier comparison script, not in this repository): SHMI_FN_OVERRIDE(inputs, climate,
# start, end) returns one row per unit with MGT_combo and SHMI (0-100); yields, N and
# everything else come from the full run. With SCORE_WB_DIR_OVERRIDE also set, the score
# is computed from the degraded workbook's records. Not part of SITE_CODE_HASH.
run_site_alt <- function(file, start, end, exclude, trt_regex, trt_parts, coarse = FALSE) {
  full <- run_site(file, start, end, exclude, trt_regex, trt_parts, coarse = FALSE, alt = FALSE)
  inp_score <- full$inputs
  if (coarse) {
    f2 <- file.path(SCORE_WB_DIR_OVERRIDE, basename(file))
    if (!file.exists(f2)) stop("No degraded workbook for ", basename(file), " in ", SCORE_WB_DIR_OVERRIDE, call. = FALSE)
    message("Scoring the degraded record: ", f2)
    inp_score <- prepare_shmi_inputs(f2, start_date_override = start, end_date_override = end,
                                     rotation_window = "events", exclude = exclude,
                                     calc_yield = TRUE, calc_n_rate = TRUE)
  }
  alt <- SHMI_FN_OVERRIDE(inp_score, full$climate, start, end)
  miss <- setdiff(unique(full$shmi$MGT_combo), alt$MGT_combo)
  # keep the unit metadata and period of the full run (corn_yields() needs them); swap the score columns
  shmi <- full$shmi %>% select(-any_of(setdiff(names(alt), "MGT_combo"))) %>% left_join(alt, by = "MGT_combo")
  if (length(miss)) stop("SHMI_FN_OVERRIDE lost units: ", paste(head(miss, 5), collapse = ", "), call. = FALSE)
  list(inputs = full$inputs, climate = full$climate, shmi = shmi,
       yield = corn_yields(full$inputs, shmi, trt_regex, trt_parts) %>%
         mutate(year = as.integer(year)) %>%
         left_join(corn_n(full$inputs), by = c("MGT_combo", "year")) %>%
         left_join(legume_credit_class(full$inputs$crop), by = c("MGT_combo", "year")))
}

# ---- the seven sites ------------------------------------------------------------
# The four Nebraska workbooks (NELITCSE, NEMLTCRS, NEMEIRR, NECCIRR) are built from the
# AgCROS tables by agcros2shmi.py (run_conversion.sh; NECCIRR's fertilizer and manure
# from supplements/NECCIRR_amendments.csv), so their treatment codes are AgCROS codes.
# Ohio and WICST come from DRIVES (hand-built workbooks).

# NELITCSE: era 1 (1986-2014). AgCROS codes = tillage (CH chisel, DK disk, MP moldboard,
# NT no-till, RT ridge, SS subsoil) + rotation (CC continuous corn, CSB / SBC the two
# phases of corn-soybean, SB continuous soybean); the 1981-85 baseline (CC, RTCC ...) and
# the 2015 redesign fall outside the period. The 600 block lacks harvests and yields in
# 1986-1995 and is excluded, as before (EXCLUDE_NEL_600).
EXCLUDE_NEL_600 <- TRUE
NEL_FILE <- "Inputs/yield_workbooks/SHMI_NELITCSE.xlsx"
nel_mgt  <- readxl::read_excel(NEL_FILE, sheet = "Mgt_Unit", skip = 3)
exclude.nel <- if (EXCLUDE_NEL_600)
  nel_mgt$MGT_combo[!is.na(nel_mgt$MGT_field) & suppressWarnings(as.integer(nel_mgt$MGT_field)) %in% 601:624] else NULL
nel <- run_site(NEL_FILE, "1986-01-01", "2014-12-31", exclude = exclude.nel,
                trt_regex = "^(CH|DK|MP|NT|RT|SS)(CC|CSB|SBC|SB)$", trt_parts = c("tillage", "rotation"))
nel$yield <- nel$yield %>%       # treatment = tillage x rotation (both phases of corn-soybean together)
  filter(!is.na(tillage)) %>%
  mutate(MGT_trt_base = paste0(tillage, ifelse(rotation %in% c("CSB", "SBC"), "CS", rotation)))
# NEMLTCRS: first era only (conventional tillage 1983-2006). AgCROS codes ROT<r><n>:
# rotation 1-7 x N level 0-2; rotation phases are separate plots of one treatment
nem <- run_site("Inputs/yield_workbooks/SHMI_NEMLTCRS.xlsx", "1983-01-01", "2006-12-31",
                trt_regex = "^ROT([0-9])([0-9])$", trt_parts = c("rotation", "N_code"))
# Ohio: Hoytville and Wooster (Wooster follows the Hoytville tillage protocol)
oh  <- run_site("Inputs/yield_workbooks/SHMI_OHTVD_Ho.xlsx", "1964-01-01", "2012-12-31",
                trt_regex = "(R[0-9]+)(T[0-9]+)(P[0-9]+)", trt_parts = c("rotation", "tillage", "phase"))
wo  <- run_site("Inputs/yield_workbooks/SHMI_OHTVD_Wo.xlsx", "1962-01-01", "2021-12-31",
                trt_regex = "(R[0-9]+)(T[0-9]+)(P[0-9]+)", trt_parts = c("rotation", "tillage", "phase"))
# WICST: cropping systems CS1-CS6 (CS3 and CS5 organic)
wi  <- run_site("Inputs/yield_workbooks/SHMI_WICST_ARS.xlsx", "1990-01-01", "2020-12-31",
                trt_regex = "(CS[0-9]+)(P[0-9]+)", trt_parts = c("rotation", "phase"))
# NEMEIRR: continuous corn, no-till (NT) vs disk (DK) x 0/50/100% residue removal
mei <- run_site("Inputs/yield_workbooks/SHMI_NEMEIRR.xlsx", "2001-01-01", "2015-12-31",
                trt_regex = "(NT|DK)([0-9]+)", trt_parts = c("tillage", "residue"))
# NECCIRR: irrigation (Full/Dfct) x amelioration (X none, C rye cover, M manure) x N x residue
cci <- run_site("Inputs/yield_workbooks/SHMI_NECCIRR.xlsx", "2010-01-01", "2018-12-31",
                trt_regex = "(Full|Dfct)([XCM])([0-9]+)([NA])",
                trt_parts = c("irrigation", "amelioration", "N_rate", "residue"))

runs <- list(NELITCSE = nel, NEMLTCRS = nem, OHTVD_Ho = oh, OHTVD_Wo = wo,
             WICST_ARS = wi, NEMEIRR = mei, NECCIRR = cci)

# ---- AgCROS sites (INCLUDE_AGCROS) ------------------------------------------------
# Period: the site's management record (first to last record year, from the units
# file). Treatment = system: AgCROS treatments that differ only by rotation phase share
# rotation, tillage, N, cover crop, residue and organic descriptors, so MGT_trt_base is
# built from those. Not scored by SHMI, so entered as "other": residue removal
# (Partial / Yes) and organic management.
agcros_runs <- list(); agcros_yield <- NULL; agcros_start <- integer()
if (INCLUDE_AGCROS) {
  # Period = the years with complete management records: planting, tillage (at sites
  # that till) and fertilizer N. AgCROS sometimes holds yields for years without them
  # (West Lafayette: yields 2003-2011, management 2004-2009); in such years tillage would
  # look absent (raising SHMI) and N would be unknown, so they are left out.
  mgt_years <- function(file) {
    rd <- function(sh) readxl::read_excel(file, sheet = sh, skip = 3)
    yrs <- function(x) unique(stats::na.omit(as.integer(format(as.Date(x), "%Y"))))
    cd <- rd("Crop_Diversity"); sd <- rd("Soil_Disturbance"); sa <- rd("Soil_Amendments")
    y  <- yrs(cd$CD_plant_date)
    ty <- yrs(sd$SD_date[!grepl("^SPRAYER", sd$SD_equip)])      # the no-till placeholder does not count
    ny <- yrs(sa$SA_date[!is.na(sa$SA_N)])
    if (length(ty)) y <- intersect(y, ty)
    if (length(ny)) y <- intersect(y, ny)
    if (!length(y)) stop(basename(file), ": no year has planting, tillage and N records")
    range(y)
  }
  run_agcros <- function(loc) {
    u <- read.csv(file.path(AGCROS_DIR, paste0("SHMI_", loc, "_units.csv")), stringsAsFactors = FALSE)
    wb <- file.path(AGCROS_DIR, paste0("SHMI_", loc, ".xlsx"))
    py <- mgt_years(wb)
    if (!is.null(AGCROS_END[[loc]])) py[2] <- min(py[2], AGCROS_END[[loc]])
    if (py[1] > py[2]) stop(loc, ": empty period (", py[1], "-", py[2], ")")
    start <- paste0(py[1], "-01-01")
    end   <- paste0(py[2], "-12-31")
    r <- run_site(file.path(AGCROS_DIR, paste0("SHMI_", loc, ".xlsx")), start, end,
                  trt_regex = "(.+)", trt_parts = "trt_code")
    # system = the AgCROS descriptors plus the treatment's own description (which separates,
    # e.g., intermittent tillage, split N and the modified planter from plain no-till)
    desc <- u %>% transmute(MGT_combo, treatment_descriptor,
                            system = paste(treatment_descriptor, rotation_descriptor, tillage_descriptor, n_treatment_descriptor,
                                           ifelse(is.na(cover_crop) | cover_crop == "", "no cover", cover_crop),
                                           paste("residue removal", residue_removal),
                                           ifelse(organic_management %in% "Yes", "organic", "conventional"), sep = " | "),
                            residue = residue_removal, organic = organic_management %in% "Yes")
    r$yield <- r$yield %>% left_join(desc, by = "MGT_combo") %>%
      mutate(MGT_trt_base = system,
             other = case_when(organic & residue %in% c("Partial", "Yes") ~ paste0(loc, "_organic_res", residue),
                               organic ~ paste0(loc, "_organic"),
                               residue %in% c("Partial", "Yes") ~ paste0(loc, "_res", residue),
                               TRUE ~ "none"))
    if (!is.null(AGCROS_OTHER[[loc]])) for (k in names(AGCROS_OTHER[[loc]])) {
      hit <- grepl(k, if (AGCROS_REFINE) r$yield$system else r$yield$treatment_descriptor, ignore.case = TRUE)
      r$yield$other[hit] <- ifelse(r$yield$other[hit] == "none", paste0(loc, "_", AGCROS_OTHER[[loc]][[k]]),
                                   paste0(r$yield$other[hit], "_", AGCROS_OTHER[[loc]][[k]]))
    }
    if (loc %in% AGCROS_CONVENTIONAL_ONLY) {
      n0 <- nrow(r$yield)
      r$yield <- r$yield %>% filter(!organic)
      message(loc, ": organic systems left out (", n0 - nrow(r$yield), " corn plot-years)")
    }
    r$start_year <- as.integer(substr(start, 1, 4))
    message(loc, ": ", n_distinct(r$yield$MGT_trt_base), " systems, ", nrow(r$yield), " corn plot-years, ",
            start, " to ", end)
    r
  }
  agcros_runs  <- setNames(lapply(AGCROS_SITES, run_agcros), AGCROS_SITES)
  runs         <- c(runs, agcros_runs)
  agcros_yield <- bind_rows(lapply(agcros_runs, `[[`, "yield"))
  agcros_start <- vapply(agcros_runs, `[[`, integer(1), "start_year")
}

# Treatment factors that SHMI does not score enter as "other" (baseline "none" at
# every site, so the term only compares plots within a site): NEMEIRR residue
# removal; NECCIRR irrigation level x residue removal; WICST organic systems (CS3, CS5).
# Ohio (Triplett-Van Doren experiment): established at Wooster in 1962 and at Hoytville in
# 1963 (Dick & Van Doren 1985). The first years had many compromised stands while equipment
# and field operations were refined, and are excluded from the yield analysis as in
# de Camargo Santos et al. (2025, Field Crops Res. 331, 109993), who analyse Wooster from
# 1964 and Hoytville from 1965. Management in those years still counts toward SHMI.
OHIO_YIELD_FROM <- c(OHTVD_Wo = 1964, OHTVD_Ho = 1965)
# Hoytville 1981: 10.75" of rain in June and cutworms; many plots replanted "to preserve
# history", yields incomplete, and NT corn "essentially eliminated from consideration"
# (Hoytville field notes, Dick, Triplett & Van Doren 2013, OSU Knowledge Bank).
OHIO_HO_DROP <- 1981
site_start <- c(NELITCSE = 1986, NEMLTCRS = 1983, OHTVD_Ho = 1963, OHTVD_Wo = 1962,
                WICST_ARS = 1990, NEMEIRR = 2001, NECCIRR = 2010, agcros_start)
corn <- bind_rows(
  nel$yield %>% mutate(other = "none"),
  nem$yield %>% mutate(other = "none") %>% filter(year <= 2006),
  oh$yield  %>% mutate(other = "none") %>% filter(year >= OHIO_YIELD_FROM[["OHTVD_Ho"]], !year %in% OHIO_HO_DROP),
  wo$yield  %>% mutate(other = "none") %>% filter(year >= OHIO_YIELD_FROM[["OHTVD_Wo"]]),
  wi$yield  %>% mutate(other = ifelse(rotation %in% c("CS3", "CS5"), "WICST_organic", "none")),
  mei$yield %>% mutate(other = ifelse(residue == "0", "none", paste0("NEMEIRR_res", residue))),
  cci$yield %>% mutate(other = ifelse(irrigation == "Full" & residue == "N", "none",
                                      paste0("NECCIRR_", irrigation, "_res", residue))),
  agcros_yield) %>%                                             # NULL unless INCLUDE_AGCROS
  mutate(site = MGT_farm,
         other = relevel(factor(other), ref = "none"),
         plot_id = factor(paste(site, MGT_field, sep = "_")),   # random effect: the plot
         site_year = factor(paste(site, year, sep = "_")),
         trt = factor(paste(site, MGT_trt_base)),
         age = unname(year - site_start[as.character(site)])) %>% # system age: years since the experiment began
  group_by(site) %>%
  mutate(SHMI_w = SHMI - mean(SHMI[!duplicated(MGT_combo)])) %>%  # SHMI minus the site mean
  ungroup()

# NEMEIRR: tillage is its only SHMI difference, crossed with residue removal. These
# columns let the tillage (SHMI) slope differ by removal level (centred, so SHMI_w
# stays the slope averaged over the levels; 0 at every other site).
mei_trt <- corn %>% filter(site == "NEMEIRR") %>% distinct(MGT_trt_base, residue)
for (r in c("50", "100")) {
  corn[[paste0("mei_res", r)]] <- ifelse(corn$site == "NEMEIRR",
                                         as.numeric(corn$residue %in% r) - mean(mei_trt$residue == r), 0)
}

# ---- N for each corn crop (see NITROGEN in the header) ---------------------------
zero_if_missing    <- c("NEMLTCRS", "WICST_ARS", "NECCIRR")
optimum_if_missing <- c("OHTVD_Ho", "OHTVD_Wo")
# AgCROS sites: N records are missing in some years (MNMOFS 2008-2009, INWLACRE 2003,
# INWLTPAC 2003, 2010-2011). A plot-year without a record takes the median recorded N of
# its own treatment, unless the treatment is unfertilized by design (its AgCROS N
# descriptor says no fertilizer / 0 N), which gets 0.
UNFERTILIZED <- "no fert|none|^0 |\\b0 (kg|N)|zero"
agcros_n <- corn %>% filter(site %in% AGCROS_SITES) %>%
  group_by(site, MGT_trt_base) %>%
  summarise(N_trt = ifelse(grepl(UNFERTILIZED, first(MGT_trt_base), ignore.case = TRUE) | all(is.na(N_rec)),
                           0, median(N_rec, na.rm = TRUE)), .groups = "drop")
corn <- corn %>% select(-any_of("N_trt")) %>%
  left_join(agcros_n, by = c("site", "MGT_trt_base")) %>%
  group_by(site_year) %>%
  # West Lafayette 2004-2007: every corn treatment gets the same N by design (220 lbs/A
  # UAN; split or not), so each plot takes the site-year median of the recorded amounts.
  # AgCROS records less on some plots (modified planter about 106 kg, split N 109 kg: only
  # one of the two splits), which would otherwise read as an N difference.
  mutate(N_rec = ifelse(site %in% AGCROS_UNIFORM_N & any(!is.na(N_rec)), median(N_rec, na.rm = TRUE), N_rec),
         N_kg_rec = case_when(!is.na(N_rec) ~ N_rec,
                              site %in% AGCROS_SITES & !is.na(N_trt) ~ N_trt,
                              site %in% zero_if_missing ~ 0,
                              site %in% optimum_if_missing ~ OPT_N_KG,
                              any(!is.na(N_rec)) ~ median(N_rec, na.rm = TRUE),
                              TRUE ~ 0)) %>%
  ungroup() %>%
  # legume credit, except where N is the assumed optimum (the optimum already
  # allows for the previous crop)
  mutate(N_leg = ifelse(is.na(N_rec) & site %in% optimum_if_missing, 0,
                        legume_credit_kg(legume_class, LEGUME_SCENARIO)),
         N_kg = N_kg_rec + N_leg,
         N100 = N_kg / 100)

# check: N by treatment (mean and highest single year)
n_cov <- corn %>% group_by(site, MGT_trt_base) %>%
  summarise(corn_years = n(), share_N_recorded = mean(!is.na(N_rec)),
            N_fert_manure_mean = mean(N_kg_rec), N_legume_mean = mean(N_leg),
            N_mean = mean(N_kg), N_max_year = max(N_kg), SHMI = mean(SHMI), .groups = "drop")
write.csv(n_cov, file.path(OUT, "N_by_treatment.csv"), row.names = FALSE)

# ---- model terms -------------------------------------------------------------------
# does a variable vary among plots within any site-year? (else site x year absorbs it)
varies <- function(d, v) any(tapply(d[[v]], d$site_year, function(x) diff(range(x))) > 1e-8)
# N response. Every site gets its own linear N term (fertilizer + manure + legume credit);
# it drops out where N does not vary within a site-year (NEMEIRR, OHTVD: one rate).
# A quadratic term is added where N was managed by design at three or more rates over a
# wide range (N_QUAD_SITES: NEMLTCRS 0/90/180; NECCIRR 0/125/200; WICST 39-224 by system).
# Elsewhere N differs only through legume credits (NELITCSE), rotation-adjusted rates
# with equal totals (PAUP), two designed rates (INWLACRE 220 vs 175 lbs/A from 2008) or
# steps of <= 51 kg tied to other treatment factors (IAAM7071: management level, stover
# removal), too little to estimate curvature.
# SHMI x N is estimated only where N was managed by design (N_INTERACT_SITES): the three
# quadratic sites and INWLACRE. B9 checks it without INWLACRE.
N_QUAD_SITES     <- c("NEMLTCRS", "NECCIRR", "WICST_ARS")
N_INTERACT_SITES <- c(N_QUAD_SITES, "INWLACRE")
if (exists("N_QUAD_SITES_OVERRIDE"))     N_QUAD_SITES     <- N_QUAD_SITES_OVERRIDE
if (exists("N_INTERACT_SITES_OVERRIDE")) N_INTERACT_SITES <- N_INTERACT_SITES_OVERRIDE
add_n_terms <- function(d) {
  for (s_ in unique(as.character(d$site))) {
    d[[paste0("N1_", s_)]] <- ifelse(d$site == s_, d$N100, 0)
    if (s_ %in% N_QUAD_SITES) d[[paste0("N2_", s_)]] <- d[[paste0("N1_", s_)]]^2
  }
  d
}
corn <- add_n_terms(corn)
n_terms <- function(d) {
  v <- grep("^N[12]_", names(d), value = TRUE)
  v[vapply(v, function(k) varies(d, k), logical(1))]
}
o_term <- function(d) if (n_distinct(d$other) > 1) "other" else character()
mei_terms <- function(d, x = "SHMI_w") {
  if (!MEI_SLOPE_BY_RESIDUE) return(character())    # one tillage slope at NEMEIRR
  v <- c("mei_res50", "mei_res100")
  v <- v[vapply(v, function(k) k %in% names(d) && any(abs(d[[k]]) > 1e-12), logical(1))]
  if (length(v)) paste0(x, ":", v) else character()
}

# ---- sites with at least MIN_LEVELS SHMI levels --------------------------------------
# a level = a distinct treatment-mean SHMI (whole points); N rate, residue removal
# and irrigation level do not add levels
lv <- corn %>% group_by(site, MGT_trt_base) %>% summarise(s = round(mean(SHMI)), .groups = "drop") %>%
  count(site, name = "n_trt") %>% left_join(
    corn %>% group_by(site, MGT_trt_base) %>% summarise(s = round(mean(SHMI)), .groups = "drop") %>%
      group_by(site) %>% summarise(n_levels = n_distinct(s)), by = "site")
print(lv)
corn <- corn %>% filter(site %in% lv$site[lv$n_levels >= MIN_LEVELS]) %>% droplevels()
# the AR(1) needs one yield per plot and year
if (any(duplicated(corn[c("plot_id", "year")]))) corn$plot_id <- factor(corn$MGT_combo)

# ---- one site's model; residual screening ---------------------------------------------
fit_site <- function(d, extra = NULL, x = "SHMI_w") {
  d <- droplevels(d)
  rhs <- c("site_year", n_terms(d), o_term(d), x, mei_terms(d, x), extra)
  lme(reformulate(rhs, "yield_kgha"), random = ~ 1 | trt/plot_id,
      correlation = corCAR1(form = ~ year | trt/plot_id), data = d, method = "REML", control = ctrl)
}
corn <- split(corn, corn$site) %>%
  map_dfr(function(ds) {
    ds <- arrange(ds, plot_id, year)
    mutate(ds, resid_n = as.numeric(residuals(fit_site(ds), type = "normalized")))
  })
excluded <- corn %>% filter(abs(resid_n) > RESID_CUT) %>%
  transmute(data = "research site", site, unit = MGT_combo, year, SHMI, yield = yield_kgha, resid = round(resid_n, 2))
corn <- corn %>% filter(abs(resid_n) <= RESID_CUT) %>% select(-resid_n) %>% droplevels()
message(nrow(excluded), " research plot-years set aside (|normalized residual| > ", RESID_CUT, ")")

# centred terms for the interactions (Fig 2, 3)
AGE_MEAN <- mean(corn$age)                               # average system age (years)
n_vary <- names(which(sapply(split(corn, corn$site), function(d) varies(d, "N100"))))
n_var_sites   <- intersect(N_INTERACT_SITES, n_vary)  # SHMI x N estimated here (Fig 2)
n_curve_sites <- intersect(N_QUAD_SITES, n_vary)      # N curves for N-equivalent, Fig 2b, N effect size
corn <- corn %>% mutate(age_c = age - AGE_MEAN,
                        N_c = ifelse(site %in% n_var_sites, N100 - N_REF / 100, 0))
message(sprintf("Average system age %.1f years; linear N at: %s; quadratic N at: %s; SHMI x N at: %s", AGE_MEAN,
                paste(n_vary, collapse = ", "), paste(n_curve_sites, collapse = ", "), paste(n_var_sites, collapse = ", ")))
# the site runs (all crops' yields and SHMI) and the corn data, for the soybean test
# (B2b_SHMI_soybean.R); written on every run, read by nothing else
saveRDS(list(runs = runs, site_start = site_start, corn = corn, AGE_MEAN = AGE_MEAN,
             n_var_sites = n_var_sites, AGCROS_DIR = AGCROS_DIR, AGCROS_SITES = AGCROS_SITES,
             AGCROS_CONVENTIONAL_ONLY = AGCROS_CONVENTIONAL_ONLY, EXCLUDE_NEL_600 = EXCLUDE_NEL_600),
        file.path(OUT, "B2_data.rds"))


################################################################################
# 2. COUNTY YIELDS (research sites) AND NCCPI HELPERS
################################################################################
# USDA-NASS county corn yield (bu/ac) for each site's county; irrigated yield for
# the irrigated sites. Cached in Outputs/yield.
SITE_COUNTY <- tribble(
  ~site,       ~state, ~county,     ~practice,
  "NELITCSE",  "NE",   "LANCASTER", "NON-IRRIGATED",  # not irr
  "NEMLTCRS",  "NE",   "SAUNDERS",  "NON-IRRIGATED",  # not irr  (2011 missing)
  "OHTVD_Ho",  "OH",   "WOOD",      "ALL PRODUCTION PRACTICES",  # not irr
  "OHTVD_Wo",  "OH",   "WAYNE",     "ALL PRODUCTION PRACTICES",  # not irr
  "WICST_ARS", "WI",   "COLUMBIA",  "ALL PRODUCTION PRACTICES",  # not irr
  "NEMEIRR",   "NE",   "SAUNDERS",  "IRRIGATED",
  "NECCIRR",   "NE",   "CLAY",      "IRRIGATED")
if (INCLUDE_AGCROS) SITE_COUNTY <- bind_rows(SITE_COUNTY, tribble(
  ~site,       ~state, ~county,      ~practice,
  "MNMOFS",    "MN",   "STEVENS",    "ALL PRODUCTION PRACTICES",
  "INWLACRE",  "IN",   "TIPPECANOE", "ALL PRODUCTION PRACTICES",
  "INWLTPAC",  "IN",   "TIPPECANOE", "ALL PRODUCTION PRACTICES",
  "IAAM7071",  "IA",   "BOONE",      "ALL PRODUCTION PRACTICES",
  "PAUP",      "PA",   "CENTRE",     "ALL PRODUCTION PRACTICES"))

cty_file <- file.path(CACHE, "nass_site_county_corn_yield.csv")
NASS_API_KEY <- Sys.getenv("NASS_API_KEY")   # set NASS_API_KEY=... in ~/.Renviron (never in the script)

# download only the county x practice series not yet in the cache (new sites are added)
cty_cached <- if (file.exists(cty_file)) read.csv(cty_file, stringsAsFactors = FALSE) else NULL
cty_need   <- distinct(SITE_COUNTY, state, county, practice)
if (!is.null(cty_cached))
  cty_need <- anti_join(cty_need, distinct(cty_cached, state, county, practice), by = c("state", "county", "practice"))
if (nrow(cty_need) > 0) {
  message("Downloading NASS county yields: ", paste(cty_need$county, collapse = ", "))
  message("Please check for missing data - SAUNDERS 2011")
  library(rnassqs)
  nassqs_auth(key = NASS_API_KEY)
  get_cty <- function(state, county, practice) {
    q <- list(source_desc = "SURVEY", sector_desc = "CROPS", commodity_desc = "CORN",
              statisticcat_desc = "YIELD", unit_desc = "BU / ACRE", util_practice_desc = "GRAIN",
              prodn_practice_desc = practice, reference_period_desc = "YEAR", domain_desc = "TOTAL",
              agg_level_desc = "COUNTY", state_alpha = state, county_name = county, year__GE = 1955)
    tryCatch(nassqs(q) %>% transmute(year = as.integer(year), bu_ac = as.numeric(gsub(",", "", Value))),
             error = function(e) { message("No NASS data: ", state, " ", county, " ", practice); NULL })
  }
  cty_raw <- pmap_dfr(cty_need, function(state, county, practice) {
    x <- get_cty(state, county, practice)
    if (is.null(x) && practice != "ALL PRODUCTION PRACTICES") {       # fall back to all practices
      message("  using ALL PRODUCTION PRACTICES for ", county, " instead")
      x <- get_cty(state, county, "ALL PRODUCTION PRACTICES")
    }
    if (is.null(x)) return(NULL)
    mutate(x, state = state, county = county, practice = practice)  # keep the requested practice label for the join
  }) %>% filter(!is.na(bu_ac)) %>% distinct(state, county, practice, year, .keep_all = TRUE)
  write.csv(bind_rows(cty_cached, cty_raw), cty_file, row.names = FALSE)
}

# Manually add Saunders 2011, missing from the NASS download (kept only if NASS lacks it)
# TODO Find the source for these...
MANUAL_CTY <- tribble(
  ~year, ~bu_ac, ~state, ~county,    ~practice,
  2011,  133.4,  "NE",   "SAUNDERS", "NON-IRRIGATED",
  2011,  181.7,  "NE",   "SAUNDERS", "IRRIGATED")

cty_exp <- read.csv(cty_file, stringsAsFactors = FALSE) %>%
  mutate(practice = ifelse(state %in% c("OH", "WI"), "ALL PRODUCTION PRACTICES", practice)) %>%
  bind_rows(MANUAL_CTY) %>%
  distinct(state, county, practice, year, .keep_all = TRUE) %>%
  inner_join(SITE_COUNTY, by = c("state", "county", "practice"), relationship = "many-to-many") %>%
  transmute(site, year = as.integer(year), county_kgha = BU_TO_KG_DM * bu_ac)

corn <- corn %>% left_join(cty_exp, by = c("site", "year"))
message(sum(!is.na(corn$county_kgha)), " of ", nrow(corn), " research plot-years have a county yield")

# NCCPI helpers (Minnesota farms): area-weighted NCCPI of the SSURGO map units around a point
NCCPI_RULE <- "NCCPI - National Commodity Crop Productivity Index (Ver 3.0)"
nccpi_of <- function(mukeys) {
  x <- get_SDA_interpretation(rulename = NCCPI_RULE, method = "Dominant Component", mukeys = unique(mukeys))
  rc <- grep("^rating_", names(x), value = TRUE)[1]
  tibble(mukey = as.character(x$mukey), nccpi = as.numeric(x[[rc]]))
}
nccpi_buffer <- function(pts_sf, ids, buffer_m) {          # area-weighted NCCPI around each point
  buf <- st_transform(st_buffer(st_transform(pts_sf, 5070), buffer_m), 4326)
  map_dfr(seq_len(nrow(buf)), function(i) {
    mu <- NULL
    for (k in 1:3) {                                         # Soil Data Access can drop a request
      mu <- tryCatch(SDA_spatialQuery(buf[i, ], what = "mupolygon", geomIntersection = TRUE), error = function(e) NULL)
      if (!is.null(mu) && nrow(mu)) break
      Sys.sleep(2)
    }
    if (is.null(mu) || !nrow(mu)) return(tibble(id = ids[i], nccpi = NA_real_))
    mu$area <- as.numeric(st_area(st_transform(mu, 5070)))
    st_drop_geometry(mu) %>% mutate(mukey = as.character(mukey)) %>%
      left_join(nccpi_of(mu$mukey), by = "mukey") %>%
      summarise(nccpi = weighted.mean(nccpi, area, na.rm = TRUE)) %>% mutate(id = ids[i], .before = 1)
  })
}

write.csv(excluded, file.path(OUT, "excluded_plot_years.csv"), row.names = FALSE)

################################################################################
# 4. FIGURE HELPERS
################################################################################
YLAB  <- "Corn yield adjusted for site, year,\nN and other treatments (kg ha\u207B\u00B9)"   # level = overall mean yield
XLAB  <- "Soil Health Management Index (SHMI)"
PANEL <- force_panelsizes(rows = unit(2.6, "in"), cols = unit(3.6, "in"))
# Manuscript figures: no title, subtitle or caption (the legend is written in the manuscript;
# panel letters are added by B7_Assemble_manuscript_figures.R). options(shmi.ms_figures = FALSE)
# keeps them for working plots.
MS <- if (isTRUE(getOption("shmi.ms_figures", TRUE))) ggplot2::labs(title = NULL, subtitle = NULL, caption = NULL) else NULL
# slope label (top-left) and p-value (bottom-right)
lab_slope <- function(est, lo, hi, what = "Yield change per 10 SHMI")
  sprintf("%s: %+.0f kg/ha (95%% CI %.0f%s%.0f)", what, est, lo, if (lo < 0) " to " else "\u2013", hi)
lab_p <- function(p, what = NULL) paste0(if (!is.null(what)) paste0(what, ": "), ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p)))
corner_labels <- function(top, bottom) list(
  annotate("text", x = -Inf, y = Inf,  hjust = -0.1, vjust = 1.1, size = 2.8, label = top),
  annotate("text", x = Inf,  y = -Inf, hjust = 1.1,  vjust = -1,  size = 2.8, label = bottom))
theme_fig <- theme_bw(base_size = 10) +
  theme(plot.caption = element_text(size = 7.5, colour = "grey30", hjust = 0), plot.caption.position = "plot")
# one row of the results table from a coefficient (per 10 SHMI)
coef_row <- function(figure, estimate, b, se, df, p, n)
  tibble(figure, estimate, effect = 10 * b, se = 10 * se, df,
         p_val = p, n) %>%
  mutate(lower = effect - qt(0.975, df) * se, upper = effect + qt(0.975, df) * se, .after = df)
results <- list()
# PLOT DATA (11 Nov 2026): everything the manuscript figures draw is collected in PD and
# saved as <OUT>/plotdata_yield_main.rds at the end of section 7. The figure scripts in
# R/figures/ draw from that file only, so a figure can be changed without rerunning this
# script. The figures this script still draws itself are kept for comparison.
PD <- list()


################################################################################
# 6. FIG 1-3: WITHIN SITES (research sites)
################################################################################
# The full within-site model:
#   yield ~ site x year + site N terms (linear; + quadratic at N_QUAD_SITES) + other
#           + SHMI_w + SHMI_w:N_c + SHMI_w:age_c
#           (+ NEMEIRR tillage x residue if MEI_SLOPE_BY_RESIDUE) + (1 | trt/plot), AR(1)
# N_c = (N - 150)/100 at N_INTERACT_SITES (0 elsewhere), age_c = system age
# minus its average. SHMI_w's coefficient is the slope at 150 kg N/ha and average age.
f_full <- reformulate(c("site_year", n_terms(corn), o_term(corn), "SHMI_w", mei_terms(corn),
                        "SHMI_w:N_c", "SHMI_w:age_c"), "yield_kgha")
m_full <- lme(f_full, random = ~ 1 | trt/plot_id, correlation = corCAR1(form = ~ year | trt/plot_id),              
              data = corn, method = "REML", control = ctrl)
tt_full <- summary(m_full)$tTable
print(tt_full[c("SHMI_w", "SHMI_w:N_c", "SHMI_w:age_c"), ])
b <- fixef(m_full); V <- vcov(m_full)
TERMS <- c("SHMI_w", "SHMI_w:N_c", "SHMI_w:age_c")
# slope (per 1 SHMI) and its SE at N kg/ha and age years
slope_at <- function(N = N_REF, age = AGE_MEAN) {
  L <- c(1, N / 100 - N_REF / 100, age - AGE_MEAN)
  list(b = sum(L * b[TERMS]), se = sqrt(drop(t(L) %*% V[TERMS, TERMS] %*% L)))
}
DF_W <- tt_full["SHMI_w", "DF"]

# y level: the overall mean yield of the research plots (site x year removes each
# site's level, so the figure places the points and lines on this common baseline;
# only the slopes are estimated). Until 30 Oct 2026 the baseline was the plots' average
# yield minus their county's. Lines pivot at the average SHMI of the data.
base_y <- mean(corn$yield_kgha)
SHMI_PIV <- mean(corn$SHMI)

# ---- Fig 1: the within-site slope -----------------------------------------------------
# points: each treatment's yields with the model's other terms removed (residual +
# fitted SHMI line: a partial-residual plot), at 150 kg N/ha and average age.
# FIG3A_POINTS = "system" (default): one point per management system at each site, its
#   mean SHMI and mean adjusted yield over all its plot-years. Treatments that differ only
#   in factors SHMI does not score (N rate, irrigation, residue removal) have the same SHMI
#   record and are merged into one point. The line is still fitted to every plot-year.
# FIG3A_POINTS = "level": plot-years grouped by SHMI rounded to a whole point (before 9 Nov
#   2026); one point can then mix systems, and one system can span several points.
if (!exists("FIG3A_POINTS")) FIG3A_POINTS <- "system"
s0 <- slope_at()
corn$y_adj <- as.numeric(resid(m_full, level = 0)) + s0$b * (corn$SHMI - SHMI_PIV) + base_y
if (FIG3A_POINTS == "system") {
  pts2 <- corn %>% group_by(site, MGT_trt_base) %>% mutate(SHMI_level = round(mean(SHMI), 1)) %>%
    group_by(site, SHMI_level) %>%
    summarise(SHMI = mean(SHMI), y = mean(y_adj), se = sd(y_adj) / sqrt(n()), n_trt = n_distinct(MGT_trt_base),
              treatments = paste(sort(unique(MGT_trt_base)), collapse = " "),
              years = n_distinct(year), plot_years = n(), .groups = "drop")
} else {
  pts2 <- corn %>% group_by(site, SHMI_level = round(SHMI)) %>%
    summarise(SHMI = mean(SHMI), y = mean(y_adj), se = sd(y_adj) / sqrt(n()), n_trt = n_distinct(MGT_trt_base),
              treatments = paste(sort(unique(MGT_trt_base)), collapse = " "),
              years = n_distinct(year), plot_years = n(), .groups = "drop")
}
write.csv(pts2, file.path(OUT, "Fig1_points.csv"), row.names = FALSE)
line_at <- function(sl, label = NULL) tibble(SHMI = seq(min(corn$SHMI), max(corn$SHMI), length.out = 200)) %>%
  mutate(fit = base_y + sl$b * (SHMI - SHMI_PIV), se = sl$se * abs(SHMI - SHMI_PIV),
         lower = fit - qt(0.975, DF_W) * se, upper = fit + qt(0.975, DF_W) * se, line = label)
ln2 <- line_at(s0)
fig2 <- ggplot() +
  geom_ribbon(data = ln2, aes(SHMI, ymin = lower, ymax = upper), fill = "darkgreen", alpha = 0.15) +
  geom_line(data = ln2, aes(SHMI, fit), colour = "darkgreen", linewidth = 1.1) +
  geom_errorbar(data = pts2, aes(SHMI, ymin = y - se, ymax = y + se), width = 0, colour = "grey35") +
  geom_point(data = pts2, aes(SHMI, y, fill = site), shape = 21, size = 2.6, colour = "black") +
  scale_fill_manual(values = site_cols, name = NULL) +
  corner_labels(lab_slope(10 * s0$b, 10 * (s0$b - qt(0.975, DF_W) * s0$se), 10 * (s0$b + qt(0.975, DF_W) * s0$se)),
                lab_p(tt_full["SHMI_w", "p-value"])) +
  labs(x = XLAB, y = YLAB, title = sprintf("Within research sites (%d sites)", n_distinct(corn$site)),
       caption = paste0("Treatments compared within each site and year. Model (annual plot yields): yield ~ site \u00D7 year\n",
                        "(which removes weather, county yield, soil and the site's level) + N per site + other treatment factors\n",
                        "+ SHMI + SHMI \u00D7 N + SHMI \u00D7 system age + (1 | treatment/plot), AR(1). Slope at ", N_REF, " kg N/ha and the average\n",
                        sprintf("system age (%.0f years). Each point: %s (\u00B1 SE). Each site's level is\n", AGE_MEAN,
                                if (FIG3A_POINTS == "system") "one management system at a site, mean over its plot-years"
                                else "plot-years with the same SHMI (whole point) at a site"),
                        "set by site \u00D7 year, so its points sit along the line; their spread around it is the within-site evidence.")) +
  theme_fig
ggsave(file.path(OUT, "Fig1_within_sites.png"), fig2 + MS + PANEL, width = 6.8, height = 4.6, dpi = 300)
results$within <- coef_row("Fig 1", sprintf("Within sites, at %d kg N/ha and average age", N_REF),
                           s0$b, s0$se, DF_W, tt_full["SHMI_w", "p-value"], n_distinct(corn$site))

# ---- Fig 2 and 3: the slope against N rate and system age ---------------------------
# y = the within-site SHMI slope (yield change per 10 SHMI) with its 95% band, x = N
# applied or system age. Plotting the slope itself (not lines of yield vs SHMI at a few
# N rates) avoids implying yield levels the model does not estimate: site x year
# removes the levels, so lines of yield vs SHMI at different N rates would all pivot
# at one point and suggest, wrongly, that unfertilized corn can out-yield fertilized
# corn at high SHMI. Points mark AGES / N_LINES; the rug shows where the data are.
plot_slope <- function(xs, at, marks, xlab, p_term, p_lab, title, cap, col, rug, label_right = FALSE) {
  sl <- map_dfr(xs, function(x) { s <- at(x); tibble(x = x, b = 10 * s$b, se = 10 * s$se) }) %>%
    mutate(lower = b - qt(0.975, DF_W) * se, upper = b + qt(0.975, DF_W) * se)
  mk <- map_dfr(marks, function(x) tibble(x = x, b = 10 * at(x)$b)) %>% mutate(lab = sprintf("%+.0f", b))
  ggplot(sl, aes(x, b)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = col, alpha = 0.18) +
    geom_line(colour = col, linewidth = 1.2) +
    geom_point(data = mk, size = 2.6, colour = col) +
    geom_text(data = mk, aes(label = lab), vjust = -1.1, size = 3) +
    geom_rug(data = tibble(x = rug), aes(x = x), inherit.aes = FALSE, sides = "b", alpha = 0.25, colour = "grey30") +
    annotate("text", x = if (label_right) Inf else -Inf, y = Inf, hjust = if (label_right) 1.05 else -0.05, vjust = 1.4, size = 2.8,
             label = paste0("Yield change per 10 SHMI, 95% band; ", lab_p(tt_full[p_term, "p-value"], p_lab))) +
    labs(x = xlab, y = "Corn yield change\nper 10 SHMI (kg ha⁻¹)", title = title, caption = cap) + theme_fig
}
n_rug <- corn %>% filter(site %in% n_var_sites) %>% distinct(site, MGT_trt_base, year, N100) %>% pull(N100) * 100
fig3 <- plot_slope(seq(1, max(corn$age), length.out = 100), function(a) slope_at(age = a), AGES,
                   "Time under management (years)", "SHMI_w:age_c", "SHMI × age",
                   "SHMI benefit and system age (research sites)",
                   paste0("Within-site SHMI slope from the Fig 1 model at each system age, at ", N_REF,
                          " kg N/ha; 95% band. Points: ", paste(AGES, collapse = ", "), " years.\n",
                          "The slope, not the yield level, is plotted: site × year removes each site's level. Rug: corn years in the data."),
                   "#238b45", unique(corn$age))
fig4 <- plot_slope(seq(N_FIG_MIN, max(N_LINES), length.out = 100), function(n) slope_at(N = n), N_LINES,
                   "N applied (kg N/ha)\nfertilizer + available manure N + legume credit", "SHMI_w:N_c", "SHMI × N",
                   "SHMI benefit and N rate (research sites)",
                   paste0("Within-site SHMI slope from the Fig 1 model at each N rate, at the average system age; 95% band. Estimated at\n",
                          "the sites with N contrasts (", paste(n_var_sites, collapse = ", "), "). The slope, not the yield level, is plotted:\n",
                          "at any SHMI, fertilized corn still out-yields unfertilized corn. Drawn from ", N_FIG_MIN, " kg N/ha: below that only\n",
                          "NEMLTCRS compared SHMI levels. Rug: N rates in the data."),
                   "#b30000", n_rug[n_rug >= N_FIG_MIN & n_rug <= max(N_LINES)], label_right = TRUE)
ggsave(file.path(OUT, "Fig3_SHMI_x_age.png"), fig3 + MS + PANEL, width = 6.8, height = 4.4, dpi = 300)
ggsave(file.path(OUT, "Fig2_SHMI_x_N.png"),   fig4 + MS + PANEL, width = 6.8, height = 4.4, dpi = 300)
# plot data: Fig. 3a (the within-site slope) and Fig. 4a,b (the slope against N and age)
slope_curve <- function(xs, at) map_dfr(xs, function(x) { s <- at(x); tibble(x = x, b = 10 * s$b, se = 10 * s$se) }) %>%
  mutate(lower = b - qt(0.975, DF_W) * se, upper = b + qt(0.975, DF_W) * se)
slope_marks <- function(marks, at) map_dfr(marks, function(x) tibble(x = x, b = 10 * at(x)$b)) %>% mutate(lab = sprintf("%+.0f", b))
PD$fig3a <- list(points = pts2, line = ln2, points_mode = FIG3A_POINTS, site_cols = site_cols,
                 slope = 10 * s0$b, lower = 10 * (s0$b - qt(0.975, DF_W) * s0$se), upper = 10 * (s0$b + qt(0.975, DF_W) * s0$se),
                 p = tt_full["SHMI_w", "p-value"], n_sites = n_distinct(corn$site), N_REF = N_REF, AGE_MEAN = AGE_MEAN)
PD$fig4ab <- list(
  N   = list(curve = slope_curve(seq(N_FIG_MIN, max(N_LINES), length.out = 100), function(n) slope_at(N = n)),
             marks = slope_marks(N_LINES, function(n) slope_at(N = n)),
             rug = n_rug[n_rug >= N_FIG_MIN & n_rug <= max(N_LINES)], p = tt_full["SHMI_w:N_c", "p-value"]),
  age = list(curve = slope_curve(seq(1, max(corn$age), length.out = 100), function(a) slope_at(age = a)),
             marks = slope_marks(AGES, function(a) slope_at(age = a)),
             rug = unique(corn$age), p = tt_full["SHMI_w:age_c", "p-value"]),
  global_slope = 10 * s0$b, N_REF = N_REF, AGE_MEAN = AGE_MEAN)

# ---- Fig 2b (alternative): lines of yield vs SHMI at the N_LINES rates, each at its own height
# Height at the average SHMI = baseline + the average N response of the N-trial sites
# (each site's own quadratic N curve, relative to N_REF, averaged over the sites with
# equal weight); slope = the within-site SHMI slope at that N rate. The heights describe
# a typical N-trial site, not any one site; the band is the slope's 95% CI.
n_eff <- function(n) mean(vapply(n_curve_sites, function(s_) {
  k1 <- paste0("N1_", s_); k2 <- paste0("N2_", s_)
  b1 <- if (k1 %in% names(b)) b[[k1]] else 0; b2 <- if (k2 %in% names(b)) b[[k2]] else 0
  b1 * (n - N_REF) / 100 + b2 * ((n / 100)^2 - (N_REF / 100)^2)
}, numeric(1)))
sl_n <- lapply(N_LINES, function(n) slope_at(N = n))
lab_n <- sprintf("%s kg N/ha (%+.0f)", N_LINES, vapply(sl_n, function(s) 10 * s$b, numeric(1)))
ln_n <- map_dfr(seq_along(N_LINES), function(k) line_at(sl_n[[k]], lab_n[k]) %>%
                  mutate(across(c(fit, lower, upper), ~ .x + n_eff(N_LINES[k])))) %>%
  mutate(line = factor(line, levels = lab_n))
cols_n <- setNames(c("#fdd49e", "#ef6548","#b30000", "#7f0000"), lab_n)
fig4b <- ggplot(ln_n, aes(SHMI, fit, colour = line, fill = line)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 1.1) +
  scale_colour_manual(values = cols_n, name = "N applied") + scale_fill_manual(values = cols_n, name = "N applied") +
  corner_labels("Yield change per 10 SHMI (kg/ha) in the legend", lab_p(tt_full["SHMI_w:N_c", "p-value"], "SHMI × N")) +
  labs(x = XLAB, y = YLAB, title = "SHMI × N (research sites)",
       caption = paste0("Within-site SHMI slope at each N rate (Fig 1 model, average system age), each line at the height of the\n",
                        "average N response of the N-rate sites (", paste(n_curve_sites, collapse = ", "), "),\n",
                        "relative to ", N_REF, " kg N/ha. Heights describe a typical N-trial site; bands are the 95% CI of the slopes.")) +
  theme_fig
ggsave(file.path(OUT, "Fig2b_SHMI_x_N_lines.png"), fig4b + MS + PANEL, width = 6.8, height = 4.4, dpi = 300)
message(sprintf("Fig 2b: average N response at %s kg N vs %d: %s kg/ha", paste(N_LINES, collapse = " / "), N_REF,
                paste(round(vapply(N_LINES, n_eff, numeric(1))), collapse = " / ")))
results$age <- bind_rows(
  coef_row("Fig 3", "SHMI x age: change in slope per year", tt_full["SHMI_w:age_c", "Value"],
           tt_full["SHMI_w:age_c", "Std.Error"], tt_full["SHMI_w:age_c", "DF"], tt_full["SHMI_w:age_c", "p-value"], n_distinct(corn$site)),
  map_dfr(AGES, function(a) { s <- slope_at(age = a)
    coef_row("Fig 3", sprintf("Slope at %d years", a), s$b, s$se, DF_W, 2 * pt(-abs(s$b / s$se), DF_W), n_distinct(corn$site)) }))
results$N <- bind_rows(
  coef_row("Fig 2", "SHMI x N: change in slope per 100 kg N/ha", tt_full["SHMI_w:N_c", "Value"],
           tt_full["SHMI_w:N_c", "Std.Error"], tt_full["SHMI_w:N_c", "DF"], tt_full["SHMI_w:N_c", "p-value"], length(n_var_sites)),
  map_dfr(sort(unique(c(N_LINES, 100))), function(n) { s <- slope_at(N = n)   # 100: the rate quoted in the talk
    coef_row("Fig 2", sprintf("Slope at %d kg N/ha", n), s$b, s$se, DF_W, 2 * pt(-abs(s$b / s$se), DF_W), length(n_var_sites)) }))

# ---- Effect sizes: what matters more within a place? (from the Fig 1 model) -------------
# Each term's effect on yield (kg/ha) over its typical range, same site and year, all
# else equal, with 95% CI:
#   SHMI    25th to 75th percentile of the research treatments' SHMI, at N_REF and the
#           average system age (the Fig 1 slope times that range)
#   N       25th to 75th percentile of N at the sites with N treatments: each site's N
#           curve over that range, averaged over the sites (at the site's mean SHMI)
#   other treatment factors   each level vs none (residue removal, irrigation, organic)
# Only within-site terms: site x year removes soil, climate and each site's level. System
# age has no row: its main effect is absorbed by site x year (only SHMI x age, Fig 3).
OTHER_LAB <- c(NEMEIRR_res50 = "Residue removal 50% (NEMEIRR)", NEMEIRR_res100 = "Residue removal 100% (NEMEIRR)",
               NECCIRR_Dfct_resN = "Deficit irrigation (NECCIRR)",   # NECCIRR baseline: full irrigation, residue kept
               NECCIRR_Full_resA = "Residue removal 100% (NECCIRR)",
               NECCIRR_Dfct_resA = "Deficit irrigation + residue removal 100% (NECCIRR)",
               WICST_organic = "Organic system (WICST)",
               MNMOFS_organic = "Organic system (MNMOFS)",
               IAAM7071_resPartial = "Residue removal, partial (IAAM7071)", IAAM7071_resYes = "Residue removal (IAAM7071)",
               PAUP_resPartial = "Residue removal, partial (PAUP)", PAUP_resYes = "Residue removal (PAUP)",
               INWLACRE_splitN = "Split N application (INWLACRE)", INWLACRE_planter = "No-till planter with row cleaners (INWLACRE)",
               INWLTPAC_splitN = "Split N application (INWLTPAC)", INWLTPAC_planter = "No-till planter with row cleaners (INWLTPAC)")
local({
  q25_75 <- function(x) unname(quantile(x, c(0.25, 0.75)))
  eff <- function(group, term, w) {                # w: named weights on fixed effects
    L <- setNames(numeric(length(b)), names(b)); L[names(w)] <- unlist(w)
    est <- sum(L * b); se <- sqrt(drop(t(L) %*% V %*% L))
    tibble(group, term, effect = est, se, df = DF_W, lower = est - qt(0.975, DF_W) * se,
           upper = est + qt(0.975, DF_W) * se, p_val = 2 * pt(-abs(est / se), DF_W))
  }
  qs <- q25_75(corn %>% group_by(site, MGT_trt_base) %>% summarise(s = mean(SHMI), .groups = "drop") %>% pull(s))
  rows <- list(eff("Management", sprintf("SHMI %.0f → %.0f", qs[1], qs[2]), list(SHMI_w = diff(qs))))
  n_cols <- intersect(paste0("N1_", n_curve_sites), names(b))   # N-rate sites only: the linear-only
  if (length(n_cols)) {                                         # sites' slopes do not extrapolate
    sites_n <- sub("^N1_", "", n_cols)
    qn <- q25_75(corn$N100[corn$site %in% sites_n])
    w <- list()
    for (s_ in sites_n) {
      w[[paste0("N1_", s_)]] <- (qn[2] - qn[1]) / length(sites_n)
      if (paste0("N2_", s_) %in% names(b)) w[[paste0("N2_", s_)]] <- (qn[2]^2 - qn[1]^2) / length(sites_n)
    }
    rows[[length(rows) + 1]] <- eff("Management", sprintf("N %.0f → %.0f kg/ha (%d sites)", 100 * qn[1], 100 * qn[2],
                                                          length(sites_n)), w)
  }
  for (k in grep("^other", names(b), value = TRUE)) {
    lv <- sub("^other", "", k)
    rows[[length(rows) + 1]] <- eff("Other treatments (vs none)",
                                    ifelse(lv %in% names(OTHER_LAB), OTHER_LAB[[lv]], gsub("_", " ", lv)), setNames(list(1), k))
  }
  es <- bind_rows(rows) %>%
    mutate(pct_yield = 100 * effect / mean(corn$yield_kgha),
           group = factor(group, levels = c("Management", "Other treatments (vs none)")))
  write.csv(es %>% mutate(across(where(is.numeric), ~ signif(.x, 3))),
            file.path(OUT, "Fig_effect_sizes.csv"), row.names = FALSE)
  print(es, n = Inf, width = Inf)
  results$effects <<- es %>% transmute(figure = "Effect sizes", estimate = paste0(group, ": ", term),
                                       effect, se, df, lower, upper, p_val, n = n_distinct(corn$site))
  es <- es %>% arrange(group, desc(effect)) %>% mutate(term = factor(term, levels = rev(unique(term))))
  p <- ggplot(es, aes(effect, term)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.25, colour = "grey40", orientation = "y") +
    geom_point(aes(fill = group), shape = 21, size = 3, colour = "black") +
    facet_grid(group ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_fill_manual(values = c("darkgreen", "grey70"), guide = "none") +
    labs(x = "Change in corn yield (kg ha⁻¹)", y = NULL,
         title = "What matters more within a site? Each term over its typical range",
         caption = paste0("From the Fig 1 model (research sites): the change in yield from the 25th to the 75th percentile of each term,\n",
                          "same site and year, all else equal, with 95% CI. SHMI at ", N_REF, " kg N/ha and average system age; N averaged over\n",
                          "the sites' N curves; other treatments vs none at their site. Soil and climate drop out with site × year.")) +
    theme_bw(base_size = 11) +
    theme(plot.title.position = "plot",
          strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1),
          strip.background = element_blank(), panel.grid.minor = element_blank(),
          plot.caption = element_text(size = 7.5, colour = "grey30", hjust = 0), plot.caption.position = "plot")
  ggsave(file.path(OUT, "Fig_effect_sizes.png"), p + MS, width = 7.5, height = 1.9 + 0.32 * nrow(es), dpi = 300)
})


################################################################################
# 6b. ROBUSTNESS: does SHMI capture the treatment effects? Leave-one-out refits
################################################################################
# (a) SHMI vs the treatments themselves (research sites). SHMI is computed from the
#     management records, so at a research site it is a one-number summary of each
#     treatment. If that one number captures most of what the treatments do to yield,
#     the index is doing its job. Four models, all fitted the same way (lme4, ML,
#     random plot, no AR(1), so their AIC are comparable):
#       no SHMI: site x year + N curves + other treatment factors
#       + SHMI (one slope)
#       + SHMI, SHMI x N, SHMI x age (the Fig 1 model)
#       treatment factor (every treatment at every site) in place of SHMI
#     "Share of treatment fit" = how much of the improvement in fit that the full
#     treatment factor achieves over the no-SHMI model is achieved by SHMI.
#     Note: some treatment contrasts are aliased with the N curves and "other" factors;
#     lme4 drops those columns, and fixed_df counts what is left.
corn$trt <- factor(paste(corn$site, corn$MGT_trt_base))
f_base <- c("site_year", n_terms(corn), o_term(corn))
fits_trt <- list(
  "No SHMI (site x year, N curves, other factors)"  = f_base,
  "+ SHMI, one slope"                                 = c(f_base, "SHMI_w", mei_terms(corn)),
  "+ SHMI, SHMI x N, SHMI x age (Fig 1 model)"        = c(f_base, "SHMI_w", mei_terms(corn), "SHMI_w:N_c", "SHMI_w:age_c"),
  "Treatment factor in place of SHMI"                 = c(f_base, "trt"))
trt_tab <- imap_dfr(fits_trt, function(rhs, nm) {
  m <- lme4::lmer(reformulate(c(rhs, "(1 | plot_id)"), "yield_kgha"), data = corn, REML = FALSE)
  tibble(model = nm, fixed_df = length(lme4::fixef(m)), logLik = as.numeric(logLik(m)), AIC = AIC(m), BIC = BIC(m))
})
dev0 <- -2 * trt_tab$logLik[1]; devT <- -2 * trt_tab$logLik[4]
trt_tab <- trt_tab %>%
  mutate(extra_df = fixed_df - fixed_df[1], dAIC = AIC - min(AIC),
         share_of_treatment_fit = round(100 * (dev0 + 2 * logLik) / (dev0 - devT), 1),
         n_treatments = n_distinct(corn$trt), n_sites = n_distinct(corn$site))
print(trt_tab, width = Inf)
write.csv(trt_tab, file.path(OUT, "SHMI_vs_treatment.csv"), row.names = FALSE)

# (b) leave one research site out: the Fig 1 slope (150 kg N/ha, average age) without each site
loo_site <- map_dfr(sort(unique(as.character(corn$site))), function(s_) {
  dd <- droplevels(filter(corn, site != s_))
  f <- reformulate(c("site_year", n_terms(dd), o_term(dd), "SHMI_w", mei_terms(dd),
                     if (any(dd$N_c != 0)) "SHMI_w:N_c", "SHMI_w:age_c"), "yield_kgha")
  m <- tryCatch(lme(f, random = ~ 1 | trt/plot_id, correlation = corCAR1(form = ~ year | trt/plot_id),
                    data = dd, method = "REML", control = ctrl), error = function(e) NULL)
  if (is.null(m)) { message("Leave-one-out without ", s_, ": model did not fit"); return(NULL) }
  tt <- summary(m)$tTable["SHMI_w", ]
  out <- coef_row("Fig 1, leave one site out", paste("Without", s_), tt[["Value"]], tt[["Std.Error"]], tt[["DF"]],
                  tt[["p-value"]], n_distinct(dd$site))
  # SHMI x N without each of the sites with N contrasts: does the interaction rest on one site?
  tN <- summary(m)$tTable
  if (s_ %in% n_var_sites && "SHMI_w:N_c" %in% rownames(tN))
    out <- bind_rows(out, coef_row("Fig 2, leave one N site out",
                                   paste("SHMI x N per 100 kg N/ha, without", s_),
                                   tN["SHMI_w:N_c", "Value"], tN["SHMI_w:N_c", "Std.Error"],
                                   tN["SHMI_w:N_c", "DF"], tN["SHMI_w:N_c", "p-value"],
                                   length(setdiff(n_var_sites, s_))))
  out
})
results$loo <- loo_site
print(results$loo, n = Inf, width = Inf)


################################################################################
# 7. FIG 4: FOREST PLOT
################################################################################
# Rows (yield change per 10 SHMI, 95% CI):
#   Within sites (all sites)       the Fig 1 slope
#   Pooled across sites            random-effects mean of the site slopes (metafor,
#                                  REML): each site counts once, whatever its size
#   each site                      its own within-site model (as fit_site); beneath
#                                  it, open circles: the N-contrast sites at 50 and 150
#                                  kg N/ha (N_SUB),
#                                  NEMEIRR at 0/50/100% residue removal and NECCIRR at
#                                  full vs deficit irrigation (each from the site's own
#                                  model with that interaction)
#   Minnesota farms                Fig 5 (between farms)
per10 <- function(m, term = "SHMI_w") {
  tt <- summary(m)$tTable[term, ]
  tibble(effect = 10 * tt[["Value"]], se = 10 * tt[["Std.Error"]], df = tt[["DF"]], p_val = tt[["p-value"]])
}
site_models <- split(corn, corn$site) %>% map(fit_site)
site_rows <- imap_dfr(site_models, ~ per10(.x) %>% mutate(label = .y, kind = "site", n = .x$dims$N))
pool <- rma(yi = effect, sei = se, data = site_rows, method = "REML", test = "knha")

# ---- heterogeneity among the site slopes, and how much N and system age explain -------
# Each site's slope (its own model, no interactions) is roughly its slope at its own mean
# N and mean system age. Moving every site to N_REF and the average age with the Fig 1
# model's SHMI x N and SHMI x age coefficients, then pooling again, shows how much of the
# spread among sites those two account for (the adjustment's own uncertainty is ignored).
# A meta-regression of the site slopes on each site's mean N asks the same for N alone.
site_mod <- corn %>% group_by(site) %>%
  summarise(N_mean = 100 * mean(N100), age_mean = mean(age), .groups = "drop")
b_full <- fixef(m_full)
site_rows2 <- site_rows %>% left_join(site_mod, by = c("label" = "site"))
site_adj <- site_rows2 %>%
  mutate(effect = effect + 10 * b_full[["SHMI_w:N_c"]] * (N_REF - N_mean) / 100 +
                           10 * b_full[["SHMI_w:age_c"]] * (AGE_MEAN - age_mean))
pool_adj <- rma(yi = effect, sei = se, data = site_adj, method = "REML", test = "knha")
mreg_N <- rma(yi = effect, sei = se, mods = ~ I((N_mean - N_REF) / 100), data = site_rows2, method = "REML", test = "knha")
het_rows <- function(m, label) tibble(
  figure = "Fig 4, heterogeneity",
  estimate = paste0(label, c(": Q (df = sites - 1)", ": I2 (%)", ": tau (between-site SD of the slope, kg/ha per 10 SHMI)")),
  effect = c(m$QE, m$I2, sqrt(m$tau2)), se = NA_real_, df = c(m$k - m$p, NA, NA),
  lower = NA_real_, upper = NA_real_, p_val = c(m$QEp, NA, NA), n = m$k)
results$het <- bind_rows(
  het_rows(pool, "Site slopes as estimated"),
  tibble(figure = "Fig 4, heterogeneity", estimate = sprintf("Pooled across sites, each site at %d kg N/ha and age %.0f y", N_REF, AGE_MEAN),
         effect = as.numeric(pool_adj$b), se = pool_adj$se, df = pool_adj$k - 1, lower = pool_adj$ci.lb, upper = pool_adj$ci.ub,
         p_val = pool_adj$pval, n = pool_adj$k),
  het_rows(pool_adj, sprintf("Site slopes at %d kg N/ha and average age", N_REF)),
  tibble(figure = "Fig 4, heterogeneity", estimate = "Meta-regression: change in site slope per 100 kg/ha higher site mean N",
         effect = mreg_N$beta[2], se = mreg_N$se[2], df = mreg_N$k - mreg_N$p, lower = mreg_N$ci.lb[2], upper = mreg_N$ci.ub[2],
         p_val = mreg_N$pval[2], n = mreg_N$k),
  tibble(figure = "Fig 4, heterogeneity", estimate = "Meta-regression on site mean N: share of between-site variance explained (R2, %)",
         effect = mreg_N$R2, se = NA_real_, df = NA_real_, lower = NA_real_, upper = NA_real_, p_val = mreg_N$QMp, n = mreg_N$k))
print(pool); print(pool_adj); print(mreg_N)
# sub-rows: NEMEIRR by residue removal
# (always shown: NEMEIRR's own model with SHMI x residue removal, whatever
# MEI_SLOPE_BY_RESIDUE says about the global model; its main row stays the one slope)
sub_mei <- if ("NEMEIRR" %in% names(site_models)) {
  # tryCatch: with other scores (coarse records, SHMI v1) the SHMI x residue model can be singular
  m <- tryCatch(fit_site(filter(corn, site == "NEMEIRR"),
                         extra = if (!MEI_SLOPE_BY_RESIDUE) c("SHMI_w:mei_res50", "SHMI_w:mei_res100")),
                error = function(e) { message("NEMEIRR SHMI x residue model not estimable: ", conditionMessage(e)); NULL })
}
sub_mei <- if (!is.null(sub_mei)) {
  m <- sub_mei
  bm <- fixef(m); Vm <- vcov(m)
  ints <- intersect(c("SHMI_w:mei_res50", "SHMI_w:mei_res100"), names(bm))
  corn %>% filter(site == "NEMEIRR") %>% distinct(residue, mei_res50, mei_res100) %>% arrange(as.numeric(residue)) %>%
    pmap_dfr(function(residue, mei_res50, mei_res100) {
      L <- setNames(numeric(length(bm)), names(bm)); L["SHMI_w"] <- 1
      L[ints] <- c(`SHMI_w:mei_res50` = mei_res50, `SHMI_w:mei_res100` = mei_res100)[ints]
      est <- 10 * sum(L * bm); se <- 10 * sqrt(drop(t(L) %*% Vm %*% L))
      tibble(parent = "NEMEIRR", label = paste0("NEMEIRR, no-till vs disk, ", residue, "% residue removed"), effect = est, se = se,
             df = Inf, p_val = 2 * pnorm(-abs(est / se)))
    })
} else NULL
# sub-rows: the slope at N_SUB (50 and 150 kg N/ha) at each N-contrast site. Not 0 N:
# NECCIRR's 0 N plots all have the same SHMI and NELITCSE has no 0 N plots, so a
# slope at 0 N would be an extrapolation at most sites
sub_n <- map_dfr(n_var_sites, function(s_) {
  m <- tryCatch(fit_site(filter(corn, site == s_), extra = "SHMI_w:N_c"), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  bm <- fixef(m); Vm <- vcov(m)
  map_dfr(N_SUB, function(n) {
    a <- n / 100 - N_REF / 100
    est <- 10 * (bm[["SHMI_w"]] + a * bm[["SHMI_w:N_c"]])
    se <- 10 * sqrt(Vm["SHMI_w", "SHMI_w"] + a^2 * Vm["SHMI_w:N_c", "SHMI_w:N_c"] + 2 * a * Vm["SHMI_w", "SHMI_w:N_c"])
    tibble(parent = s_, label = paste0(s_, ", ", n, " kg N/ha"), effect = est, se = se, df = Inf,
           p_val = 2 * pnorm(-abs(est / se)))
  })
})
# SHMI x N in each experiment's own model, wherever N varies within the experiment (Methods:
# "SHMI x N in each experiment's own model"); per 10 SHMI points per 100 kg N/ha.
nx_tab <- map_dfr(names(site_models), function(s_) {
  d <- filter(corn, site == s_)
  if (!is.finite(sd(d$N_c, na.rm = TRUE)) || sd(d$N_c, na.rm = TRUE) * 100 < 1) return(NULL)
  m <- tryCatch(fit_site(d, extra = "SHMI_w:N_c"), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  tt <- summary(m)$tTable
  if (!"SHMI_w:N_c" %in% rownames(tt)) return(NULL)
  r <- tt["SHMI_w:N_c", ]
  tibble(site = s_, N_range = paste(round(range(d$N_c, na.rm = TRUE) * 100 + N_REF), collapse = "-"),
         SHMIxN_per10_per100N = 10 * r[["Value"]], se = 10 * r[["Std.Error"]], df = r[["DF"]], p = r[["p-value"]]) %>%
    mutate(lower = SHMIxN_per10_per100N - qt(0.975, df) * se, upper = SHMIxN_per10_per100N + qt(0.975, df) * se)
})
if (nrow(nx_tab)) {
  write.csv(mutate(nx_tab, across(where(is.double), ~ signif(.x, 4))), file.path(OUT, "SHMI_x_N_by_site.csv"), row.names = FALSE)
  cat("\nSHMI x N in each experiment's own model (per 10 SHMI per 100 kg N/ha):\n"); print(as.data.frame(nx_tab), digits = 3)
}
# sub-rows: NECCIRR by irrigation (full vs deficit), from NECCIRR's own model with
# SHMI x irrigation (centred, so its main row stays the average slope). Exploratory: does
# the SHMI benefit differ when water is short?
sub_irr <- if ("NECCIRR" %in% names(site_models)) {
  dd <- filter(corn, site == "NECCIRR") %>% mutate(dfct = as.numeric(irrigation == "Dfct"))
  m_dfct <- mean(distinct(dd, MGT_trt_base, dfct)$dfct)
  dd$dfct_c <- dd$dfct - m_dfct
  m <- tryCatch(fit_site(dd, extra = "SHMI_w:dfct_c"), error = function(e) NULL)
  if (is.null(m)) NULL else {
    bm <- fixef(m); Vm <- vcov(m); k <- "SHMI_w:dfct_c"
    results$irr <- coef_row("Fig 4, NECCIRR", "SHMI x deficit irrigation: deficit minus full slope",
                            bm[[k]], sqrt(Vm[k, k]), summary(m)$tTable[k, "DF"], summary(m)$tTable[k, "p-value"], 1)
    map_dfr(c(full = -m_dfct, deficit = 1 - m_dfct), function(a) {
      est <- 10 * (bm[["SHMI_w"]] + a * bm[[k]])
      se <- 10 * sqrt(Vm["SHMI_w", "SHMI_w"] + a^2 * Vm[k, k] + 2 * a * Vm["SHMI_w", k])
      tibble(effect = est, se = se, df = Inf, p_val = 2 * pnorm(-abs(est / se)))
    }, .id = "irr") %>%
      transmute(parent = "NECCIRR", label = paste0("NECCIRR, ", irr, " irrigation"), effect, se, df, p_val)
  }
} else NULL
# ---- SHMI x other treatment factors, each experiment (7 Nov 2026) -----------------------
# For every experiment with a treatment factor SHMI does not score ("other": residue
# removal, irrigation, organic system, planter), its own model (as fit_site) with the SHMI
# slope free to differ by level: one SHMI_w x level term for each non-reference level in
# which SHMI varies. Levels in which SHMI does not vary have no slope of their own (noted).
# Reported per level: slope per 10 SHMI at the experiment's own N and age, and moved to
# N_REF and the average age by the experiment's shift in site_SHMI_ranges (slope_adj -
# slope); the SHMI x level interaction F test (the slope is the same in every level):
# F = Wald chi2 / numerator df, denominator df from lme for the interaction terms. Exploratory; the pooled
# slope in Fig. 3c assumes one slope per experiment. Written: SHMI_x_other_by_site.csv
# (Supplementary Data Table 5); Fig. 3c flags these experiments with a dagger.
OTHER_BASE <- c(NEMEIRR = "0% stover removed", NECCIRR = "full irrigation, residue kept",
                WICST_ARS = "conventional", IAAM7071 = "conventional, residue kept",
                PAUP = "conventional, residue kept", INWLACRE = "standard (conventional, residue kept)")
other_label <- function(x, s_) {
  if (x == "none") return(if (!is.na(OTHER_BASE[s_])) OTHER_BASE[[s_]] else "reference")
  x <- sub(paste0("^", s_, "_"), "", x)
  x <- sub("^res(\\d+)$", "\\1% stover removed", x)
  x <- gsub("Dfct", "deficit irrigation", gsub("Full", "full irrigation", x))
  x <- gsub("_resN", ", residue kept", gsub("_resA", ", residue removed", x))
  x <- gsub("resPartial", "residue partly removed", gsub("resYes", "residue removed", x))
  gsub("_", ", ", x)
}
shift_of <- setNames(site_adj$effect - site_rows$effect, site_rows$label)
other_tab <- imap_dfr(site_models, function(m0, s_) {
  d <- droplevels(filter(corn, site == s_))
  if (nlevels(d$other) < 2) return(NULL)
  lev <- levels(d$other)
  info <- map_dfr(lev, function(l) {
    x <- d %>% filter(other == l) %>% distinct(MGT_trt_base, .keep_all = TRUE)
    tibble(level = l, n_trt = nrow(x), plot_years = sum(d$other == l),
           SHMI_range = sprintf("%.0f-%.0f", min(x$SHMI), max(x$SHMI)),
           vary = n_distinct(round(x$SHMI, 1)) > 1)
  })
  est <- info$level[info$vary]
  out <- info %>% mutate(site = s_, factor_level = vapply(level, other_label, "", s_ = s_),
                         slope = NA_real_, se = NA_real_, lower = NA_real_, upper = NA_real_,
                         slope_adj = NA_real_, F_int = NA_real_, df1 = NA_integer_, df2 = NA_real_, p_int = NA_real_,
                         note = ifelse(vary, "", "SHMI does not vary within this level: no slope of its own"))
  if (length(est) < 2) { out$note[out$vary] <- "only one level with an SHMI contrast: no test"; return(out) }
  vn <- paste0("oth_", make.names(est[-1]))
  for (k in seq_along(vn)) d[[vn[k]]] <- as.numeric(d$other == est[k + 1])
  ex <- paste0("SHMI_w:", vn)
  m <- tryCatch(fit_site(d, extra = ex), error = function(e) { message("  ", s_, " SHMI x other: ", conditionMessage(e)); NULL })
  if (is.null(m)) { out$note[out$vary] <- "SHMI x level model not estimable"; return(out) }
  bm <- fixef(m); Vm <- vcov(m); ex <- intersect(ex, names(bm))
  dfm <- summary(m)$tTable["SHMI_w", "DF"]
  W <- tryCatch(drop(t(bm[ex]) %*% solve(Vm[ex, ex, drop = FALSE]) %*% bm[ex]), error = function(e) NA_real_)
  for (i in seq_along(est)) {
    L <- setNames(numeric(length(bm)), names(bm)); L["SHMI_w"] <- 1
    if (i > 1 && paste0("SHMI_w:", vn[i - 1]) %in% ex) L[paste0("SHMI_w:", vn[i - 1])] <- 1
    b <- 10 * sum(L * bm); se <- 10 * sqrt(drop(t(L) %*% Vm %*% L)); j <- which(out$level == est[i])
    out$slope[j] <- b; out$se[j] <- se
    out$lower[j] <- b - qt(0.975, dfm) * se; out$upper[j] <- b + qt(0.975, dfm) * se
    out$slope_adj[j] <- b + shift_of[[s_]]
  }
  df2 <- min(summary(m)$tTable[ex, "DF"])
  out$F_int <- W / length(ex); out$df1 <- length(ex); out$df2 <- df2
  out$p_int <- pf(W / length(ex), length(ex), df2, lower.tail = FALSE)
  out
}) %>% { if (nrow(.)) select(., site, factor_level, n_trt, plot_years, SHMI_range, slope, lower, upper, slope_adj,
                                 F_int, df1, df2, p_int, note) else . }
if (nrow(other_tab)) {
  write.csv(mutate(other_tab, across(where(is.double), ~ signif(.x, 4))), file.path(OUT, "SHMI_x_other_by_site.csv"), row.names = FALSE)
  cat("\nSHMI slope by level of the treatment factors SHMI does not score (per 10 SHMI):\n"); print(as.data.frame(other_tab), digits = 3)
}

subs <- bind_rows(tibble(parent = character(), label = character(), effect = numeric(), se = numeric(),
                         df = numeric(), p_val = numeric()), sub_mei, sub_n, sub_irr)
forest <- bind_rows(
  tibble(label = "Within sites (all sites)", kind = "within", effect = results$within$effect,
         lower = results$within$lower, upper = results$within$upper, p_val = results$within$p_val),
  tibble(label = "Pooled across sites", kind = "pooled", effect = as.numeric(pool$b),
         lower = pool$ci.lb, upper = pool$ci.ub, p_val = pool$pval),
  tibble(label = sprintf("Pooled, at %d kg N/ha and average age", N_REF), kind = "pooled", effect = as.numeric(pool_adj$b),
         lower = pool_adj$ci.lb, upper = pool_adj$ci.ub, p_val = pool_adj$pval),
  map_dfr(names(site_models), function(s_) bind_rows(
    site_rows %>% filter(label == s_) %>% mutate(lower = effect - qt(0.975, df) * se, upper = effect + qt(0.975, df) * se),
    subs %>% filter(parent == s_) %>% mutate(kind = "sub", lower = effect - 1.96 * se, upper = effect + 1.96 * se)) %>%
      mutate(site = s_))) %>%
  mutate(sig = case_when(p_val < 0.001 ~ "***", p_val < 0.01 ~ "**", p_val < 0.05 ~ "*", p_val < 0.1 ~ "+", TRUE ~ ""),
         row = rev(seq_len(n())), site = coalesce(site, ""))
# display labels (the CSV keeps the full labels)
SITE_LOC <- c(NECCIRR = "Clay Center, NE", NELITCSE = "Lincoln, NE", NEMLTCRS = "Mead, NE", NEMEIRR = "Mead, NE",
              OHTVD_Ho = "Hoytville, OH", OHTVD_Wo = "Wooster, OH", WICST_ARS = "Arlington, WI",
              MNMOFS = "Morris, MN", INWLACRE = "West Lafayette, IN (ACRE)", INWLTPAC = "West Lafayette, IN (TPAC)",
              IAAM7071 = "Ames, IA", PAUP = "University Park, PA")
forest <- forest %>% mutate(display = case_when(
  kind == "within" ~ "Global slope (one model, all sites)",
  kind == "pooled" & grepl("^Pooled across", label) ~ "Average of site slopes",
  kind == "pooled" ~ sprintf("Average of site slopes at %d kg N, %.0f y", N_REF, AGE_MEAN),
  kind == "site" ~ paste0(label, ifelse(label %in% names(SITE_LOC), paste0(" (", SITE_LOC[label], ")"), "")),
  kind == "sub" & grepl("kg N/ha$", label) ~ paste0("at ", sub("^[^,]+, ", "", label)),
  kind == "sub" ~ sub("^[^,]+, ", "", label),
  kind == "farm" ~ "Discovery Farms MN (between farms)",
  TRUE ~ label))
write.csv(select(forest, -row), file.path(OUT, "Fig4_forest_rows.csv"), row.names = FALSE)
rng <- diff(range(c(forest$lower, forest$upper)))
sty <- tribble(~kind,     ~shape, ~fill,        ~size,
               "within",  23,     "darkgreen",  4.5,
               "pooled",  23,     "firebrick",  4.5,
               "site",    21,     "steelblue",  3.2,
               "sub",     21,     "white",      2.2,
               "farm",    24,     FARM_COL,     3.6)
forest <- left_join(forest, sty, by = "kind")
fig5 <- ggplot(forest, aes(effect, row)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_hline(yintercept = max(forest$row) - sum(forest$kind %in% c("within", "pooled")) + 0.5, colour = "grey85") +   # summaries | sites
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.25, colour = "grey40", orientation = "y") +
  geom_point(aes(shape = shape, fill = fill, size = size), colour = "black", stroke = 0.4) +
  scale_shape_identity() + scale_fill_identity() + scale_size_identity() +
  geom_text(aes(x = upper + 0.02 * rng, label = sig), hjust = 0, size = 3) +
  scale_y_continuous(breaks = forest$row, labels = forest$display, expand = expansion(add = 0.6)) +
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.08))) +
  labs(x = "Corn yield change per 10 SHMI (kg ha\u207B\u00B9)", y = NULL,
       caption = paste0("Global slope: the Fig 1 model (one slope for all sites) at ", N_REF, " kg N/ha and average age. Average of site\n",
                        sprintf("slopes: random-effects mean, each site counted once (site slopes differ: Q = %.1f, df %d, p %s; I\u00B2 = %.0f%%);\n",
                                pool$QE, pool$k - 1, ifelse(pool$QEp < 0.001, "< 0.001", sprintf("= %.3f", pool$QEp)), pool$I2),
                        sprintf("below it, with each site moved to %d kg N/ha and %.0f years by the Fig 2-3 interactions (I\u00B2 = %.0f%%).\n",
                                N_REF, AGE_MEAN, pool_adj$I2),
                        "Site rows: each site's own model, averaged over its N rates and years; beneath, the slope at ", paste(N_SUB, collapse = " and "), "\n",
                        "kg N/ha where N varied, NEMEIRR by residue removal and NECCIRR by irrigation (exploratory).\n",
                        "*** p < 0.001, ** p < 0.01, * p < 0.05, + p < 0.1.")) +
  theme_fig + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())
ggsave(file.path(OUT, "Fig4_forest.png"),
       fig5 + MS + force_panelsizes(rows = unit(0.22 * nrow(forest), "in"), cols = unit(3.4, "in")),
       width = 6.8, height = 1.6 + 0.22 * nrow(forest), dpi = 300)


################################################################################
# 7b. ONE SLOPE ALONG THE SHMI SCALE: THE EXPERIMENTS AS A COLLECTION
################################################################################
# Site x year removes each site's level, so the global slope comes only from treatment
# differences within sites. Each experiment covers its own stretch of the SHMI scale
# (set by the treatments it was designed to compare), and the global model joins the
# stretches with one slope. Checks that one slope is fair along the whole scale, each
# added to the Fig 1 model (SHMI x N and SHMI x age kept, so a site's position on the
# scale is not mistaken for its system age):
#   (a) position: + SHMI_w x site mean SHMI. Does the slope depend on where a site's
#       stretch sits? Reported as the change in the slope (per 10 SHMI) per 10 points
#       higher site mean.
#   (b) curvature within sites: + (SHMI_w / 10)^2. Slope reported 15 points below and
#       above each site's mean.
#   (c) curvature on the absolute scale: + ((SHMI - pivot) / 10)^2 (with site x year this
#       combines (a) and (b)). Slope reported at SHMI 40, 60 and 80.
#   (d) the same question between site slopes: meta-regression of each site's own slope
#       on its mean SHMI (as the N meta-regression in section 7).
# Also written: site_SHMI_ranges.csv (each experiment's treatments, SHMI range, the part
# of SHMI that drives its range, and its slope: the basis of the Supplementary Table of
# contrasts) and Fig_SHMI_ranges.png (each site's stretch of the scale).
site_means <- corn %>% distinct(site, MGT_combo, SHMI) %>%
  group_by(site) %>% summarise(site_SHMI = mean(SHMI), .groups = "drop")
corn <- corn %>% select(-any_of(c("site_SHMI", "pos_c", "q_w", "q_abs"))) %>%
  left_join(site_means, by = "site") %>%
  mutate(pos_c = (site_SHMI - mean(site_means$site_SHMI)) / 10,
         q_w   = (SHMI_w / 10)^2,
         q_abs = ((SHMI - SHMI_PIV) / 10)^2)
f_scale <- function(extra) reformulate(c("site_year", n_terms(corn), o_term(corn), "SHMI_w", mei_terms(corn),
                                         "SHMI_w:N_c", "SHMI_w:age_c", extra), "yield_kgha")
fit_scale <- function(extra) tryCatch(
  lme(f_scale(extra), random = ~ 1 | trt/plot_id, correlation = corCAR1(form = ~ year | trt/plot_id),
      data = corn, method = "REML", control = ctrl),
  error = function(e) { message("Scale check (", extra, ") did not fit: ", conditionMessage(e)); NULL })
# slope (per 1 SHMI, at 150 kg N/ha and average age) from a model with a curvature term:
# d/dSHMI of b_w * SHMI_w + b_q * ((x - centre) / 10)^2 = b_w + 2 * b_q * (x - centre) / 100
slope_curved <- function(m, qterm, dx) {
  bm <- fixef(m); Vm <- vcov(m); L <- c(1, 2 * dx / 100)
  list(b = sum(L * bm[c("SHMI_w", qterm)]), se = sqrt(drop(t(L) %*% Vm[c("SHMI_w", qterm), c("SHMI_w", qterm)] %*% L)),
       df = summary(m)$tTable["SHMI_w", "DF"])
}
scale_rows <- list()
m_pos <- fit_scale("SHMI_w:pos_c")
if (!is.null(m_pos)) {
  tp <- summary(m_pos)$tTable
  scale_rows$pos <- bind_rows(
    coef_row("Fig 1, SHMI scale check", "(a) Change in slope per 10 points higher site mean SHMI",
             tp["SHMI_w:pos_c", "Value"], tp["SHMI_w:pos_c", "Std.Error"], tp["SHMI_w:pos_c", "DF"],
             tp["SHMI_w:pos_c", "p-value"], n_distinct(corn$site)),
    coef_row("Fig 1, SHMI scale check", "(a) Global slope with the position term (at the average site mean)",
             tp["SHMI_w", "Value"], tp["SHMI_w", "Std.Error"], tp["SHMI_w", "DF"], tp["SHMI_w", "p-value"],
             n_distinct(corn$site)))
}
m_qw <- fit_scale("q_w")
if (!is.null(m_qw)) {
  tq <- summary(m_qw)$tTable
  scale_rows$qw <- bind_rows(
    tibble(figure = "Fig 1, SHMI scale check", estimate = "(b) Curvature within sites: p for (SHMI_w/10)^2",
           effect = NA_real_, se = NA_real_, df = tq["q_w", "DF"], lower = NA_real_, upper = NA_real_,
           p_val = tq["q_w", "p-value"], n = n_distinct(corn$site)),
    map_dfr(c(-15, 15), function(dx) { s <- slope_curved(m_qw, "q_w", dx)
      coef_row("Fig 1, SHMI scale check", sprintf("(b) Slope %d points %s the site mean", abs(dx), ifelse(dx < 0, "below", "above")),
               s$b, s$se, s$df, 2 * pt(-abs(s$b / s$se), s$df), n_distinct(corn$site)) }))
}
m_qa <- fit_scale("q_abs")
if (!is.null(m_qa)) {
  ta <- summary(m_qa)$tTable
  scale_rows$qa <- bind_rows(
    tibble(figure = "Fig 1, SHMI scale check", estimate = "(c) Curvature on the SHMI scale: p for ((SHMI - pivot)/10)^2",
           effect = NA_real_, se = NA_real_, df = ta["q_abs", "DF"], lower = NA_real_, upper = NA_real_,
           p_val = ta["q_abs", "p-value"], n = n_distinct(corn$site)),
    map_dfr(c(40, 60, 80), function(x) { s <- slope_curved(m_qa, "q_abs", x - SHMI_PIV)
      coef_row("Fig 1, SHMI scale check", sprintf("(c) Slope at SHMI %d", x),
               s$b, s$se, s$df, 2 * pt(-abs(s$b / s$se), s$df), n_distinct(corn$site)) }))
}
# (a) and (c) with each site left out: does the trend along the scale rest on one
# experiment? One row per site left out (SHMI_scale_check_loo.csv): the position effect
# (a; per 10 points higher site mean), p for the curvature term (c) and the slope at
# SHMI 40, 60 and 80 from (c). The pivot of (c) stays the all-site average SHMI, and
# the centring of (a) does not change its interaction coefficient.
fit_scale_on <- function(d, extra) {
  d <- droplevels(d)
  f <- reformulate(c("site_year", n_terms(d), o_term(d), "SHMI_w", mei_terms(d),
                     if (any(d$N_c != 0)) "SHMI_w:N_c", "SHMI_w:age_c", extra), "yield_kgha")
  tryCatch(lme(f, random = ~ 1 | trt/plot_id, correlation = corCAR1(form = ~ year | trt/plot_id),
               data = d, method = "REML", control = ctrl), error = function(e) NULL)
}
ci_of <- function(b, se, df) c(b, b - qt(0.975, df) * se, b + qt(0.975, df) * se)
scale_loo <- map_dfr(c("none", sort(unique(as.character(corn$site)))), function(s_) {
  message("Scale check, leaving out: ", s_)
  dd <- if (s_ == "none") corn else filter(corn, site != s_)
  out <- tibble(left_out = s_, n_sites = n_distinct(dd$site))
  ma <- fit_scale_on(dd, "SHMI_w:pos_c")
  if (!is.null(ma)) {
    t_ <- summary(ma)$tTable["SHMI_w:pos_c", ]
    v <- ci_of(10 * t_[["Value"]], 10 * t_[["Std.Error"]], t_[["DF"]])
    out <- mutate(out, pos_effect = v[1], pos_lower = v[2], pos_upper = v[3], pos_p = t_[["p-value"]])
  }
  mc <- fit_scale_on(dd, "q_abs")
  if (!is.null(mc)) {
    out$curv_p <- summary(mc)$tTable["q_abs", "p-value"]
    for (x in c(40, 60, 80)) {
      s <- slope_curved(mc, "q_abs", x - SHMI_PIV); v <- ci_of(10 * s$b, 10 * s$se, s$df)
      out[[paste0("slope", x)]] <- v[1]; out[[paste0("slope", x, "_lower")]] <- v[2]; out[[paste0("slope", x, "_upper")]] <- v[3]
    }
  }
  out
})
write.csv(mutate(scale_loo, across(where(is.double), ~ signif(.x, 4))),
          file.path(OUT, "SHMI_scale_check_loo.csv"), row.names = FALSE)
print(scale_loo, n = Inf, width = Inf)
loo_only <- filter(scale_loo, left_out != "none")
scale_rows$loo <- tibble(
  figure = "Fig 1, SHMI scale check",
  estimate = c("(a) Leave one site out: smallest and largest position effect",
               "(a) Leave one site out: largest p for the position effect (site left out in the CSV)",
               "(c) Leave one site out: smallest and largest slope at SHMI 80"),
  effect = c(min(loo_only$pos_effect, na.rm = TRUE), NA, min(loo_only$slope80, na.rm = TRUE)),
  se = NA_real_, df = NA_real_,
  lower = NA_real_, upper = c(max(loo_only$pos_effect, na.rm = TRUE), NA, max(loo_only$slope80, na.rm = TRUE)),
  p_val = c(NA, max(loo_only$pos_p, na.rm = TRUE), NA), n = nrow(loo_only))

# (d) between site slopes, and how position relates to system age
site_rows3 <- site_rows %>% left_join(site_means, by = c("label" = "site")) %>% left_join(site_mod, by = c("label" = "site"))
mreg_pos <- rma(yi = effect, sei = se, mods = ~ I((site_SHMI - mean(site_SHMI)) / 10), data = site_rows3,
                method = "REML", test = "knha")
r_pos_age <- cor(site_rows3$site_SHMI, site_rows3$age_mean)
scale_rows$meta <- tibble(
  figure = "Fig 1, SHMI scale check",
  estimate = c("(d) Meta-regression: change in site slope per 10 points higher site mean SHMI",
               "(d) Correlation of site mean SHMI with site mean system age (r)"),
  effect = c(mreg_pos$beta[2], r_pos_age), se = c(mreg_pos$se[2], NA), df = c(mreg_pos$k - mreg_pos$p, NA),
  lower = c(mreg_pos$ci.lb[2], NA), upper = c(mreg_pos$ci.ub[2], NA), p_val = c(mreg_pos$pval[2], NA),
  n = mreg_pos$k)
results$scale <- bind_rows(scale_rows)
write.csv(mutate(results$scale, across(where(is.numeric), ~ signif(.x, 4))),
          file.path(OUT, "SHMI_scale_check.csv"), row.names = FALSE)
print(results$scale, n = Inf, width = Inf)

# ---- each experiment's treatments and stretch of the SHMI scale ----------------------
# One row per site. SHMI range over its treatments (each record counted once, as in
# SHMI_w). share_*: the part of the within-site SHMI variance that comes from each part
# of SHMI (weight x sub-index, or the amendment bonus): cov(part, SHMI) / var(SHMI) over
# the site's records; the four add to 1. The weights are recovered from the scores
# (SHMI_base = wC Cover + wD Diversity + wI InvDist), so overrides are respected.
parts_ok <- all(c("Cover", "Diversity", "InvDist", "SHMI_base", "Amend_bonus") %in% names(corn))
recs <- corn %>% distinct(site, MGT_combo, .keep_all = TRUE)
shares <- if (parts_ok) {
  wfit <- lm(SHMI_base ~ 0 + Cover + Diversity + InvDist, data = recs)
  if (summary(wfit)$r.squared < 0.999) warning("SHMI weights not recovered exactly (R2 = ", round(summary(wfit)$r.squared, 4), ")")
  w_hat <- coef(wfit)
  recs %>% mutate(p_cover = w_hat[["Cover"]] * Cover, p_div = w_hat[["Diversity"]] * Diversity,
                  p_dist = w_hat[["InvDist"]] * InvDist, p_amend = Amend_bonus) %>%
    group_by(site) %>%
    summarise(across(c(p_cover, p_div, p_dist, p_amend),
                     ~ if (var(SHMI) > 1e-8) cov(.x, SHMI) / var(SHMI) else NA_real_,
                     .names = "share_{sub('p_', '', .col)}"), .groups = "drop")
} else tibble(site = character())
site_tab <- corn %>% group_by(site) %>%
  summarise(treatments = n_distinct(MGT_trt_base), plots = n_distinct(plot_id),
            first_year = min(year), last_year = max(year), corn_plot_years = n(),
            N_min = round(100 * min(N100)), N_max = round(100 * max(N100)), .groups = "drop") %>%
  left_join(recs %>% group_by(site) %>%
              summarise(SHMI_min = min(SHMI), SHMI_max = max(SHMI), SHMI_mean = mean(SHMI), SHMI_sd = sd(SHMI),
                        SHMI_levels = n_distinct(round(SHMI)), .groups = "drop"), by = "site") %>%
  left_join(shares, by = "site") %>%
  left_join(site_rows %>% transmute(site = label, slope = effect, slope_lower = effect - qt(0.975, df) * se,
                                    slope_upper = effect + qt(0.975, df) * se), by = "site") %>%
  mutate(location = unname(SITE_LOC[site]), .after = site) %>%
  # each experiment's slope moved to N_REF kg N/ha and the average system age with the
  # Fig. 3a model's SHMI x N and SHMI x age coefficients (as in the heterogeneity section;
  # the interval keeps the site's own width and ignores the adjustment's uncertainty)
  left_join(site_adj %>% transmute(site = label, N_mean = round(N_mean), age_mean = round(age_mean, 1),
                                   slope_adj = effect), by = "site") %>%
  mutate(slope_adj_lower = slope_adj - (slope - slope_lower), slope_adj_upper = slope_adj + (slope_upper - slope)) %>%
  arrange(SHMI_mean)
# What drives each site's SHMI range (from the shares; the rule was set after seeing the
# site table, so groups are descriptive): tillage only (share_dist >= 0.8); rotation,
# legumes or manure (share_div + share_amend >= 0.4); mostly cover crops (share_cover
# >= 0.7); otherwise cover crops and tillage.
CONTRAST_LEVELS <- c("Tillage only", "Mostly cover crops", "Cover crops and tillage", "Rotation, legumes or manure")
if (parts_ok) site_tab <- site_tab %>%
  mutate(contrast = factor(case_when(share_dist >= 0.8 ~ CONTRAST_LEVELS[1],
                                     share_div + share_amend >= 0.4 ~ CONTRAST_LEVELS[4],
                                     share_cover >= 0.7 ~ CONTRAST_LEVELS[2],
                                     TRUE ~ CONTRAST_LEVELS[3]), levels = CONTRAST_LEVELS), .after = share_amend)
write.csv(mutate(site_tab, across(where(is.double), ~ round(.x, 2))), file.path(OUT, "site_SHMI_ranges.csv"), row.names = FALSE)
print(site_tab, n = Inf, width = Inf)

# ---- heterogeneity: the average with site slopes free, and what explains the spread ----
# (e) Random-slope global model: the Fig 1 model with the SHMI slope varying by site
#     (random SHMI_w | site; no random site intercept, which site x year already removes).
#     Its SHMI_w coefficient is the average within-site slope at 150 kg N/ha and average
#     age, allowing each site its own slope. CI on sites - 1 df (nlme's default df treat
#     SHMI_w as a plot-level term and would be too generous with ten sites). Also: the SD
#     of the site slopes (tau), a 95% prediction interval for the slope in a new
#     experiment, and a likelihood-ratio test against the one-slope model (REML, same
#     fixed effects; p halved for the boundary).
# (f) Meta-regressions of the site slopes (each site's own model; and moved to 150 kg
#     N/ha and average age, as site_adj in section 7) on what drives the site's SHMI
#     range: the disturbance share (continuous) and the contrast group, with group means.
#     Descriptive: ten sites, groups defined after seeing the data.
het_x <- list()
m_rs <- tryCatch(
  lme(f_full, random = list(site = pdDiag(~ SHMI_w - 1), trt = ~ 1, plot_id = ~ 1),
      correlation = corCAR1(form = ~ year | site/trt/plot_id), data = corn, method = "REML", control = ctrl),
  error = function(e) { message("Random-slope model did not fit: ", conditionMessage(e)); NULL })
if (!is.null(m_rs)) {
  k_s <- n_distinct(corn$site)
  b_rs <- fixef(m_rs)[["SHMI_w"]]; se_rs <- sqrt(vcov(m_rs)["SHMI_w", "SHMI_w"])
  tau_rs <- as.numeric(VarCorr(m_rs)["SHMI_w", "StdDev"])
  t_k <- qt(0.975, k_s - 1)
  lrt <- tryCatch(anova(m_full, m_rs), error = function(e) NULL)
  p_lrt <- if (!is.null(lrt)) 0.5 * lrt[2, "p-value"] else NA_real_
  het_x$rs <- tibble(
    figure = "Fig 1, heterogeneity",
    estimate = c("(e) Random-slope model: average within-site slope (site slopes free; df = sites - 1)",
                 "(e) Random-slope model: SD of site slopes (tau)",
                 "(e) Random-slope model: 95% prediction interval for a new experiment",
                 "(e) Random-slope model vs one slope: likelihood-ratio test (p halved)"),
    effect = c(10 * b_rs, 10 * tau_rs, 10 * b_rs, if (!is.null(lrt)) lrt[2, "L.Ratio"] else NA),
    se = c(10 * se_rs, NA, NA, NA), df = c(k_s - 1, NA, k_s - 1, 1),
    lower = c(10 * (b_rs - t_k * se_rs), NA, 10 * (b_rs - t_k * sqrt(se_rs^2 + tau_rs^2)), NA),
    upper = c(10 * (b_rs + t_k * se_rs), NA, 10 * (b_rs + t_k * sqrt(se_rs^2 + tau_rs^2)), NA),
    p_val = c(2 * pt(-abs(b_rs / se_rs), k_s - 1), NA, NA, p_lrt), n = k_s)
  # (each site's slope from this model is not reported: with ten sites and site-specific
  # N terms, it is pulled far toward the average, e.g. NELITCSE from -245 to +114)
}
if (parts_ok) {
  mr_dat <- site_rows %>% transmute(site = label, effect, se) %>%
    left_join(select(site_adj, site = label, effect_adj = effect), by = "site") %>%
    left_join(select(site_tab, site, share_dist, share_div, share_amend, contrast), by = "site") %>% droplevels()
  meta_rows <- function(yname, lab) {
    d <- mutate(mr_dat, yi = .data[[yname]])
    m_d <- rma(yi = yi, sei = se, mods = ~ share_dist, data = d, method = "REML", test = "knha")
    out <- tibble(figure = "Fig 1, heterogeneity",
                  estimate = paste0("(f) ", lab, ": change in slope from no disturbance share to all (share_dist 0 to 1)"),
                  effect = m_d$beta[2], se = m_d$se[2], df = m_d$k - m_d$p, lower = m_d$ci.lb[2], upper = m_d$ci.ub[2],
                  p_val = m_d$pval[2], n = m_d$k)
    # 9 Nov 2026: the return per point against the share of each experiment's SHMI spread that
    # came from crop diversity (Fig. 3c is ordered by it); per 10 percentage points of share.
    # Sensitivity: diversity plus amendment share (WICST's manure).
    for (mv in c("share_div", "share_divam")) {
      d$share_divam <- d$share_div + d$share_amend
      m_v <- rma(yi = yi, sei = se, mods = as.formula(paste("~ I(", mv, "* 10)")), data = d, method = "REML", test = "knha")
      out <- bind_rows(out, tibble(figure = "Fig 1, heterogeneity",
        estimate = paste0("(f) ", lab, c(": change in slope per 10 percentage points of ", ": share of between-site variance explained by "),
                          if (mv == "share_div") "diversity share" else "diversity + amendment share", c("", " (R2, %)")),
        effect = c(m_v$beta[2], m_v$R2), se = c(m_v$se[2], NA), df = c(m_v$k - m_v$p, NA),
        lower = c(m_v$ci.lb[2], NA), upper = c(m_v$ci.ub[2], NA), p_val = c(m_v$pval[2], NA), n = m_v$k))
    }
    if (n_distinct(d$contrast) > 1) {
      m_g  <- rma(yi = yi, sei = se, mods = ~ contrast, data = d, method = "REML", test = "knha")
      m_g0 <- rma(yi = yi, sei = se, mods = ~ contrast - 1, data = d, method = "REML", test = "knha")
      m_0  <- rma(yi = yi, sei = se, data = d, method = "REML", test = "knha")
      out <- bind_rows(out,
        tibble(figure = "Fig 1, heterogeneity",
               estimate = paste0("(f) ", lab, c(": contrast groups, omnibus test (F)",
                                                ": I2 (%) without and with contrast groups",
                                                ": share of between-site variance explained by contrast groups (R2, %)")),
               effect = c(m_g$QM, m_0$I2, m_g$R2), se = NA_real_,
               df = c(if (is.null(m_g$QMdf)) m_g$p - 1 else m_g$QMdf[1], NA, NA),
               lower = NA_real_, upper = c(NA, m_g$I2, NA), p_val = c(m_g$QMp, NA, NA), n = m_g$k),
        tibble(figure = "Fig 1, heterogeneity",
               estimate = paste0("(f) ", lab, ": mean slope, ", sub("^contrast", "", rownames(m_g0$beta)),
                                 " (", as.integer(table(d$contrast)[sub("^contrast", "", rownames(m_g0$beta))]), " sites)"),
               effect = as.numeric(m_g0$beta), se = m_g0$se, df = m_g0$k - m_g0$p,
               lower = m_g0$ci.lb, upper = m_g0$ci.ub, p_val = m_g0$pval, n = m_g0$k))
    }
    out
  }
  het_x$meta <- bind_rows(meta_rows("effect", "Site slopes as estimated"),
                          meta_rows("effect_adj", sprintf("Site slopes at %d kg N/ha and average age", N_REF)))
}
results$het_x <- bind_rows(het_x)
write.csv(mutate(results$het_x, across(where(is.numeric), ~ signif(.x, 4))),
          file.path(OUT, "SHMI_heterogeneity.csv"), row.names = FALSE)
print(results$het_x, n = Inf, width = Inf)
write.csv(mutate(site_tab, across(where(is.double), ~ round(.x, 2))), file.path(OUT, "site_SHMI_ranges.csv"), row.names = FALSE)

# ---- Extended Data Fig. 4: each experiment's stretch of the SHMI scale (a) and slope (b) ----
# Sites in the same order in both panels: grouped by what drives their SHMI range
# (CONTRAST_LEVELS, top to bottom: rotation, legumes or manure; cover crops and tillage;
# mostly cover crops; tillage only), by mean SHMI within a group. Colours = groups.
# (a) Fig_SHMI_ranges.png: each treatment record's SHMI (points) and the site's range
#     (bar); together the sites span the scale that no single experiment does.
# (b) EDFig7b_forest.png: the rows of Fig4_forest.png without the averages of site slopes
#     and without stars; summary rows are the global slope (Fig 1) and the random-slope
#     average (e); site labels carry each site's SHMI range.
CONTRAST_COLS <- setNames(c("#D55E00", "#E69F00", "#56B4E9", "#0072B2"), CONTRAST_LEVELS)
has_grp <- "contrast" %in% names(site_tab)
site_seq <- if (has_grp) site_tab %>% arrange(desc(as.integer(contrast)), SHMI_mean) %>% pull(site) else
  site_tab %>% arrange(SHMI_mean) %>% pull(site)
rng_dat <- site_tab %>% mutate(site = factor(site, levels = rev(site_seq)),     # first group at the top
                               grp = if (has_grp) contrast else factor("all"))
trt_pts <- recs %>% distinct(site = as.character(site), SHMI = round(SHMI, 1)) %>%
  left_join(transmute(rng_dat, site = as.character(site), grp), by = "site") %>%
  mutate(site = factor(site, levels = levels(rng_dat$site)))
grp_scale <- if (has_grp) list(scale_colour_manual(values = CONTRAST_COLS, name = NULL, drop = FALSE),
                               scale_fill_manual(values = CONTRAST_COLS, name = NULL, drop = FALSE)) else
  list(scale_colour_manual(values = "grey40", guide = "none"), scale_fill_manual(values = "grey40", guide = "none"))
fig_rng <- ggplot(rng_dat, aes(y = site)) +
  geom_segment(aes(x = SHMI_min, xend = SHMI_max, yend = site, colour = grp), linewidth = 2.2, alpha = 0.45) +
  geom_point(data = trt_pts, aes(x = SHMI, fill = grp), shape = 21, size = 2.2, colour = "black") +
  grp_scale +
  scale_y_discrete(labels = function(s) paste0(s, " (", SITE_LOC[s], ")")) +
  scale_x_continuous(limits = c(floor(min(rng_dat$SHMI_min) / 10) * 10, 100), breaks = seq(20, 100, 10)) +
  guides(colour = "none", fill = guide_legend(nrow = 2, override.aes = list(size = 3))) +
  labs(x = XLAB, y = NULL, title = "Each experiment's stretch of the SHMI scale",
       caption = paste0("Points: SHMI of each treatment record; bars: the site's range; colours: what drives the site's SHMI range.\n",
                        "The global model estimates one slope from the differences within each site; together the sites span SHMI ",
                        sprintf("%.0f-%.0f.", min(rng_dat$SHMI_min), max(rng_dat$SHMI_max)))) +
  theme_fig + theme(panel.grid.major.y = element_blank(), legend.position = "bottom")
ggsave(file.path(OUT, "Fig_SHMI_ranges.png"), fig_rng + MS, width = 6.8, height = 1.6 + 0.32 * nrow(rng_dat), dpi = 300)

# Extended Data Fig. 4a (since 30 Oct 2026): the ranges beside each experiment's own slope,
# same rows and colours; one point per experiment (no sub-rows), with the global slope
# (Fig. 3a) as a dashed line. Replaces the forest plot (EDFig7b_forest.png, still written).
if (all(c("slope", "slope_lower", "slope_upper") %in% names(rng_dat)) && requireNamespace("patchwork", quietly = TRUE)) {
  XLIM <- c(-1000, 1250); XLAB_EDGE <- 1500      # values printed right of the panel's data range
  # ED5A_COMMON (default TRUE, 6 Nov 2026): plot each experiment's slope at N_REF kg N/ha
  # and the average system age (filled, with its interval), the slope as estimated in the
  # experiment's own model shown as an open grey circle; FALSE restores the earlier panel.
  if (!exists("ED5A_COMMON")) ED5A_COMMON <- TRUE
  common <- ED5A_COMMON && "slope_adj" %in% names(rng_dat)
  sl_dat <- rng_dat %>% mutate(x  = if (common) slope_adj else slope,
                               lo = if (common) slope_adj_lower else slope_lower,
                               hi = if (common) slope_adj_upper else slope_upper,
                               lo_c = pmax(lo, XLIM[1]), hi_c = pmin(hi, XLIM[2]),
                               cut = lo < XLIM[1] | hi > XLIM[2])
  fig_sl <- ggplot(sl_dat, aes(y = site)) +
    geom_vline(xintercept = 0, colour = "grey60") +
    geom_vline(xintercept = 10 * s0$b, colour = "darkgreen", linetype = "dashed") +
    geom_errorbar(aes(xmin = lo_c, xmax = hi_c, colour = grp), width = 0, orientation = "y") +
    { if (common) geom_point(aes(x = pmin(pmax(slope, XLIM[1]), XLIM[2])), shape = 21, size = 2.2,
                             colour = "grey45", fill = "white") } +
    geom_point(aes(x = x, fill = grp), shape = 21, size = 2.6, colour = "black") +
    geom_text(aes(x = XLAB_EDGE, label = sprintf("%+.0f", x)), hjust = 1, size = 3) +
    grp_scale + guides(colour = "none", fill = "none") +
    scale_x_continuous(limits = c(XLIM[1], XLAB_EDGE), breaks = seq(-1000, 1000, 500)) +
    labs(x = if (common) expression(atop("Each experiment's own slope at 150 kg N ha"^-1*" and average age", "(kg ha"^-1*" per 10 SHMI; open: as estimated)")) else
               expression("Each experiment's own slope (kg ha"^-1*" per 10 SHMI)"), y = NULL) +
    theme_fig + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major.y = element_blank())
  fig_rs <- (fig_rng + MS + labs(title = NULL) | fig_sl) + patchwork::plot_layout(widths = c(1.15, 1), guides = "collect") &
    theme(legend.position = "bottom")
  ggsave(file.path(OUT, "Fig_SHMI_ranges_slopes.png"), fig_rs, width = 11, height = 1.9 + 0.32 * nrow(rng_dat), dpi = 300)
  if (any(sl_dat$cut)) message("Extended Data Fig. 4a: CI cut at the panel edge for ", paste(sl_dat$site[sl_dat$cut], collapse = ", "))
}

# ---- Corn-every-year systems: yield per SHMI point by route (7 Nov 2026, exploratory) -----
# Systems whose only harvested crop is corn (cover crops, living mulches and manure allowed:
# they have no yield), so rotation cannot explain the slope. Within each experiment, SHMI is
# split into the points from Inverse Disturbance (p_dist = wI x InvDist) and the rest (living
# cover, diversity and the amendment bonus: SHMI - p_dist), each centred on the experiment's
# mean over its units, as SHMI_w is. One model for all such plots, with the terms of each
# experiment's own model (site x year, N terms, other factors) and a slope for each route;
# slopes per 10 points at the data's own N and age. Written: SHMI_route_continuous_corn.csv
# (Its objects are suffixed _cc so the main model's b, V and tt are not overwritten.)
if (parts_ok && exists("runs")) {
  cc_units <- imap_dfr(runs, function(r, s_) {
    y <- r$inputs$yield
    if (is.null(y) || !nrow(y)) return(NULL)
    y %>% filter(!is.na(yield_kg_ha)) %>% distinct(MGT_combo, crop = tolower(trimws(CD_name))) %>%
      group_by(MGT_combo) %>% summarise(cc = all(grepl("^corn", crop)), .groups = "drop") %>% mutate(site = s_)
  })
  dcc <- corn %>% left_join(cc_units, by = c("site", "MGT_combo")) %>% filter(cc %in% TRUE) %>%
    mutate(p_dist = w_hat[["InvDist"]] * InvDist, p_rest = SHMI - p_dist) %>%
    group_by(site) %>%
    mutate(dist_w = p_dist - mean(p_dist[!duplicated(MGT_combo)]),
           rest_w = p_rest - mean(p_rest[!duplicated(MGT_combo)])) %>%
    filter(n_distinct(MGT_trt_base) >= 2, diff(range(SHMI)) >= 2) %>% ungroup() %>% droplevels()
  cc_sites <- dcc %>% group_by(site) %>%
    summarise(treatments = n_distinct(MGT_trt_base), plot_years = n(),
              SHMI_range = sprintf("%.0f-%.0f", min(SHMI), max(SHMI)),
              dist_range = round(diff(range(p_dist)), 1), rest_range = round(diff(range(p_rest)), 1), .groups = "drop")
  print(cc_sites)
  m_cc <- tryCatch(fit_site(dcc, x = c("dist_w", "rest_w")),
                   error = function(e) { message("corn-every-year route model failed: ", conditionMessage(e)); NULL })
  if (!is.null(m_cc)) {
    tt_cc <- summary(m_cc)$tTable; b_cc <- fixef(m_cc); V_cc <- vcov(m_cc)
    row_of <- function(k, lab) tibble(route = lab, effect = 10 * tt_cc[k, "Value"], se = 10 * tt_cc[k, "Std.Error"],
                                      df = tt_cc[k, "DF"], p = tt_cc[k, "p-value"]) %>%
      mutate(lower = effect - qt(0.975, df) * se, upper = effect + qt(0.975, df) * se)
    dlt <- b_cc[["rest_w"]] - b_cc[["dist_w"]]
    se_d <- sqrt(V_cc["rest_w", "rest_w"] + V_cc["dist_w", "dist_w"] - 2 * V_cc["rest_w", "dist_w"])
    df_d <- min(tt_cc[c("dist_w", "rest_w"), "DF"])
    route_tab <- bind_rows(row_of("dist_w", "Inverse Disturbance points (less tillage)"),
                           row_of("rest_w", "Other points (living cover, diversity, amendments)"),
                           tibble(route = "Difference (other minus disturbance)", effect = 10 * dlt, se = 10 * se_d, df = df_d,
                                  p = 2 * pt(-abs(dlt / se_d), df_d)) %>%
                             mutate(lower = effect - qt(0.975, df) * se, upper = effect + qt(0.975, df) * se)) %>%
      mutate(experiments = paste(cc_sites$site, collapse = ", "), plot_years = nrow(dcc))
    write.csv(mutate(route_tab, across(where(is.double), ~ signif(.x, 4))), file.path(OUT, "SHMI_route_continuous_corn.csv"), row.names = FALSE)
    write.csv(cc_sites, file.path(OUT, "SHMI_route_continuous_corn_sites.csv"), row.names = FALSE)
    cat("\nCorn-every-year systems: yield change per 10 SHMI points by route:\n"); print(as.data.frame(route_tab), digits = 3)
  }
}

# ---- Fig. 4 (7 Nov 2026): what the SHMI yield benefit depends on -----------------------------
# a the slope against N (Fig2_SHMI_x_N), b against system age (Fig3_SHMI_x_age), c each
# experiment's own slope at N_REF kg N/ha and the average age (slope_adj; filled: 95% CI
# excludes zero), grouped by what drove its SHMI range, with each group's mean (meta-analysis
# of the adjusted slopes) as a line and band; NEMEIRR also by share of stover removed
# (diamonds). All three on one y-axis (yield change per 10 SHMI); the dashed green line in c
# is the slope of Fig. 3a. FIG4_SAME_Y (default TRUE): one y range for a-c; experiment
# intervals beyond it are cut at the edge (arrow). FIG4_YLIM: set c(lo, hi) to choose it.
# Written: Fig4_moderators.png, lettered a-c (B7 copies it as Fig. 4). Replaces Extended Data Fig. 5a.
if (has_grp && "slope_adj" %in% names(site_tab) && requireNamespace("patchwork", quietly = TRUE)) {
  if (!exists("FIG4_SAME_Y")) FIG4_SAME_Y <- TRUE
  if (!exists("FIG3C_LEVELS")) FIG3C_LEVELS <- FALSE   # TRUE: NEMEIRR's slope by stover removal as diamonds
  oth_sites <- if (exists("other_tab") && nrow(other_tab)) unique(other_tab$site[!is.na(other_tab$slope)]) else character()
  GRP_SHORT <- setNames(c("Tillage\nonly", "Mostly\ncover crops", "Cover crops\nand tillage", "Rotation,\nlegumes, manure"),
                        CONTRAST_LEVELS)
  c_dat <- site_tab %>% filter(!is.na(contrast)) %>%
    mutate(sig = slope_adj_lower > 0 | slope_adj_upper < 0,
           lab = sprintf("%s (%.0f–%.0f)%s", site, SHMI_min, SHMI_max, ifelse(site %in% oth_sites, " \u2020", "")),
           grp = factor(GRP_SHORT[as.character(contrast)], levels = GRP_SHORT)) %>%
    arrange(contrast, slope_adj) %>% mutate(lab = factor(lab, levels = lab))
  g_mean <- results$het_x %>%
    filter(grepl(sprintf("Site slopes at %d kg N/ha and average age: mean slope, ", N_REF), estimate, fixed = TRUE)) %>%
    mutate(contrast = sub(" \\(.*$", "", sub("^.*mean slope, ", "", estimate))) %>%
    filter(contrast %in% CONTRAST_LEVELS) %>%
    transmute(grp = factor(GRP_SHORT[contrast], levels = GRP_SHORT), est = effect, lo = lower, hi = upper)
  # NEMEIRR by stover removed: its own SHMI x residue rows, moved by NEMEIRR's N and age shift
  mei <- forest %>% filter(FIG3C_LEVELS, grepl("^NEMEIRR, no-till vs disk, \\d+% residue removed", label))
  if (nrow(mei) && "NEMEIRR" %in% c_dat$site) {
    nm <- filter(c_dat, site == "NEMEIRR")
    mei <- mei %>% mutate(pct = sub("^.*disk, (\\d+)% residue removed.*$", "\\1", label),
                          y = effect + (nm$slope_adj - nm$slope), lab = nm$lab, grp = nm$grp)
  } else mei <- NULL
  # the moderator panels' bands and the experiments' point estimates set the shared range
  band <- function(xs, at) map_dfr(xs, function(x) { s <- at(x); tibble(b = 10 * s$b, se = 10 * s$se) }) %>%
    mutate(lo = b - qt(0.975, DF_W) * se, hi = b + qt(0.975, DF_W) * se)
  bands <- bind_rows(band(seq(N_FIG_MIN, max(N_LINES), length.out = 50), function(n) slope_at(N = n)),
                     band(seq(1, max(corn$age), length.out = 50), function(a) slope_at(age = a)))
  PD$fig4c <- list(sites = c_dat %>% left_join(select(site_rows2, site = label, se), by = "site"),
                   groups = g_mean, nemeirr = mei, bands = bands, global_slope = 10 * s0$b,
                   group_short = GRP_SHORT, contrast_cols = CONTRAST_COLS)
  YL <- if (exists("FIG4_YLIM")) FIG4_YLIM else {
    r <- range(c(bands$lo, bands$hi, c_dat$slope_adj, g_mean$lo, g_mean$hi, 0), na.rm = TRUE)   # NEMEIRR's diamonds may fall outside (drawn at the edge)
    c(floor(r[1] / 100) * 100 - 50, ceiling(r[2] / 100) * 100 + 50)
  }
  c_dat <- c_dat %>% mutate(lo_c = pmax(slope_adj_lower, YL[1]), hi_c = pmin(slope_adj_upper, YL[2]),
                            cut_lo = slope_adj_lower < YL[1], cut_hi = slope_adj_upper > YL[2])
  if (!is.null(mei)) mei <- mei %>% mutate(cut = y < YL[1] | y > YL[2], y = pmin(pmax(y, YL[1]), YL[2]),
                                           tag = ifelse(cut, paste0(pct, "% \u2193"), paste0(pct, "%")))   # value of a cut diamond: in the legend
  YT <- "Corn yield change\nper 10 SHMI (kg ha⁻¹)"
  GCOLS <- setNames(CONTRAST_COLS[CONTRAST_LEVELS], GRP_SHORT)
  if (!exists("FIG3C_ORDER")) FIG3C_ORDER <- "diversity"   # "diversity" (default; no fitted line since 10 Nov 2026), "scatter" or "groups" (earlier)
  if (FIG3C_ORDER == "scatter") {
    # Each experiment's slope against the share of its SHMI spread that came from crop
    # diversity, on a continuous axis, with the meta-regression of slope on that share
    # (precision-weighted, REML, Knapp-Hartung) as one straight line and its 95% band.
    # Points shaded by the disturbance share; filled if the experiment's CI excludes zero.
    dd <- c_dat %>% left_join(select(site_rows2, site = label, se), by = "site") %>%
      mutate(x = 100 * share_div)
    m_div <- rma(yi = slope_adj, sei = se, mods = ~ share_div, data = dd, method = "REML", test = "knha")
    xs <- seq(0, max(dd$share_div), length.out = 100)
    pr <- predict(m_div, newmods = xs)
    ln <- tibble(x = 100 * xs, fit = pr$pred, lo = pmax(pr$ci.lb, YL[1]), hi = pmin(pr$ci.ub, YL[2]))
    # leave-one-experiment-out: does the line depend on any one experiment?
    loo <- map_dfr(dd$site, function(s) {
      m <- rma(yi = slope_adj, sei = se, mods = ~ share_div, data = filter(dd, site != s), method = "REML", test = "knha")
      tibble(left_out = s, per10 = m$beta[2] / 10, lower = m$ci.lb[2] / 10, upper = m$ci.ub[2] / 10, p = m$pval[2])
    })
    write.csv(loo, file.path(OUT, "Fig3c_meta_loo.csv"), row.names = FALSE)
    cat(sprintf("Fig. 3c meta-regression, leave one out: %.0f to %.0f per 10 points of share, P %.3f to %.3f\n",
                min(loo$per10), max(loo$per10), min(loo$p), max(loo$p)))
    pr_s <- predict(m_div, newmods = dd$share_div)
    write.csv(dd %>% transmute(site, share_cover, share_div, share_dist, share_amend, slope_adj, slope_adj_lower, slope_adj_upper, se,
                               pred = pr_s$pred, pred_lo = pr_s$ci.lb, pred_hi = pr_s$ci.ub),
              file.path(OUT, "Fig3c_sites_by_diversity.csv"), row.names = FALSE)
    lab_m <- sprintf("%+.0f per 10 points of diversity share (95%% CI %.0f to %.0f), P = %s; R\u00B2 = %.0f%%",
                     m_div$beta[2] / 10, m_div$ci.lb[2] / 10, m_div$ci.ub[2] / 10,
                     format.pval(m_div$pval[2], digits = 2, eps = 0.001), m_div$R2)
    dd <- dd %>% mutate(y_pt = pmin(pmax(slope_adj, YL[1]), YL[2]), cut_pt = slope_adj < YL[1] | slope_adj > YL[2])
    DCOL <- c("#56B4E9", "#D55E00")
    rep_ok <- requireNamespace("ggrepel", quietly = TRUE)
    lab_layer <- if (rep_ok) ggrepel::geom_text_repel(aes(y = y_pt, label = lab), size = 2.3, colour = "grey25",
                                                      min.segment.length = 0.2, box.padding = 0.35, seed = 1, max.overlaps = Inf)
                 else geom_text(aes(y = y_pt, label = lab), size = 2.3, colour = "grey25", hjust = -0.12, vjust = -0.6)
    fig_c <- ggplot(dd, aes(x, slope_adj)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
      geom_hline(yintercept = 10 * s0$b, colour = "darkgreen", linetype = "dashed") +
      geom_ribbon(data = ln, inherit.aes = FALSE, aes(x, ymin = lo, ymax = hi), fill = "grey50", alpha = 0.18) +
      geom_line(data = ln, inherit.aes = FALSE, aes(x, fit), colour = "grey20", linewidth = 0.9) +
      geom_linerange(aes(ymin = lo_c, ymax = hi_c, colour = pmax(share_dist, 0)), linewidth = 0.6) +
      geom_point(data = filter(dd, cut_hi), aes(y = YL[2], colour = pmax(share_dist, 0)), shape = 17, size = 1.8) +
      geom_point(data = filter(dd, cut_lo), aes(y = YL[1], colour = pmax(share_dist, 0)), shape = 25, size = 1.8, fill = "grey40") +
      geom_point(aes(y = y_pt, colour = pmax(share_dist, 0), fill = ifelse(sig, pmax(share_dist, 0), NA),
                     shape = ifelse(cut_pt, "cut", "in")), size = 2.8, stroke = 1) +
      scale_shape_manual(values = c(`in` = 21, cut = 25), guide = "none") +
      lab_layer +
      annotate("text", x = 0, y = YL[2], hjust = 0, vjust = 1, size = 2.5, label = lab_m) +
      scale_x_continuous(expand = expansion(mult = c(0.03, 0.06))) +
      scale_colour_gradient(low = DCOL[1], high = DCOL[2], limits = c(0, 1), name = "Disturbance\nshare") +
      scale_fill_gradient(low = DCOL[1], high = DCOL[2], limits = c(0, 1), na.value = "white", guide = "none") +
      labs(x = "Share of the experiment's SHMI spread from Diversity (%)", y = YT) +
      theme_fig + theme(panel.grid.minor = element_blank(),
                        legend.position = c(0.99, 0.02), legend.justification = c(1, 0),
                        legend.key.height = unit(8, "pt"), legend.title = element_text(size = 7),
                        legend.text = element_text(size = 6.5), legend.background = element_blank())
  } else if (FIG3C_ORDER == "diversity") {
    # Experiments ordered by the share of their SHMI spread that came from crop diversity
    # (computed in this run, so the order follows the data); points shaded by the disturbance
    # share; line and band: meta-regression of the slope on the diversity share (precision-
    # weighted, Knapp-Hartung), evaluated at each experiment's share.
    dd <- c_dat %>% left_join(select(site_rows2, site = label, se), by = "site") %>%
      arrange(share_div, slope_adj) %>% mutate(i = row_number())
    # The meta-regression of slope on diversity share is reported in the text only (no line
    # on the categorical axis, 10 Nov 2026); it and its leave-one-out refits are written here.
    m_div <- rma(yi = slope_adj, sei = se, mods = ~ share_div, data = dd, method = "REML", test = "knha")
    loo <- map_dfr(dd$site, function(s) {
      m <- rma(yi = slope_adj, sei = se, mods = ~ share_div, data = filter(dd, site != s), method = "REML", test = "knha")
      tibble(left_out = s, per10 = m$beta[2] / 10, lower = m$ci.lb[2] / 10, upper = m$ci.ub[2] / 10, p = m$pval[2])
    })
    write.csv(loo, file.path(OUT, "Fig3c_meta_loo.csv"), row.names = FALSE)
    cat(sprintf("Fig. 3c meta-regression: %+.0f per 10 points of diversity share (%.0f to %.0f), P = %.3f, R2 = %.0f%%; leave one out %.0f to %.0f, P %.3f to %.3f\n",
                m_div$beta[2] / 10, m_div$ci.lb[2] / 10, m_div$ci.ub[2] / 10, m_div$pval[2], m_div$R2,
                min(loo$per10), max(loo$per10), min(loo$p), max(loo$p)))
    pr <- predict(m_div, newmods = dd$share_div)
    dd <- dd %>% mutate(pred = pr$pred, pred_lo = pmax(pr$ci.lb, YL[1]), pred_hi = pmin(pr$ci.ub, YL[2]),
                        xlab = sprintf("%s\n%.0f%%", as.character(lab), 100 * share_div))
    cat(sprintf("Fig. 3c order (diversity share): %s\n", paste(sprintf("%s %.0f%%", dd$site, 100 * dd$share_div), collapse = ", ")))
    lab_m <- sprintf("Slope vs diversity share: %+.0f per 10 points of share (95%% CI %.0f to %.0f), P = %s",
                     m_div$beta[2] / 10, m_div$ci.lb[2] / 10, m_div$ci.ub[2] / 10, format.pval(m_div$pval[2], digits = 2, eps = 0.001))
    write.csv(dd %>% transmute(site, share_cover, share_div, share_dist, share_amend, slope_adj, slope_adj_lower, slope_adj_upper, se,
                               pred = pr$pred, pred_lo = pr$ci.lb, pred_hi = pr$ci.ub),
              file.path(OUT, "Fig3c_sites_by_diversity.csv"), row.names = FALSE)
    DCOL <- c("#56B4E9", "#D55E00")
    fig_c <- ggplot(dd, aes(i, slope_adj)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
      geom_hline(yintercept = 10 * s0$b, colour = "darkgreen", linetype = "dashed") +
      geom_linerange(aes(ymin = lo_c, ymax = hi_c, colour = pmax(share_dist, 0)), linewidth = 0.6) +
      geom_point(data = filter(dd, cut_hi), aes(y = YL[2], colour = pmax(share_dist, 0)), shape = 17, size = 1.8) +
      geom_point(data = filter(dd, cut_lo), aes(y = YL[1], colour = pmax(share_dist, 0)), shape = 25, size = 1.8,
                 fill = "grey40") +
      geom_point(aes(colour = pmax(share_dist, 0), fill = ifelse(sig, pmax(share_dist, 0), NA)), shape = 21, size = 2.8, stroke = 1) +
      scale_x_continuous(breaks = dd$i, labels = dd$xlab, expand = expansion(add = 0.5)) +
      scale_colour_gradient(low = DCOL[1], high = DCOL[2], limits = c(0, 1), name = "Disturbance\nshare") +
      scale_fill_gradient(low = DCOL[1], high = DCOL[2], limits = c(0, 1), na.value = "white", guide = "none") +
      labs(x = "Experiment (SHMI range), ordered by the diversity share of its SHMI spread", y = YT) +
      theme_fig + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
                        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
                        legend.position = c(0.99, 0.02), legend.justification = c(1, 0),
                        legend.key.height = unit(8, "pt"), legend.title = element_text(size = 7),
                        legend.text = element_text(size = 6.5), legend.background = element_blank())
  } else {
  fig_c <- ggplot(c_dat, aes(lab, slope_adj)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
      geom_hline(yintercept = 10 * s0$b, colour = "darkgreen", linetype = "dashed") +
      geom_rect(data = g_mean, aes(xmin = -Inf, xmax = Inf, ymin = pmax(lo, YL[1]), ymax = pmin(hi, YL[2]), fill = grp),
                inherit.aes = FALSE, alpha = 0.15) +
      geom_hline(data = g_mean, aes(yintercept = est, colour = grp), linewidth = 0.9) +
      geom_linerange(aes(ymin = lo_c, ymax = hi_c, colour = grp), linewidth = 0.6) +
      geom_point(data = filter(c_dat, cut_hi), aes(y = YL[2], fill = grp), shape = 24, size = 1.8, colour = NA) +
      geom_point(data = filter(c_dat, cut_lo), aes(y = YL[1], fill = grp), shape = 25, size = 1.8, colour = NA) +
      geom_point(aes(colour = grp, fill = ifelse(sig, as.character(grp), NA)), shape = 21, size = 2.8, stroke = 1) +
      { if (!is.null(mei)) list(
          geom_point(data = mei, aes(lab, y), inherit.aes = FALSE, position = position_nudge(x = 0.32),
                     shape = 23, size = 1.9, colour = GCOLS[[1]], fill = "white", stroke = 0.7),
          geom_text(data = mei, aes(lab, y, label = tag), inherit.aes = FALSE,
                    position = position_nudge(x = 0.52), hjust = 0, size = 2.3, colour = "grey30")) } +
      facet_grid(. ~ grp, scales = "free_x", space = "free_x") +
      scale_colour_manual(values = GCOLS, guide = "none") +
      scale_fill_manual(values = GCOLS, na.value = "white", guide = "none") +
      labs(x = NULL, y = YT) +
      theme_fig + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5),
                        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
                        strip.background = element_rect(fill = "grey95"), strip.text = element_text(size = 7.5),
                        panel.spacing.x = unit(2, "pt"))
  }
  same_y <- if (FIG4_SAME_Y) coord_cartesian(ylim = YL) else NULL
  # 10 Nov 2026: a and b side by side above c (c full width); FIG4_ONE_ROW = TRUE for the old 1 x 3 row
  if (!exists("FIG4_ONE_ROW")) FIG4_ONE_ROW <- FALSE
  pc <- fig_c + (if (FIG4_SAME_Y) coord_cartesian(ylim = YL, clip = "off") else NULL)
  if (FIG4_ONE_ROW) {
    fig4_all <- patchwork::wrap_plots(fig4 + MS + same_y + labs(y = YT), fig3 + MS + same_y + labs(y = NULL),
                                      pc + labs(y = NULL), nrow = 1, widths = c(1, 1, 1.35))
    f4w <- 13; f4h <- 4.6
  } else {
    fig4_all <- (fig4 + MS + same_y + labs(y = YT) | fig3 + MS + same_y + labs(y = NULL)) / (pc + labs(y = YT)) +
      patchwork::plot_layout(heights = c(1, 1.15))
    f4w <- 9; f4h <- 8.4
  }
  fig4_all <- fig4_all + patchwork::plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = "bold", size = 14))
  ggsave(file.path(OUT, "Fig4_moderators.png"), fig4_all, width = f4w, height = f4h, dpi = 300, bg = "white")
  cat(sprintf("Fig. 4: y range %.0f to %.0f kg/ha per 10 SHMI; experiments with intervals cut: %s\n", YL[1], YL[2],
              paste(c_dat$site[c_dat$cut_lo | c_dat$cut_hi], collapse = ", ")))
}

# (b) the forest plot for Extended Data Fig. 4
fr_sum <- forest %>% filter(kind == "within") %>% select(label, kind, effect, lower, upper, display) %>%
  mutate(display = "Global slope (one model, all sites)")
if (!is.null(het_x$rs)) fr_sum <- bind_rows(fr_sum, tibble(
  label = "Random-slope average", kind = "pooled", effect = het_x$rs$effect[1], lower = het_x$rs$lower[1],
  upper = het_x$rs$upper[1], display = "Average, each site its own slope (random slopes)"))
rng_lab <- with(site_tab, setNames(sprintf("%.0f\u2013%.0f", SHMI_min, SHMI_max), site))
fr_sites <- forest %>% filter(kind %in% c("site", "sub")) %>% select(label, kind, effect, lower, upper, display, site) %>%
  mutate(.ord = seq_len(n())) %>% arrange(factor(site, levels = site_seq), .ord) %>% select(-.ord) %>%
  mutate(display = ifelse(kind == "site", paste0(display, ", SHMI ", rng_lab[site]), display)) %>%
  left_join(transmute(rng_dat, site = as.character(site), grp), by = "site")
fr <- bind_rows(fr_sum, fr_sites) %>%
  mutate(row = rev(seq_len(n())),
         fcat = case_when(kind == "within" ~ "Global slope", kind == "pooled" ~ "Random-slope average",
                          kind == "sub" ~ "Sub-row", TRUE ~ as.character(grp))) %>%
  left_join(select(sty, kind, shape, size), by = "kind")
fcols <- c("Global slope" = "darkgreen", "Random-slope average" = "firebrick", "Sub-row" = "white",
           if (has_grp) CONTRAST_COLS else c(all = "steelblue"))
write.csv(select(fr, -row), file.path(OUT, "EDFig7b_forest_rows.csv"), row.names = FALSE)
fig_fr <- ggplot(fr, aes(effect, row)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_hline(yintercept = max(fr$row) - nrow(fr_sum) + 0.5, colour = "grey85") +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.25, colour = "grey40", orientation = "y") +
  geom_point(aes(shape = shape, size = size, fill = fcat), colour = "black", stroke = 0.4) +
  scale_shape_identity() + scale_size_identity() +
  scale_fill_manual(values = fcols, breaks = if (has_grp) CONTRAST_LEVELS else NULL, name = NULL) +
  guides(fill = guide_legend(nrow = 2, override.aes = list(shape = 21, size = 3))) +
  scale_y_continuous(breaks = fr$row, labels = fr$display, expand = expansion(add = 0.6)) +
  labs(x = "Corn yield change per 10 SHMI (kg ha\u207B\u00B9)", y = NULL,
       caption = paste0("Global slope: the Fig 1 model at ", N_REF, " kg N/ha and average age. Random-slope average: the same model with the\n",
                        "slope free to differ among sites. Site rows: each site's own model; open circles: sub-rows as in Fig4_forest.png.")) +
  theme_fig + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), legend.position = "bottom")
ggsave(file.path(OUT, "EDFig7b_forest.png"),
       fig_fr + MS + force_panelsizes(rows = unit(0.22 * nrow(fr), "in"), cols = unit(3.4, "in")),
       width = 6.8, height = 2.0 + 0.22 * nrow(fr), dpi = 300)


saveRDS(PD, file.path(OUT, "plotdata_yield_main.rds"))
# the fitted main model and what B4 needs with it (11 Nov 2026), so B4 no longer reruns B2
saveRDS(list(m_full = m_full, corn = corn, n_var_sites = n_var_sites, n_curve_sites = n_curve_sites,
             N_REF = N_REF, AGE_MEAN = AGE_MEAN, OUT = OUT),
        file.path(OUT, "B2_model.rds"))
message("Plot data for the figure scripts: ", file.path(OUT, "plotdata_yield_main.rds"))


################################################################################
# 8. RESULTS TABLE
################################################################################
results$sites <- site_rows %>% transmute(figure = "Fig 4", estimate = paste("Site:", label), effect, se, df,
                                        lower = effect - qt(0.975, df) * se, upper = effect + qt(0.975, df) * se,
                                        p_val, n = 1)
results$pooled <- tibble(figure = "Fig 4", estimate = "Pooled across sites (random effects, Knapp-Hartung)",
                         effect = as.numeric(pool$b), se = pool$se, df = pool$k - 1, lower = pool$ci.lb,
                         upper = pool$ci.ub, p_val = pool$pval, n = pool$k)
res <- bind_rows(results) %>%
  arrange(factor(sub("[ ,(].*$", "", sub("^Fig ", "", figure)), levels = c("1", "2", "3", "4", "5"))) %>%  # talk order
  mutate(across(c(effect, se, lower, upper), ~ round(.x, 1)),
                                     df = round(df, 1), p_val = signif(p_val, 3))
print(res, n = Inf, width = Inf)
write.csv(res, file.path(OUT, "results_main.csv"), row.names = FALSE)
# units are kg/ha per 10 SHMI, except "SHMI x age" (per 10 SHMI per year) and
# "SHMI x N" (per 10 SHMI per 100 kg N/ha)
