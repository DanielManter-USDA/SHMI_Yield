################################################################################
# datasets.R
#
# The analysis datasets, defined once and sourced by every pipeline script
# (SEM, weight tests, calibration, summary). Each dataset is taken through the
# whole chain on its own: indicator adjustment -> SEM -> weight tests ->
# SHMI calibration. Results are then compared between datasets.
#
#   all   all units, after the exclusions in Prepare_SEM.R (woody plantations
#         and grazing lands); reference soils included
#   noMX  all, without the Mexican sites
#
# A filter works on any table with `sitecode` and `system` columns
# (semDF_prepped.rds, the SEM score files).
################################################################################

is_mx        <- function(d) grepl("^MX", as.character(d$sitecode))
is_perennial <- function(d) !is.na(d$system) & d$system == "perennial"

DATASETS <- list(
  all  = list(label = "All units (woody and grazing-land treatments excluded in Prepare_SEM.R)",
              keep  = function(d) rep(TRUE, nrow(d))),
  noMX = list(label = "Without Mexico",
              keep  = function(d) !is_mx(d))
)

# Dataset whose calibration is reported as the main result (can be changed)
PRIMARY_DATASET <- "all"

# Candidates for the primary dataset: each also gets the GAM-adjustment
# sensitivity run and its own transfer figure
CANDIDATE_DATASETS <- c("all", "noMX")

ds_filter <- function(d, key) d[DATASETS[[key]]$keep(d), , drop = FALSE]
ds_label  <- function(key) vapply(key, function(k) DATASETS[[k]]$label, character(1))
