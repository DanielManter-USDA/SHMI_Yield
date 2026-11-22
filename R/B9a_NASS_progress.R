################################################################################
# B9a_NASS_progress.R  (17 Nov 2026)
#
# Typical planting and harvest dates from USDA-NASS Crop Progress, for the coarse-record
# scenarios that give SHMI only the dates public statistics provide (B9, operation
# nassdates in degrade_records.py). For each state (and the US), crop, stage and year:
# the day of year on which the state reached 50% planted or harvested, interpolated
# linearly between the two weekly reports either side of 50%.
#
# Needs NASS_API_KEY in ~/.Renviron (never in the script). Downloads once; B9 calls this
# script only when the output is missing (delete it to download again).
#
# Output  data/cache/nass_progress_50pct.csv  (state, commodity, stage, year, doy50, weeks)
#         data/cache/nass_progress_raw.csv    (the weekly reports used)
################################################################################
suppressPackageStartupMessages({ library(dplyr); library(rnassqs) })
CACHE <- here::here("data", "cache"); dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)
key <- Sys.getenv("NASS_API_KEY")
if (!nzchar(key)) stop("NASS_API_KEY is not set: add NASS_API_KEY=... to ~/.Renviron and restart R", call. = FALSE)
nassqs_auth(key = key)

# commodities NASS reports progress for that occur in the workbooks (degrade_records.py
# maps crop names to these; wheat is split by class below)
COMM <- c("CORN", "SOYBEANS", "SORGHUM", "OATS", "BARLEY", "WHEAT", "COTTON", "SUNFLOWER",
          "PEANUTS", "RICE", "SUGARBEETS", "BEANS")
STAGE <- c(planted = "PCT PLANTED", harvested = "PCT HARVESTED")

get1 <- function(comm, unit, agg) {
  q <- list(source_desc = "SURVEY", sector_desc = "CROPS", commodity_desc = comm,
            statisticcat_desc = "PROGRESS", unit_desc = unit, agg_level_desc = agg)
  tryCatch(nassqs(q), error = function(e) { message("  none: ", comm, ", ", unit, ", ", agg); NULL })
}
raw <- bind_rows(lapply(COMM, function(cm) bind_rows(lapply(names(STAGE), function(st) {
  message("NASS Crop Progress: ", cm, ", ", STAGE[[st]])
  bind_rows(lapply(c("STATE", "NATIONAL"), function(agg) {
    d <- get1(cm, STAGE[[st]], agg)
    if (is.null(d) || !nrow(d)) return(NULL)
    d %>% transmute(state = if (agg == "NATIONAL") "US" else state_alpha, commodity_desc, class_desc,   # if, not ifelse: agg is one value
                    util_practice_desc, short_desc, stage = st, week_ending = as.Date(week_ending),
                    value = suppressWarnings(as.numeric(gsub(",", "", Value))))
  }))
}))))
if (!nrow(raw)) stop("No NASS Crop Progress records returned", call. = FALSE)
if (n_distinct(raw$state) < 5) stop("Only ", paste(unique(raw$state), collapse = ", "),
                                    " returned: the state queries failed", call. = FALSE)

# one series per crop as degrade_records.py names them: silage dropped; wheat by class;
# dry edible beans only
raw <- raw %>%
  filter(!grepl("SILAGE", util_practice_desc), !is.na(value), !is.na(week_ending)) %>%
  mutate(commodity = case_when(
    commodity_desc == "WHEAT" & grepl("WINTER", class_desc) ~ "WHEAT, WINTER",
    commodity_desc == "WHEAT" & grepl("SPRING", class_desc) ~ "WHEAT, SPRING",
    commodity_desc == "WHEAT" ~ NA_character_,
    commodity_desc == "BEANS" & grepl("DRY EDIBLE", class_desc) ~ "BEANS, DRY EDIBLE",
    commodity_desc == "BEANS" ~ NA_character_,
    commodity_desc == "COTTON" & !grepl("PIMA", class_desc) ~ "COTTON",
    commodity_desc == "COTTON" ~ NA_character_,
    TRUE ~ commodity_desc)) %>%
  filter(!is.na(commodity)) %>%
  mutate(year = as.integer(format(week_ending, "%Y"))) %>%
  distinct(state, commodity, stage, year, week_ending, .keep_all = TRUE)
write.csv(raw, file.path(CACHE, "nass_progress_raw.csv"), row.names = FALSE)

# 50% date: linear interpolation between the last report below 50% and the first at or
# above it (the calendar year of the reports, so winter wheat planting falls in the autumn
# of its sowing year, as in the workbooks); a series that starts at or above 50% takes its
# first report; one that never reaches 50% is dropped
doy50 <- raw %>% group_by(state, commodity, stage, year) %>% arrange(week_ending, .by_group = TRUE) %>%
  summarise(doy50 = {
    i <- which(value >= 50)[1]
    if (is.na(i)) NA_real_ else if (i == 1) as.numeric(format(week_ending[1], "%j")) else {
      d0 <- as.numeric(week_ending[i - 1]); d1 <- as.numeric(week_ending[i])
      d <- d0 + (50 - value[i - 1]) / (value[i] - value[i - 1]) * (d1 - d0)
      as.numeric(format(as.Date(d, origin = "1970-01-01"), "%j"))
    }
  }, weeks = n(), .groups = "drop") %>% filter(!is.na(doy50))
write.csv(doy50, file.path(CACHE, "nass_progress_50pct.csv"), row.names = FALSE)

message(sprintf("NASS 50%% dates: %d state-crop-stage-years, %d states (+US), %d-%d",
                nrow(doy50), n_distinct(doy50$state[doy50$state != "US"]), min(doy50$year), max(doy50$year)))
print(doy50 %>% group_by(commodity, stage) %>%
        summarise(states = n_distinct(state[state != "US"]), years = paste(range(year), collapse = "-"),
                  median_doy = median(doy50), .groups = "drop"), n = Inf)
