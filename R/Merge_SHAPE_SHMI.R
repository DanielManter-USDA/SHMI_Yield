################################################################################
# Merge_SHAPE_SHMI.R
#
# Averages the SHAPE app's sample scores by management unit (field_name =
# MGT_combo) and joins them to the final SHMI scores and the measured
# soil-health score (SH).
#
# Inputs
#   Outputs/SHAPE/batch_Combined_scores.csv   SHAPE app output (one row per sample)
#   Outputs/Final_plus/shmi_scores.csv        SHMI, SHMI_base, Amend_bonus, sub-indices
#   Outputs/Final_plus/us_oof_predictions.csv sitecode, system, SH
# Output
#   Outputs/SHAPE/SHAPE_by_unit_with_SHMI.csv  one row per unit
#
# SHAPE_<indicator>   mean over the unit's samples of SHAPE's Score_Mean (0-1)
# SHAPE_<ind>_sd      SD of the sample scores within the unit
# SHAPE_mean          mean of the unit's indicator scores (SOC, AC, Resp,
#                     KSSL_WAS; SHAPE returned no ACE scores)
# <indicator>_raw     mean measured value
# The 2.5% / 97.5% score bounds and the benchmarks are per sample and are not
# averaged (an average of interval bounds is not an interval).
################################################################################
suppressPackageStartupMessages({ library(dplyr); library(here) })
OUT <- function(...) here::here("Outputs", ...)

IND <- c("SOC", "AC", "ACE", "Resp", "KSSL_WAS")
shape <- read.csv(OUT("SHAPE", "batch_Combined_scores.csv"), check.names = FALSE)
scored <- IND[colSums(!is.na(shape[paste0("Score_Mean_", IND)])) > 0]
if (length(setdiff(IND, scored))) message("No SHAPE scores for: ", paste(setdiff(IND, scored), collapse = ", "))

by_unit <- shape %>%
  group_by(MGT_combo = field_name) %>%
  summarise(n_samples = n(), MAPmm = mean(MAPmm), MATC = mean(MATC), Suborder = first(Suborder),
            Texture = paste(sort(unique(na.omit(Texture))), collapse = "/"),
            across(all_of(paste0("Score_Mean_", scored)), list(m = ~ mean(.x, na.rm = TRUE), sd = ~ sd(.x, na.rm = TRUE))),
            across(all_of(IND), ~ mean(.x, na.rm = TRUE), .names = "{.col}_raw"),
            .groups = "drop") %>%
  rename_with(~ sub("^Score_Mean_(.*)_m$", "SHAPE_\\1", .x)) %>%
  rename_with(~ sub("^Score_Mean_(.*)_sd$", "SHAPE_\\1_sd", .x)) %>%
  mutate(SHAPE_mean = rowMeans(across(all_of(paste0("SHAPE_", scored))), na.rm = TRUE),
         SHAPE_n_indicators = rowSums(!is.na(across(all_of(paste0("SHAPE_", scored))))),
         across(where(is.numeric), ~ ifelse(is.nan(.x), NA_real_, .x)))

shmi <- read.csv(OUT("Final_plus", "shmi_scores.csv")) %>%
  select(MGT_combo, SHMI, SHMI_base, Amend_bonus, Cover, Diversity, InvDist, p_amend, p_animal)
sh <- read.csv(OUT("Final_plus", "us_oof_predictions.csv")) %>% select(MGT_combo, sitecode, system, SH)

merged <- by_unit %>% left_join(sh, by = "MGT_combo") %>% left_join(shmi, by = "MGT_combo") %>%
  relocate(MGT_combo, sitecode, system, n_samples, Suborder, Texture, MAPmm, MATC,
           SHMI, SHMI_base, Amend_bonus, SH, SHAPE_mean, SHAPE_n_indicators, starts_with("SHAPE_"))
if (any(is.na(merged$SHMI))) warning(sum(is.na(merged$SHMI)), " units have no SHMI score")
write.csv(merged, OUT("SHAPE", "SHAPE_by_unit_with_SHMI.csv"), row.names = FALSE)

cols <- c("SHMI", "SHMI_base", "SH", "SHAPE_mean", paste0("SHAPE_", scored))
cat("\nCorrelations, all units:\n"); print(round(cor(merged[cols], use = "pairwise")[, c("SHAPE_mean", "SH")], 2))
within <- merged %>% group_by(sitecode) %>% mutate(across(all_of(cols), ~ .x - mean(.x, na.rm = TRUE))) %>% ungroup()
cat("\nCorrelations within sites:\n"); print(round(cor(within[cols], use = "pairwise")[, c("SHAPE_mean", "SH")], 2))
message("Saved ", OUT("SHAPE", "SHAPE_by_unit_with_SHMI.csv"), " (", nrow(merged), " units)")
