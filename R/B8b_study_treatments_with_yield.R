################################################################################
# B8b_study_treatments_with_yield.R
#
# The treatment table (B8, Outputs/SHMI_study_treatments.xlsx) with each treatment's
# SHMI and its parts (sub-indices 0-100, amendment bonus, and each part's points
# toward SHMI), workbook treatment codes (MGT_trt) and mean corn and soybean yields added,
# joined on the Key column (= trt in the R outputs: site + MGT_trt_base).
#   yield_mean, yield_sd   plain mean corn yield (kg/ha dry matter) over its corn plot-years
#   yield_vs_siteyr        mean difference from the site-year mean of all corn treatments
#                          (removes year: compare treatments within an experiment with this)
#   soy_*                  the same for soybean (B2b's six experiments)
# Neither yield adjusts for N or the other treatment factors (that is the model's job).
#
# Input   Outputs/SHMI_study_treatments.xlsx (B8); Outputs/yield_main_default_total_agcros/
#         B2_data.rds (B2); data/treatment_labels.csv (from B2f, retired);
#         Outputs/yield_soybean/soy_treatment_yields.csv (B2b)
# Output  Outputs/SHMI_study_treatments_with_yield.xlsx
################################################################################
suppressPackageStartupMessages(library(dplyr))
setwd(here::here())
XL  <- "Outputs/SHMI_study_treatments.xlsx"
OUT <- "Outputs/SHMI_study_treatments_with_yield.xlsx"
for (f in c(XL, "Outputs/yield_main_default_total_agcros/B2_data.rds"))
  if (!file.exists(f)) stop(f, " not found (run B8 / B2 first)")

xl  <- readxl::read_excel(XL, sheet = "Treatments")
key <- grep("^Key", names(xl), value = TRUE)
if (!length(key)) stop(XL, " has no Key column: rebuild it with the updated python/build_study_treatments.py")
xl  <- rename(xl, trt = all_of(key))

corn <- readRDS("Outputs/yield_main_default_total_agcros/B2_data.rds")$corn
corn_trt <- corn %>% group_by(site, year) %>% mutate(dev = yield_kgha - mean(yield_kgha)) %>% ungroup() %>%
  mutate(trt = as.character(trt)) %>% group_by(trt) %>%
  summarise(MGT_trt = paste(sort(unique(as.character(MGT_trt))), collapse = ", "),
            SHMI = round(mean(SHMI[!duplicated(MGT_combo)]), 1),
            across(any_of(c("Cover", "Diversity", "InvDist", "SHMI_base", "Amend_bonus")),
                   ~ round(mean(.x[!duplicated(MGT_combo)]), 1)),
            corn_years = n_distinct(year), corn_plot_years_model = n(),
            yield_mean = round(mean(yield_kgha)), yield_sd = round(sd(yield_kgha)),
            yield_vs_siteyr = round(mean(dev)), .groups = "drop")

# Treatment labels, median recorded N and other factors were written by B2f (retired; not in the
# paper). Its output is shipped as data/treatment_labels.csv, so this table does not need B2f.
lab_f <- c("Outputs/yield_paired_budget/treatment_labels.csv", "data/treatment_labels.csv")
lab_f <- lab_f[file.exists(lab_f)][1]
if (is.na(lab_f)) warning("data/treatment_labels.csv not found: label, N_rec_median and other_factor columns left out")
lab <- if (!is.na(lab_f)) read.csv(lab_f, stringsAsFactors = FALSE) %>%
  transmute(trt = as.character(trt), N_rec_median = N_rec, other_factor = other,
            label_cover = cover, label_tillage = tillage, label_rotation = rotation, label_manure = manure) else NULL
soy_f <- "Outputs/yield_soybean/soy_treatment_yields.csv"
soy <- if (file.exists(soy_f)) read.csv(soy_f, stringsAsFactors = FALSE) %>% select(-site) else NULL

out <- xl %>% left_join(corn_trt, by = "trt")
if (!is.null(lab)) out <- left_join(out, lab, by = "trt")
if (!is.null(soy)) out <- out %>% left_join(soy, by = "trt") %>%
  mutate(MGT_trt = coalesce(MGT_trt, soy_MGT_trt), SHMI = coalesce(SHMI, soy_SHMI)) %>% select(-soy_MGT_trt, -soy_SHMI)
# each part's contribution to SHMI (points): weight x sub-index, and the amendment bonus;
# the weights are the published SHMI 2.0.0 formula weights (B2 respects overrides; these
# columns use the defaults)
W <- c(Cover = 0.687, Diversity = 0.150, InvDist = 0.163)
if (all(names(W) %in% names(out)))
  out <- out %>% mutate(pts_Cover = round(W[["Cover"]] * Cover, 1), pts_Diversity = round(W[["Diversity"]] * Diversity, 1),
                        pts_InvDist = round(W[["InvDist"]] * InvDist, 1))
out <- out %>% relocate(any_of(c("MGT_trt", "SHMI", "Cover", "Diversity", "InvDist", "SHMI_base", "Amend_bonus",
                                 "pts_Cover", "pts_Diversity", "pts_InvDist")), .after = trt) %>% rename(!!key := trt)

miss_corn <- sum(out$`Corn plot-years with yield` > 0 & is.na(out$yield_mean))
cat("Treatments:", nrow(out), "| with corn yields:", sum(!is.na(out$yield_mean)),
    "| with soybean yields:", if (is.null(soy)) 0 else sum(!is.na(out$soy_yield_mean)),
    "| corn in the table but no R match:", miss_corn, "\n")
if (miss_corn) print(out %>% filter(`Corn plot-years with yield` > 0, is.na(yield_mean)) %>% select(Site, all_of(key)))
sheets <- list(Treatments = out)
for (sh in setdiff(readxl::excel_sheets(XL), "Treatments")) sheets[[sh]] <- readxl::read_excel(XL, sheet = sh)
writexl::write_xlsx(sheets, OUT)
cat("Written:", OUT, "\n")
