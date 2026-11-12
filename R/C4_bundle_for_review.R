################################################################################
# C4_bundle_for_review.R
#
# One zip of every small results file in Outputs/, so the numbers quoted in the
# manuscript can be checked against the latest run without hunting through folders.
#   - every .csv / .txt / .log under Outputs/ up to MAX_MB each (summaries, not plot-level data)
#   - plus Outputs/Manuscript/numbered/ in full (figures and the Supplementary workbook)
# Folder structure is kept inside the zip. Nothing from Inputs/ is included
# (the raw AgCROS tables live there and must not be shared).
# Writes Outputs/review_bundle_<date>.zip and lists what it skipped for size.
################################################################################
MAX_MB <- if (exists("REVIEW_MAX_MB")) REVIEW_MAX_MB else 2
root   <- here::here()
outs   <- here::here("Outputs")
f_all  <- list.files(outs, pattern = "\\.(csv|txt|log)$", recursive = TRUE, full.names = TRUE)
f_all  <- f_all[!grepl("review_bundle_", f_all)]
mb     <- file.size(f_all) / 2^20
keep   <- f_all[mb <= MAX_MB]
num    <- list.files(here::here("Outputs", "Manuscript", "numbered"), recursive = TRUE, full.names = TRUE)
files  <- unique(c(keep, num))
root_n <- normalizePath(root, winslash = "/")
rel    <- substring(normalizePath(files, winslash = "/"), nchar(root_n) + 2)   # paths relative to the repo
zipf   <- here::here("Outputs", paste0("review_bundle_", format(Sys.Date()), ".zip"))
if (file.exists(zipf)) unlink(zipf)
old <- setwd(root); on.exit(setwd(old), add = TRUE)
if (requireNamespace("zip", quietly = TRUE)) {
  zip::zip(zipf, files = rel, mode = "mirror")
} else {
  utils::zip(zipf, files = rel, flags = "-q")   # needs a zip program on the PATH (Rtools on Windows)
}
setwd(old)
skipped <- f_all[mb > MAX_MB]
message("Review bundle: ", zipf, "\n  ", length(rel), " files, ",
        round(file.size(zipf) / 2^20, 1), " MB")
if (length(skipped)) message("  skipped (over ", MAX_MB, " MB): ",
                             paste(sub(paste0(".*Outputs/"), "", skipped), collapse = "; "))
