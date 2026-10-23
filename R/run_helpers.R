################################################################################
# run_helpers.R: helpers for run_all.Rmd
# Each analysis script runs in a fresh R session (callr) with the project root as the
# working directory; python steps run with the python found on the PATH. Every step
# writes its own log in Outputs/logs/<run time>/ and a line in steps.csv.
################################################################################

start_log <- function(root) {
  stamp <- format(Sys.time(), "%Y%m%d_%H%M")
  LOG_DIR <<- file.path(root, "Outputs", "logs", stamp)
  dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
  write.csv(data.frame(step = character(), ok = logical(), minutes = numeric(), log = character()),
            file.path(LOG_DIR, "steps.csv"), row.names = FALSE)
  invisible(LOG_DIR)
}

log_step <- function(label, ok, t0, log) {
  row <- data.frame(step = label, ok = ok, minutes = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
                    log = basename(log))
  write.table(row, file.path(LOG_DIR, "steps.csv"), sep = ",", append = TRUE, col.names = FALSE, row.names = FALSE)
  message(format(Sys.time(), "%H:%M"), "  ", if (ok) "done   " else "FAILED ", label, " (", row$minutes, " min)")
}

fail <- function(label, log) {
  tail_txt <- tryCatch(utils::tail(readLines(log, warn = FALSE), 25), error = function(e) "")
  stop("Step '", label, "' failed. Last lines of ", log, ":\n", paste(tail_txt, collapse = "\n"), call. = FALSE)
}

# run an R script from R/ in a fresh session; `vars` are set in its global environment
# first (e.g. OUT_OVERRIDE), `opts` are R options (e.g. calibration iterations)
run_r <- function(script, vars = list(), opts = list(), label = sub("\\.R$", "", script)) {
  path <- here::here("R", script)
  if (!file.exists(path)) stop("Script not found: ", path)
  log <- file.path(LOG_DIR, paste0(label, ".log"))
  t0 <- Sys.time(); message(format(t0, "%H:%M"), "  start  ", label)
  if (isTRUE(getOption("shmi.in_session", FALSE))) {
    # debugging / low-memory fallback: run in this R session (output to the console),
    # in a fresh environment so the steps still do not share objects
    message("  (in this session; log not written)")
    env <- new.env(parent = globalenv())
    old <- options(opts); on.exit(options(old), add = TRUE)
    list2env(vars, envir = env)
    ok <- tryCatch({ sys.source(path, envir = env, keep.source = FALSE); TRUE },
                   error = function(e) { message("ERROR in ", label, ": ", conditionMessage(e)); FALSE })
    rm(env); gc(verbose = FALSE)
    log_step(label, ok, t0, log)
    if (!ok) stop("Step '", label, "' failed (see the console above).", call. = FALSE)
    return(invisible(TRUE))
  }
  ok <- tryCatch({
    callr::r(function(root, path, vars, opts) {
      setwd(root)
      options(opts)
      list2env(vars, envir = globalenv())
      source(path, echo = TRUE, max.deparse.length = Inf, spaced = FALSE, encoding = "UTF-8")
      invisible(TRUE)
    }, args = list(root = here::here(), path = path, vars = vars, opts = opts),
    stdout = log, stderr = "2>&1", show = FALSE)
    TRUE
  }, error = function(e) { cat("\nRUN ERROR:", conditionMessage(e), "\n", file = log, append = TRUE); FALSE })
  log_step(label, ok, t0, log)
  if (!ok) fail(label, log)
  invisible(TRUE)
}

# the first python on the PATH that actually runs (on Windows, python3 is often a
# Microsoft Store alias that does not)
find_python <- function(py = "") {
  if (nzchar(py)) return(py)
  for (p in c("python3", "python", "py")) {
    w <- unname(Sys.which(p))
    if (nzchar(w) && identical(suppressWarnings(system2(shQuote(w), "--version", stdout = FALSE, stderr = FALSE)), 0L))
      return(w)
  }
  stop("No working python on the PATH; set params$python to the python executable.")
}

check_python <- function(py) {
  out <- suppressWarnings(system2(shQuote(py), c("-c", shQuote("import pandas, openpyxl, numpy; print('pandas', pandas.__version__)")),
                                  stdout = TRUE, stderr = TRUE))
  if (!is.null(attr(out, "status")) && attr(out, "status") != 0)
    stop("python needs pandas, openpyxl and numpy: pip install pandas openpyxl numpy\n", paste(out, collapse = "\n"))
  cat(out, sep = "\n")
}

run_py <- function(py, args, label) {
  log <- file.path(LOG_DIR, paste0(label, ".log"))
  t0 <- Sys.time(); message(format(t0, "%H:%M"), "  start  ", label)
  status <- system2(shQuote(py), shQuote(args), stdout = log, stderr = log)
  ok <- identical(as.integer(status), 0L)
  log_step(label, ok, t0, log)
  if (!ok) fail(label, log)
  cat(utils::tail(readLines(log, warn = FALSE), 3), sep = "\n")
  invisible(TRUE)
}

# ---- inputs ------------------------------------------------------------------------
CONVERTED <- c("NELITCSE", "NEMLTCRS", "NEMEIRR", "NECCIRR", "INWLACRE", "IAAM7071", "PAUP")
# the AgCROS workbooks are converted when asked (run_conversion) or when any is missing;
# otherwise the existing ones are kept, so their site runs stay cached
conv_needed <- function(p, sites = CONVERTED)
  isTRUE(p$run_conversion) || !all(file.exists(sprintf("Inputs/yield_workbooks/SHMI_%s.xlsx", sites)))
# files that only python can build, and that this run needs
python_needed <- function(p) {
  out <- character()
  if (isTRUE(p$run_napeshm_rebuild)) out <- c(out, "the NAPESHM workbook (A0)")
  if (conv_needed(p)) out <- c(out, "the AgCROS workbooks in Inputs/yield_workbooks (B0)")
  if (!file.exists("Outputs/yield/soil_c.csv")) out <- c(out, "Outputs/yield/soil_c.csv (B0)")
  if (isTRUE(p$run_site_sensitivity) && conv_needed(p, c("INWLTPAC", "MNMOFS")))
    out <- c(out, "SHMI_INWLTPAC.xlsx and SHMI_MNMOFS.xlsx for B9")
  out
}

required_inputs <- function(p) {
  conv  <- c("NELITCSE", "NEMLTCRS", "NEMEIRR", "NECCIRR", "INWLACRE", "IAAM7071", "PAUP")
  units <- c("INWLACRE", "IAAM7071", "PAUP")
  tab   <- function(f) file.path(p$agcros_tables, f)
  # the AgCROS tables are needed when a workbook must be built; INWLTPAC and MNMOFS
  # (B9 only) are converted whenever their workbooks are missing
  need_tables <- conv_needed(p) || (isTRUE(p$run_site_sensitivity) && conv_needed(p, c("INWLTPAC", "MNMOFS")))
  rows <- list(
    # Part A
    c("Inputs/napeshm_raw/mgmtdata.csv", "A0 NAPESHM rebuild", p$run_napeshm_rebuild),
    c("Inputs/napeshm_raw/treatments.csv", "A0 NAPESHM rebuild", p$run_napeshm_rebuild),
    c("Inputs/napeshm_raw/sites.csv", "A0 NAPESHM rebuild", p$run_napeshm_rebuild),
    c("Inputs/napeshm_raw/NAPESHM_test.xlsx", "A0 NAPESHM rebuild (layout)", p$run_napeshm_rebuild),
    c("Inputs/NAPESHM_rebuilt.xlsx", "A1 SHMI inputs for the calibration", p$run_calibration),
    c("Inputs/sites.csv", "A1 site climate and descriptors", p$run_calibration),
    c("Inputs/CSV Files/sites.csv", "A1 site coordinates (soils)", p$run_calibration),
    c("Inputs/CSV Files/treatments.csv", "A1 NAPESHM treatments", p$run_calibration),
    c("Inputs/CSV Files/soildata.csv", "A1 NAPESHM soil measurements", p$run_calibration),
    c("Inputs/CSV Files/plfa_ng.csv", "A1 NAPESHM PLFA", p$run_calibration),
    c("Inputs/SHAPE/batch_Combined_scores.csv", "A6 SHAPE app output (optional)", FALSE),
    # Part B
    c(tab("meassoilchem.csv"), "B0 soil C (AgCROS)", conv_needed(p) || !file.exists("Outputs/yield/soil_c.csv")),
    c(tab("mgtplanting.csv"), "B0 AgCROS conversion", need_tables),
    c(tab("measyieldnutuptake.csv"), "B0 AgCROS conversion", need_tables),
    c("Inputs/yield_workbooks/SHMI_OHTVD_Ho.xlsx", "B2 Hoytville (hand-built from DRIVES)", p$run_yield),
    c("Inputs/yield_workbooks/SHMI_OHTVD_Wo.xlsx", "B2 Wooster (hand-built from DRIVES)", p$run_yield),
    c("Inputs/yield_workbooks/SHMI_WICST_ARS.xlsx", "B2 Arlington (hand-built from DRIVES)", p$run_yield),
    c("Inputs/climate/spei03.nc", "B1 SPEIbase 3-month SPEI", p$run_yield || p$run_site_sensitivity))
  rows <- c(rows, list(c(tab("metaunits.csv"), "B0 AgCROS conversion", need_tables)))
  if (!conv_needed(p)) {
    rows <- c(rows, lapply(conv, function(s) c(sprintf("Inputs/yield_workbooks/SHMI_%s.xlsx", s), "B2 converted workbook", p$run_yield)),
              lapply(units, function(s) c(sprintf("Inputs/yield_workbooks/SHMI_%s_units.csv", s), "B2 AgCROS treatment descriptors", p$run_yield)))
  }
  d <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(d) <- c("path", "used_by", "required")
  d$required <- as.logical(d$required)
  d$present <- file.exists(d$path)
  nass <- data.frame(path = "NASS_API_KEY (~/.Renviron)", used_by = "B2 county yields (NASS Quick Stats)",
                     required = p$run_yield && !file.exists("Outputs/yield/nass_site_county_corn_yield.csv"),
                     present = nzchar(Sys.getenv("NASS_API_KEY")))
  rbind(d, nass)
}

show_inputs <- function(d) {
  d$status <- ifelse(d$present, "ok", ifelse(d$required, "MISSING", "not needed for this run"))
  print(knitr::kable(d[, c("path", "used_by", "status")]))
}

# ---- results -----------------------------------------------------------------------
show_csv <- function(f, caption) {
  if (!file.exists(f)) { cat("\n*", caption, ": not found (", f, ")*\n\n"); return(invisible()) }
  print(knitr::kable(utils::read.csv(f, check.names = FALSE), digits = 3, caption = caption))
  cat("\n")
}

# within-site SHMI slope (Fig. 3a model) and the pooled site slope for every main run
site_combination_table <- function() {
  dirs <- list.dirs("Outputs", recursive = FALSE)
  dirs <- dirs[grepl("^yield_main_", basename(dirs)) & file.exists(file.path(dirs, "results_main.csv"))]
  do.call(rbind, lapply(dirs, function(d) {
    r <- utils::read.csv(file.path(d, "results_main.csv"))
    fig1 <- r[r$figure == "Fig 1", ][1, ]
    pool <- r[r$estimate == "Pooled across sites (random effects, Knapp-Hartung)", ][1, ]
    sxn <- r[grepl("^SHMI x N: change", r$estimate), ][1, ]
    data.frame(run = sub("^yield_main_", "", basename(d)), slope_per10 = fig1$effect, lower = fig1$lower,
               upper = fig1$upper, pooled = pool$effect, pooled_lower = pool$lower, pooled_upper = pool$upper,
               sites = pool$n, SHMIxN = sxn$effect, SHMIxN_p = sxn$p_val, N_sites = sxn$n)
  }))
}
