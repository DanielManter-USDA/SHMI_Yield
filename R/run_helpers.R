################################################################################
# run_helpers.R: helpers for run_all.R (running) and report.Rmd (results)
# Each analysis script runs in a fresh R session (callr) with the project root as the
# working directory, or in this session with in_session = TRUE; python steps run with
# the python found on the PATH. Every step writes its full output to its own log in
# Outputs/logs/<run time>/ and a line in steps.csv; the console shows one line per step.
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
    # low-memory fallback: run in this R session, in a fresh environment so the steps
    # still do not share objects. Output and messages go to the step's log, as with callr.
    env <- new.env(parent = globalenv())
    old <- options(opts); on.exit(options(old), add = TRUE)
    list2env(vars, envir = env)
    con <- file(log, open = "wt")
    sink(con); sink(con, type = "message")
    ok <- tryCatch({ sys.source(path, envir = env, keep.source = FALSE); TRUE },
                   error = function(e) { cat("\nRUN ERROR:", conditionMessage(e), "\n"); FALSE },
                   finally = { sink(type = "message"); sink(); close(con) })
    rm(env); gc(verbose = FALSE)
    log_step(label, ok, t0, log)
    if (!ok) fail(label, log)
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

# The python to use: the one given (setting `python`), else the first that runs and has
# pandas, openpyxl and numpy among RETICULATE_PYTHON or RETICULATE_PYTHON_FALLBACK (set by RStudio's Global Options >
# Python), then python3, python and py on the PATH (on Windows, python3 is often a
# Microsoft Store alias that does not run). Returns NA with the reasons as an attribute.
find_python <- function(py = "") {
  cands <- if (nzchar(py)) py else
    unique(c(Sys.getenv(c("RETICULATE_PYTHON", "RETICULATE_PYTHON_FALLBACK")), unname(Sys.which(c("python3", "python", "py"))),
             # conda installs that RStudio does not see on its PATH
             Sys.getenv("CONDA_PYTHON_EXE"), path.expand(file.path("~", c("miniconda3", "anaconda3", "miniforge3"), "bin", "python"))))
  cands <- cands[nzchar(cands)]
  why <- character()
  
  for (w in cands) {
    if (!file.exists(w) && !nzchar(Sys.which(w))) { 
      why <- c(why, paste0(w, ": not found"))
      next 
    }
    
    # 1. First, check if reticulate is locked up. If it is already bound to this specific path, return it immediately.
    is_already_active <- tryCatch({
      if (reticulate::py_available(initialize = FALSE)) {
        active_path <- reticulate::py_config()$python
        if (identical(normalizePath(active_path, mustWork = FALSE), normalizePath(w, mustWork = FALSE))) return(w)
      }
      FALSE
    }, error = function(e) FALSE)
    
    # 2. To prevent reticulate initialization collision bugs, execute the module check
    # using a pure, primitive system command that uses Python's -c import flag.
    # We pass the text cleanly via single quotes to keep the shell parser stable.
    is_alive <- tryCatch({
      cmd <- paste(shQuote(w), "-c", shQuote("import pandas, openpyxl, numpy"))
      system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE) == 0L
    }, error = function(e) FALSE)
    
    if (is_alive) return(w)
    
    # 3. If primitive evaluation failed, fall back to capturing explicit module errors
    err <- tryCatch({ 
      check_python(w, quiet = TRUE)
      "" 
    }, error = function(e) conditionMessage(e))
    
    if (!nzchar(err)) return(w)
    last <- utils::tail(strsplit(err, "\n")[], 1)
    why <- c(why, paste0(w, ": cannot import pandas, openpyxl and numpy [", last, "]"))
  }
  
  if (!length(cands)) why <- "no python on the PATH and neither RETICULATE_PYTHON nor RETICULATE_PYTHON_FALLBACK set"
  structure(NA_character_, why = why)
}

check_python <- function(py, quiet = FALSE) {
  out <- suppressWarnings(system2(shQuote(py), c("-c", shQuote("import pandas, openpyxl, numpy; print('pandas', pandas.__version__)")),
                                  stdout = TRUE, stderr = TRUE))
  if (!is.null(attr(out, "status")) && attr(out, "status") != 0)
    stop("python needs pandas, openpyxl and numpy: pip install pandas openpyxl numpy\n", paste(out, collapse = "\n"))
  if (!quiet) cat(out, sep = "\n")
  invisible(out)
}

run_py <- function(py, args, label) {
  log <- file.path(LOG_DIR, paste0(label, ".log"))
  t0 <- Sys.time(); message(format(t0, "%H:%M"), "  start  ", label)
  status <- system2(shQuote(py), shQuote(args), stdout = log, stderr = log)
  ok <- identical(as.integer(status), 0L)
  log_step(label, ok, t0, log)
  if (!ok) fail(label, log)
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
  if (isTRUE(p$run_coarse) && !dir.exists("Outputs/Coarse/napeshm")) out <- c(out, "the coarse-record workbooks (B9)")
  out
}

required_inputs <- function(p) {
  conv  <- c("NELITCSE", "NEMLTCRS", "NEMEIRR", "NECCIRR", "INWLACRE", "IAAM7071", "PAUP")
  units <- c("INWLACRE", "IAAM7071", "PAUP")
  tab   <- function(f) file.path(p$agcros_tables, f)
  # the AgCROS tables are needed when a workbook must be built
  need_tables <- conv_needed(p)
  rows <- list(
    # Part A
    c("Inputs/napeshm/mgmtdata.csv", "A0 NAPESHM rebuild", p$run_napeshm_rebuild),
    c("Inputs/napeshm/treatments.csv", "A0 NAPESHM rebuild", p$run_napeshm_rebuild),
    c("Inputs/napeshm/sites.csv", "A0 NAPESHM rebuild", p$run_napeshm_rebuild),
    c("Inputs/napeshm/NAPESHM_test.xlsx", "A0 NAPESHM rebuild (layout)", p$run_napeshm_rebuild),
    c("Inputs/NAPESHM_rebuilt.xlsx", "A1 SHMI inputs for the calibration", p$run_calibration),
    c("Inputs/napeshm/sites.csv", "A1 site climate, descriptors and coordinates", p$run_calibration),
    c("Inputs/napeshm/treatments.csv", "A1 NAPESHM treatments", p$run_calibration),
    c("Inputs/napeshm/soildata.csv", "A1 NAPESHM soil measurements", p$run_calibration),
    c("Inputs/napeshm/plfa_ng.csv", "A1 NAPESHM PLFA", p$run_calibration),
    c("Inputs/SHAPE/batch_Combined_scores.csv", "A6 SHAPE app output (optional)", FALSE),
    # Part B
    c(tab("meassoilchem.csv"), "B0 soil C (AgCROS)", conv_needed(p) || !file.exists("Outputs/yield/soil_c.csv")),
    c(tab("mgtplanting.csv"), "B0 AgCROS conversion", need_tables),
    c(tab("measyieldnutuptake.csv"), "B0 AgCROS conversion", need_tables),
    c("Inputs/yield_workbooks/SHMI_OHTVD_Ho.xlsx", "B2 Hoytville (hand-built from DRIVES)", p$run_yield),
    c("Inputs/yield_workbooks/SHMI_OHTVD_Wo.xlsx", "B2 Wooster (hand-built from DRIVES)", p$run_yield),
    c("Inputs/yield_workbooks/SHMI_WICST_ARS.xlsx", "B2 Arlington (hand-built from DRIVES)", p$run_yield),
    c("Inputs/climate/spei03.nc", "B1 SPEIbase 3-month SPEI (else data/cache/site_spei.csv)", p$run_yield && !file.exists("data/cache/site_spei.csv")),
    c("Inputs/NAPESHM_rebuilt.xlsx", "B9 coarse records (NAPESHM)", isTRUE(p$run_coarse)),
    c("Inputs/napeshm/treatments.csv", "B9 coarse records (irrigation)", isTRUE(p$run_coarse)),
    c("Outputs/Final_plus/final_results.rds", "B9 coarse records (A3 output)", isTRUE(p$run_coarse) && !isTRUE(p$run_calibration)),
    c("Outputs/site_climate_monthly.rds", "B9 coarse records (A1 output)", isTRUE(p$run_coarse) && !isTRUE(p$run_calibration)))
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
  print(d[, c("path", "used_by", "status")], row.names = FALSE)
}

# ---- results -----------------------------------------------------------------------
show_csv <- function(f, caption) {
  if (!file.exists(f)) { cat("\n*", caption, ": not found (", f, ")*\n\n"); return(invisible()) }
  print(knitr::kable(utils::read.csv(f, check.names = FALSE), digits = 3, caption = caption))
  cat("\n")
}

# ---- shipped copies of derived inputs ---------------------------------------------
# Files in data/ that a step would otherwise build from restricted or very large
# sources (AgCROS tables, gSSURGO) or download with a key (NASS). Each is copied into
# Outputs/ only when missing there, so a local run keeps its own files.
SEED_FILES <- c(
  "data/cache/nass_site_county_corn_yield.csv" = "Outputs/yield/nass_site_county_corn_yield.csv",  # NASS Quick Stats (no key needed)
  "data/cache/site_soils_ssurgo.csv"           = "Outputs/site_soils_ssurgo.csv",                  # gSSURGO soil groups
  "data/soil_c.csv"                            = "Outputs/yield/soil_c.csv",                       # AgCROS soil C (derived)
  "data/SHAPE_by_unit_with_SHMI.csv"           = "Outputs/SHAPE/SHAPE_by_unit_with_SHMI.csv")      # SHAPE scores (A6 output)
seed_outputs <- function(files = SEED_FILES) {
  done <- character()
  for (src in names(files)) {
    dst <- files[[src]]
    if (file.exists(src) && !file.exists(dst)) {
      dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
      file.copy(src, dst); done <- c(done, dst)
    }
  }
  if (length(done)) cat("Copied from data/:", paste(done, collapse = ", "), "\n")
  invisible(done)
}
