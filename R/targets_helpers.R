################################################################################
# targets_helpers.R  (11 Nov 2026)
# Small helpers used by _targets.R. Each pipeline step runs one existing script (through
# run_r() in run_helpers.R, with its log in Outputs/logs/<run>/) and returns the files it
# wrote that later steps read. targets tracks those files: a step reruns only when its
# script, an input file or an upstream step's files change.
################################################################################
source(here::here("R", "run_helpers.R"))

# Run one script; the arguments after `outputs` are upstream targets, listed only so that
# targets knows the order (their values are not used). Returns the output files, which
# must all exist afterwards.
start_run <- function() {
  if (exists("LOG_DIR", envir = globalenv())) return(invisible())
  start_log(here::here())
  seed_outputs()     # shipped copies of derived inputs (data/), as run_all.R does
}
# `needs`: files the script cannot run without that no step here writes (the frozen Part A
# outputs). If one is missing, the step keeps its existing outputs when it has them (as
# run_all.R skips the Part A redraws) and otherwise stops, naming what to restore.
# SHMI_KEEP: steps (script prefixes, e.g. "B5b" or "B4,B5b") whose existing outputs are
# accepted without running them, e.g. to register work that finished when tar_make() itself
# stopped: Sys.setenv(SHMI_KEEP = "B5b"); targets::tar_make(); Sys.unsetenv("SHMI_KEEP")
kept <- function(script, outputs) {
  key <- sub("_.*$", "", script)
  keep <- trimws(strsplit(Sys.getenv("SHMI_KEEP"), ",")[[1]])
  if (key %in% keep && all(file.exists(outputs))) {
    message(script, ": existing outputs accepted (SHMI_KEEP)")
    return(TRUE)
  }
  FALSE
}
step <- function(script, outputs, ..., vars = list(), needs = character()) {
  if (kept(script, outputs)) return(outputs)
  start_run()
  miss <- needs[!file.exists(needs)]
  if (length(miss)) {
    if (all(file.exists(outputs))) {
      message(script, ": not rerun (missing ", paste(miss, collapse = ", "), "); its existing outputs are kept")
      return(outputs)
    }
    stop(script, " needs ", paste(miss, collapse = ", "),
         ". These are frozen Part A outputs: copy them back into Outputs/ (e.g. from Outputs_ref).")
  }
  run_r(script, vars = vars)
  check_outputs(script, outputs)
}

# Like step(), but for very long steps (B9, hours): reruns only when forced with
# Sys.setenv(SHMI_FORCE = "B9") (or a comma-separated list); otherwise the existing
# output files are kept and tracked.
step_slow <- function(script, outputs, ..., key, vars = list(), needs = character()) {
  forced <- key %in% trimws(strsplit(Sys.getenv("SHMI_FORCE"), ",")[[1]])
  if (forced || !all(file.exists(outputs))) return(step(script, outputs, vars = vars, needs = needs))
  message(script, ": kept the existing outputs (Sys.setenv(SHMI_FORCE = \"", key, "\") to rerun)")
  outputs
}

# A python step (same rules as step())
step_py <- function(args, label, outputs, ...) {
  start_run()
  py <- as.character(find_python(""))
  if (is.na(py)) stop("Python is needed for ", label, " (pandas, openpyxl, numpy); see Inputs/README.md")
  run_py(py, args, label)
  check_outputs(label, outputs)
}

# A figure script from R/figures/ (seconds; no analysis rerun)
figure <- function(script, outputs, ...) {
  sys.source(here::here("R", "figures", script), envir = new.env(parent = globalenv()))
  check_outputs(script, outputs)
}

check_outputs <- function(what, outputs) {
  miss <- outputs[!file.exists(outputs)]
  if (length(miss)) stop(what, " did not write: ", paste(miss, collapse = ", "))
  outputs
}

# Input files that may be absent on some machines (frozen Part A outputs, caches): the ones
# present are tracked, a missing one is reported once.
present <- function(...) {
  f <- c(...)
  if (any(!file.exists(f))) message("Not found (not tracked): ", paste(f[!file.exists(f)], collapse = ", "))
  f[file.exists(f)]
}
