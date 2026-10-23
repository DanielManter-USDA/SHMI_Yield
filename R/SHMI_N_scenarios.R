################################################################################
# SHMI_N_scenarios.R
# Is WICST's organic yield penalty really N supply? (legume credits and manure N)
#
# Uses the data of SHMI_yield_main.R (its sections 1-3). Only WICST is refitted.
#
# For every combination of
#   legume credit scenario   none / low / default / high   (SHMI_legume_credit.R)
#   manure N scenario        none / low / default / high / total  (SHMI_manure_n.R)
#   organic credit factor    1, 0.5, 0: share of the legume credit realized in the
#                            organic systems (CS3, CS5) only
# N for each corn crop is recomputed (fertilizer + available manure N + legume
# credit) and WICST's own within-site model is fitted twice by maximum likelihood:
#   with    the organic factor:  yield ~ site x year + N + N^2 + organic + SHMI
#   without the organic factor:  yield ~ site x year + N + N^2 + SHMI
#   random plot, AR(1) over years (as fit_site in SHMI_yield_main.R)
# The yields are the same in every fit; only N changes. So:
#   AIC_with (lower = better): which N accounting explains WICST's yields best
#   organic_needed = AIC_without - AIC_with: how much the organic factor still
#     improves the fit once N is counted that way (> ~4: clearly needed; near 0
#     or negative: N supply alone accounts for the organic systems' yields)
#   organic_effect: organic minus conventional at the same SHMI and N (kg/ha)
#   SHMI_slope: WICST's within-site slope (kg/ha per 10 SHMI), with the organic factor
# If organic_needed stays large in every scenario, the penalty is not an N-accounting
# problem (weeds, timing); if it falls to ~0 when the organic credit is discounted,
# the organic systems are better described as getting less N than credited.
# Caveat: lower N for the organic plots and the organic factor explain the same gap,
# so a shrinking penalty shows the gap is consistent with N supply, not that N is
# the cause.
#
# Output: <OUT>/N_scenarios_WICST.csv (OUT = the main script's output folder)
#
# PART 2: ALL RESEARCH SITES (the Fig 1 model; ten sites with INCLUDE_AGCROS)
# For every legume x manure combination, N is recomputed for every plot-year with the
# rules of SHMI_yield_main.R (a missing record is 0 N at NEMLTCRS, WICST and NECCIRR,
# the optimum at the Ohio sites; at the AgCROS sites the median recorded N of the
# plot's treatment, 0 where unfertilized by design; else the site-year median; no
# legume credit where the optimum is assumed) and the Fig 1 model is refitted by
# maximum likelihood:
#   yield ~ site x year + linear N per site (where N varies within site-years; quadratic
#           at N_QUAD_SITES) + other + SHMI_w + SHMI_w:N_c (N_INTERACT_SITES) + SHMI_w:age_c
#           + (1 | treatment/plot), AR(1)
# on exactly the same plot-years (those of the default run in memory), so the AICs
# are comparable: dAIC_best = 0 is the N accounting that explains the yields best,
# within ~2 equally good. Also reported: the SHMI slope at 150 kg N/ha and average age
# (Fig 1), SHMI x N, and the WICST organic effect.
# About 20 fits of the full model: allow 30-60 minutes. Results are written after
# each fit, so an interrupted run keeps what it has.
# Output: <OUT>/N_scenarios_all_sites.csv
################################################################################
MAIN_FILE <- if (exists("MAIN_FILE")) MAIN_FILE else here::here("R", "SHMI_yield_main.R")

# Always reloaded (sections 1-3 of SHMI_yield_main.R), so the sites and settings match
# the main script rather than whatever `corn` is left in memory; the site runs come
# from the cache.
local({
  txt <- readLines(MAIN_FILE, encoding = "UTF-8")
  stop_at <- grep("^# 4\\. FIGURE HELPERS", txt)[1] - 2
  eval(parse(text = txt[1:stop_at], encoding = "UTF-8"), envir = globalenv())
})
message("N scenarios: ", n_distinct(corn$site), " sites; output folder ", OUT)

wi0 <- corn %>% filter(site == "WICST_ARS") %>% droplevels()
need <- c(paste0("N_rec_", c("none", "low", "default", "high", "total")), "legume_class", "other", "SHMI_w")
miss <- setdiff(need, names(wi0))
if (length(miss)) stop("Not in corn: ", paste(miss, collapse = ", "))
wi0 <- wi0 %>% mutate(organic = as.numeric(other == "WICST_organic"))
cat("WICST:", nrow(wi0), "plot-years;", n_distinct(wi0$MGT_trt_base), "treatments;",
    sum(wi0$organic), "organic plot-years\n")

fit_wi <- function(d, organic = TRUE) {
  rhs <- c("site_year", "N100", "I(N100^2)", if (organic) "organic", "SHMI_w")
  lme(reformulate(rhs, "yield_kgha"), random = ~ 1 | trt/plot_id,
      correlation = corCAR1(form = ~ year | trt/plot_id), data = d, method = "ML", control = ctrl)
}

grid <- tidyr::expand_grid(legume = c("none", "low", "default", "high"),
                           manure = c("none", "low", "default", "high", "total"),
                           org_credit = c(1, 0.5, 0))
res <- purrr::pmap_dfr(grid, function(legume, manure, org_credit) {
  d <- wi0 %>%
    mutate(N_leg = legume_credit_kg(legume_class, legume) * ifelse(organic == 1, org_credit, 1),
           N100  = (coalesce(.data[[paste0("N_rec_", manure)]], 0) + N_leg) / 100)   # WICST: missing record = 0
  m1 <- tryCatch(fit_wi(d, TRUE),  error = function(e) NULL)
  m0 <- tryCatch(fit_wi(d, FALSE), error = function(e) NULL)
  if (is.null(m1) || is.null(m0)) return(tibble(legume, manure, org_credit, note = "fit failed"))
  tt <- summary(m1)$tTable
  tibble(legume, manure, org_credit,
         N_conv_mean = 100 * mean(d$N100[d$organic == 0]), N_org_mean = 100 * mean(d$N100[d$organic == 1]),
         AIC_with = AIC(m1), AIC_without = AIC(m0), organic_needed = AIC(m0) - AIC(m1),
         organic_effect = tt["organic", "Value"], organic_se = tt["organic", "Std.Error"],
         SHMI_slope = 10 * tt["SHMI_w", "Value"], SHMI_se = 10 * tt["SHMI_w", "Std.Error"])
}) %>%
  mutate(dAIC_best = AIC_with - min(AIC_with, na.rm = TRUE)) %>%    # 0 = best-fitting N accounting
  arrange(AIC_with)

write.csv(res %>% mutate(across(where(is.numeric), ~ round(.x, 1))),
          file.path(OUT, "N_scenarios_WICST.csv"), row.names = FALSE)
print(res, n = Inf, width = Inf)
cat("\nDefault accounting (default / default / 1):\n")
print(filter(res, legume == "default", manure == "default", org_credit == 1), width = Inf)


################################################################################
# PART 2: all research sites, the Fig 1 model
################################################################################
RUN_ALL_SITES <- TRUE
if (RUN_ALL_SITES) {
  if (!(LEGUME_SCENARIO == "default" && MANURE_SCENARIO %in% c("default", "total")))
    warning("corn in memory comes from a ", LEGUME_SCENARIO, "/", MANURE_SCENARIO,
            " run; part 2 is meant to start from the default run's plot-years")
  c0 <- corn %>% select(-matches("^N[12]_"))
  all_file <- file.path(OUT, "N_scenarios_all_sites.csv")
  cat("Part 2:", n_distinct(c0$site), "sites; quadratic N at", paste(N_QUAD_SITES, collapse = ", "),
      "; SHMI x N at", paste(N_INTERACT_SITES, collapse = ", "), "\n")
  grid2 <- tidyr::expand_grid(legume = c("none", "low", "default", "high"),
                              manure = c("none", "low", "default", "high", "total"))
  res2 <- list()
  for (i in seq_len(nrow(grid2))) {
    lg <- grid2$legume[i]; mn <- grid2$manure[i]
    message(sprintf("[%d/%d] legume = %s, manure = %s", i, nrow(grid2), lg, mn))
    d <- c0 %>% select(-any_of("N_trt")) %>%
      mutate(N_rec = .data[[paste0("N_rec_", mn)]])
    # AgCROS sites: a missing record takes the treatment's median (0 if unfertilized by design)
    trt_n <- d %>% filter(site %in% AGCROS_SITES) %>% group_by(site, MGT_trt_base) %>%
      summarise(N_trt = ifelse(grepl(UNFERTILIZED, first(MGT_trt_base), ignore.case = TRUE) | all(is.na(N_rec)),
                               0, median(N_rec, na.rm = TRUE)), .groups = "drop")
    d <- d %>% left_join(trt_n, by = c("site", "MGT_trt_base")) %>%
      group_by(site_year) %>%
      mutate(N_rec = ifelse(site %in% AGCROS_UNIFORM_N & any(!is.na(N_rec)), median(N_rec, na.rm = TRUE), N_rec),
             N_kg_rec = case_when(!is.na(N_rec) ~ N_rec,
                                  site %in% AGCROS_SITES & !is.na(N_trt) ~ N_trt,
                                  site %in% zero_if_missing ~ 0,
                                  site %in% optimum_if_missing ~ OPT_N_KG,
                                  any(!is.na(N_rec)) ~ median(N_rec, na.rm = TRUE),
                                  TRUE ~ 0)) %>%
      ungroup() %>%
      mutate(N_leg = ifelse(is.na(N_rec) & site %in% optimum_if_missing, 0, legume_credit_kg(legume_class, lg)),
             N100  = (N_kg_rec + N_leg) / 100)
    d <- add_n_terms(d)                       # linear N per site; quadratic at N_QUAD_SITES
    nvs <- intersect(N_INTERACT_SITES, names(which(sapply(split(d, d$site), function(x) varies(x, "N100")))))
    d <- d %>% mutate(N_c = ifelse(site %in% nvs, N100 - N_REF / 100, 0))
    f <- reformulate(c("site_year", n_terms(d), o_term(d), "SHMI_w", mei_terms(d), "SHMI_w:N_c", "SHMI_w:age_c"),
                     "yield_kgha")
    m <- tryCatch(lme(f, random = ~ 1 | trt/plot_id, correlation = corCAR1(form = ~ year | trt/plot_id),
                      data = d, method = "ML", control = ctrl), error = function(e) { message("  failed: ", conditionMessage(e)); NULL })
    if (is.null(m)) { res2[[i]] <- tibble(legume = lg, manure = mn, note = "fit failed"); next }
    tt <- summary(m)$tTable
    org <- grep("^otherWICST_organic$", rownames(tt), value = TRUE)
    res2[[i]] <- tibble(legume = lg, manure = mn, n = nrow(d), N_mean = 100 * mean(d$N100), n_N_sites = length(nvs),
                        AIC = AIC(m), logLik = as.numeric(logLik(m)),
                        SHMI_slope = 10 * tt["SHMI_w", "Value"], SHMI_se = 10 * tt["SHMI_w", "Std.Error"],
                        SHMIxN = 10 * tt["SHMI_w:N_c", "Value"],
                        WICST_organic = if (length(org)) tt[org, "Value"] else NA_real_)
    out <- bind_rows(res2) %>% mutate(dAIC_best = AIC - min(AIC, na.rm = TRUE)) %>% arrange(AIC)
    write.csv(out %>% mutate(across(where(is.numeric), ~ round(.x, 1))), all_file, row.names = FALSE)
  }
  res_all <- bind_rows(res2) %>% mutate(dAIC_best = AIC - min(AIC, na.rm = TRUE)) %>% arrange(AIC)
  print(res_all, n = Inf, width = Inf)
}
