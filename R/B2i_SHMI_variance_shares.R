################################################################################
# B2i_SHMI_variance_shares.R
#
# How SHMI varies within sites in the two datasets (Supplementary Methods, "Composition
# of SHMI variation in the two datasets"; Fig. 6a). Fig. 6a compares the principles for
# soil health in the NAPESHM systems and for corn yield in the ten experiments, so this
# checks that SHMI's parts varied similarly in both.
#
# Share of within-site SHMI variance from each part = cov(part, SHMI) / var(SHMI), where a
# part is weight x sub-index (Cover/Roots, Diversity, Inverse Disturbance) or the amendment
# bonus. The parts add up to SHMI, so the four shares sum to 1; a negative share means that
# part moved against SHMI.
#   pooled  deviations from site means, pooled over sites. Equivalent to weighting each
#           site's shares by its within-site SHMI variance and number of systems.
#   median  median of the per-site shares, over sites where SHMI varies (as in
#           Supplementary Table 2 for the experiments)
#   sd      median within-site standard deviation of SHMI, over the same sites
# NAPESHM: US systems at sites with at least two systems; all systems, and annual systems
#   only (system == "annual", as defined in A1b_Prepare_SEM.R).
# Experiments: each treatment record once, as in B2's SHMI_w and Supplementary Table 2.
#
# Input   Outputs/Final_plus/final_results.rds (A3a)
#         Outputs/yield_main_default_total_agcros/B2_data.rds (B2)
# Outputs Outputs/variance_shares/ variance_shares.csv (one row per dataset),
#         variance_shares_sites.csv (one row per site)
################################################################################
suppressPackageStartupMessages({ library(dplyr) })
setwd(here::here())
IN_SOIL  <- "Outputs/Final_plus/final_results.rds"
IN_YIELD <- "Outputs/yield_main_default_total_agcros/B2_data.rds"
OUT      <- "Outputs/variance_shares"
for (f in c(IN_SOIL, IN_YIELD)) if (!file.exists(f)) stop(f, " not found: run A3a and B2 first")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

PARTS <- c("cover", "div", "dist", "amend")
NEEDED <- c("site", "SHMI", "SHMI_base", "Cover", "Diversity", "InvDist", "Amend_bonus")

# weight x sub-index for each part. The weights are recovered from the scores
# (SHMI_base = wC Cover + wD Diversity + wI InvDist), as in B2, so overrides are respected.
add_parts <- function(d) {
  miss <- setdiff(NEEDED, names(d))
  if (length(miss)) stop("missing columns: ", paste(miss, collapse = ", "))
  wfit <- lm(SHMI_base ~ 0 + Cover + Diversity + InvDist, data = d)
  w <- coef(wfit)
  gap <- max(abs(fitted(wfit) - d$SHMI_base))
  if (gap > 0.01) warning("SHMI weights not recovered exactly (max gap ", signif(gap, 3), " points)")
  d %>% mutate(p_cover = w[["Cover"]] * Cover, p_div = w[["Diversity"]] * Diversity,
               p_dist = w[["InvDist"]] * InvDist, p_amend = Amend_bonus)
}

shares <- function(d, label) {
  d <- add_parts(d) %>% group_by(site) %>% filter(n() >= 2) %>% ungroup()
  gap2 <- max(abs(d$p_cover + d$p_div + d$p_dist + d$p_amend - d$SHMI))
  if (gap2 > 0.01) warning(label, ": parts do not add up to SHMI (max gap ", signif(gap2, 3), ")")

  # pooled over sites: deviations from site means
  dw <- d %>% group_by(site) %>%
    mutate(across(c(SHMI, all_of(paste0("p_", PARTS))), ~ .x - mean(.x))) %>% ungroup()
  vs <- sum(dw$SHMI^2)
  pooled <- vapply(PARTS, function(k) sum(dw[[paste0("p_", k)]] * dw$SHMI) / vs, numeric(1))

  # each site on its own (sites where SHMI varies)
  per_site <- d %>% group_by(site) %>% filter(var(SHMI) > 1e-8) %>%
    summarise(dataset = label, systems = n(), SHMI_sd = sd(SHMI),
              share_cover = cov(p_cover, SHMI) / var(SHMI), share_div = cov(p_div, SHMI) / var(SHMI),
              share_dist = cov(p_dist, SHMI) / var(SHMI), share_amend = cov(p_amend, SHMI) / var(SHMI),
              .groups = "drop") %>%
    relocate(dataset)

  summ <- tibble(dataset = label, systems = nrow(d), sites = n_distinct(d$site),
                 sites_varying = nrow(per_site),
                 pooled_cover = pooled[["cover"]], pooled_div = pooled[["div"]],
                 pooled_dist = pooled[["dist"]], pooled_amend = pooled[["amend"]],
                 median_cover = median(per_site$share_cover), median_div = median(per_site$share_div),
                 median_dist = median(per_site$share_dist), median_amend = median(per_site$share_amend),
                 median_within_sd = median(per_site$SHMI_sd))
  list(summary = summ, sites = per_site)
}

# ---- NAPESHM (soil health; US systems used in the calibration) ----------------------
soil <- readRDS(IN_SOIL)$us %>% mutate(site = as.character(sitecode))
if (!"system" %in% names(soil)) stop("final_results.rds$us has no `system` column")

# ---- the ten experiments (corn yield; each treatment record once) -------------------
corn <- readRDS(IN_YIELD)$corn %>% mutate(site = as.character(site))
recs <- corn %>% distinct(site, MGT_combo, .keep_all = TRUE)

res <- list(shares(soil, "NAPESHM, all systems"),
            shares(filter(soil, system == "annual"), "NAPESHM, annual systems"),
            shares(recs, "Ten experiments"))
summary_tab <- bind_rows(lapply(res, `[[`, "summary"))
site_tab    <- bind_rows(lapply(res, `[[`, "sites"))

write.csv(mutate(summary_tab, across(where(is.double), ~ round(.x, 3))),
          file.path(OUT, "variance_shares.csv"), row.names = FALSE)
write.csv(mutate(site_tab, across(where(is.double), ~ round(.x, 3))),
          file.path(OUT, "variance_shares_sites.csv"), row.names = FALSE)

cat("\nShare of within-site SHMI variance (Cover / Diversity / Disturbance / Amendments)\n")
for (i in seq_len(nrow(summary_tab))) with(summary_tab[i, ], cat(sprintf(
  "%-26s %3d systems, %2d sites | pooled %.2f / %.2f / %.2f / %.2f | median %.2f / %.2f / %.2f / %.2f | median SD %.1f\n",
  dataset, systems, sites, pooled_cover, pooled_div, pooled_dist, pooled_amend,
  median_cover, median_div, median_dist, median_amend, median_within_sd)))
