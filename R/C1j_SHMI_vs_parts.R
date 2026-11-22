################################################################################
# C1j_SHMI_vs_parts.R  (17 Nov 2026)   Table 2: SHMI against its parts
#
# Does SHMI describe soil health better than any one of its parts? Replaces the standalone
# SHMI_vs_parts_aicc.R and SHMI_vs_parts_robust.R (the free four-part combination is dropped:
# fitted to lv_SH on the same data, it is nearly SHMI's own calibration).
#   1. ΔAICc: lv_SH ~ score + (1 | site), ML fits (lme4). SHMI is charged N_TUNED extra
#      parameters for the values its calibration tuned against lv_SH (free Cover/InvDist split,
#      Diversity weight, growing-season share, livestock credit, amendment share). AICc with
#      n = systems (reported) and n = sites (same ranking).
#   2. Leave one site out: within-site R² (site-demeaned outcome and score; the slope refitted
#      without the site predicts the differences among its systems). Differences from SHMI get
#      95% intervals from a site bootstrap (B resamples).
#   3. The same for SHAPE, a soil-health score calculated without reference to the calibration.
# Corn yield (the last columns of Table 2) comes from B2d (yield_parts.csv); C2 joins the two.
#
# Reads   Outputs/SHAPE/SHAPE_by_unit_with_SHMI.csv (Part A; else data/SHAPE_by_unit_with_SHMI.csv)
# Writes  Outputs/Parts/SHMI_vs_parts.csv  (one row per outcome x score)
################################################################################
suppressPackageStartupMessages({ library(dplyr); library(lme4) })
OUT <- here::here("Outputs", "Parts"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
if (!exists("N_TUNED")) N_TUNED <- 5
if (!exists("B_BOOT"))  B_BOOT  <- 2000
set.seed(1)

f <- here::here("Outputs", "SHAPE", "SHAPE_by_unit_with_SHMI.csv")
if (!file.exists(f)) f <- here::here("data", "SHAPE_by_unit_with_SHMI.csv")
d <- read.csv(f) %>% filter(!is.na(SH), !is.na(SHMI), !is.na(SHAPE_mean)) %>% mutate(sitecode = factor(sitecode))
n_obs <- nrow(d); n_sites <- nlevels(d$sitecode)
message(sprintf("C1j: %d systems at %d sites (%s)", n_obs, n_sites, basename(dirname(f))))

SCORES <- c("SHMI" = "SHMI", "Cover/Roots" = "Cover", "Diversity" = "Diversity",
            "Inverse Disturbance" = "InvDist", "Amendments" = "p_amend")
OUTCOMES <- c(lv_SH = "SH", SHAPE = "SHAPE_mean")

# ---- 1. ΔAICc (random site intercept) ----------------------------------------------------
aicc_tab <- function(y) {
  t <- bind_rows(lapply(names(SCORES), function(nm) {
    m <- lmer(as.formula(paste(y, "~", SCORES[[nm]], "+ (1 | sitecode)")), data = d, REML = FALSE)
    ll <- logLik(m); k <- attr(ll, "df") + if (nm == "SHMI") N_TUNED else 0
    aic <- -2 * as.numeric(ll) + 2 * k
    tibble(score = nm, k = k, AIC = aic,
           AICc_obs = aic + 2 * k * (k + 1) / (n_obs - k - 1),
           AICc_sites = aic + 2 * k * (k + 1) / (n_sites - k - 1))
  }))
  t %>% mutate(dAICc = AICc_obs - AICc_obs[score == "SHMI"],
               dAICc_sites = AICc_sites - AICc_sites[score == "SHMI"])
}

# ---- 2. Leave one site out, within sites -------------------------------------------------
wd <- function(v) v - ave(v, d$sitecode)
loso_tab <- function(y) {
  Y <- wd(d[[y]]); site <- d$sitecode
  P <- sapply(names(SCORES), function(nm) {
    x <- wd(d[[SCORES[[nm]]]]); p <- numeric(length(Y))
    for (s in levels(site)) {
      tr <- site != s
      b <- sum(x[tr] * Y[tr]) / sum(x[tr]^2)        # slope through the origin (both demeaned)
      p[!tr] <- b * x[!tr]
    }
    p
  })
  r2 <- function(idx) apply(P[idx, , drop = FALSE], 2, function(p) 1 - sum((Y[idx] - p)^2) / sum(Y[idx]^2))
  est <- r2(seq_along(Y))
  rows <- split(seq_along(Y), site)
  boot <- replicate(B_BOOT, { s <- sample(levels(site), replace = TRUE); v <- r2(unlist(rows[s])); v - v[["SHMI"]] })
  tibble(score = names(SCORES), R2_site_out = est, d_R2 = est - est[["SHMI"]],
         d_lo = apply(boot, 1, quantile, 0.025), d_hi = apply(boot, 1, quantile, 0.975))
}

out <- bind_rows(lapply(names(OUTCOMES), function(o) {
  y <- OUTCOMES[[o]]
  aicc_tab(y) %>% left_join(loso_tab(y), by = "score") %>% mutate(outcome = o, .before = 1)
})) %>% mutate(n_systems = n_obs, n_sites = n_sites, n_tuned = N_TUNED, B = B_BOOT)
write.csv(mutate(out, across(where(is.double), ~ signif(.x, 5))), file.path(OUT, "SHMI_vs_parts.csv"), row.names = FALSE)
print(out %>% select(outcome, score, dAICc, dAICc_sites, R2_site_out, d_R2, d_lo, d_hi) %>%
        mutate(across(where(is.double), ~ round(.x, 3))), n = Inf, width = Inf)
message("Saved ", file.path("Outputs", "Parts", "SHMI_vs_parts.csv"))
