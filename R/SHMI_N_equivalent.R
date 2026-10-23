################################################################################
# SHMI_N_equivalent.R
#
# How much fertilizer N is 10 SHMI points worth?  Uses the fitted within-site model
# from SHMI_yield_main.R (m_full: site x year + a quadratic N curve per N-trial site
# + SHMI_w + SHMI_w:N_c + SHMI_w:age_c, AR(1)).
#
# At a reference N rate N0 (and the average system age):
#   gain(N0)      = 10 x SHMI slope at N0                       (kg grain/ha)
#   marginal(N0)  = slope of the fitted N curve at N0           (kg grain per kg N)
#   NE_margin     = gain / marginal   (kg N/ha; the N that gives the same yield at
#                   the margin, valid for small changes)
#   NE_curve      = the extra N, added on top of N0, that moves the fitted curve up by
#                   'gain' (solves the quadratic; NA when the curve plateaus first,
#                   i.e. no amount of N reproduces the SHMI gain)
# Reported for the average N-trial site (equal-weight mean of the site curves, as in
# Fig. 2b) and for each N-trial site within its own observed N range. Uncertainty:
# 5,000 draws from the fixed-effect estimates and their covariance (parametric
# simulation), so CIs carry the correlation between the SHMI and N terms.
# The old manuscript number (13 kg grain per kg N) is the average over 112-200 kg N,
# printed here for comparison.
#
# RUN: after SHMI_yield_main.R in the same R session (it uses m_full, corn,
#      n_var_sites, N_REF, AGE_MEAN, OUT); sources it if m_full is missing.
# Outputs (OUT): N_equivalent.csv, N_equivalent_by_site.csv, Fig_N_equivalent.png
################################################################################
if (!exists("m_full")) source(here::here("R", "SHMI_yield_main.R"))
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(nlme) })

N_GRID   <- seq(25, 250, by = 5)     # kg N/ha for the curves
N_REPORT <- c(50, 100, 150, 200)     # rows in the table
N_PRICE  <- 1.10                     # $ per kg N [placeholder; ~ $0.50/lb N]
NDRAW    <- 5000
set.seed(42)

b <- fixef(m_full); V <- vcov(m_full)
# the N-rate sites' curves (quadratic); linear-only sites cannot locate the curve's flattening
curve_sites <- if (exists("n_curve_sites")) n_curve_sites else n_var_sites
sites_n <- curve_sites[paste0("N1_", curve_sites) %in% names(b)]
stopifnot(length(sites_n) > 0)
keep <- c("SHMI_w", "SHMI_w:N_c", "SHMI_w:age_c", paste0("N1_", sites_n),
          intersect(paste0("N2_", sites_n), names(b)))
B <- MASS::mvrnorm(NDRAW, b[keep], V[keep, keep])          # draws x terms
B <- rbind(est = b[keep], B)                               # row 1 = point estimate
co <- function(k) if (k %in% colnames(B)) B[, k] else rep(0, nrow(B))

# curve coefficients (per N100 = N/100) for one site, or the equal-weight average
curve_coef <- function(s = NULL) {
  ss <- if (is.null(s)) sites_n else s
  list(b1 = rowMeans(do.call(cbind, lapply(ss, function(x) co(paste0("N1_", x))))),
       b2 = rowMeans(do.call(cbind, lapply(ss, function(x) co(paste0("N2_", x))))))
}
gain_at <- function(N0) 10 * (co("SHMI_w") + (N0 - N_REF) / 100 * co("SHMI_w:N_c"))  # age at mean
ne_at <- function(N0, cc) {
  g    <- gain_at(N0)
  marg <- (cc$b1 + 2 * cc$b2 * N0 / 100) / 100               # kg grain per kg N at N0
  # extra N d (kg/ha): (b2/1e4) d^2 + (b1/100 + 2 b2 N0/1e4) d - g = 0
  A <- cc$b2 / 1e4; Bq <- cc$b1 / 100 + 2 * cc$b2 * N0 / 1e4
  disc <- Bq^2 + 4 * A * g
  d <- ifelse(abs(A) < 1e-12, g / Bq,
              ifelse(disc < 0, NA, (-Bq + sqrt(pmax(disc, 0))) / (2 * A)))
  d[!is.na(d) & (d < 0 | Bq <= 0)] <- NA                     # curve already flat or falling
  data.frame(draw = seq_along(g), gain = g, marginal = marg,
             NE_margin = ifelse(marg > 0, g / marg, NA), NE_curve = d)
}
summ <- function(x) {
  est <- x[1]; x <- x[-1]
  c(est = est, lo = unname(quantile(x, 0.025, na.rm = TRUE)), hi = unname(quantile(x, 0.975, na.rm = TRUE)),
    share_NA = mean(is.na(x)))
}
tab_for <- function(cc, label, Ns) {
  bind_rows(lapply(Ns, function(N0) {
    r <- ne_at(N0, cc)
    bind_rows(lapply(c("gain", "marginal", "NE_margin", "NE_curve"), function(v)
      data.frame(curve = label, N = N0, quantity = v, t(summ(r[[v]])))))
  }))
}

# ---- average N-trial site ---------------------------------------------------------
cc_avg <- curve_coef()
tab_avg <- tab_for(cc_avg, "Average N-trial site", N_GRID)
avg_112_200 <- with(cc_avg, (b1 * (2.00 - 1.12) + b2 * (2.00^2 - 1.12^2)) / 88)
cat(sprintf("\nAverage response 112 -> 200 kg N (manuscript's 13 kg/kg): %.1f kg grain per kg N (95%% CI %.1f-%.1f)\n",
            avg_112_200[1], quantile(avg_112_200[-1], .025), quantile(avg_112_200[-1], .975)))

# ---- each N-trial site, within its own N range ------------------------------------------
rng <- corn %>% filter(site %in% sites_n) %>% group_by(site) %>%
  summarise(N_min = 100 * min(N100), N_max = 100 * max(N100), .groups = "drop")
tab_site <- bind_rows(lapply(sites_n, function(s) {
  r <- rng[rng$site == s, ]
  Ns <- N_REPORT[N_REPORT >= r$N_min & N_REPORT <= r$N_max]
  if (!length(Ns)) return(NULL)
  tab_for(curve_coef(s), s, Ns)
})) %>% left_join(rng, by = c(curve = "site"))

# ---- tables -------------------------------------------------------------------------
wide <- function(t) t %>% filter(N %in% N_REPORT) %>%
  mutate(txt = ifelse(quantity == "marginal", sprintf("%.1f (%.1f to %.1f)", est, lo, hi),
                      sprintf("%.0f (%.0f to %.0f)", est, lo, hi)),
         txt = ifelse(share_NA > 0.025, paste0(txt, sprintf(" [%.0f%% of draws beyond plateau]", 100 * share_NA)), txt)) %>%
  select(curve, N, quantity, txt) %>% pivot_wider(names_from = quantity, values_from = txt)
out_avg <- wide(tab_avg) %>%
  left_join(tab_avg %>% filter(N %in% N_REPORT, quantity == "NE_curve") %>%
              transmute(N, value_usd_ha = round(est * N_PRICE)), by = "N")
cat("\nN-equivalent of 10 SHMI points (average N-trial site; estimate and 95% CI):\n")
print(as.data.frame(out_avg), row.names = FALSE)
cat("\nBy site (within each site's observed N range):\n")
print(as.data.frame(wide(tab_site)), row.names = FALSE)
cat("\nObserved N ranges (kg/ha):\n"); print(as.data.frame(rng), row.names = FALSE)
write.csv(out_avg, file.path(OUT, "N_equivalent.csv"), row.names = FALSE)
write.csv(wide(tab_site), file.path(OUT, "N_equivalent_by_site.csv"), row.names = FALSE)

# ---- figure ---------------------------------------------------------------------------
lab <- c(gain = "a  Yield gain from 10 SHMI points (kg/ha)",
         marginal = "b  Marginal N response (kg grain per kg N)",
         NE_curve = "c  N-equivalent of 10 SHMI points (kg N/ha)")
pd <- tab_avg %>% filter(quantity %in% names(lab)) %>% mutate(panel = factor(lab[quantity], lab))
rug <- data.frame(N = 100 * corn$N100[corn$site %in% sites_n])
p <- ggplot(pd, aes(N, est)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#4292c6", alpha = 0.2) +
  geom_line(colour = "#08519c", linewidth = 0.9) +
  geom_rug(data = rug, aes(N), inherit.aes = FALSE, sides = "b", alpha = 0.05) +
  facet_wrap(~ panel, scales = "free_y", ncol = 3) +
  labs(x = "N applied (kg/ha: fertilizer + available manure N + legume credit)", y = NULL,
       caption = paste("Average of the quadratic N curves at the N-trial sites (", paste(sites_n, collapse = ", "),
                       "), at average system age; 95% bands from 5,000 draws of the fixed effects.",
                       "\nc: extra N on top of the x-axis rate that raises the fitted curve as much as 10 SHMI points; gaps = the curve plateaus first.")) +
  theme_bw(base_size = 10) + theme(strip.background = element_blank(), strip.text = element_text(hjust = 0, face = "bold"),
                                   plot.caption = element_text(hjust = 0, size = 8))
ggsave(file.path(OUT, "Fig_N_equivalent.png"), p, width = 10, height = 3.6, dpi = 300)
cat("\nSaved N_equivalent.csv, N_equivalent_by_site.csv and Fig_N_equivalent.png in", OUT, "\n")
