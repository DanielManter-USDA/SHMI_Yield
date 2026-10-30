################################################################################
# C1d_practice_additivity.R
#
# Do practices add up in MEASURED soil health? (no SHMI involved)
# Annual NAPESHM systems with all four practice labels defined; within sites (site
# fixed effects, by demeaning), lv_SH is modelled on four practice indicators and
# their six pairwise interactions:
#   cover     cover crop (grasses, mix, legume) vs none
#   organic   organic or organic + synthetic nutrients vs synthetic only
#   reduced   reduced tillage (none, minimum) vs moderate or intense
#   rotation  two or more crops vs one
# The main effects are each practice's effect on its own (from a system with none of
# them); "added up" is their sum; "together" adds the interactions. 95% intervals from
# resampling sites (each drawn site keeps its own fixed effect).
# A model without interactions is also fitted: within these experiments the practices
# are nearly uncorrelated, so the overlap seen in SHMI's parts (Fig. 1a) does not
# appear as correlated labels; non-additivity appears as negative interactions.
#
# Input   data/SHAPE_by_unit_with_SHMI.csv (lv_SH = SH), Inputs/napeshm/treatments.csv
# Output  Outputs/Practice_additivity/additivity_effects.csv, additivity_bootstrap.csv,
#         Fig_additivity.png (talk/Supplementary candidate)
# Runs in about a minute.
################################################################################
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })
OUT <- here::here("Outputs", "Practice_additivity")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
B_BOOT <- if (exists("ADD_BOOT")) ADD_BOOT else 2000
set.seed(2026)

d <- read.csv(here::here("data", "SHAPE_by_unit_with_SHMI.csv"))
trt <- read.csv(here::here("Inputs", "napeshm", "treatments.csv"), stringsAsFactors = FALSE) %>%
  mutate(MGT_combo = paste0("NAPESHM_", sitecode, "_", treatmentid),
         across(c(eu_tillage, eu_cover_crops, eu_crop_rot, eu_nutr_type), ~ tolower(trimws(.x)))) %>%
  select(MGT_combo, eu_tillage, eu_cover_crops, eu_crop_rot, eu_nutr_type)
a <- d %>% filter(system == "annual") %>% inner_join(trt, by = "MGT_combo") %>%
  transmute(sitecode, SH,
            cover    = case_when(eu_cover_crops %in% c("grasses", "mix", "legume") ~ 1, eu_cover_crops == "none" ~ 0),
            organic  = case_when(eu_nutr_type %in% c("organic", "synthetic_and_organic") ~ 1, eu_nutr_type == "synthetic" ~ 0),
            reduced  = case_when(eu_tillage %in% c("none", "minimum") ~ 1, eu_tillage %in% c("intense", "moderate") ~ 0),
            rotation = case_when(eu_crop_rot %in% c("2", "3", "greater_than_3") ~ 1, eu_crop_rot == "1" ~ 0)) %>%
  filter(if_all(c(cover, organic, reduced, rotation), ~ !is.na(.x)))
P <- c("cover", "organic", "reduced", "rotation")
PAIRS <- combn(P, 2, simplify = FALSE)
for (pq in PAIRS) a[[paste(pq, collapse = "_x_")]] <- a[[pq[1]]] * a[[pq[2]]]
COLS <- c(P, vapply(PAIRS, paste, "", collapse = "_x_"))
cat(sprintf("%d annual systems at %d sites; sites where each term varies:\n", nrow(a), n_distinct(a$sitecode)))
print(sapply(COLS, function(k) sum(tapply(a[[k]], a$sitecode, function(x) length(unique(x)) > 1))))

# within-site least squares (demeaned by site = site fixed effects)
fit <- function(df, cols, grp = df$sitecode) {
  dm <- function(x) x - ave(x, grp)
  X <- sapply(cols, function(k) dm(df[[k]])); y <- dm(df$SH)
  b <- qr.coef(qr(X), y); b[is.na(b)] <- 0; b
}
summ <- function(b) {
  three <- c("cover", "reduced", "rotation", "cover_x_reduced", "cover_x_rotation", "reduced_x_rotation")
  c(b[P], added = sum(b[P]), together = sum(b[COLS]), interactions = sum(b[setdiff(COLS, P)]),
    added3 = sum(b[c("cover", "reduced", "rotation")]), together3 = sum(b[three]))
}
est  <- summ(fit(a, COLS))
nofx <- fit(a, P)                      # no interactions
solo <- sapply(P, function(k) fit(a, k))

sites <- unique(a$sitecode)
boot <- t(replicate(B_BOOT, {
  s <- sample(sites, length(sites), replace = TRUE)
  bb <- bind_rows(lapply(seq_along(s), function(j) mutate(filter(a, sitecode == s[j]), draw = j)))
  summ(fit(bb, COLS, grp = bb$draw))
}))
write.csv(as.data.frame(boot), file.path(OUT, "additivity_bootstrap.csv"), row.names = FALSE)
q <- function(x, p) unname(quantile(x, p))
res <- tibble(term = names(est), estimate = unname(est),
              lower = apply(boot, 2, q, 0.025), upper = apply(boot, 2, q, 0.975),
              share_below_0 = colMeans(boot < 0))
res <- bind_rows(res,
  tibble(term = "ratio together/added (three practices)", estimate = est[["together3"]] / est[["added3"]],
         lower = q(boot[, "together3"] / boot[, "added3"], 0.025), upper = q(boot[, "together3"] / boot[, "added3"], 0.975)),
  tibble(term = paste0(P, " (alone, no other terms)"), estimate = unname(solo)),
  tibble(term = paste0(P, " (all four, no interactions)"), estimate = unname(nofx)))
write.csv(mutate(res, across(where(is.double), ~ round(.x, 3))), file.path(OUT, "additivity_effects.csv"), row.names = FALSE)
cat("\nMeasured soil health (lv_SH points), within sites; 95% intervals resampling sites\n")
print(as.data.frame(mutate(res, across(where(is.double), ~ round(.x, 2)))), row.names = FALSE)

# figure: each practice alone (own colour); added up = the four stacked; all four together =
# the same four colours, each shrunk in proportion so they sum to the estimate (the split of
# the overlap among practices is not estimated: ADD_SPLIT <- FALSE draws one solid bar instead)
COL <- c(cover = "#2f7fb5", organic = "#b55d1e", reduced = "#7a4e2a", rotation = "#c99310")
NAME <- c(cover = "Cover crop", organic = "Organic inputs", reduced = "Reduced tillage", rotation = "Rotation")
XL <- c(cover = "Cover\ncrop", organic = "Organic\ninputs", reduced = "Reduced\ntillage", rotation = "Rotation",
        added = "Added up,\none at a time", together = "All four\ntogether")
split <- !exists("ADD_SPLIT") || isTRUE(ADD_SPLIT)
m <- est[P]; tog <- est[["together"]]
fd <- bind_rows(
  tibble(x = P, practice = P, value = m),
  tibble(x = "added", practice = P, value = m),
  if (split) tibble(x = "together", practice = P, value = m * tog / sum(m))
  else tibble(x = "together", practice = "together", value = tog)) %>%
  mutate(x = factor(XL[x], XL), practice = factor(practice, c(rev(P), "together")))
tot <- tibble(x = factor(XL[c("added", "together")], XL), value = c(est[["added"]], tog))
g <- ggplot(fd, aes(x, value, fill = practice)) +
  geom_col(width = 0.65, colour = "white", linewidth = 0.3) +
  geom_text(data = tot, aes(x, value, label = sprintf("+%.0f", value)), inherit.aes = FALSE, vjust = -0.5, size = 4.2, fontface = "bold") +
  geom_hline(yintercept = 0, colour = "grey50") +
  scale_fill_manual(values = c(COL, together = "#2e5e2e"), breaks = P, labels = NAME[P], name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = "Change in measured soil health (lv_SH points)") +
  theme_classic(base_size = 12) + theme(legend.position = "none")
ggsave(file.path(OUT, "Fig_additivity.png"), g, width = 7.5, height = 4.2, dpi = 300)
message("Practice additivity done: ", OUT)
