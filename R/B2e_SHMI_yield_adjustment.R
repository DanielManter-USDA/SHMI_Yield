################################################################################
# B2e_SHMI_yield_adjustment.R
#
# What "within sites" means, in three steps, on the yield scale (kg/ha), for readers
# and talks. Same data and model as B2 (Fig. 3a).
#   a  Raw: each treatment's mean corn yield at each experiment against its SHMI.
#      Experiments differ in soil, climate and yield level, so they do not share a line.
#   b  Site and year removed: each plot-year minus its site-year mean (plus the overall
#      mean yield), averaged by treatment. Weather, soil and each experiment's level are
#      gone; the treatments line up.
#   c  Nitrogen, other treatment factors and system age also accounted for: the partial
#      residuals of B2's full model (yield ~ site x year + N per site + other factors +
#      SHMI + SHMI x N + SHMI x age, random treatment/plot, AR(1)) at 150 kg N/ha and the
#      average age, plus the overall mean yield; treatments with the same SHMI at a site
#      averaged (as Fig. 3a). One line, 220 kg/ha per 10 SHMI.
# Fig. 3a is panel c with a different baseline (county yield there; overall mean here).
#
# Input   Outputs/yield_main_default_total_agcros/B2_data.rds (B2)
# Output  Outputs/yield_adjustment/Fig_yield_adjustment.png/.pdf, yield_adjustment_points.csv
# Runs in a few minutes (one model fit).
################################################################################
suppressPackageStartupMessages({ library(dplyr); library(nlme); library(ggplot2); library(patchwork) })
setwd(here::here())
IN  <- "Outputs/yield_main_default_total_agcros/B2_data.rds"
OUT <- "Outputs/yield_adjustment"
if (!file.exists(IN)) stop(IN, " not found: run B2_SHMI_yield_main.R first")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
ctrl  <- lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200)
N_REF <- 150
site_cols <- c(NELITCSE = "#0072B2", NEMLTCRS = "#D55E00", OHTVD_Ho = "#56B4E9", OHTVD_Wo = "#CC79A7",
               WICST_ARS = "#009E73", NEMEIRR = "#E69F00", NECCIRR = "#8C510A", MNMOFS = "#999999",
               INWLACRE = "#000000", INWLTPAC = "#882255", IAAM7071 = "#F0E442", PAUP = "#6A3D9A")

D <- readRDS(IN); corn <- D$corn
corn <- corn %>% group_by(site) %>% mutate(SHMI_w = SHMI - mean(SHMI[!duplicated(MGT_combo)])) %>% ungroup()
GRAND <- mean(corn$yield_kgha)

# ---- the model of Fig. 3a (B2, section 6) ---------------------------------------------
varies  <- function(d, v) any(tapply(d[[v]], d$site_year, function(x) diff(range(x))) > 1e-8)
n_terms <- function(d) { v <- grep("^N[12]_", names(d), value = TRUE); v[vapply(v, function(k) varies(d, k), logical(1))] }
o_term  <- function(d) if (n_distinct(d$other) > 1) "other" else character()
f <- reformulate(c("site_year", n_terms(corn), o_term(corn), "SHMI_w", "SHMI_w:N_c", "SHMI_w:age_c"), "yield_kgha")
m <- lme(f, random = ~ 1 | trt/plot_id, correlation = corCAR1(form = ~ year | trt/plot_id),
         data = corn, method = "REML", control = ctrl)
tt <- summary(m)$tTable["SHMI_w", ]
b  <- tt[["Value"]]; se <- tt[["Std.Error"]]; DF <- tt[["DF"]]
cat(sprintf("Slope at %d kg N/ha and average age: %.0f kg/ha per 10 SHMI (95%% CI %.0f-%.0f)\n",
            N_REF, 10 * b, 10 * (b - qt(0.975, DF) * se), 10 * (b + qt(0.975, DF) * se)))
PIV <- mean(corn$SHMI)

# ---- the three steps ---------------------------------------------------------------------
corn <- corn %>% group_by(site_year) %>% mutate(y_sy = yield_kgha - mean(yield_kgha) + GRAND) %>% ungroup() %>%
  mutate(y_full = as.numeric(resid(m, level = 0)) + b * (SHMI - PIV) + GRAND)
by_trt <- function(v) corn %>% group_by(site, MGT_trt_base) %>%
  summarise(SHMI = mean(SHMI), y = mean(.data[[v]]), .groups = "drop")
pa <- by_trt("yield_kgha") %>% mutate(step = "a")
pb <- by_trt("y_sy") %>% mutate(step = "b")
pc <- corn %>% group_by(site, SHMI_level = round(SHMI)) %>%
  summarise(SHMI = mean(SHMI), y = mean(y_full), se = sd(y_full) / sqrt(n()), .groups = "drop") %>% mutate(step = "c")
write.csv(bind_rows(pa, pb, pc), file.path(OUT, "yield_adjustment_points.csv"), row.names = FALSE)
ln <- tibble(SHMI = seq(min(corn$SHMI), max(corn$SHMI), length.out = 200)) %>%
  mutate(fit = GRAND + b * (SHMI - PIV), half = qt(0.975, DF) * se * abs(SHMI - PIV))

YL <- range(c(pa$y, pb$y, pc$y)); YL <- YL + c(-0.04, 0.04) * diff(YL)    # one y scale for all three
base <- function(d, title, sub) ggplot(d, aes(SHMI, y)) +
  scale_fill_manual(values = site_cols, name = NULL) +
  coord_cartesian(ylim = YL) +
  labs(title = title, subtitle = sub, x = "Soil Health Management Index (SHMI)", y = expression("Corn yield (kg ha"^-1*")")) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11), plot.subtitle = element_text(size = 8.5, colour = "grey30"),
        panel.grid.minor = element_blank())
ga <- base(pa, "a  Raw yields", "Each treatment's mean yield at its experiment") +
  geom_point(aes(fill = site), shape = 21, size = 2.3, colour = "black")
gb <- base(pb, "b  Site and year removed", "Each yield minus its site-year mean (+ overall mean)") +
  geom_point(aes(fill = site), shape = 21, size = 2.3, colour = "black")
gc <- base(pc, "c  N, other treatments and age also removed", sprintf("Model of Fig. 3a at %d kg N/ha and average age", N_REF)) +
  geom_ribbon(data = ln, aes(SHMI, ymin = fit - half, ymax = fit + half), inherit.aes = FALSE, fill = "darkgreen", alpha = 0.15) +
  geom_line(data = ln, aes(SHMI, fit), colour = "darkgreen", linewidth = 1.1) +
  geom_errorbar(aes(ymin = y - se, ymax = y + se), width = 0, colour = "grey40") +
  geom_point(aes(fill = site), shape = 21, size = 2.3, colour = "black") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3, size = 3,
           label = sprintf("+%.0f kg/ha per 10 SHMI (95%% CI %.0f–%.0f)", 10 * b, 10 * (b - qt(0.975, DF) * se), 10 * (b + qt(0.975, DF) * se)))
fig <- (ga | gb | gc) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave(file.path(OUT, "Fig_yield_adjustment.png"), fig, width = 12, height = 4.6, dpi = 300)
tryCatch(ggsave(file.path(OUT, "Fig_yield_adjustment.pdf"), fig, width = 12, height = 4.6, device = cairo_pdf),
         error = function(e) message("PDF not written: ", conditionMessage(e)))
message("Yield adjustment figure: ", OUT)
