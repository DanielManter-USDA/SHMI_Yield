################################################################################
# F1_soil_health_model.R  (13 Nov 2026; two panels 17 Nov 2026)   Fig. 1: soil health and SHMI as systems
#   a  the soil-health model (lv_SH): the second-order SEM below
#   b  how SHMI scores a system: the path diagram among its sub-indices (C1's panel_path.rds)
# Reads:  Outputs/sem_noMX_standardized.csv (A8; lavaan standardizedSolution of the A1c model),
#         Outputs/panels/panel_path.rds (C1)
# Writes: Outputs/figures/F1_soil_health_model.png
# Draws the second-order SEM: soil health -> three factors -> 11 indicators, standardized
# loadings on the arrows (dashed: P > 0.05; * constrained equal), residual correlations
# between indicators as red curves under the boxes. Positions are computed from the table, so
# the figure follows the model if it changes. Edit the name maps below to relabel.
################################################################################
source(here::here("R", "figures", "fig_common.R"), local = TRUE)   # local: keeps the pipeline's globals intact

FACTOR_LAB <- c(SH.Bio = "Biological activity", SH.CN = "Carbon & nitrogen", SH.Phys = "Physical condition")
FACTOR_COL <- c(SH.Bio = "#cfcbf5", SH.CN = "#9fe1cb", SH.Phys = "#f9c677")
IND_SUB <- c(Cmin = "C mineralization", Nmin = "N mineralization", ARYL = "arylsulfatase", NAG = "NAGase",
             SOC = "organic C", TN = "total N", POXC = "active C", ACE = "ACE protein",
             SSI = "slaking index", WAS = "aggregate stab.", Ksat = "sat. hydr. cond.")
TOP_LAB <- "Soil health (lv_SH)"

s <- read_tab("sem")
ld <- s %>% filter(op == "=~")
latents <- unique(ld$lhs)
top <- setdiff(latents, ld$rhs)[1]                           # the second-order factor
fac <- ld %>% filter(lhs == top) %>% pull(rhs)
ind <- ld %>% filter(lhs %in% fac)
sig <- function(p) is.na(p) | p <= 0.05

# ---- positions ---------------------------------------------------------------------------
GAP <- 0.6; Y_TOP <- 3.2; Y_FAC <- 2.05; Y_IND <- 0.75; BW <- 0.94; BH <- 0.42
ind <- ind %>% group_by(lhs) %>% mutate(k = row_number()) %>% ungroup() %>%
  mutate(g = match(lhs, fac)) %>% arrange(g, k) %>%
  mutate(x = seq_len(n()) + (g - 1) * GAP, x0 = x)
fpos <- ind %>% group_by(lhs) %>% summarise(x = mean(x), .groups = "drop") %>% rename(rhs = lhs)
tx <- mean(fpos$x)
# loadings sharing a lavaan label were constrained equal (A1c): marked with *
f2 <- ld %>% filter(lhs == top) %>% left_join(fpos, by = "rhs") %>%
  mutate(lb = if ("label" %in% names(.)) coalesce(as.character(label), "") else "",
         eq = nzchar(lb) & duplicated(lb) | nzchar(lb) & duplicated(lb, fromLast = TRUE),
         lab = sprintf("%.2f%s", est.std, ifelse(eq, "*", "")))

ell <- function(id, x, y, a = 0.95, b = 0.27) {
  t <- seq(0, 2 * pi, length.out = 80); tibble(id = id, x = x + a * cos(t), y = y + b * sin(t))
}
E <- bind_rows(ell("top", tx, Y_TOP, 1.05), bind_rows(lapply(seq_len(nrow(fpos)), function(i) ell(fpos$rhs[i], fpos$x[i], Y_FAC))))
E$fill <- ifelse(E$id == "top", "#f1efe8", FACTOR_COL[E$id]); E$fill[is.na(E$fill)] <- "#dddddd"

# arrows from the top node (tx) to each factor (x from fpos); xend is set before x is overwritten
seg_top <- f2 %>% transmute(xend = x, x = tx, y = Y_TOP - 0.27, yend = Y_FAC + 0.27, lab, dashed = !sig(pvalue))
seg_ind <- ind %>% transmute(x = fpos$x[match(lhs, fpos$rhs)], y = Y_FAC - 0.27, xend = x0, yend = Y_IND + BH / 2,
                             lab = sprintf("%.2f", est.std), dashed = !sig(pvalue))
labpos <- function(d, f) d %>% mutate(lx = x + f * (xend - x), ly = y + f * (yend - y))

# residual correlations between indicators (red curves under the boxes)
rc <- s %>% filter(op == "~~", lhs != rhs, lhs %in% ind$rhs, rhs %in% ind$rhs) %>%
  mutate(x1 = ind$x[match(lhs, ind$rhs)], x2 = ind$x[match(rhs, ind$rhs)],
         x = pmin(x1, x2), xend = pmax(x1, x2), dx = xend - x, dashed = !sig(pvalue)) %>% arrange(dx)
# each curve is drawn as a parabola in data units (depth grows with span, capped), label at its lowest point
Y_RC <- Y_IND - BH / 2 - 0.03
rc <- rc %>% mutate(depth = pmin(0.12 * dx, 0.6) + 0.12, lx = (x + xend) / 2, ly = Y_RC - depth,
                    lab = sprintf("%.2f", est.std), id = row_number())
rc_path <- if (nrow(rc)) bind_rows(lapply(seq_len(nrow(rc)), function(i) {
  u <- seq(0, 1, length.out = 60)
  tibble(id = i, x = rc$x[i] + u * rc$dx[i], y = Y_RC - rc$depth[i] * 4 * u * (1 - u), dashed = rc$dashed[i])
})) else tibble()
lt <- function(d) ifelse(d, "dashed", "solid")
p <- ggplot() +
  geom_segment(data = seg_top, aes(x, y, xend = xend, yend = yend, linetype = lt(dashed)),
               colour = "grey35", arrow = arrow(length = unit(2.2, "mm"), type = "closed")) +
  geom_segment(data = seg_ind, aes(x, y, xend = xend, yend = yend, linetype = lt(dashed),
                                   colour = ifelse(dashed, "grey65", "grey35")),
               arrow = arrow(length = unit(2.2, "mm"), type = "closed")) +
  geom_label(data = labpos(seg_top, 0.5), aes(lx, ly, label = lab), size = 3.6, fontface = "bold",
             label.size = 0, fill = "white", label.padding = unit(0.1, "lines")) +
  geom_label(data = labpos(seg_ind, 0.55), aes(lx, ly, label = lab), size = 3.2, fontface = "bold",
             label.size = 0, fill = "white", label.padding = unit(0.08, "lines")) +
  { if (nrow(rc)) geom_path(data = rc_path, aes(x, y, group = id, linetype = lt(dashed), alpha = ifelse(dashed, 0.45, 1)),
                            colour = "#a0391a", linewidth = 0.6,
                            arrow = arrow(length = unit(1.6, "mm"), ends = "both", type = "closed")) } +
  { if (nrow(rc)) geom_label(data = rc, aes(lx, ly, label = lab), colour = "#a0391a", size = 3, fontface = "bold",
                             label.size = 0, fill = "white", label.padding = unit(0.06, "lines")) } +
  geom_polygon(data = E, aes(x, y, group = id, fill = fill), colour = "grey30", linewidth = 0.5) +
  geom_text(data = tibble(x = tx, y = Y_TOP, l = TOP_LAB), aes(x, y, label = l), size = 4.4, fontface = "bold") +
  geom_text(data = fpos, aes(x, Y_FAC, label = ifelse(rhs %in% names(FACTOR_LAB), FACTOR_LAB[rhs], rhs)),
            size = 4.2, fontface = "bold") +
  geom_tile(data = ind, aes(x, Y_IND), width = BW, height = BH, fill = "white", colour = "grey55", linewidth = 0.4) +
  geom_text(data = ind, aes(x, Y_IND + 0.07, label = rhs), size = 3.8, fontface = "bold") +
  geom_text(data = ind, aes(x, Y_IND - 0.11, label = ifelse(rhs %in% names(IND_SUB), IND_SUB[rhs], "")),
            size = 2.3, colour = "grey35") +
  scale_fill_identity() + scale_colour_identity() + scale_alpha_identity() +
  scale_linetype_identity() +
  coord_cartesian(xlim = range(ind$x) + c(-0.6, 0.6), ylim = c(Y_RC - max(c(rc$depth, 0)) - 0.3, Y_TOP + 0.35),
                  expand = FALSE, clip = "off") +
  theme_void()
# b: the path diagram saved by C1 (was the top of Fig. 2 until 17 Nov 2026)
pb <- enlarge(read_panel("path"), 1.3)
fig <- wrap_elements(full = enlarge(p, 1.25)) / pb + plot_layout(heights = c(5.6, 4.6))
save_fig(tag_fig(fig), "F1_soil_health_model.png", width = 12.5, height = 10.6)
