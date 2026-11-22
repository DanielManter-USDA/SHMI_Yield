################################################################################
# F5_records.R  (13 Nov 2026; was F7_records.R, renumbered 17 Nov 2026)   Fig. 5: routine records are enough
# Reads:  panel_records.rds (C1c, from B9's results)
# Writes: Outputs/figures/F5_records.png
#   a  what each simplified record keeps
#   b  soil-health performance and the corn yield slope with each record (median, 95% range)
################################################################################
source(here::here("R", "figures", "fig_common.R"), local = TRUE)   # local: keeps the pipeline's globals intact
rp <- read_panel("records")
k <- FONT_SCALE
pa <- enlarge(rp$a) + theme(axis.text.x.top = element_text(face = "bold", size = 9 * k),
                            axis.text.y = element_text(size = 10 * k, colour = "black", hjust = 1))
# b was drawn beside a, so its record names were hidden; stacked, it needs them back
pb <- enlarge(rp$b) + theme(axis.text.y = element_text(size = 10 * k, colour = "black", hjust = 1),
                            axis.ticks.y = element_line(),
                            strip.text = element_text(face = "bold", size = 9 * k, hjust = 0))
fig <- pa / pb + plot_layout(heights = c(0.8, 1))
save_fig(tag_fig(fig), "F5_records.png", width = 12.5, height = 8.5)
