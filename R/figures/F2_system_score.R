################################################################################
# F2_system_score.R  (13 Nov 2026; 17 Nov 2026: path diagram moved to Fig. 1b)
# Fig. 2: SHMI tracks measured soil health
# Reads:  panel_SHMI_vs_SH.rds (C1h)
# Writes: Outputs/figures/F2_system_score.png
#   a  all systems   b  between site means   c  within sites
################################################################################
source(here::here("R", "figures", "fig_common.R"), local = TRUE)   # local: keeps the pipeline's globals intact
sh <- lapply(read_panel("shmi_vs_sh"), enlarge)
fig <- sh$all | sh$between | sh$within
save_fig(tag_fig(fig), "F2_system_score.png", width = 12.5, height = 4.6)
