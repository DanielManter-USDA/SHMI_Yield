################################################################################
# Validate_converter_SHMI.R
# Do workbooks built by agcros2shmi.py give the same SHMI as the hand-built ones?
# For each validation site, SHMI is computed from both workbooks over the same years
# (the hand-built unit's record period), and compared field by field: each field's
# main unit (longest record) in each workbook.
# Put the hand-built workbooks in HAND_DIR and the converted ones in CONV_DIR.
################################################################################
library(SHMI); library(dplyr); library(ggplot2)
setwd(here::here())
HAND_DIR <- "Inputs/handbuilt_workbooks"                 # SHMI_<site>.xlsx, hand-built
CONV_DIR <- "Inputs/yield_workbooks"  # SHMI_<site>.xlsx from agcros2shmi.py
SITES    <- c("NELITCSE", "NEMLTCRS", "NEMEIRR", "NECCIRR")
OUT      <- "Outputs/agcros_validation"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

score <- function(path, start = NULL, end = NULL) {
  inp <- prepare_shmi_inputs(path, verbose = FALSE)
  if (!is.null(start)) {
    rb <- inp$rot_bounds
    rb$rot_start <- pmax(as.Date(rb$rot_start), start); rb$rot_end <- pmin(as.Date(rb$rot_end), end)
    rb <- rb[rb$rot_start <= rb$rot_end, ]
    rb$rot_start_yr <- as.integer(format(rb$rot_start, "%Y")); rb$rot_end_yr <- as.integer(format(rb$rot_end, "%Y"))
    inp$rot_bounds <- rb; inp$mgt <- inp$mgt[inp$mgt$MGT_combo %in% rb$MGT_combo, ]
  }
  cl <- get_shmi_climate(inp)
  build_shmi(inp, cl)$indicator_df %>%
    left_join(inp$mgt %>% select(MGT_combo, MGT_field), by = "MGT_combo") %>%
    left_join(inp$rot_bounds %>% select(MGT_combo, rot_start, rot_end), by = "MGT_combo") %>%
    mutate(field = sub("^0+", "", as.character(MGT_field)), days = as.numeric(rot_end - rot_start))
}

res <- bind_rows(lapply(SITES, function(s) {
  h <- score(file.path(HAND_DIR, paste0("SHMI_", s, ".xlsx")))
  per <- c(min(h$rot_start), max(h$rot_end))
  a <- score(file.path(CONV_DIR, paste0("SHMI_", s, ".xlsx")), per[1], per[2])
  main <- function(d) d %>% group_by(field) %>% slice_max(days, n = 1, with_ties = FALSE) %>% ungroup()
  inner_join(main(h) %>% select(field, SHMI_hand = SHMI, Cover_hand = Cover, InvDist_hand = InvDist, Diversity_hand = Diversity),
             main(a) %>% select(field, SHMI_conv = SHMI, Cover_conv = Cover, InvDist_conv = InvDist, Diversity_conv = Diversity),
             by = "field") %>% mutate(site = s)
}))
write.csv(res, file.path(OUT, "converter_vs_handbuilt_SHMI.csv"), row.names = FALSE)
print(res %>% group_by(site) %>%
        summarise(fields = n(), r = cor(SHMI_hand, SHMI_conv), mean_abs_diff = mean(abs(SHMI_hand - SHMI_conv)),
                  max_abs_diff = max(abs(SHMI_hand - SHMI_conv)), .groups = "drop"))
p <- ggplot(res, aes(SHMI_hand, SHMI_conv, colour = site)) + geom_abline(linetype = 2) + geom_point() +
  labs(x = "SHMI, hand-built workbook", y = "SHMI, converted from AgCROS") + theme_bw()
ggsave(file.path(OUT, "converter_vs_handbuilt_SHMI.png"), p, width = 5.5, height = 4.5, dpi = 300)
