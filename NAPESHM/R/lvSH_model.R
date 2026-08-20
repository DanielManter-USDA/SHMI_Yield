############################################################
# lvSH_model.R — SEM for soil health latent variable
############################################################

run_lvSH_model <- function(finalDF) {
  
  model_soil_health <- '
    SH.Bio     =~ Cmin + Nmin + ARYL + NAG
    SH.CN      =~ SOC + TN + POXC + ACE
    SH.Phys    =~ SSI + WAS + AWHC + Ksat
    SH         =~ SH.Bio + SH.CN + SH.Phys
  '
  
  fit <- sem(model_soil_health, data = finalDF,
             std.lv = TRUE, missing = "fiml")
  
  lv_scores <- lavPredict(fit, type = "lv")
  idx <- unlist(fit@Data@case.idx)
  
  latent_names <- c("SH.Bio", "SH.CN", "SH.Phys", "SH")
  for (nm in latent_names) {
    finalDF[[nm]] <- NA
    finalDF[[nm]][idx] <- lv_scores[, nm]
    finalDF[[nm]] <- minmax_scale(finalDF[[nm]])
  }
  
  saveRDS(finalDF %>% select(MGT_combo, SH, SH.Bio, SH.CN, SH.Phys),
          "sh_scores.rds")
  
  return(list(fit = fit, finalDF = finalDF))
}

############################################################
# End lvSH_model.R
############################################################
