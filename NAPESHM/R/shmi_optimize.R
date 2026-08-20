############################################################
# shmi_optimize.R — SHMI weighting, tuning, Bayesian optimization
############################################################

# This file contains:
#   - build_indices()
#   - objective_fun()
#   - Bayesian optimization (mlrMBO)
#   - build_shmi()

############################################################
# 1. Build indices for a given parameter set
############################################################

build_indices <- function(theta, sh_scores, crop_harmonized, rot_bounds, daily_dist, daily, amend, animal) {
  
  # 1. Cover
  cover_df <- compute_w.cover(
    daily      = daily,
    rot_bounds = rot_bounds,
    w_winter   = theta$w_winter,
    w_spring   = theta$w_spring,
    w_summer   = theta$w_summer,
    w_fall     = theta$w_fall
  )
  
  # 2. Diversity
  div_df <- compute_rot_diversity(
    crop_harmonized = crop_harmonized,
    daily           = daily,
    hill            = theta$hill,
    max_div         = 8
  )
  
  # 3. Disturbance
  dist_df <- compute_avg_annual_disturbance(
    daily_dist  = daily_dist,
    rot_bounds  = rot_bounds,
    transform   = "raw",
    max_stir    = theta$max_stir
  )
  
  # 4. Organic Inputs
  org_df <- compute_orginput(
    rot_bounds  = rot_bounds,
    amend       = amend,
    animal      = animal,
    w_amend     = 1,
    w_animal    = 1
  )
  
  # 5. Join everything to SH scores
  full <- sh_scores %>%
    left_join(cover_df, by = "MGT_combo") %>%
    left_join(div_df,   by = "MGT_combo") %>%
    left_join(dist_df,  by = "MGT_combo") %>%
    left_join(org_df,   by = "MGT_combo")
  
  return(full)
}


############################################################
# 2. Objective function for Bayesian optimization
############################################################

objective_fun<- function(
    w_winter,
    w_spring,
    w_summer,
    w_fall,
    hill,
    max_stir,
    w_cover,
    w_diversity,
    w_invdist,
    w_orginput,
    sh_scores
) {
  
  theta <- list(
    w_winter = w_winter,
    w_spring = w_spring,
    w_summer = w_summer,
    w_fall   = w_fall,
    hill     = hill,
    max_stir = max_stir,
    w_cover     = w_cover,
    w_diversity = w_diversity,
    w_invdist   = w_invdist,
    w_orginput  = w_orginput
  )
  
  out <- tryCatch({

    dat <- build_indices(
      theta = theta,
      sh_scores = sh_scores,
      crop_harmonized = crop_harmonized,
      rot_bounds = rot_bounds,
      daily_dist = daily_dist,
      daily = daily,
      amend = amend,
      animal = animal
    )
    
    df <- dat %>%
      dplyr::select(SH, Cover, Diversity, InvDist, Animals) %>%
      na.omit()
    
    # Normalize SHMI pillar weights
    w_sum <- w_cover + w_diversity + w_invdist + w_orginput
    w_cover     <- w_cover     / w_sum
    w_diversity <- w_diversity / w_sum
    w_invdist   <- w_invdist   / w_sum
    w_orginput  <- w_orginput  / w_sum
    
    df$SHMI <- w_cover     * df$Cover +
      w_diversity * df$Diversity +
      w_invdist   * df$InvDist +
      w_orginput  * df$Animals
    
    # Simple regression
    fit <- lm(SH ~ SHMI, data = df)
    
    preds <- predict(fit)
    y <- df$SH
    
    mae <- mean(abs(y - preds))
    
    ss_res <- sum((y - preds)^2)
    ss_tot <- sum((y - mean(y))^2)
    r2 <- 1 - ss_res / ss_tot
    
    w_sum <- w_winter + w_spring + w_summer + w_fall
    w_winter <- w_winter / w_sum
    w_spring <- w_spring / w_sum
    w_summer <- w_summer / w_sum
    w_fall   <- w_fall / w_sum
    
    log_row <- data.frame(
      w_winter = w_winter,
      w_spring = w_spring,
      w_summer = w_summer,
      w_fall   = w_fall,
      hill     = hill,
      max_stir = max_stir,
      w_cover     = w_cover,
      w_diversity = w_diversity,
      w_invdist   = w_invdist,
      w_orginput  = w_orginput,
      R2          = r2,
      MAE         = mae
    )
    
    write.table(
      log_row,
      file = "bo_log.final.csv",
      sep = ",",
      row.names = FALSE,
      col.names = FALSE,
      append = TRUE
    )
    
    list(mae = mae, r2 = r2)
    
  }, error = function(e) {
    list(mae = 1e6, r2 = -1e6)
  })
  
  list(
    Score = -out$mae,
    MAE   = out$mae,
    R2    = out$r2
  )
}

############################################################
# 3. Parameter space definition (ps <- makeParamSet(...))
############################################################

ps <- makeParamSet(
  # cover settings
  makeNumericParam("w_winter", lower = 0.1, upper = 1),
  makeNumericParam("w_spring", lower = 0.1, upper = 1),
  makeNumericParam("w_summer", lower = 0.1, upper = 1),
  makeNumericParam("w_fall",   lower = 0.1, upper = 1),
  
  # diversity settings
  makeIntegerParam("hill", lower = 1, upper = 2),
  
  # disturbance settings
  makeIntegerParam("max_stir", lower = 80, upper = 300),
  
  # SHMI settings
  makeNumericParam("w_cover",     lower = 0.1, upper = 1),
  makeNumericParam("w_diversity", lower = 0.1, upper = 1),
  makeNumericParam("w_invdist",   lower = 0.1, upper = 1),
  makeNumericParam("w_orginput",  lower = 0.1, upper = 1)
)


############################################################
# 4. Run Bayesian optimization
############################################################

run_mbo <- function(sh_scores) {
  write.table(
    data.frame(
      w_winter       = numeric(),
      w_spring       = numeric(),
      w_summer       = numeric(),
      w_fall         = numeric(),
      hill           = integer(),
      max_stir       = integer(),
      w_cover        = numeric(),
      w_diversity    = numeric(),
      w_invdist      = numeric(),
      w_orginput     = numeric(),
      R2             = numeric(),
      MAE            = numeric()
    ),
    file = "NAPESHM/data/bo_log.final.csv",
    sep = ",",
    row.names = FALSE,
    col.names = TRUE
  )
  
  lrn <- makeLearner(
    "regr.km",
    predict.type = "se",
    covtype = "matern3_2",
    nugget.estim = TRUE
  )
  
  ctrl <- makeMBOControl()
  ctrl <- setMBOControlTermination(ctrl, iters = 500)
  
  ctrl <- setMBOControlInfill(
    ctrl,
    crit = makeMBOInfillCritCB(cb.lambda = 1),
    opt.focussearch.points = 500,
    opt.focussearch.maxit = 5        # improves local refinement
  )
  
  design <- generateDesign(500, ps)   # improves early stability
  
  obj_for_mbo <- makeSingleObjectiveFunction(
    name = "shmi_tuning",
    fn = function(x) {
      
      out <- tryCatch(
        {
          args <- c(as.list(x), list(sh_scores = sh_scores))
          do.call(objective_fun, args)
        },
        error = function(e) {
          list(MAE = 1e6, R2 = -1e6)
        }
      )
      
      mae <- out$MAE
      
      # Absolute guarantee: never return NA/NaN/Inf
      if (!is.finite(mae) || is.na(mae)) {
        mae <- 1e6
      }
      
      mae
    },
    par.set = ps,
    minimize = TRUE
  )
  
  parallelStartSocket(parallel::detectCores() - 1)
  
  mbo_result <- mbo(
    fun = obj_for_mbo,
    design = design,
    learner = lrn,
    control = ctrl,
    show.info = TRUE
  )
  
  parallelStop()
  saveRDS(mbo_result, "NAPESHM/data/mbo_result.rds")
}


############################################################
# 5. Extract top models, consensus parameters, CI bands
############################################################

extract_mbo_results <- function(path = "NAPESHM/data/mbo_result.rds") {
  
  # 1. Load results
  mbo_result <- readRDS(path)
  
  # 2. Normalize weights
  res <- as.data.frame(mbo_result$opt.path) %>%
    filter(prop.type != "initdesign")
  
  bo_res <- res %>%
    arrange(y) %>%
    dplyr::select(w_winter, w_spring, w_summer, w_fall,
           hill,
           max_stir,
           w_cover, w_diversity, w_invdist, w_orginput,
           y) %>%
    mutate(
      w_sum = w_winter + w_spring + w_summer + w_fall,
      w_winter = w_winter / w_sum,
      w_spring = w_spring / w_sum,
      w_summer = w_summer / w_sum,
      w_fall   = w_fall / w_sum
    ) %>%
    mutate(
      w_sum = w_cover + w_diversity + w_invdist + w_orginput,
      w_cover     = w_cover     / w_sum,
      w_diversity = w_diversity / w_sum,
      w_invdist   = w_invdist   / w_sum,
      w_orginput  = w_orginput  / w_sum
    ) %>%
    dplyr::select(-w_sum) %>%
    rename(MAE = y)
  
  # 3. Top 100 models
  top100 <- bo_res %>%
    arrange(MAE) %>%
    slice(1:100)
  
  # 4. Consensus parameter set
  global_consensus <- top100 %>%
    dplyr::select(where(is.numeric)) %>%
    summarise(across(everything(), median))
  
  # 5. Return everything cleanly
  list(
    bo_res = bo_res,
    top100 = top100,
    global_consensus = global_consensus
  )
}

############################################################
# 6. Final SHMI assembly
############################################################

build_shmi <- function(shmi_inputs,
                       settings = NULL,
                       expert_mode = TRUE) {
  
  # --------------------------------------------------------------------------
  # 1. Official national SHMI settings (locked mode)
  # --------------------------------------------------------------------------
  official <- list(
    # cover
    w_winter = 0.130,
    w_spring = 0.129,
    w_summer = 0.513,
    w_fall   = 0.227,
    
    # diversity
    hill      = 2,
    max_div   = 8,
    
    # disturbance
    transform   = "raw",
    max_stir    = 205,
    
    # animals
    w_amend   = 1,
    w_animals = 1,
    
    # shmi weights
    w_cover    = 0.492,
    w_div      = 0.052,
    w_dist     = 0.131,
    w_ani      = 0.324
  )
  
  # --------------------------------------------------------------------------
  # 2. Determine which settings to use
  # --------------------------------------------------------------------------
  if (!expert_mode) {
    # Locked mode: always use official settings
    if (!is.null(settings)) {
      message("Note: Custom settings ignored because expert_mode = FALSE. ",
              "Using official national SHMI settings.")
    }
    settings <- official
  } else {
    # Expert mode: allow tuning but warn
    message("Expert mode enabled: SHMI scores will NOT be comparable ",
            "to the national SHMI scale.")
    # Fill missing settings with official defaults
    settings <- modifyList(official, settings)
  }
  
  # --------------------------------------------------------------------------
  # 3. Extract inputs
  # --------------------------------------------------------------------------
  rot_bounds      <- shmi_inputs$rot_bounds
  crop_harmonized <- shmi_inputs$crop_harmonized
  daily           <- shmi_inputs$daily
  daily_dist      <- shmi_inputs$daily_dist
  amend           <- shmi_inputs$amend
  animal          <- shmi_inputs$animal
  
  required <- c("rot_bounds", "crop_harmonized", "daily", "daily_dist",
                "amend", "animal")
  
  missing <- setdiff(required, names(shmi_inputs))
  if (length(missing) > 0) {
    stop("Missing required inputs in shmi_inputs: ", paste(missing, collapse = ", "))
  }
  
  # --------------------------------------------------------------------------
  # 4. Compute pillars
  # --------------------------------------------------------------------------
  
  # Cover
  cover <- compute_w.cover(
    daily       = daily,
    rot_bounds  = rot_bounds,
    w_winter    = settings$w_winter,
    w_spring    = settings$w_spring,
    w_summer    = settings$w_summer,
    w_fall      = settings$w_fall
  )
  
  # Diversity
  diversity <- compute_rot_diversity(
    crop_harmonized = crop_harmonized,
    daily           = daily,
    hill            = settings$hill,
    max_div         = settings$max_div
  )
  
  # Disturbance
  invdist <- compute_avg_annual_disturbance(
    daily_dist    = daily_dist,
    rot_bounds    = rot_bounds,
    transform     = settings$transform,
    max_stir      = settings$max_stir
  )
  
  orginput <- compute_orginput(
    rot_bounds  = rot_bounds,
    amend       = amend,
    animal      = animal,
    w_amend     = settings$w_amend,
    w_animal    = settings$w_animal
  )
  
  # --------------------------------------------------------------------------
  # 5. Combine pillars
  # --------------------------------------------------------------------------
  indicator_df <- reduce(
    list(cover, diversity, invdist, orginput),
    full_join,
    by = "MGT_combo"
  )
  
  w_sum   <- settings$w_cover + settings$w_div + settings$w_dist + settings$w_ani
  
  w_cover <- settings$w_cover / w_sum
  w_div   <- settings$w_div   / w_sum
  w_dist  <- settings$w_dist  / w_sum
  w_ani   <- settings$w_ani   / w_sum
  
  indicator_df <- indicator_df %>%
    mutate(
      SHMI = (
        w_cover * Cover +
          w_div   * Diversity +
          w_dist  * InvDist +
          w_ani   * Animals
      )
    ) %>%
    dplyr::select(MGT_combo, SHMI, Cover, Diversity, InvDist, Animals) %>%
    arrange(MGT_combo)
  
  if (!expert_mode) {
    message("SHMI computed using official national settings.")
  }
  
  # --------------------------------------------------------------------------
  # 6. Return both indicators and settings used
  # --------------------------------------------------------------------------
  return(list(
    indicator_df = indicator_df,
    settings_used = settings,
    expert_mode = expert_mode,
    shmi_version = "1.0.0",
    timestamp = Sys.time()
  ))
  
}

# ###########################################################
# # 7. Attach SH latent variables and evaluate SHMI
# ############################################################
# 
# shmi_out <- build_shmi(
#   shmi_inputs  = shmi_inputs,
#   settings     = NULL,
#   expert_mode  = FALSE
# )
# 
# # attach Soil Health latent variables
# sh_scores <- readRDS("NAPESHM/data/sh_scores.rds")
# 
# shmi_full <- shmi_out$indicator_df %>%
#   left_join(sh_scores, by = "MGT_combo")
# 
# saveRDS(shmi_out$indicator_df, "NAPESHM/data/optimum_shmi_scores.rds")


############################################################
# End shmi_optimize.R
############################################################
