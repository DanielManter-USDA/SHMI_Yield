############################################################
# setup.R — Project setup and package loading
############################################################

# Core tidyverse + utilities
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(readxl)
library(lubridate)
library(broom)

# Visualization
library(ggplot2)
library(ggh4x)
library(ggpmisc)
library(patchwork)
library(pheatmap)
library(corrplot)

# Modeling
library(lavaan)
library(nlme)
library(emmeans)
library(effectsize)
library(matrixStats)
library(quantreg)
library(glmnet)
library(scam)
library(CCA)

# SHMI Bayesian optimization
library(mlrMBO)
library(mlr)
library(ParamHelpers)
library(ranger)
library(parallelMap)
library(DiceKriging)

# Misc
library(janitor)   # used in shmi_index.R
library(scales)    # used for rescale()


# Project root
project_root <- dirname(dirname(getwd()))

# Source all module scripts
source("NAPESHM/R/utils.R")
source("NAPESHM/R/data_prep.R")
source("NAPESHM/R/lvSH_model.R")
source("NAPESHM/R/shmi_index.R")
source("NAPESHM/R/shmi_optimize.R")

############################################################
# End setup.R
############################################################
