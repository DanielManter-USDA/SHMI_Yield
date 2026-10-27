# R packages used by the analysis scripts, run_all.R and report.Rmd
pkgs <- c("callr", "here", "rmarkdown", "knitr", "dplyr", "tidyr", "tibble", "purrr", "readxl", "writexl",
          "ggplot2", "patchwork", "ggh4x", "scales", "magick", "rlang", "nlme", "lme4", "lmerTest", "metafor",
          "mgcv", "MASS", "nnls", "ranger", "lavaan", "semTools", "terra", "sf", "geodata", "prism", "daymetr",
          "soilDB", "DBI", "RSQLite", "rnassqs")
install.packages(setdiff(pkgs, rownames(installed.packages())))
# SHMI (version used in the paper: 2.0.0)
# install.packages("remotes"); remotes::install_github("DanielManter-USDA/SHMI")
# python (3.9+): pip install pandas openpyxl numpy
