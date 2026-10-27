# coarse_records

`degrade_records.py` writes a copy of an SHMI workbook with coarser management records (typical dates, cover crop yes/no, tillage class, main crop only, no amendments or grazing, a short history, date noise, misclassification). `R/B9_SHMI_coarse_records.R` runs the scenarios (step B10 of `run_all.Rmd`) and rescores them with the published weights.

`tdisc_implements.csv` and `tdisc_mapping.csv` are the T-DISC implement tables shipped with the SHMI package (USDA Office of the Chief Economist, T-DISC 1.1.1), used to classify each crop interval's tillage.

Requires python 3.9+ with openpyxl.
