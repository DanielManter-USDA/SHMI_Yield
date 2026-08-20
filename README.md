# SHMI: Soil Health Management Index

This repository contains all code, data, and documentation required to reproduce the analyses in:

**From practice to performance: Linking soil health management to yield**

The project is organized into two major experimental datasets:

- **NAPESHM** — soil health and management experiment  
- **DRIVES** — long‑term rotation, disturbance, and input experiment  

Each dataset contains its own data, analysis code, scripts, and figures to ensure clean separation and reproducibility.

---

## 🔧 How to Run the Pipeline

All NAPESHM analyses — data preparation, latent soil health modeling, SHMI construction, Bayesian optimization, and manuscript figures — are executed from a single script:

### **`NAPESHM/R/run_all.R`**

To reproduce the full workflow:

1. Open the R project (`SHMI_Yield/`)
2. Open `NAPESHM/R/run_all.R`
3. Click **Source** (or run `source("NAPESHM/R/run_all.R")`)

This script automatically:

- loads all required functions via `setup.R`
- prepares and cleans raw NAPESHM data  
- fits the Soil Health latent variable SEM  
- optionally runs Bayesian optimization  
- builds final SHMI scores  
- merges SHMI with soil health + yield data  
- generates all manuscript figures (Fig. 2a–c)

No other scripts need to be run manually.

All DRIVES analyses — data preparation, SHMI construction, yield modeling, and manuscript figures — are executed from a single script:

### **`NAPESHM/R/run_all.R`**

To reproduce the full workflow:

1. Open the R project (`SHMI_Yield/`)
2. Open `DRIVES/R/run_all.R`
3. Click **Source** (or run `source("DRIVES/R/run_all.R")`)

This script automatically:

- loads all required functions
- prepares and cleans raw DRIVES data  
- builds final SHMI scores 
- runs all yield and yield stability models
- generates all manuscript figures (Fig. 3a–c & Fig. 4a-c)

No other scripts need to be run manually.

---

## 📦 Contents

### **NAPESHM/**
Contains all data and analysis related to the NAPESHM long‑term experiment:

- `data/` — raw and processed datasets  
- `R/` — functions for lv_SH modeling, SHMI computation, optimization, yield modeling, and stability metrics  
- `figs/` — figures generated for the manuscript  

### **DRIVES/**
Contains all data and analysis related to the DRIVES long‑term experiment:

- `data/` — raw and processed datasets  
- `R/` — functions for SHMI computation, yield modeling, and stability analyses  
- `figs/` — figures generated for the manuscript  

---

## 📚 Citation

If you use SHMI or any part of this repository, please cite:

> **[add manuscript citation here]**

---

## 📬 Contact

For questions or collaboration, contact:  
**Daniel Manter**  
USDA‑ARS  
daniel.manter@usda.gov

