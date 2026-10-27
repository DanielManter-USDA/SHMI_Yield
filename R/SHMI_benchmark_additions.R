################################################################################
# Additions to SHMI_benchmark.R: a row for SHMI with the published weights held
# fixed (what a user applies). Only the score-to-soil-health line is fitted in
# each fold; the weights were fitted on all US sites, so this row is slightly
# optimistic and is reported in the caption, not plotted in Fig. 2d.
################################################################################

# 1. Paste after `equal_model <- ...` (around line 146):
fixed_model <- line_model(function(D) {        # published weights, not refitted
  base <- W_FINAL[["cover"]] * D$Cover_cal + W_FINAL[["diversity"]] * D$Diversity +
          W_FINAL[["disturbance"]] * D$InvDist
  base + CFG$amend_gap * D$p_amend * (100 - base)
})

# 2. In `cand <- list(...)`, add as the second row (SHMI must stay first):
#   list("SHMI, published weights (fixed)", "SHMI", fixed_model),

# 3. Optional check, paste after the "Check: SHMI OOF R2" block:
# The fixed-weight score should equal the published SHMI exactly.
if ("SHMI" %in% names(us)) {
  s_fixed <- environment(fixed_model)$score_fun(us)
  cat(sprintf("Fixed-weight score vs published SHMI: max abs difference %.2e\n",
              max(abs(s_fixed - us$SHMI))))
}
