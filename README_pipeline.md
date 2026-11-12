# Running the SHMI_Yield pipeline with `targets`

The pipeline in `_targets.R` replaces `run_all.R` for day-to-day runs. It knows which step depends on which, so it reruns only what is out of date.

## One-time setup

```r
install.packages(c("targets", "visNetwork"))          # visNetwork draws the graph
targets::tar_config_set(reporter_make = "timestamp")  # plain progress lines (writes _targets.yaml)
```

The second line matters. The default live progress bar can be cleared by the scripts' own progress messages, because they run inside the pipeline's R process. When that happens, `tar_make()` ends with "Cannot find progress bar".

## Everyday use

| You want to | Run |
|---|---|
| bring everything up to date (only stale steps run) | `targets::tar_make()` |
| see the pipeline and what is stale | `targets::tar_visnetwork()` |
| list what `tar_make()` would run, without running it | `targets::tar_outdated()` |
| get the files a step wrote | `targets::tar_read(B2_main)` |
| force one step to rerun | `targets::tar_invalidate(B3_stability)`, then `tar_make()` |
| run one step and what it needs | `targets::tar_make(F4_moderators)` |

Each step runs one script from `R/`. Its log goes to `Outputs/logs/<date-time>/<script>.log`, as with `run_all.R`.

## What triggers a rerun

- **The script was edited.** For example, editing `R/B3_SHMI_yield_stability.R` reruns B3 and then whatever uses its files.
- **An input file changed.** This covers the workbooks in `Inputs/yield_workbooks`, the NAPESHM files, the SPEI file, the caches in `data/`, and the frozen Part A outputs.
- **An upstream step's files changed.** For example, a new B2 run reruns B2b, B3, B4 and the others that use it.
- **The SHMI package was reinstalled.** Every scoring step reruns, as B2's site cache already does after a reinstall.
- **Editing figure scripts or `fig_common.R` in `R/figures/` redraws only the figures.** This takes seconds.

## Special cases

- **Accepting work that already finished.** If `tar_make()` stops after a step has written its files (for example, because of the progress-bar error), register that step without rerunning it:

  ```r
  Sys.setenv(SHMI_KEEP = "B5b")   # script prefixes, comma-separated: "B4,B5b"
  targets::tar_make()
  Sys.unsetenv("SHMI_KEEP")
  ```

- **Part A (the calibration) is frozen.** Its outputs are tracked as inputs. To recalibrate, use `run_all.R` with `run_calibration = TRUE`.
- **B9 (the coarse-record test) takes hours.** It reruns only when forced: `Sys.setenv(SHMI_FORCE = "B9"); targets::tar_make()`. Otherwise its existing outputs are kept.
- **Retired scripts are not run:**
  - B2e (talk figure);
  - B2g, B2h (practice pairs for yield);
  - B5 (mediation) and B0, which builds the crosswalk only B5 uses;
  - C1b;
  - the A5b and A6c redraws.

## Speed settings (set before `tar_make()`)

- **`B2D_WORKERS`** sets the number of R processes B2d uses for its model fits. The default is your cores minus one, at most 6. Use `B2D_WORKERS <- 1` to fit one after another; the results are the same.
- **`RL_FULL <- TRUE`** makes B2c also run its record-length analysis. By default B2c only makes the rolling scores that B2d uses; the analysis takes about an hour and is not in the paper.
- **B4 reads B2's saved model** (`B2_model.rds`) instead of rerunning B2.

## Adding a step

1. Add the script to the scripts block with `tar_target(s_X, sc("X.R"), format = "file")`.
2. Add a step:

   ```r
   tar_target(X_name, step("X.R", c("<files it writes that others read>"), s_X, <upstream targets>), format = "file")
   ```

3. List the new target among the upstream targets of any step that reads its files.

## Adding a figure

1. Have the analysis script save what the figure draws as `plotdata_<name>.rds` in its output folder.
2. Add the path to `PLOTDATA` in `R/figures/fig_common.R`.
3. Write `R/figures/F<n>_<name>.R`, which reads the plot data with `read_pd()` and saves with `save_fig()`.
4. Add a `figure()` target in `_targets.R`, and the script to `s_fig`.
