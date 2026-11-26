# GPR-and-MP

Code for the master's dissertation *Geopolitical risk and the transmission of
ECB monetary policy* (smooth-transition BVAR with Jarociński–Karadi (2020)
identification, plus robustness checks and state-dependent local projections).

## Running the whole project

```bash
Rscript run_pipeline.R
```

This runs every stage end to end and writes the results to `results/`,
grouped by subject. Tables are LaTeX (`.tex`), figures are PNG (300 dpi) and PDF.

| Folder | Content | Thesis |
|---|---|---|
| `results/descriptive_statistics/` | summary statistics, GPR-regime tables, correlation matrix, data plots | Tables 1–4, A6; Figures B1–B6 |
| `results/unit_root/` | ADF / Lee-Strazicich tests, Bai-Perron breaks, series and break plots, full test output | Tables A2–A3 |
| `results/shock_diagnostics/` | ADF / ARCH / exogeneity tests, Ljung-Box tables, shock correlations | Tables A4–A5; Figure A1 |
| `results/VAR/` | linear BVAR, baseline ST-BVAR, all robustness ST-BVARs, Cholesky ST-BVAR (+ IRF data as CSV) | Figures 1, 2, C1–C21 |
| `results/LP/` | state-dependent local projection (+ coefficient table) | Figure 3 |
| `results/logs/` | console log of every stage, timings, `sessionInfo()` | |

### Options

```bash
Rscript run_pipeline.R descriptive unit_root          # run only some stages
Rscript run_pipeline.R VAR --models=Figure1_STBVAR_Baseline,FigureC16_STBVAR_Threshold_VIX
Rscript run_pipeline.R --quick                        # 200/200 Gibbs draws, for testing
Rscript run_pipeline.R --cores=2                      # cores for unit root tests and BVARs
Rscript run_pipeline.R --burnin=2000 --ndraws=2000    # Gibbs sampler settings
Rscript run_pipeline.R --seed=123                     # base RNG seed
Rscript run_pipeline.R --formats=png                  # figure formats (png,pdf)
Rscript run_pipeline.R shocks --exo-lags=2            # 2-lag exogeneity regressions
Rscript run_pipeline.R --download                     # rebuild the data from the web
```

Stages: `descriptive`, `unit_root`, `shocks`, `VAR`, `LP`. The LP stage uses
the structural shock of the linear BVAR; when it is run on its own it reads the
shock saved by the last VAR run.

### Data

By default the pipeline reads `data_thesis_r0838656.csv`, the merged dataset
(`final_data`) used in the thesis. This keeps the results reproducible and
needs no internet access. `--download` rebuilds it with `Merge.R` (Eurostat,
ECB Data Portal, FRED, EA-MPD via `hfdshocks`); later data revisions can then
change the results.

### Reproducibility

Every BVAR has its own fixed seed (base `--seed`, default `20250101`), so the
Bayesian results are identical between runs and do not depend on the number of
cores. Gibbs sampler defaults follow the thesis: 4000 burn-in draws and 4000
saved draws, keeping every 4th draw.

### Requirements

R (≥ 4.3) with: dplyr, tidyr, readr, lubridate, ggplot2, patchwork, sandwich,
urca, strucchange, dynlm, FinTS, moments, xtable, knitr, kableExtra, gt,
stargazer, corrplot. `--download` additionally needs the packages loaded in
`Merge.R` (eurostat, ecb, readxl, countrycode, hfdshocks, ...).

## Structure

- `run_pipeline.R` – entry point and configuration
- `pipeline/` – the stages (`00_setup.R`, `01_data.R`, `02_…` – `06_lp.R`)
- function files used by the pipeline: `TS_functions.R`, `additional_functions.R`,
  `LeeStrazicichUnitRootTest.R`, `jk_replication_function.R`,
  `jk_SD_BVAR_functions.R`, `ST_LP_functions.R`
- data scripts: `Merge.R` and the scripts it sources
- the original interactive analysis scripts (`*_main.R`, `Desc_Stat_and_Stylized_Facts.R`,
  `Unit_root_testing_*.R`, `shock_analysis_and_purging.R`) are kept unchanged
  for reference.
