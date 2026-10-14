# Fish population synchrony in a tropical reservoir

Data and R code for the manuscript:

> **Integrating asymmetry, geography and timescale analyses to uncover fish population synchrony patterns in a tropical reservoir** Emanuel R. Cararo *et al.* --- *[Journal, year, DOI: add after publication]*

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23354278.svg)](https://doi.org/10.5281/zenodo.23354278)

## Overview

The workflow quantifies spatial synchrony of fish populations sampled at ten sites of the Madeira River reservoir (Brazil). It:

1.  builds species-by-site CPUE time series from quarterly samples and aggregates them to hydrological years (2013--2022);
2.  estimates mean spatial synchrony (Spearman cross-correlation) with bootstrap confidence intervals;
3.  tests the asymmetry of synchrony (upper- vs. lower-tail dependence) against copula-based surrogate datasets;
4.  explains pairwise synchrony with density dependence, dispersal difficulty (distance and dam barrier) and water-level synchrony using Multiple Regression on distance Matrices (MRM);
5.  relates species mean synchrony to traits, body length and trophic level with linear and phylogenetic generalised least squares (PGLS) models;
6.  describes timescale-specific synchrony with wavelet phasor mean fields and identifies synchrony modules (networks) tested against a null model.

## Repository structure

```         
.
├── run_all.R                 # runs the whole workflow in order
├── Scripts/
│   ├── 1_FUNCTIONS.R         # packages and helper functions
│   ├── 2_DATA.R              # reads and prepares all inputs
│   ├── 3_ANALYSIS.R          # all analyses
│   └── 4_VISUALIZATIONS.R    # all figures
├── Data/
│   ├── 1_CPUE.xlsx           # CPUE per site, year and month
│   ├── 2_LOCALS.csv          # sampling sites
│   ├── 3_WATER_LEVEL.csv     # water level per site and year
│   ├── 4_TRAITS.xlsx         # species traits
│   ├── fish_tree.tre         # phylogeny (Newick)
│   ├── mad_shape.rds         # Madeira River outline for maps
│   └── LICENSE               # CC BY 4.0 (data)
├── Outputs/                  # created by the scripts
├── Plots/                    # created by the scripts
├── LICENSE                   # MIT (code)
├── CITATION.cff
└── README.md
```

## Requirements

-   R ≥ 4.2 (tested with R 4.3.1 on Windows 11; see `Outputs/sessionInfo.txt` for exact package versions).
-   On Windows, Rtools may be needed to compile some packages.

Install the R packages once:

``` r
pkgs <- c(
  "tidyverse", "readxl", "wsyn", "ecodist", "lm.beta", "phytools", "ade4",
  "ape", "nlme", "psych", "ggplot2", "ggrepel", "sf", "ggeasy", "ggstatsplot",
  "patchwork", "ncf", "pracma", "copula", "matrixcalc", "Matrix", "mvtnorm",
  "stringr", "scales"
)
install.packages(setdiff(pkgs, rownames(installed.packages())))
```

## How to run

Set the working directory to the repository root (the folder containing `Scripts/` and `Data/`), then either:

``` r
source("run_all.R")
```

or, from a terminal:

``` bash
Rscript run_all.R
```

The scripts must run in order and in the same R session, because each one uses objects created by the previous one (`run_all.R` does this for you). `1_FUNCTIONS.R` clears the workspace and sets the random seed (`set.seed(263)`), so the bootstrap, surrogate and permutation results are reproducible when the whole pipeline is run from the start.

Approximate run time on a laptop: `1_FUNCTIONS.R` and `2_DATA.R` take seconds; `3_ANALYSIS.R` is the slow step (copula surrogates and modularity null models with 1000 replicates per species), from tens of minutes to a few hours depending on the machine; `4_VISUALIZATIONS.R` takes about a minute.

## Outputs

Tables (`Outputs/`; generated files and figures are not version-controlled, see `.gitignore`):

| File                                  | Content                                                                                                                                           |
|------------------------------------|------------------------------------|
| `mean_spatial_synchrony.xlsx`         | Mean Spearman synchrony per species with 95% bootstrap CI                                                                                         |
| `tail_synchrony_summary.xlsx`         | Mean total, lower-tail and upper-tail synchrony and tail association per species                                                                  |
| `tail_dependence_surrogate_test.xlsx` | Surrogate test of tail asymmetry (empirical sum, 2.5/97.5% surrogate quantiles, rank-based p)                                                     |
| `mrm_results.xlsx`                    | MRM standardised coefficients (± SE) for density dependence, dispersal difficulty and water-level synchrony, MRM F-test p-value (`f_stat`) and R² |
| `data_pgls.xlsx`                      | Species-level data used in the trait and PGLS models                                                                                              |
| `trait_models_summary.txt`            | PCoA eigenvalues, linear model, PGLS (Brownian, Pagel's λ, OU) summaries and AIC                                                                  |
| `wpmf_results.rds`                    | Wavelet phasor mean field objects (`wsyn::wpmf`) per species                                                                                      |
| `modularity_null_model.xlsx`          | Number of modules, observed modularity and null-model statistics per species                                                                      |
| `sessionInfo.txt`                     | R and package versions used in the last run                                                                                                       |

Figures (`Plots/`):

| File                                        | Content                                         |
|------------------------------------|------------------------------------|
| `Forest_plot_mean_synchrony.tiff`           | Mean spatial synchrony per species              |
| `Forest_plot_coefficients.tiff`             | MRM coefficients and R² per species             |
| `PCoA_plots.tiff`                           | PCoA of species traits                          |
| `Supplementary_LinearModel_Diagnostics.png` | Diagnostic plots of the trait linear model      |
| `all_wavelets_grouped.pdf`                  | Wavelet phasor mean field magnitude per species |
| `Synchrony_networks.pdf`                    | Map of sites and synchrony modules per species  |

## Input data

**`1_CPUE.xlsx`** --- one row per site, year and month (400 rows: 10 sites × 40 quarterly samples). Columns: `Local` (site code), `Year`, `Month` (3, 6, 9, 12), then one column per species (`Genus.species`), CPUE values.

**`2_LOCALS.csv`** (semicolon-separated) --- `Local`, `ID`, `Habitat` (river / tributary), `Dam` (upstream / downstream of the dam), `lat`, `long` (decimal degrees).

**`3_WATER_LEVEL.csv`** (semicolon-separated) --- `HydroYear` (2013--2022), then one column per site (`Abun`, `Kar`, `Mad1`, `Mad2`, `Mad4`, `Mad5`, `Mpar`, `Raul`, `Slou`, `Ssim`) with the annual water-level metric.

**`4_TRAITS.xlsx`** --- `Species`, `Body length`, `Trophic level`, `Diet`, `Habitat`, `Migratory`, `Parental care`, `Fertilization`, `References` (38 species).

**`fish_tree.tre`** --- Newick phylogeny with tips named `Genus_species`; all species in the analyses match the tips.

**`mad_shape.rds`** --- `sf` polygon of the Madeira River used as the base map (read with `readRDS()`).

## Citation

If you use these data or code, please cite the article and this archive:

-   Cararo, E. R. *et al.* (2026). Integrating asymmetry, geography and timescale analyses to uncover fish population synchrony patterns in a tropical reservoir. *JOURNAL*. DOI: *[add]*
-   Cararo, E. R. (2026). Data and code for "Integrating asymmetry, geography and timescale analyses to uncover fish population synchrony patterns in a tropical reservoir" (v1.0.1). Zenodo. DOI: [10.5281/zenodo.23354278](https://doi.org/10.5281/zenodo.23354278)

See also `CITATION.cff`.

## License

-   **Code** (`Scripts/`, `run_all.R`): [MIT License](LICENSE).
-   **Data** (`Data/`): [Creative Commons Attribution 4.0 International (CC BY 4.0)](Data/LICENSE). You may reuse the data, including commercially, as long as you give appropriate credit by citing the article and this archive (see Citation). Trait values in `4_TRAITS.xlsx` come from published sources listed in its `References` column; please cite those sources too.

## Contact

Emanuel R. Cararo --- [emanuelrampanellicararo\@gmail.com](mailto:emanuelrampanellicararo@gmail.com){.email}
