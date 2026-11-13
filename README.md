# Turtle trip-limit effects in the Hawaii shallow-set longline fishery

A Bayesian, right-censored negative binomial analysis of sea turtle interactions in the Hawaii shallow-set longline fishery under trip-level interaction limits. The main question is how many leatherback and loggerhead interactions were not encountered because trips that reached their limit ended early.

## Background

The fishery has operated under two management regimes for sea turtle interactions:

- **Pre period (hard cap):** fleet-wide interaction caps. Trip-level counts are observed in full.
- **Post period (trip limits):** each trip ends when it reaches **2 leatherback** or **5 loggerhead** interactions.

Under trip limits, the count on a capped trip is a lower bound. Interactions that would have happened if the trip had continued are never observed. This project treats those trips as right-censored observations, estimates the underlying count distribution, and predicts the excess over the limit for each capped trip.

## Approach

Trip-level counts $X$ for each species are modelled as negative binomial with mean $\mu$ and shape $k$, fitted in [`brms`](https://paul-buerkner.github.io/brms/):

```r
count_fit | cens(censored) ~ block            # common shape (main model)
count_fit | cens(censored) ~ block, shape ~ block   # period-specific shape
count_fit | cens(censored) ~ 1                # Post period only
```

Capped trips enter the likelihood as $P(X \ge L)$, encoded as `count_fit = L - 1` with right censoring. The excess for a capped trip is $X - L$ given $X \ge L$. Totals over all capped trips are computed by convolution for each posterior draw.

## Workflow

`bayes_trip_limit_effect.R` runs the analysis in these sections:

| Section | Content |
|:---|:---|
| 0 | Setup, data loading and data checks |
| 1 | Model fitting and posterior-extraction functions |
| 2 | Model fits (common shape by default) |
| 3 | Convergence diagnostics (R-hat, divergences, bulk and tail ESS) |
| 4 | Prior predictive checks |
| 5 | Posterior predictive checks (hard-cap tail, proportions by period) and LOO |
| 6 | Expected and simulated excess over the limit for capped trips |
| 7 | Shape-prior sensitivity: gamma(2, 1), gamma(2, 0.5), gamma(2, 0.25) |
| 8 | Maximum likelihood comparison using Post-period data only |
| 9 | Projection to future years: binomial number of capped trips plus truncated NB excess |
| 10 | Four-panel figure: observed vs fitted distributions, and posterior predictive total excess |

## Data

The script expects an `.Rdata` file containing `trip_turtles_cens`, in long format with one row per trip and species:

| Column | Description |
|:---|:---|
| `year`, `trip` | Trip identifiers |
| `n_sets` | Number of sets on the trip |
| `species` | `leatherback` or `loggerhead` |
| `count` | Observed interactions |
| `block` | `Pre` (hard cap) or `Post` (trip limits) |
| `limit` | Trip limit for the species |
| `censored` | `none` or `right` |
| `count_fit` | Response used in the model (`count - 1` for censored rows) |

The data are not included in this repository.

## Outputs

| File | Content |
|:---|:---|
| `fits.rds` | Fitted `brms` models by species and model type |
| `tab_diagnostics.csv` | Convergence diagnostics |
| `tab_saved_by_model.csv` | Excess per capped trip and total excess, with 95% intervals |
| `tab_shape_prior_sensitivity.csv` | Excess estimates under alternative shape priors |
| `tab_ppc_hardcap_tail.csv` | Posterior predictive check of the Pre-period upper tail |
| `tab_mle_post_only.csv` | Post-only maximum likelihood estimates |
| `fig_turtle_encounters_and_avoided.png/.pdf` | Four-panel summary figure |

## Current results (common-shape model)

| Species | Capped trips | Excess per capped trip | Total excess, mean (95% interval) |
|:---|---:|---:|---:|
| Leatherback | 5 | 0.12 | 0.59 (0–3) |
| Loggerhead | 3 | 1.35 | 4.05 (0–12) |

All models converged (R-hat ≤ 1.006, no divergences, tail ESS > 2,200). Loggerhead estimates are insensitive to the shape prior. Leatherback estimates shift modestly with the prior (total 0.49–0.69).

For loggerheads, the posterior predictive distribution of the total excess has its mode above zero. [`excess_distribution_note.md`](excess_distribution_note.md) explains why this follows from the model.

## Known issues and next steps

- The common-shape model under-predicts the Pre-period loggerhead tail (observed maximum 11 against a predicted median maximum of 5). Because the excess depends directly on the tail, this matters for the loggerhead estimate.
- The period-specific shape model and its LOO comparison are in the script but switched off (`fit_types <- "common"`). Running them is the main remaining check.
- The Post-only MLE implies a lighter loggerhead tail (per-trip excess 0.78) than the common-shape model.

## Requirements

R with `brms` (and a working Stan backend), `posterior`, `dplyr`, `tidyr`, `tibble`, `ggplot2` and `patchwork`.

## Usage

1. Set `data_file` and `out_dir` at the top of `bayes_trip_limit_effect.R`.
2. Choose model types in `fit_types` (Section 2).
3. Run the script section by section and check the printed output at each stage.

## Contact

Robert Ahrens, NOAA Pacific Islands Fisheries Science Center (PIFSC)
