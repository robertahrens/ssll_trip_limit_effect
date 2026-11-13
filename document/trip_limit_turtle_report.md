# An evaluation of the impact of trip limits on sea turtle interaction in the Hawaiʻi shallow-set longline

**Robert Ahrens**<sup>1</sup>

<sup>1</sup> Pacific Islands Fisheries Science Center  
National Marine Fisheries Service  
1845 Wasp Boulevard  
Honolulu, HI 96818

NOAA Data Report [DR-##-#]

[Month] [Year]

<p align="center"><img src="media/doc_seal.png" alt="U.S. Department of Commerce seal" width="110"></p>

U.S. Department of Commerce  
National Oceanic and Atmospheric Administration  
National Marine Fisheries Service  
Pacific Islands Fisheries Science Center

---

# Acknowledgements

Mention throughout this document to trade names or commercial companies is for identification purposes only and does not imply endorsement by the National Marine Fisheries Service, NOAA.

Edited by [First name] [Last name]

**Recommended citation**

Doe J.S., Smith C.R. (Year). *Title in sentence case* (Technical Memoranda Series, publication number). Pacific Islands Fisheries Science Center. doi: hyperlink

*Delete this text. Please go to File \> Properties \> Title \> insert the title of this document \> save*

**Copies of this report are available from**

Pacific Islands Fisheries Science Center  
National Marine Fisheries Service  
National Oceanic and Atmospheric Administration  
1845 Wasp Boulevard, Building #176  
Honolulu, Hawaiʻi 96818

**Or online at**

[<u>https://repository.library.noaa.gov/</u>](https://repository.library.noaa.gov/)

# Table of Contents

- [Acknowledgements](#acknowledgements)
- [List of Tables](#list-of-tables)
- [List of Figures](#list-of-figures)
- [Executive Summary](#executive-summary)
- [Introduction](#introduction)
- [Methods](#methods)
- [Output and Interpretation](#output-and-interpretation)
  - [Interactions avoided](#interactions-avoided)
  - [Sensitivity analyses](#sensitivity-analyses)
  - [Extrapolation to the removal of trip limits in future years](#extrapolation-to-the-removal-of-trip-limits-in-future-years)
  - [Assumptions and limitations](#assumptions-and-limitations)
- [Literature Cited](#literature-cited)
- [Appendix](#appendix)
  - [Prior and Posterior Predictive Checks](#prior-and-posterior-predictive-checks)
  - [R function to generate random excess values](#r-function-to-generate-random-excess-values)

# List of Tables

- [Table 1. Posterior summaries for the censored negative binomial models. CrI, credible interval; ESS, effective sample size.](#table-1)
- [Table 2. Estimated interactions avoided by the trip limit, conditional on trips having reached the limit.](#table-2)
- [Table 3. Sensitivity of the estimated excess interactions to the model and prior. Shape values are posterior means.](#table-3)
- [Table 4. Prior predictive summaries compared with the observed data. L is the trip limit. Each row reports the proportion of simulated trips with no interactions, at or above the limit, and with ten or more, and quantiles of the simulated maximum count.](#table-4)
- [Table 5. Posterior predictive check of the hard-cap (Pre) tail for the common-shape models, 1,076 trips.](#table-5)
- [Table 6. Sensitivity of the loggerhead conclusions to the four high hard-cap trips (maximum likelihood fits to all trips, with the indicated Pre trips removed).](#table-6)
- [Table 7. Dispersion and total excess estimated from the trip-limit period alone, compared with the common-shape models.](#table-7)

# List of Figures

- [Figure 1. Left: observed proportion of trips by number of turtles encountered per trip (bars) and posterior model fit with 95% credible intervals (points), by regulatory period. In the trip-limit period the bin at the limit is censored (≥ limit). Right: posterior predictive distribution of the total number of turtles not encountered across trips that reached the limit, with 95% credible intervals across posterior draws.](#figure-1)

# Executive Summary

This report estimates how many sea turtle interactions were avoided by the trip limit introduced in 2021, which caps leatherback interactions at 2 and loggerhead interactions at 5 per trip. Trips that reached the limit have censored counts, because the true number of encounters is unknown. A right-censored Bayesian negative binomial model was fit to trip-level counts for 2005–2020 (before the limits) and 2021–2025 (trip-limit period). For the capped trips, the model estimates how many encounters would have occurred beyond the limit.

**Key findings**

- **Few trips reached the limit.** Five leatherback trips and three loggerhead trips (of 294 trip-limit-period trips) reached the limit, and none recorded more than the limit.

- **Leatherback:** about 0.6 additional interactions (0.12 per capped trip; predictive 95% interval 0–3) would have occurred across the capped trips without the limit, or roughly 0.1 per year. This ranged from 0.49 to 0.69 across shape priors and from 0.50 to 0.76 when only trip-limit-period data were used.

- **Loggerhead:** about 4.1 additional interactions across the 3 capped trips (1.35 per capped trip; predictive 95% interval 0–12) using shared time period shape model. A model estimated from trip-limit-period data alone gave 2.2 (0–7), and a maximum likelihood fit gave 2.35. The shared shape parameter model is considered the most reasonable alternative given the few trips that had higher numbers of interactions in the pre-trip-limit period.

- **Effect on encounter rates.** Mean encounter rates per trip were higher in the trip-limit period for both species: about 2.9 times higher for loggerheads (95% CrI 2.1–3.9) and 1.5 times higher for leatherbacks (1.0–2.2).

**Confidence and limitations**

The estimates are small in absolute terms and uncertain, because they rest on only five and three capped trips and on the assumed negative binomial form for the unobserved tail. No leatherback trip exceeded two in the trip-limit period, and only one earlier trip recorded more (three), so the leatherback dispersion is not identified by the data and the estimate depends on the prior. The loggerhead estimate is sensitive to whether dispersion is shared across periods, and the pooled shape is held down by four high trips before the limits (counts of 6 to 11). The calculation also assumes that the trip-limit-period encounter rate is what would have applied without the limit, so any avoidance behavior by operators would lead to underestimating the turtles saved.

# Introduction

The Hawaiʻi-based shallow-set longline (SSLL) fishery primarily targets swordfish (*Xiphias gladius*) and is managed under the Fishery Ecosystem Plan for Pelagic Fisheries of the Western Pacific Region, developed by the Western Pacific Fishery Management Council and implemented by the National Marine Fisheries Service (NMFS). Incidental interactions with protected sea turtles, particularly North Pacific loggerhead (*Caretta caretta*) and western Pacific leatherback (*Dermochelys coriacea*) sea turtles, have been central to management of the fishery. Following closure of the shallow-set sector in 2001, the fishery reopened in 2004 under a suite of measures intended to reduce incidental take. These included mandatory large circle hooks, fish rather than squid bait, annual fleet-wide limits on loggerhead and leatherback interactions, and 100% observer coverage (69 FR 17329). Because every trip is observed, recorded interactions and their dates constitute a near-complete census of fleet-wide interactions.

The annual interaction limits, commonly termed *hard caps*, were initially set at 17 loggerhead and 16 leatherback interactions per calendar year. Attainment of either limit required closure of the SSLL fishery for the remainder of the year, prohibiting further shallow-set fishing north of the Equator by vessels holding Hawaiʻi longline limited access permits. The limits were tied to incidental take statements issued under Section 7 of the Endangered Species Act (ESA). Their function, however, differs from that of anticipated take levels. Anticipated take describes the expected number of interactions associated with fishing. Hard caps instead act as operational constraints that can terminate fishing before the end of the year.

The limits have been revised several times in response to new biological opinions, rulemaking, and litigation. In 2010, the loggerhead limit was increased from 17 to 46 interactions, and the annual limit on shallow sets was removed. A subsequent court order reinstated the 17-interaction loggerhead limit in 2011. In November 2012, NMFS increased the limits to 34 loggerhead and 26 leatherback interactions, consistent with the incidental take statement of the 2012 biological opinion (77 FR 60637). Following further litigation, the fishery was closed in May 2018 pursuant to a federal court order. The loggerhead limit was then reduced from 34 to 17 interactions effective January 1, 2019 (83 FR 21939; 83 FR 49484).

Between 2005 and 2019, the fishery closed before the end of the calendar year on four occasions:

- March 20, 2006, upon reaching the loggerhead limit of 17 interactions.

- November 18, 2011, upon reaching the leatherback limit of 16 interactions.

- May 4, 2018, pursuant to a federal court order rather than attainment of an annual limit.

- March 19, 2019, upon reaching the reinstated loggerhead limit of 17 interactions.

In each case, the closure ended shallow-set fishing for the rest of the year. It therefore eliminated further opportunities for interactions with both species, not only the species whose limit was reached.

The management framework changed substantially in 2020. Following a new ESA consultation and associated rulemaking, NMFS made three changes, effective September 17, 2020 (85 FR 57988):

- It removed the annual fleet-wide loggerhead limit.

- It reduced the annual leatherback limit from 26 to 16 interactions.

- It established individual vessel trip limits of five loggerhead and two leatherback interactions.

Under this framework, a vessel reaching a trip limit must cease fishing for the remainder of that trip, but loggerhead interactions no longer trigger a fleet-wide closure. The annual leatherback limit remains, and its attainment closes the fishery for the remainder of the calendar year.

This analysis considers the impact of the trip limit mitigation measure on the number of interactions with leatherback and loggerhead sea turtles in the SSLL. It estimates the number of turtles not encountered because of the trip limits using right censored distribution fits.

# Methods

Only five leatherback and three loggerhead capped trips were observed since the establishment of trip limits. Assuming a common dispersion parameter in both the pre-trip limit and post trip limit periods as well as no between-year variation within periods, a species specific negative binomial model was estimated and used to extrapolate potential encounters beyond the trip limits and the total number of leatherbacks and loggerheads not encountered because trips had reached the limit. The analysis covered 445 trips extracted from the Pacific Islands Region Observer Program System (Pacific Islands Regional Office, 2026). Trips were assigned to a regulatory period by trip year: 2005–2020 (Pre trip-limit period) and 2021–2025 (Post trip-limit period). Under the trip limit, a trip was restricted to a maximum of *L* = 2 leatherback and *L* = 5 loggerhead interactions. In the trip-limit period, a trip with a count equal to *L* was treated as right-censored, because the regulation prevented the true number of encounters from being observed. Leatherback interactions were rare: 91.3% of trips had none in the pre period and 87.1% in the post. Only 1 trip in the pre period recorded 3. Loggerhead interactions were more frequent (87% and 66.3% of trips with none; mean 0.197 per trip in the pre period and 0.561 in the post period), and the pre period included counts up to 11.

The number of interactions ($y_{i}$) that would have occurred on trip *i* were modeled as:

$$y_{i}\sim \mathrm{NegBin}\left( \mu_{i},\phi \right),\qquad \log\left( \mu_{i} \right) = \beta_{0} + \beta_{1} \cdot \mathrm{Period}_{i}$$

$\mathrm{Period}_{i}$ is an indicator equal to 1 for trips in the trip-limit period and 0 otherwise, and $\phi$ is the negative binomial shape (dispersion) parameter. A single $\phi$ was shared across both periods. An uncensored trip contributes its negative binomial probability mass, $f(y_{i}|\mu_{i},\ \phi)$ to the likelihood. A right-censored trip contributes the probability that the true count was at least the limit:

$$\Pr\left( y_{i} \geq L \right) = 1 - F(L - 1|\mu_{i},\ \phi)$$

where *F* is the negative binomial cumulative distribution function. The model was fitted using the *brms* (Bürkner, 2017) package censoring term cens(). For discrete families, brms evaluates the censored contribution with the Stan (Carpenter et al., 2017) complementary CDF, which gives $\Pr(Y > y)$ rather than $\Pr(Y \geq y)$. The response supplied for censored trips was therefore set to *L* − 1 (verified in the Stan data as bounds of 1 for leatherback and 4 for loggerhead). Uncensored trips retained their observed counts.

A preliminary model that added a random intercept for trip year produced divergent transitions (16 for leatherback and 12 for loggerhead) with a poorly identified year-level standard deviation. Because each year falls entirely within one period and only a small number of trips were capped, the year-level variance was not well informed by the data, and the random effect was omitted. The no-random-effects model sampled without divergent transitions and is the model reported here.

Weakly informative priors were used: *β*<sub>0</sub> ~ Normal(0, 1.5), *β*<sub>1</sub> ~ Normal(0, 1), and $\phi$ ~ Gamma(2, 0.5), where the regression priors are on the log scale. Models were fitted in R using the brms package (Bürkner 2017), which interfaces with Stan and its No-U-Turn sampler. Each model used four chains of 2000 iterations (1000 warm-up), giving 4000 post-warm-up draws, with adapt_delta = 0.95. Convergence was assessed with the potential scale reduction factor (R-hat), bulk and tail effective sample sizes (ESS) (Vehtari et al., 2021), and the number of divergent transitions.

Prior predictive checks compared the number of zero trips, the proportion of trips at or above the limit, and the maximum count implied by the priors with the observed data. The initial leatherback priors implied that 45% of trips had no interactions, 36% reached the limit and 8% had ten or more, against observed values of 89%, 1.3% and 0%, so we also fitted tightened priors (Intercept ~ Normal(−2, 1), period effect ~ Normal(0, 0.5), shape ~ Gamma(2, 1)), which implied 84% zeros, 4.1% at or above the limit and 0.04% at ten or more. Sensitivity to the shape prior was assessed for both species by refitting with Gamma(2, 0.25), Gamma(2, 0.5) and Gamma(2, 1) priors. For loggerheads, an intercept prior centered on the observed mean (Normal(−0.5, 1)) was also examined.

Posterior predictive checks focused on the hard-cap period, where counts are uncensored and the tail of the distribution can be compared directly with data. We computed the posterior predictive probability that the maximum count in the hard-cap period was at least as large as observed, and compared the observed and predicted proportions of hard-cap trips at or above the limit.

The number of additional interactions that would have occurred on the capped trips had the limit not been in place is conditional on $Y \geq L$. For posterior draw *m*, with trip-limit-period mean $\mu_{m} = \exp(\beta_{0,m} + \beta_{1,m})$ and shape $\phi_{m}$, the expected excess over the limit per capped trip is

$$\delta_{m} = \frac{\sum_{x \geq L}(x - L)\,f(x \mid \mu_{m}, \phi_{m})}{\sum_{x \geq L} f(x \mid \mu_{m}, \phi_{m})}$$

and the expected total avoided across the *n*<sub>c</sub> capped trips is $A_{m} = n_{c}\delta_{m}$. Sums were evaluated numerically to $x = 500$. To express the additional uncertainty from the random number of encounters on each capped trip, the total was also simulated directly. For each posterior draw, *n*<sub>c</sub> counts were drawn from the negative binomial distribution truncated below at *L* using the inverse-CDF method (*U* ~ Uniform[*F*(*L* − 1), 1]; *X* = *F*<sup>−1</sup>(*U*)), and the total excess expressed as $T_{m} = \sum_{j}(X_{j} - L)$. The median and 2.5th–97.5th percentiles of *T* form the predictive interval.

Because the expected excess depends strongly on the tail of the distribution, the Bayesian result was compared with a maximum likelihood fit of the same censored negative binomial likelihood to trip-limit-period counts only, so that the dispersion parameter was estimated without information from the pre-trip limit period. The conditional expectation was evaluated at the maximum likelihood estimates.

# Output and Interpretation

All parameters had R-hat ≤ 1.01, and the smallest tail ESS was 2253 (Table 1). Neither model produced divergent transitions. For leatherbacks, the mean number of interactions per trip was approximately 0.15 in the trip-limit period compared with 0.10 before the trip limits (period effect 0.40, 95% CrI 0.02–0.77; rate ratio 1.5, 95% CrI 1.0–2.2), so the credible interval for the period effect only marginally excluded zero. The leatherback shape was not well identified by the data (posterior mean 2.3, 95% CrI 0.56–7.32): its posterior mean moved with the prior mean across the three shape priors (1.5, 2.3 and 4.2 for prior means of 2, 4 and 8), whereas the mean per trip in the trip-limit period was stable at 0.15.

<a id="table-1"></a>**Table 1.** Posterior summaries for the censored negative binomial models. CrI, credible interval; ESS, effective sample size.

| **Species / parameter** | **Mean** | **SD** | **95% CrI**  | **Min. tail ESS** |
|-------------------------|----------|--------|--------------|-------------------|
| Leatherback             |          |        |              |                   |
| β₀ (intercept)          | −2.33    | 0.10   | −2.53, −2.14 | 2828              |
| β₁ (trip-limit period)  | 0.40     | 0.19   | 0.01, 0.76   | 2251              |
| φ (shape)               | 2.29     | 1.83   | 0.56, 7.35   | 2780              |
| Loggerhead              |          |        |              |                   |
| β₀ (intercept)          | -1.61    | 0.09   | −1.78, −1.45 | 2437              |
| β₁ (trip-limit period)  | 1.04     | 0.15   | 0.74, 1.34   | 2429              |
| φ (shape)               | 0.35     | 0.05   | 0.26, 0.47   | 2417              |

Intercept and period coefficients are on the log scale. Pre trip limit period is the reference level.

For loggerheads, the mean number of interactions per trip was approximately 0.20 before the trip limits and 0.57 in the trip-limit period (period effect 1.05, 95% CrI 0.75–1.35; rate ratio 2.9, 95% CrI 2.1–3.9). Under the common-shape model, the shape was well identified (0.35, 95% CrI 0.26–0.46) and was essentially unchanged across the three shape priors (posterior means 0.347–0.351). In the posterior predictive check of the pre-trip-limit tail, the probability of a maximum at least as large as the observed 11 was 0.004, and the predicted proportion of trips at or above five was 0.12% (95% interval 0–0.37%) compared with 0.37% observed, indicating that the common-shape model underpredicts the pre-trip-limit tail. This result depends on four trips with counts of 6 to 11 (see the Appendix).

## Interactions avoided

Across the five capped leatherback trips, the expected number of additional interactions had no limit been in place was 0.59 (about 0.12 per capped trip; Table 2 and Figure 1). Across the three capped loggerhead trips, it was 4.12 (about 1.35 per capped trip). The predictive intervals are wide (leatherback 0–2; loggerhead 0–12), reflecting the small number of capped trips and the heavy tail of the fitted loggerhead distribution.

<a id="table-2"></a>**Table 2.** Estimated interactions avoided by the trip limit, conditional on trips having reached the limit.

| **Species** | **Limit (L)** | **Capped trips (n<sub>c</sub>)** | **Excess per capped trip, δ** | **Expected total, A** | **Predictive median** | **Predictive 95% interval** |
|-------------|---------------|----------------------------------|-------------------------------|-----------------------|-----------------------|-----------------------------|
| Leatherback | 2             | 5                                | 0.117                         | 0.587                 | 0                     | 0–3                         |
| Loggerhead  | 5             | 3                                | 1.351                         | 4.115                 | 3                     | 0–12                        |

δ and A are posterior means of the exact conditional expectations. The predictive median and interval are from simulation of the total excess across capped trips.

<a id="figure-1"></a>

![Four-panel figure showing frequency distributions for leatherback (upper panels) and loggerhead (lower panels) sea turtles. The upper-left panel shows trip level interaction pre and post trip limits in orange and blue declining quickly from 0 to a max of 3. The upper-right panel shows a right-skewed distribution for expected turtles not encountered that declines from 0 to a max of 6. The lower-left panel shows trip level interaction pre and post trip limits in orange and blue declining quickly from 0 to a max of 11. The lower-right panel shows a right-skewed distribution for expected turtles not encountered that increases from 0 to 2 and then declines to a max of 20.](figures/fig1_encounters_and_avoided.png)

**Figure 1.** Left: observed proportion of trips by number of turtles encountered per trip (bars) and posterior model fit with 95% credible intervals (points), by regulatory period. In the trip-limit period the bin at the limit is censored (≥ limit). Right: posterior predictive distribution of the total number of turtles not encountered across trips that reached the limit, with 95% credible intervals across posterior draws.

## Sensitivity analyses

The shape prior was varied alone (Gamma(2, 1), Gamma(2, 0.5) and Gamma(2, 0.25), with prior means of 2, 4 and 8) and the estimates recomputed (Table 3). The leatherback total ranged from 0.49 to 0.69, and the posterior mean shape moved with the prior (1.5 to 4.2), consistent with dispersion that is only weakly identified. The common-shape loggerhead estimate was insensitive to the prior (4.04–4.07, shape about 0.35), showing that the data fix the shape of the pooled model, although that value is held down by the four extreme Pre trips (see Appendix).

<a id="table-3"></a>**Table 3.** Sensitivity of the estimated excess interactions to the model and prior. Shape values are posterior means.

| **Model and prior**             | **Shape** | **Excess per capped trip** | **Expected total** |
|---------------------------------|-----------|----------------------------|--------------------|
| ***Leatherback***               |           |                            |                    |
| Shape ~ Gamma(2, 1)             | 1.5       | 0.137                      | 0.687              |
| Shape ~ Gamma(2, 0.5) (initial) | 2.28      | 0.117                      | 0.587              |
| Shape ~ Gamma(2, 0.25)          | 4.21      | 0.0983                     | 0.497              |
| ***Loggerhead***                |           |                            |                    |
| Shape ~ Gamma(2, 1)             | 0.347     | 1.35                       | 4.06               |
| Shape ~ Gamma(2, 0.5) (initial) | 0.35      | 1.36                       | 4.07               |
| Shape ~ Gamma(2, 0.25)          | 0.351     | 1.35                       | 4.04               |

## Extrapolation to the removal of trip limits in future years

Estimating the interactions that would no longer be avoided in a future year, after the trip limits are removed, depends on the expected number of trips in that year and on how many of those trips would reach the limit. The simplest approach, ceteris paribus, applies the expectation for the period assessed, which is 0.12 leatherbacks and 0.81 loggerheads per year (0.59 / 5 and 4.05 / 5, where 5 is the number of years, 2021–2025). The intervals for the period average are 0 to 0.6 per year for leatherback and 0 to 2.4 per year for loggerhead (from the predictive intervals of 0 to 3 and 0 to 12 for the totals).

Incorporating the uncertainty in the number of trips that would reach the limit requires a two-stage simulation in which that number is random. In the first stage, the number of capped trips is drawn as $n_{c}\sim Binomial(n,p)$, where *n* is the number of trips in the year and *p* is the probability that a trip reaches the limit. In the empirical approach, *p* is estimated from the observed fraction of trip-limit-period trips that reached the limit (5 of 294 for leatherback and 3 of 294 for loggerhead), using a conjugate Beta(1, 1) prior so that *p* ~ Beta(6, 290) for leatherback and Beta(4, 292) for loggerhead. Alternatively, *p* can be computed for each posterior draw as Pr(*Y* ≥ *L*) from the fitted negative binomial, which keeps *p* consistent with the model. In the second stage, for each of the $n_{c}$ capped trips, a count is drawn from the negative binomial, using the posterior draws of its parameters, truncated below at *L* by the inverse-CDF method, and the excess over the limit is summed to give the total for the year.

## Assumptions and limitations

The estimates are conditional on the negative binomial form, which determines the extrapolated tail beyond the limit, and they rest on only five (leatherback) and three (loggerhead) capped trips. No leatherback trip exceeded two in the trip-limit period, and only one trip before the limits recorded more (three), so the leatherback tail beyond the limit, and its dispersion, are inferred almost entirely from the model and the prior.

A single dispersion parameter was assumed for both periods. For loggerheads, the common-shape model underpredicts the pre-trip-limit tail (probability 0.004 of a maximum as large as the observed 11), but this depends on four trips with counts of 6 to 11. An exploratory maximum likelihood model with period-specific shape gave heavier tails before the limits (shape 0.20) than during them (0.67) and a lower total (2.35). The evidence for period-specific dispersion is strong when all trips are used (likelihood ratio 16.5) but borderline when the four extreme trips are removed (likelihood ratio 3.3, P = 0.07), in which case the pooled shape rises from 0.33 to 0.50 and the common-shape total falls from 4.2 to 3.0. The common-shape estimate of 4.05 should therefore be read as the upper end of a range of about 2.2 to 4.1 across the trip-limit-only, maximum likelihood and common-shape results.

The analysis also assumes that the encounter rate observed under the trip limit is the rate that would have applied without it. If operators changed behavior to avoid reaching the limit, the counterfactual is underestimated. Mean encounter rates differed between periods (about 1.5 times higher for leatherbacks and 2.9 times higher for loggerheads after the limits), so the periods are not exchangeable, and assuming the same dispersion in both is a strong assumption. Finally, because reaching the limit may end a trip, encounters with the other species after the limit was reached may also have been truncated, which is not accounted for here.

# Literature Cited

Bürkner, P.-C. (2017). brms: An R package for Bayesian multilevel models using Stan. *Journal of statistical software*, *80*, 1-28. <https://doi.org/10.18637/jss.v080.i01>

Carpenter, B., Gelman, A., Hoffman, M. D., Lee, D., Goodrich, B., Betancourt, M., Brubaker, M., Guo, J., Li, P., & Riddell, A. (2017). Stan: A probabilistic programming language. *Journal of statistical software*, *76*, 1-32. <https://doi.org/10.18637/jss.v076.i01>

Pacific Islands Regional Office. (2026). *Longline Observer Data System* <https://www.fisheries.noaa.gov/inport/item/9027>

Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Bürkner, P. (2021). Rank-normalization, folding, and localization: An improved R-hat for assessing convergence of MCMC. *Bayesian Analysis*, *16*(2), 667–718. <https://doi.org/10.1214/20-BA1221>

# Appendix

## Prior and Posterior Predictive Checks

The number of turtles not encountered on capped trips, depends on the upper tail of the negative binomial distribution above the trip limit. That part of the distribution is never observed in the trip-limit period, because the regulation censors counts at the limit. The prior predictive checks ask whether the priors imply plausible counts before the data are seen. The posterior predictive checks use the Pre period, in which counts are uncensored and the tail can be compared directly with data. All model checks used the 1,370 trips from 2005–2025 (1,076 Pre and 294 Post).

### Prior predictive checks

Prior predictive draws were generated by sampling the intercept, period effect and shape from their priors and simulating a negative binomial count for each trip, with the period indicator set from the data. Three statistics were compared with the observed data: the proportion of trips with no interactions, the proportion at or above the limit, and the maximum count across all trips (Table 4). Three prior sets were examined: the initial priors (intercept ~ Normal(0, 1.5), period effect ~ Normal(0, 1), shape ~ Gamma(2, 0.5)); tightened leatherback priors (intercept ~ Normal(−2, 1), period effect ~ Normal(0, 0.5), shape ~ Gamma(2, 1)); and centered loggerhead priors (intercept ~ Normal(−0.5, 1), period effect ~ Normal(0, 0.5), shape ~ Gamma(2, 0.5)).

<a id="table-4"></a>**Table 4.** Prior predictive summaries compared with the observed data. L is the trip limit. Each row reports the proportion of simulated trips with no interactions, at or above the limit, and with ten or more, and quantiles of the simulated maximum count.

| **Species and prior set** | **P(zero)** | **P(≥ L)** | **P(≥ 10)** | **Maximum count: median / 95% / 99%** |
|---------------------------|-------------|------------|-------------|---------------------------------------|
| ***Leatherback (L = 2)*** |             |            |             |                                       |
| Observed (1,370 trips)    | 0.904       | 0.010      | 0.000       | 3 (observed maximum)                  |
| Initial priors            | 0.441       | 0.374      | 0.079       | 9 / 89 / 270                          |
| Tightened priors          | 0.839       | 0.040      | 0.0003      | 3 / 9 / 18                            |
| ***Loggerhead (L = 5)***  |             |            |             |                                       |
| Observed (1,370 trips)    | 0.826       | 0.005      | 0.001       | 11 (observed maximum)                 |
| Initial priors            | 0.445       | 0.158      | 0.073       | 9 / 84 / 281                          |
| Centered priors           | 0.553       | 0.045      | 0.009       | 6 / 22 / 46                           |

Proportions are averaged over prior draws and trips. Maxima are taken over the 1,370 trips in each simulated data set. Simulation error is about ±0.005 for proportions and larger for the upper quantiles of the maximum.

The initial priors were diffuse on the count scale. For leatherbacks they implied that 37% of trips would reach the limit and that a typical data set would contain a trip with nine leatherbacks, whereas 1.0% of trips reached the limit and the largest observed count was three. The same priors implied a 99th percentile maximum of about 270. The tightened leatherback priors reproduced the observed zero and limit proportions closely and gave a plausible maximum (median 3). For loggerheads, centering the intercept on the observed mean reduced the implied proportion at or above the limit from 15.8% to 4.5%, which is still above the observed 0.5% but not implausible, and the data dominate the intercept and period effect with 1,370 trips.

#### Posterior predictive check of the hard-cap tail

For the common-shape models, we generated posterior predictive counts for each of the 1,076 Pre trips and compared the maximum count and the proportion of trips at or above the limit with the observed values (Table 5). Because posterior predictive draws are uncensored, the check is valid only where counts are fully observed, so it was applied to the Pre period alone.

<a id="table-5"></a>**Table 5.** Posterior predictive check of the hard-cap (Pre) tail for the common-shape models, 1,076 trips.

| **Species** | **Observed maximum** | **Predicted maximum (median)** | **P(predicted max ≥ observed)** | **Observed % at or above L** | **Predicted % at or above L (95% interval)** |
|-------------|----------------------|--------------------------------|---------------------------------|------------------------------|----------------------------------------------|
| Leatherback | 3                    | 2                              | 0.406                           | 0.84                         | 0.69 (0.19–1.39)                             |
| Loggerhead  | 11                   | 5                              | 0.004                           | 0.37                         | 0.12 (0.00–0.37)                             |

The predicted maximum is the maximum count over the 1,076 Pre trips in each posterior predictive data set. The P-value is the proportion of posterior draws in which the predicted maximum was at least as large as the observed maximum.

The leatherback model reproduced both statistics: the observed maximum of three was well within the predictive distribution (P = 0.40), and the observed proportion of trips at or above two lay inside the predicted interval. This check has limited power for leatherbacks, because the observed maximum comes from a single trip. For loggerheads, the probability of a predicted maximum as large as the observed 11 was 0.004, the median predicted maximum was 5, and the observed proportion of trips at or above five (0.37%, four trips) lay at the upper limit of the predicted interval. This result depends on four Pre trips with counts of 6, 8, 9 and 11, which are the only Pre trips above four, and the maximum is a post hoc statistic that a single extreme trip can drive. The posterior check conclusion depends on those trips (Table 6). Removing the trip with 11 alone left strong evidence that dispersion differs between periods (likelihood ratio 11.9, P = 0.0006), but removing all four trips reduced it to a borderline result (likelihood ratio 3.3, P = 0.07) and raised the pooled shape from 0.33 to 0.50. The common shape of about 0.35 is held down by four extreme hard-cap trips, and the data cannot determine whether the Pre tail is genuinely heavier or whether these trips are unusual events. The check was not applied to the trip-limit period, where counts at the limit are censored and the tail above the limit is unobserved.

<a id="table-6"></a>**Table 6.** Sensitivity of the loggerhead conclusions to the four high hard-cap trips (maximum likelihood fits to all trips, with the indicated Pre trips removed).

| **Pre trips removed**          | **Common shape** | **Common-shape total** | **Likelihood ratio: common vs period-specific shape** | **Trip-limit-only total** |
|--------------------------------|------------------|------------------------|-------------------------------------------------------|---------------------------|
| None                           | 0.33             | 4.17                   | 16.5 (P \< 0.0001)                                    | 2.35                      |
| Trip with 11                   | 0.37             | 3.79                   | 11.9 (P = 0.0006)                                     | 2.35                      |
| Four trips with 6, 8, 9 and 11 | 0.50             | 2.97                   | 3.3 (P = 0.07)                                        | 2.35                      |

The common-shape total is the expected number of turtles above the limit across the 3 capped trips from the single-shape model. The trip-limit-only total uses trip-limit data only and is therefore unchanged by removing Pre trips. The likelihood ratio compares one shape with a period-specific shape (1 degree of freedom).

#### Identifiability of dispersion from trip-limit data alone

Because the common-shape loggerhead model was in tension with the Pre tail and with the trip-limit data, we estimated the loggerhead model from trip-limit-period data only, with no period effect and no information from the hard-cap period. In this model, shape is informed by the counts below the limit and by the number of trips at or above it, with the censored likelihood accounting for the limit. The likelihood profile for shape (Table 7) shows that the trip-limit data identify the loggerhead dispersion reasonably well (95% interval 0.43–1.15) but do not identify the leatherback dispersion, for which the interval extends from 0.27 to unbounded values. The common-shape value for loggerheads (0.35) lies outside the profile interval of the trip-limit data (log-likelihood 3.9 below the maximum against a cutoff of 1.9), whereas the pooled shape of 0.50 obtained without the four extreme Pre trips lies inside it.

<a id="table-7"></a>**Table 7.** Dispersion and total excess estimated from the trip-limit period alone, compared with the common-shape models.

| **Quantity (trip-limit period, 294 trips)**            | **Leatherback**          | **Loggerhead**               |
|--------------------------------------------------------|--------------------------|------------------------------|
| Counts of 0, 1, 2, …, and capped (≥ L)                 | 256, 33; 5 capped        | 195, 59, 23, 11, 3; 3 capped |
| Maximum likelihood shape (95% profile interval)        | 0.93 (0.27 to unbounded) | 0.66 (0.43–1.15)             |
| Maximum likelihood total excess (profile range)        | 0.78 (0.25–2.05)         | 2.37 (1.53–3.35)             |
| Post-only Bayesian: posterior mean shape               | 3.42                     | 0.79                         |
| Post-only Bayesian: total excess                       | 0.50                     | 2.18                         |
| Common-shape Bayesian: shape, total excess             | 2.28, 0.59               | 0.35, 4.05                   |
| Log-likelihood at the common-shape value minus maximum | —                        | −3.9 (95% cutoff −1.9)       |

Profile intervals are likelihood-based. “Total excess” is the expected number of turtles above the limit across the capped trips (5 leatherback, 3 loggerhead).

These results led to different structures for the two species. For leatherbacks, the trip-limit data cannot identify the dispersion, and the hard-cap data contain usable tail information (nine Pre trips at or above the limit, one with three) and passed the tail check, so the common-shape model was retained. For loggerheads, the pooled shape is held down by four extreme Pre trips and is in tension with the trip-limit data, so the trip-limit-only model was used as the primary model and the common-shape model as an upper sensitivity bound. In both cases, the estimate of the number not encountered rests on the negative binomial form above the limit, which none of the checks can verify, since counts above the limit are never observed under the trip limit.

## R function to generate random excess values

```r
rand_excess <- function(n_draws, n_trips, limit, mu, shape, p = NULL, seed = NULL) {
  # n_trips: vector of trip counts, one per year (name it by year)
  # mu, shape: posterior draws of the NB parameters
  # p probability of a trip hitting the limit
  # limit is the trip cap
  if (!is.null(seed)) set.seed(seed)
  p_below <- pnbinom(limit - 1, mu = mu, size = shape) # F(L - 1)
  nc <- rbinom(n_draws, size = n_trips, prob = p) # capped trips per draw
  excess <- vapply(nc, function(k) {
    if (k == 0) return(0)
    u <- pmin(runif(k, p_below, 1), 1 - 1e-12) # U ~ Uniform[F(L-1), 1]
    sum(qnbinom(u, mu = mu, size = shape) - limit) # X = F^-1(U); excess = X - L
  }, numeric(1))
  return(excess)
}

#leatherbacks
rand_excess(n_draws = 10, n_trips = 60, limit = 2, mu = 0.147, shape = 2.28, p = 6/290, seed = NULL)

#loggerheads
rand_excess(n_draws = 10, n_trips = 60, limit = 5, mu = 0.574, shape = 0.344, p = 4/292, seed = NULL)
```
