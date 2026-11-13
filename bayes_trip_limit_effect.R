# =============================================================================
# Turtle trip-limit analysis: right-censored negative binomial workflow
# -----------------------------------------------------------------------------
# Order of the script
#   0. Setup, data and checks
#   1. Model fitting functions (common shape, period-specific shape, Post-only)
#   2. Fit models
#   3. Convergence diagnostics
#   4. Prior predictive checks
#   5. Posterior predictive checks and LOO
#   6. Turtles not encountered (excess over the limit) for the capped trips
#   7. Shape-prior sensitivity
#   8. Maximum likelihood comparison
#   9. Projection to a future year (binomial + truncated NB)
#  10. Four-panel figure
#
# Data: long format, one row per trip x species, with columns
#   year, trip, n_sets, species, count, block (Pre/Post), limit,
#   censored ("none"/"right"), count_fit (= count - 1 for censored rows)
#
# NOTE: written without being run against your objects. Run section by section
# and check the printed output against the earlier results.
# =============================================================================

library(dplyr)
library(tidyr)
library(tibble)
library(brms)
library(ggplot2)
library(patchwork)

# ---- 0. Setup ---------------------------------------------------------------
data_file <- "/Users/robert.ahrens/Documents/Projects/Turtles/trip_limit_effect/bayes/trip_turtles_cens.Rdata"          # <- edit path
out_dir   <- "/Users/robert.ahrens/Documents/Projects/Turtles/trip_limit_effect/bayes/output"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

load(data_file)                                  # loads trip_turtles_cens
dat <- trip_turtles_cens %>%
  mutate(block = factor(block, levels = c("Pre", "Post")))   # Pre = reference level

limits    <- c(leatherback = 2, loggerhead = 5)
species_v <- names(limits)

set.seed(123)
ITER <- 2000; CHAINS <- 4; CORES <- 4

# Data checks: counts per block, capped trips, censoring bounds
dat %>% filter(species == "leatherback") %>%
  group_by(block) %>%
  summarise(first_year = min(year), last_year = max(year), n_trips = n())
table(dat$species, dat$censored)
dat %>% filter(censored == "right") %>% count(species, count, count_fit)   # bounds: L - 1
stopifnot(all(dat$count_fit[dat$censored == "none"] == dat$count[dat$censored == "none"]))

obs_summary <- dat %>%
  group_by(species, block) %>%
  summarise(n = n(), mean = mean(count), p_zero = mean(count == 0),
            p_ge_limit = mean(count >= first(limit)), max = max(count), .groups = "drop")
obs_summary

n_capped <- sapply(species_v, function(sp) sum(dat$species == sp & dat$censored == "right"))
n_post_trips <- sapply(species_v, function(sp) sum(dat$species == sp & dat$block == "Post"))
post_years <- sort(unique(dat$year[dat$block == "Post"]))

# ---- 1. Model fitting functions ----------------------------------------------
# type = "common"       : count_fit | cens(censored) ~ block            (one shape)
# type = "period_shape" : ... ~ block, shape ~ block                    (shape by period)
# type = "post_only"    : Post data only, count_fit | cens(censored) ~ 1
fit_nb <- function(sp_data, type = c("common", "period_shape", "post_only"),
                   shape_prior = "gamma(2, 0.5)", sample_prior = "no",
                   intercept_prior = "normal(0, 1.5)", b_prior = "normal(0, 1)") {
  type <- match.arg(type)

  if (type == "common") {
    form <- bf(count_fit | cens(censored) ~ block)
    pri  <- c(prior_string(intercept_prior, class = "Intercept"),
              prior_string(b_prior, class = "b"),
              prior_string(shape_prior, class = "shape"))
  } else if (type == "period_shape") {
    form <- bf(count_fit | cens(censored) ~ block, shape ~ block)
    pri  <- c(prior_string(intercept_prior, class = "Intercept"),
              prior_string(b_prior, class = "b"),
              prior_string("normal(0, 1)", class = "Intercept", dpar = "shape"),
              prior_string("normal(0, 1)", class = "b", dpar = "shape"))
  } else {
    sp_data <- filter(sp_data, block == "Post")
    form <- bf(count_fit | cens(censored) ~ 1)
    pri  <- c(prior_string(intercept_prior, class = "Intercept"),
              prior_string(shape_prior, class = "shape"))
  }

  brm(form, data = sp_data, family = negbinomial(), prior = pri,
      sample_prior = sample_prior,
      control = list(adapt_delta = 0.95),
      chains = CHAINS, cores = CORES, iter = ITER, refresh = 0)
}

# Trip-limit-period (Post) parameters, one row per posterior draw
trip_limit_draws <- function(fit, type) {
  d <- posterior::as_draws_df(fit)
  b0 <- d$b_Intercept
  if (type == "post_only") return(tibble(mu = exp(b0), shape = d$shape))
  b1 <- d[[grep("^b_block", names(d), value = TRUE)]]
  if (type == "common") {
    tibble(mu = exp(b0 + b1), shape = d$shape)
  } else {
    g0 <- d$b_shape_Intercept
    g1 <- d[[grep("^b_shape_block", names(d), value = TRUE)]]
    tibble(mu = exp(b0 + b1), shape = exp(g0 + g1))
  }
}

# Both periods (needed for the figure)
period_draws <- function(fit, type) {
  stopifnot(type %in% c("common", "period_shape"))
  d <- posterior::as_draws_df(fit)
  b0 <- d$b_Intercept
  b1 <- d[[grep("^b_block", names(d), value = TRUE)]]
  if (type == "common") {
    tibble(mu_pre = exp(b0), mu_post = exp(b0 + b1), shape_pre = d$shape, shape_post = d$shape)
  } else {
    g0 <- d$b_shape_Intercept
    g1 <- d[[grep("^b_shape_block", names(d), value = TRUE)]]
    tibble(mu_pre = exp(b0), mu_post = exp(b0 + b1),
           shape_pre = exp(g0), shape_post = exp(g0 + g1))
  }
}

# ---- 2. Fit models -------------------------------------------------------------
fit_types <- c("common", "post_only", "period_shape")   # drop any you do not want
fit_types <- "common"   # drop any you do not want

fits <- lapply(species_v, function(sp) {
  sp_data <- filter(dat, species == sp)
  setNames(lapply(fit_types, function(tp) fit_nb(sp_data, tp)), fit_types)
})
names(fits) <- species_v
saveRDS(fits, file.path(out_dir, "fits.rds"))

# ---- 3. Convergence diagnostics --------------------------------------------------
check_fit <- function(fit, label) {
  np  <- nuts_params(fit)
  sm  <- posterior::summarise_draws(posterior::as_draws_array(fit), "rhat", "ess_bulk", "ess_tail")
  tibble(model = label,
         n_divergent = sum(np$Parameter == "divergent__" & np$Value == 1),
         max_rhat = max(sm$rhat, na.rm = TRUE),
         min_ess_bulk = min(sm$ess_bulk, na.rm = TRUE),
         min_ess_tail = min(sm$ess_tail, na.rm = TRUE))
}
diag_tab <- bind_rows(lapply(species_v, function(sp)
  bind_rows(lapply(fit_types, function(tp) check_fit(fits[[sp]][[tp]], paste(sp, tp))))))
diag_tab        # want R-hat <= 1.01, no divergences, tail ESS well above 400

lapply(species_v, function(sp) summary(fits[[sp]]$common))

# ---- 4. Prior predictive checks ---------------------------------------------------
prior_check <- function(sp_data, type, limit, ...) {
  fp <- fit_nb(sp_data, type, sample_prior = "only", ...)
  y  <- posterior_predict(fp)
  mx <- apply(y, 1, max)
  tibble(p_zero = mean(y == 0), p_ge_limit = mean(y >= limit), p_ge_10 = mean(y >= 10),
         max_median = quantile(mx, .5), max_95 = quantile(mx, .95), max_99 = quantile(mx, .99))
}
prior_tab <- bind_rows(lapply(species_v, function(sp)
  bind_cols(species = sp, prior_check(filter(dat, species == sp), "common", limits[[sp]]))))
prior_tab       # compare with observed p_zero, p_ge_limit and max in obs_summary

# ---- 5. Posterior predictive checks and LOO ------------------------------------------
# Hard-cap (Pre) tail: uncensored counts, so it can be compared directly with data
ppc_tail <- function(fit, sp_data, limit) {
  idx  <- which(sp_data$block == "Pre")
  obs  <- sp_data$count[idx]
  yrep <- posterior_predict(fit)[, idx]
  mx   <- apply(yrep, 1, max)
  pge  <- rowMeans(yrep >= limit)
  tibble(n_trips = length(idx), obs_max = max(obs), pred_max_median = median(mx),
         p_max_ge_obs = mean(mx >= max(obs)),
         obs_prop_ge_L = mean(obs >= limit), pred_prop_ge_L = mean(pge),
         pred_lower = quantile(pge, .025), pred_upper = quantile(pge, .975))
}
ppc_tab <- bind_rows(lapply(species_v, function(sp)
  bind_rows(lapply(setdiff(fit_types, "post_only"), function(tp)
    bind_cols(species = sp, model = tp,
              ppc_tail(fits[[sp]][[tp]], filter(dat, species == sp), limits[[sp]]))))))
ppc_tab

# Observed vs predicted proportions by block, with the limit bin pooled (>= L)
ppc_props <- function(fit, sp_data, limit) {
  yrep <- posterior_predict(fit)
  bind_rows(lapply(levels(sp_data$block), function(b) {
    idx <- which(sp_data$block == b)
    obs <- pmin(sp_data$count[idx], limit)
    rep <- pmin(yrep[, idx, drop = FALSE], limit)
    tibble(block = b, value = 0:limit,
           observed  = sapply(0:limit, function(k) mean(obs == k)),
           pred_mean = sapply(0:limit, function(k) mean(rowMeans(rep == k))),
           lower = sapply(0:limit, function(k) quantile(rowMeans(rep == k), .025)),
           upper = sapply(0:limit, function(k) quantile(rowMeans(rep == k), .975)))
  }))
}
lapply(species_v, function(sp) ppc_props(fits[[sp]]$common, filter(dat, species == sp), limits[[sp]]))

# Does shape differ by period? (common vs period-specific)
lapply(species_v, function(sp) loo(fits[[sp]]$common, fits[[sp]]$period_shape))

# ---- 6. Turtles not encountered for the capped trips ---------------------------------------
# Expected excess per capped trip, E[X - L | X >= L], one value per posterior draw
excess_per_trip <- function(mu, shape, limit, xmax = 2000) {
  x <- limit:xmax
  vapply(seq_along(mu), function(i) {
    f <- dnbinom(x, mu = mu[i], size = shape[i])
    sum((x - limit) * f) / sum(f)
  }, numeric(1))
}

# Total excess over n_capped trips, drawn from the NB truncated below at L
# (inverse CDF: U ~ Uniform[F(L-1), 1]; X = F^-1(U); excess = X - L)
sim_total_saved <- function(mu, shape, limit, n_capped) {
  F_Lm1 <- pnbinom(limit - 1, mu = mu, size = shape)
  vapply(seq_along(mu), function(m) {
    u <- pmin(runif(n_capped, F_Lm1[m], 1), 1 - 1e-12)
    sum(qnbinom(u, mu = mu[m], size = shape[m]) - limit)
  }, numeric(1))
}

saved_summary <- function(fit, type, species, limit, n_c, label = type) {
  dr    <- trip_limit_draws(fit, type)
  delta <- excess_per_trip(dr$mu, dr$shape, limit)
  sim   <- sim_total_saved(dr$mu, dr$shape, limit, n_c)
  tibble(species = species, model = label, n_capped = n_c,
         mu_post = mean(dr$mu), shape_mean = mean(dr$shape), shape_median = median(dr$shape),
         p_reach_limit = mean(pnbinom(limit - 1, mu = dr$mu, size = dr$shape, lower.tail = FALSE)),
         excess_per_trip = mean(delta), total_exact = n_c * mean(delta),
         total_sim_mean = mean(sim), total_sim_median = median(sim),
         lower_95 = quantile(sim, .025, names = FALSE), upper_95 = quantile(sim, .975, names = FALSE))
}

saved_tab <- bind_rows(lapply(species_v, function(sp)
  bind_rows(lapply(fit_types, function(tp)
    saved_summary(fits[[sp]][[tp]], tp, sp, limits[[sp]], n_capped[[sp]])))))
saved_tab        # exact and simulated totals should agree up to Monte Carlo noise

# ---- 7. Shape-prior sensitivity (shape prior varied ONLY) -----------------------------------
shape_priors <- c("gamma(2, 1)", "gamma(2, 0.5)", "gamma(2, 0.25)")
sens_tab <- bind_rows(lapply(species_v, function(sp) {
  sp_data <- filter(dat, species == sp)
  bind_rows(lapply(shape_priors, function(sp_prior) {
    f <- fit_nb(sp_data, "common", shape_prior = sp_prior)
    bind_cols(shape_prior = sp_prior,
              saved_summary(f, "common", sp, limits[[sp]], n_capped[[sp]])) %>%
      bind_cols(n_div = check_fit(f, sp_prior)$n_divergent)
  }))
}))
sens_tab         # is the estimate (and the shape posterior) driven by the prior?

# ---- 8. Maximum likelihood comparison (Post data only) --------------------------------------
mle_post_only <- function(y_fit, censored, limit) {
  nll <- function(th) {
    mu <- exp(th[1]); size <- exp(th[2])
    -sum(ifelse(censored, pnbinom(y_fit, mu = mu, size = size, lower.tail = FALSE, log.p = TRUE),
                          dnbinom(y_fit, mu = mu, size = size, log = TRUE)))
  }
  o <- optim(c(log(mean(y_fit) + 0.1), 0), nll)
  mu <- exp(o$par[1]); size <- exp(o$par[2])
  d <- excess_per_trip(mu, size, limit)
  c(mu = mu, shape = size, excess_per_trip = d)
}
mle_tab <- bind_rows(lapply(species_v, function(sp) {
  s <- filter(dat, species == sp, block == "Post")
  r <- mle_post_only(s$count_fit, s$censored == "right", limits[[sp]])
  tibble(species = sp, mu = r[["mu"]], shape = r[["shape"]],
         excess_per_trip = r[["excess_per_trip"]], total = n_capped[[sp]] * r[["excess_per_trip"]])
}))
mle_tab

# ---- 9. Projection to a future year (binomial + truncated NB) ----------------------------------
# Stage 1: n_capped ~ Binomial(n_trips, p), p = P(trip reaches the limit)
# Stage 2: excess for each capped trip from the NB truncated below at L
p_capped_from_nb <- function(limit, mu, shape) {
  pnbinom(limit - 1, mu = mu, size = shape, lower.tail = FALSE)
}

# Empirical p: Beta posterior from the observed capped fraction in the Post period
p_capped_beta <- function(n_capped, n_trips, n_draws, prior = c(1, 1)) {
  rbeta(n_draws, prior[1] + n_capped, prior[2] + n_trips - n_capped)
}

rand_excess <- function(n_draws, n_trips, limit, mu, shape, p, seed = NULL) {
  # n_trips: trips in ONE year; p: binomial probability a trip is capped
  # mu, shape, p: scalars or vectors of posterior draws (length n_draws)
  if (!is.null(seed)) set.seed(seed)
  stopifnot(length(n_trips) == 1, !missing(p))
  mu <- rep_len(mu, n_draws); shape <- rep_len(shape, n_draws); p <- rep_len(p, n_draws)
  F_Lm1 <- pnbinom(limit - 1, mu = mu, size = shape)
  nc <- rbinom(n_draws, size = n_trips, prob = p)
  excess <- vapply(seq_len(n_draws), function(m) {
    if (nc[m] == 0) return(0)
    u <- pmin(runif(nc[m], F_Lm1[m], 1), 1 - 1e-12)
    sum(qnbinom(u, mu = mu[m], size = shape[m]) - limit)
  }, numeric(1))
  data.frame(n_capped = nc, excess = excess)
}

# Several years at once, sharing posterior draws so totals keep the right uncertainty
rand_excess_years <- function(n_trips, limit, mu, shape, p) {
  # n_trips: named vector of trips per year
  M <- length(mu)
  res <- lapply(n_trips, function(n) rand_excess(M, n, limit, mu, shape, p))
  excess <- sapply(res, `[[`, "excess")
  list(excess = excess, n_capped = sapply(res, `[[`, "n_capped"), total = rowSums(excess))
}

# Example: observed Post-period trip counts per year, using the common-shape fit
project_species <- function(sp, type = "common", p_source = c("nb", "binomial")) {
  p_source <- match.arg(p_source)
  dr <- trip_limit_draws(fits[[sp]][[type]], type)
  nt <- dat %>% filter(species == sp, block == "Post") %>% count(year) %>%
    { setNames(.$n, .$year) }
  p <- if (p_source == "nb") p_capped_from_nb(limits[[sp]], dr$mu, dr$shape)
       else p_capped_beta(n_capped[[sp]], n_post_trips[[sp]], nrow(dr))
  sim <- rand_excess_years(nt, limits[[sp]], dr$mu, dr$shape, p)
  list(sim = sim,
       by_year = tibble(year = as.integer(names(nt)), n_trips = nt,
                        mean = colMeans(sim$excess),
                        lower_95 = apply(sim$excess, 2, quantile, .025),
                        upper_95 = apply(sim$excess, 2, quantile, .975)),
       per_year_avg = quantile(sim$total / length(nt), c(.025, .5, .975)))
}
proj_lb <- project_species("leatherback"); proj_lh <- project_species("loggerhead")
proj_lb$by_year; proj_lb$per_year_avg
proj_lh$by_year; proj_lh$per_year_avg

# Check: the expected annual excess should equal n_trips * p * delta
# (compare mean(rand_excess(...)$excess) with the analytic value)

# ---- 10. Four-panel figure ---------------------------------------------------------------------
# Left: observed proportions by period (bars) and posterior fit with 95% CrI (points)
# Right: posterior predictive distribution of turtles not encountered (capped trips)
fig_type <- "common"     # or "period_shape"

cols <- c(Pre = "#0072B2", Post = "#D55E00")
theme_fig <- theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 11),
        panel.grid.minor = element_blank(),
        plot.background  = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA))

left_panel <- function(fit, sp, L, title) {
  g <- period_draws(fit, fig_type)
  kmax <- max(dat$count[dat$species == sp & dat$block == "Pre"])
  mus    <- list(Pre = g$mu_pre, Post = g$mu_post)
  shapes <- list(Pre = g$shape_pre, Post = g$shape_post)

  fit_df <- bind_rows(lapply(c("Pre", "Post"), function(b) {
    kk <- if (b == "Pre") 0:kmax else 0:L            # Post stops at the limit bin
    pm <- sapply(kk, function(k) {
      if (b == "Post" && k == L) pnbinom(L - 1, mu = mus[[b]], size = shapes[[b]], lower.tail = FALSE)
      else dnbinom(k, mu = mus[[b]], size = shapes[[b]])
    })
    tibble(block = b, count = kk, mean = colMeans(pm),
           lower = apply(pm, 2, quantile, .025), upper = apply(pm, 2, quantile, .975))
  }))

  obs_df <- dat %>% filter(species == sp) %>%
    mutate(count = if_else(block == "Post", pmin(count, L), count)) %>%
    count(block, count) %>% group_by(block) %>% mutate(prop = n / sum(n)) %>% ungroup() %>%
    select(block, count, prop)
  obs_df <- fit_df %>% select(block, count) %>%
    left_join(obs_df, by = c("block", "count")) %>% mutate(prop = coalesce(prop, 0))

  dodge <- position_dodge(width = 0.8)
  ggplot() +
    geom_col(data = obs_df, aes(count, prop, fill = block), position = dodge, width = 0.7, alpha = 0.55) +
    geom_pointrange(data = fit_df, aes(count, mean, ymin = lower, ymax = upper, colour = block),
                    position = dodge, size = 0.35, linewidth = 0.6) +
    geom_vline(xintercept = L + 0.5, linetype = "dashed") +
    annotate("text", x = L + 0.5, y = Inf, label = paste0("limit = ", L),
             angle = 90, vjust = -0.4, hjust = 1.1, size = 3, fontface = "italic") +
    scale_fill_manual(values = cols, name = "Period") +
    scale_colour_manual(values = cols, name = "Period") +
    scale_x_continuous(breaks = 0:kmax) +
    coord_cartesian(xlim = c(-0.6, kmax + 0.8)) +
    labs(x = "Turtles encountered per trip", y = "Proportion of trips", title = title) +
    theme_fig
}

# P(T = k) for the total excess over n_c capped trips, by repeated convolution
total_excess_pmf <- function(mu, phi, L, n_c, K = 80) {
  f <- dnbinom(L:(L + K), mu = mu, size = phi)
  q <- f / sum(f)                                   # distribution of X - L given X >= L
  tot <- q
  if (n_c > 1) for (j in 2:n_c) tot <- convolve(tot, rev(q), type = "open")
  pmax(tot, 0)
}

right_panel <- function(fit, L, n_c, Tmax, title, K = 80) {
  dr  <- trip_limit_draws(fit, fig_type)
  pmf <- t(vapply(seq_len(nrow(dr)),
                  function(i) total_excess_pmf(dr$mu[i], dr$shape[i], L, n_c, K),
                  numeric(n_c * K + 1)))
  mean_pmf <- colMeans(pmf)
  cdf  <- cumsum(mean_pmf)
  q_lo <- which(cdf >= 0.025)[1] - 1
  q_hi <- which(cdf >= 0.975)[1] - 1
  exp_tot <- sum((0:(n_c * K)) * mean_pmf)
  p_tail  <- sum(mean_pmf[(Tmax + 2):length(mean_pmf)])

  df <- tibble(k = 0:Tmax, mean = mean_pmf[1:(Tmax + 1)],
               lower = apply(pmf[, 1:(Tmax + 1)], 2, quantile, .025),
               upper = apply(pmf[, 1:(Tmax + 1)], 2, quantile, .975))
  lab <- sprintf("Expected total = %.2f\n95%% predictive interval: %d\u2013%d\nP(T > %d) = %.3f",
                 exp_tot, q_lo, q_hi, Tmax, p_tail)

  ggplot(df, aes(k, mean)) +
    geom_col(fill = cols[["Post"]], alpha = 0.55, width = 0.7) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25, linewidth = 0.5) +
    annotate("text", x = Inf, y = Inf, label = lab, hjust = 1.05, vjust = 1.3, size = 3.2) +
    scale_x_continuous(breaks = seq(0, Tmax, by = if (Tmax > 10) 2 else 1)) +
    labs(x = "Total turtles not encountered (excess over limit, all capped trips)",
         y = "Posterior predictive probability", title = title) +
    theme_fig
}

p_a <- left_panel(fits$leatherback[[fig_type]], "leatherback", limits[["leatherback"]],
                  "Leatherback: trips by number encountered")
p_b <- right_panel(fits$leatherback[[fig_type]], limits[["leatherback"]], n_capped[["leatherback"]],
                   Tmax = 6, title = sprintf("Leatherback: turtles not encountered (%d capped trips)", n_capped[["leatherback"]]))
p_c <- left_panel(fits$loggerhead[[fig_type]], "loggerhead", limits[["loggerhead"]],
                  "Loggerhead: trips by number encountered")
p_d <- right_panel(fits$loggerhead[[fig_type]], limits[["loggerhead"]], n_capped[["loggerhead"]],
                   Tmax = 20, title = sprintf("Loggerhead: turtles not encountered (%d capped trips)", n_capped[["loggerhead"]]))

fig <- (p_a | p_b) / (p_c | p_d) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")",
                  theme = theme(plot.background = element_rect(fill = "transparent", color = NA))) &
  theme(legend.position = "bottom")
fig

ggsave(file.path(out_dir, "fig_turtle_encounters_and_avoided.png"), fig,
       width = 11, height = 8, dpi = 300, bg = "transparent")
ggsave(file.path(out_dir, "fig_turtle_encounters_and_avoided.pdf"), fig,
       width = 11, height = 8, bg = "transparent")

# ---- Save result tables --------------------------------------------------------------------------
write.csv(diag_tab,  file.path(out_dir, "tab_diagnostics.csv"),   row.names = FALSE)
write.csv(saved_tab, file.path(out_dir, "tab_saved_by_model.csv"), row.names = FALSE)
write.csv(sens_tab,  file.path(out_dir, "tab_shape_prior_sensitivity.csv"), row.names = FALSE)
write.csv(ppc_tab,   file.path(out_dir, "tab_ppc_hardcap_tail.csv"), row.names = FALSE)
write.csv(mle_tab,   file.path(out_dir, "tab_mle_post_only.csv"), row.names = FALSE)
