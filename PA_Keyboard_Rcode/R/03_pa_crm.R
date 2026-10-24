# =============================================================================
# PA-CRM: adult-information borrowing comparator
# =============================================================================
# This module implements the PA-CRM comparison design. Adult doses are mapped to
# the pediatric effective-dose scale; adult and pediatric DLT data are combined
# using a normalized adult power prior and Bayesian model averaging (BMA).
#
# Function guide:
#   pa_crm_standardized_doses()       Construct pediatric/adult effective doses.
#   make_pa_crm_quadrature(),
#     make_pa_crm_basis()             Precompute numerical integration quantities.
#   bernoulli_count_log_likelihood() Calculate the binomial log-likelihood.
#   pa_crm_fit()                      Fit candidate models and calculate BMA weights.
#   simulate_pa_crm_pediatric()       Simulate PA-CRM pediatric dose assignment.
# =============================================================================

# =============================================================================
# PA-CRM model with deterministic quadrature and model averaging
# =============================================================================
pa_crm_standardized_doses <- function(config) {
  child_x          <- logit(config$pa_crm_skeleton) - config$pa_crm_intercept
  adult_dose_mg_m2 <- config$adult_doses_mg / config$bsa_conversion
  adult_x          <- linear_interpolate_extrapolate(
    config$pediatric_doses_mg_m2,
    child_x,
    adult_dose_mg_m2
  )
  list(
    child_x          = child_x,
    adult_x          = adult_x,
    adult_dose_mg_m2 = adult_dose_mg_m2
  )
}

make_pa_crm_quadrature <- function(config) {
  cache_key <- paste(
    "pa_crm_quad",
    paste(config$pa_crm_skeleton, collapse  = ","),
    config$pa_crm_intercept,
    paste(config$adult_doses_mg, collapse  = ","),
    paste(config$pediatric_doses_mg_m2, collapse  = ","),
    config$bsa_conversion,
    config$target,
    config$pa_crm_quad_n_per_side,
    paste(config$pa_crm_alpha_bounds, collapse  = ","),
    sep  = "|"
  )
  if (exists(cache_key, envir  = .PAK_NUMERIC_CACHE, inherits  = FALSE)) {
    return(get(cache_key, envir  = .PAK_NUMERIC_CACHE, inherits  = FALSE))
  }
  standardized    <- pa_crm_standardized_doses(config)
  x1              <- standardized$child_x[1L]
  threshold_ratio <- (logit(config$target) - config$pa_crm_intercept) / x1
  if (!is.finite(threshold_ratio) || threshold_ratio <= 0) {
    stop("PA-CRM lowest-dose safety threshold is not well-defined.", call.  = FALSE)
  }
  threshold <- log(threshold_ratio)
  bounds    <- config$pa_crm_alpha_bounds
  if (threshold <= bounds[1L] || threshold >= bounds[2L]) {
    q <- map_legendre_interval(bounds[1L], bounds[2L],
      2L * config$pa_crm_quad_n_per_side)
  } else {
    q_lower <- map_legendre_interval(bounds[1L], threshold,
      config$pa_crm_quad_n_per_side)
    q_upper <- map_legendre_interval(threshold, bounds[2L],
      config$pa_crm_quad_n_per_side)
    q <- list(nodes  = c(q_lower$nodes, q_upper$nodes),
      weights  = c(q_lower$weights, q_upper$weights))
  }
  answer <- list(
    alpha                         = q$nodes,
    integration_weight            = q$weights,
    lowest_safety_alpha_threshold = threshold,
    standardized                  = standardized
  )
  assign(cache_key, answer, envir  = .PAK_NUMERIC_CACHE)
  answer
}

make_pa_crm_basis <- function(config) {
  cache_key <- paste(
    "pa_crm_basis", config_signature(config),
    format(config$pa_crm_alpha_mean, digits  = 16),
    format(config$pa_crm_alpha_variance, digits  = 16), sep  = "|"
  )
  if (exists(cache_key, envir  = .PAK_NUMERIC_CACHE, inherits  = FALSE)) {
    return(get(cache_key, envir  = .PAK_NUMERIC_CACHE, inherits  = FALSE))
  }
  quadrature            <- make_pa_crm_quadrature(config)
  alpha                 <- quadrature$alpha
  exp_alpha             <- exp(alpha)
  pediatric_probability <- inv_logit(
    config$pa_crm_intercept +
    outer(exp_alpha, quadrature$standardized$child_x, `*`)
  )
  adult_probability <- inv_logit(
    config$pa_crm_intercept +
    outer(exp_alpha, quadrature$standardized$adult_x, `*`)
  )
  log_base_weight <- log(quadrature$integration_weight) + stats::dnorm(
    alpha,
    mean = config$pa_crm_alpha_mean,
    sd   = sqrt(config$pa_crm_alpha_variance),
    log  = TRUE
  )
  answer <- c(quadrature, list(
      pediatric_probability = pediatric_probability,
      adult_probability     = adult_probability,
      log_base_weight       = log_base_weight
  ))
  assign(cache_key, answer, envir  = .PAK_NUMERIC_CACHE)
  answer
}

bernoulli_count_log_likelihood <- function(probability_matrix, x, n) {
  p <- clamp(probability_matrix, 1e-15, 1 - 1e-15)
  x <- as.numeric(x)
  n <- as.numeric(n)
  as.numeric(log(p) %*% x + log1p(-p) %*% (n - x))
}

pa_crm_fit <- function(adult_n, adult_x, pediatric_n, pediatric_x, config) {
  if (length(adult_n) != 5L || length(adult_x) != 5L ||
    length(pediatric_n) != 5L || length(pediatric_x) != 5L) {
    stop("PA-CRM count vectors must all have length five.", call.  = FALSE)
  }
  basis                 <- make_pa_crm_basis(config)
  alpha                 <- basis$alpha
  standardized          <- basis$standardized
  pediatric_probability <- basis$pediatric_probability
  adult_probability     <- basis$adult_probability
  pediatric_ll          <- bernoulli_count_log_likelihood(
    pediatric_probability, pediatric_x, pediatric_n
  )
  adult_ll        <- bernoulli_count_log_likelihood(adult_probability, adult_x, adult_n)
  log_base_weight <- basis$log_base_weight

  n_model                           <- length(config$pa_crm_delta)
  log_joint_evidence                <- numeric(n_model)
  log_adult_power_normalizer        <- numeric(n_model)
  log_model_marginal                <- numeric(n_model)
  node_weights                      <- vector("list", n_model)
  model_pediatric_mean              <- matrix(NA_real_, nrow  = n_model, ncol  = 5L)
  model_lowest_overdose_probability <- numeric(n_model)
  model_alpha_mean                  <- numeric(n_model)
  for (m in seq_len(n_model)) {
    log_adult_kernel              <- log_base_weight + config$pa_crm_delta[m] * adult_ll
    log_adult_power_normalizer[m] <- log_sum_exp(log_adult_kernel)
    log_kernel                    <- log_base_weight + pediatric_ll +
    config$pa_crm_delta[m] * adult_ll
    log_joint_evidence[m] <- log_sum_exp(log_kernel)
    # Equation (8) of Li and Yuan: marginal pediatric likelihood under the
    # normalized adult power prior, not the unnormalized joint integral.
    log_model_marginal[m] <- log_joint_evidence[m] -
    log_adult_power_normalizer[m]
    weight                    <- exp(log_kernel - log_joint_evidence[m])
    node_weights[[m]]         <- weight
    model_pediatric_mean[m, ] <- colSums(pediatric_probability * weight)
    # The quadrature is split exactly at the P1 safety threshold.
    model_lowest_overdose_probability[m] <- sum(
      (pediatric_probability[, 1L] > config$target) * weight
    )
    model_alpha_mean[m] <- sum(alpha * weight)
  }
  log_model_weight            <- log(config$pa_crm_delta_prior) + log_model_marginal
  model_weight                <- exp(log_model_weight - log_sum_exp(log_model_weight))
  pediatric_mean              <- colSums(model_pediatric_mean * model_weight)
  lowest_overdose_probability <- sum(
    model_lowest_overdose_probability * model_weight
  )
  overdose_probability <- c(lowest_overdose_probability, rep(NA_real_, 4L))
  alpha_mean           <- sum(model_alpha_mean * model_weight)
  distance             <- abs(pediatric_mean - config$target)
  recommended          <- which(distance <= min(distance) + 1e-12)[1L]
  names(model_weight) <- paste0("delta_", format(config$pa_crm_delta, trim  = TRUE))
  list(
    model_weight                = model_weight,
    log_evidence                = log_model_marginal,
    log_joint_evidence          = log_joint_evidence,
    log_adult_power_normalizer  = log_adult_power_normalizer,
    pediatric_mean              = pediatric_mean,
    overdose_probability        = overdose_probability,
    lowest_overdose_probability = lowest_overdose_probability,
    alpha_mean                  = alpha_mean,
    model_alpha_mean            = model_alpha_mean,
    recommended                 = as.integer(recommended),
    child_x                     = standardized$child_x,
    adult_x                     = standardized$adult_x,
    adult_dose_mg_m2            = standardized$adult_dose_mg_m2
  )
}


# =============================================================================
# Pediatric PA-CRM trial simulation
# =============================================================================
simulate_pa_crm_pediatric <- function(
  scenario, replication, potentials,
  adult_trial, config
) {
  method               <- "PA_CRM"
  n_dose               <- 5L
  n                    <- integer(n_dose)
  x                    <- integer(n_dose)
  log_auc_sum          <- numeric(n_dose)
  auc_n                <- integer(n_dose)
  current              <- 1L
  early_stop           <- FALSE
  stop_reason          <- "pediatric_max_n"
  interim_rows         <- list()
  path_rows            <- list()
  diagnostic_rows      <- list()
  max_cohorts          <- config$pediatric_max_n %/% config$cohort_size
  final_adult_snapshot <- NULL
  final_fit            <- NULL

  for (cohort in seq_len(max_cohorts)) {
    dose     <- current
    enrolled <- enroll_pediatric_cohort(
      dose, n, x, log_auc_sum, auc_n, scenario, potentials, config
    )
    n              <- enrolled$n
    x              <- enrolled$x
    log_auc_sum    <- enrolled$log_auc_sum
    auc_n          <- enrolled$auc_n
    adult_snapshot <- adult_snapshot_for_pediatric_interim(
      adult_trial, cohort, config
    )
    final_adult_snapshot <- adult_snapshot
    fit                  <- pa_crm_fit(adult_snapshot$n, adult_snapshot$x, n, x, config)
    final_fit            <- fit
    reached_maximum      <- sum(n) >= config$pediatric_max_n
    safety_stop          <- fit$lowest_overdose_probability > config$pa_crm_safety_cutoff

    if (safety_stop) {
      actual_action <- "STOP_SAFETY"
      next_dose     <- NA_integer_
      early_stop    <- TRUE
      stop_reason   <- "pa_crm_lowest_dose_overdose"
    } else if (reached_maximum) {
      actual_action <- "END_MAX_N"
      next_dose     <- NA_integer_
      stop_reason   <- "pediatric_max_n"
    } else {
      direction     <- sign(fit$recommended - dose)
      next_dose     <- as.integer(clamp(dose + direction, 1L, n_dose))
      actual_action <- if (next_dose > dose) "E" else if (next_dose < dose) "D" else "S"
    }

    interim_rows[[cohort]] <- data.frame(
      config_id              = config$config_id,
      scenario               = scenario$id,
      method                 = method,
      replication            = replication,
      pediatric_cohort       = cohort,
      pediatric_total_n      = sum(n),
      current_dose           = dose,
      current_n              = n[dose],
      current_dlt            = x[dose],
      cohort_dlt             = enrolled$cohort_dlt,
      current_gm_auc         = geometric_mean_from_logs(log_auc_sum, auc_n)[dose],
      adult_requested_cohort = adult_snapshot$requested_adult_cohort,
      adult_available_cohort = adult_snapshot$available_adult_cohort,
      adult_total_n          = adult_snapshot$total_n,
      ordinary_action        = if (fit$recommended > dose) "E" else if (fit$recommended < dose) "D" else "S",
      actual_action          = actual_action,
      next_dose              = next_dose,
      safety_trigger_dose    = if (safety_stop) 1L else NA_integer_,
      early_stop             = early_stop,
      matched_adult_dose     = NA_integer_,
      match_distance         = NA_real_,
      s_a                    = NA_real_,
      w_post                 = NA_real_,
      borrowed_ess_elir      = NA_real_,
      stringsAsFactors       = FALSE
    )
    path_rows[[cohort]] <- data.frame(
      config_id        = config$config_id,
      scenario         = scenario$id,
      method           = method,
      replication      = replication,
      pediatric_cohort = cohort,
      assigned_dose    = dose,
      cohort_dlt       = enrolled$cohort_dlt,
      actual_action    = actual_action,
      next_dose        = next_dose,
      stringsAsFactors = FALSE
    )
    diagnostic_rows[[cohort]] <- data.frame(
      config_id                   = config$config_id,
      scenario                    = scenario$id,
      method                      = method,
      replication                 = replication,
      pediatric_cohort            = cohort,
      current_dose                = dose,
      raw_recommended_dose        = fit$recommended,
      next_dose                   = next_dose,
      alpha_mean                  = fit$alpha_mean,
      delta_0_weight              = fit$model_weight[1L],
      delta_0_5_weight            = fit$model_weight[2L],
      delta_1_weight              = fit$model_weight[3L],
      lowest_overdose_probability = fit$lowest_overdose_probability,
      predicted_p1                = fit$pediatric_mean[1L],
      predicted_p2                = fit$pediatric_mean[2L],
      predicted_p3                = fit$pediatric_mean[3L],
      predicted_p4                = fit$pediatric_mean[4L],
      predicted_p5                = fit$pediatric_mean[5L],
      stringsAsFactors            = FALSE
    )

    if (early_stop || reached_maximum) break
    current <- next_dose
  }

  selected <- if (early_stop) NA_integer_ else final_fit$recommended
  common   <- make_trial_summary_and_dose_tables(
    method, scenario, replication, n, x, log_auc_sum, auc_n,
    selected, stop_reason, early_stop, final_adult_snapshot,
    rep(FALSE, n_dose), config
  )
  common$dose$posterior_median            <- NA_real_
  common$dose$isotonic_estimate           <- NA_real_
  common$dose$final_pava_weight           <- NA_real_
  common$dose$final_borrowed_ess_positive <- NA_real_
  common$dose$final_matched_adult_dose    <- NA_integer_
  c(common, list(
      interim         = rbind_fill(interim_rows),
      path            = rbind_fill(path_rows),
      matching        = data.frame(),
      borrowing       = data.frame(),
      final_borrowing = data.frame(),
      pa_crm          = rbind_fill(diagnostic_rows),
      failure         = data.frame()
  ))
}

