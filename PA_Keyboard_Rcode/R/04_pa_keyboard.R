# =============================================================================
# PA-Keyboard: exposure-guided local adult-information borrowing
# =============================================================================
# This module provides numerical calculations and output helpers for the
# PA-Keyboard workflows defined in run_analysis.R. It fits local exposure-toxicity
# models, computes OVL, limits adult ESS, and performs final MTD selection.
#
# Reading order:
#   1. Configure and validate the local exposure-toxicity model.
#   2. Recover observed AUC/DLT pairs and evaluate the likelihood on a theta grid.
#   3. Normalize pediatric/adult posteriors and integrate their overlap.
#   4. Convert OVL to an operational weight and cap the adult effective sample size.
#   5. Update the beta posterior and select the final dose by weighted PAVA.
#   6. Assemble interim and final diagnostic records.
#
# Execution relationship: run_analysis.R Section 4a calls these numerical
# functions for local borrowing; Section 4b repeats that workflow after each
# pediatric cohort. This module supplies calculations, not the simulation loop.
#
# Function guide:
#   make_pa_keyboard_config(),
#     validate_pa_keyboard_config()        Specify local borrowing parameters.
#   extract_observed_auc_dlt()             Retrieve matched patient-level data.
#   pa_keyboard_loglik_on_theta_grid(),
#     pa_keyboard_normalize_log_density()  Evaluate exposure-toxicity posteriors.
#   pa_keyboard_posterior_ovl()            Calculate posterior compatibility.
#   pa_keyboard_operational_delta()        Apply the prescribed OVL convention.
#   pa_keyboard_effective_posterior()      Update weighted beta posterior counts.
#   make_pa_keyboard_base_config(),
#     make_pa_keyboard_variants(),
#     validate_pa_keyboard_variant_config()
#                                         Define supported borrowing settings.
#   pa_keyboard_apply_ess_cap()            Limit borrowed adult sample size.
#   select_pa_keyboard_mtd()               Apply posterior medians and weighted PAVA.
#   record_pa_keyboard_interim()           Record interim diagnostic tables.
#   assemble_pa_keyboard_trial()           Assemble trial and dose-level results.
#   empty_pa_keyboard_final_fit()          Represent an early-stopped trial.
#   initialize_pa_keyboard_borrowing()      Initialize the pediatric-only result.
#   complete_pa_keyboard_borrowing()        Attach borrowing diagnostics.
# =============================================================================

# =============================================================================
# Exposure-toxicity overlap model and compatibility borrowing
# =============================================================================

# Configure numerical integration, operational OVL and final dose-selection rules.
# Baseline overrides describe the trial; pa_keyboard_overrides describe borrowing.
# The returned configuration is checked before any posterior calculation.
make_pa_keyboard_config <- function(
  baseline_overrides    = list(),
  pa_keyboard_overrides = list()
) {

  if (!"output_dir" %in% names(baseline_overrides)) {
    baseline_overrides$output_dir <- file.path(
      getwd(), "PA_Keyboard_results"
    )
  }

  config <- make_default_config(baseline_overrides)

  overlap_defaults <- list(
    pa_keyboard_prior_id            = "P4_q_beta11",
    pa_keyboard_grid_n              = 4096L,
    pa_keyboard_prior_grid_quantile = 0.999999,
    pa_keyboard_delta_digits        = 2L,
    pa_keyboard_delta_min           = 0.01,
    # Require more than one pediatric cohort before evaluating OVL borrowing.
    pa_keyboard_min_pediatric_n       = 2L * config$cohort_size,
    pa_keyboard_final_min_pediatric_n = 2L * config$cohort_size,
    # Combine pediatric information with the OVL-weighted adult sample size.
    pa_keyboard_final_weight  = "nP_plus_delta_nA"
  )

  unknown <- setdiff(names(pa_keyboard_overrides), names(overlap_defaults))
  if (length(unknown)) {
    stopf("Unknown PA-Keyboard config field(s): %s", paste(unknown, collapse  = ", "))
  }
  overlap_defaults[names(pa_keyboard_overrides)] <- pa_keyboard_overrides
  config[names(overlap_defaults)] <- overlap_defaults

  validate_pa_keyboard_config(config)
}


# Check numerical-grid resolution, prior truncation, borrowing gates and weights.
# Validation rejects unsupported values; it does not estimate model parameters.
validate_pa_keyboard_config <- function(config) {
  # Validate the shared simulation configuration.
  validate_config(config)

  assert_true(
    identical(config$pa_keyboard_prior_id, "P4_q_beta11"),
    "The local exposure-toxicity model requires the reference Beta(1,1) prior."
  )
  assert_true(
    is.numeric(config$pa_keyboard_grid_n) && length(config$pa_keyboard_grid_n) == 1L &&
    config$pa_keyboard_grid_n >= 100L &&
    config$pa_keyboard_grid_n == as.integer(config$pa_keyboard_grid_n),
    "pa_keyboard_grid_n must be an integer >= 100."
  )
  assert_scalar_numeric(
    config$pa_keyboard_prior_grid_quantile,
    "pa_keyboard_prior_grid_quantile", 0.90, 1, TRUE, TRUE
  )
  assert_true(
    is.numeric(config$pa_keyboard_delta_digits) &&
    config$pa_keyboard_delta_digits == as.integer(config$pa_keyboard_delta_digits) &&
    config$pa_keyboard_delta_digits >= 0,
    "pa_keyboard_delta_digits must be a nonnegative integer."
  )
  assert_scalar_numeric(config$pa_keyboard_delta_min, "pa_keyboard_delta_min", 0, 1)

  for (nm in c("pa_keyboard_min_pediatric_n", "pa_keyboard_final_min_pediatric_n")) {
    value <- config[[nm]]
    assert_true(
      is.numeric(value) && length(value) == 1L && value >= config$cohort_size &&
      value == as.integer(value) && value %% config$cohort_size == 0L,
      paste0(nm, " must be a positive cohort-compatible integer.")
    )
  }

  assert_true(
    config$pa_keyboard_final_weight %in% c("nP_plus_delta_nA", "nP_only"),
    "pa_keyboard_final_weight must be 'nP_plus_delta_nA' or 'nP_only'."
  )
  config
}


# Recover patient-level observations at one population/dose combination.
# Only the first n_at_dose potential outcomes are used, matching the patients
# already enrolled at that dose. This reconstruction draws no new random numbers.
# Return binary toxicity, AUC and log(AUC) vectors for the local likelihood.
extract_observed_auc_dlt <- function(
  population  = c("pediatric", "adult"),
  dose,
  n_at_dose,
  scenario,
  potentials,
  config
) {

  population <- match.arg(population)
  dose       <- as.integer(dose)
  n_at_dose  <- as.integer(n_at_dose)

  if (n_at_dose <= 0L) {
    return(list(
        toxicity = integer(),
        auc      = numeric(),
        log_auc  = numeric()
    ))
  }

  idx <- seq_len(n_at_dose)

  if (population == "pediatric") {
    truth     <- scenario$pediatric_dlt[dose]
    auc_truth <- scenario$pediatric_auc[dose]
    dlt_u     <- potentials$pediatric_dlt_u[idx, dose]
    auc_z     <- potentials$pediatric_auc_z[idx, dose]
  } else {
    truth     <- scenario$adult_dlt[dose]
    auc_truth <- scenario$adult_auc[dose]
    dlt_u     <- potentials$adult_dlt_u[idx, dose]
    auc_z     <- potentials$adult_auc_z[idx, dose]
  }

  toxicity <- as.integer(dlt_u < truth)
  log_auc  <- log(auc_truth) + config$auc_sigma_log * auc_z
  auc      <- exp(log_auc)

  list(
    toxicity = toxicity,
    auc      = auc,
    log_auc  = log_auc
  )
}


# Evaluate the Bernoulli log likelihood for q(AUC; theta) = 1 - exp(-theta * AUC).
# DLT observations contribute log(1 - exp(-theta * AUC)); non-DLT observations
# contribute -theta * AUC. The returned vector has one value per theta grid point.
# The model is fitted locally to the current pediatric or matched adult dose.
pa_keyboard_loglik_on_theta_grid <- function(toxicity, auc, theta_grid) {
  if (length(toxicity) != length(auc)) {
    stop("toxicity and AUC vectors have different lengths.", call.  = FALSE)
  }
  if (!length(auc)) return(rep(0, length(theta_grid)))
  if (any(!is.finite(auc)) || any(auc <= 0)) {
    stop("All AUC observations must be finite and >0.", call.  = FALSE)
  }
  toxicity <- as.integer(toxicity)
  if (any(!toxicity %in% c(0L, 1L))) {
    stop("toxicity must contain only 0/1.", call.  = FALSE)
  }

  loglik <- numeric(length(theta_grid))

  # -expm1(-z) evaluates 1 - exp(-z) accurately when z is close to zero.
  auc_tox <- auc[toxicity == 1L]
  if (length(auc_tox)) {
    z      <- outer(theta_grid, auc_tox, `*`)
    loglik <- loglik + rowSums(log(-expm1(-z)))
  }

  auc_notox <- auc[toxicity == 0L]
  if (length(auc_notox)) {
    loglik <- loglik - theta_grid * sum(auc_notox)
  }

  loglik
}


# Convert grid log densities to a density integrating to one by the midpoint rule.
# Subtracting the largest finite log density prevents numerical overflow; the
# constant cancels during normalization. dx is the common grid-cell width.
pa_keyboard_normalize_log_density <- function(log_density, dx) {
  finite <- is.finite(log_density)
  if (!any(finite)) {
    stop("PA-Keyboard posterior has no finite grid values.", call.  = FALSE)
  }
  shift        <- max(log_density[finite])
  dens         <- numeric(length(log_density))
  dens[finite] <- exp(log_density[finite] - shift)
  z            <- sum(dens) * dx
  if (!is.finite(z) || z <= 0) {
    stop("Failed to normalize PA-Keyboard posterior density.", call.  = FALSE)
  }
  dens / z
}

# Posterior overlap of the two dose-local exposure-toxicity models.
# Input: patient-level DLT/AUC for the current pediatric and matched adult dose.
# Output: overlap weight, reference exposure A*, and arithmetic mean exposures.
pa_keyboard_posterior_ovl <- function(
  toxicity_p,
  auc_p,
  toxicity_a,
  auc_a,
  config
) {

  if (!length(auc_p) || !length(auc_a)) {
    return(list(
        delta_raw  = 0,
        a_star     = NA_real_,
        mean_auc_p = if (length(auc_p)) mean(auc_p) else NA_real_,
        mean_auc_a = if (length(auc_a)) mean(auc_a) else NA_real_
    ))
  }

  if (length(toxicity_p) != length(auc_p) ||
    length(toxicity_a) != length(auc_a)) {
    stop("PA-Keyboard toxicity/AUC vector lengths are inconsistent.", call.  = FALSE)
  }

  # Define one common reference exposure from the two arithmetic mean AUCs.
  # These means differ from the geometric means used for dose matching.
  mean_auc_p <- mean(auc_p)
  mean_auc_a <- mean(auc_a)
  a_star     <- (mean_auc_p + mean_auc_a) / 2

  if (!is.finite(a_star) || a_star <= 0) {
    stop("Invalid PA-Keyboard common reference exposure A*.", call.  = FALSE)
  }

  # A common q(A*) ~ Beta(1,1) prior corresponds to theta ~ Exp(rate=A*).
  theta_max <- stats::qexp(
    config$pa_keyboard_prior_grid_quantile,
    rate  = a_star
  )
  # Use equal-width midpoint cells on the prior-quantile-truncated theta range.
  n_grid     <- as.integer(config$pa_keyboard_grid_n)
  dx         <- theta_max / n_grid
  theta_grid <- (seq_len(n_grid) - 0.5) * dx

  logprior <- stats::dexp(theta_grid, rate  = a_star, log  = TRUE)

  # Fit each population separately on the same grid with the same prior.
  dens_p <- pa_keyboard_normalize_log_density(
    pa_keyboard_loglik_on_theta_grid(toxicity_p, auc_p, theta_grid) + logprior,
    dx
  )
  dens_a <- pa_keyboard_normalize_log_density(
    pa_keyboard_loglik_on_theta_grid(toxicity_a, auc_a, theta_grid) + logprior,
    dx
  )

  # Integrate the pointwise smaller density: OVL = integral min(fP, fA) dtheta.
  # The clamp removes numerical excursions outside the probability range.
  ovl <- sum(pmin(dens_p, dens_a)) * dx
  ovl <- max(0, min(1, ovl))

  list(
    delta_raw  = ovl,
    a_star     = a_star,
    mean_auc_p = mean_auc_p,
    mean_auc_a = mean_auc_a
  )
}


# Convert numerical overlap to the weight used by the dose-finding algorithm.
# Round first, then set interim weights at or below the configured minimum to zero.
# Final assessment uses the rounded overlap without the interim cutoff. Both
# paths restrict the result to [0, 1]; the adult ESS cap is applied separately.
pa_keyboard_operational_delta <- function(delta_raw, config, final  = FALSE) {
  delta <- round(delta_raw, digits  = config$pa_keyboard_delta_digits)
  # Interim borrowing applies the configured minimum; final OVL is only rounded.
  if (!isTRUE(final) && delta <= config$pa_keyboard_delta_min) delta <- 0
  max(0, min(1, delta))
}


# Update a Beta(1, 1) prior using pediatric counts plus weighted adult counts.
# Effective DLTs are xP + delta*xA; effective total count is nP + delta*nA.
# The posterior is Beta(1 + x_eff, 1 + n_eff - x_eff). Its one-component list
# uses the interface expected by posterior_key_probabilities() in 02_keyboard.R.
# With delta = 0, the result is exactly the pediatric-only beta posterior.
pa_keyboard_effective_posterior <- function(
  pediatric_x,
  pediatric_n,
  adult_x = 0,
  adult_n = 0,
  delta   = 0
) {

  x_eff <- pediatric_x + delta * adult_x
  n_eff <- pediatric_n + delta * adult_n
  list(
    weights           = 1,
    a                 = 1 + x_eff,
    b                 = 1 + n_eff - x_eff,
    x_eff             = x_eff,
    n_eff             = n_eff,
    effective_adult_n = delta * adult_n
  )
}


# =============================================================================
# Borrowing configurations, adult ESS cap and final dose selection
# =============================================================================

# Combine baseline trial settings with the exposure-toxicity borrowing settings.
# This configuration is the starting point for the prespecified borrowing variants.
make_pa_keyboard_base_config <- function(
  baseline_overrides = list(),
  output_dir         = NULL
) {

  if (is.null(output_dir)) {
    output_dir <- file.path(
      getwd(),
      "PA_Keyboard_results"
    )
  }

  if (!"output_dir" %in% names(baseline_overrides)) {
    baseline_overrides$output_dir <- output_dir
  }

  cfg <- make_pa_keyboard_config(
    baseline_overrides    = baseline_overrides,
    pa_keyboard_overrides = list(
      pa_keyboard_min_pediatric_n       = 6L,
      pa_keyboard_final_min_pediatric_n = 6L,
      pa_keyboard_final_weight          = "nP_plus_delta_nA"
    )
  )

  # Add the ET-OVL effective-sample-size fields.
  cfg$pa_keyboard_ess_cap      <- Inf
  cfg$pa_keyboard_diag_variant <- "PAK_original"
  cfg$pa_keyboard_diag_label   <- PA_KEYBOARD_VARIANT_LABELS[["PAK_original"]]
  cfg$config_id                <- "PAK_original"
  cfg
}


# Construct borrowing variants from a common base configuration.
# Each variant sets its pediatric sample-size gate, adult ESS cap and output labels.
# The formal-study setup selects the capped variant before setting study metadata.
make_pa_keyboard_variants <- function(base_config) {
  base_config <- validate_pa_keyboard_variant_config(base_config)

  original <- base_config
  original$pa_keyboard_min_pediatric_n <- 6L
  original$pa_keyboard_final_min_pediatric_n <- 6L
  original$pa_keyboard_ess_cap <- Inf
  original$pa_keyboard_diag_variant <- "PAK_original"
  original$pa_keyboard_diag_label <- PA_KEYBOARD_VARIANT_LABELS[["PAK_original"]]
  original$config_id <- "PAK_original"

  early                                   <- base_config
  early$pa_keyboard_min_pediatric_n       <- 3L
  early$pa_keyboard_final_min_pediatric_n <- 3L
  early$pa_keyboard_ess_cap               <- Inf
  early$pa_keyboard_diag_variant          <- "PAK_early"
  early$pa_keyboard_diag_label            <- PA_KEYBOARD_VARIANT_LABELS[["PAK_early"]]
  early$config_id                         <- "PAK_early"

  capped                                   <- base_config
  capped$pa_keyboard_min_pediatric_n       <- 6L
  capped$pa_keyboard_final_min_pediatric_n <- 6L
  capped$pa_keyboard_ess_cap               <- 6
  capped$pa_keyboard_diag_variant          <- "PAK_capped"
  capped$pa_keyboard_diag_label            <- PA_KEYBOARD_VARIANT_LABELS[["PAK_capped"]]
  capped$config_id                         <- "PAK_capped"

  out <- list(
    PAK_original = validate_pa_keyboard_variant_config(original),
    PAK_early    = validate_pa_keyboard_variant_config(early),
    PAK_capped   = validate_pa_keyboard_variant_config(capped)
  )
  out[PA_KEYBOARD_VARIANT_ORDER]
}


# Validate the general PA-Keyboard settings and the adult ESS cap.
# A finite cap must be positive; an infinite cap represents uncapped borrowing.
validate_pa_keyboard_variant_config <- function(config) {
  # Validate shared settings before checking variant-specific fields.
  known                          <- config
  known$pa_keyboard_ess_cap      <- NULL
  known$pa_keyboard_diag_variant <- NULL
  known$pa_keyboard_diag_label   <- NULL
  known                          <- validate_pa_keyboard_config(known)

  assert_true(
    is.numeric(config$pa_keyboard_ess_cap) && length(config$pa_keyboard_ess_cap) == 1L &&
    (is.infinite(config$pa_keyboard_ess_cap) || config$pa_keyboard_ess_cap > 0),
    "pa_keyboard_ess_cap must be a positive number or Inf."
  )
  assert_true(
    config$pa_keyboard_diag_variant %in% PA_KEYBOARD_VARIANT_ORDER,
    "Unknown PA-Keyboard diagnostic variant."
  )
  config
}


# Limit the realized adult contribution while preserving the observed adult DLT rate.
# Uncapped adult ESS is delta_ovl * adult_n. When it exceeds cap, replace the
# weight by cap/adult_n. Return both weights/ESS values and the cap indicator.
# No adult patients or a zero overlap weight implies zero adult contribution.
pa_keyboard_apply_ess_cap <- function(
  delta_ovl,
  adult_n,
  cap  = Inf
) {

  delta_ovl <- max(0, min(1, as.numeric(delta_ovl)))
  adult_n   <- as.numeric(adult_n)

  if (!is.finite(adult_n) || adult_n <= 0 || delta_ovl <= 0) {
    return(list(
        delta_applied              = 0,
        effective_adult_n_uncapped = 0,
        effective_adult_n          = 0,
        cap_triggered              = FALSE
    ))
  }

  eff_uncapped <- delta_ovl * adult_n

  if (is.infinite(cap) || eff_uncapped <= cap + 1e-12) {
    return(list(
        delta_applied              = delta_ovl,
        effective_adult_n_uncapped = eff_uncapped,
        effective_adult_n          = eff_uncapped,
        cap_triggered              = FALSE
    ))
  }

  delta_applied <- cap / adult_n
  delta_applied <- max(0, min(delta_ovl, delta_applied))

  list(
    delta_applied              = delta_applied,
    effective_adult_n_uncapped = eff_uncapped,
    effective_adult_n          = delta_applied * adult_n,
    cap_triggered              = TRUE
  )
}


# Compute final dose estimates and select among tried, non-eliminated doses.
# For each eligible dose, call the Section 4a borrowing workflow with final = TRUE,
# then calculate the median of its resulting beta posterior. Thus final estimates
# follow the implemented borrowing rule; they are not automatically pediatric-only.
# The configuration chooses PAVA weights: nP_only uses pediatric sample size;
# nP_plus_delta_nA also includes the realized adult ESS. Fit nondecreasing estimates
# and select the dose closest to target using the shared dose tie-breaking rule.
# Return the selected index, dose estimates, weights and final borrowing diagnostics.
select_pa_keyboard_mtd <- function(
  n,
  x,
  log_auc_sum,
  auc_n,
  eliminated,
  adult_snapshot,
  scenario,
  potentials,
  config,
  replication
) {

  n_dose   <- length(n)
  tried    <- which(n > 0L)
  eligible <- tried[!eliminated[tried]]

  posterior_median           <- rep(NA_real_, n_dose)
  isotonic                   <- rep(NA_real_, n_dose)
  pava_weight                <- rep(NA_real_, n_dose)
  delta_raw                  <- rep(0, n_dose)
  delta_ovl                  <- rep(0, n_dose)
  delta_applied              <- rep(0, n_dose)
  matched                    <- rep(NA_integer_, n_dose)
  effective_adult_n_uncapped <- rep(0, n_dose)
  effective_adult_n          <- rep(0, n_dose)
  cap_triggered              <- rep(FALSE, n_dose)
  diagnostics                <- vector("list", n_dose)

  if (!length(eligible)) {
    return(list(
        selected                   = NA_integer_,
        posterior_median           = posterior_median,
        isotonic                   = isotonic,
        pava_weight                = pava_weight,
        delta_raw                  = delta_raw,
        delta_ovl                  = delta_ovl,
        delta_applied              = delta_applied,
        matched_adult_dose         = matched,
        effective_adult_n_uncapped = effective_adult_n_uncapped,
        effective_adult_n          = effective_adult_n,
        cap_triggered              = cap_triggered,
        diagnostics                = data.frame()
    ))
  }

  gm_auc <- geometric_mean_from_logs(log_auc_sum, auc_n)

  # Final estimates are computed only for visited doses that remain admissible.
  for (dose in eligible) {
    # Match AUC, evaluate overlap, cap ESS, and update the local posterior.
    borrow <- pa_keyboard_borrow_at_dose(
      pediatric_dose   = dose,
      pediatric_n      = n[dose],
      pediatric_x      = x[dose],
      pediatric_gm_auc = gm_auc[dose],
      adult_snapshot   = adult_snapshot,
      scenario         = scenario,
      potentials       = potentials,
      config           = config,
      final            = TRUE
    )

    matched[dose]                    <- borrow$matched_adult_dose
    delta_raw[dose]                  <- borrow$delta_raw
    delta_ovl[dose]                  <- borrow$delta_ovl
    delta_applied[dose]              <- borrow$delta_applied
    effective_adult_n_uncapped[dose] <- borrow$effective_adult_n_uncapped
    effective_adult_n[dose]          <- borrow$effective_adult_n
    cap_triggered[dose]              <- borrow$cap_triggered

    posterior_median[dose] <- stats::qbeta(
      0.5,
      shape1 = borrow$posterior$a,
      shape2 = borrow$posterior$b
    )

    # The PAVA weight option affects isotonic weighting, not the posterior median.
    pava_weight[dose] <- if (config$pa_keyboard_final_weight == "nP_only") {
      n[dose]
    } else {
      n[dose] + effective_adult_n[dose]
    }

    diagnostics[[dose]] <- data.frame(
      config_id          = config$config_id,
      variant            = config$pa_keyboard_diag_variant,
      variant_label      = config$pa_keyboard_diag_label,
      scenario           = scenario$id,
      method             = PA_KEYBOARD_METHOD,
      replication        = replication,
      pediatric_dose     = dose,
      eligible           = TRUE,
      pediatric_n        = n[dose],
      pediatric_x        = x[dose],
      pediatric_gm_auc   = gm_auc[dose],
      match_available    = isTRUE(borrow$match$available),
      ovl_evaluated      = isTRUE(borrow$ovl_evaluated),
      matched_adult_dose = borrow$matched_adult_dose,
      match_distance     = if (isTRUE(borrow$match$available)) {
        borrow$match$distance
      } else {
        NA_real_
      },
      adult_n  = if (!is.na(borrow$matched_adult_dose)) {
        adult_snapshot$n[borrow$matched_adult_dose]
      } else {
        0L
      },
      adult_x  = if (!is.na(borrow$matched_adult_dose)) {
        adult_snapshot$x[borrow$matched_adult_dose]
      } else {
        0L
      },
      delta_raw                  = borrow$delta_raw,
      delta_ovl                  = borrow$delta_ovl,
      delta_applied              = borrow$delta_applied,
      effective_adult_n_uncapped = borrow$effective_adult_n_uncapped,
      effective_adult_n          = borrow$effective_adult_n,
      cap_triggered              = borrow$cap_triggered,
      a_star                     = borrow$a_star,
      posterior_median           = posterior_median[dose],
      pava_weight                = pava_weight[dose],
      stringsAsFactors           = FALSE
    )
  }

  # Enforce nondecreasing toxicity across eligible doses before target selection.
  isotonic[eligible] <- pava_weighted(
    posterior_median[eligible],
    weights  = pava_weight[eligible]
  )

  selected <- select_closest_target_dose(
    eligible,
    isotonic,
    config$target,
    tolerance  = 1e-12
  )

  diagnostics <- rbind_fill(diagnostics)
  if (nrow(diagnostics)) {
    diagnostics$selected          <- diagnostics$pediatric_dose == selected
    diagnostics$isotonic_estimate <- isotonic[diagnostics$pediatric_dose]
  }

  list(
    selected                   = as.integer(selected),
    posterior_median           = posterior_median,
    isotonic                   = isotonic,
    pava_weight                = pava_weight,
    delta_raw                  = delta_raw,
    delta_ovl                  = delta_ovl,
    delta_applied              = delta_applied,
    matched_adult_dose         = matched,
    effective_adult_n_uncapped = effective_adult_n_uncapped,
    effective_adult_n          = effective_adult_n,
    cap_triggered              = cap_triggered,
    diagnostics                = diagnostics
  )
}

# Collect diagnostics for one completed pediatric interim. This helper records
# an already computed decision; it does not update the posterior or assign doses.
record_pa_keyboard_interim <- function(
  scenario, config, dose, borrow, method, replication, cohort, child_gm,
  adult_snapshot, n, x, decision, actual_action, next_dose, safety, early_stop,
  enrolled
) {
  oracle_global <- oracle_auc_match(
    scenario$pediatric_auc[dose],
    scenario$adult_auc,
    config$adult_doses_mg,
    config$auc_caliper_log,
    config$auc_tie_tolerance
  )
  oracle_eligible <- oracle_auc_match(
    scenario$pediatric_auc[dose],
    scenario$adult_auc,
    config$adult_doses_mg,
    config$auc_caliper_log,
    config$auc_tie_tolerance,
    eligible  = borrow$match$eligible
  )

  matching <- data.frame(
    config_id          = config$config_id,
    variant            = config$pa_keyboard_diag_variant,
    variant_label      = config$pa_keyboard_diag_label,
    scenario           = scenario$id,
    method             = method,
    replication        = replication,
    pediatric_cohort   = cohort,
    pediatric_dose     = dose,
    pediatric_gm_auc   = child_gm[dose],
    match_available    = isTRUE(borrow$match$available),
    ovl_evaluated      = isTRUE(borrow$ovl_evaluated),
    matched_adult_dose = borrow$matched_adult_dose,
    match_distance     = if (isTRUE(borrow$match$available)) {
      borrow$match$distance
    } else {
      NA_real_
    },
    oracle_global                  = oracle_global,
    oracle_eligible                = oracle_eligible,
    correct_given_oracle_available = if (!is.na(oracle_eligible)) {
      identical(
        as.integer(borrow$matched_adult_dose),
        as.integer(oracle_eligible)
      )
    } else {
      NA
    },
    stringsAsFactors  = FALSE
  )

  # Record whether the posterior overlap coefficient was evaluated.
  borrowing <- data.frame(
    config_id          = config$config_id,
    variant            = config$pa_keyboard_diag_variant,
    variant_label      = config$pa_keyboard_diag_label,
    scenario           = scenario$id,
    method             = method,
    replication        = replication,
    pediatric_cohort   = cohort,
    pediatric_dose     = dose,
    pediatric_n        = n[dose],
    pediatric_x        = x[dose],
    match_available    = isTRUE(borrow$match$available),
    ovl_evaluated      = isTRUE(borrow$ovl_evaluated),
    matched_adult_dose = borrow$matched_adult_dose,
    adult_n            = if (!is.na(borrow$matched_adult_dose)) {
      adult_snapshot$n[borrow$matched_adult_dose]
    } else {
      0L
    },
    adult_x  = if (!is.na(borrow$matched_adult_dose)) {
      adult_snapshot$x[borrow$matched_adult_dose]
    } else {
      0L
    },
    delta_raw                  = borrow$delta_raw,
    delta_ovl                  = borrow$delta_ovl,
    delta_applied              = borrow$delta_applied,
    effective_adult_n_uncapped = borrow$effective_adult_n_uncapped,
    effective_adult_n          = borrow$effective_adult_n,
    cap_triggered              = borrow$cap_triggered,
    a_star                     = borrow$a_star,
    mean_auc_p_arith           = borrow$mean_auc_p_arith,
    mean_auc_a_arith           = borrow$mean_auc_a_arith,
    stringsAsFactors           = FALSE
  )

  interim <- data.frame(
    config_id              = config$config_id,
    variant                = config$pa_keyboard_diag_variant,
    variant_label          = config$pa_keyboard_diag_label,
    scenario               = scenario$id,
    method                 = method,
    replication            = replication,
    pediatric_cohort       = cohort,
    current_dose           = dose,
    pediatric_n_at_dose    = n[dose],
    pediatric_x_at_dose    = x[dose],
    adult_total_n          = adult_snapshot$total_n,
    adult_requested_cohort = adult_snapshot$requested_adult_cohort,
    adult_available_cohort = adult_snapshot$available_adult_cohort,
    matched_adult_dose     = borrow$matched_adult_dose,
    match_distance         = if (isTRUE(borrow$match$available)) {
      borrow$match$distance
    } else {
      NA_real_
    },
    ovl_evaluated              = isTRUE(borrow$ovl_evaluated),
    delta_ovl                  = borrow$delta_ovl,
    delta_applied              = borrow$delta_applied,
    effective_adult_n_uncapped = borrow$effective_adult_n_uncapped,
    effective_adult_n          = borrow$effective_adult_n,
    cap_triggered              = borrow$cap_triggered,
    raw_action                 = decision$action,
    actual_action              = actual_action,
    next_dose                  = next_dose,
    safety_trigger_dose        = safety$trigger_dose,
    early_stop                 = early_stop,
    stringsAsFactors           = FALSE
  )

  path <- data.frame(
    config_id        = config$config_id,
    variant          = config$pa_keyboard_diag_variant,
    variant_label    = config$pa_keyboard_diag_label,
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

  list(matching  = matching, borrowing  = borrowing, interim  = interim, path  = path)
}


# Combine final dose estimates with the trial summary and recorded interim tables.
# The selected dose is missing for early-stopped trials. Output fields follow the
# common schema consumed by the operating-characteristic functions in 06_results.R.
assemble_pa_keyboard_trial <- function(
  scenario, replication, config, method, n, x, log_auc_sum, auc_n, early_stop,
  stop_reason, final_adult_snapshot, eliminated, final_fit, interim_rows, path_rows,
  matching_rows, borrowing_rows
) {
  selected <- if (early_stop) NA_integer_ else final_fit$selected

  common <- make_trial_summary_and_dose_tables(
    method               = method,
    scenario             = scenario,
    replication          = replication,
    n                    = n,
    x                    = x,
    log_auc_sum          = log_auc_sum,
    auc_n                = auc_n,
    selected_mtd         = selected,
    stop_reason          = stop_reason,
    early_stop           = early_stop,
    final_adult_snapshot = final_adult_snapshot,
    eliminated           = eliminated,
    config               = config
  )

  # Label trial and dose records with the analysis variant.
  common$trial$variant       <- config$pa_keyboard_diag_variant
  common$trial$variant_label <- config$pa_keyboard_diag_label
  common$dose$variant        <- config$pa_keyboard_diag_variant
  common$dose$variant_label  <- config$pa_keyboard_diag_label

  common$dose$posterior_median <- final_fit$posterior_median
  common$dose$isotonic_estimate <- final_fit$isotonic
  common$dose$final_pava_weight <- final_fit$pava_weight
  common$dose$final_borrowed_ess_positive <- NA_real_
  common$dose$final_matched_adult_dose <- final_fit$matched_adult_dose
  common$dose$final_pa_keyboard_delta_raw <- final_fit$delta_raw
  common$dose$final_pa_keyboard_delta_ovl <- final_fit$delta_ovl
  common$dose$final_pa_keyboard_delta_applied <- final_fit$delta_applied
  common$dose$final_pa_keyboard_effective_adult_n_uncapped <-
  final_fit$effective_adult_n_uncapped
  common$dose$final_pa_keyboard_effective_adult_n <- final_fit$effective_adult_n
  common$dose$final_pa_keyboard_cap_triggered <- final_fit$cap_triggered

  c(
    common,
    list(
      interim         = rbind_fill(interim_rows),
      path            = rbind_fill(path_rows),
      matching        = rbind_fill(matching_rows),
      borrowing       = rbind_fill(borrowing_rows),
      final_borrowing = final_fit$diagnostics,
      pa_crm          = data.frame(),
      failure         = data.frame()
    )
  )
}


# Represent a trial stopped without an MTD recommendation.
# Dose estimates are missing, borrowing contributions are zero, and the final
# borrowing table is empty; vector lengths match the pediatric dose grid.
empty_pa_keyboard_final_fit <- function(n_dose) {
  list(
    selected                   = NA_integer_,
    posterior_median           = rep(NA_real_, n_dose),
    isotonic                   = rep(NA_real_, n_dose),
    pava_weight                = rep(NA_real_, n_dose),
    delta_raw                  = rep(0, n_dose),
    delta_ovl                  = rep(0, n_dose),
    delta_applied              = rep(0, n_dose),
    matched_adult_dose         = rep(NA_integer_, n_dose),
    effective_adult_n_uncapped = rep(0, n_dose),
    effective_adult_n          = rep(0, n_dose),
    cap_triggered              = rep(FALSE, n_dose),
    diagnostics                = data.frame()
  )
}


# Create the default result before borrowing eligibility is checked.
# Record the AUC match but initialize all adult contributions to zero and retain
# the pediatric-only posterior. Failed eligibility checks return this object.
initialize_pa_keyboard_borrowing <- function(match_result, pediatric_x, pediatric_n) {
  list(
    match              = match_result,
    matched_adult_dose = if (isTRUE(match_result$available)) {
      as.integer(match_result$match)
    } else {
      NA_integer_
    },
    ovl_evaluated              = FALSE,
    delta_raw                  = 0,
    delta_ovl                  = 0,
    delta_applied              = 0,
    effective_adult_n_uncapped = 0,
    effective_adult_n          = 0,
    cap_triggered              = FALSE,
    a_star                     = NA_real_,
    mean_auc_p_arith           = NA_real_,
    mean_auc_a_arith           = NA_real_,
    posterior                  = child_only_posterior(pediatric_x, pediatric_n)
  )

}


# Attach evaluated overlap, capped weights, adult ESS and the updated posterior.
# This helper only packages computed quantities; it does not recompute borrowing.
complete_pa_keyboard_borrowing <- function(out, fit, delta_ovl, capped, posterior) {
  out$ovl_evaluated              <- TRUE
  out$delta_raw                  <- fit$delta_raw
  out$delta_ovl                  <- delta_ovl
  out$delta_applied              <- capped$delta_applied
  out$effective_adult_n_uncapped <- capped$effective_adult_n_uncapped
  out$effective_adult_n          <- capped$effective_adult_n
  out$cap_triggered              <- capped$cap_triggered
  out$a_star                     <- fit$a_star
  out$mean_auc_p_arith           <- fit$mean_auc_p
  out$mean_auc_a_arith           <- fit$mean_auc_a
  out$posterior                  <- posterior
  out
}

