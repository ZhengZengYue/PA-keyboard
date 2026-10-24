# =============================================================================
# Original Keyboard design, dose matching, and reference trial
# =============================================================================
# This module implements the original Keyboard decision rules, pediatric-only
# safety elimination, AUC dose matching, and simulations for the adult reference
# trial and pediatric-only Keyboard comparator.
#
# Function guide:
#   make_keyboard_keys()              Partition the toxicity probability interval.
#   beta_key_probabilities(),
#     posterior_key_probabilities()   Calculate probabilities for each key.
#   choose_strongest_key(),
#     keyboard_decision()             Determine escalate/stay/de-escalate decisions.
#   child_only_posterior()            Update the pediatric-only beta posterior.
#   pediatric_overdose_probability(),
#     update_keyboard_elimination()   Apply the pediatric safety rule.
#   pava_weighted(),
#     select_closest_target_dose(),
#     select_keyboard_mtd()           Select the MTD using isotonic estimates.
#   next_keyboard_dose()              Apply one-level dose movement restrictions.
#   choose_auc_candidate(),
#     select_auc_match(),
#     oracle_auc_match()              Identify observed and truth-level AUC matches.
#   simulate_adult_trial()            Generate the concurrent adult Keyboard trial.
#   adult_snapshot_for_pediatric_interim()
#                                    Retrieve concurrent adult data; exported as
#                                    BASELINE_ADULT_SNAPSHOT_FOR_PEDIATRIC_INTERIM
#                                    before the sensitivity wrapper is loaded.
#   simulate_keyboard_pediatric()     Simulate the pediatric-only comparator.
# =============================================================================


# Construct target-width intervals on [0, 1], truncating intervals at the boundaries.
# relation = -1, 0, 1 marks below-target, target and above-target keys, respectively.
make_keyboard_keys <- function(config) {
  cache_key <- sprintf("keyboard_keys_%.16g_%.16g", config$target_lower,
    config$target_upper)
  if (exists(cache_key, envir  = .PAK_NUMERIC_CACHE, inherits  = FALSE)) {
    return(get(cache_key, envir  = .PAK_NUMERIC_CACHE, inherits  = FALSE))
  }
  width       <- config$target_upper - config$target_lower
  lower_parts <- list()
  upper_edge  <- config$target_lower
  while (upper_edge > 1e-14) {
    lower_edge                              <- max(0, upper_edge - width)
    lower_parts[[length(lower_parts) + 1L]] <- data.frame(
      lower            = lower_edge,
      upper            = upper_edge,
      relation         = -1L,
      stringsAsFactors = FALSE
    )
    upper_edge <- lower_edge
  }
  if (length(lower_parts)) lower_parts <- rev(lower_parts)
  target_part <- data.frame(
    lower            = config$target_lower,
    upper            = config$target_upper,
    relation         = 0L,
    stringsAsFactors = FALSE
  )
  upper_parts <- list()
  lower_edge  <- config$target_upper
  while (lower_edge < 1 - 1e-14) {
    upper_edge                              <- min(1, lower_edge + width)
    upper_parts[[length(upper_parts) + 1L]] <- data.frame(
      lower            = lower_edge,
      upper            = upper_edge,
      relation         = 1L,
      stringsAsFactors = FALSE
    )
    lower_edge <- upper_edge
  }
  keys           <- do.call(rbind, c(lower_parts, list(target_part), upper_parts))
  keys$key_id    <- seq_len(nrow(keys))
  keys$midpoint  <- (keys$lower + keys$upper) / 2
  keys$is_target <- keys$relation == 0L
  rownames(keys) <- NULL
  assign(cache_key, keys, envir  = .PAK_NUMERIC_CACHE)
  keys
}


# Integrate one beta posterior over each key using differences of its CDF.
beta_key_probabilities <- function(a, b, keys) {
  stats::pbeta(keys$upper, shape1  = a, shape2  = b) -
  stats::pbeta(keys$lower, shape1  = a, shape2  = b)
}


# Average key probabilities across the supplied beta-mixture components.
# A pediatric-only or effective-count posterior is represented by one component.
posterior_key_probabilities <- function(posterior, keys) {
  assert_true(length(posterior$weights) == length(posterior$a) &&
    length(posterior$a) == length(posterior$b),
    "Posterior mixture components are inconsistent.")
  probabilities <- numeric(nrow(keys))
  for (k in seq_along(posterior$weights)) {
    probabilities <- probabilities + posterior$weights[k] *
    beta_key_probabilities(posterior$a[k], posterior$b[k], keys)
  }
  probabilities
}


# Choose the key with greatest posterior mass. Ties prefer the target key, then
# the key closest to it, then the higher-toxicity key. E/D/S denote escalate,
# de-escalate and stay; the rightmost tie rule therefore favors caution.
choose_strongest_key <- function(probabilities, keys, tolerance  = 1e-14) {
  max_probability <- max(probabilities)
  tied            <- which(abs(probabilities - max_probability) <= tolerance)
  if (length(tied) > 1L) {
    target_tie <- tied[keys$is_target[tied]]
    if (length(target_tie)) {
      chosen <- target_tie[1L]
    } else {
      distance_to_target <- abs(keys$midpoint[tied] - mean(c(
            keys$lower[keys$is_target], keys$upper[keys$is_target]
      )))
      closest <- tied[abs(distance_to_target - min(distance_to_target)) <= tolerance]
      # An exact remaining tie is resolved conservatively toward the right key.
      chosen <- closest[which.max(keys$relation[closest])]
    }
  } else {
    chosen <- tied
  }
  relation <- keys$relation[chosen]
  action   <- if (relation < 0L) "E" else if (relation > 0L) "D" else "S"
  list(
    key_id       = keys$key_id[chosen],
    relation     = relation,
    action       = action,
    probability  = probabilities[chosen],
    tied_key_ids = paste(keys$key_id[tied], collapse  = ",")
  )
}


# Update a Beta(1, 1) prior with x DLTs among n pediatric patients.
child_only_posterior <- function(x, n) {
  list(weights  = 1, a  = 1 + x, b  = 1 + n - x)
}


# Translate the supplied posterior into a key-based action. If posterior is NULL,
# use pediatric counts alone. Safety elimination is applied separately.
keyboard_decision <- function(x, n, config, posterior  = NULL) {
  if (is.null(posterior)) posterior <- child_only_posterior(x, n)
  keys                       <- make_keyboard_keys(config)
  probabilities              <- posterior_key_probabilities(posterior, keys)
  decision                   <- choose_strongest_key(probabilities, keys)
  decision$key_probabilities <- probabilities
  decision
}


# Compute Pr(p > target | x, n) under the pediatric-only Beta(1+x, 1+n-x) posterior.
pediatric_overdose_probability <- function(x, n, target  = 0.30) {
  1 - stats::pbeta(target, shape1  = 1 + x, shape2  = 1 + n - x)
}


# Evaluate the safety rule only at doses with at least safety_min_n patients.
# Eliminate the lowest triggering dose and all higher doses permanently. The trial
# stops if the lowest dose is eliminated. Inputs contain no borrowed adult counts.
update_keyboard_elimination <- function(n, x, eliminated, config) {
  n_dose        <- length(n)
  probabilities <- vapply(seq_len(n_dose), function(j) {
      if (n[j] < config$safety_min_n) return(NA_real_)
      pediatric_overdose_probability(x[j], n[j], config$target)
    }, numeric(1))
  trigger_candidates <- which(
    !eliminated & n >= config$safety_min_n &
    is.finite(probabilities) & probabilities > config$safety_cutoff
  )
  trigger <- if (length(trigger_candidates)) min(trigger_candidates) else NA_integer_
  if (!is.na(trigger)) eliminated[seq.int(trigger, n_dose)] <- TRUE
  list(
    eliminated           = eliminated,
    trigger_dose         = trigger,
    stop                 = isTRUE(eliminated[1L]),
    overdose_probability = probabilities
  )
}


# Pool adjacent blocks whenever estimated toxicity decreases with increasing dose.
# Replace a violating pair by its weighted mean, then recheck the preceding block.
# Expand the final block levels back to a fitted value for every input dose.
pava_weighted <- function(y, weights  = rep(1, length(y)), tolerance  = 1e-14) {
  if (length(y) != length(weights) || !length(y)) {
    stop("y and weights must be nonempty vectors of equal length.", call.  = FALSE)
  }
  if (any(!is.finite(y)) || any(!is.finite(weights)) || any(weights <= 0)) {
    stop("PAVA inputs must be finite and weights must be positive.", call.  = FALSE)
  }
  levels        <- as.numeric(y)
  block_weights <- as.numeric(weights)
  starts        <- seq_along(y)
  ends          <- seq_along(y)
  i             <- 1L
  while (i < length(levels)) {
    if (levels[i] <= levels[i + 1L] + tolerance) {
      i <- i + 1L
    } else {
      merged_weight <- block_weights[i] + block_weights[i + 1L]
      levels[i]     <- (levels[i] * block_weights[i] +
        levels[i + 1L] * block_weights[i + 1L]) / merged_weight
      block_weights[i] <- merged_weight
      ends[i]          <- ends[i + 1L]
      remove           <- i + 1L
      levels           <- levels[-remove]
      block_weights    <- block_weights[-remove]
      starts           <- starts[-remove]
      ends             <- ends[-remove]
      i                <- max(1L, i - 1L)
    }
  }
  fitted <- numeric(length(y))
  for (k in seq_along(levels)) fitted[starts[k]:ends[k]] <- levels[k]
  fitted
}


# Select among eligible doses by absolute distance from the target. For a tie,
# use the higher dose when the first tied estimate is below target; otherwise
# use the lower dose. Return NA when no eligible dose remains.
select_closest_target_dose <- function(
  eligible, isotonic, target,
  tolerance  = 1e-12
) {
  if (!length(eligible)) return(NA_integer_)
  distance <- abs(isotonic[eligible] - target)
  tied     <- eligible[abs(distance - min(distance)) <= tolerance]
  if (length(tied) == 1L) return(as.integer(tied))
  tied_estimate <- isotonic[tied[1L]]
  if (tied_estimate < target - tolerance) {
    as.integer(max(tied))
  } else {
    # Above or exactly at target, resolve toward the lower pediatric dose.
    as.integer(min(tied))
  }
}


# Calculate pediatric-only beta posterior medians at tried doses and fit weighted
# PAVA with pediatric sample sizes. Select only among tried, non-eliminated doses.
select_keyboard_mtd <- function(n, x, eliminated, config) {
  tried     <- which(n > 0L)
  eligible  <- tried[!eliminated[tried]]
  estimates <- rep(NA_real_, length(n))
  isotonic  <- rep(NA_real_, length(n))
  if (!length(eligible)) {
    return(list(
        selected              = NA_integer_,
        posterior_median      = estimates,
        isotonic              = isotonic,
        pava_weight           = ifelse(n > 0L, as.numeric(n), NA_real_),
        borrowed_ess_positive = ifelse(n > 0L, 0, NA_real_),
        matched_adult_dose    = rep(NA_integer_, length(n)),
        diagnostics           = data.frame()
    ))
  }
  # Fit all tried doses; the final choice is restricted to non-eliminated doses.
  estimates[tried] <- stats::qbeta(
    0.5,
    shape1 = 1 + x[tried],
    shape2 = 1 + n[tried] - x[tried]
  )
  isotonic[tried] <- pava_weighted(estimates[tried], weights  = n[tried])
  selected        <- select_closest_target_dose(
    eligible, isotonic, config$target, tolerance  = 1e-12
  )
  list(
    selected              = as.integer(selected),
    posterior_median      = estimates,
    isotonic              = isotonic,
    pava_weight           = ifelse(n > 0L, as.numeric(n), NA_real_),
    borrowed_ess_positive = ifelse(n > 0L, 0, NA_real_),
    matched_adult_dose    = rep(NA_integer_, length(n)),
    diagnostics           = data.frame()
  )
}


# Apply a one-level E/D move or stay within the dose grid. If the candidate has
# been eliminated, use the highest remaining lower dose, or NA if none exists.
next_keyboard_dose <- function(current_dose, action, eliminated) {
  n_dose <- length(eliminated)
  if (action == "E") {
    candidate <- min(n_dose, current_dose + 1L)
  } else if (action == "D") {
    candidate <- max(1L, current_dose - 1L)
  } else {
    candidate <- current_dose
  }
  if (!eliminated[candidate]) return(candidate)
  lower <- which(seq_len(n_dose) < candidate & !eliminated)
  if (length(lower)) return(max(lower))
  NA_integer_
}


# Resolve minimum log-AUC-distance ties by higher adult AUC, then higher nominal
# adult dose. tolerance determines numerical ties in distance and log AUC.
choose_auc_candidate <- function(
  candidate_indices, distances, adult_auc,
  adult_nominal_dose, tolerance  = 1e-12
) {
  minimum_distance <- min(distances[candidate_indices])
  tied             <- candidate_indices[
    abs(distances[candidate_indices] - minimum_distance) <= tolerance
  ]
  if (length(tied) > 1L) {
    highest_auc <- max(adult_auc[tied])
    tied        <- tied[abs(log(adult_auc[tied]) - log(highest_auc)) <= tolerance]
  }
  if (length(tied) > 1L) tied <- tied[which.max(adult_nominal_dose[tied])]
  as.integer(tied[1L])
}


# Match pediatric geometric-mean AUC to adult doses with observations and valid AUC.
# Choose the closest absolute log ratio, then require it to be within the caliper.
# Return the matched index, distance, candidate set, availability and failure reason.
select_auc_match <- function(
  pediatric_gm, adult_gm, adult_n,
  adult_nominal_dose, caliper, tolerance  = 1e-12
) {
  if (!is.finite(pediatric_gm) || pediatric_gm <= 0) {
    return(list(
        match     = NA_integer_,
        distance  = NA_real_,
        available = FALSE,
        eligible  = integer(),
        reason    = "pediatric_auc_unavailable"
    ))
  }
  eligible <- which(adult_n > 0L & is.finite(adult_gm) & adult_gm > 0)
  if (!length(eligible)) {
    return(list(
        match     = NA_integer_,
        distance  = NA_real_,
        available = FALSE,
        eligible  = eligible,
        reason    = "no_eligible_adult_dose"
    ))
  }
  distances           <- rep(Inf, length(adult_gm))
  distances[eligible] <- abs(log(pediatric_gm / adult_gm[eligible]))
  chosen              <- choose_auc_candidate(
    eligible, distances, adult_gm, adult_nominal_dose, tolerance
  )
  if (distances[chosen] > caliper + tolerance) {
    return(list(
        match     = NA_integer_,
        distance  = distances[chosen],
        available = FALSE,
        eligible  = eligible,
        reason    = "outside_caliper"
    ))
  }
  list(
    match     = chosen,
    distance  = distances[chosen],
    available = TRUE,
    eligible  = eligible,
    reason    = "matched"
  )
}


# Apply the matching rule to true scenario AUCs for performance evaluation only.
# This oracle is a diagnostic reference and is not used to allocate patients.
oracle_auc_match <- function(
  pediatric_auc, adult_auc, adult_nominal_dose,
  caliper, tolerance  = 1e-12,
  eligible  = seq_along(adult_auc)
) {
  eligible <- eligible[is.finite(adult_auc[eligible]) & adult_auc[eligible] > 0]
  if (!length(eligible)) return(NA_integer_)
  distances           <- rep(Inf, length(adult_auc))
  distances[eligible] <- abs(log(pediatric_auc / adult_auc[eligible]))
  chosen              <- choose_auc_candidate(
    eligible, distances, adult_auc, adult_nominal_dose, tolerance
  )
  if (distances[chosen] > caliper + tolerance) NA_integer_ else chosen
}


# Simulate one adult Keyboard trial using the shared potential-outcome matrices.
# Save cumulative snapshots after each cohort so pediatric trials can access only
# the adult data available at their corresponding interim analysis.
simulate_adult_trial <- function(scenario, potentials, config) {
  n_dose      <- length(config$adult_doses_mg)
  n           <- integer(n_dose)
  x           <- integer(n_dose)
  log_auc_sum <- numeric(n_dose)
  auc_n       <- integer(n_dose)
  eliminated  <- rep(FALSE, n_dose)
  current     <- 1L
  stopped     <- FALSE
  stop_reason <- "adult_max_n"
  snapshots   <- list()
  path_rows   <- list()
  max_cohorts <- config$adult_max_n %/% config$cohort_size

  for (cohort in seq_len(max_cohorts)) {
    dose    <- current
    indices <- seq.int(n[dose] + 1L, n[dose] + config$cohort_size)
    if (max(indices) > nrow(potentials$adult_dlt_u)) {
      stop("Adult potential outcome matrix was exhausted.", call.  = FALSE)
    }
    cohort_dlt     <- potentials$adult_dlt_u[indices, dose] < scenario$adult_dlt[dose]
    cohort_log_auc <- log(scenario$adult_auc[dose]) +
    config$auc_sigma_log * potentials$adult_auc_z[indices, dose]
    n[dose]           <- n[dose] + config$cohort_size
    x[dose]           <- x[dose] + sum(cohort_dlt)
    log_auc_sum[dose] <- log_auc_sum[dose] + sum(cohort_log_auc)
    auc_n[dose]       <- auc_n[dose] + config$cohort_size

    safety            <- update_keyboard_elimination(n, x, eliminated, config)
    eliminated        <- safety$eliminated
    ordinary_decision <- keyboard_decision(x[dose], n[dose], config)
    if (safety$stop) {
      action      <- "STOP_SAFETY"
      next_dose   <- NA_integer_
      stopped     <- TRUE
      stop_reason <- "adult_lowest_dose_overdose"
    } else if (eliminated[dose]) {
      action    <- "D_SAFETY"
      next_dose <- next_keyboard_dose(dose, "D", eliminated)
    } else {
      action    <- ordinary_decision$action
      next_dose <- next_keyboard_dose(dose, action, eliminated)
    }

    # Retain the complete cumulative state for later pediatric interim lookups.
    snapshot <- list(
      cohort       = cohort,
      total_n      = sum(n),
      n            = n,
      x            = x,
      log_auc_sum  = log_auc_sum,
      auc_n        = auc_n,
      gm_auc       = geometric_mean_from_logs(log_auc_sum, auc_n),
      eliminated   = eliminated,
      current_dose = dose,
      next_dose    = next_dose,
      stopped      = stopped,
      stop_reason  = stop_reason
    )
    snapshots[[cohort]] <- snapshot
    path_rows[[cohort]] <- data.frame(
      adult_cohort           = cohort,
      adult_dose             = dose,
      adult_total_n          = sum(n),
      cohort_dlt             = sum(cohort_dlt),
      cumulative_n_at_dose   = n[dose],
      cumulative_dlt_at_dose = x[dose],
      action                 = action,
      next_dose              = next_dose,
      safety_trigger_dose    = safety$trigger_dose,
      stopped                = stopped,
      stringsAsFactors       = FALSE
    )
    if (stopped || sum(n) >= config$adult_max_n) break
    if (is.na(next_dose)) {
      stopped     <- TRUE
      stop_reason <- "adult_no_assignable_dose"
      break
    }
    current <- next_dose
  }

  list(
    n           = n,
    x           = x,
    log_auc_sum = log_auc_sum,
    auc_n       = auc_n,
    eliminated  = eliminated,
    stopped     = stopped,
    stop_reason = stop_reason,
    snapshots   = snapshots,
    path        = rbind_fill(path_rows)
  )
}


# Retrieve the concurrent adult snapshot after the adult lead-in and current
# pediatric cohort. If the adult trial has ended, use its last completed snapshot.
adult_snapshot_for_pediatric_interim <- function(
  adult_trial, pediatric_cohort,
  config
) {
  if (!length(adult_trial$snapshots)) {
    stop("Adult trial has no completed cohort.", call.  = FALSE)
  }
  lead_cohorts                    <- config$adult_lead_n %/% config$cohort_size
  requested_adult_cohort          <- lead_cohorts + pediatric_cohort
  available_index                 <- min(requested_adult_cohort, length(adult_trial$snapshots))
  snapshot                        <- adult_trial$snapshots[[available_index]]
  snapshot$requested_adult_cohort <- requested_adult_cohort
  snapshot$available_adult_cohort <- available_index
  snapshot$adult_lead_completed   <- snapshot$total_n >= config$adult_lead_n
  snapshot
}

# Preserve the concurrent-data implementation for the availability wrapper
# in 05_simulation.R, which also supports completed-adult sensitivity scenarios.
BASELINE_ADULT_SNAPSHOT_FOR_PEDIATRIC_INTERIM <- adult_snapshot_for_pediatric_interim


# =============================================================================
# Pediatric reference trial using the original Keyboard design
# =============================================================================

# Run the pediatric-only Keyboard comparator with the common trial-output schema.
# Adult snapshots are recorded for comparison, but adult data do not enter its
# posterior, safety rule or final MTD selection.
simulate_keyboard_pediatric <- function(
  method, scenario, replication,
  potentials, adult_trial, config
) {
  if (!identical(method, "Keyboard")) {
    stop("This standalone Keyboard engine accepts method = 'Keyboard' only.",
      call.  = FALSE)
  }
  n_dose               <- 5L
  n                    <- integer(n_dose)
  x                    <- integer(n_dose)
  log_auc_sum          <- numeric(n_dose)
  auc_n                <- integer(n_dose)
  eliminated           <- rep(FALSE, n_dose)
  current              <- 1L
  early_stop           <- FALSE
  stop_reason          <- "pediatric_max_n"
  interim_rows         <- list()
  path_rows            <- list()
  matching_rows        <- list()
  borrowing_rows       <- list()
  max_cohorts          <- config$pediatric_max_n %/% config$cohort_size
  final_adult_snapshot <- NULL

  for (cohort in seq_len(max_cohorts)) {
    dose     <- current
    enrolled <- enroll_pediatric_cohort(
      dose, n, x, log_auc_sum, auc_n, scenario, potentials, config
    )
    n              <- enrolled$n
    x              <- enrolled$x
    log_auc_sum    <- enrolled$log_auc_sum
    auc_n          <- enrolled$auc_n
    child_gm       <- geometric_mean_from_logs(log_auc_sum, auc_n)
    adult_snapshot <- adult_snapshot_for_pediatric_interim(
      adult_trial, cohort, config
    )
    final_adult_snapshot <- adult_snapshot

    # The standard Keyboard design does not borrow adult information.
    borrowing <- list(s_a  = 0, w_post  = 0, borrowed_ess  = 0)
    posterior <- child_only_posterior(x[dose], n[dose])

    ordinary_decision <- keyboard_decision(
      x[dose], n[dose], config, posterior  = posterior
    )
    safety     <- update_keyboard_elimination(n, x, eliminated, config)
    eliminated <- safety$eliminated

    reached_maximum <- sum(n) >= config$pediatric_max_n
    if (safety$stop) {
      action      <- "STOP_SAFETY"
      next_dose   <- NA_integer_
      early_stop  <- TRUE
      stop_reason <- "pediatric_lowest_dose_overdose"
    } else if (reached_maximum) {
      action      <- "END_MAX_N"
      next_dose   <- NA_integer_
      stop_reason <- "pediatric_max_n"
    } else if (eliminated[dose]) {
      action    <- "D_SAFETY"
      next_dose <- next_keyboard_dose(dose, "D", eliminated)
    } else {
      action    <- ordinary_decision$action
      next_dose <- next_keyboard_dose(dose, action, eliminated)
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
      current_gm_auc         = child_gm[dose],
      adult_requested_cohort = adult_snapshot$requested_adult_cohort,
      adult_available_cohort = adult_snapshot$available_adult_cohort,
      adult_total_n          = adult_snapshot$total_n,
      ordinary_action        = ordinary_decision$action,
      actual_action          = action,
      next_dose              = next_dose,
      safety_trigger_dose    = safety$trigger_dose,
      early_stop             = early_stop,
      matched_adult_dose     = NA_integer_,
      match_distance         = NA_real_,
      s_a                    = borrowing$s_a,
      w_post                 = borrowing$w_post,
      borrowed_ess_elir      = borrowing$borrowed_ess,
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
      actual_action    = action,
      next_dose        = next_dose,
      stringsAsFactors = FALSE
    )

    if (early_stop || reached_maximum) break
    if (is.na(next_dose)) {
      early_stop  <- TRUE
      stop_reason <- "pediatric_no_assignable_dose"
      break
    }
    current <- next_dose
  }

  final_selection <- if (early_stop) {
    list(
      selected              = NA_integer_,
      posterior_median      = rep(NA_real_, n_dose),
      isotonic              = rep(NA_real_, n_dose),
      pava_weight           = rep(NA_real_, n_dose),
      borrowed_ess_positive = rep(NA_real_, n_dose),
      matched_adult_dose    = rep(NA_integer_, n_dose),
      diagnostics           = data.frame()
    )
  } else {
    select_keyboard_mtd(n, x, eliminated, config)
  }
  common <- make_trial_summary_and_dose_tables(
    method, scenario, replication, n, x, log_auc_sum, auc_n,
    final_selection$selected, stop_reason, early_stop, final_adult_snapshot,
    eliminated, config
  )
  common$dose$posterior_median            <- final_selection$posterior_median
  common$dose$isotonic_estimate           <- final_selection$isotonic
  common$dose$final_pava_weight           <- final_selection$pava_weight
  common$dose$final_borrowed_ess_positive <- final_selection$borrowed_ess_positive
  common$dose$final_matched_adult_dose    <- final_selection$matched_adult_dose
  c(common, list(
      interim         = rbind_fill(interim_rows),
      path            = rbind_fill(path_rows),
      matching        = rbind_fill(matching_rows),
      borrowing       = rbind_fill(borrowing_rows),
      final_borrowing = data.frame(),
      pa_crm          = data.frame(),
      failure         = data.frame()
  ))
}

