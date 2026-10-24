# =============================================================================
# Simulation engine: trial data, scenario execution, and sensitivity analyses
# =============================================================================
# This module generates reproducible patient-level potential outcomes shared
# across designs and executes the formal simulation study. It also controls
# adult-data availability, sensitivity variants, parallel execution, and the
# scenario-level results saved for reproducibility.
#
# Function guide:
#   seed_from_keys(), set_reproducible_seed(),
#     make_potential_outcomes()           Generate independent reproducible streams.
#   enroll_pediatric_cohort(),
#     make_trial_summary_and_dose_tables() Update trial observations and outputs.
#   empty_trial_output(),
#     combine_trial_outputs()             Combine trial data across replications.
#   adult_snapshot_for_pediatric_interim()
#                                         Handle concurrent/completed adult trials.
#   make_formal_config(), validate_formal_config(),
#     formal_to_core_config()             Configure the primary study.
#   make_formal_sensitivity_variants(),
#     validate_formal_sensitivity_config(),
#     sensitivity_to_core_config(),
#     sensitivity_methods_for_variant(),
#     formal_sensitivity_plan_table()    Configure/tabulate sensitivity analyses.
#   simulate_formal_replication_safe()   Handle replication failures.
#   run_formal_scenario(), run_formal_scenario_checkpointed()
#                                         Run and checkpoint each scenario.
#   assemble_formal_main_results(), run_formal_main()
#                                         Summarize and export main results.
#   simulate_formal_sensitivity_replication(),
#     simulate_formal_sensitivity_replication_safe(),
#     run_formal_sensitivity_variant_scenario(),
#     run_formal_sensitivity()            Execute sensitivity analyses.
#   run_formal_full_analysis()            Combine main and sensitivity results.
#   formal_replication_signature()      Identify the selected replication workflow.
#   formal_config_signature(),
#     formal_sensitivity_signature()      Validate checkpoint compatibility.
#   run_study_method()                  Isolate method failures and attach adult paths.
#   run_sensitivity_scenario_checkpointed()
#                                        Run or restore one sensitivity scenario.
#   assemble_formal_sensitivity_results()
#                                        Summarize sensitivity variants.
#   Other helpers                       Summarize failures, metadata and outputs.
# =============================================================================

# =============================================================================
# Reproducible potential outcomes and shared trial functions
# =============================================================================
seed_from_keys <- function(master_seed, scenario_index, replication, domain) {
  modulus <- 2147483646
  if (scenario_index < 1L || replication < 1L || replication >= 1000000L ||
    domain < 1L || domain > 9L) {
    stop("Invalid keyed RNG indices.", call.  = FALSE)
  }
  # Nonoverlapping decimal fields give a collision-free map for this study:
  # scenario | replication (up to 999,999) | domain (1...9).
  keyed_offset <- (as.numeric(scenario_index) - 1) * 10000000 +
  (as.numeric(replication) - 1) * 10 + as.numeric(domain)
  seed <- (as.numeric(master_seed) + keyed_offset) %% modulus
  as.integer(seed + 1)
}

set_reproducible_seed <- function(seed) {
  set.seed(
    seed,
    kind        = "L'Ecuyer-CMRG",
    normal.kind = "Inversion",
    sample.kind = "Rejection"
  )
  invisible(seed)
}

make_potential_outcomes <- function(scenario_index, replication, config) {
  adult_rows     <- config$adult_max_n
  pediatric_rows <- config$pediatric_max_n
  n_dose         <- 5L

  set_reproducible_seed(seed_from_keys(
      config$master_seed, scenario_index, replication, 1L
  ))
  adult_dlt_u <- matrix(stats::runif(adult_rows * n_dose),
    nrow  = adult_rows, ncol  = n_dose)
  set_reproducible_seed(seed_from_keys(
      config$master_seed, scenario_index, replication, 2L
  ))
  adult_auc_z <- matrix(stats::rnorm(adult_rows * n_dose),
    nrow  = adult_rows, ncol  = n_dose)
  set_reproducible_seed(seed_from_keys(
      config$master_seed, scenario_index, replication, 3L
  ))
  pediatric_dlt_u <- matrix(stats::runif(pediatric_rows * n_dose),
    nrow  = pediatric_rows, ncol  = n_dose)
  set_reproducible_seed(seed_from_keys(
      config$master_seed, scenario_index, replication, 4L
  ))
  pediatric_auc_z <- matrix(stats::rnorm(pediatric_rows * n_dose),
    nrow  = pediatric_rows, ncol  = n_dose)
  list(
    adult_dlt_u     = adult_dlt_u,
    adult_auc_z     = adult_auc_z,
    pediatric_dlt_u = pediatric_dlt_u,
    pediatric_auc_z = pediatric_auc_z
  )
}


enroll_pediatric_cohort <- function(
  dose, n, x, log_auc_sum, auc_n,
  scenario, potentials, config
) {
  indices <- seq.int(n[dose] + 1L, n[dose] + config$cohort_size)
  if (max(indices) > nrow(potentials$pediatric_dlt_u)) {
    stop("Pediatric potential outcome matrix was exhausted.", call.  = FALSE)
  }
  cohort_dlt <- potentials$pediatric_dlt_u[indices, dose] <
  scenario$pediatric_dlt[dose]
  cohort_log_auc <- log(scenario$pediatric_auc[dose]) +
  config$auc_sigma_log * potentials$pediatric_auc_z[indices, dose]
  n[dose]           <- n[dose] + config$cohort_size
  x[dose]           <- x[dose] + sum(cohort_dlt)
  log_auc_sum[dose] <- log_auc_sum[dose] + sum(cohort_log_auc)
  auc_n[dose]       <- auc_n[dose] + config$cohort_size
  list(
    n              = n,
    x              = x,
    log_auc_sum    = log_auc_sum,
    auc_n          = auc_n,
    cohort_dlt     = sum(cohort_dlt),
    cohort_log_auc = cohort_log_auc
  )
}

make_trial_summary_and_dose_tables <- function(
  method, scenario, replication, n, x, log_auc_sum, auc_n,
  selected_mtd, stop_reason, early_stop, final_adult_snapshot,
  eliminated  = rep(FALSE, length(n)), config
) {
  total_n           <- sum(n)
  overtoxic         <- scenario$pediatric_dlt >= 0.35
  true_mtd          <- scenario$true_mtd
  correct_selection <- if (is.na(true_mtd)) {
    NA
  } else {
    !is.na(selected_mtd) && selected_mtd == true_mtd
  }
  correct_no_mtd <- if (is.na(true_mtd)) is.na(selected_mtd) else NA
  summary        <- data.frame(
    config_id                  = config$config_id,
    scenario                   = scenario$id,
    method                     = method,
    replication                = replication,
    selected_mtd               = selected_mtd,
    true_mtd                   = true_mtd,
    correct_selection          = correct_selection,
    correct_no_mtd             = correct_no_mtd,
    selected_no_mtd            = is.na(selected_mtd),
    early_stop                 = early_stop,
    stop_reason                = stop_reason,
    pediatric_n                = total_n,
    total_dlt                  = sum(x),
    overtoxic_n                = sum(n[overtoxic]),
    overtoxic_allocation       = if (total_n > 0) sum(n[overtoxic]) / total_n else NA_real_,
    true_mtd_n                 = if (!is.na(true_mtd)) n[true_mtd] else NA_integer_,
    true_mtd_allocation        = if (!is.na(true_mtd) && total_n > 0) n[true_mtd] / total_n else NA_real_,
    adult_n_at_pediatric_end   = final_adult_snapshot$total_n,
    adult_dlt_at_pediatric_end = sum(final_adult_snapshot$x),
    stringsAsFactors           = FALSE
  )
  dose <- data.frame(
    config_id        = config$config_id,
    scenario         = scenario$id,
    method           = method,
    replication      = replication,
    pediatric_dose   = seq_along(n),
    dose_mg_m2       = config$pediatric_doses_mg_m2,
    true_dlt         = scenario$pediatric_dlt,
    true_auc         = scenario$pediatric_auc,
    n                = n,
    dlt              = x,
    allocation       = if (total_n > 0) n / total_n else NA_real_,
    gm_auc           = geometric_mean_from_logs(log_auc_sum, auc_n),
    eliminated       = eliminated,
    selected         = if (is.na(selected_mtd)) rep(FALSE, length(n)) else seq_along(n) == selected_mtd,
    true_mtd         = if (is.na(true_mtd)) rep(FALSE, length(n)) else seq_along(n) == true_mtd,
    overtoxic        = overtoxic,
    stringsAsFactors = FALSE
  )
  adult <- data.frame(
    config_id        = config$config_id,
    scenario         = scenario$id,
    method           = method,
    replication      = replication,
    adult_dose       = seq_along(final_adult_snapshot$n),
    dose_mg          = config$adult_doses_mg,
    true_dlt         = scenario$adult_dlt,
    true_auc         = scenario$adult_auc,
    n                = final_adult_snapshot$n,
    dlt              = final_adult_snapshot$x,
    gm_auc           = final_adult_snapshot$gm_auc,
    stringsAsFactors = FALSE
  )
  list(trial  = summary, dose  = dose, adult  = adult)
}


empty_trial_output <- function() {
  list(
    trial           = data.frame(),
    dose            = data.frame(),
    adult           = data.frame(),
    interim         = data.frame(),
    path            = data.frame(),
    matching        = data.frame(),
    borrowing       = data.frame(),
    final_borrowing = data.frame(),
    pa_crm          = data.frame(),
    adult_path      = data.frame(),
    failure         = data.frame()
  )
}

combine_trial_outputs <- function(outputs) {
  fields <- names(empty_trial_output())
  setNames(lapply(fields, function(field) {
        rbind_fill(lapply(outputs, `[[`, field))
    }), fields)
}


# =============================================================================
# Adult information available at pediatric interim analyses
# =============================================================================
adult_snapshot_for_pediatric_interim <- function(
  adult_trial,
  pediatric_cohort,
  config
) {

  completed_mode <- isTRUE(config$formal_adult_completed)

  if (!completed_mode) {
    return(
      BASELINE_ADULT_SNAPSHOT_FOR_PEDIATRIC_INTERIM(
        adult_trial,
        pediatric_cohort,
        config
      )
    )
  }

  if (!length(adult_trial$snapshots)) {
    stop("Adult trial has no completed cohort.", call.  = FALSE)
  }

  available_index                 <- length(adult_trial$snapshots)
  snapshot                        <- adult_trial$snapshots[[available_index]]
  snapshot$requested_adult_cohort <- available_index
  snapshot$available_adult_cohort <- available_index
  snapshot$adult_lead_completed   <- TRUE
  snapshot
}


# =============================================================================
# Main and sensitivity simulation workflows
# =============================================================================
make_formal_config <- function(
  overrides  = list(),
  output_dir = NULL
) {

  if (is.null(output_dir)) {
    output_dir <- file.path(
      getwd(),
      "PA_keyboard_ETOVL_formal_v03_results"
    )
  }

  if (!"output_dir" %in% names(overrides)) {
    overrides$output_dir <- output_dir
  }

  # Construct the exposure-guided borrowing configuration.
  base <- make_pa_keyboard_base_config(
    baseline_overrides = overrides,
    output_dir         = output_dir
  )

  # Apply the prespecified adult ESS cap.
  cfg <- make_pa_keyboard_variants(base)[["PAK_capped"]]

  # Assign the formal analysis label.
  cfg$config_id                <- "MAIN"
  cfg$pa_keyboard_diag_variant <- "Formal_ETOVL"
  cfg$pa_keyboard_diag_label   <- paste0(
    "Formal PA-Keyboard ET-OVL: nP>=6, adult ESS cap=6"
  )

  # Primary study settings: pediatric borrowing gate and adult ESS cap.
  cfg$pa_keyboard_min_pediatric_n       <- 6L
  cfg$pa_keyboard_final_min_pediatric_n <- 6L
  cfg$pa_keyboard_ess_cap               <- 6

  # Adult information is concurrent in the primary analysis.
  cfg$formal_adult_completed <- FALSE
  cfg$formal_adult_data_mode <- "Concurrent: adult lead-in = 6, then same cohort speed"

  validate_formal_config(cfg)
}

validate_formal_config <- function(config) {
  # Validate the shared fields before checking the formal design parameters.
  known                          <- config
  known$pa_keyboard_ess_cap      <- NULL
  known$pa_keyboard_diag_variant <- NULL
  known$pa_keyboard_diag_label   <- NULL
  known                          <- validate_pa_keyboard_config(known)

  assert_true(
    identical(as.integer(config$pa_keyboard_min_pediatric_n), 6L),
    "Formal ET-OVL borrowing must start at pediatric n >= 6."
  )
  assert_true(
    identical(as.integer(config$pa_keyboard_final_min_pediatric_n), 6L),
    "Formal ET-OVL final borrowing must require pediatric n >= 6."
  )
  assert_true(
    isTRUE(all.equal(as.numeric(config$pa_keyboard_ess_cap), 6)),
    "Formal ET-OVL adult effective-sample-size cap must equal 6."
  )
  assert_true(
    identical(as.integer(config$safety_min_n), 3L),
    "Keyboard pediatric safety elimination must remain available from n >= 3."
  )
  assert_true(
    is.logical(config$formal_adult_completed) &&
    length(config$formal_adult_completed) == 1L &&
    !is.na(config$formal_adult_completed),
    "formal_adult_completed must be TRUE/FALSE."
  )
  config
}

formal_to_core_config <- function(config) {
  config                        <- validate_formal_config(config)
  core                          <- config
  core$pa_keyboard_diag_variant <- "PAK_capped"
  core$pa_keyboard_diag_label   <- PA_KEYBOARD_VARIANT_LABELS[["PAK_capped"]]
  core
}

validate_formal_sensitivity_config <- function(config) {
  # Validate the shared fields before checking sensitivity parameters.
  known                          <- config
  known$pa_keyboard_ess_cap      <- NULL
  known$pa_keyboard_diag_variant <- NULL
  known$pa_keyboard_diag_label   <- NULL
  known                          <- validate_pa_keyboard_config(known)

  assert_true(
    config$config_id %in% FORMAL_SENSITIVITY_PLAN$variant_id,
    "Unknown formal sensitivity variant."
  )
  assert_true(
    is.numeric(config$pa_keyboard_ess_cap) && length(config$pa_keyboard_ess_cap) == 1L &&
    is.finite(config$pa_keyboard_ess_cap) && config$pa_keyboard_ess_cap > 0,
    "Sensitivity pa_keyboard_ess_cap must be finite and > 0."
  )
  assert_true(
    config$pa_keyboard_min_pediatric_n %in% c(3L, 6L) &&
    config$pa_keyboard_final_min_pediatric_n %in% c(3L, 6L),
    "Sensitivity borrowing gate must be pediatric n=3 or n=6."
  )
  assert_true(
    is.logical(config$formal_adult_completed) &&
    length(config$formal_adult_completed) == 1L &&
    !is.na(config$formal_adult_completed),
    "formal_adult_completed must be TRUE/FALSE."
  )
  config
}

sensitivity_to_core_config <- function(config) {
  config                        <- validate_formal_sensitivity_config(config)
  core                          <- config
  core$pa_keyboard_diag_variant <- "PAK_capped"
  core$pa_keyboard_diag_label   <- PA_KEYBOARD_VARIANT_LABELS[["PAK_capped"]]
  validate_pa_keyboard_variant_config(core)
}

make_formal_sensitivity_variants <- function(
  base_config  = make_formal_config()
) {

  base_config <- validate_formal_config(base_config)

  make_one <- function(id, group, setting) {
    cfg                            <- base_config
    cfg$config_id                  <- id
    cfg$formal_sensitivity_group   <- group
    cfg$formal_sensitivity_setting <- setting
    cfg
  }

  out <- list()

  out$CAL16                 <- make_one("CAL16", "AUC caliper", "0.16")
  out$CAL16$auc_caliper_log <- 0.16

  out$CAL28                 <- make_one("CAL28", "AUC caliper", "0.28")
  out$CAL28$auc_caliper_log <- 0.28

  out$LEAD9              <- make_one("LEAD9", "Adult information", "Lead 9")
  out$LEAD9$adult_lead_n <- 9L
  out$LEAD9$formal_adult_data_mode <-
  "Concurrent: adult lead-in = 9, then same cohort speed"

  out$LEAD12              <- make_one("LEAD12", "Adult information", "Lead 12")
  out$LEAD12$adult_lead_n <- 12L
  out$LEAD12$formal_adult_data_mode <-
  "Concurrent: adult lead-in = 12, then same cohort speed"

  out$ADULT_COMPLETED <- make_one(
    "ADULT_COMPLETED",
    "Adult information",
    "Adult completed"
  )
  out$ADULT_COMPLETED$formal_adult_completed <- TRUE
  out$ADULT_COMPLETED$formal_adult_data_mode <-
  "Adult trial completed (or early-stopped) before pediatric start"

  out$CAP3                     <- make_one("CAP3", "Adult ESS cap", "3")
  out$CAP3$pa_keyboard_ess_cap <- 3

  out$CAP9                     <- make_one("CAP9", "Adult ESS cap", "9")
  out$CAP9$pa_keyboard_ess_cap <- 9

  out$GATE3 <- make_one("GATE3", "Borrowing gate", "nP >= 3")
  out$GATE3$pa_keyboard_min_pediatric_n <- 3L
  out$GATE3$pa_keyboard_final_min_pediatric_n <- 3L

  # Keep the order exactly as pre-specified in FORMAL_SENSITIVITY_PLAN.
  out <- out[FORMAL_SENSITIVITY_PLAN$variant_id]
  lapply(out, validate_formal_sensitivity_config)
}

sensitivity_methods_for_variant <- function(variant_id) {
  if (variant_id %in% c("LEAD9", "LEAD12", "ADULT_COMPLETED")) {
    # Adult information availability changes both borrowing designs.
    return(c("PA-Keyboard", "PA-CRM"))
  }

  # Caliper / ET-ESS cap / pediatric borrowing gate are PA-Keyboard-specific.
  "PA-Keyboard"
}

formal_sensitivity_plan_table <- function() {
  data.frame(
    Group  = c(
      "AUC caliper", "AUC caliper", "AUC caliper",
      "Adult information", "Adult information", "Adult information", "Adult information",
      "Adult ESS cap", "Adult ESS cap", "Adult ESS cap",
      "Borrowing gate", "Borrowing gate"
    ),
    Setting  = c(
      "0.16", "0.22 (MAIN)", "0.28",
      "Lead 6 (MAIN)", "Lead 9", "Lead 12", "Adult completed",
      "3", "6 (MAIN)", "9",
      "nP >= 3", "nP >= 6 (MAIN)"
    ),
    Scenarios        = paste(FORMAL_SENSITIVITY_SCENARIOS, collapse  = ", "),
    Analysis         = "One-factor-at-a-time",
    stringsAsFactors = FALSE
  )
}

relabel_method_in_output <- function(result, new_method) {
  for (nm in names(result)) {
    obj <- result[[nm]]
    if (is.data.frame(obj) && "method" %in% names(obj)) {
      obj$method   <- new_method
      result[[nm]] <- obj
    }
  }
  result
}

make_formal_failure_output <- function(
  scenario_id,
  method,
  replication,
  error
) {

  z         <- empty_trial_output()
  z$failure <- data.frame(
    config_id        = "MAIN",
    scenario         = scenario_id,
    method           = method,
    replication      = replication,
    error_class      = class(error)[1L],
    error_message    = conditionMessage(error),
    stringsAsFactors = FALSE
  )
  z
}

attach_adult_path <- function(
  result,
  adult_trial,
  scenario_id,
  method,
  replication,
  config
) {

  result$adult_path <- transform(
    adult_trial$path,
    config_id   = config$config_id,
    scenario    = scenario_id,
    method      = method,
    replication = replication
  )
  result
}

simulate_formal_replication_safe <- function(
  scenario,
  replication,
  config,
  replication_function  = simulate_formal_replication
) {

  tryCatch(
    replication_function(scenario, replication, config),
    error  = function(e) {
      if (isTRUE(config$fail_fast)) stop(e)
      scenario_id <- if (is.character(scenario)) scenario else scenario$id
      setNames(
        lapply(
          FORMAL_METHODS,
          function(method) {
            make_formal_failure_output(
              scenario_id, method, replication, e
            )
          }
        ),
        FORMAL_METHODS
      )
    }
  )
}

combine_formal_outputs <- function(replication_outputs) {
  by_method <- lapply(FORMAL_METHODS, function(method) {
      combine_trial_outputs(lapply(replication_outputs, `[[`, method))
  })
  names(by_method) <- FORMAL_METHODS

  list(
    by_method = by_method,
    combined  = combine_trial_outputs(by_method)
  )
}


formal_config_signature <- function(config, n_sim  = config$n_sim) {
  paste(
    FORMAL_FILE_VERSION,
    paste(FORMAL_METHODS, collapse  = ","),
    paste(FORMAL_SCENARIOS, collapse  = ","),
    paste0("target=", config$target),
    paste0("cohort=", config$cohort_size),
    paste0("adultN=", config$adult_max_n),
    paste0("pedN=", config$pediatric_max_n),
    paste0("lead=", config$adult_lead_n),
    paste0("adult_completed=", isTRUE(config$formal_adult_completed)),
    paste0("auc_sigma=", config$auc_sigma_log),
    paste0("auc_caliper=", config$auc_caliper_log),
    paste0("borrow_min=", config$pa_keyboard_min_pediatric_n),
    paste0("final_min=", config$pa_keyboard_final_min_pediatric_n),
    paste0("ess_cap=", config$pa_keyboard_ess_cap),
    paste0("safety_min=", config$safety_min_n),
    paste0("safety_cut=", config$safety_cutoff),
    paste0("nsim=", as.integer(n_sim)),
    paste0("seed=", config$master_seed),
    paste0("all_parameters=", config_signature(config)),
    sep  = "|"
  )
}

detect_formal_script_path <- function() {
  # The entrypoint path identifies the project root for worker module loading.
  path <- get0("PAK_ENTRY_SCRIPT", envir  = .GlobalEnv, inherits  = FALSE)
  if (is.character(path) && length(path) == 1L && file.exists(path)) {
    return(normalizePath(path, mustWork  = TRUE))
  }
  candidate <- file.path(getwd(), "run_analysis.R")
  if (file.exists(candidate)) return(normalizePath(candidate, mustWork  = TRUE))
  stop("Cannot find run_analysis.R; source the project entrypoint first.",
    call.  = FALSE)
}

run_formal_scenario <- function(
  scenario_id,
  config,
  n_sim                = config$n_sim,
  n_workers            = config$n_workers,
  script_path          = detect_formal_script_path(),
  replication_function = simulate_formal_replication
) {

  config    <- validate_formal_config(config)
  scenarios <- make_scenarios()
  if (!scenario_id %in% FORMAL_SCENARIOS) {
    stopf("Scenario %s is not in FORMAL_SCENARIOS.", scenario_id)
  }
  scenario <- scenarios[[scenario_id]]

  n_sim     <- as.integer(n_sim)
  n_workers <- as.integer(max(1L, min(n_workers, n_sim)))

  message("============================================================")
  message("Formal ET-OVL scenario: ", scenario_id)
  message("Methods: ", paste(FORMAL_METHODS, collapse  = ", "))
  message("N simulations: ", n_sim)
  message("Workers: ", n_workers)
  message("Borrowing gate: pediatric n >= 6")
  message("Adult borrowed ESS cap: 6")
  message("============================================================")

  t0 <- proc.time()[["elapsed"]]

  if (n_workers == 1L) {
    rep_outputs <- lapply(seq_len(n_sim), function(r) {
        simulate_formal_replication_safe(
          scenario             = scenario,
          replication          = r,
          config               = config,
          replication_function = replication_function
        )
    })
  } else {
    cl <- parallel::makePSOCKcluster(n_workers)
    on.exit(parallel::stopCluster(cl), add  = TRUE)

    pak_initialize_workers(cl, dirname(script_path))

    rep_outputs <- parallel::parLapplyLB(
      cl,
      seq_len(n_sim),
      function(r, scenario, config, replication_function) {
        simulate_formal_replication_safe(
          scenario             = scenario,
          replication          = r,
          config               = config,
          replication_function = replication_function
        )
      },
      scenario             = scenario,
      config               = config,
      replication_function = replication_function
    )
  }

  runtime_sec <- proc.time()[["elapsed"]] - t0
  combined    <- combine_formal_outputs(rep_outputs)
  raw         <- combined$combined

  result <- list(
    scenario              = scenario_id,
    n_sim                 = n_sim,
    n_workers             = n_workers,
    config_signature      = formal_config_signature(config, n_sim),
    replication_signature = formal_replication_signature(replication_function),
    runtime_sec           = runtime_sec,
    main                  = make_formal_main_table(raw, runtime_sec),
    dose                  = make_formal_dose_table(raw),
    matching              = make_formal_matching_table(raw),
    final_borrowing       = make_formal_final_borrowing_table(raw),
    online_borrowing      = make_formal_online_borrowing_table(raw),
    failures              = raw$failure
  )

  message(
    "Finished ", scenario_id,
    " in ", round(runtime_sec / 60, 3), " min."
  )
  result
}

run_formal_scenario_checkpointed <- function(
  scenario_id, config, n_sim, n_workers, resume  = TRUE,
  script_path          = detect_formal_script_path(),
  replication_function = simulate_formal_replication
) {
  dir.create(config$output_dir, recursive  = TRUE, showWarnings  = FALSE)
  checkpoint_dir <- file.path(config$output_dir, "checkpoints")
  dir.create(checkpoint_dir, recursive  = TRUE, showWarnings  = FALSE)
  checkpoint <- file.path(checkpoint_dir, paste0(scenario_id, "_formal_summary_v03.rds"))
  if (isTRUE(resume) && file.exists(checkpoint)) {
    cached       <- readRDS(checkpoint)
    expected_sig <- formal_config_signature(config, n_sim)
    if (!is.null(cached$config_signature) &&
      identical(cached$config_signature, expected_sig) &&
      identical(cached$replication_signature,
        formal_replication_signature(replication_function))) {
      message("Loading completed checkpoint: ", scenario_id)
      return(cached)
    }
    message("Ignoring stale checkpoint (config/run mismatch): ", scenario_id)
  }
  result <- run_formal_scenario(
    scenario_id          = scenario_id,
    config               = config,
    n_sim                = n_sim,
    n_workers            = n_workers,
    script_path          = script_path,
    replication_function = replication_function
  )
  saveRDS(result, checkpoint, compress  = "gzip")
  gc(verbose  = FALSE)
  result
}

assemble_formal_main_results <- function(
  pieces, config, scenario_ids, n_sim, n_workers,
  output_xlsx      = file.path(config$output_dir, "PA_keyboard_ETOVL_formal_results.xlsx"),
  save_compact_rds = TRUE
) {
  main     <- rbind_fill(lapply(pieces, `[[`, "main"))
  dose     <- rbind_fill(lapply(pieces, `[[`, "dose"))
  matching <- rbind_fill(lapply(pieces, `[[`, "matching"))
  final    <- rbind_fill(lapply(pieces, `[[`, "final_borrowing"))
  online   <- rbind_fill(lapply(pieces, `[[`, "online_borrowing"))

  main <- sort_formal_rows(main)

  if (nrow(dose)) {
    dose_num <- as.integer(sub("^P", "", dose$Dose))
    dose     <- dose[
      order(
        match(dose$Scenario, FORMAL_SCENARIOS),
        match(dose$Method, FORMAL_METHODS),
        dose_num
      ),
      ,
      drop  = FALSE
    ]
  }

  if (nrow(matching)) {
    matching_dose <- as.integer(sub("^P", "", matching$Dose))
    matching      <- matching[
      order(match(matching$Scenario, FORMAL_SCENARIOS), matching_dose),
      ,
      drop  = FALSE
    ]
  }

  if (nrow(final)) {
    final_dose <- as.integer(sub("^P", "", final$Dose))
    final      <- final[
      order(match(final$Scenario, FORMAL_SCENARIOS), final_dose),
      ,
      drop  = FALSE
    ]
  }

  if (nrow(online)) {
    online <- online[
      order(match(online$Scenario, FORMAL_SCENARIOS)),
      ,
      drop  = FALSE
    ]
  }

  write_formal_workbook(
    config           = config,
    main_table       = main,
    dose_table       = dose,
    matching_table   = matching,
    final_borrowing  = final,
    online_borrowing = online,
    path             = output_xlsx
  )

  compact <- list(
    metadata  = list(
      code_version   = FORMAL_FILE_VERSION,
      formal_version = FORMAL_VERSION,
      protocol_date  = FORMAL_PROTOCOL_DATE,
      methods        = FORMAL_METHODS,
      scenarios      = scenario_ids,
      n_sim          = n_sim,
      n_workers      = n_workers,
      completed_at   = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      output_xlsx    = output_xlsx
    ),
    config           = config,
    main             = main,
    dose             = dose,
    matching         = matching,
    final_borrowing  = final,
    online_borrowing = online
  )

  if (isTRUE(save_compact_rds)) {
    saveRDS(
      compact,
      file.path(
        config$output_dir,
        "PA_keyboard_ETOVL_formal_main_compact_v03.rds"
      ),
      compress  = "gzip"
    )
  }

  message("============================================================")
  message("FORMAL MAIN ANALYSIS COMPLETE")
  message("Single XLSX: ", normalizePath(output_xlsx, mustWork  = FALSE))
  message(
    "Sheets: 00_Design / 01_Main / 02_Dose / 03_Borrowing / 04_Sensitivity"
  )
  message("04_Sensitivity contains a placeholder until sensitivity is run.")
  message("============================================================")

  invisible(compact)
}

run_formal_main <- function(
  config           = make_formal_config(),
  scenario_ids     = FORMAL_SCENARIOS,
  n_sim            = config$n_sim,
  n_workers        = config$n_workers,
  output_xlsx      = file.path(config$output_dir, "PA_keyboard_ETOVL_formal_results.xlsx"),
  save_compact_rds = TRUE,
  resume           = TRUE,
  script_path      = detect_formal_script_path()
) {
  config       <- validate_formal_config(config)
  scenario_ids <- intersect(FORMAL_SCENARIOS, scenario_ids)
  if (!length(scenario_ids)) stop("No valid formal scenarios requested.", call.  = FALSE)
  pieces <- setNames(vector("list", length(scenario_ids)), scenario_ids)
  for (id in scenario_ids) {
    pieces[[id]] <- run_formal_scenario_checkpointed(
      scenario_id = id,
      config      = config,
      n_sim       = n_sim,
      n_workers   = n_workers,
      resume      = resume,
      script_path = script_path
    )
  }
  assemble_formal_main_results(
    pieces           = pieces,
    config           = config,
    scenario_ids     = scenario_ids,
    n_sim            = n_sim,
    n_workers        = n_workers,
    output_xlsx      = output_xlsx,
    save_compact_rds = save_compact_rds
  )
}

combine_outputs_for_methods <- function(replication_outputs, methods) {
  by_method <- lapply(methods, function(method) {
      combine_trial_outputs(lapply(replication_outputs, `[[`, method))
  })
  names(by_method) <- methods

  list(
    by_method = by_method,
    combined  = combine_trial_outputs(by_method)
  )
}

make_sensitivity_failure_output <- function(
  config,
  scenario_id,
  method,
  replication,
  error
) {

  z         <- empty_trial_output()
  z$failure <- data.frame(
    config_id        = config$config_id,
    scenario         = scenario_id,
    method           = method,
    replication      = replication,
    error_class      = class(error)[1L],
    error_message    = conditionMessage(error),
    stringsAsFactors = FALSE
  )
  z
}

simulate_formal_sensitivity_replication <- function(
  scenario,
  replication,
  config,
  methods  = sensitivity_methods_for_variant(config$config_id)
) {

  config    <- validate_formal_sensitivity_config(config)
  scenarios <- make_scenarios()

  if (is.character(scenario)) scenario <- scenarios[[scenario]]
  if (is.null(scenario) || !(scenario$id %in% FORMAL_SENSITIVITY_SCENARIOS)) {
    stop("Unknown/non-sensitivity formal scenario.", call.  = FALSE)
  }

  unknown <- setdiff(methods, c("PA-Keyboard", "PA-CRM"))
  if (length(unknown)) {
    stopf("Unknown sensitivity method(s): %s", paste(unknown, collapse  = ", "))
  }

  scenario_index <- match(scenario$id, names(scenarios))
  potentials     <- make_potential_outcomes(
    scenario_index,
    replication,
    config
  )
  adult_trial <- simulate_adult_trial(
    scenario,
    potentials,
    config
  )

  run_one_method <- function(method) {
    tryCatch(
      {
        result <- if (method == "PA-Keyboard") {
          core_cfg <- sensitivity_to_core_config(config)
          z        <- simulate_pa_keyboard_pediatric(
            scenario    = scenario,
            replication = replication,
            potentials  = potentials,
            adult_trial = adult_trial,
            config      = core_cfg
          )
          relabel_method_in_output(z, "PA-Keyboard")
        } else {
          z <- simulate_pa_crm_pediatric(
            scenario,
            replication,
            potentials,
            adult_trial,
            config
          )
          relabel_method_in_output(z, "PA-CRM")
        }

        attach_adult_path(
          result,
          adult_trial,
          scenario$id,
          method,
          replication,
          config
        )
      },
      error  = function(e) {
        if (isTRUE(config$fail_fast)) stop(e)
        make_sensitivity_failure_output(
          config,
          scenario$id,
          method,
          replication,
          e
        )
      }
    )
  }

  setNames(lapply(methods, run_one_method), methods)
}

simulate_formal_sensitivity_replication_safe <- function(
  scenario,
  replication,
  config,
  methods  = sensitivity_methods_for_variant(config$config_id)
) {

  tryCatch(
    simulate_formal_sensitivity_replication(
      scenario,
      replication,
      config,
      methods
    ),
    error  = function(e) {
      if (isTRUE(config$fail_fast)) stop(e)
      scenario_id <- if (is.character(scenario)) scenario else scenario$id
      setNames(
        lapply(methods, function(method) {
            make_sensitivity_failure_output(
              config,
              scenario_id,
              method,
              replication,
              e
            )
        }),
        methods
      )
    }
  )
}

formal_sensitivity_signature <- function(config, n_sim) {
  paste(
    FORMAL_FILE_VERSION,
    config$config_id,
    paste0("n_sim=", n_sim),
    paste0("seed=", config$master_seed),
    paste0("lead=", config$adult_lead_n),
    paste0("adult_completed=", isTRUE(config$formal_adult_completed)),
    paste0("caliper=", format(config$auc_caliper_log, digits  = 16)),
    paste0("borrow_min=", config$pa_keyboard_min_pediatric_n),
    paste0("final_min=", config$pa_keyboard_final_min_pediatric_n),
    paste0("cap=", format(config$pa_keyboard_ess_cap, digits  = 16)),
    paste0("all_parameters=", config_signature(config)),
    sep  = "|"
  )
}

run_formal_sensitivity_variant_scenario <- function(
  scenario_id,
  config,
  n_sim       = config$n_sim,
  n_workers   = config$n_workers,
  methods     = sensitivity_methods_for_variant(config$config_id),
  script_path = detect_formal_script_path()
) {

  config <- validate_formal_sensitivity_config(config)
  if (!scenario_id %in% FORMAL_SENSITIVITY_SCENARIOS) {
    stopf("Scenario %s is not a formal sensitivity scenario.", scenario_id)
  }

  scenarios <- make_scenarios()
  scenario  <- scenarios[[scenario_id]]
  n_sim     <- as.integer(n_sim)
  n_workers <- as.integer(max(1L, min(n_workers, n_sim)))

  message("============================================================")
  message("Sensitivity: ", config$config_id, " | ", scenario_id)
  message("Group: ", config$formal_sensitivity_group)
  message("Setting: ", config$formal_sensitivity_setting)
  message("Methods: ", paste(methods, collapse  = ", "))
  message("N simulations: ", n_sim)
  message("Workers: ", n_workers)
  message("Adult data mode: ", config$formal_adult_data_mode)
  message("============================================================")

  t0 <- proc.time()[["elapsed"]]

  if (n_workers == 1L) {
    rep_outputs <- lapply(seq_len(n_sim), function(r) {
        simulate_formal_sensitivity_replication_safe(
          scenario    = scenario,
          replication = r,
          config      = config,
          methods     = methods
        )
    })
  } else {
    cl <- parallel::makePSOCKcluster(n_workers)
    on.exit(parallel::stopCluster(cl), add  = TRUE)

    pak_initialize_workers(cl, dirname(script_path))

    rep_outputs <- parallel::parLapplyLB(
      cl,
      seq_len(n_sim),
      function(r, scenario, config, methods) {
        simulate_formal_sensitivity_replication_safe(
          scenario    = scenario,
          replication = r,
          config      = config,
          methods     = methods
        )
      },
      scenario = scenario,
      config   = config,
      methods  = methods
    )
  }

  runtime_sec <- proc.time()[["elapsed"]] - t0
  combined    <- combine_outputs_for_methods(rep_outputs, methods)
  raw         <- combined$combined

  list(
    variant_id       = config$config_id,
    group            = config$formal_sensitivity_group,
    setting          = config$formal_sensitivity_setting,
    scenario         = scenario_id,
    methods          = methods,
    n_sim            = n_sim,
    n_workers        = n_workers,
    runtime_sec      = runtime_sec,
    config_signature = formal_sensitivity_signature(config, n_sim),
    main             = make_formal_main_table(raw, runtime_sec),
    online_borrowing = make_formal_online_borrowing_table(raw),
    failures         = raw$failure
  )
}


run_formal_sensitivity <- function(
  base_config  = make_formal_config(),
  main_results,
  scenario_ids    = FORMAL_SENSITIVITY_SCENARIOS,
  n_sim           = base_config$n_sim,
  n_workers       = base_config$n_workers,
  resume          = TRUE,
  update_workbook = TRUE,
  output_xlsx     = file.path(
    base_config$output_dir,
    "PA_keyboard_ETOVL_formal_results.xlsx"
  ),
  script_path  = detect_formal_script_path()
) {

  base_config  <- validate_formal_config(base_config)
  scenario_ids <- intersect(FORMAL_SENSITIVITY_SCENARIOS, scenario_ids)
  if (!length(scenario_ids)) {
    stop("No valid sensitivity scenarios requested.", call.  = FALSE)
  }
  if (missing(main_results) || is.null(main_results$main)) {
    stop(
      "main_results must be the compact object returned by run_formal_main().",
      call.  = FALSE
    )
  }

  variants       <- make_formal_sensitivity_variants(base_config)
  checkpoint_dir <- file.path(
    base_config$output_dir,
    "sensitivity_checkpoints"
  )
  dir.create(checkpoint_dir, recursive  = TRUE, showWarnings  = FALSE)

  pieces <- list()

  for (variant_id in names(variants)) {
    cfg     <- variants[[variant_id]]
    methods <- sensitivity_methods_for_variant(variant_id)

    for (scenario_id in scenario_ids) {
      key        <- paste(variant_id, scenario_id, sep  = "__")
      checkpoint <- file.path(
        checkpoint_dir,
        paste0(key, "_summary.rds")
      )

      expected_sig <- formal_sensitivity_signature(cfg, n_sim)
      if (isTRUE(resume) && file.exists(checkpoint)) {
        cached <- readRDS(checkpoint)
        if (
          !is.null(cached$config_signature) &&
          identical(cached$config_signature, expected_sig)
        ) {
          message("Loading sensitivity checkpoint: ", key)
          pieces[[key]] <- cached
          next
        }
      }

      ans <- run_formal_sensitivity_variant_scenario(
        scenario_id = scenario_id,
        config      = cfg,
        n_sim       = n_sim,
        n_workers   = n_workers,
        methods     = methods,
        script_path = script_path
      )
      pieces[[key]] <- ans
      saveRDS(ans, checkpoint, compress  = "gzip")
      gc(verbose  = FALSE)
    }
  }

  assemble_formal_sensitivity_results(
    pieces, base_config, main_results, scenario_ids, n_sim, n_workers,
    variants, update_workbook, output_xlsx
  )
}

assemble_formal_sensitivity_results <- function(
  pieces, base_config, main_results, scenario_ids, n_sim, n_workers, variants,
  update_workbook, output_xlsx
) {
  variant_rows <- rbind_fill(lapply(pieces, function(x) {
        add_sensitivity_metadata(
          x$main,
          group            = x$group,
          variant_id       = x$variant_id,
          setting          = x$setting,
          reference        = FALSE,
          online_borrowing = x$online_borrowing
        )
  }))

  reference_rows <- make_sensitivity_reference_rows(main_results)
  table          <- sort_sensitivity_rows(
    rbind_fill(list(reference_rows, variant_rows))
  )

  sensitivity <- list(
    metadata  = list(
      code_version = FORMAL_FILE_VERSION,
      scenarios    = scenario_ids,
      n_sim        = n_sim,
      n_workers    = n_workers,
      variants     = names(variants),
      completed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    ),
    plan   = formal_sensitivity_plan_table(),
    table  = table,
    pieces = pieces
  )

  saveRDS(
    sensitivity,
    file.path(
      base_config$output_dir,
      "PA_keyboard_ETOVL_formal_sensitivity_compact_v03.rds"
    ),
    compress  = "gzip"
  )

  if (isTRUE(update_workbook)) {
    write_formal_workbook(
      config            = base_config,
      main_table        = main_results$main,
      dose_table        = main_results$dose,
      matching_table    = main_results$matching,
      final_borrowing   = main_results$final_borrowing,
      online_borrowing  = main_results$online_borrowing,
      sensitivity_table = table,
      path              = output_xlsx
    )
  }

  invisible(sensitivity)
}

run_formal_full_analysis <- function(
  config                = make_formal_config(),
  main_scenarios        = FORMAL_SCENARIOS,
  sensitivity_scenarios = FORMAL_SENSITIVITY_SCENARIOS,
  n_sim                 = config$n_sim,
  n_workers             = config$n_workers,
  resume                = TRUE,
  output_xlsx           = file.path(
    config$output_dir,
    "PA_keyboard_ETOVL_formal_results.xlsx"
  ),
  script_path  = detect_formal_script_path()
) {

  config <- validate_formal_config(config)

  main <- run_formal_main(
    config           = config,
    scenario_ids     = main_scenarios,
    n_sim            = n_sim,
    n_workers        = n_workers,
    output_xlsx      = output_xlsx,
    save_compact_rds = TRUE,
    resume           = resume,
    script_path      = script_path
  )

  sensitivity <- run_formal_sensitivity(
    base_config     = config,
    main_results    = main,
    scenario_ids    = sensitivity_scenarios,
    n_sim           = n_sim,
    n_workers       = n_workers,
    resume          = resume,
    update_workbook = TRUE,
    output_xlsx     = output_xlsx,
    script_path     = script_path
  )

  out <- list(
    metadata  = list(
      code_version = FORMAL_FILE_VERSION,
      completed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      output_xlsx  = output_xlsx
    ),
    main        = main,
    sensitivity = sensitivity
  )

  saveRDS(
    out,
    file.path(
      config$output_dir,
      "PA_keyboard_ETOVL_formal_full_compact_v03.rds"
    ),
    compress  = "gzip"
  )

  message("============================================================")
  message("FORMAL MAIN + SENSITIVITY ANALYSES COMPLETE")
  message("Single XLSX: ", normalizePath(output_xlsx, mustWork  = FALSE))
  message(
    "Sheets: 00_Design / 01_Main / 02_Dose / 03_Borrowing / 04_Sensitivity"
  )
  message("============================================================")

  invisible(out)
}


run_study_method <- function(
  method, expression, scenario, replication,
  adult_trial, config
) {
  tryCatch({
      result <- relabel_method_in_output(force(expression), method)
      attach_adult_path(result, adult_trial, scenario$id, method, replication, config)
    }, error  = function(e) {
      if (isTRUE(config$fail_fast)) stop(e)
      make_formal_failure_output(scenario$id, method, replication, e)
  })
}

run_sensitivity_scenario_checkpointed <- function(
  scenario_id, config, n_sim,
  n_workers, resume  = TRUE, script_path  = detect_formal_script_path()
) {
  checkpoint_dir <- file.path(config$output_dir, "sensitivity_checkpoints")
  dir.create(checkpoint_dir, recursive  = TRUE, showWarnings  = FALSE)
  key          <- paste(config$config_id, scenario_id, sep  = "__")
  checkpoint   <- file.path(checkpoint_dir, paste0(key, "_summary.rds"))
  expected_sig <- formal_sensitivity_signature(config, n_sim)
  if (isTRUE(resume) && file.exists(checkpoint)) {
    cached <- readRDS(checkpoint)
    if (!is.null(cached$config_signature) &&
      identical(cached$config_signature, expected_sig)) return(cached)
  }
  ans <- run_formal_sensitivity_variant_scenario(
    scenario_id, config, n_sim, n_workers,
    sensitivity_methods_for_variant(config$config_id), script_path)
  saveRDS(ans, checkpoint, compress  = "gzip")
  ans
}

formal_replication_signature <- function(replication_function) {
  # Store the callback definition as well as parameters in checkpoint metadata.
  # If external objects used by a custom callback change, rerun with resume = FALSE.
  stopifnot(is.function(replication_function))
  list(formals  = formals(replication_function), body  = body(replication_function))
}
