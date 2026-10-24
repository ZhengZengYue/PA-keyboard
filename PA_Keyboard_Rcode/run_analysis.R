# =============================================================================
# PA-Keyboard: simulation study and reproducibility
# =============================================================================
# Purpose: reproduce the main study, sensitivity analyses, tables and Figures 1-2.
# Requirements: keep R/ beside this file; install openxlsx, dplyr, tidyr, ggplot2.
# Execution: open PA_Keyboard.Rproj, open this file, edit settings, click Source.
# Alternatively set project_dir to the full extracted project directory.
# Definitions in Sections 4-5 do not run simulations. Section 6 calls them for
# every replication in each scenario; Section 8 uses the same PA-Keyboard flow.
# Outputs: compact RDS records, checkpoints, XLSX tables, TIFFs and session info.
# To regenerate outputs without simulation, use render_results.R.
# =============================================================================

# 1. User settings ----------------------------------------------------------
project_dir <- "D:/thesis_and_paper/杨悦1/Biometrical Journal/PA_Keyboard_BJ_reproducibility_v6/PA_Keyboard_BJ_reproducibility_v6"
output_dir  <- NULL                 # NULL: results/ under project_dir; or full path
n_sim       <- 5000L                      # Replications per scenario
n_workers   <- 4L                   # Number of parallel worker processes
resume      <- TRUE                     # Reuse checkpoints with matching configurations

# 2. Load study modules -----------------------------------------------------
source(file.path(project_dir, "R", "00_setup.R"), encoding  = "UTF-8")
# Validate the project layout and resolve the output directory.
paths <- pak_initialize_paths(project_dir, output_dir)
# Load method, simulation, summary and plotting functions from R/.
pak_load_modules(paths)
# Workers use the entry script path to locate R/; they do not execute it.
PAK_ENTRY_SCRIPT <- paths$entry_script
output_dir       <- paths$output_dir
pak_check_packages()

# 3. Configure the study ----------------------------------------------------
# config contains the common trial design, borrowing rules and simulation settings.
# The keyed master seed is available as config$master_seed.
config <- make_formal_config(
  overrides  = list(
    n_sim       = n_sim,
    n_workers   = n_workers,
    write_excel = TRUE,
    verbose     = TRUE
  ),
  output_dir  = output_dir
)
config <- validate_formal_config(config)

# scenarios contains the true dose-specific DLT probabilities and AUC profiles.
scenarios             <- make_scenarios()
main_scenarios        <- c("E0", "E1", "E2", "E3", "M1", "M2", "L1", "L2", "S1", "S2")
sensitivity_scenarios <- FORMAL_SENSITIVITY_SCENARIOS

# 4a. Define the current-dose matching and borrowing workflow ----------------
# Section 4b calls this function after each pediatric cohort at its current dose.
# It returns the updated posterior and matching/borrowing diagnostics. Final-dose
# assessment also uses this function, with final = TRUE and its specified gate.
pa_keyboard_borrow_at_dose <- function(
  pediatric_dose,
  pediatric_n,
  pediatric_x,
  pediatric_gm_auc,
  adult_snapshot,
  scenario,
  potentials,
  config,
  final  = FALSE
) {

  config <- validate_pa_keyboard_variant_config(config)

  # Select the pediatric sample-size gate for interim or final assessment.
  min_n <- if (isTRUE(final)) {
    config$pa_keyboard_final_min_pediatric_n
  } else {
    config$pa_keyboard_min_pediatric_n
  }

  # Step 1: Identify a comparable adult dose from cumulative geometric-mean AUC.
  match_result <- select_auc_match(
    pediatric_gm       = pediatric_gm_auc,
    adult_gm           = adult_snapshot$gm_auc,
    adult_n            = adult_snapshot$n,
    adult_nominal_dose = config$adult_doses_mg,
    caliper            = config$auc_caliper_log,
    tolerance          = config$auc_tie_tolerance
  )

  # Start with zero adult borrowing and a pediatric-only posterior.
  # If a borrowing condition fails, this result is returned unchanged.
  out <- initialize_pa_keyboard_borrowing(match_result, pediatric_x, pediatric_n)

  # Step 2: Require sufficient pediatric data and an acceptable AUC match.
  if (pediatric_n < min_n || !isTRUE(match_result$available)) {
    return(out)
  }

  j <- out$matched_adult_dose
  if (is.na(j) || adult_snapshot$n[j] <= 0L) return(out)

  # Retrieve observed AUC and DLT pairs at the current and matched doses.
  ped_dat <- extract_observed_auc_dlt(
    "pediatric", pediatric_dose, pediatric_n, scenario, potentials, config
  )
  adult_dat <- extract_observed_auc_dlt(
    "adult", j, adult_snapshot$n[j], scenario, potentials, config
  )

  # Step 3: Assess local exposure-toxicity compatibility by posterior overlap.
  fit <- pa_keyboard_posterior_ovl(
    toxicity_p = ped_dat$toxicity,
    auc_p      = ped_dat$auc,
    toxicity_a = adult_dat$toxicity,
    auc_a      = adult_dat$auc,
    config     = config
  )

  # Step 4: Apply the prespecified operational overlap convention.
  delta_ovl <- pa_keyboard_operational_delta(fit$delta_raw, config, final  = final)

  # Step 5: Apply the adult ESS limit to the borrowing weight.
  capped <- pa_keyboard_apply_ess_cap(
    delta_ovl = delta_ovl,
    adult_n   = adult_snapshot$n[j],
    cap       = config$pa_keyboard_ess_cap
  )

  # Step 6: Update the pediatric posterior using weighted adult counts.
  posterior <- pa_keyboard_effective_posterior(
    pediatric_x = pediatric_x,
    pediatric_n = pediatric_n,
    adult_x     = adult_snapshot$x[j],
    adult_n     = adult_snapshot$n[j],
    delta       = capped$delta_applied
  )

  # Attach raw/applied weights, effective adult sample sizes and the posterior.
  complete_pa_keyboard_borrowing(out, fit, delta_ovl, capped, posterior)
}

# 4b. Define the sequential PA-Keyboard trial -------------------------------
# One call represents one pediatric trial. The cohort loop updates only the
# currently treated dose; safety uses pediatric DLT observations alone.
# The returned trial, dose and interim records feed the study summaries.
simulate_pa_keyboard_pediatric <- function(
  scenario,
  replication,
  potentials,
  adult_trial,
  config
) {

  config <- validate_pa_keyboard_variant_config(config)
  method <- PA_KEYBOARD_METHOD
  n_dose <- length(config$pediatric_doses_mg_m2)

  # Dose-level state: n = treated patients; x = observed DLTs.
  # AUC sums and counts provide the cumulative geometric mean at each dose.
  n           <- integer(n_dose)
  x           <- integer(n_dose)
  log_auc_sum <- numeric(n_dose)
  auc_n       <- integer(n_dose)
  eliminated  <- rep(FALSE, n_dose)
  current     <- 1L
  early_stop  <- FALSE
  stop_reason <- "pediatric_max_n"

  # Collect interim decisions, patient allocation, matching and borrowing records.
  interim_rows         <- list()
  path_rows            <- list()
  matching_rows        <- list()
  borrowing_rows       <- list()
  max_cohorts          <- config$pediatric_max_n %/% config$cohort_size
  final_adult_snapshot <- NULL

  for (cohort in seq_len(max_cohorts)) {
    dose <- current

    # Enroll the next pediatric cohort and update its observed outcomes.
    enrolled <- enroll_pediatric_cohort(
      dose, n, x, log_auc_sum, auc_n,
      scenario, potentials, config
    )
    n           <- enrolled$n
    x           <- enrolled$x
    log_auc_sum <- enrolled$log_auc_sum
    auc_n       <- enrolled$auc_n

    # Restrict adult information to what is available at this pediatric interim.
    child_gm       <- geometric_mean_from_logs(log_auc_sum, auc_n)
    adult_snapshot <- adult_snapshot_for_pediatric_interim(
      adult_trial, cohort, config
    )
    final_adult_snapshot <- adult_snapshot

    # Apply the matching-to-posterior workflow defined in Section 4a.
    borrow <- pa_keyboard_borrow_at_dose(
      pediatric_dose   = dose,
      pediatric_n      = n[dose],
      pediatric_x      = x[dose],
      pediatric_gm_auc = child_gm[dose],
      adult_snapshot   = adult_snapshot,
      scenario         = scenario,
      potentials       = potentials,
      config           = config,
      final            = FALSE
    )

    # Use the borrowed posterior for the current-dose Keyboard decision.
    decision <- keyboard_decision(
      x         = x[dose],
      n         = n[dose],
      config    = config,
      posterior = borrow$posterior
    )

    # Permanent dose elimination is based on pediatric DLT data only.
    safety     <- update_keyboard_elimination(n, x, eliminated, config)
    eliminated <- safety$eliminated

    # Safety stopping takes precedence over sample-size completion and movement.
    reached_maximum <- sum(n) >= config$pediatric_max_n
    if (safety$stop) {
      actual_action <- "STOP_SAFETY"
      next_dose     <- NA_integer_
      early_stop    <- TRUE
      stop_reason   <- "pediatric_lowest_dose_overdose"
    } else if (reached_maximum) {
      actual_action <- "END_MAX_N"
      next_dose     <- NA_integer_
      stop_reason   <- "pediatric_max_n"
    } else if (eliminated[dose]) {
      actual_action <- "D_SAFETY"
      next_dose     <- next_keyboard_dose(dose, "D", eliminated)
    } else {
      actual_action <- decision$action
      next_dose     <- next_keyboard_dose(dose, decision$action, eliminated)
    }

    # Record matching, borrowing and dose decisions for study summaries.
    rows <- record_pa_keyboard_interim(
      scenario, config, dose, borrow, method, replication, cohort, child_gm,
      adult_snapshot, n, x, decision, actual_action, next_dose, safety,
      early_stop, enrolled
    )
    matching_rows[[cohort]]  <- rows$matching
    borrowing_rows[[cohort]] <- rows$borrowing
    interim_rows[[cohort]]   <- rows$interim
    path_rows[[cohort]]      <- rows$path

    # End the trial when required; otherwise assign the next cohort.
    if (early_stop || sum(n) >= config$pediatric_max_n) break
    if (is.na(next_dose)) {
      early_stop  <- TRUE
      stop_reason <- "pediatric_no_assignable_dose"
      break
    }
    current <- next_dose
  }

  if (is.null(final_adult_snapshot)) {
    final_adult_snapshot <- adult_snapshot_for_pediatric_interim(
      adult_trial, 1L, config
    )
  }

  # Early-stopped trials have no selected MTD. Other trials apply the final
  # dose-selection rule implemented in select_pa_keyboard_mtd() (R/04_pa_keyboard.R).
  final_fit <- if (early_stop) {
    empty_pa_keyboard_final_fit(n_dose)
  } else {
    select_pa_keyboard_mtd(
      n              = n,
      x              = x,
      log_auc_sum    = log_auc_sum,
      auc_n          = auc_n,
      eliminated     = eliminated,
      adult_snapshot = final_adult_snapshot,
      scenario       = scenario,
      potentials     = potentials,
      config         = config,
      replication    = replication
    )
  }

  # Combine final dose estimates with the accumulated interim records.
  assemble_pa_keyboard_trial(
    scenario, replication, config, method, n, x, log_auc_sum, auc_n,
    early_stop, stop_reason, final_adult_snapshot, eliminated, final_fit,
    interim_rows, path_rows, matching_rows, borrowing_rows
  )
}

# 5. Define one three-design simulation replication ------------------------
# Section 6 passes this function to the scenario runner as replication_function.
# Each call generates one shared dataset and returns one result per design.
simulate_formal_replication <- function(scenario, replication, config) {
  config    <- validate_formal_config(config)
  scenarios <- make_scenarios()

  if (is.character(scenario)) scenario <- scenarios[[scenario]]
  if (is.null(scenario) || !(scenario$id %in% FORMAL_SCENARIOS)) {
    stop("Unknown/non-primary formal scenario.", call.  = FALSE)
  }

  # Use the fixed scenario ordering to identify reproducible random streams.
  scenario_index <- match(scenario$id, names(scenarios))

  # Generate shared adult/pediatric potential outcomes for all three designs.
  # Each of the four random domains uses its prespecified keyed seed.
  potentials <- make_potential_outcomes(scenario_index, replication, config)

  # Simulate the concurrent adult Keyboard trial once per replication.
  adult_trial <- simulate_adult_trial(scenario, potentials, config)

  # run_study_method() (R/05_simulation.R) runs the supplied trial expression,
  # records method-specific failures, and attaches the adult allocation path.
  # Each design uses the same potentials and adult trial within this replication.
  keyboard <- run_study_method(
    method     = "Keyboard",
    expression = simulate_keyboard_pediatric(
      "Keyboard", scenario, replication, potentials, adult_trial, config
    ),
    scenario    = scenario,
    replication = replication,
    adult_trial = adult_trial,
    config      = config
  )

  # formal_to_core_config() maps the formal-study configuration to the
  # PA-Keyboard variant fields required by the workflow in Section 4b.
  pa_keyboard <- run_study_method(
    method     = "PA-Keyboard",
    expression = simulate_pa_keyboard_pediatric(
      scenario, replication, potentials, adult_trial, formal_to_core_config(config)
    ),
    scenario    = scenario,
    replication = replication,
    adult_trial = adult_trial,
    config      = config
  )

  pa_crm <- run_study_method(
    method     = "PA-CRM",
    expression = simulate_pa_crm_pediatric(
      scenario, replication, potentials, adult_trial, config
    ),
    scenario    = scenario,
    replication = replication,
    adult_trial = adult_trial,
    config      = config
  )

  # Preserve the design labels used by the scenario-level summary functions.
  list("Keyboard"  = keyboard, "PA-Keyboard"  = pa_keyboard, "PA-CRM"  = pa_crm)
}

# 6. Run the main study across scenarios ------------------------------------
# run_formal_scenario_checkpointed() (R/05_simulation.R) repeats the Section 5
# workflow n_sim times per scenario, using serial or parallel execution.
# It saves scenario results and, with resume = TRUE, reuses a saved result only
# when both its configuration and replication-function definition match.
# main_pieces contains one scenario result, including summary tables, per entry.
main_pieces <- vector("list", length(main_scenarios))
names(main_pieces) <- main_scenarios

for (scenario_id in main_scenarios) {
  main_pieces[[scenario_id]] <- run_formal_scenario_checkpointed(
    scenario_id          = scenario_id,
    config               = config,
    replication_function = simulate_formal_replication,
    n_sim                = n_sim,
    n_workers            = n_workers,
    resume               = resume,
    script_path          = paths$entry_script
  )
}

# 7. Combine the main-study summaries ---------------------------------------
# Each scenario runner has summarized selection, allocation, safety, matching
# and borrowing. Assemble these tables across scenarios, write the results
# workbook, and save the compact main-study record for subsequent analyses.
raw_workbook <- file.path(output_dir, "PA_Keyboard_raw_results.xlsx")

main_results <- assemble_formal_main_results(
  pieces       = main_pieces,
  config       = config,
  scenario_ids = main_scenarios,
  n_sim        = n_sim,
  n_workers    = n_workers,
  output_xlsx  = raw_workbook
)

# 8. Run sensitivity analyses ----------------------------------------------
# sensitivity_plan describes the one-factor-at-a-time comparisons; variants
# contains their executable configurations. The settings vary the AUC caliper,
# adult ESS cap, pediatric borrowing gate and adult-information availability.
sensitivity_plan <- formal_sensitivity_plan_table()
variants         <- make_formal_sensitivity_variants(config)

sensitivity_pieces <- list()

# The checkpoint runner uses the methods specified for each variant and reuses
# the PA-Keyboard workflow in Sections 4a-4b. Results are indexed by both variant
# and scenario; compatible completed sensitivity runs can be resumed.
for (variant_id in names(variants)) {
  for (scenario_id in sensitivity_scenarios) {
    key <- paste(variant_id, scenario_id, sep  = "__")

    sensitivity_pieces[[key]] <- run_sensitivity_scenario_checkpointed(
      scenario_id = scenario_id,
      config      = variants[[variant_id]],
      n_sim       = n_sim,
      n_workers   = n_workers,
      resume      = resume,
      script_path = paths$entry_script
    )
  }
}

# Combine variant summaries with the main-study reference rows and update
# the sensitivity sheet in the results workbook.
sensitivity_results <- assemble_formal_sensitivity_results(
  pieces          = sensitivity_pieces,
  base_config     = config,
  main_results    = main_results,
  scenario_ids    = sensitivity_scenarios,
  n_sim           = n_sim,
  n_workers       = n_workers,
  variants        = variants,
  update_workbook = TRUE,
  output_xlsx     = raw_workbook
)

# 9. Save the complete analysis --------------------------------------------
# Keep main and sensitivity results together for independent table/figure
# regeneration with render_results.R; no new simulations are needed there.
analysis <- list(
  metadata  = list(
    code_version = FORMAL_FILE_VERSION,
    completed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    output_xlsx  = raw_workbook
  ),
  main        = main_results,
  sensitivity = sensitivity_results
)

saveRDS(
  analysis,
  file     = file.path(output_dir, "PA_keyboard_ETOVL_formal_full_compact_v03.rds"),
  compress = "gzip"
)

# 10. Generate publication tables and figures -------------------------------
# Define paths together so the publication outputs are easy to locate.
manuscript_file <- file.path(output_dir, "PA_Keyboard_manuscript_tables.xlsx")
figure1_file    <- file.path(output_dir, "Figure1_scenario_truth.tiff")
figure2_file    <- file.path(output_dir, "Figure2_dose_finding_performance.tiff")

# Tables use the complete analysis. Figure 1 displays scenario truths;
# Figure 2 summarizes dose-finding performance from the simulation results.
manuscript <- create_manuscript_tables(analysis, manuscript_file)
figure1    <- create_figure1(figure1_file)
figure2    <- create_figure2(analysis, figure2_file)

# Record the R version and loaded packages for reproducibility.
writeLines(
  capture.output(utils::sessionInfo()),
  con  = file.path(output_dir, "session_info.txt")
)
message("Analysis outputs: ", output_dir)
