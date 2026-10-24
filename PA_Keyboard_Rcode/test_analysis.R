# =============================================================================
# PA-Keyboard: algorithm and reproducibility checks
# =============================================================================
# Run this file with Source in RStudio or with Rscript test_analysis.R.
# Tests cover the formal settings, borrowing cap, simulation summaries,
# deterministic random streams, and a small three-design simulation.
# An optional parallel check can be requested with Rscript test_analysis.R --parallel.
# These small tests do not reproduce the manuscript's 5,000-replication values.
#
# Function guide:
#   test_replication_dispatch()        Verify callback dispatch and checkpoint identity.
#   test_workflow_loading()            Verify definitions-only loading.
#   unit_test_formal_config()          Verify primary design and adult ESS cap.
#   unit_test_formal_sensitivity_config()
#                                      Verify sensitivity-analysis settings.
#   unit_test_formal_borrowing_metrics()
#                                      Verify denominator conventions.
#   unit_test_formal()                 Run core checks and a small trial simulation.
#   test_seed_reproducibility()        Verify fixed RNG streams across repeats.
#   test_parallel_reproducibility()    Compare serial and PSOCK outputs.
# =============================================================================

# Open PA_Keyboard.Rproj first, or set a full project directory below.
test_project_dir <- "."
check_parallel   <- FALSE  # TRUE: also compare one and two PSOCK workers
source(file.path(test_project_dir, "R", "00_setup.R"), encoding  = "UTF-8")
paths <- pak_initialize_paths(test_project_dir, file.path(tempdir(), "PAK_tests"))
pak_load_modules(paths)
pak_load_workflows(paths$project_dir)
PAK_ENTRY_SCRIPT <- paths$entry_script

test_workflow_loading <- function() {
  definitions <- new.env(parent  = globalenv())
  pak_load_workflows(paths$project_dir, envir  = definitions)
  stopifnot(setequal(ls(definitions), c(
        "pa_keyboard_borrow_at_dose", "simulate_pa_keyboard_pediatric",
        "simulate_formal_replication")))
  message("PASS: workflow loading defines functions without running analysis.")
  invisible(TRUE)
}

# A deliberate sentinel error verifies that a supplied callback is actually
# invoked, rather than silently replaced by the default replication workflow.
test_replication_dispatch <- function(n_workers  = 1L) {
  cfg <- make_formal_config(
    overrides  = list(
      n_sim       = 2L,
      n_workers   = n_workers,
      fail_fast   = TRUE,
      write_excel = FALSE
    ),
    output_dir  = tempfile("PAK_callback_")
  )
  sentinel <- function(scenario, replication, config) {
    stop("PAK_CALLBACK_INVOKED", call.  = FALSE)
  }
  error_text <- tryCatch({
      run_formal_scenario_checkpointed(
        scenario_id          = "E1",
        config               = cfg,
        replication_function = sentinel,
        n_sim                = 2L,
        n_workers            = n_workers,
        resume               = FALSE
      )
      "NO_ERROR"
    }, error  = conditionMessage)
  stopifnot(grepl("PAK_CALLBACK_INVOKED", error_text, fixed  = TRUE))

  # Callback definitions are part of checkpoint compatibility.
  other <- function(scenario, replication, config) NULL
  stopifnot(!identical(formal_replication_signature(sentinel),
      formal_replication_signature(other)))
  message("PASS: supplied replication callback invoked with ", n_workers, " worker(s).")
  invisible(TRUE)
}

unit_test_formal_config <- function() {
  cfg <- make_formal_config(
    overrides  = list(
      n_sim       = 5L,
      n_workers   = 1L,
      write_excel = FALSE
    )
  )

  stopifnot(
    cfg$pa_keyboard_min_pediatric_n == 6L,
    cfg$pa_keyboard_final_min_pediatric_n == 6L,
    cfg$pa_keyboard_ess_cap == 6,
    cfg$safety_min_n == 3L
  )

  cap_check <- pa_keyboard_apply_ess_cap(
    delta_ovl = 0.8,
    adult_n   = 15,
    cap       = cfg$pa_keyboard_ess_cap
  )
  stopifnot(
    abs(cap_check$effective_adult_n - 6) < 1e-12,
    abs(cap_check$delta_applied - 0.4) < 1e-12,
    isTRUE(cap_check$cap_triggered)
  )

  message("PASS: formal config + ESS-cap checks.")
  invisible(TRUE)
}

unit_test_formal_sensitivity_config <- function() {
  cfg <- make_formal_config(
    overrides  = list(
      n_sim       = 3L,
      n_workers   = 1L,
      write_excel = FALSE
    ),
    output_dir  = file.path(tempdir(), "ETOVL_formal_sensitivity_test")
  )

  variants <- make_formal_sensitivity_variants(cfg)
  stopifnot(
    identical(names(variants), FORMAL_SENSITIVITY_PLAN$variant_id),
    variants$CAL16$auc_caliper_log == 0.16,
    variants$CAL28$auc_caliper_log == 0.28,
    variants$LEAD9$adult_lead_n == 9L,
    variants$LEAD12$adult_lead_n == 12L,
    isTRUE(variants$ADULT_COMPLETED$formal_adult_completed),
    variants$CAP3$pa_keyboard_ess_cap == 3,
    variants$CAP9$pa_keyboard_ess_cap == 9,
    variants$GATE3$pa_keyboard_min_pediatric_n == 3L,
    variants$GATE3$pa_keyboard_final_min_pediatric_n == 3L
  )

  # Completed-adult mode must expose the final available adult snapshot from
  # the first pediatric interim onward.
  scenarios   <- make_scenarios()
  potentials  <- make_potential_outcomes(2L, 1L, variants$ADULT_COMPLETED)
  adult_trial <- simulate_adult_trial(
    scenarios$E1,
    potentials,
    variants$ADULT_COMPLETED
  )
  snap <- adult_snapshot_for_pediatric_interim(
    adult_trial,
    pediatric_cohort = 1L,
    config           = variants$ADULT_COMPLETED
  )
  stopifnot(
    snap$available_adult_cohort == length(adult_trial$snapshots),
    isTRUE(snap$adult_lead_completed)
  )

  message("PASS: formal sensitivity configuration checks.")
  invisible(TRUE)
}

unit_test_formal_borrowing_metrics <- function() {
  toy_dose <- data.frame(
    scenario                            = rep("E1", 6),
    method                              = rep("PA-Keyboard", 6),
    replication                         = rep(1:3, each  = 2),
    pediatric_dose                      = rep(1:2, times  = 3),
    n                                   = c(6, 0, 6, 3, 6, 0),
    final_pava_weight                   = c(9, NA, 9, NA, 9, NA),
    final_matched_adult_dose            = c(2, NA, 2, NA, NA, NA),
    final_pa_keyboard_delta_raw         = c(0.5, 0, 0.6, 0, 0, 0),
    final_pa_keyboard_delta_ovl         = c(0.5, 0, 0.6, 0, 0, 0),
    final_pa_keyboard_delta_applied     = c(0.5, 0, 0.5, 0, 0, 0),
    final_pa_keyboard_effective_adult_n = c(3, 0, 3, 0, 0, 0),
    final_pa_keyboard_cap_triggered     = c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE),
    stringsAsFactors                    = FALSE
  )

  z  <- make_formal_final_borrowing_table(list(dose  = toy_dose))
  z1 <- z[z$Dose == "P1", , drop  = FALSE]
  z2 <- z[z$Dose == "P2", , drop  = FALSE]

  stopifnot(
    nrow(z1) == 1L,
    abs(z1$`Dose visited (%)` - 100) < 1e-12,
    abs(z1$`Borrow activated (%) | all trials` - (200 / 3)) < 1e-10,
    abs(z1$`Mean realized adult ESS | all trials` - 2) < 1e-12,
    abs(z1$`Median realized adult ESS | all trials` - 3) < 1e-12,
    abs(z1$`Median adult ESS | activated` - 3) < 1e-12,
    abs(z2$`Mean realized adult ESS | all trials`) < 1e-12,
    abs(z2$`Median realized adult ESS | all trials`) < 1e-12
  )

  toy_match <- data.frame(
    scenario                       = rep("E1", 4),
    method                         = rep("PA-Keyboard", 4),
    replication                    = 1:4,
    pediatric_dose                 = rep(2L, 4),
    oracle_eligible                = c(3L, 3L, 3L, NA_integer_),
    match_available                = c(TRUE, TRUE, FALSE, TRUE),
    correct_given_oracle_available = c(TRUE, FALSE, FALSE, NA),
    match_distance                 = c(0.01, 0.02, NA, 0.03),
    stringsAsFactors               = FALSE
  )

  m <- make_formal_matching_table(list(matching  = toy_match))
  stopifnot(
    nrow(m) == 1L,
    abs(
      m$`Correct AUC match (%) | oracle and sample match available` - 50
    ) < 1e-10,
    abs(
      m$`Sample match available (%) | oracle available` - (200 / 3)
    ) < 1e-10
  )

  message(
    "PASS: unconditional realized-ESS + conditional matching denominator checks."
  )
  invisible(TRUE)
}

unit_test_formal <- function() {
  unit_test_formal_config()
  unit_test_formal_sensitivity_config()
  unit_test_formal_borrowing_metrics()

  cfg <- make_formal_config(
    overrides  = list(
      n_sim       = 3L,
      n_workers   = 1L,
      write_excel = FALSE,
      fail_fast   = TRUE
    ),
    output_dir  = file.path(tempdir(), "ETOVL_formal_test")
  )

  ans <- run_formal_scenario(
    "E1",
    config      = cfg,
    n_sim       = 3L,
    n_workers   = 1L,
    script_path = detect_formal_script_path()
  )

  stopifnot(
    setequal(unique(ans$main$Method), FORMAL_METHODS),
    !"Proposed_PA_Keyboard" %in% ans$main$Method,
    nrow(ans$dose) > 0,
    nrow(ans$matching) > 0,
    nrow(ans$online_borrowing) > 0,
    nrow(ans$final_borrowing) > 0,
    "Mean realized adult ESS | all trials" %in% names(ans$final_borrowing)
  )

  message("PASS: formal three-method smoke test.")
  invisible(TRUE)
}


# Verify determinism of every random domain without altering the study seed.
test_seed_reproducibility <- function() {
  cfg <- make_formal_config(overrides  = list(n_sim  = 3L, n_workers  = 1L))
  x   <- make_potential_outcomes(2L, 1L, cfg)
  y   <- make_potential_outcomes(2L, 1L, cfg)
  z   <- make_potential_outcomes(2L, 2L, cfg)
  stopifnot(identical(x, y), !identical(x, z))
  stopifnot(length(unique(vapply(1:4, function(i)
          seed_from_keys(cfg$master_seed, 2L, 1L, i), integer(1)))) == 4L)
  message("PASS: repeated seeds, replication separation and four RNG domains.")
  invisible(TRUE)
}

# Optional: same scenario, identical outputs using serial and PSOCK parallel.
test_parallel_reproducibility <- function() {
  cfg <- make_formal_config(overrides  = list(
      n_sim       = 3L,
      n_workers   = 2L,
      fail_fast   = TRUE,
      write_excel = FALSE
    ), output_dir  = file.path(tempdir(), "PAK_reproducibility_test"))
  a <- run_formal_scenario("E1", cfg, n_sim  = 3L, n_workers  = 1L)
  b <- run_formal_scenario("E1", cfg, n_sim  = 3L, n_workers  = 2L)
  for (nm in c("main", "dose", "matching", "final_borrowing", "online_borrowing")) {
    x <- a[[nm]]
    y <- b[[nm]]
    if (nm == "main") {
      x$`Scenario runtime (min)` <- NULL
      y$`Scenario runtime (min)` <- NULL
    }
    stopifnot(isTRUE(all.equal(x, y, tolerance  = 0, check.attributes  = TRUE)))
  }
  message("PASS: identical results with 1 and 2 workers.")
  invisible(TRUE)
}

test_workflow_loading()
test_replication_dispatch()
test_seed_reproducibility()
unit_test_formal()
if (isTRUE(check_parallel) || "--parallel" %in% commandArgs(trailingOnly  = TRUE)) {
  test_replication_dispatch(n_workers  = 2L)
  test_parallel_reproducibility()
}
message("All requested tests passed.")
