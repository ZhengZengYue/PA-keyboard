# =============================================================================
# Operating characteristics, manuscript tables, and result workbooks
# =============================================================================
# This module summarizes simulated dose selection, patient allocation, safety,
# and borrowing outcomes. It creates the manuscript and supplementary tables
# and coordinates generation of the complete publication outputs.
#
# Function guide:
#   split_rows(), group_values()          Group simulation results.
#   summarize_trial_table(),
#     summarize_dose_table()             Aggregate trial/dose operating measures.
#   make_formal_main_table(),
#     make_formal_dose_table()           Tabulate simulation performance.
#   make_formal_matching_table(),
#     make_formal_online_borrowing_table(),
#     make_formal_final_borrowing_table() Summarize AUC matching and borrowed ESS.
#   formal_config_table(),
#     formal_scenario_table()           Tabulate simulation specifications.
#   write_formal_workbook()              Write full simulation summary workbook.
#   create_manuscript_tables()           Write Tables 1-4 and Tables B-1 to B-5.
#   run_publication_analysis()           Programmatic complete-analysis wrapper.
#   rebuild_publication_outputs()        Recreate tables/figures from saved RDS.
#   Other formatting helpers             Prepare tables and Excel styles.
# =============================================================================

# =============================================================================
# Summaries of simulated operating characteristics
# =============================================================================
split_rows <- function(data, columns) {
  if (!nrow(data)) return(list())
  key <- interaction(data[columns], drop  = TRUE, lex.order  = TRUE, sep  = "\r")
  split(data, key, drop  = TRUE)
}

group_values <- function(data, columns) {
  as.list(data[1L, columns, drop  = FALSE])
}

summarize_trial_table <- function(trial) {
  groups <- split_rows(trial, c("config_id", "scenario", "method"))
  rbind_fill(lapply(groups, function(d) {
        key <- group_values(d, c("config_id", "scenario", "method"))
        data.frame(
          config_id                                = key$config_id,
          scenario                                 = key$scenario,
          method                                   = key$method,
          successful_trials                        = nrow(d),
          PCS                                      = safe_mean(d$correct_selection),
          PCS_mcse                                 = safe_mcse(d$correct_selection),
          correct_no_mtd_probability               = safe_mean(d$correct_no_mtd),
          correct_no_mtd_mcse                      = safe_mcse(d$correct_no_mtd),
          selected_no_mtd_probability              = safe_mean(d$selected_no_mtd),
          early_stop_probability                   = safe_mean(d$early_stop),
          early_stop_mcse                          = safe_mcse(d$early_stop),
          mean_pediatric_n                         = safe_mean(d$pediatric_n),
          mean_total_dlt                           = safe_mean(d$total_dlt),
          mean_total_dlt_mcse                      = safe_mcse(d$total_dlt),
          mean_overtoxic_allocation                = safe_mean(d$overtoxic_allocation),
          mean_overtoxic_allocation_mcse           = safe_mcse(d$overtoxic_allocation),
          probability_overtoxic_allocation_ge_0_60 = safe_mean(
            d$overtoxic_allocation >= 0.60
          ),
          mean_true_mtd_allocation      = safe_mean(d$true_mtd_allocation),
          mean_adult_n_at_pediatric_end = safe_mean(d$adult_n_at_pediatric_end),
          stringsAsFactors              = FALSE
        )
  }))
}

summarize_dose_table <- function(dose) {
  groups <- split_rows(dose, c("config_id", "scenario", "method", "pediatric_dose"))
  rbind_fill(lapply(groups, function(d) {
        key             <- group_values(d, c("config_id", "scenario", "method", "pediatric_dose"))
        final_evaluated <- is.finite(d$final_pava_weight)
        data.frame(
          config_id                        = key$config_id,
          scenario                         = key$scenario,
          method                           = key$method,
          pediatric_dose                   = key$pediatric_dose,
          dose_mg_m2                       = d$dose_mg_m2[1L],
          true_dlt                         = d$true_dlt[1L],
          true_auc                         = d$true_auc[1L],
          mean_n                           = safe_mean(d$n),
          mean_allocation                  = safe_mean(d$allocation),
          mean_dlt                         = safe_mean(d$dlt),
          selection_probability            = safe_mean(d$selected),
          elimination_probability          = safe_mean(d$eliminated),
          mean_final_posterior_median      = safe_mean(d$posterior_median),
          mean_final_isotonic_estimate     = safe_mean(d$isotonic_estimate),
          mean_final_pava_weight           = safe_mean(d$final_pava_weight),
          mean_final_borrowed_ess_positive = safe_mean(
            d$final_borrowed_ess_positive
          ),
          final_match_rate  = if (any(final_evaluated)) {
            mean(!is.na(d$final_matched_adult_dose[final_evaluated]))
          } else NA_real_,
          overtoxic        = d$overtoxic[1L],
          stringsAsFactors = FALSE
        )
  }))
}


pct <- function(x) {
  ifelse(is.na(x), NA_real_, 100 * x)
}

sort_formal_rows <- function(df) {
  if (!nrow(df)) return(df)
  s <- match(df$Scenario, FORMAL_SCENARIOS)
  m <- match(df$Method, FORMAL_METHODS)
  df[order(s, m), , drop  = FALSE]
}

make_formal_main_table <- function(raw, runtime_sec  = NA_real_) {
  tab <- summarize_trial_table(raw$trial)
  if (!nrow(tab)) return(data.frame())

  out <- data.frame(
    Scenario = tab$scenario,
    Method   = tab$method,
    `Successful trials`  = tab$successful_trials,
    `PCS (%)`  = pct(tab$PCS),
    `Correct no-MTD (%)`  = pct(tab$correct_no_mtd_probability),
    `Early stop (%)`  = pct(tab$early_stop_probability),
    `Patients at true MTD (%)`  = pct(tab$mean_true_mtd_allocation),
    `Patients at overly toxic doses (%)`  = pct(tab$mean_overtoxic_allocation),
    `P(overtoxic allocation >=60%)`  = 
    pct(tab$probability_overtoxic_allocation_ge_0_60),
    `Mean pediatric N`  = tab$mean_pediatric_n,
    `Mean DLT`  = tab$mean_total_dlt,
    `Mean adult N at pediatric end`  = tab$mean_adult_n_at_pediatric_end,
    check.names      = FALSE,
    stringsAsFactors = FALSE
  )

  failures <- if (nrow(raw$failure)) {
    aggregate(
      replication ~ scenario + method,
      raw$failure,
      FUN  = length
    )
  } else {
    data.frame(
      scenario    = character(),
      method      = character(),
      replication = integer()
    )
  }
  names(failures)[names(failures) == "replication"] <- "Failures"

  key_out      <- paste(out$Scenario, out$Method, sep  = "||")
  key_fail     <- paste(failures$scenario, failures$method, sep  = "||")
  idx          <- match(key_out, key_fail)
  out$Failures <- ifelse(is.na(idx), 0L, failures$Failures[idx])
  out$`Scenario runtime (min)` <- runtime_sec / 60

  sort_formal_rows(out)
}

make_formal_dose_table <- function(raw) {
  tab <- summarize_dose_table(raw$dose)
  if (!nrow(tab)) return(data.frame())

  out <- data.frame(
    Scenario = tab$scenario,
    Method   = tab$method,
    Dose     = paste0("P", tab$pediatric_dose),
    `Dose (mg/m2)`  = tab$dose_mg_m2,
    `True pediatric DLT`  = tab$true_dlt,
    `True pediatric AUC`  = tab$true_auc,
    `Selection (%)`  = pct(tab$selection_probability),
    `Mean N treated`  = tab$mean_n,
    `Allocation (%)`  = pct(tab$mean_allocation),
    `Mean DLT at dose`  = tab$mean_dlt,
    `Elimination (%)`  = pct(tab$elimination_probability),
    Overtoxic        = tab$overtoxic,
    check.names      = FALSE,
    stringsAsFactors = FALSE
  )

  s <- match(out$Scenario, FORMAL_SCENARIOS)
  m <- match(out$Method, FORMAL_METHODS)
  d <- as.integer(sub("^P", "", out$Dose))
  out[order(s, m, d), , drop  = FALSE]
}

formal_safe_quantile <- function(x, prob) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(stats::quantile(
      x,
      probs = prob,
      names = FALSE,
      type  = 7,
      na.rm = TRUE
  ))
}

make_formal_matching_table <- function(raw) {
  d <- raw$matching
  if (!is.data.frame(d) || !nrow(d)) return(data.frame())

  d <- d[d$method == "PA-Keyboard", , drop  = FALSE]
  if (!nrow(d)) return(data.frame())

  key <- interaction(
    d$scenario,
    d$pediatric_dose,
    drop      = TRUE,
    lex.order = TRUE
  )
  groups <- split(d, key, drop  = TRUE)

  out <- rbind_fill(lapply(groups, function(z) {
        oracle_available <- !is.na(z$oracle_eligible)
        sample_match     <- !is.na(z$match_available) & z$match_available
        correct          <- !is.na(z$correct_given_oracle_available) &
        z$correct_given_oracle_available
        correct_evaluable <- oracle_available & sample_match
        matched_rows      <- sample_match & is.finite(z$match_distance)

        data.frame(
          Scenario = z$scenario[1L],
          Dose     = paste0("P", z$pediatric_dose[1L]),
          `N online matching evaluations`  = nrow(z),
          `N trials represented`  = length(unique(z$replication)),
          `Oracle comparable opportunity (%)`  = 100 * mean(oracle_available),
          `Sample AUC match available (%)`  = 100 * mean(sample_match),
          `Correct AUC match (%) | oracle and sample match available`  = 
          if (any(correct_evaluable)) {
            100 * mean(correct[correct_evaluable])
          } else {
            NA_real_
          },
          `Sample match available (%) | oracle available`  = 
          if (any(oracle_available)) {
            100 * mean(sample_match[oracle_available])
          } else {
            NA_real_
          },
          `Mean absolute log-AUC distance | matched`  = if (any(matched_rows)) {
            mean(z$match_distance[matched_rows])
          } else {
            NA_real_
          },
          check.names      = FALSE,
          stringsAsFactors = FALSE
        )
  }))

  dose_num <- as.integer(sub("^P", "", out$Dose))
  out[
    order(match(out$Scenario, FORMAL_SCENARIOS), dose_num),
    ,
    drop  = FALSE
  ]
}

make_formal_online_borrowing_table <- function(raw) {
  d <- raw$borrowing
  if (!is.data.frame(d) || !nrow(d)) return(data.frame())

  d <- d[d$method == "PA-Keyboard", , drop  = FALSE]
  if (!nrow(d)) return(data.frame())

  groups <- split(d, d$scenario, drop  = TRUE)
  out    <- rbind_fill(lapply(groups, function(z) {
        evaluated <- !is.na(z$ovl_evaluated) & z$ovl_evaluated

        data.frame(
          Scenario  = z$scenario[1L],
          `N online borrowing evaluations`  = nrow(z),
          `OVL evaluated (%)`  = 100 * mean(evaluated),
          `Mean raw OVL | evaluated`  = if (any(evaluated)) {
            mean(z$delta_raw[evaluated])
          } else {
            NA_real_
          },
          `Mean operational OVL delta | evaluated`  = if (any(evaluated)) {
            mean(z$delta_ovl[evaluated])
          } else {
            NA_real_
          },
          `Mean applied delta | evaluated`  = if (any(evaluated)) {
            mean(z$delta_applied[evaluated])
          } else {
            NA_real_
          },
          `Mean effective adult ESS | evaluated`  = if (any(evaluated)) {
            mean(z$effective_adult_n[evaluated])
          } else {
            NA_real_
          },
          `Cap triggered (%) | evaluated`  = if (any(evaluated)) {
            100 * mean(z$cap_triggered[evaluated])
          } else {
            NA_real_
          },
          check.names      = FALSE,
          stringsAsFactors = FALSE
        )
  }))

  out[order(match(out$Scenario, FORMAL_SCENARIOS)), , drop  = FALSE]
}

make_formal_final_borrowing_table <- function(raw) {
  d <- raw$dose
  if (!is.data.frame(d) || !nrow(d)) return(data.frame())

  d <- d[d$method == "PA-Keyboard", , drop  = FALSE]
  if (!nrow(d)) return(data.frame())

  required <- c(
    "scenario", "replication", "pediatric_dose", "n",
    "final_pava_weight",
    "final_matched_adult_dose",
    "final_pa_keyboard_delta_raw",
    "final_pa_keyboard_delta_ovl",
    "final_pa_keyboard_delta_applied",
    "final_pa_keyboard_effective_adult_n",
    "final_pa_keyboard_cap_triggered"
  )
  missing_cols <- setdiff(required, names(d))
  if (length(missing_cols)) {
    stop(
      paste0(
        "Cannot construct realized borrowing table; raw dose output is missing: ",
        paste(missing_cols, collapse  = ", ")
      ),
      call.  = FALSE
    )
  }

  # A non-borrowing state contributes an adult ESS of zero.
  d$realized_ess <- as.numeric(d$final_pa_keyboard_effective_adult_n)
  d$realized_ess[!is.finite(d$realized_ess)] <- 0
  d$realized_ess[d$realized_ess < 0] <- 0

  d$visited               <- is.finite(d$n) & d$n > 0
  d$final_eligible        <- is.finite(d$final_pava_weight)
  d$final_match_available <- !is.na(d$final_matched_adult_dose)
  d$borrow_activated      <- d$realized_ess > 1e-12

  key <- interaction(
    d$scenario,
    d$pediatric_dose,
    drop      = TRUE,
    lex.order = TRUE
  )
  groups <- split(d, key, drop  = TRUE)

  out <- rbind_fill(lapply(groups, function(z) {
        visited    <- z$visited
        eligible   <- z$final_eligible
        activated  <- z$borrow_activated
        active_ess <- z$realized_ess[activated]

        data.frame(
          Scenario = z$scenario[1L],
          Dose     = paste0("P", z$pediatric_dose[1L]),
          `N successful trials`  = length(unique(z$replication)),
          `Dose visited (%)`  = 100 * mean(visited),
          `Final eligible (%) | all trials`  = 100 * mean(eligible),
          `Final AUC match available (%) | visited`  = if (any(visited)) {
            100 * mean(z$final_match_available[visited])
          } else {
            NA_real_
          },
          `Borrow activated (%) | visited`  = if (any(visited)) {
            100 * mean(activated[visited])
          } else {
            NA_real_
          },
          `Borrow activated (%) | all trials`  = 100 * mean(activated),
          `Mean realized adult ESS | all trials`  = 
          mean(z$realized_ess),
          `Median realized adult ESS | all trials`  = 
          formal_safe_quantile(z$realized_ess, 0.50),
          `Q25 realized adult ESS | all trials`  = 
          formal_safe_quantile(z$realized_ess, 0.25),
          `Q75 realized adult ESS | all trials`  = 
          formal_safe_quantile(z$realized_ess, 0.75),
          `Mean adult ESS | activated`  = if (any(activated)) {
            mean(active_ess)
          } else {
            NA_real_
          },
          `Median adult ESS | activated`  = if (any(activated)) {
            formal_safe_quantile(active_ess, 0.50)
          } else {
            NA_real_
          },
          `Q25 adult ESS | activated`  = if (any(activated)) {
            formal_safe_quantile(active_ess, 0.25)
          } else {
            NA_real_
          },
          `Q75 adult ESS | activated`  = if (any(activated)) {
            formal_safe_quantile(active_ess, 0.75)
          } else {
            NA_real_
          },
          `Cap triggered (%) | activated`  = if (any(activated)) {
            100 * mean(
              as.logical(z$final_pa_keyboard_cap_triggered[activated]),
              na.rm  = TRUE
            )
          } else {
            NA_real_
          },
          check.names      = FALSE,
          stringsAsFactors = FALSE
        )
  }))

  dose_num <- as.integer(sub("^P", "", out$Dose))
  out[
    order(match(out$Scenario, FORMAL_SCENARIOS), dose_num),
    ,
    drop  = FALSE
  ]
}

formal_config_table <- function(config) {
  values <- c(
    Code_version                    = FORMAL_FILE_VERSION,
    Protocol_date                   = FORMAL_PROTOCOL_DATE,
    Methods                         = paste(FORMAL_METHODS, collapse  = ", "),
    Scenarios                       = paste(FORMAL_SCENARIOS, collapse  = ", "),
    Target_DLT                      = config$target,
    Target_key                      = paste0("(", config$target_lower, ", ", config$target_upper, ")"),
    Cohort_size                     = config$cohort_size,
    Adult_max_N                     = config$adult_max_n,
    Pediatric_max_N                 = config$pediatric_max_n,
    Adult_lead_N                    = config$adult_lead_n,
    Adult_data_mode                 = config$formal_adult_data_mode,
    Keyboard_online_prior           = "Beta(1,1)",
    Keyboard_safety_min_N           = config$safety_min_n,
    Keyboard_safety_cutoff          = config$safety_cutoff,
    AUC_DGP                         = "Direct log-normal patient-level AUC",
    AUC_sigma_log                   = config$auc_sigma_log,
    AUC_match_metric                = "absolute log-ratio of cumulative geometric means",
    AUC_caliper_log                 = config$auc_caliper_log,
    ET_model                        = "P(DLT=1|AUC,theta)=1-exp(-theta*AUC)",
    ET_prior                        = "q(A*)~Beta(1,1); theta~Exponential(rate=A*)",
    Borrowing_weight                = "posterior overlap coefficient (OVL)",
    Borrowing_start_pediatric_N     = config$pa_keyboard_min_pediatric_n,
    Final_borrowing_min_pediatric_N = config$pa_keyboard_final_min_pediatric_n,
    Maximum_borrowed_adult_ESS      = config$pa_keyboard_ess_cap,
    Realized_ESS_reporting          = 
    "All successful trials x dose; no borrowing contributes ESS=0",
    Conditional_ESS_reporting  = 
    "ESS distribution among final realized ESS > 0",
    Correct_AUC_match_reporting  = 
    paste0(
      "Online matching evaluations; denominator = evaluations with ",
      "both oracle and sample AUC match available"
    ),
    Hard_safety_uses_adult_data = "No; pediatric-only",
    Final_candidate_rule        = "tried + non-eliminated pediatric doses",
    N_sim_per_scenario          = config$n_sim,
    N_workers                   = config$n_workers,
    Master_seed                 = config$master_seed
  )

  data.frame(
    Parameter        = names(values),
    Value            = unname(values),
    stringsAsFactors = FALSE
  )
}

formal_scenario_table <- function() {
  scenarios <- make_scenarios()
  tab       <- scenario_truth_table(scenarios)
  tab       <- tab[tab$scenario %in% FORMAL_SCENARIOS, , drop  = FALSE]
  tab       <- tab[match(FORMAL_SCENARIOS, tab$scenario), , drop  = FALSE]
  tab
}

write_formal_workbook <- function(
  config,
  main_table,
  dose_table,
  matching_table,
  final_borrowing,
  online_borrowing,
  sensitivity_table  = data.frame(),
  path
) {

  if (!requireNamespace("openxlsx", quietly  = TRUE)) {
    stop(
      "Package 'openxlsx' is required for the formal XLSX output.",
      call.  = FALSE
    )
  }

  dir.create(dirname(path), recursive  = TRUE, showWarnings  = FALSE)

  wb <- openxlsx::createWorkbook()

  # Simple restrained styles.
  title_style <- openxlsx::createStyle(
    textDecoration = "bold",
    fontSize       = 12
  )
  header_style <- openxlsx::createStyle(
    textDecoration = "bold",
    halign         = "center",
    valign         = "center",
    border         = "Bottom"
  )
  pct_style <- openxlsx::createStyle(numFmt  = "0.00")
  num_style <- openxlsx::createStyle(numFmt  = "0.000")

  # --------------------------------------------------------------------------
  # 00_Design: two compact tables on one sheet.
  # --------------------------------------------------------------------------
  openxlsx::addWorksheet(wb, "00_Design")
  cfg  <- formal_config_table(config)
  scen <- formal_scenario_table()

  openxlsx::writeData(wb, "00_Design", "Formal configuration", startRow  = 1)
  openxlsx::addStyle(wb, "00_Design", title_style, rows  = 1, cols  = 1)
  openxlsx::writeData(wb, "00_Design", cfg, startRow  = 2, headerStyle  = header_style)

  scen_row <- nrow(cfg) + 5L
  openxlsx::writeData(wb, "00_Design", "Scenario truths", startRow  = scen_row)
  openxlsx::addStyle(wb, "00_Design", title_style, rows  = scen_row, cols  = 1)
  openxlsx::writeData(
    wb, "00_Design", scen,
    startRow    = scen_row + 1L,
    headerStyle = header_style
  )

  sens_plan     <- formal_sensitivity_plan_table()
  sens_plan_row <- scen_row + nrow(scen) + 4L
  openxlsx::writeData(
    wb, "00_Design", "Pre-specified OFAT sensitivity plan",
    startRow  = sens_plan_row
  )
  openxlsx::addStyle(
    wb, "00_Design", title_style,
    rows  = sens_plan_row, cols  = 1
  )
  openxlsx::writeData(
    wb, "00_Design", sens_plan,
    startRow    = sens_plan_row + 1L,
    headerStyle = header_style
  )
  openxlsx::setColWidths(
    wb, "00_Design",
    cols   = 1:max(2, ncol(scen), ncol(sens_plan)),
    widths = "auto"
  )

  # --------------------------------------------------------------------------
  # 01_Main
  # --------------------------------------------------------------------------
  openxlsx::addWorksheet(wb, "01_Main")
  openxlsx::writeData(wb, "01_Main", main_table, headerStyle  = header_style)
  openxlsx::freezePane(wb, "01_Main", firstRow  = TRUE)
  openxlsx::setColWidths(wb, "01_Main", cols  = seq_len(ncol(main_table)), widths  = "auto")

  # --------------------------------------------------------------------------
  # 02_Dose
  # --------------------------------------------------------------------------
  openxlsx::addWorksheet(wb, "02_Dose")
  openxlsx::writeData(wb, "02_Dose", dose_table, headerStyle  = header_style)
  openxlsx::freezePane(wb, "02_Dose", firstRow  = TRUE)
  openxlsx::setColWidths(wb, "02_Dose", cols  = seq_len(ncol(dose_table)), widths  = "auto")

  # --------------------------------------------------------------------------
  # 03_Borrowing
  #   A) AUC matching performance: WHERE to borrow
  #   B) Final realized ESS: HOW MUCH was actually borrowed
  #      - ALL successful replications are in the denominator
  #      - no borrowing = ESS 0
  #   C) Compact online OVL diagnostics: supplementary mechanism description
  # --------------------------------------------------------------------------
  openxlsx::addWorksheet(wb, "03_Borrowing")

  row_a <- 1L
  openxlsx::writeData(
    wb, "03_Borrowing",
    "A. Dose-local AUC matching performance (online evaluations)",
    startRow  = row_a
  )
  openxlsx::addStyle(wb, "03_Borrowing", title_style, rows  = row_a, cols  = 1)
  if (nrow(matching_table)) {
    openxlsx::writeData(
      wb, "03_Borrowing", matching_table,
      startRow    = row_a + 1L,
      headerStyle = header_style
    )
  }

  row_b <- row_a + max(nrow(matching_table), 1L) + 4L
  openxlsx::writeData(
    wb, "03_Borrowing",
    paste0(
      "B. Final realized borrowing by pediatric dose ",
      "(all successful trials; no borrowing = ESS 0)"
    ),
    startRow  = row_b
  )
  openxlsx::addStyle(wb, "03_Borrowing", title_style, rows  = row_b, cols  = 1)
  if (nrow(final_borrowing)) {
    openxlsx::writeData(
      wb, "03_Borrowing", final_borrowing,
      startRow    = row_b + 1L,
      headerStyle = header_style
    )
  }

  row_c <- row_b + max(nrow(final_borrowing), 1L) + 4L
  openxlsx::writeData(
    wb, "03_Borrowing",
    paste0(
      "C. Supplementary online ET-OVL diagnostics ",
      "(conditional on OVL evaluation where indicated)"
    ),
    startRow  = row_c
  )
  openxlsx::addStyle(wb, "03_Borrowing", title_style, rows  = row_c, cols  = 1)
  if (nrow(online_borrowing)) {
    openxlsx::writeData(
      wb, "03_Borrowing", online_borrowing,
      startRow    = row_c + 1L,
      headerStyle = header_style
    )
  }

  max_borrow_cols <- max(
    ncol(matching_table),
    ncol(final_borrowing),
    ncol(online_borrowing),
    1L
  )
  openxlsx::setColWidths(
    wb, "03_Borrowing",
    cols   = seq_len(max_borrow_cols),
    widths = "auto"
  )
  openxlsx::freezePane(wb, "03_Borrowing", firstRow  = TRUE)

  # --------------------------------------------------------------------------
  # 04_Sensitivity: all OFAT settings on ONE compact sheet.
  # --------------------------------------------------------------------------
  openxlsx::addWorksheet(wb, "04_Sensitivity")
  if (nrow(sensitivity_table)) {
    openxlsx::writeData(
      wb, "04_Sensitivity", sensitivity_table,
      headerStyle  = header_style
    )
    openxlsx::freezePane(wb, "04_Sensitivity", firstRow  = TRUE)
    openxlsx::setColWidths(
      wb, "04_Sensitivity",
      cols   = seq_len(ncol(sensitivity_table)),
      widths = "auto"
    )
  } else {
    openxlsx::writeData(
      wb, "04_Sensitivity",
      data.frame(
        Note  = paste0(
          "Sensitivity analyses have not yet been run. Use ",
          "run_formal_full_analysis() or run_formal_sensitivity()."
        ),
        stringsAsFactors  = FALSE
      )
    )
  }

  openxlsx::saveWorkbook(wb, path, overwrite  = TRUE)
  invisible(path)
}


add_sensitivity_metadata <- function(
  main_table,
  group,
  variant_id,
  setting,
  reference        = FALSE,
  online_borrowing = data.frame()
) {

  if (!nrow(main_table)) return(data.frame())

  out <- main_table
  out <- out[, setdiff(names(out), "Scenario runtime (min)"), drop  = FALSE]

  out <- cbind(
    data.frame(
      Group            = group,
      Variant          = variant_id,
      Setting          = setting,
      Reference        = if (reference) "Yes" else "No",
      stringsAsFactors = FALSE
    ),
    out
  )

  out$`Mean applied delta | evaluated` <- NA_real_
  out$`Mean effective adult ESS | evaluated` <- NA_real_
  out$`Cap triggered (%) | evaluated` <- NA_real_

  if (nrow(online_borrowing)) {
    key_b <- paste(online_borrowing$Scenario, "PA-Keyboard", sep  = "||")
    key_o <- paste(out$Scenario, out$Method, sep  = "||")
    idx   <- match(key_o, key_b)
    hit   <- !is.na(idx)
    out$`Mean applied delta | evaluated`[hit] <-
    online_borrowing$`Mean applied delta | evaluated`[idx[hit]]
    out$`Mean effective adult ESS | evaluated`[hit] <-
    online_borrowing$`Mean effective adult ESS | evaluated`[idx[hit]]
    out$`Cap triggered (%) | evaluated`[hit] <-
    online_borrowing$`Cap triggered (%) | evaluated`[idx[hit]]
  }

  out
}

make_sensitivity_reference_rows <- function(main_results) {
  main_table <- main_results$main
  online     <- main_results$online_borrowing
  if (!nrow(main_table)) return(data.frame())

  rows <- list()

  add_reference <- function(group, setting, methods) {
    z <- main_table[
      main_table$Scenario %in% FORMAL_SENSITIVITY_SCENARIOS &
      main_table$Method %in% methods,
      ,
      drop  = FALSE
    ]
    add_sensitivity_metadata(
      z,
      group            = group,
      variant_id       = "MAIN",
      setting          = setting,
      reference        = TRUE,
      online_borrowing = online
    )
  }

  rows[["AUC caliper"]] <- add_reference(
    "AUC caliper", "0.22 (MAIN)", "PA-Keyboard"
  )
  rows[["Adult information"]] <- add_reference(
    "Adult information", "Lead 6 (MAIN)", c("PA-Keyboard", "PA-CRM")
  )
  rows[["Adult ESS cap"]] <- add_reference(
    "Adult ESS cap", "6 (MAIN)", "PA-Keyboard"
  )
  rows[["Borrowing gate"]] <- add_reference(
    "Borrowing gate", "nP >= 6 (MAIN)", "PA-Keyboard"
  )

  rbind_fill(rows)
}

sort_sensitivity_rows <- function(df) {
  if (!nrow(df)) return(df)

  group_order <- c(
    "AUC caliper",
    "Adult information",
    "Adult ESS cap",
    "Borrowing gate"
  )

  setting_order <- c(
    "0.16", "0.22 (MAIN)", "0.28",
    "Lead 6 (MAIN)", "Lead 9", "Lead 12", "Adult completed",
    "3", "6 (MAIN)", "9",
    "nP >= 3", "nP >= 6 (MAIN)"
  )

  df[
    order(
      match(df$Group, group_order),
      match(df$Setting, setting_order),
      match(df$Scenario, FORMAL_SENSITIVITY_SCENARIOS),
      match(df$Method, FORMAL_METHODS)
    ),
    ,
    drop  = FALSE
  ]
}


# =============================================================================
# Formatted manuscript and supplementary tables
# =============================================================================
create_manuscript_tables <- function(analysis, output_file) {
  required_pkgs <- c("dplyr", "tidyr", "openxlsx")
  missing_pkgs  <- required_pkgs[
    !vapply(required_pkgs, requireNamespace, logical(1), quietly  = TRUE)
  ]
  if (length(missing_pkgs) > 0L) {
    stop(
      "Required packages are unavailable: ",
      paste(missing_pkgs, collapse  = ", "),
      call.  = FALSE
    )
  }

  suppressPackageStartupMessages({
      library(dplyr)
      library(tidyr)
      library(openxlsx)
  })
  # ------------------------------------------------------------
  # 1. Manuscript constants and display order
  # ------------------------------------------------------------

  scenario_order <- c(
    "E0", "E1", "E2", "E3",
    "M1", "M2",
    "L1", "L2",
    "S1", "S2"
  )

  main_scenarios <- c(
    "E0", "E1", "E2", "E3",
    "M1", "M2",
    "L1", "L2"
  )

  safety_scenarios      <- c("S1", "S2")
  sensitivity_scenarios <- c("E1", "E3", "L2", "S2")

  method_order           <- c("Keyboard", "PA-Keyboard", "PA-CRM")
  borrowing_method_order <- c("PA-Keyboard", "PA-CRM")
  dose_order             <- paste0("P", 1:5)

  # Adult nominal doses are part of the manuscript trial specification
  # but are not stored as a standalone vector in the formal result workbook.
  adult_nominal_doses_mg <- c(60, 85, 110, 150, 180)

  # Display wording used in Supplementary Table B-1.
  purpose_map <- c(
    E0 = "Exact exposure and toxicity alignment",
    E1 = "Systematic exposure shift",
    E2 = "Larger exposure shift",
    E3 = "Irregular and many-to-one exposure mapping",
    M1 = "Pediatric MTD at the lowest dose",
    M2 = "Pediatric MTD at the highest dose",
    L1 = "Local incompatibility below the target dose",
    L2 = "Local incompatibility at the target dose",
    S1 = "No acceptable pediatric dose; compatible matched adult information",
    S2 = "No acceptable pediatric dose; lower adult toxicity at matched exposures"
  )


  # ------------------------------------------------------------
  # 2. Helper functions
  # ------------------------------------------------------------

  as_num <- function(x) {
    suppressWarnings(as.numeric(x))
  }

  get_design_value <- function(design_df, parameter) {
    out <- design_df$Value[design_df$Parameter == parameter]
    if (length(out) == 0L) {
      stop("Design parameter not found: ", parameter)
    }
    out[1]
  }

  fmt_pct <- function(x) {
    x <- as_num(x)
    ifelse(is.na(x), "\u2014", sprintf("%.1f", x))
  }

  fmt_ess <- function(x) {
    x <- as_num(x)
    ifelse(
      is.na(x),
      "\u2014",
      ifelse(x > 0 & x < 0.01, "<0.01", sprintf("%.2f", x))
    )
  }

  fmt_vector <- function(x, digits  = 2L, star_index  = NA_integer_) {
    vals <- as_num(strsplit(as.character(x), ",", fixed  = TRUE)[[1]])
    out  <- sprintf(paste0("%.", digits, "f"), vals)

    if (!is.na(star_index) &&
      star_index >= 1L &&
      star_index <= length(out)) {
      out[star_index] <- paste0(out[star_index], "*")
    }

    paste0("(", paste(out, collapse  = ", "), ")")
  }

  # Use the in-memory analysis objects produced by the simulation.
  if (is.null(analysis$main) || is.null(analysis$sensitivity)) {
    stop("analysis must be returned by run_formal_full_analysis().", call.  = FALSE)
  }

  config              <- analysis$main$config
  design              <- formal_config_table(config)
  scenario_truth      <- formal_scenario_table()
  main                <- analysis$main$main
  dose                <- analysis$main$dose
  matching_raw        <- analysis$main$matching
  borrowing_final_raw <- analysis$main$final_borrowing
  sensitivity         <- analysis$sensitivity$table

  dir.create(dirname(output_file), recursive  = TRUE, showWarnings  = FALSE)
  # ------------------------------------------------------------
  # 4. Main manuscript Table 1
  #    Basic settings for the simulation study
  # ------------------------------------------------------------

  pediatric_nominal_doses <- dose %>%
  filter(Scenario == "E0", Method == "Keyboard") %>%
  arrange(factor(Dose, levels  = dose_order)) %>%
  distinct(Dose, `Dose (mg/m2)`) %>%
  pull(`Dose (mg/m2)`)

  table1 <- tibble(
    Setting  = c(
      "Pediatric doses",
      "Adult doses",
      "Target DLT probability",
      "Target key",
      "Cohort size",
      "Maximum pediatric sample size",
      "Maximum adult sample size",
      "Adult lead-in"
    ),
    Value  = c(
      paste0(
        paste(as.integer(pediatric_nominal_doses), collapse  = ", "),
        " mg/m^2"
      ),
      paste0(
        paste(adult_nominal_doses_mg, collapse  = ", "),
        " mg"
      ),
      sprintf("%.2f", as_num(get_design_value(design, "Target_DLT"))),
      as.character(get_design_value(design, "Target_key")),
      as.character(as.integer(as_num(get_design_value(design, "Cohort_size")))),
      as.character(as.integer(as_num(get_design_value(design, "Pediatric_max_N")))),
      as.character(as.integer(as_num(get_design_value(design, "Adult_max_N")))),
      paste0(
        as.integer(as_num(get_design_value(design, "Adult_lead_N"))),
        " patients"
      )
    )
  )


  # ------------------------------------------------------------
  # 5. Main manuscript Table 2
  #    PCS, MTD allocation, overtoxic allocation
  # ------------------------------------------------------------

  table2_long <- bind_rows(
    main %>%
    filter(Scenario %in% main_scenarios) %>%
    transmute(
      `Operating characteristic`  = "PCS (%)",
      Method,
      Scenario,
      Value  = as_num(`PCS (%)`)
    ),

    main %>%
    filter(Scenario %in% main_scenarios) %>%
    transmute(
      `Operating characteristic`  = "MTD allocation (%)",
      Method,
      Scenario,
      Value  = as_num(`Patients at true MTD (%)`)
    ),

    main %>%
    filter(Scenario %in% main_scenarios) %>%
    transmute(
      `Operating characteristic`  = "Overtoxic allocation (%)",
      Method,
      Scenario,
      Value  = as_num(`Patients at overly toxic doses (%)`)
    )
  )

  table2 <- table2_long %>%
  mutate(
    `Operating characteristic`  = factor(
      `Operating characteristic`,
      levels  = c(
        "PCS (%)",
        "MTD allocation (%)",
        "Overtoxic allocation (%)"
      )
    ),
    Method   = factor(Method, levels  = method_order),
    Scenario = factor(Scenario, levels  = main_scenarios)
  ) %>%
  arrange(`Operating characteristic`, Method, Scenario) %>%
  pivot_wider(
    names_from  = Scenario,
    values_from = Value
  ) %>%
  arrange(`Operating characteristic`, Method) %>%
  mutate(
    across(
      all_of(main_scenarios),
      ~ round(as_num(.x), 1)
    ),
    Method  = as.character(Method),
    `Operating characteristic`  = 
    as.character(`Operating characteristic`)
  ) %>%
  select(
    `Operating characteristic`,
    Method,
    all_of(main_scenarios)
  )


  # ------------------------------------------------------------
  # 6. Main manuscript Table 3
  #    Safety operating characteristics in S1 and S2
  # ------------------------------------------------------------

  p1_allocation <- dose %>%
  filter(
    Scenario %in% safety_scenarios,
    Dose == "P1"
  ) %>%
  transmute(
    Scenario,
    Method,
    P1_allocation  = as_num(`Allocation (%)`)
  )

  table3_long <- bind_rows(
    main %>%
    filter(Scenario %in% safety_scenarios) %>%
    transmute(
      `Operating characteristic`  = "No-MTD selection (%)",
      Scenario,
      Method,
      Value  = as_num(`Correct no-MTD (%)`)
    ),

    main %>%
    filter(Scenario %in% safety_scenarios) %>%
    transmute(
      `Operating characteristic`  = "Early stopping (%)",
      Scenario,
      Method,
      Value  = as_num(`Early stop (%)`)
    ),

    p1_allocation %>%
    transmute(
      `Operating characteristic`  = "Patients at P1 (%)",
      Scenario,
      Method,
      Value  = P1_allocation
    ),

    p1_allocation %>%
    transmute(
      `Operating characteristic`  = "Patients above P1 (%)",
      Scenario,
      Method,
      Value  = 100 - P1_allocation
    )
  )

  table3 <- table3_long %>%
  mutate(
    `Operating characteristic`  = factor(
      `Operating characteristic`,
      levels  = c(
        "No-MTD selection (%)",
        "Early stopping (%)",
        "Patients at P1 (%)",
        "Patients above P1 (%)"
      )
    ),
    Method   = factor(Method, levels  = method_order),
    Scenario = factor(Scenario, levels  = safety_scenarios)
  ) %>%
  arrange(`Operating characteristic`, Method, Scenario) %>%
  pivot_wider(
    names_from  = Scenario,
    values_from = Value
  ) %>%
  arrange(`Operating characteristic`, Method) %>%
  mutate(
    across(
      all_of(safety_scenarios),
      ~ round(as_num(.x), 1)
    ),
    Method  = as.character(Method),
    `Operating characteristic`  = 
    as.character(`Operating characteristic`)
  ) %>%
  select(
    `Operating characteristic`,
    Method,
    all_of(safety_scenarios)
  )


  # ------------------------------------------------------------
  # 7. Main manuscript Table 4
  #    Adult-information availability sensitivity
  # ------------------------------------------------------------

  adult_info <- sensitivity %>%
  filter(
    Group == "Adult information",
    Scenario %in% sensitivity_scenarios,
    Method %in% borrowing_method_order
  )

  table4_long <- bind_rows(
    adult_info %>%
    filter(Scenario %in% c("E1", "E3", "L2")) %>%
    transmute(
      Scenario,
      `Operating characteristic`  = "PCS (%)",
      Method,
      Setting,
      Value  = as_num(`PCS (%)`)
    ),

    adult_info %>%
    filter(Scenario %in% c("E1", "E3", "L2")) %>%
    transmute(
      Scenario,
      `Operating characteristic`  = "Overtoxic allocation (%)",
      Method,
      Setting,
      Value  = as_num(`Patients at overly toxic doses (%)`)
    ),

    adult_info %>%
    filter(Scenario == "S2") %>%
    transmute(
      Scenario,
      `Operating characteristic`  = "Early stopping (%)",
      Method,
      Setting,
      Value  = as_num(`Early stop (%)`)
    )
  )

  table4 <- table4_long %>%
  mutate(
    Scenario  = factor(
      Scenario,
      levels  = c("E1", "E3", "L2", "S2")
    ),
    `Operating characteristic`  = factor(
      `Operating characteristic`,
      levels  = c(
        "PCS (%)",
        "Overtoxic allocation (%)",
        "Early stopping (%)"
      )
    ),
    Method  = factor(Method, levels  = borrowing_method_order)
  ) %>%
  arrange(
    Scenario,
    `Operating characteristic`,
    Method
  ) %>%
  pivot_wider(
    names_from  = Setting,
    values_from = Value
  ) %>%
  mutate(
    across(
      any_of(c(
          "Lead 6 (MAIN)",
          "Lead 9",
          "Lead 12",
          "Adult completed"
      )),
      ~ round(as_num(.x), 1)
    ),
    Scenario = as.character(Scenario),
    Method   = as.character(Method),
    `Operating characteristic`  = 
    as.character(`Operating characteristic`)
  ) %>%
  select(
    Scenario,
    `Operating characteristic`,
    Method,
    `Lead 6 (MAIN)`,
    `Lead 9`,
    `Lead 12`,
    `Adult completed`
  )


  # ------------------------------------------------------------
  # 8. Supplementary Table B-1
  #    Simulation scenarios
  # ------------------------------------------------------------

  table_b1 <- scenario_truth %>%
  mutate(
    scenario  = factor(scenario, levels  = scenario_order)
  ) %>%
  arrange(scenario) %>%
  rowwise() %>%
  mutate(
    Purpose  = unname(purpose_map[as.character(scenario)]),
    `DLT probabilities`  = paste0(
      "Pediatric: ",
      fmt_vector(
        pediatric_dlt,
        digits     = 2L,
        star_index = ifelse(
          is.na(true_mtd),
          NA_integer_,
          as.integer(true_mtd)
        )
      ),
      "\nAdult: ",
      fmt_vector(adult_dlt, digits  = 2L)
    ),
    AUC  = paste0(
      "Pediatric: ",
      fmt_vector(pediatric_auc, digits  = 2L),
      "\nAdult: ",
      fmt_vector(adult_auc, digits  = 2L)
    )
  ) %>%
  ungroup() %>%
  transmute(
    Scenario  = as.character(scenario),
    Purpose,
    `DLT probabilities`,
    AUC
  )


  # ------------------------------------------------------------
  # 9. Supplementary Table B-2
  #    Dose-specific operating characteristics
  # ------------------------------------------------------------

  dose_core <- dose %>%
  filter(
    Scenario %in% scenario_order,
    Method %in% method_order,
    Dose %in% dose_order
  ) %>%
  mutate(
    Scenario = factor(Scenario, levels  = scenario_order),
    Method   = factor(Method, levels  = method_order),
    Dose     = factor(Dose, levels  = dose_order)
  )

  b2_selection <- dose_core %>%
  transmute(
    Scenario,
    Method,
    Dose,
    Value  = as_num(`Selection (%)`)
  ) %>%
  pivot_wider(
    names_from   = Dose,
    values_from  = Value,
    names_prefix = "Sel_"
  )

  b2_no_mtd <- dose_core %>%
  group_by(Scenario, Method) %>%
  summarise(
    No_MTD  = max(
      0,
      100 - sum(as_num(`Selection (%)`), na.rm  = TRUE)
    ),
    .groups  = "drop"
  )

  b2_allocation <- dose_core %>%
  transmute(
    Scenario,
    Method,
    Dose,
    Value  = as_num(`Allocation (%)`)
  ) %>%
  pivot_wider(
    names_from   = Dose,
    values_from  = Value,
    names_prefix = "Alloc_"
  )

  table_b2 <- b2_selection %>%
  left_join(
    b2_no_mtd,
    by  = c("Scenario", "Method")
  ) %>%
  left_join(
    b2_allocation,
    by  = c("Scenario", "Method")
  ) %>%
  arrange(Scenario, Method) %>%
  mutate(
    across(
      matches("^(Sel_|Alloc_)|^No_MTD$"),
      ~ round(as_num(.x), 1)
    ),
    Scenario = as.character(Scenario),
    Method   = as.character(Method)
  ) %>%
  select(
    Scenario,
    Method,
    Sel_P1, Sel_P2, Sel_P3, Sel_P4, Sel_P5,
    No_MTD,
    Alloc_P1, Alloc_P2, Alloc_P3, Alloc_P4, Alloc_P5
  )


  # ------------------------------------------------------------
  # 10. Supplementary Table B-3
  #     Dose-specific realized borrowing
  # ------------------------------------------------------------

  table_b3 <- borrowing_final_raw %>%
  transmute(
    Scenario  = as.character(Scenario),
    `Pediatric dose`  = as.character(Dose),
    `Dose visited (%)`  = 
    fmt_pct(`Dose visited (%)`),
    `Final AUC match available (%)`  = 
    fmt_pct(`Final AUC match available (%) | visited`),
    `Borrowing activated (%)`  = 
    fmt_pct(`Borrow activated (%) | all trials`),
    `Mean realized adult ESS`  = 
    fmt_ess(`Mean realized adult ESS | all trials`)
  ) %>%
  mutate(
    Scenario  = factor(Scenario, levels  = scenario_order),
    `Pediatric dose`  = factor(
      `Pediatric dose`,
      levels  = dose_order
    )
  ) %>%
  arrange(Scenario, `Pediatric dose`) %>%
  mutate(
    Scenario  = as.character(Scenario),
    `Pediatric dose`  = as.character(`Pediatric dose`)
  )


  # ------------------------------------------------------------
  # 11. Supplementary Table B-4
  #     Dose-specific AUC matching performance
  # ------------------------------------------------------------

  matching_core <- matching_raw %>%
  transmute(
    Scenario                     = as.character(Scenario),
    Dose                         = as.character(Dose),
    Match_available_given_oracle = 
    as_num(`Sample match available (%) | oracle available`),
    Correct_match_given_availability  = 
    as_num(`Correct AUC match (%) | oracle and sample match available`)
  )

  matching_grid <- expand_grid(
    Scenario = scenario_order,
    Dose     = dose_order
  ) %>%
  left_join(
    matching_core,
    by  = c("Scenario", "Dose")
  )

  b4_match_available <- matching_grid %>%
  transmute(
    `Matching measure`  = "Match available, given oracle (%)",
    Scenario,
    Dose,
    Value  = fmt_pct(Match_available_given_oracle)
  )

  b4_correct_match <- matching_grid %>%
  transmute(
    `Matching measure`  = "Correct match, given availability (%)",
    Scenario,
    Dose,
    Value  = fmt_pct(Correct_match_given_availability)
  )

  table_b4 <- bind_rows(
    b4_match_available,
    b4_correct_match
  ) %>%
  mutate(
    `Matching measure`  = factor(
      `Matching measure`,
      levels  = c(
        "Match available, given oracle (%)",
        "Correct match, given availability (%)"
      )
    ),
    Scenario  = factor(
      Scenario,
      levels  = scenario_order
    ),
    Dose  = factor(
      Dose,
      levels  = dose_order
    )
  ) %>%
  arrange(
    `Matching measure`,
    Scenario,
    Dose
  ) %>%
  pivot_wider(
    names_from  = Dose,
    values_from = Value
  ) %>%
  mutate(
    `Matching measure`  = 
    as.character(`Matching measure`),
    Scenario  = as.character(Scenario)
  ) %>%
  select(
    `Matching measure`,
    Scenario,
    P1, P2, P3, P4, P5
  )


  # ------------------------------------------------------------
  # 12. Supplementary Table B-5
  #     Full sensitivity analyses
  # ------------------------------------------------------------

  b5_long <- bind_rows(
    sensitivity %>%
    filter(Scenario %in% c("E1", "E3", "L2")) %>%
    transmute(
      Group,
      Setting,
      Method,
      Scenario,
      Metric = "PCS",
      Value  = as_num(`PCS (%)`)
    ),

    sensitivity %>%
    filter(Scenario %in% c("E1", "E3", "L2")) %>%
    transmute(
      Group,
      Setting,
      Method,
      Scenario,
      Metric = "Overtoxic",
      Value  = as_num(`Patients at overly toxic doses (%)`)
    ),

    sensitivity %>%
    filter(Scenario == "S2") %>%
    transmute(
      Group,
      Setting,
      Method,
      Scenario,
      Metric = "EarlyStop",
      Value  = as_num(`Early stop (%)`)
    )
  ) %>%
  mutate(
    Key  = paste(Scenario, Metric, sep  = "_")
  )

  table_b5 <- b5_long %>%
  select(
    Group,
    Setting,
    Method,
    Key,
    Value
  ) %>%
  distinct() %>%
  pivot_wider(
    names_from  = Key,
    values_from = Value
  ) %>%
  mutate(
    `Sensitivity factor`  = recode(
      Group,
      "AUC caliper"  = "AUC caliper",
      "Adult information"  = 
      "Adult data available before pediatric trial initiation",
      "Adult ESS cap"  = 
      "Maximum borrowed adult ESS",
      "Borrowing gate"  = 
      "Borrowing gate"
    ),
    group_order  = match(
      Group,
      c(
        "AUC caliper",
        "Adult information",
        "Adult ESS cap",
        "Borrowing gate"
      )
    ),
    setting_order  = case_when(
      Group == "AUC caliper" &
      Setting == "0.16" ~ 1,
      Group == "AUC caliper" &
      Setting == "0.22 (MAIN)" ~ 2,
      Group == "AUC caliper" &
      Setting == "0.28" ~ 3,

      Group == "Adult information" &
      Setting == "Lead 6 (MAIN)" ~ 1,
      Group == "Adult information" &
      Setting == "Lead 9" ~ 2,
      Group == "Adult information" &
      Setting == "Lead 12" ~ 3,
      Group == "Adult information" &
      Setting == "Adult completed" ~ 4,

      Group == "Adult ESS cap" &
      Setting == "3" ~ 1,
      Group == "Adult ESS cap" &
      Setting == "6 (MAIN)" ~ 2,
      Group == "Adult ESS cap" &
      Setting == "9" ~ 3,

      Group == "Borrowing gate" &
      Setting == "nP >= 3" ~ 1,
      Group == "Borrowing gate" &
      Setting == "nP >= 6 (MAIN)" ~ 2,

      TRUE ~ 99
    ),
    method_order_num  = match(
      Method,
      borrowing_method_order
    )
  ) %>%
  arrange(
    group_order,
    setting_order,
    method_order_num
  ) %>%
  transmute(
    `Sensitivity factor`,
    Setting,
    Method,
    E1_PCS       = round(as_num(E1_PCS), 1),
    E1_Overtoxic = round(as_num(E1_Overtoxic), 1),
    E3_PCS       = round(as_num(E3_PCS), 1),
    E3_Overtoxic = round(as_num(E3_Overtoxic), 1),
    L2_PCS       = round(as_num(L2_PCS), 1),
    L2_Overtoxic = round(as_num(L2_Overtoxic), 1),
    S2_EarlyStop = round(as_num(S2_EarlyStop), 1)
  )


  # ------------------------------------------------------------
  # 13. Workbook styles
  # ------------------------------------------------------------

  wb <- createWorkbook()

  title_style <- createStyle(
    fontSize       = 12,
    textDecoration = "bold",
    halign         = "left",
    valign         = "center"
  )

  header_style <- createStyle(
    textDecoration = "bold",
    fgFill         = "#E9EDF2",
    halign         = "center",
    valign         = "center",
    wrapText       = TRUE,
    border         = "TopBottomLeftRight",
    borderColour   = "#A6A6A6",
    borderStyle    = "thin"
  )

  subheader_style <- createStyle(
    textDecoration = "bold",
    fgFill         = "#F3F5F7",
    halign         = "center",
    valign         = "center",
    wrapText       = TRUE,
    border         = "TopBottomLeftRight",
    borderColour   = "#BFBFBF",
    borderStyle    = "thin"
  )

  body_style <- createStyle(
    valign       = "center",
    halign       = "center",
    border       = "Bottom",
    borderColour = "#E0E0E0",
    borderStyle  = "hair"
  )

  body_left_style <- createStyle(
    valign       = "center",
    halign       = "left",
    wrapText     = TRUE,
    border       = "Bottom",
    borderColour = "#E0E0E0",
    borderStyle  = "hair"
  )

  note_style <- createStyle(
    fontSize       = 9,
    textDecoration = "italic",
    halign         = "left",
    valign         = "top",
    wrapText       = TRUE
  )


  # ------------------------------------------------------------
  # 14. Generic writing helpers
  # ------------------------------------------------------------

  write_simple_table <- function(
    wb,
    sheet,
    title,
    df,
    note      = NULL,
    widths    = NULL,
    left_cols = NULL
  ) {
    addWorksheet(wb, sheet, gridLines  = FALSE)

    n_cols <- ncol(df)

    writeData(
      wb,
      sheet,
      title,
      startRow = 1,
      startCol = 1,
      colNames = FALSE
    )
    mergeCells(
      wb,
      sheet,
      cols = 1:n_cols,
      rows = 1
    )
    addStyle(
      wb,
      sheet,
      title_style,
      rows       = 1,
      cols       = 1:n_cols,
      gridExpand = TRUE
    )
    setRowHeights(wb, sheet, rows  = 1, heights  = 22)

    writeData(
      wb,
      sheet,
      df,
      startRow    = 2,
      startCol    = 1,
      colNames    = TRUE,
      rowNames    = FALSE,
      headerStyle = header_style,
      keepNA      = FALSE
    )

    if (nrow(df) > 0L) {
      data_rows <- 3:(nrow(df) + 2L)

      addStyle(
        wb,
        sheet,
        body_style,
        rows       = data_rows,
        cols       = 1:n_cols,
        gridExpand = TRUE,
        stack      = TRUE
      )

      if (!is.null(left_cols) && length(left_cols) > 0L) {
        addStyle(
          wb,
          sheet,
          body_left_style,
          rows       = data_rows,
          cols       = left_cols,
          gridExpand = TRUE,
          stack      = TRUE
        )
      }
    }

    if (!is.null(note) && nzchar(note)) {
      note_row <- nrow(df) + 4L

      writeData(
        wb,
        sheet,
        note,
        startRow = note_row,
        startCol = 1,
        colNames = FALSE
      )
      mergeCells(
        wb,
        sheet,
        cols = 1:n_cols,
        rows = note_row
      )
      addStyle(
        wb,
        sheet,
        note_style,
        rows       = note_row,
        cols       = 1:n_cols,
        gridExpand = TRUE
      )
      setRowHeights(
        wb,
        sheet,
        rows    = note_row,
        heights = 45
      )
    }

    freezePane(
      wb,
      sheet,
      firstActiveRow  = 3
    )

    if (is.null(widths)) {
      setColWidths(
        wb,
        sheet,
        cols   = 1:n_cols,
        widths = "auto"
      )
    } else {
      setColWidths(
        wb,
        sheet,
        cols   = 1:n_cols,
        widths = widths
      )
    }
  }


  write_table_b2 <- function(wb, df) {
    sheet <- "Supp_Table_B2"
    title <- "Table B-2. Dose-specific operating characteristics"

    addWorksheet(wb, sheet, gridLines  = FALSE)

    writeData(
      wb, sheet, title,
      startRow = 1, startCol  = 1,
      colNames = FALSE
    )
    mergeCells(wb, sheet, cols  = 1:13, rows  = 1)
    addStyle(
      wb, sheet, title_style,
      rows       = 1, cols  = 1:13,
      gridExpand = TRUE
    )

    # Row 2: grouped headers
    writeData(wb, sheet, "Scenario", startRow  = 2, startCol  = 1, colNames  = FALSE)
    writeData(wb, sheet, "Method", startRow  = 2, startCol  = 2, colNames  = FALSE)
    writeData(wb, sheet, "Selection (%)", startRow  = 2, startCol  = 3, colNames  = FALSE)
    writeData(wb, sheet, "Allocation (%)", startRow  = 2, startCol  = 9, colNames  = FALSE)

    mergeCells(wb, sheet, cols  = 1, rows  = 2:3)
    mergeCells(wb, sheet, cols  = 2, rows  = 2:3)
    mergeCells(wb, sheet, cols  = 3:8, rows  = 2)
    mergeCells(wb, sheet, cols  = 9:13, rows  = 2)

    # Row 3: subheaders
    writeData(
      wb,
      sheet,
      t(c(
          "P1", "P2", "P3", "P4", "P5", "No MTD",
          "P1", "P2", "P3", "P4", "P5"
      )),
      startRow = 3,
      startCol = 3,
      colNames = FALSE
    )

    addStyle(
      wb, sheet, header_style,
      rows       = 2, cols  = 1:13,
      gridExpand = TRUE
    )
    addStyle(
      wb, sheet, subheader_style,
      rows       = 3, cols  = 1:13,
      gridExpand = TRUE,
      stack      = TRUE
    )

    writeData(
      wb,
      sheet,
      df,
      startRow = 4,
      startCol = 1,
      colNames = FALSE,
      rowNames = FALSE,
      keepNA   = FALSE
    )

    if (nrow(df) > 0L) {
      data_rows <- 4:(nrow(df) + 3L)

      addStyle(
        wb, sheet, body_style,
        rows       = data_rows,
        cols       = 1:13,
        gridExpand = TRUE
      )
      addStyle(
        wb, sheet, body_left_style,
        rows       = data_rows,
        cols       = 2,
        gridExpand = TRUE,
        stack      = TRUE
      )
    }

    note <- paste0(
      "Note: Values are percentages and are rounded to one decimal place. ",
      "Selection (%) gives the probability that each pediatric dose was selected as the MTD; ",
      "No MTD indicates that no MTD was selected. Allocation (%) gives the mean proportion ",
      "of pediatric patients treated at each dose. Bold values in the selection block indicate ",
      "the correct selection outcome (the true pediatric MTD for E0-E3, M1-M2, and L1-L2; ",
      "No MTD for S1-S2). Bold values in the allocation block indicate the true pediatric MTD."
    )

    note_row <- nrow(df) + 5L
    writeData(
      wb, sheet, note,
      startRow = note_row,
      startCol = 1,
      colNames = FALSE
    )
    mergeCells(
      wb, sheet,
      cols = 1:13,
      rows = note_row
    )
    addStyle(
      wb, sheet, note_style,
      rows       = note_row,
      cols       = 1:13,
      gridExpand = TRUE
    )
    setRowHeights(
      wb, sheet,
      rows    = note_row,
      heights = 60
    )

    # Bold correct-selection and true-MTD allocation cells.
    true_mtd_map <- c(
      E0 = 3,
      E1 = 3,
      E2 = 3,
      E3 = 3,
      M1 = 1,
      M2 = 5,
      L1 = 3,
      L2 = 3
    )

    for (r in seq_len(nrow(df))) {
      excel_row <- r + 3L
      sc        <- df$Scenario[r]

      if (sc %in% names(true_mtd_map)) {
        mtd <- true_mtd_map[[sc]]

        # Selection columns C:G correspond to P1:P5.
        sel_col <- 2L + mtd

        # Allocation columns I:M correspond to P1:P5.
        alloc_col <- 8L + mtd

        addStyle(
          wb, sheet,
          createStyle(textDecoration  = "bold"),
          rows       = excel_row,
          cols       = c(sel_col, alloc_col),
          gridExpand = TRUE,
          stack      = TRUE
        )
      } else if (sc %in% c("S1", "S2")) {
        # No-MTD is column H.
        addStyle(
          wb, sheet,
          createStyle(textDecoration  = "bold"),
          rows  = excel_row,
          cols  = 8,
          stack = TRUE
        )
      }
    }

    setColWidths(
      wb,
      sheet,
      cols   = 1:13,
      widths = c(
        9, 15,
        8, 8, 8, 8, 8, 10,
        8, 8, 8, 8, 8
      )
    )

    freezePane(
      wb,
      sheet,
      firstActiveRow  = 4
    )
  }


  write_table_b5 <- function(wb, df) {
    sheet <- "Supp_Table_B5"
    title <- "Table B-5. Results of the sensitivity analyses"

    addWorksheet(wb, sheet, gridLines  = FALSE)

    writeData(
      wb, sheet, title,
      startRow = 1, startCol  = 1,
      colNames = FALSE
    )
    mergeCells(wb, sheet, cols  = 1:10, rows  = 1)
    addStyle(
      wb, sheet, title_style,
      rows       = 1, cols  = 1:10,
      gridExpand = TRUE
    )

    # Top-level headers
    writeData(wb, sheet, "Sensitivity factor", startRow  = 2, startCol  = 1, colNames  = FALSE)
    writeData(wb, sheet, "Setting", startRow  = 2, startCol  = 2, colNames  = FALSE)
    writeData(wb, sheet, "Method", startRow  = 2, startCol  = 3, colNames  = FALSE)

    writeData(wb, sheet, "E1", startRow  = 2, startCol  = 4, colNames  = FALSE)
    writeData(wb, sheet, "E3", startRow  = 2, startCol  = 6, colNames  = FALSE)
    writeData(wb, sheet, "L2", startRow  = 2, startCol  = 8, colNames  = FALSE)
    writeData(wb, sheet, "S2", startRow  = 2, startCol  = 10, colNames  = FALSE)

    mergeCells(wb, sheet, cols  = 1, rows  = 2:3)
    mergeCells(wb, sheet, cols  = 2, rows  = 2:3)
    mergeCells(wb, sheet, cols  = 3, rows  = 2:3)
    mergeCells(wb, sheet, cols  = 4:5, rows  = 2)
    mergeCells(wb, sheet, cols  = 6:7, rows  = 2)
    mergeCells(wb, sheet, cols  = 8:9, rows  = 2)

    writeData(
      wb,
      sheet,
      t(c(
          "PCS (%)", "Overtoxic allocation (%)",
          "PCS (%)", "Overtoxic allocation (%)",
          "PCS (%)", "Overtoxic allocation (%)",
          "Early stopping (%)"
      )),
      startRow = 3,
      startCol = 4,
      colNames = FALSE
    )

    addStyle(
      wb, sheet, header_style,
      rows       = 2, cols  = 1:10,
      gridExpand = TRUE
    )
    addStyle(
      wb, sheet, subheader_style,
      rows       = 3, cols  = 1:10,
      gridExpand = TRUE,
      stack      = TRUE
    )

    writeData(
      wb,
      sheet,
      df,
      startRow = 4,
      startCol = 1,
      colNames = FALSE,
      rowNames = FALSE,
      keepNA   = FALSE
    )

    if (nrow(df) > 0L) {
      data_rows <- 4:(nrow(df) + 3L)

      addStyle(
        wb, sheet, body_style,
        rows       = data_rows,
        cols       = 1:10,
        gridExpand = TRUE
      )

      addStyle(
        wb, sheet, body_left_style,
        rows       = data_rows,
        cols       = 1:3,
        gridExpand = TRUE,
        stack      = TRUE
      )
    }

    note <- paste0(
      "Note: All values are percentages and are rounded to one decimal place. ",
      "Each design parameter was varied separately, while all other parameters were kept ",
      "at their primary-analysis values. For adult data availability, PA-Keyboard and PA-CRM ",
      "were both evaluated because this factor affects the information available to both borrowing ",
      "designs; the remaining sensitivity factors are specific to PA-Keyboard. ",
      "\"Lead 6\", \"Lead 9\", and \"Lead 12\" indicate the number of adult patients enrolled ",
      "before pediatric trial initiation. \"Adult completed\" indicates that the adult trial was ",
      "completed before pediatric enrollment began."
    )

    note_row <- nrow(df) + 5L
    writeData(
      wb, sheet, note,
      startRow = note_row,
      startCol = 1,
      colNames = FALSE
    )
    mergeCells(
      wb, sheet,
      cols = 1:10,
      rows = note_row
    )
    addStyle(
      wb, sheet, note_style,
      rows       = note_row,
      cols       = 1:10,
      gridExpand = TRUE
    )
    setRowHeights(
      wb, sheet,
      rows    = note_row,
      heights = 70
    )

    setColWidths(
      wb,
      sheet,
      cols   = 1:10,
      widths = c(
        32, 18, 15,
        10, 18,
        10, 18,
        10, 18,
        14
      )
    )

    freezePane(
      wb,
      sheet,
      firstActiveRow  = 4
    )
  }


  # ------------------------------------------------------------
  # 15. Write all nine tables
  # ------------------------------------------------------------

  write_simple_table(
    wb        = wb,
    sheet     = "Main_Table1",
    title     = "Table 1. Basic settings for the simulation study",
    df        = table1,
    widths    = c(32, 45),
    left_cols = 1:2
  )

  write_simple_table(
    wb    = wb,
    sheet = "Main_Table2",
    title = paste0(
      "Table 2. Operating characteristics of PA-Keyboard, Keyboard, ",
      "and PA-CRM across simulation scenarios"
    ),
    df        = table2,
    widths    = c(28, 15, rep(9, length(main_scenarios))),
    left_cols = c(1, 2)
  )

  write_simple_table(
    wb    = wb,
    sheet = "Main_Table3",
    title = paste0(
      "Table 3. Safety operating characteristics of PA-Keyboard, ",
      "Keyboard, and PA-CRM in scenarios S1 and S2"
    ),
    df        = table3,
    widths    = c(28, 15, 10, 10),
    left_cols = c(1, 2)
  )

  write_simple_table(
    wb    = wb,
    sheet = "Main_Table4",
    title = paste0(
      "Table 4. Sensitivity of PA-Keyboard and PA-CRM to the amount ",
      "of adult data available before pediatric trial initiation"
    ),
    df   = table4,
    note = paste0(
      "Note: Lead 6 is the primary-analysis setting. Adult completed indicates ",
      "that the adult trial was completed before pediatric enrollment began."
    ),
    widths    = c(10, 26, 15, 16, 12, 12, 16),
    left_cols = c(1, 2, 3)
  )

  write_simple_table(
    wb    = wb,
    sheet = "Supp_Table_B1",
    title = "Table B-1. Simulation scenarios",
    df    = table_b1,
    note  = paste0(
      "Note: Pediatric and adult values correspond to dose levels P1-P5 and A1-A5, ",
      "respectively. Values marked with an asterisk indicate the true pediatric MTD. ",
      "Scenarios S1 and S2 have no acceptable pediatric MTD."
    ),
    widths    = c(9, 36, 46, 46),
    left_cols = 1:4
  )

  # Increase row height for wrapped scenario cells.
  setRowHeights(
    wb,
    "Supp_Table_B1",
    rows    = 3:(nrow(table_b1) + 2L),
    heights = 42
  )

  write_table_b2(
    wb = wb,
    df = table_b2
  )

  write_simple_table(
    wb    = wb,
    sheet = "Supp_Table_B3",
    title = "Table B-3. Dose-specific realized borrowing of PA-Keyboard",
    df    = table_b3,
    note  = paste0(
      "Notes: Dose visited (%) is the percentage of successful simulation replications ",
      "in which at least one pediatric patient was treated at the corresponding dose. ",
      "Final AUC match available (%) is calculated among replications in which that ",
      "pediatric dose was visited. Borrowing activated (%) is calculated over all successful ",
      "simulation replications; borrowing is considered activated when the realized adult ESS ",
      "is greater than zero. Mean realized adult ESS is calculated over all successful simulation ",
      "replications. Replications in which the dose was not visited or no adult information was ",
      "borrowed contribute an ESS of zero. The maximum borrowed adult ESS is 6."
    ),
    widths    = c(9, 12, 14, 24, 20, 20),
    left_cols = c(1, 2)
  )

  write_simple_table(
    wb    = wb,
    sheet = "Supp_Table_B4",
    title = "Table B-4. Dose-specific AUC matching performance of PA-Keyboard",
    df    = table_b4,
    note  = paste0(
      "Notes: Match available, given oracle (%) is the percentage of online evaluations ",
      "in which the sample-based AUC matching procedure identified an acceptable adult match, ",
      "conditional on the population-level oracle dose being available. Correct match, given ",
      "availability (%) is the percentage of evaluations in which the sample-based matched adult ",
      "dose coincided with the oracle dose, conditional on both the oracle dose and a sample-based ",
      "AUC match being available. Evaluations with no sample-based match are excluded from the latter ",
      "denominator. \"\u2014\" indicates that the corresponding conditional measure was not defined ",
      "because no oracle matching opportunity was available."
    ),
    widths    = c(34, 10, rep(9, 5)),
    left_cols = c(1, 2)
  )

  write_table_b5(
    wb = wb,
    df = table_b5
  )


  # ------------------------------------------------------------
  # 16. Final workbook settings and export
  # ------------------------------------------------------------

  # Uniform tab colour is intentionally omitted to keep the workbook simple.
  # Activate the first worksheet when opening the file.
  activeSheet(wb) <- 1

  saveWorkbook(
    wb,
    output_file,
    overwrite  = TRUE
  )

  message("Manuscript tables saved to: ", normalizePath(output_file, winslash  = "/", mustWork  = FALSE))

  invisible(list(
      path   = output_file,
      tables = list(
        Main_Table1   = table1,
        Main_Table2   = table2,
        Main_Table3   = table3,
        Main_Table4   = table4,
        Supp_Table_B1 = table_b1,
        Supp_Table_B2 = table_b2,
        Supp_Table_B3 = table_b3,
        Supp_Table_B4 = table_b4,
        Supp_Table_B5 = table_b5
      )
  ))
}

# =============================================================================
# Publication analysis and output regeneration
# =============================================================================
run_publication_analysis <- function(
  output_dir = file.path(PAK_PROJECT_DIR, "results"),
  config     = NULL,
  n_sim      = NULL,
  n_workers  = NULL,
  resume     = TRUE
) {

  required_pkgs <- c("openxlsx", "dplyr", "tidyr", "ggplot2")
  missing_pkgs  <- required_pkgs[
    !vapply(required_pkgs, requireNamespace, logical(1), quietly  = TRUE)
  ]
  if (length(missing_pkgs) > 0L) {
    stop(
      "Required packages are unavailable: ",
      paste(missing_pkgs, collapse  = ", "),
      call.  = FALSE
    )
  }

  dir.create(output_dir, recursive  = TRUE, showWarnings  = FALSE)

  if (is.null(config)) {
    config <- make_formal_config(
      overrides  = list(write_excel  = TRUE, verbose  = TRUE),
      output_dir = output_dir
    )
  } else {
    config$output_dir <- output_dir
  }

  if (!is.null(n_sim)) config$n_sim <- as.integer(n_sim)
  if (!is.null(n_workers)) config$n_workers <- as.integer(n_workers)
  config    <- validate_formal_config(config)
  n_sim     <- config$n_sim
  n_workers <- config$n_workers

  output_paths <- list(
    raw_results       = file.path(output_dir, "PA_Keyboard_raw_results.xlsx"),
    manuscript_tables = file.path(output_dir, "PA_Keyboard_manuscript_tables.xlsx"),
    figure1           = file.path(output_dir, "Figure1_scenario_truth.tiff"),
    figure2           = file.path(output_dir, "Figure2_dose_finding_performance.tiff")
  )

  analysis <- run_formal_full_analysis(
    config                = config,
    main_scenarios        = FORMAL_SCENARIOS,
    sensitivity_scenarios = FORMAL_SENSITIVITY_SCENARIOS,
    n_sim                 = n_sim,
    n_workers             = n_workers,
    resume                = resume,
    output_xlsx           = output_paths$raw_results,
    script_path           = detect_formal_script_path()
  )

  manuscript <- create_manuscript_tables(
    analysis    = analysis,
    output_file = output_paths$manuscript_tables
  )
  figure1 <- create_figure1(output_paths$figure1)
  figure2 <- create_figure2(analysis, output_paths$figure2)

  writeLines(capture.output(utils::sessionInfo()),
    file.path(output_dir, "session_info.txt"))
  message("Publication analysis completed.")
  message("Output directory: ", normalizePath(output_dir, winslash  = "/", mustWork  = FALSE))

  invisible(list(
      analysis   = analysis,
      manuscript = manuscript,
      figures    = list(figure1  = figure1, figure2  = figure2),
      files      = output_paths
  ))
}


# Recreate the manuscript and supplementary tables and figures directly from
# a completed simulation's saved RDS summaries.
rebuild_publication_outputs <- function(
  output_dir   = file.path(PAK_PROJECT_DIR, "results"),
  analysis_rds = file.path(output_dir,
    "PA_keyboard_ETOVL_formal_full_compact_v03.rds")
) {
  if (!file.exists(analysis_rds)) {
    stop("Saved full analysis not found: ", analysis_rds, call.  = FALSE)
  }
  analysis <- readRDS(analysis_rds)
  if (is.null(analysis$main$config) || is.null(analysis$sensitivity$table)) {
    stop("RDS does not contain complete main/sensitivity results.", call.  = FALSE)
  }
  dir.create(output_dir, recursive  = TRUE, showWarnings  = FALSE)
  paths <- list(
    raw     = file.path(output_dir, "PA_Keyboard_raw_results.xlsx"),
    tables  = file.path(output_dir, "PA_Keyboard_manuscript_tables.xlsx"),
    figure1 = file.path(output_dir, "Figure1_scenario_truth.tiff"),
    figure2 = file.path(output_dir, "Figure2_dose_finding_performance.tiff")
  )
  main <- analysis$main
  write_formal_workbook(
    config            = main$config,
    main_table        = main$main,
    dose_table        = main$dose,
    matching_table    = main$matching,
    final_borrowing   = main$final_borrowing,
    online_borrowing  = main$online_borrowing,
    sensitivity_table = analysis$sensitivity$table,
    path              = paths$raw
  )
  create_manuscript_tables(analysis, paths$tables)
  create_figure1(paths$figure1)
  create_figure2(analysis, paths$figure2)
  writeLines(capture.output(utils::sessionInfo()),
    file.path(output_dir, "session_info.txt"))
  invisible(paths)
}

