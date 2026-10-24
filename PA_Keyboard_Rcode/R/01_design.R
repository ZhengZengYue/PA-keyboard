# =============================================================================
# PA-Keyboard simulation: configuration and study scenarios
# =============================================================================
# This module defines the simulation settings, prespecified pediatric and adult
# dose-toxicity/AUC profiles, and numerical utilities shared by the trial designs.
# Formal scenario IDs, method labels, and sensitivity settings are defined here.
#
# Function guide:
#   stopf(), assert_*()                Input validation and informative errors.
#   clamp(), logit(), inv_logit(),
#     log_sum_exp()                    Stable numerical transformations.
#   safe_mean(), safe_mcse(),
#     geometric_mean_from_logs()      Descriptive summaries and AUC aggregation.
#   linear_interpolate_extrapolate()  Map converted doses to CRM effective doses.
#   rbind_fill()                      Combine simulation outputs with shared columns.
#   simple_text_hash(),
#     config_signature()              Identify simulation configuration settings.
#   gauss_legendre_01(),
#     map_legendre_interval()         Numerical integration nodes and weights.
#   make_default_config(),
#     validate_config()               Define and validate trial/simulation parameters.
#   make_scenarios(),
#     scenario_truth_table()          Specify and tabulate simulation scenarios.
# =============================================================================

# Version identifiers used in configuration and PA-CRM cache signatures.
PAK_VERSION       <- "0.3.0"
PAK_FILE_VERSION  <- "v03"
PAK_PROTOCOL_DATE <- "2026-08-21"

# Numeric cache used by Keyboard keys and deterministic PA-CRM quadrature.
.PAK_NUMERIC_CACHE <- new.env(parent  = emptyenv())

# =============================================================================
# Simulation utilities, defaults and prespecified scenarios
# =============================================================================
stopf <- function(fmt, ...) {
  stop(sprintf(fmt, ...), call.  = FALSE)
}

assert_true <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call.  = FALSE)
  invisible(TRUE)
}

assert_scalar_numeric <- function(
  x, name, lower  = -Inf, upper  = Inf,
  lower_open  = FALSE, upper_open  = FALSE
) {
  ok <- is.numeric(x) && length(x) == 1L && is.finite(x)
  if (ok) {
    ok <- if (lower_open) x > lower else x >= lower
    ok <- ok && if (upper_open) x < upper else x <= upper
  }
  if (!ok) {
    stopf("%s must be one finite numeric value in the required range.", name)
  }
  invisible(TRUE)
}

assert_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stopf("%s must be TRUE or FALSE.", name)
  }
  invisible(TRUE)
}

clamp <- function(x, lower, upper) {
  pmin(pmax(x, lower), upper)
}

logit <- function(p) {
  p <- clamp(p, .Machine$double.eps, 1 - .Machine$double.eps)
  log(p) - log1p(-p)
}

inv_logit <- function(x) {
  # Start from x so matrix/array dimensions and names are preserved.
  out       <- x
  pos       <- x >= 0
  out[pos]  <- 1 / (1 + exp(-x[pos]))
  ex        <- exp(x[!pos])
  out[!pos] <- ex / (1 + ex)
  out
}

log_sum_exp <- function(x) {
  if (!length(x)) return(-Inf)
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}

safe_mean <- function(x) {
  if (!length(x) || all(is.na(x))) NA_real_ else mean(x, na.rm  = TRUE)
}

safe_mcse <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  stats::sd(x) / sqrt(length(x))
}

geometric_mean_from_logs <- function(log_sum, n) {
  ifelse(n > 0, exp(log_sum / n), NA_real_)
}

linear_interpolate_extrapolate <- function(x, y, xout) {
  if (length(x) != length(y) || length(x) < 2L || is.unsorted(x)) {
    stop("x and y must be ordered vectors of equal length >= 2.", call.  = FALSE)
  }
  ans   <- stats::approx(x, y, xout  = xout, rule  = 2, ties  = "ordered")$y
  below <- xout < x[1L]
  above <- xout > x[length(x)]
  if (any(below)) {
    slope      <- (y[2L] - y[1L]) / (x[2L] - x[1L])
    ans[below] <- y[1L] + slope * (xout[below] - x[1L])
  }
  if (any(above)) {
    k          <- length(x)
    slope      <- (y[k] - y[k - 1L]) / (x[k] - x[k - 1L])
    ans[above] <- y[k] + slope * (xout[above] - x[k])
  }
  ans
}

rbind_fill <- function(data_frames) {
  data_frames <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L, data_frames)
  if (!length(data_frames)) return(data.frame())
  all_names <- unique(unlist(lapply(data_frames, names), use.names  = FALSE))
  aligned   <- lapply(data_frames, function(d) {
      missing_names <- setdiff(all_names, names(d))
      for (nm in missing_names) d[[nm]] <- NA
      d[all_names]
  })
  out <- do.call(rbind, aligned)
  rownames(out) <- NULL
  out
}

simple_text_hash <- function(text) {
  bytes   <- as.integer(charToRaw(enc2utf8(text)))
  h       <- 2166136261
  modulus <- 4294967291
  for (b in bytes) h <- (h * 16777619 + b) %% modulus
  sprintf("%010.0f", h)
}

config_signature <- function(config) {
  x                        <- config
  x$output_dir             <- NULL
  x$n_workers              <- NULL
  x$enable_checkpoint      <- NULL
  x$resume_from_checkpoint <- NULL
  x$checkpoint_every       <- NULL
  text                     <- paste(capture.output(dput(x)), collapse  = "")
  simple_text_hash(text)
}

gauss_legendre_01 <- function(n  = 64L) {
  n   <- as.integer(n)
  key <- paste0("legendre01_", n)
  if (exists(key, envir  = .PAK_NUMERIC_CACHE, inherits  = FALSE)) {
    return(get(key, envir  = .PAK_NUMERIC_CACHE, inherits  = FALSE))
  }
  if (n < 2L) stop("Gauss-Legendre order must be at least 2.", call.  = FALSE)
  i                        <- seq_len(n - 1L)
  off                      <- i / sqrt(4 * i^2 - 1)
  jacobi                   <- matrix(0, nrow  = n, ncol  = n)
  jacobi[cbind(i, i + 1L)] <- off
  jacobi[cbind(i + 1L, i)] <- off
  eig                      <- eigen(jacobi, symmetric  = TRUE)
  ord                      <- order(eig$values)
  nodes                    <- (eig$values[ord] + 1) / 2
  weights                  <- eig$vectors[1L, ord]^2
  ans                      <- list(nodes  = nodes, weights  = weights)
  assign(key, ans, envir  = .PAK_NUMERIC_CACHE)
  ans
}

map_legendre_interval <- function(a, b, n  = 48L) {
  q <- gauss_legendre_01(n)
  list(
    nodes   = a + (b - a) * q$nodes,
    weights = (b - a) * q$weights
  )
}

make_default_config <- function(overrides  = list()) {
  config <- list(
    script_version                   = PAK_VERSION,
    file_version                     = PAK_FILE_VERSION,
    protocol_date                    = PAK_PROTOCOL_DATE,
    config_id                        = "MAIN",
    target                           = 0.30,
    target_lower                     = 0.25,
    target_upper                     = 0.35,
    cohort_size                      = 3L,
    adult_max_n                      = 30L,
    pediatric_max_n                  = 30L,
    adult_doses_mg                   = c(60, 85, 110, 150, 180),
    pediatric_doses_mg_m2            = c(35, 50, 65, 85, 105),
    bsa_conversion                   = 1.73,
    auc_sigma_log                    = 0.15,
    auc_caliper_log                  = 0.22,
    auc_tie_tolerance                = 1e-12,
    adult_lead_n                     = 6L,
    adult_reached_guard              = FALSE,
    safety_cutoff                    = 0.95,
    safety_min_n                     = 3L,
    mixture_w0                       = 0.50,
    adult_component_parameterization = "proportional_count_compression",
    elir_cap                         = 6,
    # Nodes per logit-space interval; 13 intervals are used for boundary accuracy.
    elir_quad_n            = 16L,
    elir_grid_n            = 41L,
    elir_root_tolerance    = 1e-6,
    elir_cap_tolerance     = 1e-6,
    enable_elir_cache      = TRUE,
    pa_crm_skeleton        = c(0.09, 0.17, 0.30, 0.45, 0.60),
    pa_crm_intercept       = 3,
    pa_crm_alpha_mean      = 0,
    pa_crm_alpha_variance  = 2,
    pa_crm_delta           = c(0, 0.5, 1),
    pa_crm_delta_prior     = rep(1 / 3, 3),
    pa_crm_safety_cutoff   = 0.90,
    pa_crm_quad_n_per_side = 80L,
    pa_crm_alpha_bounds    = c(-10, 10),
    n_sim                  = 5000L,
    n_workers              = 36L,
    master_seed            = 24082026L,
    enable_checkpoint      = FALSE,
    checkpoint_every       = 500L,
    resume_from_checkpoint = FALSE,
    write_excel            = TRUE,
    output_dir             = file.path(getwd(), "PA_keyboard_results_v03"),
    keep_interim_raw       = TRUE,
    fail_fast              = FALSE,
    verbose                = TRUE
  )
  if (length(overrides)) {
    unknown <- setdiff(names(overrides), names(config))
    if (length(unknown)) stopf("Unknown config field(s): %s", paste(unknown, collapse  = ", "))
    config[names(overrides)] <- overrides
  }
  validate_config(config)
}

validate_config <- function(config) {
  assert_scalar_numeric(config$target, "target", 0, 1, TRUE, TRUE)
  assert_scalar_numeric(config$target_lower, "target_lower", 0, config$target, FALSE, TRUE)
  assert_scalar_numeric(config$target_upper, "target_upper", config$target, 1, TRUE, FALSE)
  assert_true(abs((config$target - config$target_lower) -
      (config$target_upper - config$target)) < 1e-12,
    "The locked target key must be symmetric around target.")
  for (nm in c("cohort_size", "adult_max_n", "pediatric_max_n", "adult_lead_n",
      "n_sim", "n_workers", "checkpoint_every", "elir_quad_n",
      "elir_grid_n", "pa_crm_quad_n_per_side")) {
    value <- config[[nm]]
    if (!is.numeric(value) || length(value) != 1L || value < 1 || value != as.integer(value)) {
      stopf("%s must be a positive integer.", nm)
    }
  }
  assert_scalar_numeric(config$master_seed, "master_seed", 1, 2147483645)
  assert_true(config$master_seed == as.integer(config$master_seed),
    "master_seed must be an integer.")
  assert_true(config$adult_max_n %% config$cohort_size == 0L,
    "adult_max_n must be divisible by cohort_size.")
  assert_true(config$pediatric_max_n %% config$cohort_size == 0L,
    "pediatric_max_n must be divisible by cohort_size.")
  assert_true(config$adult_lead_n %% config$cohort_size == 0L,
    "adult_lead_n must be divisible by cohort_size.")
  assert_true(config$adult_lead_n < config$adult_max_n,
    "adult_lead_n must be smaller than adult_max_n.")
  assert_true(config$n_sim < 1000000L,
    "n_sim must be below 1,000,000 for the collision-free keyed RNG map.")
  assert_true(length(config$adult_doses_mg) == 5L &&
    all(diff(config$adult_doses_mg) > 0),
    "The locked adult dose vector must contain five increasing doses.")
  assert_true(length(config$pediatric_doses_mg_m2) == 5L &&
    all(diff(config$pediatric_doses_mg_m2) > 0),
    "The locked pediatric dose vector must contain five increasing doses.")
  assert_scalar_numeric(config$auc_sigma_log, "auc_sigma_log", 0, Inf, TRUE, FALSE)
  assert_scalar_numeric(config$auc_caliper_log, "auc_caliper_log", 0, Inf, TRUE, FALSE)
  assert_scalar_numeric(config$mixture_w0, "mixture_w0", 0, 1, TRUE, TRUE)
  assert_true(
    identical(config$adult_component_parameterization,
      "proportional_count_compression"),
    paste(
      "v03 requires adult_component_parameterization =",
      "'proportional_count_compression'."
    )
  )
  assert_scalar_numeric(config$elir_cap, "elir_cap", 0, Inf, FALSE, FALSE)
  assert_scalar_numeric(config$safety_cutoff, "safety_cutoff", 0, 1, TRUE, TRUE)
  assert_scalar_numeric(config$pa_crm_safety_cutoff, "pa_crm_safety_cutoff", 0, 1, TRUE, TRUE)
  assert_true(length(config$pa_crm_skeleton) == 5L &&
    all(diff(config$pa_crm_skeleton) > 0),
    "PA-CRM skeleton must contain five increasing probabilities.")
  assert_true(length(config$pa_crm_delta) == length(config$pa_crm_delta_prior),
    "PA-CRM delta values and prior probabilities must have equal length.")
  assert_true(abs(sum(config$pa_crm_delta_prior) - 1) < 1e-12,
    "PA-CRM delta prior probabilities must sum to one.")
  for (nm in c("adult_reached_guard", "enable_elir_cache", "enable_checkpoint",
      "resume_from_checkpoint", "write_excel", "keep_interim_raw",
      "fail_fast", "verbose")) assert_flag(config[[nm]], nm)
  assert_true(!config$adult_reached_guard,
    "adult_reached_guard is locked OFF in v03.")
  if (config$resume_from_checkpoint && !config$enable_checkpoint) {
    stop("resume_from_checkpoint requires enable_checkpoint = TRUE.", call.  = FALSE)
  }
  config
}

make_scenarios <- function() {
  pediatric_auc  <- c(0.50, 0.70, 1.00, 1.35, 1.80)
  pediatric_base <- c(0.09, 0.17, 0.30, 0.45, 0.60)
  e1_auc         <- c(0.35, 0.50, 0.70, 1.00, 1.35)

  make_one <- function(
    id, pediatric_dlt, adult_dlt, adult_auc,
    family, purpose, primary  = TRUE, true_mtd  = NA_integer_
  ) {
    if (is.na(true_mtd) && pediatric_dlt[1L] < 0.35) {
      true_mtd <- which.min(abs(pediatric_dlt - 0.30))
    }
    list(
      id            = id,
      family        = family,
      purpose       = purpose,
      primary       = primary,
      pediatric_dlt = pediatric_dlt,
      adult_dlt     = adult_dlt,
      pediatric_auc = pediatric_auc,
      adult_auc     = adult_auc,
      true_mtd      = if (length(true_mtd) && !is.na(true_mtd)) as.integer(true_mtd) else NA_integer_,
      no_mtd        = is.na(true_mtd)
    )
  }

  scenarios <- list(
    E0  = make_one("E0", pediatric_base, pediatric_base,
      c(0.50, 0.70, 1.00, 1.35, 1.80), "Exposure aligned",
      "Complete exposure and toxicity agreement"),
    E1  = make_one("E1", pediatric_base, c(0.04, 0.09, 0.17, 0.30, 0.45),
      e1_auc, "Exposure aligned", "Systematic exposure shift"),
    E2  = make_one("E2", pediatric_base, c(0.02, 0.04, 0.09, 0.17, 0.30),
      c(0.25, 0.35, 0.50, 0.70, 1.00), "Exposure aligned",
      "Larger exposure shift"),
    E3  = make_one("E3", pediatric_base, c(0.09, 0.13, 0.17, 0.37, 0.60),
      c(0.50, 0.60, 0.70, 1.17, 1.80), "Exposure aligned",
      "Irregular and many-to-one exposure mapping"),
    M1  = make_one("M1", c(0.30, 0.45, 0.60, 0.72, 0.82),
      c(0.17, 0.30, 0.45, 0.60, 0.72), e1_auc, "MTD location",
      "True pediatric MTD at P1", true_mtd  = 1L),
    M2  = make_one("M2", c(0.03, 0.06, 0.10, 0.17, 0.30),
      c(0.01, 0.03, 0.06, 0.10, 0.17), e1_auc, "MTD location",
      "True pediatric MTD at P5; P5 has no exposure match", true_mtd  = 5L),
    L1  = make_one("L1", pediatric_base, c(0.04, 0.09, 0.27, 0.30, 0.45),
      e1_auc, "Local conflict", "Conflict at the P2-A3 exposure match"),
    L2  = make_one("L2", pediatric_base, c(0.04, 0.09, 0.17, 0.20, 0.45),
      e1_auc, "Local conflict", "Conflict at the target P3-A4 match"),
    S1  = make_one("S1", c(0.35, 0.45, 0.55, 0.65, 0.75),
      c(0.25, 0.35, 0.45, 0.55, 0.65), e1_auc, "No MTD safety",
      "No acceptable pediatric dose; compatible adult match",
      true_mtd  = NA_integer_),
    S2  = make_one("S2", c(0.35, 0.45, 0.55, 0.65, 0.75),
      c(0.15, 0.25, 0.35, 0.45, 0.55), e1_auc, "No MTD safety",
      "No acceptable pediatric dose; adult match is 0.10 safer",
      true_mtd  = NA_integer_),
    G1  = make_one("G1", pediatric_base, c(0.02, 0.05, 0.10, 0.20, 0.35),
      e1_auc, "Global conflict", "Pediatric toxicity systematically higher",
      primary  = FALSE),
    G2  = make_one("G2", pediatric_base, c(0.12, 0.20, 0.27, 0.40, 0.55),
      e1_auc, "Global conflict", "Pediatric toxicity systematically lower",
      primary  = FALSE)
  )
  scenarios
}

scenario_truth_table <- function(scenarios  = make_scenarios()) {
  rows <- lapply(scenarios, function(s) {
      data.frame(
        scenario         = s$id,
        family           = s$family,
        primary          = s$primary,
        purpose          = s$purpose,
        true_mtd         = s$true_mtd,
        no_mtd           = s$no_mtd,
        pediatric_dlt    = paste(s$pediatric_dlt, collapse  = ","),
        adult_dlt        = paste(s$adult_dlt, collapse  = ","),
        pediatric_auc    = paste(s$pediatric_auc, collapse  = ","),
        adult_auc        = paste(s$adult_auc, collapse  = ","),
        stringsAsFactors = FALSE
      )
  })
  do.call(rbind, rows)
}


# -----------------------------------------------------------------------------
# Formal analysis constants (final matchconditional specification)
# -----------------------------------------------------------------------------
FORMAL_VERSION       <- "3.0.0"
FORMAL_FILE_VERSION  <- "PA_keyboard_ETOVL_formal_main_v03_borrowing_metrics_matchconditional"
FORMAL_PROTOCOL_DATE <- "2026-08-26"

FORMAL_SCENARIOS <- c(
  "E0", "E1", "E2", "E3",
  "M1", "M2", "L1", "L2",
  "S1", "S2"
)
FORMAL_METHODS               <- c("Keyboard", "PA-Keyboard", "PA-CRM")
FORMAL_METHOD_ORDER          <- setNames(seq_along(FORMAL_METHODS), FORMAL_METHODS)
FORMAL_SENSITIVITY_SCENARIOS <- c("E1", "E3", "L2", "S2")

FORMAL_SENSITIVITY_PLAN <- data.frame(
  group  = c(
    "AUC caliper", "AUC caliper",
    "Adult information", "Adult information", "Adult information",
    "Adult ESS cap", "Adult ESS cap",
    "Borrowing gate"
  ),
  variant_id  = c(
    "CAL16", "CAL28",
    "LEAD9", "LEAD12", "ADULT_COMPLETED",
    "CAP3", "CAP9",
    "GATE3"
  ),
  setting  = c(
    "AUC caliper log = 0.16",
    "AUC caliper log = 0.28",
    "Adult lead-in = 9 patients",
    "Adult lead-in = 12 patients",
    "Adult trial completed before pediatric start",
    "Maximum borrowed adult ESS = 3",
    "Maximum borrowed adult ESS = 9",
    "Borrowing starts at pediatric n = 3"
  ),
  stringsAsFactors  = FALSE
)

# Internal identifiers used by the ET-OVL simulation engine.
PA_KEYBOARD_METHOD         <- "PA-Keyboard"
PA_KEYBOARD_VARIANT_ORDER  <- c("PAK_original", "PAK_early", "PAK_capped")
PA_KEYBOARD_VARIANT_LABELS <- c(
  PAK_original = "PA-Keyboard original: nmin=6, uncapped delta*nA",
  PAK_early    = "PA-Keyboard early: nmin=3, uncapped delta*nA",
  PAK_capped   = "PA-Keyboard capped: nmin=6, max effective adult N=6"
)

