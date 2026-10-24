# =============================================================================
# Project setup and output paths
# =============================================================================
# This module validates the project layout, resolves user-selected paths, and
# loads the statistical methods in their required order. It does not run a
# simulation or modify the working directory.
#
# Functions:
#   pak_initialize_paths()     Validate project and results directories.
#   pak_load_modules()         Load the seven statistical study modules.
#   pak_check_packages()       Check packages required for publication outputs.
#   pak_load_workflows()       Load only workflow definitions for tests.
#   pak_initialize_workers()   Load modules and export workflows to PSOCK workers.
# =============================================================================

pak_initialize_paths <- function(project_dir, output_dir  = NULL) {
  if (!is.character(project_dir) || length(project_dir) != 1L ||
    is.na(project_dir) || !nzchar(project_dir)) {
    stop("Set 'project_dir' to the extracted project folder, or open PA_Keyboard.Rproj before running the script.", call.  = FALSE)
  }
  root         <- normalizePath(project_dir, winslash  = "/", mustWork  = TRUE)
  module_files <- file.path(root, "R", sprintf("%02d_%s.R", 1:7,
      c("design", "keyboard", "pa_crm", "pa_keyboard", "simulation", "results", "figures")))
  if (!all(file.exists(module_files))) {
    stop("Cannot find the seven R/ study modules under: ", root,
      ". Extract the complete ZIP, or correct 'project_dir'.", call.  = FALSE)
  }
  if (is.null(output_dir)) output_dir <- file.path(root, "results")
  if (!is.character(output_dir) || length(output_dir) != 1L ||
    is.na(output_dir) || !nzchar(output_dir)) {
    stop("'output_dir' must be NULL or a nonempty directory path.", call.  = FALSE)
  }
  # Relative output directories are interpreted from the project root.
  # Recognize Unix, home-relative, Windows drive, and network-share paths.
  is_absolute <- startsWith(output_dir, "/") || startsWith(output_dir, "~") ||
  grepl("^[A-Za-z]:", output_dir) || startsWith(output_dir, "\\\\")
  if (!is_absolute) output_dir <- file.path(root, output_dir)
  output_dir <- path.expand(output_dir)
  dir.create(output_dir, recursive  = TRUE, showWarnings  = FALSE)
  if (!dir.exists(output_dir)) stop("Cannot create output directory: ", output_dir, call.  = FALSE)
  list(
    project_dir  = root,
    output_dir   = normalizePath(output_dir, winslash  = "/", mustWork  = TRUE),
    entry_script = file.path(root, "run_analysis.R"),
    module_files = module_files
  )
}

pak_load_modules <- function(paths) {
  for (file in paths$module_files) source(file, local  = FALSE, encoding  = "UTF-8")
  invisible(paths$module_files)
}

pak_check_packages <- function() {
  required <- c("openxlsx", "dplyr", "tidyr", "ggplot2")
  missing  <- required[!vapply(required, requireNamespace, logical(1L), quietly  = TRUE)]
  if (length(missing)) {
    stop("Required R packages not installed: ", paste(missing, collapse  = ", "),
      ". See README.txt.", call.  = FALSE)
  }
  invisible(TRUE)
}

# Read function definitions only; no settings, package checks or analysis run.
pak_load_workflows <- function(project_dir, envir  = .GlobalEnv) {
  expressions <- parse(file.path(project_dir, "run_analysis.R"), encoding  = "UTF-8")
  required    <- c("pa_keyboard_borrow_at_dose", "simulate_pa_keyboard_pediatric",
    "simulate_formal_replication")
  found <- character()
  for (expr in expressions) {
    if (is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
      is.symbol(expr[[2L]]) && is.call(expr[[3L]]) &&
      identical(expr[[3L]][[1L]], as.name("function"))) {
      name <- as.character(expr[[2L]])
      if (name %in% required) {
        eval(expr, envir  = envir)
        found <- c(found, name)
      }
    }
  }
  if (!setequal(found, required)) stop("Missing study workflow definitions.")
  invisible(found)
}

pak_initialize_workers <- function(cl, project_dir) {
  parallel::clusterCall(cl, function(root) {
      source(file.path(root, "R", "00_setup.R"), encoding  = "UTF-8")
      paths <- pak_initialize_paths(root, file.path(tempdir(), "PAK_worker"))
      pak_load_modules(paths)
      NULL
    }, normalizePath(project_dir, mustWork  = TRUE))
  parallel::clusterExport(cl,
    c("pa_keyboard_borrow_at_dose", "simulate_pa_keyboard_pediatric",
      "simulate_formal_replication"), envir  = .GlobalEnv)
  invisible(NULL)
}
