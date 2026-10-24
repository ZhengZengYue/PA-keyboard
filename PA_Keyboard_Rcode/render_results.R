# =============================================================================
# Regenerate publication outputs from a completed compact analysis
# =============================================================================
# Open PA_Keyboard.Rproj and click Source. No simulations are performed.
# Set output_dir to the directory containing the completed analysis RDS.
project_dir <- "."
output_dir  <- NULL
source(file.path(project_dir, "R", "00_setup.R"), encoding  = "UTF-8")
paths <- pak_initialize_paths(project_dir, output_dir)
pak_load_modules(paths)
pak_check_packages()
publication_results <- rebuild_publication_outputs(output_dir  = paths$output_dir)
