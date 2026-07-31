this_file <- sys.frame(1)$ofile
if (is.null(this_file)) {
  stop(
    "Cargue este archivo con source('R/load_endogenous_model.R').",
    call. = FALSE
  )
}
project_dir <- normalizePath(file.path(dirname(this_file), ".."))
source(
  file.path(project_dir, "R", "load_validated_model.R"),
  local = .GlobalEnv
)
source(
  file.path(project_dir, "R", "endogenous_choice.R"),
  local = .GlobalEnv
)
source(
  file.path(project_dir, "R", "calibration_becerra2026.R"),
  local = .GlobalEnv
)
