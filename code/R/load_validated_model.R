this_file <- sys.frame(1)$ofile
if (is.null(this_file)) {
  stop(
    "Cargue este archivo con source('R/load_validated_model.R').",
    call. = FALSE
  )
}
project_dir <- normalizePath(file.path(dirname(this_file), ".."))

model_files <- c(
  "utils.R",
  "parameters.R",
  "density.R",
  "households.R",
  "economy.R",
  "steady_state.R",
  "refined_integration.R",
  "transition.R",
  "validated_api.R"
)

for (model_file in model_files) {
  source(
    file.path(project_dir, "R", model_file),
    local = .GlobalEnv
  )
}
