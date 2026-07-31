model_files <- c(
  "utils.R",
  "parameters.R",
  "density.R",
  "households.R",
  "economy.R",
  "steady_state.R"
)

model_root <- normalizePath(
  file.path(dirname(sys.frame(1)$ofile %||% "R/load_model.R"), ".."),
  mustWork = FALSE
)

for (model_file in model_files) {
  source(file.path(model_root, "R", model_file), local = .GlobalEnv)
}
