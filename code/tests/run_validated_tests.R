args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- if (length(file_arg) == 1L) {
  dirname(normalizePath(file_arg))
} else {
  getwd()
}
project_dir <- normalizePath(file.path(script_dir, ".."))
source(file.path(project_dir, "R", "load_validated_model.R"))

grid <- make_type_grid(501L)
stopifnot(abs(grid$raw_mass - 1) < 1e-6)
stopifnot(abs(sum(grid$weights) - 1) < 1e-12)

compete <- default_parameters("compete")
compete_ss <- solve_steady_state_validated(
  compete,
  grid = grid,
  tolerance = 1e-9,
  residual_tolerance = 2e-6
)
stopifnot(compete_ss$validated)
stopifnot(all(compete_ss$cohort$micro$consumption_young > 0))
stopifnot(all(compete_ss$cohort$micro$consumption_old > 0))
stopifnot(all(compete_ss$cohort$micro$voluntary_saving >= 0))
stopifnot(abs(
  compete_ss$identities[["supply_minus_income"]]
) < 1e-7)

pillars <- default_parameters("pillars")
pillars$payg_income_threshold <- 1.0
pillars_ss <- solve_steady_state_validated(
  pillars,
  grid = grid,
  tolerance = 1e-9,
  residual_tolerance = 2e-6
)
stopifnot(pillars_ss$validated)
stopifnot(abs(
  pillars_ss$identities[["supply_minus_income"]]
) < 1e-7)

no_reform <- solve_transition_validated(
  initial_param = compete,
  final_param = compete,
  periods = 12L,
  initial_solution = compete_ss,
  final_solution = compete_ss,
  grid = grid,
  tolerance = 1e-8,
  max_iterations = 500L
)
stopifnot(no_reform$validated)
stopifnot(no_reform$terminal_gap < 1e-5)
stopifnot(no_reform$max_interior_budget_residual < 2e-6)

cat("Todas las pruebas validadas pasaron.\n")
