args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- if (length(file_arg) == 1L) {
  dirname(normalizePath(file_arg))
} else {
  getwd()
}
project_dir <- normalizePath(file.path(script_dir, ".."))
source(file.path(project_dir, "R", "load_endogenous_model.R"))

grid <- make_type_grid(501L)
param <- default_endogenous_parameters()
eligibility <- payg_eligibility_probability(grid$i, param)
stopifnot(all(diff(eligibility) > 0))
stopifnot(all(eligibility > 0), all(eligibility < 1))

coarse <- solve_steady_state_endogenous(
  param,
  grid = make_type_grid(101L),
  old_weight = 0.90,
  tolerance = 1e-6,
  residual_tolerance = 1e-4
)
solution <- solve_steady_state_endogenous(
  param,
  initial = as.list(coarse$state),
  grid = grid,
  old_weight = 0.90,
  tolerance = 1e-8,
  residual_tolerance = 2e-6
)
stopifnot(solution$validated)
shares <- c(
  informal = solution$cohort$informal_share,
  capitalization = solution$cohort$capitalization_share,
  payg = solution$cohort$payg_share
)
stopifnot(abs(sum(shares) - 1) < 1e-8)
stopifnot(all(shares > 0.01))
stopifnot(solution$cohort$diagnostics$ordered_three_choice)
stopifnot(all(solution$cohort$micro$consumption_young > 0))
stopifnot(all(solution$cohort$micro$consumption_old > 0))
stopifnot(all(solution$cohort$micro$voluntary_saving >= 0))
stopifnot(abs(
  solution$identities[["supply_minus_income"]]
) < 1e-7)
stopifnot(max(abs(solution$fixed_point_residuals)) < 2e-6)
roots <- solution$cohort$integration$roots
stopifnot(length(roots) == 2L)
stopifnot(roots[1L] > 0, roots[2L] < 1, roots[1L] < roots[2L])

cat("Todas las pruebas de eleccion endogena pasaron.\n")
cat(sprintf(
  "Participaciones: informal=%.4f, capitalizacion=%.4f, PAYG=%.4f\n",
  shares[["informal"]],
  shares[["capitalization"]],
  shares[["payg"]]
))
cat(sprintf("Umbrales: i1=%.6f, i2=%.6f\n", roots[1L], roots[2L]))
