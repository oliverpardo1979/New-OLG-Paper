args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[1]))
} else {
  normalizePath("code/tests/test_model.R")
}
paper_dir <- normalizePath(file.path(dirname(script_path), "..", ".."))

source(file.path(paper_dir, "code", "R", "model.R"))
source(file.path(paper_dir, "code", "R", "solver.R"))

par <- default_parameters()
d0 <- 0.150170218
reference_mix <- c(0.45, 0.35, 0.20)
fit <- calibrate_parameters(par, d0, reference_mix)
stopifnot(isTRUE(fit$converged))
par <- fit$parameters

eq <- solve_stationary_equilibrium(d0, reference_mix, par)
stopifnot(isTRUE(eq$converged))
stopifnot(max(abs(eq$residuals)) < 1e-5)
stopifnot(abs(eq$informality - 0.554) < 2e-4)
stopifnot(abs(eq$social_security_share - 0.057) < 2e-4)
stopifnot(all(c(eq$tax_consumption, eq$tax_payroll, eq$tax_capital) >= 0))
stopifnot(abs(par$beta * (1 + par$bond_rate) - 1) < 1e-12)
stopifnot(abs(par$old_age_outlay - 0.426) < 0.001)
stopifnot(abs(sum(reference_mix) - 1) < 1e-12)
revenue_check <-
  eq$tax_consumption * (
    eq$consumption_formal +
      par$informal_vat_coverage * eq$price_informal *
        eq$consumption_informal
  ) +
  eq$tax_payroll * eq$wage_formal * eq$labor_formal +
  eq$tax_capital * eq$rental_rate * eq$capital
stopifnot(abs(revenue_check - eq$required_revenue) < 1e-7)
resource_check <-
  eq$consumption_formal + par$delta * eq$capital +
  eq$public_absorption + eq$collection_cost - eq$output_formal
stopifnot(abs(resource_check) < 1e-6)

grid <- evaluate_tax_grid(d0, par, step = 0.25)
stopifnot(any(grid$converged))
stopifnot(sum(grid$is_optimum) == 1)

frictionless <- par
frictionless$collection_cost_consumption <- 0
frictionless$collection_cost_payroll <- 0
frictionless$collection_cost_capital <- 0
frictionless_grid <- evaluate_tax_grid(d0, frictionless, step = 0.25)
frictionless_best <- optimal_grid_row(frictionless_grid)
stopifnot(abs(frictionless_best$mix_consumption - 1) < 1e-12)

demographics <- read.csv(file.path(paper_dir, "data", "demographic_path.csv"))
short_path <- solve_demographic_path(
  demographics[c(1, 14, 27), ],
  as.numeric(optimal_grid_row(grid)[1, c(
    "mix_consumption", "mix_payroll", "mix_capital"
  )]),
  par
)
stopifnot(nrow(short_path) == 3)
stopifnot(all(diff(short_path$social_security_share) > 0))

cat("All aging-tax-mix model tests passed.\n")
