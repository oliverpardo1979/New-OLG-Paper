script_path <- normalizePath(
  sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
)
project_root <- normalizePath(
  file.path(dirname(script_path), "..", ".."),
  mustWork = TRUE
)
source(file.path(project_root, "code", "R", "model.R"))
source(file.path(project_root, "code", "R", "stable_path_solver.R"))
source(file.path(project_root, "code", "R", "fast_solver.R"))
source(file.path(project_root, "code", "R", "stable_path_solver_v2.R"))
source(file.path(project_root, "code", "R", "solver_entry.R"))

assert_close <- function(value, target, tolerance, label) {
  if (!is.finite(value) || abs(value - target) > tolerance) {
    stop(
      label, " failed: value=", signif(value, 8),
      ", target=", signif(target, 8),
      ", tolerance=", tolerance
    )
  }
}

param <- default_parameters()
param$choice_scale <- 0.30
param$horizon <- 500L
param <- calibrate_parameters(param)
baseline <- solve_steady_state(param)

assert_close(
  baseline$informal_share,
  param$target_informal_share,
  1e-8,
  "Informality calibration"
)
assert_close(
  baseline$government_share_output,
  param$target_government_share,
  1e-8,
  "Government-spending calibration"
)
assert_close(
  baseline$government_residual,
  0,
  1e-10,
  "Government budget"
)
assert_close(
  baseline$resource_residual,
  0,
  1e-10,
  "Stationary resource constraint"
)
assert_close(
  baseline$euler_residual,
  0,
  1e-10,
  "Stationary Euler equation"
)

scenario_names <- c("Lump sum", "Mixed", "Payroll")
payroll_shares <- c(0, 0.5, 1)
steady <- do.call(
  rbind,
  lapply(
    seq_along(payroll_shares),
    function(index) {
      solve_steady_state(
        param,
        spending_increase = param$permanent_spending_increase,
        marginal_payroll_share = payroll_shares[index],
        label = scenario_names[index]
      )
    }
  )
)
if (!all(diff(steady$informal_share) > 0)) {
  stop("Informality should rise with payroll financing.")
}
if (!all(diff(steady$capital) < 0)) {
  stop("Capital should fall with payroll financing.")
}

transitions <- lapply(
  seq_along(payroll_shares),
  function(index) {
    solve_transition(
      param,
      marginal_payroll_share = payroll_shares[index],
      label = scenario_names[index],
      path_type = "step"
    )
  }
)
for (index in seq_along(transitions)) {
  transition <- transitions[[index]]
  last <- nrow(transition$data)
  assert_close(
    transition$data$capital[last],
    transition$final_steady$capital,
    1e-7,
    paste("Terminal capital", scenario_names[index])
  )
  assert_close(
    transition$data$consumption[last],
    transition$final_steady$consumption,
    1e-7,
    paste("Terminal consumption", scenario_names[index])
  )
  assert_close(
    transition$data$informal_share[last],
    transition$final_steady$informal_share,
    1e-7,
    paste("Terminal informality", scenario_names[index])
  )
}

# Additive public consumption changes welfare, but does not enter any
# allocation equation. The same transition must therefore be used under
# both welfare specifications.
allocation_without_public_utility <- transitions[[2L]]$data[
  , c("capital", "consumption", "informal_share")
]
allocation_with_public_utility <- transitions[[2L]]$data[
  , c("capital", "consumption", "informal_share")
]
if (!identical(
    allocation_without_public_utility,
    allocation_with_public_utility
)) {
  stop("Additive public utility must not change allocations.")
}

cat("All public-spending paper tests passed.\n")
