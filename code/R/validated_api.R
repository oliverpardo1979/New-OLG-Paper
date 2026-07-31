solve_steady_state_validated <- function(
    param = default_parameters(),
    initial = NULL,
    grid = make_type_grid(1001L),
    old_weight = 0.97,
    tolerance = 1e-9,
    residual_tolerance = 1e-5,
    max_iterations = 6000L,
    verbose = FALSE
) {
  solution <- solve_steady_state(
    param = param,
    initial = initial,
    grid = grid,
    old_weight = old_weight,
    tolerance = tolerance,
    max_iterations = max_iterations,
    verbose = verbose
  )

  max_residual <- max(abs(solution$fixed_point_residuals))
  solution$validated <- solution$converged &&
    max_residual < residual_tolerance
  solution$residual_tolerance <- residual_tolerance

  if (!solution$validated) {
    warning(sprintf(
      paste(
        "El estado estacionario no superó la validación:",
        "converged=%s, max_residual=%.3e."
      ),
      solution$converged,
      max_residual
    ))
  }
  solution
}

solve_transition_validated <- function(...) {
  transition <- solve_transition(...)
  final <- transition$final_solution
  last <- nrow(transition$path)

  transition$path$consumption_young[last] <-
    final$cohort$consumption_young
  transition$path$consumption_old[last] <-
    final$cohort$consumption_old_next
  transition$path$informal_share[last] <-
    final$cohort$informal_share
  transition$path$payg_contributions[last] <-
    final$cohort$payg_contributions
  transition$path$payg_benefits[last] <-
    final$cohort$payg_benefits_next
  transition$path$budget_residual[last] <-
    final$budget$residual

  interior_budget <- transition$path$budget_residual[
    2:(last - 1L)
  ]
  transition$max_interior_budget_residual <- max(
    abs(interior_budget),
    na.rm = TRUE
  )
  transition$validated <- transition$converged_internal &&
    transition$terminal_consistent
  transition
}
