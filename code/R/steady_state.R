initial_steady_state_guess <- function(param, grid, cutoff = 0.468) {
  informal <- as.numeric(grid$i <= cutoff)
  formal <- 1 - informal

  list(
    k = 0.05,
    formal_labor = weighted_sum(
      formal * param$base_formal * exp(param$rho_formal * grid$i),
      grid$weights
    ),
    informal_labor = weighted_sum(
      informal * param$base_informal *
        exp(param$rho_informal * grid$i),
      grid$weights
    ),
    consumption_tax = param$tau_consumption
  )
}

solve_steady_state <- function(
    param = default_parameters(),
    initial = NULL,
    grid = make_type_grid(),
    old_weight = 0.97,
    tolerance = 1e-7,
    max_iterations = 3000L,
    verbose = FALSE
) {
  validate_parameters(param)
  assert_scalar(old_weight, "old_weight", 0, 0.999999)
  assert_scalar(tolerance, "tolerance", .Machine$double.eps, Inf)

  if (is.null(initial)) {
    initial <- initial_steady_state_guess(param, grid)
  }

  state <- c(
    k = initial$k,
    formal_labor = initial$formal_labor,
    informal_labor = initial$informal_labor,
    consumption_tax = initial$consumption_tax
  )
  converged <- FALSE
  gap <- Inf

  for (iteration in seq_len(max_iterations)) {
    iter_param <- copy_list(param)
    iter_param$tau_consumption <- state[["consumption_tax"]]
    prices <- factor_prices(
      state[["k"]],
      state[["formal_labor"]],
      state[["consumption_tax"]],
      iter_param
    )
    cohort <- evaluate_cohort(iter_param, prices, grid)
    budget <- government_budget(
      iter_param,
      prices,
      state[["k"]],
      cohort,
      cohort$consumption_old_next,
      cohort$payg_benefits_next
    )

    target <- c(
      k = cohort$capital_next,
      formal_labor = cohort$formal_labor,
      informal_labor = cohort$informal_labor,
      consumption_tax = budget$implied_tax
    )

    if (target[["k"]] <= 0 ||
        target[["formal_labor"]] <= 0 ||
        target[["consumption_tax"]] <= -0.999999 ||
        any(!is.finite(target))) {
      stop_model(
        "La iteración %s produjo un estado económicamente inválido.",
        iteration
      )
    }

    updated <- old_weight * state + (1 - old_weight) * target
    gap <- max_relative_gap(updated, state)
    state <- updated

    if (verbose && (iteration == 1L || iteration %% 50L == 0L)) {
      message(sprintf(
        paste(
          "iter=%d gap=%.3e k=%.6f nf=%.6f",
          "ns=%.6f tauc=%.6f informalidad=%.4f"
        ),
        iteration,
        gap,
        state[["k"]],
        state[["formal_labor"]],
        state[["informal_labor"]],
        state[["consumption_tax"]],
        cohort$informal_share
      ))
    }

    if (gap < tolerance) {
      converged <- TRUE
      break
    }
  }

  final_param <- copy_list(param)
  final_param$tau_consumption <- state[["consumption_tax"]]
  prices <- factor_prices(
    state[["k"]],
    state[["formal_labor"]],
    state[["consumption_tax"]],
    final_param
  )
  cohort <- evaluate_cohort(final_param, prices, grid)
  budget <- government_budget(
    final_param,
    prices,
    state[["k"]],
    cohort,
    cohort$consumption_old_next,
    cohort$payg_benefits_next
  )
  identities <- steady_state_identities(
    final_param,
    prices,
    state[["k"]],
    cohort
  )

  fixed_point_residuals <- c(
    capital = cohort$capital_next - state[["k"]],
    formal_labor = cohort$formal_labor - state[["formal_labor"]],
    informal_labor = cohort$informal_labor - state[["informal_labor"]],
    consumption_tax = budget$implied_tax -
      state[["consumption_tax"]],
    government_budget = budget$residual
  )

  structure(
    list(
      converged = converged,
      iterations = iteration,
      iteration_gap = gap,
      state = state,
      param = final_param,
      prices = prices,
      cohort = cohort,
      budget = budget,
      identities = identities,
      fixed_point_residuals = fixed_point_residuals,
      grid = grid
    ),
    class = "pension_steady_state"
  )
}

print.pension_steady_state <- function(x, ...) {
  cat("Estado estacionario del modelo pensional\n")
  cat(sprintf("  Convergió: %s (%d iteraciones)\n",
    ifelse(x$converged, "sí", "no"), x$iterations))
  cat(sprintf("  Capital: %.8f\n", x$state[["k"]]))
  cat(sprintf("  Trabajo formal: %.8f\n",
    x$state[["formal_labor"]]))
  cat(sprintf("  Trabajo informal: %.8f\n",
    x$state[["informal_labor"]]))
  cat(sprintf("  Informalidad: %.4f%%\n",
    100 * x$cohort$informal_share))
  cat(sprintf("  IVA: %.4f%%\n",
    100 * x$state[["consumption_tax"]]))
  cat(sprintf("  Retorno: %.4f%%\n", 100 * x$prices$r))
  cat(sprintf("  Cruces formal/informal: %d\n",
    x$cohort$diagnostics$number_of_switches))
  cat("  Máximo residuo de punto fijo: ")
  cat(sprintf("%.3e\n", max(abs(x$fixed_point_residuals))))
  invisible(x)
}
