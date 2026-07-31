realize_initial_old_cohort <- function(
    initial_solution,
    new_param,
    realized_return,
    realized_consumption_price
) {
  old_param <- initial_solution$param
  old_prices <- initial_solution$prices
  intervals <- initial_solution$cohort$integration
  left_i <- intervals$left_i
  right_i <- intervals$right_i
  formal <- intervals$formal
  integrate_piecewise <- make_interval_integrator(left_i, right_i)

  old_left <- household_outcomes(
    left_i,
    formal,
    old_param,
    old_prices
  )
  old_right <- household_outcomes(
    right_i,
    formal,
    old_param,
    old_prices
  )

  realized_prices <- list(r = realized_return)

  realized_consumption <- function(i, old_outcome) {
    wage <- individual_wage(i, formal, old_param, old_prices)
    pension <- pension_accounts(
      i,
      formal,
      wage,
      old_param,
      realized_prices
    )

    # Derechos contributivos y saldos de ahorro quedan adquiridos con las
    # reglas antiguas. Las pensiones básicas vigentes se pagan con la política
    # nueva desde el primer periodo de la reforma.
    benefit <-
      pension$funded_benefit +
      pension$payg_benefit +
      (1 + realized_return) *
        old_param$universal_basic_saving / (1 + old_param$g) +
      (1 + realized_return) *
        old_param$conditional_basic_saving *
        as.numeric(!formal) / (1 + old_param$g) +
      (1 - old_param$m) *
        new_param$universal_basic_pension +
      (1 - old_param$m) *
        new_param$conditional_basic_pension *
        as.numeric(!formal)

    list(
      consumption_old = (
        (1 + realized_return) * old_outcome$voluntary_saving +
          benefit
      ) / ((1 - old_param$m) * realized_consumption_price),
      payg_benefit = pension$payg_benefit
    )
  }

  left_realized <- realized_consumption(left_i, old_left)
  right_realized <- realized_consumption(right_i, old_right)

  list(
    consumption_old = (1 - old_param$m) / (1 + old_param$n) *
      integrate_piecewise(
        left_realized$consumption_old,
        right_realized$consumption_old
      ),
    payg_benefits = integrate_piecewise(
      left_realized$payg_benefit,
      right_realized$payg_benefit
    ) / (1 + old_param$n),
    informal_share = initial_solution$cohort$informal_share
  )
}

solve_transition <- function(
    initial_param,
    final_param,
    periods = 40L,
    initial_solution = NULL,
    final_solution = NULL,
    grid = make_type_grid(1001L),
    old_weight = 0.95,
    tolerance = 1e-6,
    terminal_tolerance = 2e-3,
    max_iterations = 2000L,
    verbose = FALSE
) {
  validate_parameters(initial_param)
  validate_parameters(final_param)
  if (periods < 6L) {
    stop_model("La transición requiere al menos seis periodos.")
  }

  fixed_demography <- c("m", "n", "g")
  changed <- fixed_demography[vapply(
    fixed_demography,
    function(name) initial_param[[name]] != final_param[[name]],
    logical(1)
  )]
  if (length(changed) > 0L) {
    stop_model(
      "Esta versión mantiene fijos durante la reforma: %s.",
      paste(changed, collapse = ", ")
    )
  }

  if (is.null(initial_solution)) {
    initial_solution <- solve_steady_state(
      initial_param,
      grid = grid
    )
  }
  if (is.null(final_solution)) {
    final_solution <- solve_steady_state(
      final_param,
      grid = grid
    )
  }
  if (!initial_solution$converged || !final_solution$converged) {
    stop_model(
      "Los dos estados estacionarios deben converger antes de la transición."
    )
  }

  horizon <- periods + 1L
  k <- c(
    initial_solution$state[["k"]],
    seq(
      initial_solution$state[["k"]],
      final_solution$state[["k"]],
      length.out = periods
    )
  )
  formal_labor <- c(
    initial_solution$state[["formal_labor"]],
    seq(
      initial_solution$state[["formal_labor"]],
      final_solution$state[["formal_labor"]],
      length.out = periods
    )
  )
  informal_labor <- c(
    initial_solution$state[["informal_labor"]],
    seq(
      initial_solution$state[["informal_labor"]],
      final_solution$state[["informal_labor"]],
      length.out = periods
    )
  )
  consumption_tax <- c(
    initial_solution$state[["consumption_tax"]],
    seq(
      initial_solution$state[["consumption_tax"]],
      final_solution$state[["consumption_tax"]],
      length.out = periods
    )
  )

  internal_gap <- Inf
  terminal_gap <- Inf
  converged_internal <- FALSE
  cohort_path <- vector("list", horizon)

  for (iteration in seq_len(max_iterations)) {
    old_k <- k
    old_formal <- formal_labor
    old_informal <- informal_labor
    old_tax <- consumption_tax

    price_path <- vector("list", horizon)
    for (t in seq_len(horizon)) {
      period_param <- if (t == 1L) initial_param else final_param
      period_param <- copy_list(period_param)
      period_param$tau_consumption <- consumption_tax[t]
      price_path[[t]] <- factor_prices(
        k[t],
        formal_labor[t],
        consumption_tax[t],
        period_param
      )
    }

    consumption_young <- rep(NA_real_, horizon)
    consumption_old <- rep(NA_real_, horizon)
    payg_contributions <- rep(NA_real_, horizon)
    payg_benefits <- rep(NA_real_, horizon)
    informal_share <- rep(NA_real_, horizon)
    old_informal_share <- rep(NA_real_, horizon)
    implied_k <- k
    implied_formal <- formal_labor
    implied_informal <- informal_labor
    implied_tax <- consumption_tax
    budget_residual <- rep(NA_real_, horizon)

    consumption_young[1L] <-
      initial_solution$cohort$consumption_young
    consumption_old[1L] <-
      initial_solution$cohort$consumption_old_next
    payg_contributions[1L] <-
      initial_solution$cohort$payg_contributions
    payg_benefits[1L] <-
      initial_solution$cohort$payg_benefits_next
    informal_share[1L] <-
      initial_solution$cohort$informal_share
    old_informal_share[1L] <-
      initial_solution$cohort$informal_share

    surprised_old <- realize_initial_old_cohort(
      initial_solution,
      final_param,
      price_path[[2L]]$r,
      price_path[[2L]]$p
    )
    consumption_old[2L] <- surprised_old$consumption_old
    payg_benefits[2L] <- surprised_old$payg_benefits
    old_informal_share[2L] <- surprised_old$informal_share

    for (t in 2:periods) {
      decision_prices <- list(
        r = price_path[[t + 1L]]$r,
        wage_formal = price_path[[t]]$wage_formal,
        wage_informal = (1 + consumption_tax[t]) *
          final_param$A_informal,
        p = 1 + consumption_tax[t],
        q = 1 + consumption_tax[t + 1L]
      )

      period_param <- copy_list(final_param)
      period_param$tau_consumption <- consumption_tax[t]
      cohort <- evaluate_cohort(
        period_param,
        decision_prices,
        grid
      )
      cohort_path[[t]] <- cohort

      implied_k[t + 1L] <- cohort$capital_next
      implied_formal[t] <- cohort$formal_labor
      implied_informal[t] <- cohort$informal_labor
      consumption_young[t] <- cohort$consumption_young
      consumption_old[t + 1L] <- cohort$consumption_old_next
      payg_contributions[t] <- cohort$payg_contributions
      payg_benefits[t + 1L] <- cohort$payg_benefits_next
      informal_share[t] <- cohort$informal_share
      old_informal_share[t + 1L] <- cohort$informal_share
    }

    for (t in 2:periods) {
      period_param <- copy_list(final_param)
      period_param$tau_consumption <- consumption_tax[t]
      budget_cohort <- cohort_path[[t]]
      budget <- government_budget(
        period_param,
        price_path[[t]],
        k[t],
        budget_cohort,
        consumption_old[t],
        payg_benefits[t],
        old_informal_share[t]
      )
      implied_tax[t] <- budget$implied_tax
      budget_residual[t] <- budget$residual
    }

    # k[1] es el estado inicial. k[2] fue predeterminado por la cohorte que
    # ahorró antes de la reforma. El último estado se fija como condición
    # terminal y su discrepancia se reporta, no se oculta.
    if (periods >= 3L) {
      update_k <- 3:periods
      k[update_k] <- old_weight * k[update_k] +
        (1 - old_weight) * implied_k[update_k]
    }
    update_flow <- 2:periods
    formal_labor[update_flow] <-
      old_weight * formal_labor[update_flow] +
      (1 - old_weight) * implied_formal[update_flow]
    informal_labor[update_flow] <-
      old_weight * informal_labor[update_flow] +
      (1 - old_weight) * implied_informal[update_flow]
    consumption_tax[update_flow] <-
      old_weight * consumption_tax[update_flow] +
      (1 - old_weight) * implied_tax[update_flow]

    internal_gap <- max(
      max_relative_gap(k[3:periods], old_k[3:periods]),
      max_relative_gap(
        formal_labor[2:periods],
        old_formal[2:periods]
      ),
      max_relative_gap(
        informal_labor[2:periods],
        old_informal[2:periods]
      ),
      max_relative_gap(
        consumption_tax[2:periods],
        old_tax[2:periods]
      )
    )
    terminal_gap <- abs(
      implied_k[horizon] - final_solution$state[["k"]]
    ) / final_solution$state[["k"]]

    if (verbose && (iteration == 1L || iteration %% 50L == 0L)) {
      message(sprintf(
        "iter=%d gap=%.3e brecha_terminal=%.3e",
        iteration,
        internal_gap,
        terminal_gap
      ))
    }

    if (internal_gap < tolerance) {
      converged_internal <- TRUE
      break
    }
  }

  path <- data.frame(
    period = seq_len(horizon),
    capital = k,
    formal_labor = formal_labor,
    informal_labor = informal_labor,
    consumption_tax = consumption_tax,
    interest_rate = vapply(price_path, `[[`, numeric(1), "r"),
    formal_wage = vapply(
      price_path,
      `[[`,
      numeric(1),
      "wage_formal"
    ),
    informal_wage = vapply(
      price_path,
      `[[`,
      numeric(1),
      "wage_informal"
    ),
    consumption_young = consumption_young,
    consumption_old = consumption_old,
    informal_share = informal_share,
    payg_contributions = payg_contributions,
    payg_benefits = payg_benefits,
    budget_residual = budget_residual
  )

  structure(
    list(
      converged_internal = converged_internal,
      terminal_consistent = terminal_gap < terminal_tolerance,
      iterations = iteration,
      internal_gap = internal_gap,
      terminal_gap = terminal_gap,
      path = path,
      initial_solution = initial_solution,
      final_solution = final_solution,
      timing = paste(
        "La reforma comienza en t=2.",
        "Los derechos contributivos de la cohorte vieja quedan adquiridos;",
        "las pensiones básicas nuevas se pagan desde t=2."
      )
    ),
    class = "pension_transition"
  )
}

print.pension_transition <- function(x, ...) {
  cat("Transición del modelo pensional\n")
  cat(sprintf(
    "  Punto fijo interior: %s (%d iteraciones)\n",
    ifelse(x$converged_internal, "sí", "no"),
    x$iterations
  ))
  cat(sprintf(
    "  Brecha terminal del capital: %.3e\n",
    x$terminal_gap
  ))
  cat(sprintf(
    "  Condición terminal: %s\n",
    ifelse(x$terminal_consistent, "satisfecha", "no satisfecha")
  ))
  invisible(x)
}
