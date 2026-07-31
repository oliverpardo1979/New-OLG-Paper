# Este archivo reemplaza la agregación sobre nodos por una integración por
# intervalos. Debe cargarse después de economy.R.

make_type_grid_from_nodes <- function(i) {
  i <- sort(unique(c(0, i, 1)))
  if (length(i) < 3L || any(diff(i) <= 0)) {
    stop_model("Los nodos de integración no son válidos.")
  }

  density <- type_density(i)
  dx <- diff(i)
  raw_mass <- sum(dx * (density[-length(i)] + density[-1L]) / 2)

  list(
    i = i,
    density = density,
    raw_mass = raw_mass
  )
}

utility_difference <- function(i, param, prices) {
  formal <- household_outcomes(i, TRUE, param, prices)
  informal <- household_outcomes(i, FALSE, param, prices)
  formal$decision_utility - informal$decision_utility
}

find_choice_roots <- function(micro, param, prices) {
  changes <- which(diff(as.integer(micro$formal)) != 0)
  if (length(changes) == 0L) {
    return(numeric(0))
  }

  vapply(changes, function(index) {
    lower <- micro$i[index]
    upper <- micro$i[index + 1L]

    if (param$regime == "compete" &&
        lower <= param$payg_type_threshold &&
        upper >= param$payg_type_threshold) {
      return(param$payg_type_threshold)
    }

    tryCatch(
      uniroot(
        function(x) utility_difference(x, param, prices),
        interval = c(lower, upper),
        tol = 1e-11
      )$root,
      error = function(e) (lower + upper) / 2
    )
  }, numeric(1))
}

prepare_choice_intervals <- function(param, prices, grid) {
  base_nodes <- grid$i
  if (param$regime == "compete") {
    base_nodes <- c(base_nodes, param$payg_type_threshold)
  }

  diagnostic_grid <- make_type_grid_from_nodes(base_nodes)
  diagnostic_micro <- choose_sector(diagnostic_grid$i, param, prices)
  roots <- find_choice_roots(diagnostic_micro, param, prices)
  integration_grid <- make_type_grid_from_nodes(c(diagnostic_grid$i, roots))

  left_i <- integration_grid$i[-length(integration_grid$i)]
  right_i <- integration_grid$i[-1L]
  midpoint_i <- (left_i + right_i) / 2
  interval_formal <- choose_sector(
    midpoint_i,
    param,
    prices
  )$formal

  list(
    diagnostic_micro = diagnostic_micro,
    roots = roots,
    nodes = integration_grid$i,
    left_i = left_i,
    right_i = right_i,
    interval_formal = interval_formal
  )
}

make_interval_integrator <- function(left_i, right_i) {
  dx <- right_i - left_i
  left_density <- type_density(left_i)
  right_density <- type_density(right_i)
  mass <- sum(dx * (left_density + right_density) / 2)

  function(left_value, right_value) {
    sum(dx * (
      left_density * left_value + right_density * right_value
    ) / 2) / mass
  }
}

evaluate_cohort <- function(param, prices, grid) {
  intervals <- prepare_choice_intervals(param, prices, grid)
  left_i <- intervals$left_i
  right_i <- intervals$right_i
  formal <- intervals$interval_formal
  integrate_piecewise <- make_interval_integrator(left_i, right_i)

  left <- household_outcomes(left_i, formal, param, prices)
  right <- household_outcomes(right_i, formal, param, prices)

  formal_mass <- integrate_piecewise(
    as.numeric(formal),
    as.numeric(formal)
  )
  informal_share <- 1 - formal_mass

  formal_labor <- integrate_piecewise(
    as.numeric(formal) * param$base_formal *
      exp(param$rho_formal * left_i),
    as.numeric(formal) * param$base_formal *
      exp(param$rho_formal * right_i)
  )
  informal_labor <- integrate_piecewise(
    as.numeric(!formal) * param$base_informal *
      exp(param$rho_informal * left_i),
    as.numeric(!formal) * param$base_informal *
      exp(param$rho_informal * right_i)
  )

  voluntary_capital <- integrate_piecewise(
    left$voluntary_saving,
    right$voluntary_saving
  ) / (1 + param$n)
  funded_capital <- integrate_piecewise(
    left$funded_contribution,
    right$funded_contribution
  ) / ((1 + param$n) * (1 + param$g))
  public_saving_capital <- (
    param$universal_basic_saving +
      param$conditional_basic_saving * informal_share
  ) / ((1 + param$n) * (1 + param$g))

  list(
    micro = intervals$diagnostic_micro,
    diagnostics = formality_diagnostics(intervals$diagnostic_micro),
    informal_share = informal_share,
    formal_labor = formal_labor,
    informal_labor = informal_labor,
    capital_next = voluntary_capital + funded_capital +
      public_saving_capital,
    voluntary_capital = voluntary_capital,
    funded_capital = funded_capital,
    public_saving_capital = public_saving_capital,
    consumption_young = integrate_piecewise(
      left$consumption_young,
      right$consumption_young
    ),
    consumption_old_next = (1 - param$m) / (1 + param$n) *
      integrate_piecewise(
        left$consumption_old,
        right$consumption_old
      ),
    payg_contributions = integrate_piecewise(
      left$payg_contribution,
      right$payg_contribution
    ),
    payg_benefits_next = integrate_piecewise(
      left$payg_benefit,
      right$payg_benefit
    ) / (1 + param$n),
    decision_welfare = integrate_piecewise(
      left$decision_utility,
      right$decision_utility
    ),
    experienced_welfare = integrate_piecewise(
      left$experienced_utility,
      right$experienced_utility
    ),
    integration = list(
      nodes = intervals$nodes,
      left_i = left_i,
      right_i = right_i,
      formal = formal,
      roots = intervals$roots
    )
  )
}
