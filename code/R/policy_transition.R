# Reforma pensional: transicion desde competencia de regimenes hacia pilares.
# Un periodo del modelo representa una generacion de aproximadamente 40 anos.

make_law2381_parameters <- function(
    initial_solution,
    pillar_threshold_smlmv = 2.30,
    solidarity_benefit_smlmv = 218846 / 1160000,
    solidarity_target_mass = 0.126,
    payg_benefit_multiplier = 1.00
) {
  if (!inherits(initial_solution, "pension_steady_state") ||
      !initial_solution$validated) {
    stop_model("Se requiere un estado estacionario inicial validado.")
  }
  if (initial_solution$param$regime != "compete") {
    stop_model("El escenario inicial debe representar competencia de regimenes.")
  }
  assert_scalar(
    pillar_threshold_smlmv,
    "pillar_threshold_smlmv",
    1,
    Inf
  )
  assert_scalar(
    solidarity_benefit_smlmv,
    "solidarity_benefit_smlmv",
    0,
    Inf
  )
  assert_scalar(
    solidarity_target_mass,
    "solidarity_target_mass",
    0,
    1
  )
  assert_scalar(
    payg_benefit_multiplier,
    "payg_benefit_multiplier",
    0,
    Inf
  )

  roots <- initial_solution$cohort$integration$roots
  entry_type <- if (length(roots) > 0L) {
    roots[1L]
  } else {
    initial_solution$cohort$diagnostics$first_formal_type
  }
  minimum_wage <- individual_wage(
    entry_type,
    TRUE,
    initial_solution$param,
    initial_solution$prices
  )
  solidarity_cutoff <- if (solidarity_target_mass == 0) {
    0
  } else if (solidarity_target_mass == 1) {
    1
  } else {
    uniroot(
      function(cutoff) type_cdf(cutoff) - solidarity_target_mass,
      interval = c(0, 1),
      tol = 1e-10
    )$root
  }

  param <- copy_list(initial_solution$param)
  param$regime <- "pillars"
  # Articulo 23: 13.2 puntos del excedente van a la cuenta individual.
  param$funded_account_rate <- 0.132
  param$pillar_replacement_intercept <- 0.655
  param$pillar_replacement_slope <- 0.005
  param$minimum_wage_model_units <- minimum_wage
  param$pillar_threshold_smlmv <- pillar_threshold_smlmv
  param$solidarity_benefit_smlmv <- solidarity_benefit_smlmv
  param$solidarity_target_mass <- solidarity_target_mass
  param$solidarity_type_cutoff <- solidarity_cutoff
  param$payg_income_threshold <- pillar_threshold_smlmv * minimum_wage
  param$conditional_basic_pension <-
    solidarity_benefit_smlmv * minimum_wage
  param$payg_benefit_multiplier <- payg_benefit_multiplier
  # El piso de la prestacion contributiva se normaliza a un SMLMV.
  param$payg_minimum_benefit <- minimum_wage
  validate_endogenous_parameters(param)
  param
}

realize_initial_old_endogenous <- function(
    initial_solution,
    new_param,
    realized_return,
    realized_consumption_price
) {
  old_param <- initial_solution$param
  old_prices <- initial_solution$prices
  base_nodes <- initial_solution$cohort$integration$nodes
  integration_grid <- make_type_grid_from_nodes(
    c(base_nodes, new_param$solidarity_type_cutoff)
  )
  left_i <- integration_grid$i[-length(integration_grid$i)]
  right_i <- integration_grid$i[-1L]
  midpoint_i <- (left_i + right_i) / 2
  choice <- choose_endogenous_regime(
    midpoint_i,
    old_param,
    old_prices
  )$choice
  formal <- choice != "informal"
  solidarity_interval <- !formal &
    midpoint_i <= new_param$solidarity_type_cutoff
  integrate_piecewise <- make_interval_integrator(left_i, right_i)
  old_left <- household_choice_outcomes(
    left_i,
    choice,
    old_param,
    old_prices
  )
  old_right <- household_choice_outcomes(
    right_i,
    choice,
    old_param,
    old_prices
  )
  realized_prices <- list(r = realized_return)

  realized <- function(i, old_outcome) {
    wage <- individual_wage(i, formal, old_param, old_prices)
    pension <- pension_choice_accounts(
      i,
      choice,
      wage,
      old_param,
      realized_prices
    )
    benefit <-
      pension$funded_benefit + pension$payg_benefit +
      (1 + realized_return) * old_param$universal_basic_saving /
        (1 + old_param$g) +
      (1 + realized_return) * old_param$conditional_basic_saving *
        as.numeric(!formal) / (1 + old_param$g) +
      (1 - old_param$m) * new_param$universal_basic_pension +
      (1 - old_param$m) * new_param$conditional_basic_pension *
        as.numeric(solidarity_interval)
    list(
      consumption_old = (
        (1 + realized_return) * old_outcome$voluntary_saving +
          benefit
      ) / ((1 - old_param$m) * realized_consumption_price),
      payg_benefit = pension$payg_benefit
    )
  }

  left_realized <- realized(left_i, old_left)
  right_realized <- realized(right_i, old_right)
  beneficiary_indicator <- as.numeric(solidarity_interval)
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
    informal_share = initial_solution$cohort$informal_share,
    solidarity_beneficiary_share = integrate_piecewise(
      beneficiary_indicator,
      beneficiary_indicator
    )
  )
}

tax_to_fiscal_state <- function(tax, tax_upper) {
  if (any(!is.finite(tax)) || any(tax <= 0) ||
      any(tax >= tax_upper)) {
    stop_model("El impuesto debe estar estrictamente entre cero y tax_upper.")
  }
  qlogis(tax / tax_upper)
}

fiscal_state_to_tax <- function(state, tax_upper) {
  tax_upper * plogis(state)
}

make_fiscal_rule <- function(
    final_solution,
    initial_debt = 0,
    debt_target = 0,
    tax_anchor_weight = 0.05,
    tax_upper = 0.90,
    domestic_debt_share = 0.25
) {
  if (!inherits(final_solution, "pension_steady_state") ||
      !final_solution$validated) {
    stop_model("La regla fiscal requiere un estado final validado.")
  }
  assert_scalar(initial_debt, "initial_debt", -Inf, Inf)
  assert_scalar(debt_target, "debt_target", -Inf, Inf)
  assert_scalar(tax_anchor_weight, "tax_anchor_weight", 0, 0.999999)
  assert_scalar(tax_upper, "tax_upper", 0, Inf)
  assert_scalar(domestic_debt_share, "domestic_debt_share", 0, 1)
  target_tax <- final_solution$state[["consumption_tax"]]
  if (target_tax <= 0 || target_tax >= tax_upper) {
    stop_model(
      "tax_upper debe exceder el impuesto del estado estacionario final."
    )
  }
  list(
    initial_debt = initial_debt,
    debt_target = debt_target,
    tax_anchor_weight = tax_anchor_weight,
    tax_upper = tax_upper,
    domestic_debt_share = domestic_debt_share,
    target_tax = target_tax,
    target_state = tax_to_fiscal_state(target_tax, tax_upper)
  )
}

government_debt_next <- function(
    current_debt,
    primary_deficit,
    interest_rate,
    param
) {
  growth <- (1 + param$n) * (1 + param$g)
  ((1 + interest_rate) * current_debt + primary_deficit) / growth
}

fiscal_rule_tax_current <- function(
    current_debt,
    interest_rate,
    budget,
    param,
    fiscal_rule
) {
  growth <- (1 + param$n) * (1 + param$g)
  stabilizing_tax <- (
    budget$spending - budget$non_consumption_revenue +
      (1 + interest_rate) * current_debt -
      growth * fiscal_rule$debt_target
  ) / budget$consumption_tax_base
  tax_epsilon <- 1e-8
  bounded_stabilizing_tax <- min(
    max(stabilizing_tax, tax_epsilon),
    fiscal_rule$tax_upper - tax_epsilon
  )
  desired_state <- tax_to_fiscal_state(
    bounded_stabilizing_tax,
    fiscal_rule$tax_upper
  )
  current_state <-
    fiscal_rule$tax_anchor_weight *
      fiscal_rule$target_state +
    (1 - fiscal_rule$tax_anchor_weight) * desired_state
  fiscal_state_to_tax(current_state, fiscal_rule$tax_upper)
}

solve_transition_endogenous <- function(
    initial_param,
    final_param,
    periods = 12L,
    initial_solution = NULL,
    final_solution = NULL,
    initial_path = NULL,
    grid = make_type_grid(501L),
    old_weight = 0.92,
    fiscal_old_weight = 0.95,
    tolerance = 1e-7,
    terminal_tolerance = 2e-3,
    fiscal_tolerance = 2e-5,
    fiscal_rule = NULL,
    max_iterations = 3000L,
    verbose = FALSE
) {
  validate_endogenous_parameters(initial_param)
  validate_endogenous_parameters(final_param)
  if (initial_param$regime != "compete" ||
      final_param$regime != "pillars") {
  assert_scalar(
    fiscal_old_weight,
    "fiscal_old_weight",
    0,
    0.999999
  )
    stop_model("La transicion debe ir de compete a pillars.")
  }
  if (periods < 6L) {
    stop_model("La transicion requiere al menos seis generaciones.")
  }
  fixed_demography <- c("m", "n", "g")
  changed <- fixed_demography[vapply(
    fixed_demography,
    function(name) initial_param[[name]] != final_param[[name]],
    logical(1)
  )]
  if (length(changed) > 0L) {
    stop_model(
      "La reforma mantiene fijos durante la transicion: %s.",
      paste(changed, collapse = ", ")
    )
  }

  if (is.null(initial_solution)) {
    initial_solution <- solve_steady_state_endogenous(
      initial_param,
      grid = grid
    )
  }
  if (is.null(final_solution)) {
    final_solution <- solve_steady_state_endogenous(
      final_param,
      initial = as.list(initial_solution$state),
      grid = grid,
      old_weight = 0.90
    )
  }
  if (!initial_solution$validated || !final_solution$validated) {
    stop_model("Los dos estados estacionarios deben estar validados.")
  }
  if (is.null(fiscal_rule)) {
    fiscal_rule <- make_fiscal_rule(final_solution)
  }
  required_rule <- c(
    "initial_debt", "debt_target", "tax_anchor_weight",
    "tax_upper", "domestic_debt_share", "target_tax", "target_state"
  )
  if (!all(required_rule %in% names(fiscal_rule))) {
    stop_model("La especificacion de la regla fiscal esta incompleta.")
  }
  initial_tax <- initial_solution$state[["consumption_tax"]]
  if (initial_tax <= 0 || initial_tax >= fiscal_rule$tax_upper) {
    stop_model("El impuesto inicial no pertenece al dominio de la regla.")
  }

  horizon <- periods + 1L
  interpolate_state <- function(initial, final) {
    c(initial, seq(initial, final, length.out = periods))
  }
  k <- interpolate_state(
    initial_solution$state[["k"]],
    final_solution$state[["k"]]
  )
  formal_labor <- interpolate_state(
    initial_solution$state[["formal_labor"]],
    final_solution$state[["formal_labor"]]
  )
  informal_labor <- interpolate_state(
    initial_solution$state[["informal_labor"]],
    final_solution$state[["informal_labor"]]
  )
  consumption_tax <- interpolate_state(
    initial_tax,
    final_solution$state[["consumption_tax"]]
  )

  consumption_tax[2L] <- initial_tax
  public_debt <- c(
    fiscal_rule$initial_debt,
    seq(
      fiscal_rule$initial_debt,
      fiscal_rule$debt_target,
      length.out = periods
    )
  )
  if (!is.null(initial_path)) {
    required_path <- c(
      "capital", "formal_labor", "informal_labor",
      "consumption_tax", "public_debt"
    )
    if (!is.data.frame(initial_path) ||
        nrow(initial_path) != horizon ||
        !all(required_path %in% names(initial_path))) {
      stop_model(
        "initial_path debe contener %d filas y las cinco sendas requeridas.",
        horizon
      )
    }
    k <- initial_path$capital
    formal_labor <- initial_path$formal_labor
    informal_labor <- initial_path$informal_labor
    consumption_tax <- initial_path$consumption_tax
    public_debt <- initial_path$public_debt
    if (any(!is.finite(k)) || any(k <= 0) ||
        any(!is.finite(consumption_tax)) ||
        any(consumption_tax <= 0) ||
        any(consumption_tax >= fiscal_rule$tax_upper)) {
      stop_model("initial_path queda fuera del dominio economico.")
    }
  }
  evaluate_paths <- function(
      k_path,
      formal_path,
      informal_path,
      tax_path
  ) {
    price_path <- vector("list", horizon)
    for (t in seq_len(horizon)) {
      period_param <- if (t == 1L) initial_param else final_param
      period_param <- copy_list(period_param)
      period_param$tau_consumption <- tax_path[t]
      price_path[[t]] <- factor_prices(
        k_path[t],
        formal_path[t],
        tax_path[t],
        period_param
      )
    }

    names_to_track <- c(
      "consumption_young", "consumption_old", "informal_share",
      "formal_share", "solidarity_beneficiary_share",
      "funded_only_share", "payg_only_share",
      "mixed_share", "payg_contributions", "payg_benefits",
      "decision_welfare", "experienced_welfare", "budget_residual",
      "primary_balance", "primary_deficit", "debt_identity_residual",
      "fiscal_rule_residual"
    )
    tracked <- setNames(
      lapply(names_to_track, function(x) rep(NA_real_, horizon)),
      names_to_track
    )
    implied_k <- k_path
    implied_formal <- formal_path
    implied_informal <- informal_path
    implied_tax <- tax_path
    implied_debt <- rep(fiscal_rule$debt_target, horizon)
    implied_debt[1:2] <- fiscal_rule$initial_debt
    implied_tax[1L] <- initial_tax
    gross_assets_next <- k_path
    cohort_path <- vector("list", horizon)
    cohort_path[[1L]] <- initial_solution$cohort
    cohort_path[[horizon]] <- final_solution$cohort

    set_endpoint <- function(index, solution) {
      tracked$consumption_young[index] <<-
        solution$cohort$consumption_young
      tracked$consumption_old[index] <<-
        solution$cohort$consumption_old_next
      tracked$informal_share[index] <<-
        solution$cohort$informal_share
      tracked$formal_share[index] <<- solution$cohort$formal_share
      tracked$solidarity_beneficiary_share[index] <<-
        solution$cohort$solidarity_beneficiary_share
      tracked$funded_only_share[index] <<-
        solution$cohort$funded_only_share
      tracked$payg_only_share[index] <<-
        solution$cohort$payg_only_share
      tracked$mixed_share[index] <<- solution$cohort$mixed_share
      tracked$payg_contributions[index] <<-
        solution$cohort$payg_contributions
      tracked$payg_benefits[index] <<-
        solution$cohort$payg_benefits_next
      tracked$decision_welfare[index] <<-
        solution$cohort$decision_welfare
      tracked$experienced_welfare[index] <<-
        solution$cohort$experienced_welfare
      tracked$budget_residual[index] <<- solution$budget$residual
    }
    set_endpoint(1L, initial_solution)
    set_endpoint(horizon, final_solution)

    surprised_old <- realize_initial_old_endogenous(
      initial_solution,
      final_param,
      price_path[[2L]]$r,
      price_path[[2L]]$p
    )
    tracked$consumption_old[2L] <- surprised_old$consumption_old
    tracked$payg_benefits[2L] <- surprised_old$payg_benefits
    old_solidarity_share <- rep(NA_real_, horizon)
    old_solidarity_share[1L] <-
      initial_solution$cohort$solidarity_beneficiary_share
    old_solidarity_share[2L] <-
      surprised_old$solidarity_beneficiary_share

    for (t in 2:periods) {
      decision_prices <- list(
        r = price_path[[t + 1L]]$r,
        wage_formal = price_path[[t]]$wage_formal,
        wage_informal = (1 + tax_path[t]) * final_param$A_informal,
        p = 1 + tax_path[t],
        q = 1 + tax_path[t + 1L]
      )
      period_param <- copy_list(final_param)
      period_param$tau_consumption <- tax_path[t]
      cohort <- evaluate_endogenous_cohort(
        period_param,
        decision_prices,
        grid
      )
      cohort_path[[t]] <- cohort
      gross_assets_next[t + 1L] <- cohort$capital_next
      implied_formal[t] <- cohort$formal_labor
      implied_informal[t] <- cohort$informal_labor
      tracked$consumption_young[t] <- cohort$consumption_young
      tracked$consumption_old[t + 1L] <- cohort$consumption_old_next
      tracked$informal_share[t] <- cohort$informal_share
      tracked$formal_share[t] <- cohort$formal_share
      tracked$solidarity_beneficiary_share[t] <-
        cohort$solidarity_beneficiary_share
      tracked$funded_only_share[t] <- cohort$funded_only_share
      tracked$payg_only_share[t] <- cohort$payg_only_share
      tracked$mixed_share[t] <- cohort$mixed_share
      tracked$payg_contributions[t] <- cohort$payg_contributions
      tracked$payg_benefits[t + 1L] <- cohort$payg_benefits_next
      tracked$decision_welfare[t] <- cohort$decision_welfare
      tracked$experienced_welfare[t] <- cohort$experienced_welfare
      old_solidarity_share[t + 1L] <- cohort$solidarity_beneficiary_share
    }

    for (t in 2:periods) {
      period_param <- copy_list(final_param)
      period_param$tau_consumption <- tax_path[t]
      budget <- government_budget(
        period_param,
        price_path[[t]],
        k_path[t],
        cohort_path[[t]],
        tracked$consumption_old[t],
        tracked$payg_benefits[t],
        old_solidarity_share[t]
      )
      implied_tax[t] <- fiscal_rule_tax_current(
        implied_debt[t],
        price_path[[t]]$r,
        budget,
        final_param,
        fiscal_rule
      )
      tracked$primary_balance[t] <-
        implied_tax[t] * budget$consumption_tax_base +
        budget$non_consumption_revenue - budget$spending
      tracked$primary_deficit[t] <- -tracked$primary_balance[t]
      tracked$budget_residual[t] <- tracked$primary_balance[t]
      implied_debt[t + 1L] <- government_debt_next(
        implied_debt[t],
        tracked$primary_deficit[t],
        price_path[[t]]$r,
        final_param
      )
      implied_k[t + 1L] <- gross_assets_next[t + 1L] -
        fiscal_rule$domestic_debt_share * implied_debt[t + 1L]
      if (!is.finite(implied_k[t + 1L])) {
        stop_model("El capital fisico implicito no es finito.")
      }
      growth <- (1 + final_param$n) * (1 + final_param$g)
      tracked$debt_identity_residual[t] <-
        growth * implied_debt[t + 1L] -
        ((1 + price_path[[t]]$r) * implied_debt[t] +
          tracked$primary_deficit[t])
      tracked$fiscal_rule_residual[t] <-
        tax_path[t] - implied_tax[t]
    }

    list(
      price_path = price_path,
      cohort_path = cohort_path,
      tracked = tracked,
      implied_k = implied_k,
      implied_formal = implied_formal,
      implied_informal = implied_informal,
      implied_tax = implied_tax,
      implied_debt = implied_debt,
      gross_assets_next = gross_assets_next
    )
  }

  converged_internal <- FALSE
  internal_gap <- Inf
  terminal_gap <- Inf
  terminal_gaps <- function(evaluated) {
    c(
      capital = abs(
        evaluated$implied_k[horizon] -
          final_solution$state[["k"]]
      ) / final_solution$state[["k"]],
      debt = abs(
        evaluated$implied_debt[horizon] - fiscal_rule$debt_target
      ) / final_solution$state[["k"]],
      tax = abs(
        evaluated$implied_tax[periods] - fiscal_rule$target_tax
      ) / fiscal_rule$target_tax
    )
  }

  for (iteration in seq_len(max_iterations)) {
    old_k <- k
    old_formal <- formal_labor
    old_informal <- informal_labor
    old_tax <- consumption_tax
    old_debt <- public_debt
    evaluated <- evaluate_paths(
      k,
      formal_labor,
      informal_labor,
      consumption_tax
    )

    if (periods >= 3L) {
      update_k <- 3:periods
      k[update_k] <- old_weight * k[update_k] +
        (1 - old_weight) * evaluated$implied_k[update_k]
    }
    if (any(!is.finite(k[update_k])) || any(k[update_k] <= 0)) {
      stop_model(
        "La iteracion produjo capital no positivo; ajuste la regla fiscal."
      )
    }
    update_flow <- 2:periods
    formal_labor[update_flow] <-
      old_weight * formal_labor[update_flow] +
      (1 - old_weight) * evaluated$implied_formal[update_flow]
    informal_labor[update_flow] <-
      old_weight * informal_labor[update_flow] +
      (1 - old_weight) * evaluated$implied_informal[update_flow]
    update_fiscal <- 3:periods
    public_debt[update_fiscal] <-
      fiscal_old_weight * public_debt[update_fiscal] +
      (1 - fiscal_old_weight) * evaluated$implied_debt[update_fiscal]
    update_tax <- 2:periods
    consumption_tax[update_tax] <-
      fiscal_old_weight * consumption_tax[update_tax] +
      (1 - fiscal_old_weight) * evaluated$implied_tax[update_tax]

    internal_gap_components <- c(
      capital = max_relative_gap(k[3:periods], old_k[3:periods]),
      max_relative_gap(
        formal_labor[2:periods],
        old_formal[2:periods]
      ) |> setNames("formal_labor"),
      max_relative_gap(
        informal_labor[2:periods],
        old_informal[2:periods]
      ) |> setNames("informal_labor"),
      max_relative_gap(
        consumption_tax[2:periods],
        old_tax[2:periods]
      ) |> setNames("consumption_tax"),
      public_debt = max(
        abs(public_debt[3:periods] - old_debt[3:periods])
      ) / final_solution$state[["k"]]
    )
    internal_gap <- max(internal_gap_components)
    terminal_gap <- max(terminal_gaps(evaluated))
    fiscal_gap_iteration <- max(
      abs(evaluated$tracked$fiscal_rule_residual[2:periods]),
      na.rm = TRUE
    )
    if (verbose && (iteration == 1L || iteration %% 50L == 0L)) {
      message(sprintf(
        "iter=%d gap=%.3e (%s) brecha_terminal=%.3e brecha_fiscal=%.3e",
        iteration,
        internal_gap,
        names(which.max(internal_gap_components)),
        terminal_gap,
        fiscal_gap_iteration
      ))
    }
    if (internal_gap < tolerance &&
        terminal_gap < terminal_tolerance &&
        fiscal_gap_iteration < fiscal_tolerance) {
      converged_internal <- TRUE
      break
    }
  }

  evaluated <- evaluate_paths(
    k,
    formal_labor,
    informal_labor,
    consumption_tax
  )
  terminal_gap_components <- terminal_gaps(evaluated)
  terminal_gap <- max(terminal_gap_components)
  tracked <- evaluated$tracked
  max_primary_deficit <- max(
    abs(tracked$primary_deficit[2:periods]),
    na.rm = TRUE
  )
  max_debt_identity_residual <- max(
    abs(tracked$debt_identity_residual[2:periods]),
    na.rm = TRUE
  )
  max_fiscal_rule_residual <- max(
    abs(tracked$fiscal_rule_residual[2:periods]),
    na.rm = TRUE
  )
  tracked$primary_balance[c(1L, horizon)] <- 0
  tracked$primary_deficit[c(1L, horizon)] <- 0
  tracked$debt_identity_residual[c(1L, horizon)] <- 0
  tracked$fiscal_rule_residual[c(1L, horizon)] <- 0
  formal_output <- k^final_param$alpha *
    formal_labor^(1 - final_param$alpha)
  informal_output <- final_param$A_informal * informal_labor
  total_output <- formal_output + informal_output
  path <- data.frame(
    cohort = 0:periods,
    capital = k,
    public_debt = evaluated$implied_debt,
    output = total_output,
    debt_to_output = evaluated$implied_debt / total_output,
    formal_labor = formal_labor,
    informal_labor = informal_labor,
    consumption_tax = consumption_tax,
    interest_rate = vapply(
      evaluated$price_path,
      `[[`,
      numeric(1),
      "r"
    ),
    formal_wage = vapply(
      evaluated$price_path,
      `[[`,
      numeric(1),
      "wage_formal"
    ),
    informal_wage = vapply(
      evaluated$price_path,
      `[[`,
      numeric(1),
      "wage_informal"
    ),
    consumption_young = tracked$consumption_young,
    consumption_old = tracked$consumption_old,
    informal_share = tracked$informal_share,
    formal_share = tracked$formal_share,
    solidarity_beneficiary_share = tracked$solidarity_beneficiary_share,
    funded_only_share = tracked$funded_only_share,
    payg_only_share = tracked$payg_only_share,
    mixed_share = tracked$mixed_share,
    payg_contributions = tracked$payg_contributions,
    payg_benefits = tracked$payg_benefits,
    decision_welfare = tracked$decision_welfare,
    primary_balance = tracked$primary_balance,
    primary_deficit = tracked$primary_deficit,
    debt_identity_residual = tracked$debt_identity_residual,
    fiscal_rule_residual = tracked$fiscal_rule_residual,
    experienced_welfare = tracked$experienced_welfare,
    budget_residual = tracked$budget_residual,
    stringsAsFactors = FALSE
  )

  terminal_consistent <- terminal_gap < terminal_tolerance
  fiscal_consistent <-
    max_debt_identity_residual < fiscal_tolerance &&
    max_fiscal_rule_residual < fiscal_tolerance
  dynamic_outcome <- if (
      converged_internal && terminal_consistent && fiscal_consistent
  ) {
    "steady_state"
  } else {
    "unresolved"
  }
  structure(
    list(
      converged_internal = converged_internal,
      terminal_consistent = terminal_consistent,
      fiscal_consistent = fiscal_consistent,
      dynamic_outcome = dynamic_outcome,
      validated = converged_internal &&
        terminal_consistent && fiscal_consistent,
      iterations = iteration,
      internal_gap = internal_gap,
      terminal_gap = terminal_gap,
      terminal_gap_components = terminal_gap_components,
      max_primary_deficit = max_primary_deficit,
      max_debt_identity_residual = max_debt_identity_residual,
      max_fiscal_rule_residual = max_fiscal_rule_residual,
      path = path,
      price_path = evaluated$price_path,
      cohort_path = evaluated$cohort_path,
      initial_solution = initial_solution,
      final_solution = final_solution,
      internal_gap_components = internal_gap_components,
      fiscal_rule = fiscal_rule,
      grid = grid,
      timing = paste(
        "La reforma comienza con la cohorte 1.",
        "La cohorte vieja conserva derechos contributivos de Ley 100.",
        "Cada periodo representa cerca de 40 anos."
      )
    ),
    class = "pension_policy_transition"
  )
}

weighted_quantile <- function(x, weights, probabilities) {
  valid <- is.finite(x) & is.finite(weights) & weights >= 0
  x <- x[valid]
  weights <- weights[valid]
  if (length(x) == 0L || sum(weights) <= 0) {
    return(rep(NA_real_, length(probabilities)))
  }
  order_x <- order(x)
  x <- x[order_x]
  weights <- weights[order_x] / sum(weights)
  cumulative <- cumsum(weights)
  vapply(probabilities, function(probability) {
    x[which(cumulative >= probability)[1L]]
  }, numeric(1))
}

weighted_gini <- function(x, weights) {
  valid <- is.finite(x) & is.finite(weights) & weights >= 0 & x >= 0
  x <- x[valid]
  weights <- weights[valid]
  if (length(x) == 0L || sum(weights) <= 0 ||
      sum(x * weights) <= 0) {
    return(NA_real_)
  }
  order_x <- order(x)
  x <- x[order_x]
  weights <- weights[order_x] / sum(weights)
  income_share <- x * weights / sum(x * weights)
  cumulative_income <- cumsum(income_share)
  previous_income <- c(0, cumulative_income[-length(cumulative_income)])
  1 - sum(weights * (cumulative_income + previous_income))
}

experienced_consumption_equivalent <- function(
    policy_utility,
    baseline_utility,
    param,
    old_only = FALSE
) {
  valid <- is.finite(policy_utility) & is.finite(baseline_utility)
  result <- rep(NA_real_, length(policy_utility))
  if (param$sigma == 1) {
    scale_weight <- if (old_only) {
      1
    } else {
      1 + (1 - param$m) * param$delta
    }
    result[valid] <- exp(
      (policy_utility[valid] - baseline_utility[valid]) /
        scale_weight
    ) - 1
  } else {
    ratio <- policy_utility[valid] / baseline_utility[valid]
    valid_ratio <- is.finite(ratio) & ratio > 0
    converted <- rep(NA_real_, length(ratio))
    converted[valid_ratio] <-
      ratio[valid_ratio]^(1 / (1 - param$sigma)) - 1
    result[valid] <- converted
  }
  result
}

build_policy_microdata <- function(transition) {
  if (!inherits(transition, "pension_policy_transition")) {
    stop_model("Se requiere una transicion pensional endogena.")
  }
  grid <- transition$grid
  initial <- transition$initial_solution
  final <- transition$final_solution
  baseline <- choose_endogenous_regime(
    grid$i,
    initial$param,
    initial$prices
  )

  make_young_cohort <- function(cohort) {
    if (cohort == 0L) {
      outcome <- baseline
      policy_label <- "Law 100"
    } else if (cohort == max(transition$path$cohort)) {
      outcome <- choose_endogenous_regime(
        grid$i,
        final$param,
        final$prices
      )
      policy_label <- "Law 2381 pillars"
    } else {
      t <- cohort + 1L
      tax_path <- transition$path$consumption_tax
      period_param <- copy_list(final$param)
      period_param$tau_consumption <- tax_path[t]
      decision_prices <- list(
        r = transition$price_path[[t + 1L]]$r,
        wage_formal = transition$price_path[[t]]$wage_formal,
        wage_informal = (1 + tax_path[t]) * final$param$A_informal,
        p = 1 + tax_path[t],
        q = 1 + tax_path[t + 1L]
      )
      outcome <- choose_endogenous_regime(
        grid$i,
        period_param,
        decision_prices
      )
      policy_label <- "Law 2381 pillars"
    }
    cev <- experienced_consumption_equivalent(
      outcome$experienced_utility,
      baseline$experienced_utility,
      initial$param
    )
    data.frame(
      cohort = cohort,
      life_stage = "young_lifetime",
      policy = policy_label,
      i = grid$i,
      weight = grid$weights,
      choice = outcome$choice,
      pension_class = outcome$pension_class,
      formal = outcome$formal,
      solidarity_eligible = outcome$solidarity_eligible,
      wage = outcome$wage,
      consumption_young = outcome$consumption_young,
      consumption_old = outcome$consumption_old,
      experienced_utility = outcome$experienced_utility,
      baseline_utility = baseline$experienced_utility,
      welfare_cev = cev,
      stringsAsFactors = FALSE
    )
  }

  young <- do.call(rbind, lapply(
    0:max(transition$path$cohort),
    make_young_cohort
  ))

  realized_return <- transition$price_path[[2L]]$r
  realized_price <- transition$price_path[[2L]]$p
  old_pension <- pension_choice_accounts(
    grid$i,
    baseline$choice,
    baseline$wage,
    initial$param,
    list(r = realized_return)
  )
  old_solidarity_eligible <- !baseline$formal &
    grid$i <= final$param$solidarity_type_cutoff
  old_benefit <-
    old_pension$funded_benefit + old_pension$payg_benefit +
    (1 + realized_return) * initial$param$universal_basic_saving /
      (1 + initial$param$g) +
    (1 + realized_return) * initial$param$conditional_basic_saving *
      as.numeric(!baseline$formal) / (1 + initial$param$g) +
    (1 - initial$param$m) * final$param$universal_basic_pension +
    (1 - initial$param$m) * final$param$conditional_basic_pension *
      as.numeric(old_solidarity_eligible)
  old_consumption <- (
    (1 + realized_return) * baseline$voluntary_saving + old_benefit
  ) / ((1 - initial$param$m) * realized_price)
  old_scaled <- (1 + initial$param$g) * old_consumption /
    (1 - initial$param$m)
  baseline_old_scaled <- (1 + initial$param$g) *
    baseline$consumption_old / (1 - initial$param$m)
  old_utility <- crra_utility(old_scaled, initial$param$sigma)
  baseline_old_utility <- crra_utility(
    baseline_old_scaled,
    initial$param$sigma
  )
  old <- data.frame(
    cohort = -1L,
    life_stage = "old_at_reform",
    policy = "Law 2381 pillars",
    i = grid$i,
    weight = grid$weights,
    choice = baseline$choice,
    pension_class = baseline$pension_class,
    formal = baseline$formal,
    wage = baseline$wage,
    solidarity_eligible = old_solidarity_eligible,
    consumption_young = NA_real_,
    consumption_old = old_consumption,
    experienced_utility = old_utility,
    baseline_utility = baseline_old_utility,
    welfare_cev = experienced_consumption_equivalent(
      old_utility,
      baseline_old_utility,
      initial$param,
      old_only = TRUE
    ),
    stringsAsFactors = FALSE
  )
  rownames(young) <- NULL
  rbind(old, young)
}

summarize_policy_distributions <- function(micro, param) {
  groups <- split(
    micro,
    interaction(micro$cohort, micro$life_stage, drop = TRUE)
  )
  summary <- lapply(groups, function(group) {
    weights <- group$weight / sum(group$weight)
    wage_quantiles <- weighted_quantile(
      group$wage,
      weights,
      c(0.10, 0.50, 0.90)
    )
    welfare_quantiles <- weighted_quantile(
      group$welfare_cev,
      weights,
      c(0.10, 0.50, 0.90)
    )
    aggregate_cev <- experienced_consumption_equivalent(
      weighted_sum(group$experienced_utility, weights),
      weighted_sum(group$baseline_utility, weights),
      param,
      old_only = group$life_stage[1L] == "old_at_reform"
    )
    data.frame(
      cohort = group$cohort[1L],
      life_stage = group$life_stage[1L],
      policy = group$policy[1L],
      mean_wage = weighted_sum(group$wage, weights),
      wage_p10 = wage_quantiles[1L],
      wage_p50 = wage_quantiles[2L],
      wage_p90 = wage_quantiles[3L],
      wage_gini = weighted_gini(group$wage, weights),
      aggregate_welfare_cev = aggregate_cev,
      mean_individual_welfare_cev = weighted_sum(
        group$welfare_cev,
        weights
      ),
      welfare_p10 = welfare_quantiles[1L],
      welfare_p50 = welfare_quantiles[2L],
      welfare_p90 = welfare_quantiles[3L],
      winner_share = sum(weights[group$welfare_cev > 1e-8]),
      loser_share = sum(weights[group$welfare_cev < -1e-8]),
      informal_share = sum(weights[group$pension_class == "informal"]),
      funded_only_share = sum(
        weights[group$pension_class == "funded_only"]
      ),
      payg_only_share = sum(
        weights[group$pension_class == "payg_only"]
      ),
      mixed_share = sum(weights[group$pension_class == "mixed"]),
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, summary)
  result <- result[order(result$cohort), ]
  rownames(result) <- NULL
  result
}

print.pension_policy_transition <- function(x, ...) {
  cat("Transicion Ley 100 a sistema de pilares\n")
  cat(sprintf(
    "  Punto fijo interior: %s (%d iteraciones)\n",
    ifelse(x$converged_internal, "si", "no"),
    x$iterations
  ))
  cat(sprintf("  Brecha terminal conjunta: %.3e\n", x$terminal_gap))
  cat(sprintf("  Resultado dinamico: %s\n", x$dynamic_outcome))
  cat(sprintf("  Identidades fiscales: %s\n",
              ifelse(x$fiscal_consistent, "si", "no")))
  cat(sprintf(
    "  Maximo residuo de deuda: %.3e\n",
    x$max_debt_identity_residual
  ))
  cat(sprintf(
    "  Maximo residuo de regla fiscal: %.3e\n",
    x$max_fiscal_rule_residual
  ))
  cat(sprintf("  Validada: %s\n", ifelse(x$validated, "si", "no")))
  invisible(x)
}
