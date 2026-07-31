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

solve_transition_endogenous <- function(
    initial_param,
    final_param,
    periods = 12L,
    initial_solution = NULL,
    final_solution = NULL,
    grid = make_type_grid(501L),
    old_weight = 0.92,
    tolerance = 1e-7,
    terminal_tolerance = 2e-3,
    budget_tolerance = 2e-5,
    cycle_tolerance = 1e-2,
    boundary_buffer = 3L,
    max_iterations = 3000L,
    verbose = FALSE
) {
  validate_endogenous_parameters(initial_param)
  validate_endogenous_parameters(final_param)
  if (initial_param$regime != "compete" ||
      final_param$regime != "pillars") {
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
  if (boundary_buffer < 2L ||
      periods < boundary_buffer + 6L) {
    stop_model(
      "Se requieren al menos seis cohortes interiores y un buffer terminal de dos."
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
    initial_solution$state[["consumption_tax"]],
    final_solution$state[["consumption_tax"]]
  )

  evaluate_paths <- function(k_path, formal_path, informal_path, tax_path) {
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
      "decision_welfare", "experienced_welfare", "budget_residual"
    )
    tracked <- setNames(
      lapply(names_to_track, function(x) rep(NA_real_, horizon)),
      names_to_track
    )
    implied_k <- k_path
    implied_formal <- formal_path
    implied_informal <- informal_path
    implied_tax <- tax_path
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
      implied_k[t + 1L] <- cohort$capital_next
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
      implied_tax[t] <- budget$implied_tax
      tracked$budget_residual[t] <- budget$residual
    }

    list(
      price_path = price_path,
      cohort_path = cohort_path,
      tracked = tracked,
      implied_k = implied_k,
      implied_formal = implied_formal,
      implied_informal = implied_informal,
      implied_tax = implied_tax
    )
  }

  converged_internal <- FALSE
  internal_gap <- Inf
  terminal_gap <- Inf
  for (iteration in seq_len(max_iterations)) {
    old_k <- k
    old_formal <- formal_labor
    old_informal <- informal_labor
    old_tax <- consumption_tax
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
    update_flow <- 2:periods
    formal_labor[update_flow] <-
      old_weight * formal_labor[update_flow] +
      (1 - old_weight) * evaluated$implied_formal[update_flow]
    informal_labor[update_flow] <-
      old_weight * informal_labor[update_flow] +
      (1 - old_weight) * evaluated$implied_informal[update_flow]
    consumption_tax[update_flow] <-
      old_weight * consumption_tax[update_flow] +
      (1 - old_weight) * evaluated$implied_tax[update_flow]

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
      evaluated$implied_k[horizon] -
        final_solution$state[["k"]]
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

  evaluated <- evaluate_paths(
    k,
    formal_labor,
    informal_labor,
    consumption_tax
  )
  terminal_gap <- abs(
    evaluated$implied_k[horizon] - final_solution$state[["k"]]
  ) / final_solution$state[["k"]]
  tracked <- evaluated$tracked
  max_budget_residual <- max(
    abs(tracked$budget_residual[2:periods]),
    na.rm = TRUE
  )
  path <- data.frame(
    cohort = 0:periods,
    capital = k,
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
    experienced_welfare = tracked$experienced_welfare,
    budget_residual = tracked$budget_residual,
    stringsAsFactors = FALSE
  )

  cycle_variables <- c(
    "capital", "formal_labor", "informal_labor", "consumption_tax"
  )
  cycle_end <- periods - boundary_buffer
  late_cohorts <- (cycle_end - 1L):cycle_end
  early_cohorts <- late_cohorts - 2L
  cycle_gap <- max(vapply(
    cycle_variables,
    function(variable) {
      max_relative_gap(
        path[path$cohort %in% late_cohorts, variable],
        path[path$cohort %in% early_cohorts, variable]
      )
    },
    numeric(1)
  ))
  phase_gap <- max(vapply(
    cycle_variables,
    function(variable) max_relative_gap(
      path[path$cohort == late_cohorts[2L], variable],
      path[path$cohort == late_cohorts[1L], variable]
    ),
    numeric(1)
  ))
  terminal_consistent <- terminal_gap < terminal_tolerance
  cycle_detected <- cycle_gap < cycle_tolerance &&
    phase_gap > 5 * cycle_tolerance
  dynamic_outcome <- if (terminal_consistent) {
    "steady_state"
  } else if (cycle_detected) {
    "two_cycle"
  } else {
    "unresolved"
  }
  structure(
    list(
      converged_internal = converged_internal,
      terminal_consistent = terminal_consistent,
      cycle_detected = cycle_detected,
      dynamic_outcome = dynamic_outcome,
      cycle_gap = cycle_gap,
      phase_gap = phase_gap,
      validated = converged_internal &&
        (terminal_consistent || cycle_detected) &&
        max_budget_residual < budget_tolerance,
      iterations = iteration,
      internal_gap = internal_gap,
      terminal_gap = terminal_gap,
      max_budget_residual = max_budget_residual,
      path = path,
      price_path = evaluated$price_path,
      cohort_path = evaluated$cohort_path,
      initial_solution = initial_solution,
      final_solution = final_solution,
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
  cat(sprintf("  Brecha terminal del capital: %.3e\n", x$terminal_gap))
  cat(sprintf("  Resultado dinamico: %s\n", x$dynamic_outcome))
  cat(sprintf(
    "  Brecha de repeticion a dos periodos: %.3e\n",
    x$cycle_gap
  ))
  cat(sprintf(
    "  Maximo residuo fiscal interior: %.3e\n",
    x$max_budget_residual
  ))
  cat(sprintf("  Validada: %s\n", ifelse(x$validated, "si", "no")))
  invisible(x)
}
