# Eleccion endogena entre informalidad, capitalizacion y PAYG.
# Se carga despues de load_validated_model.R y conserva los benchmarks previos.

default_endogenous_parameters <- function() {
  param <- default_parameters("compete")
  param$choice_model <- "endogenous"
  param$payg_eligibility_intercept <- -9.75
  param$payg_eligibility_slope <- 15.00
  param$payg_refund_rate <- 0.50
  param$payg_benefit_multiplier <- 2.80
  param$payg_minimum_benefit <- 0.05
  param$payg_maximum_benefit <- 100.00
  param
}

validate_endogenous_parameters <- function(param) {
  validate_parameters(param)
  required <- c(
    "choice_model", "payg_eligibility_intercept",
    "payg_eligibility_slope", "payg_refund_rate",
    "payg_benefit_multiplier", "payg_minimum_benefit",
    "payg_maximum_benefit"
  )
  missing <- setdiff(required, names(param))
  if (length(missing) > 0L) {
    stop_model(
      "Faltan parametros de eleccion endogena: %s.",
      paste(missing, collapse = ", ")
    )
  }
  if (!identical(param$choice_model, "endogenous")) {
    stop_model("choice_model debe ser endogenous.")
  }
  assert_scalar(
    param$payg_eligibility_intercept,
    "payg_eligibility_intercept"
  )
  assert_scalar(
    param$payg_eligibility_slope,
    "payg_eligibility_slope",
    0,
    Inf
  )
  assert_scalar(param$payg_refund_rate, "payg_refund_rate", 0, 1)
  assert_scalar(
    param$payg_benefit_multiplier,
    "payg_benefit_multiplier",
    0,
    Inf
  )
  assert_scalar(
    param$payg_minimum_benefit,
    "payg_minimum_benefit",
    0,
    Inf
  )
  assert_scalar(
    param$payg_maximum_benefit,
    "payg_maximum_benefit",
    param$payg_minimum_benefit,
    Inf
  )
  invisible(TRUE)
}

payg_eligibility_probability <- function(i, param) {
  plogis(
    param$payg_eligibility_intercept +
      param$payg_eligibility_slope * i
  )
}

pension_choice_accounts <- function(
    i,
    choice,
    wage,
    param,
    prices
) {
  n_obs <- length(i)
  if (length(choice) == 1L) {
    choice <- rep(choice, n_obs)
  }
  if (length(choice) != n_obs) {
    stop_model("i y choice deben tener la misma longitud.")
  }
  valid <- c("informal", "capitalization", "payg")
  if (any(!choice %in% valid)) {
    stop_model("Alternativa pensional no valida.")
  }

  funded <- choice == "capitalization"
  in_payg <- choice == "payg"
  formal <- funded | in_payg
  covered_wage <- wage * as.numeric(formal)
  funded_contribution <- param$tau_pension * covered_wage *
    as.numeric(funded)
  payg_contribution <- param$tau_pension * covered_wage *
    as.numeric(in_payg)

  eligibility <- payg_eligibility_probability(i, param) *
    as.numeric(in_payg)
  eligible_pension_base <- pmin(
    param$payg_maximum_benefit,
    pmax(
      param$payg_minimum_benefit,
      param$payg_benefit_multiplier *
        param$replacement_rate * covered_wage
    )
  )
  eligible_benefit <- (1 - param$m) * eligible_pension_base /
    (1 + param$g)
  refund_benefit <- param$payg_refund_rate *
    (1 + prices$r) * payg_contribution / (1 + param$g)

  list(
    funded_contribution = funded_contribution,
    payg_contribution = payg_contribution,
    funded_benefit = (1 + prices$r) * funded_contribution /
      (1 + param$g),
    payg_benefit = as.numeric(in_payg) * (
      eligibility * eligible_benefit +
        (1 - eligibility) * refund_benefit
    ),
    eligibility_probability = eligibility,
    eligible_payg_benefit = eligible_benefit * as.numeric(in_payg),
    refund_benefit = refund_benefit * as.numeric(in_payg)
  )
}

household_choice_outcomes <- function(i, choice, param, prices) {
  validate_endogenous_parameters(param)
  n_obs <- length(i)
  if (length(choice) == 1L) {
    choice <- rep(choice, n_obs)
  }
  formal <- choice != "informal"
  wage <- individual_wage(i, formal, param, prices)
  pension <- pension_choice_accounts(
    i,
    choice,
    wage,
    param,
    prices
  )

  labor_taxes <- (param$tau_labor + param$tau_pension) *
    wage * as.numeric(formal) - param$transfer_young
  lump_sum_benefit <-
    (1 - param$m) * param$universal_basic_pension +
    (1 + prices$r) * param$universal_basic_saving / (1 + param$g) +
    (1 - param$m) * param$conditional_basic_pension *
      as.numeric(!formal) +
    (1 + prices$r) * param$conditional_basic_saving *
      as.numeric(!formal) / (1 + param$g)
  benefit <- lump_sum_benefit + pension$funded_benefit +
    pension$payg_benefit

  euler_term <- (
    (1 + prices$r) * param$delta * param$beta *
      prices$p / prices$q
  )^(1 / param$sigma) * prices$q / prices$p
  savings_numerator <- euler_term * (wage - labor_taxes) -
    (1 + param$g) * benefit / (1 - param$m)
  savings_denominator <- (1 + param$g) *
    (euler_term + (1 + prices$r) / (1 - param$m))

  voluntary_saving <- pmax(0, savings_numerator / savings_denominator)
  consumption_young <- (
    wage - labor_taxes - (1 + param$g) * voluntary_saving
  ) / prices$p
  consumption_old <- (
    (1 + prices$r) * voluntary_saving + benefit
  ) / ((1 - param$m) * prices$q)
  current_utility <- crra_utility(consumption_young, param$sigma)
  scaled_old_consumption <- (1 + param$g) * consumption_old /
    (1 - param$m)
  old_utility <- crra_utility(scaled_old_consumption, param$sigma)

  data.frame(
    i = i,
    choice = choice,
    formal = formal,
    wage = wage,
    labor_taxes = labor_taxes,
    voluntary_saving = voluntary_saving,
    consumption_young = consumption_young,
    consumption_old = consumption_old,
    decision_utility = current_utility +
      (1 - param$m) * param$delta * param$beta * old_utility,
    experienced_utility = current_utility +
      (1 - param$m) * param$delta * old_utility,
    funded_contribution = pension$funded_contribution,
    payg_contribution = pension$payg_contribution,
    funded_benefit = pension$funded_benefit,
    payg_benefit = pension$payg_benefit,
    eligibility_probability = pension$eligibility_probability,
    eligible_payg_benefit = pension$eligible_payg_benefit,
    refund_benefit = pension$refund_benefit,
    stringsAsFactors = FALSE
  )
}

choose_endogenous_regime <- function(i, param, prices) {
  alternatives <- c("informal", "capitalization", "payg")
  outcomes <- lapply(
    alternatives,
    function(choice) household_choice_outcomes(
      i,
      choice,
      param,
      prices
    )
  )
  utilities <- do.call(
    cbind,
    lapply(outcomes, function(x) x$decision_utility)
  )
  if (any(!is.finite(utilities))) {
    stop_model("Alguna utilidad de las tres alternativas no es finita.")
  }

  chosen_index <- max.col(utilities, ties.method = "first")
  selected <- outcomes[[1L]]
  for (index in 2:3) {
    use <- chosen_index == index
    selected[use, ] <- outcomes[[index]][use, ]
  }
  selected$utility_informal <- utilities[, 1L]
  selected$utility_capitalization <- utilities[, 2L]
  selected$utility_payg <- utilities[, 3L]
  selected$utility_difference_formal_informal <-
    pmax(utilities[, 2L], utilities[, 3L]) - utilities[, 1L]
  selected$utility_difference_payg_capitalization <-
    utilities[, 3L] - utilities[, 2L]
  rownames(selected) <- NULL
  selected
}

endogenous_choice_diagnostics <- function(micro) {
  changes <- which(micro$choice[-1L] !=
    micro$choice[-nrow(micro)])
  sequence <- rle(micro$choice)$values
  list(
    number_of_switches = length(changes),
    switch_intervals = if (length(changes) == 0L) {
      matrix(numeric(0), ncol = 2L)
    } else {
      cbind(micro$i[changes], micro$i[changes + 1L])
    },
    choice_sequence = sequence,
    ordered_three_choice = identical(
      sequence,
      c("informal", "capitalization", "payg")
    ),
    monotone_single_crossing = length(changes) <= 2L,
    first_formal_type = if (any(micro$formal)) {
      min(micro$i[micro$formal])
    } else {
      1
    }
  )
}

endogenous_utility_difference <- function(
    i,
    left_choice,
    right_choice,
    param,
    prices
) {
  left <- household_choice_outcomes(
    i,
    left_choice,
    param,
    prices
  )$decision_utility
  right <- household_choice_outcomes(
    i,
    right_choice,
    param,
    prices
  )$decision_utility
  left - right
}

prepare_endogenous_intervals <- function(param, prices, grid) {
  diagnostic_grid <- make_type_grid_from_nodes(grid$i)
  diagnostic_micro <- choose_endogenous_regime(
    diagnostic_grid$i,
    param,
    prices
  )
  changes <- which(
    diagnostic_micro$choice[-1L] !=
      diagnostic_micro$choice[-nrow(diagnostic_micro)]
  )
  roots <- vapply(changes, function(index) {
    lower <- diagnostic_micro$i[index]
    upper <- diagnostic_micro$i[index + 1L]
    left_choice <- diagnostic_micro$choice[index]
    right_choice <- diagnostic_micro$choice[index + 1L]
    tryCatch(
      uniroot(
        function(x) endogenous_utility_difference(
          x,
          left_choice,
          right_choice,
          param,
          prices
        ),
        interval = c(lower, upper),
        tol = 1e-11
      )$root,
      error = function(e) (lower + upper) / 2
    )
  }, numeric(1))

  integration_grid <- make_type_grid_from_nodes(
    c(diagnostic_grid$i, roots)
  )
  left_i <- integration_grid$i[-length(integration_grid$i)]
  right_i <- integration_grid$i[-1L]
  midpoint_i <- (left_i + right_i) / 2
  interval_choice <- choose_endogenous_regime(
    midpoint_i,
    param,
    prices
  )$choice

  list(
    diagnostic_micro = diagnostic_micro,
    roots = roots,
    nodes = integration_grid$i,
    left_i = left_i,
    right_i = right_i,
    interval_choice = interval_choice
  )
}

evaluate_endogenous_cohort <- function(param, prices, grid) {
  intervals <- prepare_endogenous_intervals(param, prices, grid)
  left_i <- intervals$left_i
  right_i <- intervals$right_i
  choice <- intervals$interval_choice
  formal <- choice != "informal"
  integrate_piecewise <- make_interval_integrator(left_i, right_i)
  left <- household_choice_outcomes(left_i, choice, param, prices)
  right <- household_choice_outcomes(right_i, choice, param, prices)

  integrate_indicator <- function(label) {
    indicator <- as.numeric(choice == label)
    integrate_piecewise(indicator, indicator)
  }
  informal_share <- integrate_indicator("informal")
  capitalization_share <- integrate_indicator("capitalization")
  payg_share <- integrate_indicator("payg")

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
    diagnostics = endogenous_choice_diagnostics(
      intervals$diagnostic_micro
    ),
    informal_share = informal_share,
    formal_share = capitalization_share + payg_share,
    capitalization_share = capitalization_share,
    payg_share = payg_share,
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
      choice = choice,
      formal = formal,
      roots = intervals$roots
    )
  )
}

solve_steady_state_endogenous <- function(
    param = default_endogenous_parameters(),
    initial = NULL,
    grid = make_type_grid(1001L),
    old_weight = 0.97,
    tolerance = 1e-9,
    residual_tolerance = 1e-5,
    max_iterations = 6000L,
    verbose = FALSE
) {
  validate_endogenous_parameters(param)
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
    cohort <- evaluate_endogenous_cohort(iter_param, prices, grid)
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
        "La iteracion %s produjo un estado endogeno invalido.",
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
          "tauc=%.6f I=%.4f C=%.4f P=%.4f"
        ),
        iteration,
        gap,
        state[["k"]],
        state[["formal_labor"]],
        state[["consumption_tax"]],
        cohort$informal_share,
        cohort$capitalization_share,
        cohort$payg_share
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
  cohort <- evaluate_endogenous_cohort(final_param, prices, grid)
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
    informal_labor = cohort$informal_labor -
      state[["informal_labor"]],
    consumption_tax = budget$implied_tax -
      state[["consumption_tax"]],
    government_budget = budget$residual
  )
  max_residual <- max(abs(fixed_point_residuals))

  structure(
    list(
      converged = converged,
      validated = converged && max_residual < residual_tolerance,
      iterations = iteration,
      iteration_gap = gap,
      state = state,
      param = final_param,
      prices = prices,
      cohort = cohort,
      budget = budget,
      identities = identities,
      fixed_point_residuals = fixed_point_residuals,
      residual_tolerance = residual_tolerance,
      grid = grid
    ),
    class = "pension_steady_state"
  )
}
