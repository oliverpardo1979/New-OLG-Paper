# Eleccion endogena entre informalidad, capitalizacion y PAYG.
# Se carga despues de load_validated_model.R y conserva los benchmarks previos.

default_endogenous_parameters <- function() {
  param <- default_parameters("compete")

  # Calibracion colombiana basada en Becerra (2025), CEDE Pension Model.
  # Un periodo del modelo equivale a 40 anos.
  param$tau_pension <- 0.160
  param$replacement_rate <- 0.650
  param$A_informal <- 0.08215

  param$choice_model <- "endogenous"
  param$payg_eligibility_intercept <- -8.79
  param$payg_eligibility_slope <- 15.00
  param$payg_refund_rate <- 1.00
  param$payg_refund_annual_return <- 0.00
  param$funded_account_rate <- 0.115
  param$payg_notional_account_rate <- 0.130
  param$payg_benefit_multiplier <- 2.90
  param$payg_minimum_benefit <- 0.05
  param$payg_maximum_benefit <- 100.00
  param$minimum_wage_model_units <- 0.50
  param$pillar_threshold_smlmv <- 2.30
  param$pillar_replacement_intercept <- 0.655
  param$pillar_replacement_slope <- 0.005
  param$solidarity_benefit_smlmv <- 218846 / 1160000
  param$solidarity_target_mass <- 0.00
  param$solidarity_type_cutoff <- 0.00

  # Momentos objetivo. La cobertura contributiva aproxima la probabilidad
  # media de elegibilidad entre formales; no es una identidad contable.
  param$target_formal_share <- 0.485
  param$target_payg_share_among_formal <- 0.100
  param$target_contributory_coverage <- 0.250
  param$target_annual_real_return <- 0.040
  param
}

validate_endogenous_parameters <- function(param) {
  validate_parameters(param)
  required <- c(
    "choice_model", "payg_eligibility_intercept",
    "payg_eligibility_slope", "payg_refund_rate",
    "payg_refund_annual_return", "funded_account_rate",
    "payg_notional_account_rate",
    "payg_benefit_multiplier", "payg_minimum_benefit",
    "payg_maximum_benefit", "minimum_wage_model_units",
    "pillar_threshold_smlmv", "pillar_replacement_intercept",
    "pillar_replacement_slope", "solidarity_benefit_smlmv",
    "solidarity_target_mass", "solidarity_type_cutoff",
    "target_formal_share",
    "target_payg_share_among_formal",
    "target_contributory_coverage",
    "target_annual_real_return"
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
    param$payg_refund_annual_return,
    "payg_refund_annual_return",
    -0.999999,
    Inf
  )
  assert_scalar(
    param$funded_account_rate,
    "funded_account_rate",
    0,
    param$tau_pension
  )
  assert_scalar(
    param$payg_notional_account_rate,
    "payg_notional_account_rate",
    0,
    param$tau_pension
  )
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
  assert_scalar(
    param$minimum_wage_model_units,
    "minimum_wage_model_units",
    .Machine$double.eps,
    Inf
  )
  assert_scalar(
    param$pillar_threshold_smlmv,
    "pillar_threshold_smlmv",
    1,
    Inf
  )
  assert_scalar(
    param$pillar_replacement_intercept,
    "pillar_replacement_intercept",
    0,
    1
  )
  assert_scalar(
    param$pillar_replacement_slope,
    "pillar_replacement_slope",
    0,
    Inf
  )
  assert_scalar(
    param$solidarity_benefit_smlmv,
    "solidarity_benefit_smlmv",
    0,
    Inf
  )
  assert_scalar(
    param$solidarity_target_mass,
    "solidarity_target_mass",
    0,
    1
  )
  assert_scalar(
    param$solidarity_type_cutoff,
    "solidarity_type_cutoff",
    0,
    1
  )
  assert_scalar(param$target_formal_share, "target_formal_share", 0, 1)
  assert_scalar(
    param$target_payg_share_among_formal,
    "target_payg_share_among_formal",
    0,
    1
  )
  assert_scalar(
    param$target_contributory_coverage,
    "target_contributory_coverage",
    0,
    param$target_formal_share
  )
  assert_scalar(
    param$target_annual_real_return,
    "target_annual_real_return",
    -0.999999,
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
  valid <- switch(
    param$regime,
    compete = c("informal", "capitalization", "payg"),
    pillars = c("informal", "pillars"),
    funded = c("informal", "capitalization")
  )
  if (any(!choice %in% valid)) {
    stop_model("Alternativa pensional no valida.")
  }

  formal <- choice != "informal"
  covered_wage <- wage * as.numeric(formal)
  if (param$regime == "compete") {
    funded_base <- covered_wage * as.numeric(
      choice == "capitalization"
    )
    payg_base <- covered_wage * as.numeric(choice == "payg")
  } else if (param$regime == "pillars") {
    funded_base <- pmax(
      0,
      covered_wage - param$payg_income_threshold
    )
    payg_base <- pmin(
      covered_wage,
      param$payg_income_threshold
    )
  } else {
    funded_base <- covered_wage
    payg_base <- rep(0, n_obs)
  }
  has_funded <- funded_base > 0
  has_payg <- payg_base > 0
  funded_contribution <- param$funded_account_rate * funded_base
  funded_public_contribution <- (
    param$tau_pension - param$funded_account_rate
  ) * funded_base
  payg_contribution <- param$tau_pension * payg_base
  payg_notional_contribution <-
    param$payg_notional_account_rate * payg_base

  eligibility <- payg_eligibility_probability(i, param) *
    as.numeric(has_payg)
  replacement_rate <- if (identical(param$regime, "pillars")) {
    wage_smlmv <- payg_base / param$minimum_wage_model_units
    pmax(
      0,
      param$pillar_replacement_intercept -
        param$pillar_replacement_slope * wage_smlmv
    )
  } else {
    rep(param$replacement_rate, n_obs)
  }
  eligible_pension_base <- pmin(
    param$payg_maximum_benefit,
    param$payg_benefit_multiplier * pmax(
      param$payg_minimum_benefit,
      replacement_rate * payg_base
    )
  )
  eligible_benefit <- (1 - param$m) * eligible_pension_base /
    (1 + param$g)
  refund_benefit <- param$payg_refund_rate *
    (1 + param$payg_refund_annual_return)^40 *
    payg_notional_contribution / (1 + param$g)

  list(
    funded_contribution = funded_contribution,
    funded_public_contribution = funded_public_contribution,
    payg_contribution = payg_contribution,
    payg_notional_contribution = payg_notional_contribution,
    funded_base = funded_base,
    payg_base = payg_base,
    funded_benefit = (1 + prices$r) * funded_contribution /
      (1 + param$g),
    payg_benefit = as.numeric(has_payg) * (
      eligibility * eligible_benefit +
        (1 - eligibility) * refund_benefit
    ),
    eligibility_probability = eligibility,
    eligible_payg_benefit = eligible_benefit *
      as.numeric(has_payg),
    refund_benefit = refund_benefit * as.numeric(has_payg),
    has_funded_component = has_funded,
    has_payg_component = has_payg
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
  solidarity_eligible <- !formal &
    i <= param$solidarity_type_cutoff
  lump_sum_benefit <-
    (1 - param$m) * param$universal_basic_pension +
    (1 + prices$r) * param$universal_basic_saving / (1 + param$g) +
    (1 - param$m) * param$conditional_basic_pension *
      as.numeric(solidarity_eligible) +
    (1 + prices$r) * param$conditional_basic_saving *
      as.numeric(solidarity_eligible) / (1 + param$g)
  benefit <- lump_sum_benefit + pension$funded_benefit +
    pension$payg_benefit
  pension_class <- ifelse(
    !formal,
    "informal",
    if (param$regime == "compete") {
      ifelse(
        choice == "capitalization",
        "funded_only",
        "payg_only"
      )
    } else if (param$regime == "pillars") {
      ifelse(
        pension$has_funded_component,
        "mixed",
        "payg_only"
      )
    } else {
      rep("funded_only", n_obs)
    }
  )

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
    funded_public_contribution = pension$funded_public_contribution,
    payg_contribution = pension$payg_contribution,
    payg_notional_contribution = pension$payg_notional_contribution,
    funded_base = pension$funded_base,
    payg_base = pension$payg_base,
    funded_benefit = pension$funded_benefit,
    payg_benefit = pension$payg_benefit,
    eligibility_probability = pension$eligibility_probability,
    eligible_payg_benefit = pension$eligible_payg_benefit,
    refund_benefit = pension$refund_benefit,
    solidarity_eligible = solidarity_eligible,
    pension_class = pension_class,
    stringsAsFactors = FALSE
  )
}

choose_endogenous_regime <- function(i, param, prices) {
  alternatives <- switch(
    param$regime,
    compete = c("informal", "capitalization", "payg"),
    pillars = c("informal", "pillars"),
    funded = c("informal", "capitalization")
  )
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
    stop_model("Alguna utilidad de las alternativas no es finita.")
  }

  chosen_index <- max.col(utilities, ties.method = "first")
  selected <- outcomes[[1L]]
  for (index in seq_along(alternatives)[-1L]) {
    use <- chosen_index == index
    selected[use, ] <- outcomes[[index]][use, ]
  }
  utility_by_choice <- setNames(
    lapply(seq_along(alternatives), function(index) {
      utilities[, index]
    }),
    paste0("utility_", alternatives)
  )
  for (name in names(utility_by_choice)) {
    selected[[name]] <- utility_by_choice[[name]]
  }
  if (param$regime == "pillars") {
    selected$utility_capitalization <- NA_real_
    selected$utility_payg <- NA_real_
  } else if (param$regime == "funded") {
    selected$utility_payg <- NA_real_
    selected$utility_pillars <- NA_real_
  }
  selected$utility_difference_formal_informal <-
    apply(utilities[, -1L, drop = FALSE], 1L, max) -
      utilities[, 1L]
  selected$utility_difference_payg_capitalization <- if (
      param$regime == "compete"
  ) {
    utilities[, 3L] - utilities[, 2L]
  } else {
    rep(NA_real_, nrow(utilities))
  }
  selected$utility_difference_pillars_informal <- if (
      param$regime == "pillars"
  ) {
    utilities[, 2L] - utilities[, 1L]
  } else {
    rep(NA_real_, nrow(utilities))
  }
  selected$utility_difference_funded_informal <- if (
      param$regime == "funded"
  ) {
    utilities[, 2L] - utilities[, 1L]
  } else {
    rep(NA_real_, nrow(utilities))
  }
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
    ordered_pillars_choice = identical(
      sequence,
      c("informal", "pillars")
    ),
    ordered_funded_choice = identical(
      sequence,
      c("informal", "capitalization")
    ),
    ordered_policy_choice = identical(
      sequence,
      c("informal", "capitalization", "payg")
    ) || identical(sequence, c("informal", "pillars")) ||
      identical(sequence, c("informal", "capitalization")),
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

  pension_split_root <- numeric(0)
  if (param$regime == "pillars" &&
      param$rho_formal > 0 &&
      param$payg_income_threshold > 0) {
    candidate <- log(
      param$payg_income_threshold /
        (prices$wage_formal * param$base_formal)
    ) / param$rho_formal
    if (is.finite(candidate) && candidate > 0 && candidate < 1) {
      pension_split_root <- candidate
    }
  }
  solidarity_root <- if (
      param$solidarity_type_cutoff > 0 &&
      param$solidarity_type_cutoff < 1
  ) {
    param$solidarity_type_cutoff
  } else numeric(0)

  integration_grid <- make_type_grid_from_nodes(
    c(diagnostic_grid$i, roots, pension_split_root,
      solidarity_root)
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
    pension_split_root = pension_split_root,
    solidarity_root = solidarity_root,
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
  midpoint_i <- (left_i + right_i) / 2
  choice <- intervals$interval_choice
  formal <- choice != "informal"
  integrate_piecewise <- make_interval_integrator(left_i, right_i)
  left <- household_choice_outcomes(left_i, choice, param, prices)
  right <- household_choice_outcomes(right_i, choice, param, prices)
  midpoint <- household_choice_outcomes(
    midpoint_i,
    choice,
    param,
    prices
  )

  integrate_choice <- function(label) {
    indicator <- as.numeric(choice == label)
    integrate_piecewise(indicator, indicator)
  }
  integrate_class <- function(label) {
    indicator <- as.numeric(midpoint$pension_class == label)
    integrate_piecewise(indicator, indicator)
  }
  informal_share <- integrate_choice("informal")
  formal_share <- integrate_piecewise(
    as.numeric(formal),
    as.numeric(formal)
  )
  capitalization_share <- integrate_choice("capitalization")
  payg_share <- integrate_choice("payg")
  funded_only_share <- integrate_class("funded_only")
  payg_only_share <- integrate_class("payg_only")
  mixed_share <- integrate_class("mixed")
  solidarity_beneficiary_share <- integrate_piecewise(
    as.numeric(midpoint$solidarity_eligible),
    as.numeric(midpoint$solidarity_eligible)
  )
  funded_covered_share <- funded_only_share + mixed_share
  payg_covered_share <- payg_only_share + mixed_share
  formal_eligibility_mass <- integrate_piecewise(
    as.numeric(formal) * payg_eligibility_probability(left_i, param),
    as.numeric(formal) * payg_eligibility_probability(right_i, param)
  )

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
      param$conditional_basic_saving * solidarity_beneficiary_share
  ) / ((1 + param$n) * (1 + param$g))

  list(
    micro = intervals$diagnostic_micro,
    diagnostics = endogenous_choice_diagnostics(
      intervals$diagnostic_micro
    ),
    informal_share = informal_share,
    formal_share = formal_share,
    capitalization_share = capitalization_share,
    payg_share = payg_share,
    funded_only_share = funded_only_share,
    payg_only_share = payg_only_share,
    mixed_share = mixed_share,
    funded_covered_share = funded_covered_share,
    solidarity_beneficiary_share = solidarity_beneficiary_share,
    payg_covered_share = payg_covered_share,
    payg_share_among_formal = switch(
      param$regime,
      compete = payg_share / formal_share,
      pillars = payg_only_share / formal_share,
      funded = 0
    ),
    formal_eligibility_probability = formal_eligibility_mass /
      formal_share,
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
    funded_public_contributions = integrate_piecewise(
      left$funded_public_contribution,
      right$funded_public_contribution
    ),
    pension_public_revenue = integrate_piecewise(
      left$funded_public_contribution + left$payg_contribution,
      right$funded_public_contribution + right$payg_contribution
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
      pension_class = midpoint$pension_class,
      formal = formal,
      roots = intervals$roots,
      pension_split_root = intervals$pension_split_root,
      solidarity_root = intervals$solidarity_root
    )
  )
}

solve_steady_state_endogenous <- function(
    param = default_endogenous_parameters(),
    initial = NULL,
    grid = make_type_grid(1001L),
    fiscal_closure = NULL,
    old_weight = 0.97,
    tolerance = 1e-9,
    residual_tolerance = 1e-5,
    max_iterations = 6000L,
    verbose = FALSE
) {
  validate_endogenous_parameters(param)
  fiscal_closure <- normalize_fiscal_steady_state_closure(
    fiscal_closure)
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
    fiscal <- steady_state_fiscal_targets(
      iter_param,
      state[["k"]],
      cohort,
      fiscal_closure
    )
    target <- c(
      k = cohort$capital_next - fiscal$domestic_debt_share *
        fiscal$public_debt,
      formal_labor = cohort$formal_labor,
      informal_labor = cohort$informal_labor,
      consumption_tax = budget$implied_tax +
        fiscal$primary_balance / budget$consumption_tax_base
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
  fiscal <- steady_state_fiscal_targets(
    final_param,
    state[["k"]],
    cohort,
    fiscal_closure
  )
  identities <- steady_state_identities(
    final_param,
    prices,
    state[["k"]],
    cohort
  )
  fixed_point_residuals <- c(
    capital = cohort$capital_next - fiscal$domestic_debt_share *
      fiscal$public_debt - state[["k"]],
    formal_labor = cohort$formal_labor - state[["formal_labor"]],
    informal_labor = cohort$informal_labor -
      state[["informal_labor"]],
    consumption_tax = budget$implied_tax +
      fiscal$primary_balance / budget$consumption_tax_base -
      state[["consumption_tax"]],
    government_budget = budget$residual - fiscal$primary_balance
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
      grid = grid,
      fiscal = fiscal
    ),
    class = "pension_steady_state"
  )
}
