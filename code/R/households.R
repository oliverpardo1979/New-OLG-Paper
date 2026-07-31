factor_prices <- function(k, formal_labor, consumption_tax, param) {
  validate_parameters(param)
  assert_scalar(k, "k", .Machine$double.eps, Inf)
  assert_scalar(formal_labor, "formal_labor", .Machine$double.eps, Inf)
  assert_scalar(consumption_tax, "consumption_tax", -0.999999, Inf)

  mpk_net <- param$alpha * (formal_labor / k)^(1 - param$alpha) -
    param$depreciation

  list(
    r = (1 - param$tau_capital) * mpk_net,
    wage_formal = ((1 - param$alpha) *
      (k / formal_labor)^param$alpha) / (1 + param$tau_payroll),
    wage_informal = (1 + consumption_tax) * param$A_informal,
    p = 1 + consumption_tax,
    q = 1 + consumption_tax
  )
}

individual_wage <- function(i, formal, param, prices) {
  ifelse(
    formal,
    prices$wage_formal * param$base_formal * exp(param$rho_formal * i),
    prices$wage_informal * param$base_informal *
      exp(param$rho_informal * i)
  )
}

pension_accounts <- function(i, formal, wage, param, prices) {
  covered_wage <- wage * as.numeric(formal)

  if (param$regime == "compete") {
    in_payg <- i >= param$payg_type_threshold
    funded_base <- covered_wage * (!in_payg)
    payg_base <- covered_wage * in_payg
  } else {
    funded_base <- pmax(0, covered_wage - param$payg_income_threshold)
    payg_base <- pmin(covered_wage, param$payg_income_threshold)
  }

  funded_contribution <- param$tau_pension * funded_base
  payg_contribution <- param$tau_pension * payg_base

  list(
    funded_contribution = funded_contribution,
    payg_contribution = payg_contribution,
    funded_benefit = (1 + prices$r) * funded_contribution / (1 + param$g),
    payg_benefit = (1 - param$m) * param$replacement_rate *
      payg_base / (1 + param$g)
  )
}

household_outcomes <- function(i, formal, param, prices) {
  validate_parameters(param)
  n_obs <- length(i)
  if (length(formal) == 1L) {
    formal <- rep(formal, n_obs)
  }
  if (length(formal) != n_obs) {
    stop_model("`i` y `formal` deben tener la misma longitud.")
  }

  wage <- individual_wage(i, formal, param, prices)
  pension <- pension_accounts(i, formal, wage, param, prices)
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

  decision_utility <- current_utility +
    (1 - param$m) * param$delta * param$beta * old_utility
  experienced_utility <- current_utility +
    (1 - param$m) * param$delta * old_utility

  data.frame(
    i = i,
    formal = as.logical(formal),
    wage = wage,
    labor_taxes = labor_taxes,
    voluntary_saving = voluntary_saving,
    consumption_young = consumption_young,
    consumption_old = consumption_old,
    decision_utility = decision_utility,
    experienced_utility = experienced_utility,
    funded_contribution = pension$funded_contribution,
    payg_contribution = pension$payg_contribution,
    funded_benefit = pension$funded_benefit,
    payg_benefit = pension$payg_benefit
  )
}

choose_sector <- function(i, param, prices) {
  informal <- household_outcomes(i, FALSE, param, prices)
  formal <- household_outcomes(i, TRUE, param, prices)
  difference <- formal$decision_utility - informal$decision_utility

  if (any(!is.finite(difference))) {
    stop_model(
      paste(
        "La utilidad formal menos la informal no es finita.",
        "Revise positividad del consumo y parámetros."
      )
    )
  }

  choose_formal <- difference >= 0
  selected <- informal
  selected[choose_formal, ] <- formal[choose_formal, ]
  selected$utility_difference <- difference
  selected
}

formality_diagnostics <- function(micro) {
  formal <- as.integer(micro$formal)
  changes <- which(diff(formal) != 0)

  list(
    number_of_switches = length(changes),
    switch_intervals = if (length(changes) == 0L) {
      matrix(numeric(0), ncol = 2L)
    } else {
      cbind(micro$i[changes], micro$i[changes + 1L])
    },
    monotone_single_crossing = length(changes) <= 1L &&
      !any(diff(formal) < 0),
    first_formal_type = if (any(micro$formal)) {
      min(micro$i[micro$formal])
    } else {
      1
    }
  )
}
