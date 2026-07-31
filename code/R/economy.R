evaluate_cohort <- function(param, prices, grid) {
  micro <- choose_sector(grid$i, param, prices)
  weights <- grid$weights
  formal <- as.numeric(micro$formal)
  informal <- 1 - formal

  effective_formal_labor <- param$base_formal *
    exp(param$rho_formal * grid$i)
  effective_informal_labor <- param$base_informal *
    exp(param$rho_informal * grid$i)

  informal_share <- weighted_sum(informal, weights)
  formal_labor <- weighted_sum(
    formal * effective_formal_labor,
    weights
  )
  informal_labor <- weighted_sum(
    informal * effective_informal_labor,
    weights
  )

  voluntary_capital <- weighted_sum(
    micro$voluntary_saving,
    weights
  ) / (1 + param$n)
  funded_capital <- weighted_sum(
    micro$funded_contribution,
    weights
  ) / ((1 + param$n) * (1 + param$g))
  public_saving_capital <- (
    param$universal_basic_saving +
      param$conditional_basic_saving * informal_share
  ) / ((1 + param$n) * (1 + param$g))

  list(
    micro = micro,
    diagnostics = formality_diagnostics(micro),
    informal_share = informal_share,
    formal_labor = formal_labor,
    informal_labor = informal_labor,
    capital_next = voluntary_capital + funded_capital +
      public_saving_capital,
    voluntary_capital = voluntary_capital,
    funded_capital = funded_capital,
    public_saving_capital = public_saving_capital,
    consumption_young = weighted_sum(
      micro$consumption_young,
      weights
    ),
    consumption_old_next = (1 - param$m) / (1 + param$n) *
      weighted_sum(micro$consumption_old, weights),
    payg_contributions = weighted_sum(
      micro$payg_contribution,
      weights
    ),
    payg_benefits_next = weighted_sum(
      micro$payg_benefit,
      weights
    ) / (1 + param$n),
    decision_welfare = weighted_sum(
      micro$decision_utility,
      weights
    ),
    experienced_welfare = weighted_sum(
      micro$experienced_utility,
      weights
    )
  )
}

government_budget <- function(
    param,
    prices,
    k,
    cohort,
    consumption_old,
    payg_benefits,
    old_informal_share = cohort$informal_share
) {
  capital_tax_base <- (
    param$alpha * (cohort$formal_labor / k)^(1 - param$alpha) -
      param$depreciation
  ) * k

  spending <-
    param$government_spending +
    param$universal_basic_saving +
    param$conditional_basic_saving * cohort$informal_share +
    (1 - param$m) / (1 + param$n) * (
      param$universal_basic_pension +
        param$conditional_basic_pension * old_informal_share
    ) +
    payg_benefits +
    param$transfer_young

  pension_public_revenue <- if (!is.null(cohort$pension_public_revenue)) {
    cohort$pension_public_revenue
  } else {
    cohort$payg_contributions
  }

  non_consumption_revenue <-
    (param$tau_labor + param$tau_payroll) *
      prices$wage_formal * cohort$formal_labor +
    pension_public_revenue +
    param$tau_capital * capital_tax_base

  consumption_tax_base <- cohort$consumption_young +
    consumption_old - param$A_informal * cohort$informal_labor

  if (!is.finite(consumption_tax_base) ||
      consumption_tax_base <= 0) {
    stop_model("La base del impuesto al consumo no es positiva.")
  }

  implied_tax <- (spending - non_consumption_revenue) /
    consumption_tax_base
  actual_revenue <- param$tau_consumption * consumption_tax_base +
    non_consumption_revenue

  list(
    spending = spending,
    non_consumption_revenue = non_consumption_revenue,
    consumption_tax_base = consumption_tax_base,
    implied_tax = implied_tax,
    residual = actual_revenue - spending
  )
}

steady_state_identities <- function(param, prices, k, cohort) {
  output_formal <- k^param$alpha *
    cohort$formal_labor^(1 - param$alpha)
  output_informal <- param$A_informal * cohort$informal_labor
  investment <- (
    (1 + param$g) * (1 + param$n) -
      (1 - param$depreciation)
  ) * k

  supply_gdp <- output_formal + output_informal
  demand_gdp <- investment + cohort$consumption_young +
    cohort$consumption_old_next + param$government_spending
  income_gdp <- (1 + param$tau_payroll) *
    prices$wage_formal * cohort$formal_labor +
    (
      prices$r / (1 - param$tau_capital) +
        param$depreciation
    ) * k +
    output_informal

  c(
    supply_gdp = supply_gdp,
    demand_gdp = demand_gdp,
    income_gdp = income_gdp,
    supply_minus_demand = supply_gdp - demand_gdp,
    supply_minus_income = supply_gdp - income_gdp
  )
}
