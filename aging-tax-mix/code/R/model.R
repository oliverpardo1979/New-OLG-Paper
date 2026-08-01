# Core equations for "Financing Population Aging".
#
# The paper studies conditional stationary equilibria.  Each observation in the
# demographic exercise is a steady state associated with that year's DANE
# old-age dependency ratio.  It is not a perfect-foresight transition.

logistic <- function(x) {
  1 / (1 + exp(-pmax(pmin(x, 35), -35)))
}

crra_utility <- function(c, sigma = 2) {
  if (!is.finite(c) || c <= 0) return(-Inf)
  if (abs(sigma - 1) < 1e-10) return(log(c))
  (c^(1 - sigma) - 1) / (1 - sigma)
}

ces_quantity <- function(c_formal, c_informal, omega, eta) {
  if (c_formal <= 0 || c_informal <= 0) return(NA_real_)
  if (abs(eta - 1) < 1e-8) {
    return(c_formal^omega * c_informal^(1 - omega))
  }
  rho <- (eta - 1) / eta
  (omega^(1 / eta) * c_formal^rho +
     (1 - omega)^(1 / eta) * c_informal^rho)^(1 / rho)
}

ces_price_index <- function(price_formal, price_informal, omega, eta) {
  if (abs(eta - 1) < 1e-8) {
    return((price_formal / omega)^omega *
             (price_informal / (1 - omega))^(1 - omega))
  }
  (omega * price_formal^(1 - eta) +
     (1 - omega) * price_informal^(1 - eta))^(1 / (1 - eta))
}

default_parameters <- function() {
  list(
    alpha = 0.33,
    delta = 0.06,
    beta = 1 / 1.04,
    sigma = 2,
    A_formal = 1,
    A_informal = 0.60,
    omega_formal = 0.65,
    eta = 1.50,
    informal_vat_coverage = 0.15,
    choice_scale = 0.30,
    formal_cost = 0,
    old_age_outlay = 0.40,
    resource_share = 0.25,
    bond_rate = 0.04,
    output_growth = 0.01,
    debt_anchor = 0.55,
    fiscal_feedback = 0.10,
    public_good_weight = 0,
    collection_cost_consumption = 2,
    collection_cost_payroll = 2,
    collection_cost_capital = 2,
    max_tax_consumption = 0.45,
    max_tax_payroll = 0.40,
    max_tax_capital = 0.50
  )
}

validate_mix <- function(mix) {
  mix <- as.numeric(mix)
  if (length(mix) != 3 || any(!is.finite(mix)) || any(mix < -1e-10)) {
    stop("Tax mix must contain three non-negative finite shares.")
  }
  if (abs(sum(mix) - 1) > 1e-8) {
    stop("Tax mix shares must sum to one.")
  }
  names(mix) <- c("consumption", "payroll", "capital")
  mix
}

steady_components <- function(z, dependency_ratio, mix, par,
                              debt_ratio = par$debt_anchor) {
  mix <- validate_mix(mix)
  capital <- exp(z[1])
  labor_formal <- logistic(z[2])
  labor_informal <- 1 - labor_formal
  price_informal <- exp(z[3])

  output_formal <- par$A_formal * capital^par$alpha *
    labor_formal^(1 - par$alpha)
  output_informal <- par$A_informal * labor_informal
  wage_formal <- (1 - par$alpha) * output_formal / labor_formal
  rental_rate <- par$alpha * output_formal / capital
  wage_informal <- price_informal * par$A_informal
  nominal_output <- output_formal + price_informal * output_informal

  social_security <- par$old_age_outlay * dependency_ratio
  public_absorption <- par$resource_share * social_security
  consumption_formal <- exp(z[4])
  consumption_informal <- output_informal

  primary_surplus_ratio <- (par$bond_rate - par$output_growth) *
    debt_ratio + par$fiscal_feedback * (debt_ratio - par$debt_anchor)
  required_revenue <- social_security +
    primary_surplus_ratio * nominal_output

  base_consumption <- consumption_formal +
    par$informal_vat_coverage * price_informal * consumption_informal
  base_payroll <- wage_formal * labor_formal
  base_capital <- rental_rate * capital

  tax_consumption <- mix["consumption"] * required_revenue /
    pmax(base_consumption, 1e-12)
  tax_payroll <- mix["payroll"] * required_revenue /
    pmax(base_payroll, 1e-12)
  tax_capital <- mix["capital"] * required_revenue /
    pmax(base_capital, 1e-12)
  collection_cost <- 0.5 * (
    par$collection_cost_consumption * tax_consumption^2 *
      base_consumption +
      par$collection_cost_payroll * tax_payroll^2 * base_payroll +
      par$collection_cost_capital * tax_capital^2 * base_capital
  )

  consumer_price_formal <- 1 + tax_consumption
  consumer_price_informal <- price_informal *
    (1 + par$informal_vat_coverage * tax_consumption)
  price_index <- ces_price_index(
    consumer_price_formal,
    consumer_price_informal,
    par$omega_formal,
    par$eta
  )
  composite_consumption <- ces_quantity(
    consumption_formal,
    consumption_informal,
    par$omega_formal,
    par$eta
  )

  net_formal_wage <- (1 - tax_payroll) * wage_formal
  formal_choice_index <- (
    log(pmax(net_formal_wage, 1e-12)) -
      log(pmax(wage_informal, 1e-12)) -
      par$formal_cost
  ) / par$choice_scale

  demand_log_ratio <- log(par$omega_formal / (1 - par$omega_formal)) -
    par$eta * log(consumer_price_formal / consumer_price_informal)
  observed_log_ratio <- log(pmax(consumption_formal, 1e-12) /
                              pmax(consumption_informal, 1e-12))
  gross_after_tax_return <- 1 - par$delta +
    (1 - tax_capital) * rental_rate

  residuals <- c(
    labor_choice = qlogis(labor_formal) - formal_choice_index,
    ces_demand = observed_log_ratio - demand_log_ratio,
    capital_euler = log(pmax(par$beta * gross_after_tax_return, 1e-12)),
    formal_goods = (
      consumption_formal + par$delta * capital +
        public_absorption + collection_cost - output_formal
    ) / pmax(output_formal, 1e-12)
  )

  utility_private <- crra_utility(composite_consumption, par$sigma)
  utility_public <- if (par$public_good_weight > 0 &&
                        public_absorption > 0) {
    par$public_good_weight * log(public_absorption)
  } else {
    0
  }

  list(
    z = z,
    residuals = residuals,
    capital = capital,
    labor_formal = labor_formal,
    labor_informal = labor_informal,
    informality = labor_informal,
    price_informal = price_informal,
    output_formal = output_formal,
    output_informal = output_informal,
    nominal_output = nominal_output,
    wage_formal = wage_formal,
    wage_informal = wage_informal,
    rental_rate = rental_rate,
    dependency_ratio = dependency_ratio,
    social_security = social_security,
    social_security_share = social_security / nominal_output,
    public_absorption = public_absorption,
    collection_cost = unname(collection_cost),
    collection_cost_share = unname(collection_cost / nominal_output),
    consumption_formal = consumption_formal,
    consumption_informal = consumption_informal,
    composite_consumption = composite_consumption,
    consumption_price_index = price_index,
    primary_surplus_ratio = primary_surplus_ratio,
    required_revenue = required_revenue,
    revenue_share = required_revenue / nominal_output,
    tax_consumption = unname(tax_consumption),
    tax_payroll = unname(tax_payroll),
    tax_capital = unname(tax_capital),
    mix_consumption = unname(mix["consumption"]),
    mix_payroll = unname(mix["payroll"]),
    mix_capital = unname(mix["capital"]),
    gross_after_tax_return = gross_after_tax_return,
    flow_utility = utility_private + utility_public,
    lifetime_welfare = (utility_private + utility_public) / (1 - par$beta)
  )
}

equilibrium_penalty <- function(x, par) {
  penalty <- 0
  lower_positive <- c(
    x$consumption_formal,
    x$consumption_informal,
    x$composite_consumption,
    x$nominal_output,
    1 - x$tax_payroll,
    x$gross_after_tax_return
  )
  if (any(!is.finite(lower_positive))) return(1e8)
  if (any(lower_positive <= 0)) {
    penalty <- penalty + 1e5 +
      1e5 * sum(pmax(-lower_positive, 0)^2)
  }
  tax_rates <- c(
    x$tax_consumption,
    x$tax_payroll,
    x$tax_capital
  )
  tax_caps <- c(
    par$max_tax_consumption,
    par$max_tax_payroll,
    par$max_tax_capital
  )
  if (any(!is.finite(tax_rates)) || any(tax_rates < 0)) return(1e8)
  penalty + 1e4 * sum(pmax(tax_rates - tax_caps, 0)^2)
}
