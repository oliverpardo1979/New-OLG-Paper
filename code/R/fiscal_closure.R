make_fiscal_steady_state_closure <- function(
    debt_to_annual_gdp = 0,
    primary_balance_to_gdp = 0,
    years_per_period = 40,
    domestic_debt_share = 0
) {
  assert_scalar(
    debt_to_annual_gdp,
    "debt_to_annual_gdp",
    0,
    Inf
  )
  assert_scalar(
    primary_balance_to_gdp,
    "primary_balance_to_gdp",
    -Inf,
    Inf
  )
  assert_scalar(years_per_period, "years_per_period", 1, Inf)
  assert_scalar(
    domestic_debt_share,
    "domestic_debt_share",
    0,
    1
  )
  if (debt_to_annual_gdp == 0 && primary_balance_to_gdp != 0) {
    stop_model(
      "Un estado sin deuda requiere balance primario objetivo igual a cero."
    )
  }
  list(
    debt_to_annual_gdp = debt_to_annual_gdp,
    primary_balance_to_gdp = primary_balance_to_gdp,
    years_per_period = years_per_period,
    domestic_debt_share = domestic_debt_share
  )
}

normalize_fiscal_steady_state_closure <- function(fiscal_closure) {
  if (is.null(fiscal_closure)) {
    return(make_fiscal_steady_state_closure())
  }
  required <- c(
    "debt_to_annual_gdp", "primary_balance_to_gdp",
    "years_per_period", "domestic_debt_share"
  )
  if (!is.list(fiscal_closure) ||
      !all(required %in% names(fiscal_closure))) {
    stop_model("La especificacion fiscal estacionaria esta incompleta.")
  }
  do.call(
    make_fiscal_steady_state_closure,
    fiscal_closure[required]
  )
}

steady_state_fiscal_targets <- function(
    param,
    k,
    cohort,
    fiscal_closure = NULL
) {
  closure <- normalize_fiscal_steady_state_closure(fiscal_closure)
  output <- k^param$alpha *
    cohort$formal_labor^(1 - param$alpha) +
    param$A_informal * cohort$informal_labor
  public_debt <- closure$debt_to_annual_gdp * output /
    closure$years_per_period
  primary_balance <- closure$primary_balance_to_gdp * output
  growth <- (1 + param$n) * (1 + param$g)
  gross_debt_return <- if (public_debt > 0) {
    growth + primary_balance / public_debt
  } else {
    growth
  }
  if (!is.finite(gross_debt_return) || gross_debt_return <= 0) {
    stop_model("El retorno fiscal implicito no es positivo.")
  }
  list(
    output = output,
    public_debt = public_debt,
    debt_to_annual_gdp = if (output > 0) {
      closure$years_per_period * public_debt / output
    } else NA_real_,
    primary_balance = primary_balance,
    primary_balance_to_gdp = primary_balance / output,
    years_per_period = closure$years_per_period,
    domestic_debt_share = closure$domestic_debt_share,
    gross_debt_return = gross_debt_return,
    debt_interest_rate = gross_debt_return - 1,
    annual_debt_interest_rate =
      gross_debt_return^(1 / closure$years_per_period) - 1
  )
}
