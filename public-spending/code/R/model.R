logistic <- function(x) {
  1 / (1 + exp(-pmax(-700, pmin(700, x))))
}

crra_utility <- function(consumption, sigma) {
  if (any(!is.finite(consumption)) || any(consumption <= 0)) {
    stop("Consumption must be finite and strictly positive.")
  }
  if (abs(sigma - 1) < 1e-10) {
    return(log(consumption))
  }
  (consumption^(1 - sigma) - 1) / (1 - sigma)
}

default_parameters <- function() {
  list(
    alpha = 0.33,
    depreciation = 0.06,
    beta = 0.96,
    sigma = 2,
    A_informal = 0.60,
    choice_scale = 0.30,
    target_informal_share = 0.554,
    target_government_share = 0.1472,
    baseline_payroll_finance_share = 0.50,
    permanent_spending_increase = 0.02,
    public_utility_weight = 0.15,
    horizon = 120L,
    gradual_speed = 0.20
  )
}

validate_parameters <- function(param) {
  required <- c(
    "alpha", "depreciation", "beta", "sigma", "A_informal",
    "choice_scale", "target_informal_share",
    "target_government_share", "baseline_payroll_finance_share",
    "permanent_spending_increase", "public_utility_weight",
    "horizon", "gradual_speed"
  )
  missing <- setdiff(required, names(param))
  if (length(missing) > 0L) {
    stop("Missing parameters: ", paste(missing, collapse = ", "))
  }
  if (param$alpha <= 0 || param$alpha >= 1) {
    stop("alpha must lie strictly between zero and one.")
  }
  if (param$beta <= 0 || param$beta >= 1) {
    stop("beta must lie strictly between zero and one.")
  }
  if (param$depreciation <= 0 || param$depreciation >= 1) {
    stop("depreciation must lie strictly between zero and one.")
  }
  if (param$target_informal_share <= 0 ||
      param$target_informal_share >= 1) {
    stop("target_informal_share must lie strictly between zero and one.")
  }
  if (param$target_government_share <= 0 ||
      param$target_government_share >= 1) {
    stop("target_government_share must lie strictly between zero and one.")
  }
  invisible(TRUE)
}

steady_capital_per_formal_worker <- function(param) {
  target_net_return <- 1 / param$beta - 1
  (
    param$alpha /
      (target_net_return + param$depreciation)
  )^(1 / (1 - param$alpha))
}

calibrate_parameters <- function(param = default_parameters()) {
  validate_parameters(param)

  formal_share <- 1 - param$target_informal_share
  capital_per_formal <- steady_capital_per_formal_worker(param)
  capital <- capital_per_formal * formal_share
  formal_output <- capital^param$alpha *
    formal_share^(1 - param$alpha)
  informal_output <- param$A_informal * (1 - formal_share)
  output <- formal_output + informal_output
  government_spending <- param$target_government_share * output
  wage_formal <- (1 - param$alpha) *
    (capital / formal_share)^param$alpha
  payroll_requirement <- param$baseline_payroll_finance_share *
    government_spending
  payroll_tax <- payroll_requirement / (wage_formal * formal_share)

  if (payroll_tax >= 0.95) {
    stop("The calibration implies an inadmissibly high payroll tax.")
  }

  target_log_odds <- qlogis(formal_share)
  formality_cost <- log((1 - payroll_tax) * wage_formal) -
    log(param$A_informal) -
    param$choice_scale * target_log_odds

  param$capital_per_formal <- capital_per_formal
  param$baseline_formal_share <- formal_share
  param$baseline_capital <- capital
  param$baseline_output <- output
  param$baseline_government_spending <- government_spending
  param$baseline_payroll_requirement <- payroll_requirement
  param$formality_cost <- formality_cost
  param$baseline_payroll_tax <- payroll_tax
  param
}

labor_market_residual <- function(
    formal_share,
    capital,
    payroll_requirement,
    param
) {
  if (formal_share <= 0 || formal_share >= 1 || capital <= 0) {
    return(NA_real_)
  }
  wage_formal <- (1 - param$alpha) *
    (capital / formal_share)^param$alpha
  formal_wage_bill <- wage_formal * formal_share
  payroll_tax <- payroll_requirement / formal_wage_bill
  if (!is.finite(payroll_tax) || payroll_tax >= 0.999999) {
    return(NA_real_)
  }
  index <- (
    log((1 - payroll_tax) * wage_formal) -
      log(param$A_informal) -
      param$formality_cost
  ) / param$choice_scale
  qlogis(formal_share) - index
}

solve_formal_share <- function(
    capital,
    payroll_requirement,
    param,
    guess = param$baseline_formal_share
) {
  if (!is.finite(capital) || capital <= 0) {
    stop("Capital must be finite and strictly positive.")
  }
  grid <- seq(1e-4, 0.9999, length.out = 800L)
  residuals <- vapply(
    grid,
    labor_market_residual,
    numeric(1),
    capital = capital,
    payroll_requirement = payroll_requirement,
    param = param
  )
  valid <- is.finite(residuals)
  crossings <- which(
    valid[-length(valid)] & valid[-1L] &
      residuals[-length(residuals)] * residuals[-1L] <= 0
  )

  if (length(crossings) > 0L) {
    candidates <- vapply(
      crossings,
      function(index) {
        uniroot(
          labor_market_residual,
          lower = grid[index],
          upper = grid[index + 1L],
          capital = capital,
          payroll_requirement = payroll_requirement,
          param = param,
          tol = 1e-12
        )$root
      },
      numeric(1)
    )
    return(candidates[which.min(abs(candidates - guess))])
  }

  lower <- min(grid[valid])
  upper <- max(grid[valid])
  fallback <- optimize(
    function(formal_share) {
      value <- labor_market_residual(
        formal_share,
        capital,
        payroll_requirement,
        param
      )
      if (!is.finite(value)) return(1e12)
      value^2
    },
    interval = c(lower, upper),
    tol = 1e-12
  )
  if (fallback$objective > 1e-10) {
    stop("Could not solve the formal-informal allocation.")
  }
  fallback$minimum
}

evaluate_economy <- function(
    capital,
    government_spending,
    payroll_requirement,
    param,
    formal_guess = param$baseline_formal_share
) {
  formal_share <- solve_formal_share(
    capital,
    payroll_requirement,
    param,
    guess = formal_guess
  )
  informal_share <- 1 - formal_share
  formal_output <- capital^param$alpha *
    formal_share^(1 - param$alpha)
  informal_output <- param$A_informal * informal_share
  output <- formal_output + informal_output
  wage_formal <- (1 - param$alpha) *
    (capital / formal_share)^param$alpha
  formal_wage_bill <- wage_formal * formal_share
  payroll_tax <- payroll_requirement / formal_wage_bill
  net_return <- param$alpha * formal_output / capital -
    param$depreciation
  lump_sum_revenue <- government_spending - payroll_requirement

  list(
    capital = capital,
    formal_share = formal_share,
    informal_share = informal_share,
    formal_output = formal_output,
    informal_output = informal_output,
    output = output,
    wage_formal = wage_formal,
    wage_informal = param$A_informal,
    formal_wage_bill = formal_wage_bill,
    payroll_tax = payroll_tax,
    net_return = net_return,
    government_spending = government_spending,
    payroll_revenue = payroll_tax * formal_wage_bill,
    lump_sum_revenue = lump_sum_revenue,
    government_residual = payroll_tax * formal_wage_bill +
      lump_sum_revenue - government_spending
  )
}

scenario_levels <- function(
    param,
    spending_increase = param$permanent_spending_increase,
    marginal_payroll_share = 0
) {
  spending_change <- spending_increase * param$baseline_output
  list(
    government_spending =
      param$baseline_government_spending + spending_change,
    payroll_requirement =
      param$baseline_payroll_requirement +
      marginal_payroll_share * spending_change
  )
}

solve_steady_state <- function(
    param,
    spending_increase = 0,
    marginal_payroll_share = 0,
    label = "Baseline"
) {
  levels <- scenario_levels(
    param,
    spending_increase = spending_increase,
    marginal_payroll_share = marginal_payroll_share
  )
  capital_per_formal <- steady_capital_per_formal_worker(param)

  formal_residual <- function(formal_share) {
    capital <- capital_per_formal * formal_share
    labor_market_residual(
      formal_share,
      capital,
      levels$payroll_requirement,
      param
    )
  }

  grid <- seq(1e-4, 0.9999, length.out = 800L)
  residuals <- vapply(grid, formal_residual, numeric(1))
  valid <- is.finite(residuals)
  crossings <- which(
    valid[-length(valid)] & valid[-1L] &
      residuals[-length(residuals)] * residuals[-1L] <= 0
  )
  if (length(crossings) == 0L) {
    stop("No stationary formal-sector allocation was found.")
  }
  candidate_roots <- vapply(
    crossings,
    function(index) {
      uniroot(
        formal_residual,
        lower = grid[index],
        upper = grid[index + 1L],
        tol = 1e-12
      )$root
    },
    numeric(1)
  )
  formal_share <- candidate_roots[
    which.min(abs(candidate_roots - param$baseline_formal_share))
  ]
  capital <- capital_per_formal * formal_share
  economy <- evaluate_economy(
    capital,
    levels$government_spending,
    levels$payroll_requirement,
    param,
    formal_guess = formal_share
  )
  consumption <- economy$output -
    levels$government_spending -
    param$depreciation * capital
  if (!is.finite(consumption) || consumption <= 0) {
    stop("The stationary allocation has non-positive consumption.")
  }

  data.frame(
    scenario = label,
    spending_increase_initial_gdp = spending_increase,
    marginal_payroll_share = marginal_payroll_share,
    capital = capital,
    output = economy$output,
    formal_output = economy$formal_output,
    informal_output = economy$informal_output,
    consumption = consumption,
    investment = param$depreciation * capital,
    government_spending = levels$government_spending,
    government_share_output =
      levels$government_spending / economy$output,
    formal_share = economy$formal_share,
    informal_share = economy$informal_share,
    payroll_tax = economy$payroll_tax,
    net_return = economy$net_return,
    wage_formal = economy$wage_formal,
    wage_informal = economy$wage_informal,
    government_residual = economy$government_residual,
    resource_residual = economy$output - consumption -
      levels$government_spending - param$depreciation * capital,
    euler_residual = param$beta * (1 + economy$net_return) - 1
  )
}

spending_path <- function(
    param,
    path = c("step", "gradual"),
    horizon = param$horizon
) {
  path <- match.arg(path)
  time <- 0:horizon
  if (path == "step") {
    completion <- as.numeric(time > 0)
  } else {
    completion <- 1 - exp(-param$gradual_speed * time)
    completion[1L] <- 0
  }
  spending_change <- param$permanent_spending_increase *
    param$baseline_output
  data.frame(
    time = time,
    completion = completion,
    government_spending =
      param$baseline_government_spending +
      completion * spending_change
  )
}

transition_given_initial_consumption <- function(
    initial_consumption,
    param,
    marginal_payroll_share,
    path_type = "step"
) {
  policy <- spending_path(param, path_type)
  n_periods <- nrow(policy)
  capital <- rep(NA_real_, n_periods)
  consumption <- rep(NA_real_, n_periods)
  formal_share <- rep(NA_real_, n_periods)
  informal_share <- rep(NA_real_, n_periods)
  output <- rep(NA_real_, n_periods)
  formal_output <- rep(NA_real_, n_periods)
  informal_output <- rep(NA_real_, n_periods)
  payroll_tax <- rep(NA_real_, n_periods)
  net_return <- rep(NA_real_, n_periods)
  wage_formal <- rep(NA_real_, n_periods)
  payroll_requirement <- param$baseline_payroll_requirement +
    marginal_payroll_share *
      (policy$government_spending -
        param$baseline_government_spending)

  baseline <- solve_steady_state(param)
  capital[1L] <- baseline$capital
  consumption[1L] <- baseline$consumption
  baseline_economy <- evaluate_economy(
    capital[1L],
    policy$government_spending[1L],
    payroll_requirement[1L],
    param
  )
  formal_share[1L] <- baseline_economy$formal_share
  informal_share[1L] <- baseline_economy$informal_share
  output[1L] <- baseline_economy$output
  formal_output[1L] <- baseline_economy$formal_output
  informal_output[1L] <- baseline_economy$informal_output
  payroll_tax[1L] <- baseline_economy$payroll_tax
  net_return[1L] <- baseline_economy$net_return
  wage_formal[1L] <- baseline_economy$wage_formal

  consumption[2L] <- initial_consumption
  capital[2L] <- capital[1L]
  previous_formal_share <- formal_share[1L]

  for (index in 2:(n_periods - 1L)) {
    economy <- tryCatch(
      evaluate_economy(
        capital[index],
        policy$government_spending[index],
        payroll_requirement[index],
        param,
        formal_guess = previous_formal_share
      ),
      error = function(error) NULL
    )
    if (is.null(economy)) {
      return(list(valid = FALSE, terminal_gap = -1e6))
    }
    formal_share[index] <- economy$formal_share
    informal_share[index] <- economy$informal_share
    output[index] <- economy$output
    formal_output[index] <- economy$formal_output
    informal_output[index] <- economy$informal_output
    payroll_tax[index] <- economy$payroll_tax
    net_return[index] <- economy$net_return
    wage_formal[index] <- economy$wage_formal
    previous_formal_share <- economy$formal_share

    next_capital <- (1 - param$depreciation) * capital[index] +
      economy$output -
      consumption[index] -
      policy$government_spending[index]
    if (!is.finite(next_capital) || next_capital <= 1e-8) {
      return(list(valid = FALSE, terminal_gap = -1e6))
    }
    capital[index + 1L] <- next_capital

    next_economy <- tryCatch(
      evaluate_economy(
        next_capital,
        policy$government_spending[index + 1L],
        payroll_requirement[index + 1L],
        param,
        formal_guess = previous_formal_share
      ),
      error = function(error) NULL
    )
    if (is.null(next_economy) ||
        !is.finite(1 + next_economy$net_return) ||
        1 + next_economy$net_return <= 0) {
      return(list(valid = FALSE, terminal_gap = -1e6))
    }
    consumption[index + 1L] <- consumption[index] * (
      param$beta * (1 + next_economy$net_return)
    )^(1 / param$sigma)
  }

  last <- n_periods
  terminal_economy <- tryCatch(
    evaluate_economy(
      capital[last],
      policy$government_spending[last],
      payroll_requirement[last],
      param,
      formal_guess = previous_formal_share
    ),
    error = function(error) NULL
  )
  if (is.null(terminal_economy)) {
    return(list(valid = FALSE, terminal_gap = -1e6))
  }
  formal_share[last] <- terminal_economy$formal_share
  informal_share[last] <- terminal_economy$informal_share
  output[last] <- terminal_economy$output
  formal_output[last] <- terminal_economy$formal_output
  informal_output[last] <- terminal_economy$informal_output
  payroll_tax[last] <- terminal_economy$payroll_tax
  net_return[last] <- terminal_economy$net_return
  wage_formal[last] <- terminal_economy$wage_formal

  final_steady <- solve_steady_state(
    param,
    spending_increase = param$permanent_spending_increase,
    marginal_payroll_share = marginal_payroll_share
  )
  terminal_gap <- capital[last] - final_steady$capital

  data <- data.frame(
    time = policy$time,
    capital = capital,
    output = output,
    formal_output = formal_output,
    informal_output = informal_output,
    consumption = consumption,
    government_spending = policy$government_spending,
    government_share_output = policy$government_spending / output,
    formal_share = formal_share,
    informal_share = informal_share,
    payroll_tax = payroll_tax,
    net_return = net_return,
    wage_formal = wage_formal,
    payroll_requirement = payroll_requirement
  )
  list(
    valid = TRUE,
    terminal_gap = terminal_gap,
    data = data,
    final_steady = final_steady
  )
}

solve_transition <- function(
    param,
    marginal_payroll_share,
    label,
    path_type = c("step", "gradual")
) {
  path_type <- match.arg(path_type)
  final_steady <- solve_steady_state(
    param,
    spending_increase = param$permanent_spending_increase,
    marginal_payroll_share = marginal_payroll_share
  )
  policy <- spending_path(param, path_type)
  initial_state <- evaluate_economy(
    param$baseline_capital,
    policy$government_spending[2L],
    param$baseline_payroll_requirement +
      marginal_payroll_share *
        (policy$government_spending[2L] -
          param$baseline_government_spending),
    param
  )
  maximum_consumption <- (1 - param$depreciation) *
    param$baseline_capital +
    initial_state$output -
    policy$government_spending[2L]
  grid <- seq(
    0.05 * final_steady$consumption,
    0.995 * maximum_consumption,
    length.out = 300L
  )
  gaps <- vapply(
    grid,
    function(candidate) {
      transition_given_initial_consumption(
        candidate,
        param,
        marginal_payroll_share,
        path_type
      )$terminal_gap
    },
    numeric(1)
  )
  crossings <- which(
    is.finite(gaps[-length(gaps)]) &
      is.finite(gaps[-1L]) &
      gaps[-length(gaps)] * gaps[-1L] <= 0
  )
  if (length(crossings) == 0L) {
    stop("Could not bracket the transition saddle path for ", label, ".")
  }
  root_index <- crossings[which.min(
    abs(grid[crossings] - final_steady$consumption)
  )]
  initial_consumption <- uniroot(
    function(candidate) {
      transition_given_initial_consumption(
        candidate,
        param,
        marginal_payroll_share,
        path_type
      )$terminal_gap
    },
    lower = grid[root_index],
    upper = grid[root_index + 1L],
    tol = 1e-11
  )$root
  solution <- transition_given_initial_consumption(
    initial_consumption,
    param,
    marginal_payroll_share,
    path_type
  )
  if (!solution$valid) {
    stop("The transition solution is not valid for ", label, ".")
  }
  solution$data$scenario <- label
  solution$data$path_type <- path_type
  solution$data$capital_gap_final <-
    solution$data$capital - final_steady$capital
  solution$data$consumption_gap_final <-
    solution$data$consumption - final_steady$consumption
  solution
}

discounted_welfare <- function(
    transition,
    param,
    public_utility_weight = 0
) {
  data <- transition$data
  post_announcement <- data$time >= 1
  consumption <- data$consumption[post_announcement]
  spending <- data$government_spending[post_announcement]
  periods <- seq_along(consumption) - 1L
  flow <- crra_utility(consumption, param$sigma) +
    public_utility_weight * log(spending)
  terminal_flow <- crra_utility(
    transition$final_steady$consumption,
    param$sigma
  ) + public_utility_weight *
    log(transition$final_steady$government_spending)
  sum(param$beta^periods * flow) +
    param$beta^length(flow) * terminal_flow / (1 - param$beta)
}

baseline_welfare <- function(param, public_utility_weight = 0) {
  baseline <- solve_steady_state(param)
  flow <- crra_utility(baseline$consumption, param$sigma) +
    public_utility_weight *
      log(baseline$government_spending)
  flow / (1 - param$beta)
}

consumption_compensation <- function(
    transition,
    param,
    public_utility_weight = 0
) {
  target <- baseline_welfare(param, public_utility_weight)
  data <- transition$data
  post_announcement <- data$time >= 1
  consumption <- data$consumption[post_announcement]
  spending <- data$government_spending[post_announcement]
  periods <- seq_along(consumption) - 1L

  welfare_with_scale <- function(scale) {
    flow <- crra_utility(scale * consumption, param$sigma) +
      public_utility_weight * log(spending)
    terminal_flow <- crra_utility(
      scale * transition$final_steady$consumption,
      param$sigma
    ) + public_utility_weight *
      log(transition$final_steady$government_spending)
    sum(param$beta^periods * flow) +
      param$beta^length(flow) * terminal_flow / (1 - param$beta)
  }

  objective <- function(log_scale) {
    welfare_with_scale(exp(log_scale)) - target
  }
  bracket <- c(-4, 4)
  if (objective(bracket[1L]) * objective(bracket[2L]) > 0) {
    stop("Could not bracket the welfare compensation.")
  }
  scale <- exp(uniroot(objective, bracket, tol = 1e-11)$root)
  100 * (scale - 1)
}
