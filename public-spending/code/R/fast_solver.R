# Faster local roots for calculations along a known equilibrium branch.

solve_formal_share <- function(
    capital,
    payroll_requirement,
    param,
    guess = param$baseline_formal_share
) {
  if (!is.finite(capital) || capital <= 0) {
    stop("Capital must be finite and strictly positive.")
  }
  tax_floor <- if (payroll_requirement > 0) {
    (
      payroll_requirement /
        ((1 - param$alpha) * capital^param$alpha * 0.999)
    )^(1 / (1 - param$alpha))
  } else {
    1e-5
  }
  lower <- max(1e-5, 1.00001 * tax_floor)
  if (!is.finite(lower) || lower >= 0.99999) {
    stop("The payroll requirement leaves no feasible formal tax base.")
  }

  local_lower <- max(lower, 0.65 * guess)
  local_upper <- min(0.99999, 1.35 * guess)
  local_grid <- seq(local_lower, local_upper, length.out = 24L)
  full_grid <- seq(lower, 0.99999, length.out = 48L)
  grid <- sort(unique(c(local_grid, full_grid)))
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
  if (length(crossings) == 0L) {
    stop("Could not solve the formal-informal allocation.")
  }
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
        tol = 1e-10
      )$root
    },
    numeric(1)
  )
  candidates[which.min(abs(candidates - guess))]
}

inverse_transition_map <- function(
    capital_next,
    consumption_next,
    government_spending,
    payroll_requirement,
    param,
    direction,
    formal_guess = param$baseline_formal_share
) {
  next_economy <- evaluate_economy(
    capital_next,
    government_spending,
    payroll_requirement,
    param,
    formal_guess = formal_guess
  )
  euler_growth <- (
    param$beta * (1 + next_economy$net_return)
  )^(1 / param$sigma)
  consumption <- consumption_next / euler_growth

  resource_residual <- function(capital) {
    economy <- tryCatch(
      evaluate_economy(
        capital,
        government_spending,
        payroll_requirement,
        param,
        formal_guess = next_economy$formal_share
      ),
      error = function(error) NULL
    )
    if (is.null(economy)) return(NA_real_)
    (1 - param$depreciation) * capital +
      economy$output -
      consumption -
      government_spending -
      capital_next
  }

  center <- capital_next
  span <- 0.01 * max(1, center)
  roots <- numeric(0)
  for (attempt in seq_len(18L)) {
    lower <- max(1e-7, center - span)
    upper <- center + span
    grid <- seq(lower, upper, length.out = 12L)
    residuals <- vapply(grid, resource_residual, numeric(1))
    valid <- is.finite(residuals)
    crossings <- which(
      valid[-length(valid)] & valid[-1L] &
        residuals[-length(residuals)] * residuals[-1L] <= 0
    )
    if (length(crossings) > 0L) {
      roots <- vapply(
        crossings,
        function(index) {
          uniroot(
            resource_residual,
            lower = grid[index],
            upper = grid[index + 1L],
            tol = 1e-10
          )$root
        },
        numeric(1)
      )
      break
    }
    span <- 1.65 * span
  }
  if (length(roots) == 0L) {
    stop("The inverse transition map has no admissible capital root.")
  }
  directional <- direction * (roots - capital_next) >= -1e-7
  if (any(directional)) roots <- roots[directional]
  capital <- roots[which.min(abs(roots - capital_next))]
  economy <- evaluate_economy(
    capital,
    government_spending,
    payroll_requirement,
    param,
    formal_guess = next_economy$formal_share
  )
  list(
    capital = capital,
    consumption = consumption,
    economy = economy
  )
}
