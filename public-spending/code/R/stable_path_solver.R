transition_map <- function(
    capital,
    consumption,
    government_spending,
    payroll_requirement,
    param,
    formal_guess = param$baseline_formal_share
) {
  current <- evaluate_economy(
    capital,
    government_spending,
    payroll_requirement,
    param,
    formal_guess = formal_guess
  )
  capital_next <- (1 - param$depreciation) * capital +
    current$output - consumption - government_spending
  if (!is.finite(capital_next) || capital_next <= 0) {
    stop("The transition map produced non-positive capital.")
  }
  next_economy <- evaluate_economy(
    capital_next,
    government_spending,
    payroll_requirement,
    param,
    formal_guess = current$formal_share
  )
  consumption_next <- consumption * (
    param$beta * (1 + next_economy$net_return)
  )^(1 / param$sigma)
  c(capital = capital_next, consumption = consumption_next)
}

transition_jacobian <- function(
    steady,
    government_spending,
    payroll_requirement,
    param
) {
  state <- c(
    capital = steady$capital,
    consumption = steady$consumption
  )
  jacobian <- matrix(NA_real_, nrow = 2L, ncol = 2L)
  for (column in seq_along(state)) {
    step <- 1e-5 * max(1, abs(state[column]))
    upper <- state
    lower <- state
    upper[column] <- upper[column] + step
    lower[column] <- lower[column] - step
    mapped_upper <- transition_map(
      upper["capital"],
      upper["consumption"],
      government_spending,
      payroll_requirement,
      param
    )
    mapped_lower <- transition_map(
      lower["capital"],
      lower["consumption"],
      government_spending,
      payroll_requirement,
      param
    )
    jacobian[, column] <- (mapped_upper - mapped_lower) / (2 * step)
  }
  dimnames(jacobian) <- list(
    c("capital_next", "consumption_next"),
    c("capital", "consumption")
  )
  jacobian
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

  reference <- max(
    param$baseline_capital,
    capital_next,
    param$capital_per_formal
  )
  grid <- exp(seq(
    log(max(1e-6, 0.03 * min(param$baseline_capital, capital_next))),
    log(3 * reference),
    length.out = 120L
  ))
  residuals <- vapply(grid, resource_residual, numeric(1))
  valid <- is.finite(residuals)
  crossings <- which(
    valid[-length(valid)] & valid[-1L] &
      residuals[-length(residuals)] * residuals[-1L] <= 0
  )
  if (length(crossings) == 0L) {
    stop("The inverse transition map has no admissible capital root.")
  }
  roots <- vapply(
    crossings,
    function(index) {
      uniroot(
        resource_residual,
        lower = grid[index],
        upper = grid[index + 1L],
        tol = 1e-11
      )$root
    },
    numeric(1)
  )
  directional <- direction * (roots - capital_next) >= -1e-7
  if (any(directional)) {
    roots <- roots[directional]
  }
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

solve_step_transition_stable <- function(
    param,
    marginal_payroll_share,
    label
) {
  baseline <- solve_steady_state(param)
  final <- solve_steady_state(
    param,
    spending_increase = param$permanent_spending_increase,
    marginal_payroll_share = marginal_payroll_share,
    label = label
  )
  levels <- scenario_levels(
    param,
    spending_increase = param$permanent_spending_increase,
    marginal_payroll_share = marginal_payroll_share
  )

  if (abs(baseline$capital - final$capital) < 1e-10) {
    post_time <- 1:param$horizon
    post <- data.frame(
      time = post_time,
      capital = final$capital,
      consumption = final$consumption
    )
  } else {
    jacobian <- transition_jacobian(
      final,
      levels$government_spending,
      levels$payroll_requirement,
      param
    )
    decomposition <- eigen(jacobian)
    stable_index <- which(abs(decomposition$values) < 1)
    if (length(stable_index) != 1L) {
      stop("The stationary system does not have a unique stable root.")
    }
    stable_vector <- Re(decomposition$vectors[, stable_index])
    stable_vector <- stable_vector / stable_vector[1L]
    direction <- sign(baseline$capital - final$capital)

    trace_from_terminal <- function(amplitude, steps, collect = FALSE) {
      state <- list(
        capital = final$capital + direction * amplitude,
        consumption = final$consumption +
          direction * amplitude * stable_vector[2L]
      )
      states <- if (collect) list(state) else NULL
      for (step in seq_len(steps)) {
        state <- inverse_transition_map(
          state$capital,
          state$consumption,
          levels$government_spending,
          levels$payroll_requirement,
          param,
          direction
        )
        if (collect) states[[length(states) + 1L]] <- state
      }
      if (collect) {
        return(states)
      }
      state$capital
    }

    amplitude <- 1e-8 * max(1, final$capital)
    state <- list(
      capital = final$capital + direction * amplitude,
      consumption = final$consumption +
        direction * amplitude * stable_vector[2L]
    )
    steps <- 0L
    repeat {
      steps <- steps + 1L
      state <- inverse_transition_map(
        state$capital,
        state$consumption,
        levels$government_spending,
        levels$payroll_requirement,
        param,
        direction
      )
      crossed <- direction *
        (state$capital - baseline$capital) >= 0
      if (crossed || steps >= param$horizon - 2L) break
    }
    if (!crossed) {
      stop("The stable manifold did not reach inherited capital.")
    }

    amplitude_gap <- function(log_amplitude) {
      trace_from_terminal(exp(log_amplitude), steps) -
        baseline$capital
    }
    upper_log <- log(amplitude)
    lower_log <- upper_log - 2
    lower_gap <- amplitude_gap(lower_log)
    upper_gap <- amplitude_gap(upper_log)
    while (
      is.finite(lower_gap) && is.finite(upper_gap) &&
        lower_gap * upper_gap > 0 &&
        lower_log > log(.Machine$double.xmin) + 10
    ) {
      lower_log <- lower_log - 2
      lower_gap <- amplitude_gap(lower_log)
    }
    if (!is.finite(lower_gap) || !is.finite(upper_gap) ||
        lower_gap * upper_gap > 0) {
      stop("Could not bracket the stable-manifold displacement.")
    }
    root_amplitude <- exp(uniroot(
      amplitude_gap,
      lower = lower_log,
      upper = upper_log,
      tol = 1e-11
    )$root)
    backward_states <- trace_from_terminal(
      root_amplitude,
      steps,
      collect = TRUE
    )
    forward_states <- rev(backward_states)
    post <- data.frame(
      time = seq_along(forward_states),
      capital = vapply(
        forward_states,
        function(state) state$capital,
        numeric(1)
      ),
      consumption = vapply(
        forward_states,
        function(state) state$consumption,
        numeric(1)
      )
    )
    if (nrow(post) < param$horizon) {
      extension <- data.frame(
        time = (nrow(post) + 1L):param$horizon,
        capital = final$capital,
        consumption = final$consumption
      )
      post <- rbind(post, extension)
    } else {
      post <- post[seq_len(param$horizon), ]
    }
  }

  pre <- data.frame(
    time = 0L,
    capital = baseline$capital,
    consumption = baseline$consumption
  )
  path <- rbind(pre, post)
  path$government_spending <- c(
    baseline$government_spending,
    rep(levels$government_spending, nrow(path) - 1L)
  )
  path$payroll_requirement <- c(
    param$baseline_payroll_requirement,
    rep(levels$payroll_requirement, nrow(path) - 1L)
  )

  evaluated <- lapply(
    seq_len(nrow(path)),
    function(index) {
      evaluate_economy(
        path$capital[index],
        path$government_spending[index],
        path$payroll_requirement[index],
        param
      )
    }
  )
  path$output <- vapply(evaluated, `[[`, numeric(1), "output")
  path$formal_output <- vapply(
    evaluated,
    `[[`,
    numeric(1),
    "formal_output"
  )
  path$informal_output <- vapply(
    evaluated,
    `[[`,
    numeric(1),
    "informal_output"
  )
  path$formal_share <- vapply(
    evaluated,
    `[[`,
    numeric(1),
    "formal_share"
  )
  path$informal_share <- vapply(
    evaluated,
    `[[`,
    numeric(1),
    "informal_share"
  )
  path$payroll_tax <- vapply(
    evaluated,
    `[[`,
    numeric(1),
    "payroll_tax"
  )
  path$net_return <- vapply(
    evaluated,
    `[[`,
    numeric(1),
    "net_return"
  )
  path$wage_formal <- vapply(
    evaluated,
    `[[`,
    numeric(1),
    "wage_formal"
  )
  path$government_share_output <-
    path$government_spending / path$output
  path$scenario <- label
  path$path_type <- "step"
  path$capital_gap_final <- path$capital - final$capital
  path$consumption_gap_final <-
    path$consumption - final$consumption

  list(
    valid = TRUE,
    terminal_gap = tail(path$capital_gap_final, 1L),
    data = path,
    final_steady = final
  )
}

solve_transition <- function(
    param,
    marginal_payroll_share,
    label,
    path_type = c("step")
) {
  path_type <- match.arg(path_type)
  solve_step_transition_stable(
    param,
    marginal_payroll_share,
    label
  )
}
