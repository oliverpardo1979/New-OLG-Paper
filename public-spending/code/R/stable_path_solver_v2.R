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
    post <- data.frame(
      time = 1:param$horizon,
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
    reference_amplitude <- 1e-4 * max(1, final$capital)

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
      if (collect) return(states)
      state$capital
    }

    state <- list(
      capital = final$capital +
        direction * reference_amplitude,
      consumption = final$consumption +
        direction * reference_amplitude * stable_vector[2L]
    )
    steps <- 0L
    max_steps <- max(800L, param$horizon)
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
      if (crossed || steps >= max_steps) break
    }
    if (!crossed) {
      stop("The stable manifold did not reach inherited capital.")
    }

    amplitude_gap <- function(amplitude) {
      trace_from_terminal(amplitude, steps) -
        baseline$capital
    }
    upper <- reference_amplitude
    upper_gap <- amplitude_gap(upper)
    lower <- upper / 2
    lower_gap <- amplitude_gap(lower)
    while (
      is.finite(lower_gap) && is.finite(upper_gap) &&
        lower_gap * upper_gap > 0 &&
        lower > 1e-12
    ) {
      lower <- lower / 2
      lower_gap <- amplitude_gap(lower)
    }
    if (!is.finite(lower_gap) || !is.finite(upper_gap) ||
        lower_gap * upper_gap > 0) {
      stop("Could not bracket the stable-manifold displacement.")
    }
    root_amplitude <- uniroot(
      amplitude_gap,
      lower = lower,
      upper = upper,
      tol = 1e-8
    )$root
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
    } else if (nrow(post) > param$horizon) {
      param$horizon <- nrow(post)
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
  fields <- c(
    "output", "formal_output", "informal_output", "formal_share",
    "informal_share", "payroll_tax", "net_return", "wage_formal"
  )
  for (field in fields) {
    path[[field]] <- vapply(evaluated, `[[`, numeric(1), field)
  }
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
    final_steady = final,
    stable_eigenvalue = Re(
      decomposition$values[stable_index] %||% NA_real_
    )
  )
}

`%||%` <- function(left, right) {
  if (length(left) == 0L || is.null(left)) right else left
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
