solve_stationary_equilibrium <- function(dependency_ratio, mix, par,
                                         debt_ratio = par$debt_anchor,
                                         start = NULL,
                                         tolerance = 1e-9) {
  mix <- validate_mix(mix)
  supplied_start <- !is.null(start)
  if (is.null(start)) {
    start <- c(log(1.7), qlogis(0.45), log(0.85), log(0.70))
  }

  objective <- function(z) {
    x <- steady_components(z, dependency_ratio, mix, par, debt_ratio)
    if (any(!is.finite(x$residuals))) return(1e12)
    sum(x$residuals^2) + equilibrium_penalty(x, par)
  }

  starts <- if (supplied_start) {
    list(start)
  } else {
    list(
      start,
      c(log(1.0), qlogis(0.45), log(0.70), log(0.60)),
      c(log(2.5), qlogis(0.55), log(1.00), log(0.85)),
      c(log(0.7), qlogis(0.35), log(1.20), log(0.45))
    )
  }
  candidates <- lapply(starts, function(s) {
    tryCatch(
      optim(s, objective, method = "BFGS",
            control = list(maxit = 800, reltol = 1e-12)),
      error = function(e) NULL
    )
  })
  values <- vapply(candidates, function(x) {
    if (is.null(x)) Inf else x$value
  }, numeric(1))
  best <- candidates[[which.min(values)]]

  if (is.null(best) || !is.finite(best$value) || best$value > 1e-5) {
    nm <- tryCatch(
      optim(start, objective, method = "Nelder-Mead",
            control = list(maxit = 3000, reltol = 1e-12)),
      error = function(e) NULL
    )
    if (!is.null(nm) && is.finite(nm$value) &&
        (is.null(best) || nm$value < best$value)) {
      best <- nm
    }
  }

  if (is.null(best)) {
    return(list(converged = FALSE, objective = Inf))
  }
  result <- steady_components(
    best$par, dependency_ratio, mix, par, debt_ratio
  )
  result$converged <- is.finite(best$value) &&
    best$value < tolerance &&
    equilibrium_penalty(result, par) < 1e-8
  result$objective <- best$value
  result$optimizer_code <- best$convergence
  result
}

calibrate_parameters <- function(par,
                                 dependency_ratio,
                                 reference_mix = c(0.45, 0.35, 0.20),
                                 target_informality = 0.554,
                                 target_outlay_share = 0.057,
                                 debt_ratio = par$debt_anchor) {
  calibration_objective <- function(theta) {
    trial <- par
    trial$old_age_outlay <- exp(theta[1])
    trial$formal_cost <- theta[2]
    eq <- solve_stationary_equilibrium(
      dependency_ratio, reference_mix, trial, debt_ratio
    )
    if (!isTRUE(eq$converged)) return(1e5 + eq$objective)
    500 * (eq$informality - target_informality)^2 +
      500 * (eq$social_security_share - target_outlay_share)^2
  }

  fit <- optim(
    c(log(par$old_age_outlay), par$formal_cost),
    calibration_objective,
    method = "Nelder-Mead",
    control = list(maxit = 500, reltol = 1e-11)
  )
  par$old_age_outlay <- exp(fit$par[1])
  par$formal_cost <- fit$par[2]
  eq <- solve_stationary_equilibrium(
    dependency_ratio, reference_mix, par, debt_ratio
  )
  list(
    parameters = par,
    equilibrium = eq,
    objective = fit$value,
    converged = isTRUE(eq$converged) &&
      abs(eq$informality - target_informality) < 2e-4 &&
      abs(eq$social_security_share - target_outlay_share) < 2e-4
  )
}

simplex_grid <- function(step = 0.05) {
  n <- round(1 / step)
  if (abs(n * step - 1) > 1e-10) {
    stop("The grid step must divide one exactly.")
  }
  out <- vector("list", (n + 1) * (n + 2) / 2)
  index <- 1
  for (i in 0:n) {
    for (j in 0:(n - i)) {
      out[[index]] <- c(
        consumption = i / n,
        payroll = j / n,
        capital = (n - i - j) / n
      )
      index <- index + 1
    }
  }
  do.call(rbind, out)
}

equilibrium_row <- function(eq) {
  fields <- c(
    "dependency_ratio", "capital", "labor_formal", "informality",
    "price_informal", "nominal_output", "wage_formal", "wage_informal",
    "rental_rate", "social_security", "social_security_share",
    "public_absorption", "collection_cost", "collection_cost_share",
    "consumption_formal", "consumption_informal",
    "composite_consumption", "consumption_price_index",
    "primary_surplus_ratio", "required_revenue", "revenue_share",
    "tax_consumption", "tax_payroll", "tax_capital",
    "mix_consumption", "mix_payroll", "mix_capital",
    "flow_utility", "lifetime_welfare", "objective"
  )
  values <- vapply(fields, function(nm) {
    value <- eq[[nm]]
    if (is.null(value)) NA_real_ else as.numeric(value)
  }, numeric(1))
  as.data.frame(as.list(values), check.names = FALSE)
}

evaluate_tax_grid <- function(dependency_ratio, par,
                              debt_ratio = par$debt_anchor,
                              step = 0.05) {
  mixes <- simplex_grid(step)
  rows <- vector("list", nrow(mixes))
  start <- NULL
  for (i in seq_len(nrow(mixes))) {
    eq <- solve_stationary_equilibrium(
      dependency_ratio,
      mixes[i, ],
      par,
      debt_ratio,
      start = start
    )
    if (isTRUE(eq$converged)) start <- eq$z
    row <- equilibrium_row(eq)
    row$converged <- isTRUE(eq$converged)
    rows[[i]] <- row
  }
  out <- do.call(rbind, rows)
  out$welfare_loss_pct <- NA_real_
  feasible <- out$converged & is.finite(out$lifetime_welfare)
  if (any(feasible)) {
    best <- max(out$lifetime_welfare[feasible])
    c_best <- max(out$composite_consumption[feasible])
    out$welfare_loss_pct[feasible] <-
      100 * (c_best - out$composite_consumption[feasible]) / c_best
    out$is_optimum <- FALSE
    out$is_optimum[which.max(ifelse(feasible, out$lifetime_welfare, -Inf))] <-
      TRUE
  } else {
    out$is_optimum <- FALSE
  }
  out
}

optimal_grid_row <- function(grid) {
  feasible <- grid$converged & is.finite(grid$lifetime_welfare)
  if (!any(feasible)) stop("No feasible tax mix found.")
  grid[which.max(ifelse(feasible, grid$lifetime_welfare, -Inf)), ,
       drop = FALSE]
}

solve_demographic_path <- function(demographics, mix, par,
                                   debt_ratio = par$debt_anchor) {
  rows <- vector("list", nrow(demographics))
  start <- NULL
  for (i in seq_len(nrow(demographics))) {
    eq <- solve_stationary_equilibrium(
      demographics$old_age_dependency[i],
      mix,
      par,
      debt_ratio,
      start
    )
    if (!isTRUE(eq$converged)) {
      stop("Conditional steady state failed in year ",
           demographics$year[i])
    }
    start <- eq$z
    row <- equilibrium_row(eq)
    row$year <- demographics$year[i]
    rows[[i]] <- row
  }
  out <- do.call(rbind, rows)
  out[, c("year", setdiff(names(out), "year"))]
}
