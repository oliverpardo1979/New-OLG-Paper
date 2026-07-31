# Calibracion indirecta con momentos del CEDE Pension Model (Becerra, 2026).
#
# Los dos parametros libres son la productividad informal y el intercepto
# de elegibilidad PAYG. Los demas valores se mantienen fijos.

becerra2026_targets <- function(param = default_endogenous_parameters()) {
  c(
    formal_share = param$target_formal_share,
    formal_eligibility_probability =
      param$target_contributory_coverage /
        param$target_formal_share
  )
}

decode_becerra2026_parameters <- function(theta, param) {
  if (!is.numeric(theta) || length(theta) != 2L ||
      any(!is.finite(theta))) {
    stop_model("theta debe contener dos valores finitos.")
  }
  candidate <- copy_list(param)
  candidate$A_informal <- exp(theta[1L])
  candidate$payg_eligibility_intercept <- theta[2L]
  candidate
}

encode_becerra2026_parameters <- function(param) {
  c(
    log(param$A_informal),
    param$payg_eligibility_intercept
  )
}

becerra2026_moments <- function(solution) {
  c(
    formal_share = solution$cohort$formal_share,
    formal_eligibility_probability =
      solution$cohort$formal_eligibility_probability
  )
}

becerra2026_validation <- function(solution, param = solution$param) {
  data.frame(
    moment = c("payg_share_among_formal", "annual_real_return"),
    target = c(
      param$target_payg_share_among_formal,
      param$target_annual_real_return
    ),
    model = c(
      solution$cohort$payg_share_among_formal,
      (1 + solution$prices$r)^(1 / 40) - 1
    )
  )
}

becerra2026_moment_residuals <- function(solution, targets) {
  scales <- c(
    formal_share = 0.01,
    formal_eligibility_probability = 0.02
  )
  (becerra2026_moments(solution) - targets) / scales
}

calibrate_becerra2026 <- function(
    param = default_endogenous_parameters(),
    grid_size = 101L,
    maxit = 80L,
    trace = FALSE
) {
  validate_endogenous_parameters(param)
  targets <- becerra2026_targets(param)
  grid <- make_type_grid(grid_size)
  last_state <- NULL
  evaluations <- 0L

  objective <- function(theta) {
    evaluations <<- evaluations + 1L
    candidate <- decode_becerra2026_parameters(theta, param)
    result <- tryCatch(
      solve_steady_state_endogenous(
        candidate,
        initial = last_state,
        grid = grid,
        old_weight = 0.90,
        tolerance = 3e-6,
        residual_tolerance = 5e-4,
        max_iterations = 2500L
      ),
      error = function(e) NULL
    )
    if (is.null(result) || !result$converged ||
        any(!is.finite(result$state))) {
      return(1e8)
    }
    last_state <<- as.list(result$state)
    residuals <- becerra2026_moment_residuals(result, targets)
    value <- sum(residuals^2)
    if (trace) {
      message(sprintf(
        "eval=%d objective=%.6f residuals=[%s]",
        evaluations,
        value,
        paste(sprintf("%.4f", residuals), collapse = ", ")
      ))
    }
    value
  }

  fit <- optim(
    encode_becerra2026_parameters(param),
    objective,
    method = "Nelder-Mead",
    control = list(maxit = maxit, reltol = 1e-7)
  )
  calibrated <- decode_becerra2026_parameters(fit$par, param)
  coarse <- solve_steady_state_endogenous(
    calibrated,
    grid = make_type_grid(101L),
    old_weight = 0.90,
    tolerance = 1e-6,
    residual_tolerance = 1e-4
  )
  solution <- solve_steady_state_endogenous(
    calibrated,
    initial = as.list(coarse$state),
    grid = make_type_grid(501L),
    old_weight = 0.90,
    tolerance = 1e-8,
    residual_tolerance = 2e-6
  )

  list(
    parameters = calibrated,
    solution = solution,
    targets = targets,
    moments = becerra2026_moments(solution),
    validation = becerra2026_validation(solution, calibrated),
    residuals = becerra2026_moments(solution) - targets,
    objective = fit$value,
    convergence = fit$convergence,
    evaluations = evaluations
  )
}
