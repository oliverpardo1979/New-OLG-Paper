solve_step_transition_stable_v2 <- solve_step_transition_stable

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
  if (abs(baseline$capital - final$capital) >= 1e-10) {
    return(solve_step_transition_stable_v2(
      param,
      marginal_payroll_share,
      label
    ))
  }

  levels <- scenario_levels(
    param,
    spending_increase = param$permanent_spending_increase,
    marginal_payroll_share = marginal_payroll_share
  )
  path <- data.frame(
    time = 0:param$horizon,
    capital = c(
      baseline$capital,
      rep(final$capital, param$horizon)
    ),
    consumption = c(
      baseline$consumption,
      rep(final$consumption, param$horizon)
    ),
    government_spending = c(
      baseline$government_spending,
      rep(levels$government_spending, param$horizon)
    ),
    payroll_requirement = c(
      param$baseline_payroll_requirement,
      rep(levels$payroll_requirement, param$horizon)
    )
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
    terminal_gap = 0,
    data = path,
    final_steady = final,
    stable_eigenvalue = NA_real_
  )
}
