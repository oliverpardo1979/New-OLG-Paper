args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- if (length(file_arg) == 1L) {
  dirname(normalizePath(file_arg))
} else {
  getwd()
}
project_dir <- normalizePath(file.path(script_dir, ".."))
source(file.path(project_dir, "R", "load_endogenous_model.R"))

fiscal_closure <- make_fiscal_steady_state_closure(
  debt_to_annual_gdp = 0.55,
  primary_balance_to_gdp = 0.002,
  years_per_period = 40,
  domestic_debt_share = 0
)
grid <- make_type_grid(101L)
initial_param <- default_endogenous_parameters()
initial <- solve_steady_state_endogenous(
  initial_param,
  grid = grid,
  old_weight = 0.90,
  tolerance = 1e-6,
  residual_tolerance = 1e-4,
  fiscal_closure = fiscal_closure
)
stopifnot(initial$validated)

pillar_param <- make_law2381_parameters(initial)
funded_param <- make_fully_funded_parameters(initial)
stopifnot(
  pillar_param$regime == "pillars",
  funded_param$regime == "funded",
  abs(
    pillar_param$payg_income_threshold /
      pillar_param$minimum_wage_model_units - 2.30
  ) < 1e-12,
  abs(
    pillar_param$conditional_basic_pension /
      pillar_param$minimum_wage_model_units - 218846 / 1160000
  ) < 1e-12,
  abs(pillar_param$funded_account_rate - 0.132) < 1e-12,
  abs(pillar_param$payg_benefit_multiplier - 1) < 1e-12,
  abs(funded_param$funded_account_rate - funded_param$tau_pension) < 1e-12,
  abs(funded_param$payg_notional_account_rate) < 1e-12
)

solve_reform_state <- function(param, seed) {
  solve_steady_state_endogenous(
    param,
    initial = as.list(seed$state),
    grid = grid,
    old_weight = 0.90,
    tolerance = 1e-6,
    residual_tolerance = 1e-4,
    fiscal_closure = fiscal_closure
  )
}
pillar <- solve_reform_state(pillar_param, initial)
funded <- solve_reform_state(funded_param, initial)
stopifnot(pillar$validated, funded$validated)

pillar_shares <- c(
  informal = pillar$cohort$informal_share,
  payg_only = pillar$cohort$payg_only_share,
  mixed = pillar$cohort$mixed_share
)
funded_shares <- c(
  informal = funded$cohort$informal_share,
  funded_only = funded$cohort$funded_only_share
)
stopifnot(
  abs(sum(pillar_shares) - 1) < 1e-7,
  all(pillar_shares > 0),
  pillar$cohort$funded_only_share < 1e-10,
  pillar$cohort$payg_covered_share > pillar$cohort$payg_only_share,
  length(pillar$cohort$integration$pension_split_root) == 1L,
  abs(sum(funded_shares) - 1) < 1e-7,
  all(funded_shares > 0),
  funded$cohort$payg_only_share < 1e-10,
  funded$cohort$mixed_share < 1e-10,
  all(pillar$cohort$micro$consumption_young > 0),
  all(pillar$cohort$micro$consumption_old > 0),
  all(funded$cohort$micro$consumption_young > 0),
  all(funded$cohort$micro$consumption_old > 0)
)

steady_states <- list(initial, pillar, funded)
stopifnot(
  all(vapply(
    steady_states,
    function(x) abs(x$fiscal$debt_to_annual_gdp - 0.55) < 1e-10,
    logical(1)
  )),
  all(vapply(
    steady_states,
    function(x) abs(x$fiscal$primary_balance_to_gdp - 0.002) < 1e-10,
    logical(1)
  ))
)

solve_reform_transition <- function(final_param, final_solution) {
  rule <- make_fiscal_rule(
    final_solution,
    initial_debt = initial$fiscal$public_debt,
    tax_lower = -0.20,
    tax_upper = 0.90,
    debt_persistence = 0.25,
    tax_persistence = 0.70
  )
  transition <- solve_transition_endogenous(
    initial_param = initial_param,
    final_param = final_param,
    periods = 20L,
    initial_solution = initial,
    final_solution = final_solution,
    grid = grid,
    old_weight = 0.75,
    fiscal_old_weight = 0.50,
    tolerance = 1e-4,
    terminal_tolerance = 1e-3,
    fiscal_tolerance = 2e-4,
    fiscal_rule = rule,
    max_iterations = 400L
  )
  list(transition = transition, rule = rule)
}
pillar_run <- solve_reform_transition(pillar_param, pillar)
funded_run <- solve_reform_transition(funded_param, funded)

count_direction_reversals <- function(x, tolerance = 1e-8) {
  changes <- diff(x)
  directions <- sign(changes[abs(changes) > tolerance])
  if (length(directions) < 2L) return(0L)
  sum(directions[-1L] != directions[-length(directions)])
}
validate_transition <- function(run) {
  transition <- run$transition
  rule <- run$rule
  path <- transition$path
  growth <- (1 + transition$final_solution$param$n) *
    (1 + transition$final_solution$param$g)
  rows <- 2:(nrow(path) - 1L)
  independent_debt_residual <-
    growth * path$public_debt[rows + 1L] -
    ((1 + rule$debt_interest_rate) * path$public_debt[rows] -
       path$primary_balance[rows])
  reversals <- vapply(
    path[c("capital", "public_debt", "consumption_tax", "informal_share")],
    count_direction_reversals,
    integer(1)
  )
  stopifnot(
    transition$validated,
    transition$dynamic_outcome == "steady_state",
    transition$terminal_consistent,
    transition$fiscal_consistent,
    transition$max_debt_identity_residual < 2e-4,
    transition$max_fiscal_rule_residual < 2e-4,
    max(abs(independent_debt_residual)) < 2e-4,
    all(is.finite(path$capital)),
    all(path$capital > 0),
    all(is.finite(path$public_debt)),
    all(path$consumption_tax > rule$tax_lower),
    all(path$consumption_tax < rule$tax_upper),
    abs(tail(path$debt_to_output, 1L) - 0.55) < 1e-5,
    abs(path$fiscal_adjustment[1L]) < 1e-12,
    abs(tail(path$fiscal_adjustment, 1L)) < 1e-12,
    max(reversals) <= 2L
  )
  micro <- build_policy_microdata(transition)
  distribution <- summarize_policy_distributions(micro, initial_param)
  stopifnot(
    all(is.finite(micro$wage)),
    all(micro$wage > 0),
    all(is.finite(micro$welfare_cev)),
    all(is.finite(distribution$wage_gini)),
    all(distribution$winner_share >= 0 & distribution$winner_share <= 1),
    all(distribution$loser_share >= 0 & distribution$loser_share <= 1)
  )
  invisible(list(micro = micro, distribution = distribution))
}
pillar_distribution <- validate_transition(pillar_run)
funded_distribution <- validate_transition(funded_run)

cat("Todas las pruebas de reforma, regla fiscal y distribucion pasaron.\n")
print(pillar_run$transition)
print(funded_run$transition)
print(pillar_shares)
print(funded_shares)