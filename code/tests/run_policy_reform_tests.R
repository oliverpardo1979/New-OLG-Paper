args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- if (length(file_arg) == 1L) {
  dirname(normalizePath(file_arg))
} else {
  getwd()
}
project_dir <- normalizePath(file.path(script_dir, ".."))
source(file.path(project_dir, "R", "load_endogenous_model.R"))

initial_param <- default_endogenous_parameters()
initial_coarse <- solve_steady_state_endogenous(
  initial_param,
  grid = make_type_grid(101L),
  old_weight = 0.90,
  tolerance = 1e-6,
  residual_tolerance = 1e-4
)
initial <- solve_steady_state_endogenous(
  initial_param,
  initial = as.list(initial_coarse$state),
  grid = make_type_grid(101L),
  old_weight = 0.90,
  tolerance = 2e-7,
  residual_tolerance = 1e-5
)
stopifnot(initial$validated)

final_param <- make_law2381_parameters(initial)
stopifnot(
  final_param$regime == "pillars",
  abs(
    final_param$payg_income_threshold /
      final_param$minimum_wage_model_units - 2.30
  ) < 1e-12,
  abs(
    final_param$conditional_basic_pension /
      final_param$minimum_wage_model_units - 218846 / 1160000
  ) < 1e-12,
  abs(final_param$funded_account_rate - 0.132) < 1e-12,
  abs(final_param$payg_benefit_multiplier - 1) < 1e-12
)
final_coarse <- solve_steady_state_endogenous(
  final_param,
  initial = as.list(initial$state),
  grid = make_type_grid(101L),
  old_weight = 0.90,
  tolerance = 1e-6,
  residual_tolerance = 1e-4
)
final <- solve_steady_state_endogenous(
  final_param,
  initial = as.list(final_coarse$state),
  grid = make_type_grid(101L),
  old_weight = 0.90,
  tolerance = 2e-7,
  residual_tolerance = 1e-5
)
stopifnot(final$validated)
shares <- c(
  informal = final$cohort$informal_share,
  payg_only = final$cohort$payg_only_share,
  mixed = final$cohort$mixed_share
)
print(c(
  shares,
  formal_entry = final$cohort$integration$roots[1L],
  pension_split = final$cohort$integration$pension_split_root,
  consumption_tax = final$state[["consumption_tax"]]
))
stopifnot(
  abs(sum(shares) - 1) < 1e-7,
  all(shares > 0),
  final$cohort$funded_only_share < 1e-10,
  final$cohort$payg_covered_share > final$cohort$payg_only_share,
  length(final$cohort$integration$pension_split_root) == 1L,
  all(final$cohort$micro$consumption_young > 0),
  all(final$cohort$micro$consumption_old > 0)
)

fiscal_rule <- make_fiscal_rule(
  final,
  initial_debt = 0,
  debt_target = 0,
  tax_anchor_weight = 0.05,
  tax_upper = 0.90,
  domestic_debt_share = 0.25
)
saved_transition_file <- file.path(
  project_dir, "..", "results", "policy_transition.rds"
)
initial_path <- NULL
if (file.exists(saved_transition_file)) {
  saved_transition <- readRDS(saved_transition_file)
  if (inherits(saved_transition, "pension_policy_transition") &&
      nrow(saved_transition$path) == 21L &&
      "public_debt" %in% names(saved_transition$path)) {
    initial_path <- saved_transition$path
  }
}

transition <- solve_transition_endogenous(
  initial_param = initial_param,
  final_param = final_param,
  periods = 20L,
  initial_solution = initial,
  final_solution = final,
  initial_path = initial_path,
  grid = make_type_grid(101L),
  old_weight = 0.90,
  fiscal_old_weight = 0.95,
  tolerance = 1e-4,
  terminal_tolerance = 1e-3,
  fiscal_tolerance = 2e-4,
  fiscal_rule = fiscal_rule,
  max_iterations = 1200L
)
print(c(
  converged_internal = transition$converged_internal,
  terminal_consistent = transition$terminal_consistent,
  validated = transition$validated,
  fiscal_consistent = transition$fiscal_consistent,
  iterations = transition$iterations,
  internal_gap = transition$internal_gap,
  terminal_gap = transition$terminal_gap,
  max_debt_identity_residual = transition$max_debt_identity_residual,
  max_fiscal_rule_residual = transition$max_fiscal_rule_residual
))
growth <- (1 + final_param$n) * (1 + final_param$g)
fiscal_rows <- 2:(nrow(transition$path) - 1L)
independent_debt_residual <-
  growth * transition$path$public_debt[fiscal_rows + 1L] -
  ((1 + transition$path$interest_rate[fiscal_rows]) *
    transition$path$public_debt[fiscal_rows] -
    transition$path$primary_balance[fiscal_rows])

stopifnot(
  transition$validated,
  transition$dynamic_outcome == "steady_state",
  transition$terminal_consistent,
  transition$fiscal_consistent,
  transition$max_debt_identity_residual < 2e-4,
  transition$max_fiscal_rule_residual < 2e-4,
  max(abs(independent_debt_residual)) < 2e-4,
  all(is.finite(transition$path$capital)),
  all(is.finite(transition$path$public_debt)),
  all(transition$path$consumption_tax > 0),
  all(transition$path$consumption_tax < fiscal_rule$tax_upper),
  abs(tail(transition$path$public_debt, 1L)) < 1e-5,
  all(transition$path$capital > 0),
  all(is.finite(transition$path$formal_share)),
  abs(
    transition$path$informal_share[1L] +
      transition$path$funded_only_share[1L] +
      transition$path$payg_only_share[1L] - 1
  ) < 1e-7,
  abs(
    tail(transition$path$informal_share, 1L) +
      tail(transition$path$payg_only_share, 1L) +
      tail(transition$path$mixed_share, 1L) - 1
  ) < 1e-7
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

cat("Todas las pruebas de reforma y distribucion pasaron.\n")
print(transition)
print(shares)
print(distribution)
