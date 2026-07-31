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

transition <- solve_transition_endogenous(
  initial_param = initial_param,
  final_param = final_param,
  periods = 16L,
  initial_solution = initial,
  final_solution = final,
  grid = make_type_grid(101L),
  old_weight = 0.85,
  tolerance = 2e-5,
  terminal_tolerance = 1e-2,
  budget_tolerance = 2e-4,
  max_iterations = 600L
)
print(c(
  converged_internal = transition$converged_internal,
  terminal_consistent = transition$terminal_consistent,
  validated = transition$validated,
  cycle_detected = transition$cycle_detected,
  cycle_gap = transition$cycle_gap,
  phase_gap = transition$phase_gap,
  iterations = transition$iterations,
  internal_gap = transition$internal_gap,
  terminal_gap = transition$terminal_gap,
  max_budget_residual = transition$max_budget_residual
))
stopifnot(
  transition$validated,
  transition$dynamic_outcome == "two_cycle",
  transition$cycle_gap < 1e-2,
  !transition$terminal_consistent,
  all(is.finite(transition$path$capital)),
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
