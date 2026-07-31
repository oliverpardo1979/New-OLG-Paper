args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- if (length(file_arg) == 1L) {
  dirname(normalizePath(file_arg))
} else {
  getwd()
}
project_dir <- normalizePath(file.path(script_dir, ".."))
source(file.path(project_dir, "R", "load_endogenous_model.R"))

param <- default_endogenous_parameters()
coarse <- solve_steady_state_endogenous(
  param,
  grid = make_type_grid(101L),
  old_weight = 0.90,
  tolerance = 1e-6,
  residual_tolerance = 1e-4
)
solution <- solve_steady_state_endogenous(
  param,
  initial = as.list(coarse$state),
  grid = make_type_grid(501L),
  old_weight = 0.90,
  tolerance = 1e-8,
  residual_tolerance = 2e-6
)
if (!solution$validated) {
  stop("El equilibrio refinado no supero la validacion.")
}

output_dir <- if (
    basename(project_dir) == "code" &&
      dir.exists(file.path(project_dir, "..", "results"))
) {
  normalizePath(file.path(project_dir, "..", "results"))
} else {
  file.path(project_dir, "outputs")
}
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
roots <- solution$cohort$integration$roots
summary <- data.frame(
  outcome = c(
    "capital_stock", "informal_share", "capitalization_share",
    "payg_share", "formal_labor", "informal_labor",
    "consumption_tax", "period_return", "annualized_return",
    "formal_wage", "young_consumption", "old_consumption",
    "voluntary_capital", "funded_capital",
    "payg_contributions", "payg_benefits",
    "decision_welfare", "experienced_welfare",
    "first_threshold", "second_threshold",
    "max_fixed_point_residual"
  ),
  value = c(
    solution$state[["k"]],
    solution$cohort$informal_share,
    solution$cohort$capitalization_share,
    solution$cohort$payg_share,
    solution$state[["formal_labor"]],
    solution$state[["informal_labor"]],
    solution$state[["consumption_tax"]],
    solution$prices$r,
    (1 + solution$prices$r)^(1 / 40) - 1,
    solution$prices$wage_formal,
    solution$cohort$consumption_young,
    solution$cohort$consumption_old_next,
    solution$cohort$voluntary_capital,
    solution$cohort$funded_capital,
    solution$cohort$payg_contributions,
    solution$cohort$payg_benefits_next,
    solution$cohort$decision_welfare,
    solution$cohort$experienced_welfare,
    roots[1L],
    roots[2L],
    max(abs(solution$fixed_point_residuals))
  )
)
write.csv(
  summary,
  file.path(output_dir, "endogenous_steady_state.csv"),
  row.names = FALSE
)
write.csv(
  solution$cohort$micro,
  file.path(output_dir, "endogenous_micro_choices.csv"),
  row.names = FALSE
)
saveRDS(solution, file.path(output_dir, "endogenous_solution.rds"))

calibration_targets <- becerra2026_targets(param)
calibration_moments <- becerra2026_moments(solution)
calibration_validation <- becerra2026_validation(solution, param)
calibration_comparison <- rbind(
  data.frame(
    use = "calibrated",
    moment = names(calibration_targets),
    target = unname(calibration_targets),
    model = unname(calibration_moments),
    stringsAsFactors = FALSE
  ),
  data.frame(
    use = "external_validation",
    moment = calibration_validation$moment,
    target = calibration_validation$target,
    model = calibration_validation$model,
    stringsAsFactors = FALSE
  )
)
calibration_comparison$gap <-
  calibration_comparison$model - calibration_comparison$target
write.csv(
  calibration_comparison,
  file.path(output_dir, "becerra2026_moment_comparison.csv"),
  row.names = FALSE
)

parameter_provenance <- data.frame(
  parameter = c(
    "tau_pension", "funded_account_rate",
    "payg_notional_account_rate", "payg_refund_rate",
    "payg_refund_annual_return", "replacement_rate",
    "A_informal", "payg_eligibility_intercept",
    "payg_eligibility_slope", "payg_benefit_multiplier",
    "payg_minimum_benefit", "g"
  ),
  value = c(
    param$tau_pension, param$funded_account_rate,
    param$payg_notional_account_rate, param$payg_refund_rate,
    param$payg_refund_annual_return, param$replacement_rate,
    param$A_informal, param$payg_eligibility_intercept,
    param$payg_eligibility_slope, param$payg_benefit_multiplier,
    param$payg_minimum_benefit, param$g
  ),
  status = c(
    rep("direct_legal_input", 6L),
    rep("calibrated_to_comparable_moment", 2L),
    "retained_shape_normalization", "selected_for_interior_branch",
    "retained_internal_units", "retained_structural_parameter"
  ),
  note = c(
    "Total statutory contribution rate",
    "Credited to the individual retirement account",
    "Credited to the PAYG notional account",
    "Full nominal refund when ineligible",
    "Refund adjusted only for inflation",
    "Approximate statutory rate at one minimum wage",
    "Matches the long-run formality anchor",
    "Matches the contributory-eligibility proxy",
    "The source does not identify a logistic slope",
    "The source does not identify the 40-year conversion",
    "Minimum-wage units are not mapped to efficiency wages",
    "Wage growth is not balanced-growth productivity"
  ),
  stringsAsFactors = FALSE
)
write.csv(
  parameter_provenance,
  file.path(output_dir, "becerra2026_parameter_provenance.csv"),
  row.names = FALSE
)

sensitivity_values <- c(-8.90, -8.79, -8.70)
sensitivity <- do.call(rbind, lapply(
  sensitivity_values,
  function(intercept) {
    alternative <- copy_list(param)
    alternative$payg_eligibility_intercept <- intercept
    result <- solve_steady_state_endogenous(
      alternative,
      initial = as.list(solution$state),
      grid = make_type_grid(201L),
      old_weight = 0.90,
      tolerance = 2e-7,
      residual_tolerance = 1e-5
    )
    data.frame(
      eligibility_intercept = intercept,
      converged = result$converged,
      validated = result$validated,
      informal_share = result$cohort$informal_share,
      capitalization_share = result$cohort$capitalization_share,
      payg_share = result$cohort$payg_share,
      consumption_tax = result$state[["consumption_tax"]],
      first_threshold = result$cohort$integration$roots[1L],
      second_threshold = result$cohort$integration$roots[2L],
      max_residual = max(abs(result$fixed_point_residuals))
    )
  }
))
write.csv(
  sensitivity,
  file.path(output_dir, "endogenous_sensitivity.csv"),
  row.names = FALSE
)

print(solution)
print(summary)
print(calibration_comparison)
print(parameter_provenance)
print(sensitivity)
