args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- if (length(file_arg) == 1L) {
  dirname(normalizePath(file_arg))
} else {
  getwd()
}
project_dir <- normalizePath(file.path(script_dir, ".."))
source(file.path(project_dir, "R", "load_endogenous_model.R"))

repository_dir <- normalizePath(file.path(project_dir, ".."))
results_dir <- file.path(repository_dir, "results")
figures_dir <- file.path(repository_dir, "figures")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# Long-run legal anchor: 55 percent of annual GDP and a 0.2 percent
# primary structural surplus. Public debt is a fiscal liability held outside
# the household portfolio in the central experiment.
fiscal_closure <- make_fiscal_steady_state_closure(
  debt_to_annual_gdp = 0.55, primary_balance_to_gdp = 0.002,
  years_per_period = 40, domestic_debt_share = 0
)

initial_param <- default_endogenous_parameters()
initial_coarse <- solve_steady_state_endogenous(
  initial_param,
  grid = make_type_grid(101L),
  old_weight = 0.90,
  tolerance = 1e-6,
  residual_tolerance = 1e-4,
  fiscal_closure = fiscal_closure
)
initial <- solve_steady_state_endogenous(
  initial_param,
  initial = as.list(initial_coarse$state),
  grid = make_type_grid(201L),
  old_weight = 0.90,
  tolerance = 2e-7,
  residual_tolerance = 1e-5,
  fiscal_closure = fiscal_closure
)
if (!initial$validated) stop("El estado inicial de Ley 100 no fue validado.")

pillar_param <- make_law2381_parameters(initial)
pillar_coarse <- solve_steady_state_endogenous(
  pillar_param,
  initial = as.list(initial$state),
  grid = make_type_grid(101L),
  old_weight = 0.90,
  tolerance = 1e-6,
  residual_tolerance = 1e-4,
  fiscal_closure = fiscal_closure
)
pillar <- solve_steady_state_endogenous(
  pillar_param,
  initial = as.list(pillar_coarse$state),
  grid = make_type_grid(201L),
  old_weight = 0.90,
  tolerance = 2e-7,
  residual_tolerance = 1e-5,
  fiscal_closure = fiscal_closure
)
if (!pillar$validated) stop("El estado contrafactual de pilares no fue validado.")

stress_param <- make_law2381_parameters(
  initial,
  payg_benefit_multiplier = initial_param$payg_benefit_multiplier
)
stress <- solve_steady_state_endogenous(
  stress_param,
  initial = as.list(pillar$state),
  grid = make_type_grid(101L),
  old_weight = 0.92,
  tolerance = 1e-6,
  residual_tolerance = 1e-4,
  fiscal_closure = fiscal_closure
)
if (!stress$validated) stop("El estado de estres chi=2.9 no fue validado.")

funded_param <- make_fully_funded_parameters(initial)
funded_coarse <- solve_steady_state_endogenous(
  funded_param,
  initial = as.list(initial$state),
  grid = make_type_grid(101L),
  old_weight = 0.90,
  tolerance = 1e-6,
  residual_tolerance = 1e-4,
  fiscal_closure = fiscal_closure
)
funded <- solve_steady_state_endogenous(
  funded_param,
  initial = as.list(funded_coarse$state),
  grid = make_type_grid(201L),
  old_weight = 0.90,
  tolerance = 2e-7,
  residual_tolerance = 1e-5,
  fiscal_closure = fiscal_closure
)
if (!funded$validated) stop("El estado totalmente fondeado no fue validado.")
saveRDS(
  list(
    initial_param = initial_param, initial = initial,
    pillar_param = pillar_param, pillar = pillar,
    stress_param = stress_param, stress = stress,
    funded_param = funded_param, funded = funded,
    fiscal_closure = fiscal_closure
  ),
  file.path(results_dir, "policy_steady_states.rds")
)

fiscal_rule <- make_fiscal_rule(
  pillar,
  initial_debt = initial$fiscal$public_debt,
  tax_lower = -0.20,
  tax_upper = 0.90
)

transition_seed <- NULL
saved_transition_file <- file.path(results_dir, "policy_transition.rds")
if (file.exists(saved_transition_file)) {
  saved_transition <- readRDS(saved_transition_file)
  if (inherits(saved_transition, "pension_policy_transition") &&
      nrow(saved_transition$path) == 21L &&
      "public_debt" %in% names(saved_transition$path) &&
      !is.null(saved_transition$fiscal_rule$target_debt_to_annual_gdp) &&
      abs(saved_transition$fiscal_rule$target_debt_to_annual_gdp - 0.55) <
        1e-10) {
    transition_seed <- saved_transition$path
  }
}

transition <- solve_transition_endogenous(
  initial_param = initial_param,
  final_param = pillar_param,
  periods = 20L,
  initial_solution = initial,
  final_solution = pillar,
  initial_path = transition_seed,
  grid = make_type_grid(101L),
  old_weight = 0.75,
  fiscal_old_weight = 0.50,
  tolerance = 1e-4,
  terminal_tolerance = 1e-3,
  fiscal_tolerance = 2e-4,
  fiscal_rule = fiscal_rule,
  max_iterations = 400L
)
if (!transition$validated || transition$dynamic_outcome != "steady_state") {
  stop("La trayectoria central no converge al estado estacionario validado.")
}

funded_fiscal_rule <- make_fiscal_rule(
  funded,
  initial_debt = initial$fiscal$public_debt,
  tax_lower = -0.20,
  tax_upper = 0.90
)
funded_transition_seed <- NULL
saved_funded_transition_file <- file.path(
  results_dir, "policy_funded_transition.rds"
)
if (file.exists(saved_funded_transition_file)) {
  saved_funded_transition <- readRDS(saved_funded_transition_file)
  if (inherits(saved_funded_transition, "pension_policy_transition") &&
      nrow(saved_funded_transition$path) == 21L &&
      saved_funded_transition$final_solution$param$regime == "funded" &&
      abs(
        saved_funded_transition$fiscal_rule$target_debt_to_annual_gdp - 0.55
      ) < 1e-10) {
    funded_transition_seed <- saved_funded_transition$path
  }
}
funded_transition <- solve_transition_endogenous(
  initial_param = initial_param,
  final_param = funded_param,
  periods = 20L,
  initial_solution = initial,
  final_solution = funded,
  initial_path = funded_transition_seed,
  grid = make_type_grid(101L),
  old_weight = 0.75,
  fiscal_old_weight = 0.50,
  tolerance = 1e-4,
  terminal_tolerance = 1e-3,
  fiscal_tolerance = 2e-4,
  fiscal_rule = funded_fiscal_rule,
  max_iterations = 400L
)
if (!funded_transition$validated ||
    funded_transition$dynamic_outcome != "steady_state") {
  stop("La trayectoria fondeada no converge al estado estacionario validado.")
}

funded_path <- funded_transition$path
funded_path$payg_balance <-
  funded_path$payg_contributions - funded_path$payg_benefits



reporting_cutoff <- max(transition$path$cohort)
funded_path$reported <- funded_path$cohort <= reporting_cutoff
funded_path$phase <- ifelse(
  funded_path$cohort == 0L,
  "Law 100 benchmark",
  ifelse(
    funded_path$cohort == 1L,
    "first reform cohort",
    ifelse(
      funded_path$cohort == reporting_cutoff,
      "funded steady state",
      "transition"
    )
  )
)
path <- transition$path
path$reported <- path$cohort <= reporting_cutoff
path$phase <- ifelse(
  path$cohort == 0L,
  "Law 100 benchmark",
  ifelse(
    path$cohort == 1L,
    "first reform cohort",
    ifelse(path$cohort == reporting_cutoff, "pillar steady state", "transition")
  )
)
path$payg_balance <- path$payg_contributions - path$payg_benefits

micro <- build_policy_microdata(transition)
distribution <- summarize_policy_distributions(micro, initial_param)
distribution$reported <- distribution$cohort <= reporting_cutoff

funded_micro <- build_policy_microdata(funded_transition)
funded_distribution <- summarize_policy_distributions(
  funded_micro, initial_param
)
funded_distribution$reported <-
  funded_distribution$cohort <= reporting_cutoff

scenario_row <- function(label, solution, chi) {
  data.frame(
    scenario = label,
    benefit_multiplier = chi,
    capital = solution$state[["k"]],
    output = solution$fiscal$output,
    consumption_tax = solution$state[["consumption_tax"]],
    public_debt_to_annual_gdp = solution$fiscal$debt_to_annual_gdp,
    primary_balance_to_gdp = solution$fiscal$primary_balance_to_gdp,
    gross_household_assets = solution$cohort$capital_next,
    voluntary_capital = solution$cohort$voluntary_capital,
    funded_capital = solution$cohort$funded_capital,
    funded_share_of_assets = solution$cohort$funded_capital /
      solution$cohort$capital_next,
    government_spending = solution$budget$spending,
    non_consumption_revenue = solution$budget$non_consumption_revenue,
    consumption_tax_base = solution$budget$consumption_tax_base,
    public_pension_revenue = solution$cohort$pension_public_revenue,
    formal_wage = solution$prices$wage_formal,
    informal_share = solution$cohort$informal_share,
    funded_only_share = solution$cohort$funded_only_share,
    payg_only_share = solution$cohort$payg_only_share,
    mixed_share = solution$cohort$mixed_share,
    solidarity_beneficiary_share =
      solution$cohort$solidarity_beneficiary_share,
    payg_contributions = solution$cohort$payg_contributions,
    payg_benefits = solution$cohort$payg_benefits_next,
    payg_balance = solution$cohort$payg_contributions - solution$cohort$payg_benefits_next,
    experienced_welfare = solution$cohort$experienced_welfare,
    max_residual = max(abs(solution$fixed_point_residuals)),
    stringsAsFactors = FALSE
  )
}
scenario_summary <- rbind(
  scenario_row("Law 100 competition", initial, initial_param$payg_benefit_multiplier),
  scenario_row("Law 2381 pillars - central", pillar, 1.0),
  scenario_row("Law 2381 pillars - chi 2.9 stress", stress, 2.9),
  scenario_row("Fully funded formal system", funded, 0)
)

parameter_crosswalk <- data.frame(
  parameter = c(
    "total contribution rate", "PAYG common-fund credit",
    "funded credit above threshold", "fully funded account credit",
    "pillar threshold in SMLMV", "replacement intercept",
    "replacement slope per SMLMV",
    "minimum contributory benefit in SMLMV",
    "solidarity benefit in SMLMV", "solidarity target mass",
    "PAYG long-period multiplier", "net-debt anchor / annual GDP",
    "structural primary balance / GDP", "years per model period",
    "public-debt persistence", "fiscal-wedge persistence",
    "domestic debt share", "fiscal-wedge lower bound",
    "fiscal-wedge upper bound"
  ),
  value = c(
    pillar_param$tau_pension, pillar_param$payg_notional_account_rate,
    pillar_param$funded_account_rate, funded_param$funded_account_rate,
    pillar_param$pillar_threshold_smlmv,
    pillar_param$pillar_replacement_intercept,
    pillar_param$pillar_replacement_slope, 1,
    pillar_param$solidarity_benefit_smlmv,
    pillar_param$solidarity_target_mass,
    pillar_param$payg_benefit_multiplier,
    fiscal_rule$target_debt_to_annual_gdp,
    fiscal_rule$target_primary_balance_to_gdp,
    fiscal_rule$years_per_period, fiscal_rule$debt_persistence,
    fiscal_rule$tax_persistence, fiscal_rule$domestic_debt_share,
    fiscal_rule$tax_lower, fiscal_rule$tax_upper
  ),
  source = c(
    rep("Law 2381, art. 23", 3L), "counterfactual design",
    "Law 2381, arts. 3 and 19", rep("Law 2381, art. 32", 3L),
    "DANE 2023 line / 2023 SMLMV",
    "DANE 2023 extreme-poverty incidence",
    "structural normalization; chi=2.9 is stress",
    "Law 2155 of 2021, art. 60", "Law 2155 of 2021, art. 60",
    "model timing", rep("reduced-form transition calibration", 2L),
    "central fiscal-incidence assumption",
    "numerical policy domain", "numerical policy domain"
  ),
  status = c(
    rep("direct legal input", 3L), "counterfactual",
    rep("direct legal input", 4L), "calibration ratio",
    "calibration proxy", "not identified", rep("direct legal input", 2L),
    "model normalization", rep("provisional calibration", 2L),
    "identifying assumption", rep("numerical bound", 2L)
  ),
  stringsAsFactors = FALSE
)

count_direction_reversals <- function(x, tolerance = 1e-8) {
  changes <- diff(x)
  directions <- sign(changes[abs(changes) > tolerance])
  if (length(directions) < 2L) return(0L)
  sum(directions[-1L] != directions[-length(directions)])
}
max_path_reversals <- function(candidate_path) {
  variables <- c(
    "capital", "public_debt", "consumption_tax", "informal_share"
  )
  max(vapply(
    candidate_path[variables], count_direction_reversals, integer(1)
  ))
}

pillar_reversals <- max_path_reversals(path)
funded_reversals <- max_path_reversals(funded_path)
validation <- data.frame(
  criterion = c(
    "initial steady state", "pillar steady state", "stress steady state",
    "funded steady state", "pillar internal fixed point",
    "pillar terminal consistency", "pillar debt-identity residual",
    "pillar fiscal-rule residual", "pillar dynamic outcome",
    "funded internal fixed point", "funded terminal consistency",
    "funded debt-identity residual", "funded fiscal-rule residual",
    "funded dynamic outcome", "pillar maximum direction reversals",
    "funded maximum direction reversals"
  ),
  value = c(
    initial$validated, pillar$validated, stress$validated, funded$validated,
    transition$converged_internal, transition$terminal_gap,
    transition$max_debt_identity_residual,
    transition$max_fiscal_rule_residual,
    transition$dynamic_outcome == "steady_state",
    funded_transition$converged_internal, funded_transition$terminal_gap,
    funded_transition$max_debt_identity_residual,
    funded_transition$max_fiscal_rule_residual,
    funded_transition$dynamic_outcome == "steady_state",
    pillar_reversals, funded_reversals
  ),
  threshold = c(
    1, 1, 1, 1, 1, 1e-3, 2e-4, 2e-4, 1,
    1, 1e-3, 2e-4, 2e-4, 1, 2, 2
  ),
  passed = c(
    initial$validated, pillar$validated, stress$validated, funded$validated,
    transition$converged_internal, transition$terminal_consistent,
    transition$max_debt_identity_residual < 2e-4,
    transition$max_fiscal_rule_residual < 2e-4,
    transition$dynamic_outcome == "steady_state",
    funded_transition$converged_internal,
    funded_transition$terminal_consistent,
    funded_transition$max_debt_identity_residual < 2e-4,
    funded_transition$max_fiscal_rule_residual < 2e-4,
    funded_transition$dynamic_outcome == "steady_state",
    pillar_reversals <= 2L, funded_reversals <= 2L
  ),
  stringsAsFactors = FALSE
)

write.csv(path, file.path(results_dir, "policy_transition_path.csv"), row.names = FALSE)
write.csv(micro, file.path(results_dir, "policy_microdata.csv"), row.names = FALSE)
write.csv(
  distribution,
  file.path(results_dir, "policy_distribution_summary.csv"),
  row.names = FALSE
)
write.csv(
  funded_path,
  file.path(results_dir, "policy_funded_transition_path.csv"),
  row.names = FALSE
)
write.csv(
  funded_micro,
  file.path(results_dir, "policy_funded_microdata.csv"),
  row.names = FALSE
)
write.csv(
  funded_distribution,
  file.path(results_dir, "policy_funded_distribution_summary.csv"),
  row.names = FALSE
)
write.csv(
  scenario_summary,
  file.path(results_dir, "policy_scenario_summary.csv"),
  row.names = FALSE
)
write.csv(
  parameter_crosswalk,
  file.path(results_dir, "policy_parameter_crosswalk.csv"),
  row.names = FALSE
)
write.csv(validation, file.path(results_dir, "policy_validation.csv"), row.names = FALSE)
saveRDS(transition, file.path(results_dir, "policy_transition.rds"))
saveRDS(
  funded_transition,
  file.path(results_dir, "policy_funded_transition.rds")
)

colors <- c(
  navy = "#264653", orange = "#E76F51", teal = "#2A9D8F",
  gold = "#E9C46A", blue = "#457B9D", gray = "#6C757D"
)
reported_path <- path[path$reported, ]
reported_distribution <- distribution[
  distribution$reported & distribution$life_stage == "young_lifetime",
]
reported_funded_path <- funded_path[funded_path$reported, ]
reported_funded_distribution <- funded_distribution[
  funded_distribution$reported &
    funded_distribution$life_stage == "young_lifetime",
]

png(
  file.path(figures_dir, "policy_transition_trajectories.png"),
  width = 1800, height = 1200, res = 180
)
par(mfrow = c(2, 3), mar = c(4, 4, 2.5, 1), oma = c(0, 0, 2, 0))
plot(reported_path$cohort, reported_path$capital, type = "o", pch = 16,
     col = colors[["navy"]], xlab = "Cohort", ylab = "Capital",
     main = "Capital stock")
matplot(
  reported_path$cohort,
  100 * reported_path[, c("informal_share", "payg_only_share", "mixed_share")],
  type = "l", lty = c(1, 2, 3), lwd = 2,
  col = c(colors[["orange"]], colors[["blue"]], colors[["teal"]]),
  xlab = "Cohort", ylab = "Percent", main = "Pension and labor status"
)
legend("right", c("Informal", "PAYG only", "PAYG + funded"),
       col = c(colors[["orange"]], colors[["blue"]], colors[["teal"]]),
       lty = c(1, 2, 3), bty = "n", cex = 0.8)
plot(reported_path$cohort, 100 * reported_path$consumption_tax,
     type = "o", pch = 16, col = colors[["orange"]],
     xlab = "Cohort", ylab = "Percent", main = "Consumption tax")
plot(reported_path$cohort, 100 * reported_path$debt_to_output,
     type = "o", pch = 17, lty = 2, col = colors[["blue"]],
     xlab = "Cohort", ylab = "Percent of annual GDP", main = "Public debt")
abline(h = 55, col = colors[["gray"]], lty = 3)
plot(
  reported_path$cohort,
  100 * reported_path$fiscal_adjustment / reported_path$output,
  type = "h", lwd = 4,
  col = ifelse(
    reported_path$fiscal_adjustment >= 0,
    colors[["teal"]], colors[["orange"]]
  ),
  xlab = "Cohort", ylab = "Percent of output",
  main = "Additional fiscal adjustment"
)
abline(h = 0, col = colors[["gray"]])
plot(reported_distribution$cohort,
     100 * reported_distribution$aggregate_welfare_cev,
     type = "o", pch = 16, col = colors[["navy"]],
     xlab = "Cohort", ylab = "Consumption-equivalent percent",
     main = "Aggregate lifetime welfare")
abline(h = 0, col = colors[["gray"]])
mtext("Transition from Law 100 competition to Law 2381 pillars", outer = TRUE, cex = 1.15)
dev.off()

png(
  file.path(figures_dir, "policy_funded_transition_trajectories.png"),
  width = 1800, height = 1200, res = 180
)
par(mfrow = c(2, 3), mar = c(4, 4, 2.5, 1), oma = c(0, 0, 2, 0))
plot(
  reported_funded_path$cohort, reported_funded_path$capital,
  type = "o", pch = 16, col = colors[["navy"]],
  xlab = "Cohort", ylab = "Capital", main = "Capital stock"
)
matplot(
  reported_funded_path$cohort,
  100 * reported_funded_path[, c("informal_share", "funded_only_share")],
  type = "l", lty = c(1, 2), lwd = 2,
  col = c(colors[["orange"]], colors[["teal"]]),
  xlab = "Cohort", ylab = "Percent", main = "Pension and labor status"
)
legend(
  "right", c("Informal", "Funded formal"),
  col = c(colors[["orange"]], colors[["teal"]]),
  lty = c(1, 2), bty = "n", cex = 0.8
)
plot(
  reported_funded_path$cohort,
  100 * reported_funded_path$consumption_tax,
  type = "o", pch = 16, col = colors[["orange"]],
  xlab = "Cohort", ylab = "Percent", main = "Fiscal wedge / dividend"
)
abline(h = 0, col = colors[["gray"]])
plot(
  reported_funded_path$cohort,
  100 * reported_funded_path$debt_to_output,
  type = "o", pch = 17, lty = 2, col = colors[["blue"]],
  xlab = "Cohort", ylab = "Percent of annual GDP", main = "Public debt"
)
abline(h = 55, col = colors[["gray"]], lty = 3)
plot(
  reported_funded_path$cohort,
  100 * reported_funded_path$fiscal_adjustment /
    reported_funded_path$output,
  type = "h", lwd = 4,
  col = ifelse(
    reported_funded_path$fiscal_adjustment >= 0,
    colors[["teal"]], colors[["orange"]]
  ),
  xlab = "Cohort", ylab = "Percent of output",
  main = "Additional fiscal adjustment"
)
abline(h = 0, col = colors[["gray"]])
plot(
  reported_funded_distribution$cohort,
  100 * reported_funded_distribution$aggregate_welfare_cev,
  type = "o", pch = 16, col = colors[["navy"]],
  xlab = "Cohort", ylab = "Consumption-equivalent percent",
  main = "Aggregate lifetime welfare"
)
abline(h = 0, col = colors[["gray"]])
mtext(
  "Transition from Law 100 competition to a fully funded formal system",
  outer = TRUE, cex = 1.15
)
dev.off()

png(
  file.path(figures_dir, "policy_reform_comparison.png"),
  width = 1600, height = 1100, res = 180
)
par(mfrow = c(2, 2), mar = c(4, 4, 2.8, 1))
plot(
  reported_path$cohort, reported_path$capital,
  type = "l", lwd = 2.5, col = colors[["orange"]],
  ylim = range(reported_path$capital, reported_funded_path$capital),
  xlab = "Cohort", ylab = "Capital", main = "Capital stock"
)
lines(
  reported_funded_path$cohort, reported_funded_path$capital,
  lwd = 2.5, col = colors[["teal"]], lty = 2
)
legend(
  "topleft", c("Pillars", "Fully funded"),
  col = c(colors[["orange"]], colors[["teal"]]),
  lty = c(1, 2), lwd = 2.5, bty = "n"
)
plot(
  reported_path$cohort, 100 * reported_path$consumption_tax,
  type = "l", lwd = 2.5, col = colors[["orange"]],
  ylim = range(
    100 * reported_path$consumption_tax,
    100 * reported_funded_path$consumption_tax
  ),
  xlab = "Cohort", ylab = "Percent", main = "Fiscal wedge / dividend"
)
lines(
  reported_funded_path$cohort,
  100 * reported_funded_path$consumption_tax,
  lwd = 2.5, col = colors[["teal"]], lty = 2
)
abline(h = 0, col = colors[["gray"]])
plot(
  reported_path$cohort, 100 * reported_path$informal_share,
  type = "l", lwd = 2.5, col = colors[["orange"]],
  ylim = range(
    100 * reported_path$informal_share,
    100 * reported_funded_path$informal_share
  ),
  xlab = "Cohort", ylab = "Percent", main = "Informal employment"
)
lines(
  reported_funded_path$cohort,
  100 * reported_funded_path$informal_share,
  lwd = 2.5, col = colors[["teal"]], lty = 2
)
plot(
  reported_distribution$cohort,
  100 * reported_distribution$aggregate_welfare_cev,
  type = "l", lwd = 2.5, col = colors[["orange"]],
  ylim = range(
    100 * reported_distribution$aggregate_welfare_cev,
    100 * reported_funded_distribution$aggregate_welfare_cev
  ),
  xlab = "Cohort", ylab = "Consumption-equivalent percent",
  main = "Aggregate lifetime welfare"
)
lines(
  reported_funded_distribution$cohort,
  100 * reported_funded_distribution$aggregate_welfare_cev,
  lwd = 2.5, col = colors[["teal"]], lty = 2
)
abline(h = 0, col = colors[["gray"]])
dev.off()

stationary_pillar_micro <- micro[
  micro$life_stage == "young_lifetime" & micro$cohort == reporting_cutoff,
]
stationary_funded_micro <- funded_micro[
  funded_micro$life_stage == "young_lifetime" &
    funded_micro$cohort == reporting_cutoff,
]
png(
  file.path(figures_dir, "policy_stationary_distribution_comparison.png"),
  width = 1700, height = 800, res = 180
)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
plot(
  NA,
  xlim = range(stationary_pillar_micro$wage, stationary_funded_micro$wage),
  ylim = c(0, 1), xlab = "Labor earnings (model units)",
  ylab = "Cumulative population share",
  main = "Stationary earnings distribution"
)
for (entry in list(
  list(data = stationary_pillar_micro, color = colors[["orange"]], lty = 1),
  list(data = stationary_funded_micro, color = colors[["teal"]], lty = 2)
)) {
  ordered <- order(entry$data$wage)
  lines(
    entry$data$wage[ordered],
    cumsum(entry$data$weight[ordered]) / sum(entry$data$weight),
    col = entry$color, lty = entry$lty, lwd = 2.5
  )
}
legend(
  "bottomright", c("Pillars", "Fully funded"),
  col = c(colors[["orange"]], colors[["teal"]]),
  lty = c(1, 2), lwd = 2.5, bty = "n"
)
plot(
  stationary_pillar_micro$i,
  100 * stationary_pillar_micro$welfare_cev,
  type = "l", lwd = 2.5, col = colors[["orange"]],
  ylim = range(
    100 * stationary_pillar_micro$welfare_cev,
    100 * stationary_funded_micro$welfare_cev
  ),
  xlab = "Worker type", ylab = "Consumption-equivalent percent",
  main = "Stationary intragenerational welfare"
)
lines(
  stationary_funded_micro$i,
  100 * stationary_funded_micro$welfare_cev,
  lwd = 2.5, col = colors[["teal"]], lty = 2
)
abline(h = 0, col = colors[["gray"]])
dev.off()

selected_cohorts <- c(0L, 1L, 2L, 5L, 10L, 20L)
selected_micro <- micro[
  micro$life_stage == "young_lifetime" & micro$cohort %in% selected_cohorts,
]
png(file.path(figures_dir, "policy_wage_distribution.png"),
    width = 1500, height = 1000, res = 180)
plot(NA, xlim = range(selected_micro$wage), ylim = c(0, 1),
     xlab = "Labor earnings (model units)", ylab = "Cumulative population share",
     main = "Distribution of labor earnings across cohorts")
line_colors <- c(colors[["gray"]], colors[["navy"]], colors[["orange"]],
                 colors[["teal"]], colors[["gold"]], colors[["blue"]])
for (index in seq_along(selected_cohorts)) {
  group <- selected_micro[selected_micro$cohort == selected_cohorts[index], ]
  order_group <- order(group$wage)
  lines(group$wage[order_group],
        cumsum(group$weight[order_group]) / sum(group$weight),
        col = line_colors[index], lwd = 2, lty = index)
}
legend("bottomright", paste("Cohort", selected_cohorts), col = line_colors,
       lty = seq_along(selected_cohorts), lwd = 2, bty = "n", cex = 0.82)
dev.off()

welfare_cohorts <- c(1L, 2L, 5L, 10L, 20L)
welfare_micro <- micro[
  micro$life_stage == "young_lifetime" & micro$cohort %in% welfare_cohorts,
]
png(file.path(figures_dir, "policy_welfare_distribution.png"),
    width = 1500, height = 1000, res = 180)
plot(NA, xlim = c(0, 1), ylim = range(100 * welfare_micro$welfare_cev),
     xlab = "Worker type", ylab = "Consumption-equivalent percent",
     main = "Intragenerational welfare effects")
abline(h = 0, col = colors[["gray"]])
for (index in seq_along(welfare_cohorts)) {
  group <- welfare_micro[welfare_micro$cohort == welfare_cohorts[index], ]
  lines(group$i, 100 * group$welfare_cev,
        col = line_colors[index + 1L], lwd = 2, lty = index)
}
legend("bottomright", paste("Cohort", welfare_cohorts),
       col = line_colors[2:6], lty = seq_along(welfare_cohorts),
       lwd = 2, bty = "n", cex = 0.85)
dev.off()

intergen <- distribution[
  distribution$reported & distribution$cohort >= -1L,
]
png(file.path(figures_dir, "policy_intergenerational_welfare.png"),
    width = 1700, height = 800, res = 180)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
plot(intergen$cohort, 100 * intergen$aggregate_welfare_cev,
     type = "o", pch = 16, col = colors[["navy"]],
     xlab = "Cohort (-1 is old at reform)",
     ylab = "Consumption-equivalent percent",
     main = "Aggregate and distributional welfare")
polygon(
  c(intergen$cohort, rev(intergen$cohort)),
  100 * c(intergen$welfare_p10, rev(intergen$welfare_p90)),
  col = adjustcolor(colors[["blue"]], alpha.f = 0.18), border = NA
)
lines(intergen$cohort, 100 * intergen$aggregate_welfare_cev,
      type = "o", pch = 16, col = colors[["navy"]])
abline(h = 0, col = colors[["gray"]])
plot(intergen$cohort, 100 * intergen$winner_share,
     type = "o", pch = 16, col = colors[["teal"]], ylim = c(0, 100),
     xlab = "Cohort (-1 is old at reform)", ylab = "Percent",
     main = "Share with a welfare gain")
abline(h = 50, col = colors[["gray"]], lty = 2)
dev.off()

print(initial)
print(pillar)
print(stress)
print(transition)
print(scenario_summary)
print(distribution[distribution$reported, ])
print(validation)
