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
  grid = make_type_grid(201L),
  old_weight = 0.90,
  tolerance = 2e-7,
  residual_tolerance = 1e-5
)
if (!initial$validated) stop("El estado inicial de Ley 100 no fue validado.")

pillar_param <- make_law2381_parameters(initial)
pillar_coarse <- solve_steady_state_endogenous(
  pillar_param,
  initial = as.list(initial$state),
  grid = make_type_grid(101L),
  old_weight = 0.90,
  tolerance = 1e-6,
  residual_tolerance = 1e-4
)
pillar <- solve_steady_state_endogenous(
  pillar_param,
  initial = as.list(pillar_coarse$state),
  grid = make_type_grid(201L),
  old_weight = 0.90,
  tolerance = 2e-7,
  residual_tolerance = 1e-5
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
  residual_tolerance = 1e-4
)
if (!stress$validated) stop("El estado de estres chi=2.9 no fue validado.")

transition <- solve_transition_endogenous(
  initial_param = initial_param,
  final_param = pillar_param,
  periods = 16L,
  initial_solution = initial,
  final_solution = pillar,
  grid = make_type_grid(101L),
  old_weight = 0.85,
  tolerance = 2e-5,
  terminal_tolerance = 1e-2,
  budget_tolerance = 2e-4,
  cycle_tolerance = 1e-2,
  boundary_buffer = 3L,
  max_iterations = 600L
)
if (!transition$validated || transition$dynamic_outcome != "two_cycle") {
  stop("La trayectoria central no produjo el ciclo de dos periodos validado.")
}

reporting_cutoff <- max(transition$path$cohort) - 3L
path <- transition$path
path$reported <- path$cohort <= reporting_cutoff
path$phase <- ifelse(
  path$cohort == 0L,
  "Law 100 benchmark",
  ifelse(
    path$cohort == 1L,
    "first reform cohort",
    ifelse(path$cohort %% 2L == 0L, "high-tax phase", "low-tax phase")
  )
)
path$payg_balance <- path$payg_contributions - path$payg_benefits

micro <- build_policy_microdata(transition)
distribution <- summarize_policy_distributions(micro, initial_param)
distribution$reported <- distribution$cohort <= reporting_cutoff

scenario_row <- function(label, solution, chi) {
  data.frame(
    scenario = label,
    benefit_multiplier = chi,
    capital = solution$state[["k"]],
    consumption_tax = solution$state[["consumption_tax"]],
    formal_wage = solution$prices$wage_formal,
    informal_share = solution$cohort$informal_share,
    funded_only_share = solution$cohort$funded_only_share,
    payg_only_share = solution$cohort$payg_only_share,
    mixed_share = solution$cohort$mixed_share,
    solidarity_beneficiary_share =
      solution$cohort$solidarity_beneficiary_share,
    payg_contributions = solution$cohort$payg_contributions,
    payg_benefits = solution$cohort$payg_benefits_next,
    experienced_welfare = solution$cohort$experienced_welfare,
    max_residual = max(abs(solution$fixed_point_residuals)),
    stringsAsFactors = FALSE
  )
}
scenario_summary <- rbind(
  scenario_row("Law 100 competition", initial, initial_param$payg_benefit_multiplier),
  scenario_row("Law 2381 pillars - central", pillar, 1.0),
  scenario_row("Law 2381 pillars - chi 2.9 stress", stress, 2.9)
)

parameter_crosswalk <- data.frame(
  parameter = c(
    "total contribution rate", "PAYG common-fund credit",
    "funded credit above threshold", "pillar threshold in SMLMV",
    "replacement intercept", "replacement slope per SMLMV",
    "minimum contributory benefit in SMLMV",
    "solidarity benefit in SMLMV", "solidarity target mass",
    "PAYG long-period multiplier"
  ),
  value = c(
    pillar_param$tau_pension, pillar_param$payg_notional_account_rate,
    pillar_param$funded_account_rate, pillar_param$pillar_threshold_smlmv,
    pillar_param$pillar_replacement_intercept,
    pillar_param$pillar_replacement_slope, 1,
    pillar_param$solidarity_benefit_smlmv,
    pillar_param$solidarity_target_mass,
    pillar_param$payg_benefit_multiplier
  ),
  source = c(
    "Law 2381, art. 23", "Law 2381, art. 23",
    "Law 2381, art. 23", "Law 2381, arts. 3 and 19",
    "Law 2381, art. 32", "Law 2381, art. 32",
    "Law 2381, art. 32", "DANE 2023 line / 2023 SMLMV",
    "DANE 2023 extreme-poverty incidence",
    "structural normalization; chi=2.9 is stress"
  ),
  status = c(
    rep("direct legal input", 7L),
    "calibration ratio", "calibration proxy", "not identified"
  ),
  stringsAsFactors = FALSE
)

validation <- data.frame(
  criterion = c(
    "initial steady state", "pillar steady state", "stress steady state",
    "internal fixed point", "two-cycle repetition",
    "terminal steady-state consistency", "maximum fiscal residual"
  ),
  value = c(
    initial$validated, pillar$validated, stress$validated,
    transition$converged_internal, transition$cycle_gap,
    transition$terminal_gap, transition$max_budget_residual
  ),
  threshold = c(
    1, 1, 1, 1, 1e-2, 1e-2, 2e-4
  ),
  passed = c(
    initial$validated, pillar$validated, stress$validated,
    transition$converged_internal, transition$cycle_gap < 1e-2,
    transition$terminal_consistent,
    transition$max_budget_residual < 2e-4
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

colors <- c(
  navy = "#264653", orange = "#E76F51", teal = "#2A9D8F",
  gold = "#E9C46A", blue = "#457B9D", gray = "#6C757D"
)
reported_path <- path[path$reported, ]
reported_distribution <- distribution[
  distribution$reported & distribution$life_stage == "young_lifetime",
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
plot(reported_path$cohort, reported_path$formal_wage,
     type = "o", pch = 16, col = colors[["blue"]],
     xlab = "Cohort", ylab = "Model units", main = "Formal wage")
plot(reported_path$cohort, reported_path$payg_balance,
     type = "h", lwd = 4, col = ifelse(
       reported_path$payg_balance >= 0, colors[["teal"]], colors[["orange"]]
     ), xlab = "Cohort", ylab = "Contributions minus benefits",
     main = "PAYG cash balance")
abline(h = 0, col = colors[["gray"]])
plot(reported_distribution$cohort,
     100 * reported_distribution$aggregate_welfare_cev,
     type = "o", pch = 16, col = colors[["navy"]],
     xlab = "Cohort", ylab = "Consumption-equivalent percent",
     main = "Aggregate lifetime welfare")
abline(h = 0, col = colors[["gray"]])
mtext("Transition from Law 100 competition to Law 2381 pillars", outer = TRUE, cex = 1.15)
dev.off()

selected_cohorts <- c(0L, 1L, 2L, 3L, 12L, 13L)
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

welfare_cohorts <- c(1L, 2L, 3L, 12L, 13L)
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
