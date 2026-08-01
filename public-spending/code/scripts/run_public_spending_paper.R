script_path <- normalizePath(
  sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
)
project_root <- normalizePath(
  file.path(dirname(script_path), "..", ".."),
  mustWork = TRUE
)
source(file.path(project_root, "code", "R", "model.R"))
source(file.path(project_root, "code", "R", "stable_path_solver.R"))
source(file.path(project_root, "code", "R", "fast_solver.R"))
source(file.path(project_root, "code", "R", "stable_path_solver_v2.R"))
source(file.path(project_root, "code", "R", "solver_entry.R"))

results_dir <- file.path(project_root, "results")
figures_dir <- file.path(project_root, "figures")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

param <- default_parameters()
param$choice_scale <- 0.30
param$horizon <- 500L
param <- calibrate_parameters(param)
baseline <- solve_steady_state(param, label = "Benchmark")
scenario_map <- data.frame(
  scenario = c(
    "Lump-sum financing",
    "Mixed financing",
    "Payroll financing"
  ),
  marginal_payroll_share = c(0, 0.5, 1)
)

steady_scenarios <- lapply(
  seq_len(nrow(scenario_map)),
  function(index) {
    solve_steady_state(
      param,
      spending_increase = param$permanent_spending_increase,
      marginal_payroll_share =
        scenario_map$marginal_payroll_share[index],
      label = scenario_map$scenario[index]
    )
  }
)
steady <- do.call(rbind, c(list(baseline), steady_scenarios))
rownames(steady) <- NULL

transitions <- lapply(
  seq_len(nrow(scenario_map)),
  function(index) {
    solve_transition(
      param,
      marginal_payroll_share =
        scenario_map$marginal_payroll_share[index],
      label = scenario_map$scenario[index],
      path_type = "step"
    )
  }
)
transition_data <- do.call(
  rbind,
  lapply(transitions, function(solution) solution$data)
)
rownames(transition_data) <- NULL

welfare <- do.call(
  rbind,
  lapply(
    seq_along(transitions),
    function(index) {
      data.frame(
        scenario = scenario_map$scenario[index],
        public_utility_weight = c(
          0,
          param$public_utility_weight
        ),
        consumption_compensation_percent = c(
          consumption_compensation(
            transitions[[index]],
            param,
            public_utility_weight = 0
          ),
          consumption_compensation(
            transitions[[index]],
            param,
            public_utility_weight =
              param$public_utility_weight
          )
        )
      )
    }
  )
)

calibration <- data.frame(
  parameter = c(
    "Capital share in formal production",
    "Annual depreciation",
    "Annual discount factor",
    "CRRA coefficient",
    "Informal productivity",
    "Discrete-choice scale",
    "Benchmark informality",
    "Benchmark government consumption / GDP",
    "Baseline payroll-financed share of G",
    "Permanent increase in G / initial GDP"
  ),
  symbol = c(
    "alpha", "delta", "beta", "sigma", "A_I", "mu",
    "ell_I", "G/Y", "lambda_0", "Delta G/Y_0"
  ),
  value = c(
    param$alpha,
    param$depreciation,
    param$beta,
    param$sigma,
    param$A_informal,
    param$choice_scale,
    param$target_informal_share,
    param$target_government_share,
    param$baseline_payroll_finance_share,
    param$permanent_spending_increase
  ),
  status = c(
    rep("Standard / provisional", 6),
    "DANE target",
    "World Bank target",
    "Provisional policy split",
    "Policy experiment"
  )
)

write.csv(
  steady,
  file.path(results_dir, "steady_state_comparison.csv"),
  row.names = FALSE
)
write.csv(
  transition_data,
  file.path(results_dir, "transition_paths.csv"),
  row.names = FALSE
)
write.csv(
  welfare,
  file.path(results_dir, "welfare_comparison.csv"),
  row.names = FALSE
)
write.csv(
  calibration,
  file.path(results_dir, "calibration.csv"),
  row.names = FALSE
)

baseline_row <- steady[steady$scenario == "Benchmark", ]
policy_rows <- steady[steady$scenario != "Benchmark", ]
policy_rows$capital_change_percent <- 100 * (
  policy_rows$capital / baseline_row$capital - 1
)
policy_rows$output_change_percent <- 100 * (
  policy_rows$output / baseline_row$output - 1
)
policy_rows$consumption_change_percent <- 100 * (
  policy_rows$consumption / baseline_row$consumption - 1
)
policy_rows$informality_change_pp <- 100 * (
  policy_rows$informal_share - baseline_row$informal_share
)
policy_rows$payroll_tax_change_pp <- 100 * (
  policy_rows$payroll_tax - baseline_row$payroll_tax
)

table_lines <- c(
  "\\begin{table}[!htbp]",
  "\\centering",
  "\\caption{Permanent public-spending increase: stationary effects}",
  "\\label{tab:steady_results}",
  "\\begin{threeparttable}",
  "\\begin{tabular}{lrrrrr}",
  "\\toprule",
  paste0(
    "Financing of additional spending & $\\Delta K$ & $\\Delta Y$",
    " & $\\Delta C$ & $\\Delta$ informality & $\\Delta\\tau_f$ \\\\"
  ),
  " & (\\%) & (\\%) & (\\%) & (p.p.) & (p.p.) \\\\",
  "\\midrule",
  vapply(
    seq_len(nrow(policy_rows)),
    function(index) {
      paste0(
        policy_rows$scenario[index], " & ",
        sprintf("%.2f", policy_rows$capital_change_percent[index]),
        " & ",
        sprintf("%.2f", policy_rows$output_change_percent[index]),
        " & ",
        sprintf("%.2f", policy_rows$consumption_change_percent[index]),
        " & ",
        sprintf("%.2f", policy_rows$informality_change_pp[index]),
        " & ",
        sprintf("%.2f", policy_rows$payroll_tax_change_pp[index]),
        " \\\\"
      )
    },
    character(1)
  ),
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}[flushleft]\\footnotesize",
  paste0(
    "\\item Notes: The increase in government purchases equals two percent",
    " of benchmark GDP. All scenarios share the same initial steady state.",
    " Lump-sum financing means that none of the additional spending is",
    " charged to the formal wage bill; mixed and payroll financing charge",
    " 50 and 100 percent, respectively."
  ),
  "\\end{tablenotes}",
  "\\end{threeparttable}",
  "\\end{table}"
)
writeLines(
  table_lines,
  file.path(results_dir, "steady_state_table.tex")
)

find_policy <- function(label) {
  policy_rows[policy_rows$scenario == label, ]
}
payroll_row <- find_policy("Payroll financing")
lump_row <- find_policy("Lump-sum financing")
mixed_row <- find_policy("Mixed financing")
macro_lines <- c(
  paste0(
    "\\newcommand{\\BenchmarkInformality}{",
    sprintf("%.1f", 100 * baseline_row$informal_share),
    "\\%}"
  ),
  paste0(
    "\\newcommand{\\BenchmarkGovernmentShare}{",
    sprintf("%.1f", 100 * baseline_row$government_share_output),
    "\\%}"
  ),
  paste0(
    "\\newcommand{\\LumpCapitalChange}{",
    sprintf("%.1f", lump_row$capital_change_percent),
    "\\%}"
  ),
  paste0(
    "\\newcommand{\\PayrollCapitalChange}{",
    sprintf("%.1f", payroll_row$capital_change_percent),
    "\\%}"
  ),
  paste0(
    "\\newcommand{\\PayrollInformalityChange}{",
    sprintf("%.1f", payroll_row$informality_change_pp),
    " percentage points}"
  ),
  paste0(
    "\\newcommand{\\MixedInformalityChange}{",
    sprintf("%.1f", mixed_row$informality_change_pp),
    " percentage points}"
  )
)
writeLines(
  macro_lines,
  file.path(results_dir, "summary_values.tex")
)

scenario_colors <- c(
  "Lump-sum financing" = "#2A6FBB",
  "Mixed financing" = "#E69F00",
  "Payroll financing" = "#C43C39"
)

png(
  file.path(figures_dir, "transition_financing.png"),
  width = 1800,
  height = 1100,
  res = 180
)
par(mfrow = c(2, 3), mar = c(4, 4.2, 2.6, 1), las = 1)
panels <- list(
  c("capital", "Capital stock"),
  c("informal_share", "Informal employment share"),
  c("payroll_tax", "Formal payroll tax"),
  c("output", "Total output"),
  c("consumption", "Private consumption"),
  c("government_share_output", "Government spending / output")
)
for (panel in panels) {
  variable <- panel[1L]
  label <- panel[2L]
  ranges <- range(transition_data[[variable]], finite = TRUE)
  plot(
    NA,
    xlim = range(transition_data$time),
    ylim = ranges,
    xlab = "Year after announcement",
    ylab = label,
    main = label
  )
  for (scenario in names(scenario_colors)) {
    selection <- transition_data$scenario == scenario
    lines(
      transition_data$time[selection],
      transition_data[[variable]][selection],
      col = scenario_colors[[scenario]],
      lwd = 2
    )
  }
  abline(
    h = baseline_row[[variable]],
    col = "gray65",
    lty = 3
  )
}
legend(
  "topright",
  legend = names(scenario_colors),
  col = scenario_colors,
  lwd = 2,
  bty = "n",
  cex = 0.75
)
dev.off()

png(
  file.path(figures_dir, "steady_state_comparison.png"),
  width = 1500,
  height = 750,
  res = 180
)
par(mfrow = c(1, 3), mar = c(7, 4.2, 2.6, 1), las = 1)
bar_labels <- c("Lump sum", "Mixed", "Payroll")
barplot(
  policy_rows$capital_change_percent,
  names.arg = bar_labels,
  col = unname(scenario_colors),
  ylab = "Percent",
  main = "Capital",
  las = 2
)
abline(h = 0, col = "gray50")
barplot(
  policy_rows$informality_change_pp,
  names.arg = bar_labels,
  col = unname(scenario_colors),
  ylab = "Percentage points",
  main = "Informality",
  las = 2
)
abline(h = 0, col = "gray50")
barplot(
  policy_rows$consumption_change_percent,
  names.arg = bar_labels,
  col = unname(scenario_colors),
  ylab = "Percent",
  main = "Private consumption",
  las = 2
)
abline(h = 0, col = "gray50")
dev.off()

cat("Public-spending paper analysis completed.\n")
print(steady)
print(welfare)
