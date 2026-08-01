args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[1]))
} else {
  normalizePath("code/scripts/run_paper.R")
}
paper_dir <- normalizePath(file.path(dirname(script_path), "..", ".."))

source(file.path(paper_dir, "code", "R", "model.R"))
source(file.path(paper_dir, "code", "R", "solver.R"))

results_dir <- file.path(paper_dir, "results")
figures_dir <- file.path(paper_dir, "figures")
output_dir <- file.path(paper_dir, "output", "pdf")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

demographics <- read.csv(
  file.path(paper_dir, "data", "demographic_path.csv")
)
reference_mix <- c(consumption = 0.45, payroll = 0.35, capital = 0.20)

calibration <- calibrate_parameters(
  default_parameters(),
  demographics$old_age_dependency[1],
  reference_mix = reference_mix
)
if (!isTRUE(calibration$converged)) {
  stop("Calibration did not match the two targets.")
}
par <- calibration$parameters

evaluate_mixes <- function(mixes, dependency_ratio, scenario_par,
                           debt_ratio = scenario_par$debt_anchor) {
  rows <- vector("list", nrow(mixes))
  start <- NULL
  for (i in seq_len(nrow(mixes))) {
    eq <- solve_stationary_equilibrium(
      dependency_ratio, mixes[i, ], scenario_par, debt_ratio, start
    )
    if (isTRUE(eq$converged)) start <- eq$z
    row <- equilibrium_row(eq)
    row$converged <- isTRUE(eq$converged)
    rows[[i]] <- row
  }
  out <- do.call(rbind, rows)
  out$is_optimum <- FALSE
  out$welfare_loss_pct <- NA_real_
  feasible <- out$converged & is.finite(out$lifetime_welfare)
  if (any(feasible)) {
    best_i <- which.max(ifelse(feasible, out$lifetime_welfare, -Inf))
    out$is_optimum[best_i] <- TRUE
    best_c <- out$composite_consumption[best_i]
    out$welfare_loss_pct[feasible] <-
      100 * (best_c - out$composite_consumption[feasible]) / best_c
  }
  out
}

refined_policy_search <- function(dependency_ratio, scenario_par,
                                  debt_ratio = scenario_par$debt_anchor) {
  coarse_mixes <- simplex_grid(0.10)
  coarse <- evaluate_mixes(
    coarse_mixes, dependency_ratio, scenario_par, debt_ratio
  )
  coarse_best <- optimal_grid_row(coarse)
  center <- as.numeric(coarse_best[1, c(
    "mix_consumption", "mix_payroll", "mix_capital"
  )])
  fine_mixes <- simplex_grid(0.025)
  keep <- apply(
    abs(sweep(fine_mixes, 2, center, "-")),
    1,
    max
  ) <= 0.125 + 1e-10
  fine <- evaluate_mixes(
    fine_mixes[keep, , drop = FALSE],
    dependency_ratio,
    scenario_par,
    debt_ratio
  )
  combined <- rbind(coarse, fine)
  key <- paste(
    round(combined$mix_consumption, 6),
    round(combined$mix_payroll, 6),
    round(combined$mix_capital, 6),
    sep = ":"
  )
  combined <- combined[!duplicated(key), , drop = FALSE]
  feasible <- combined$converged & is.finite(combined$lifetime_welfare)
  best_i <- which.max(ifelse(feasible, combined$lifetime_welfare, -Inf))
  combined$is_optimum <- FALSE
  combined$is_optimum[best_i] <- TRUE
  best_c <- combined$composite_consumption[best_i]
  combined$welfare_loss_pct <- NA_real_
  combined$welfare_loss_pct[feasible] <-
    100 * (best_c - combined$composite_consumption[feasible]) / best_c
  list(grid = combined, optimum = combined[best_i, , drop = FALSE])
}

policy_years <- c(2024, 2050, 2070)
policy_searches <- vector("list", length(policy_years))
names(policy_searches) <- as.character(policy_years)
for (year in policy_years) {
  message("Optimizing tax mix for ", year, "...")
  d <- demographics$old_age_dependency[demographics$year == year]
  policy_searches[[as.character(year)]] <- refined_policy_search(d, par)
}

policy_grids <- do.call(rbind, lapply(names(policy_searches), function(year) {
  out <- policy_searches[[year]]$grid
  out$year <- as.integer(year)
  out$debt_anchor <- par$debt_anchor
  out
}))
policy_optima <- do.call(rbind, lapply(names(policy_searches), function(year) {
  out <- policy_searches[[year]]$optimum
  out$year <- as.integer(year)
  out$debt_anchor <- par$debt_anchor
  out
}))
rownames(policy_grids) <- NULL
rownames(policy_optima) <- NULL

debt_anchors <- c(0.40, 0.55, 0.70)
debt_searches <- vector("list", length(debt_anchors))
names(debt_searches) <- format(debt_anchors, nsmall = 2)
d2050 <- demographics$old_age_dependency[demographics$year == 2050]
for (anchor in debt_anchors) {
  message("Optimizing 2050 tax mix at debt anchor ", anchor, "...")
  scenario_par <- par
  scenario_par$debt_anchor <- anchor
  debt_searches[[format(anchor, nsmall = 2)]] <-
    refined_policy_search(d2050, scenario_par, anchor)
}
debt_optima <- do.call(rbind, lapply(names(debt_searches), function(anchor) {
  out <- debt_searches[[anchor]]$optimum
  out$debt_anchor <- as.numeric(anchor)
  out$year <- 2050
  out
}))
rownames(debt_optima) <- NULL

frictionless_par <- par
frictionless_par$collection_cost_consumption <- 0
frictionless_par$collection_cost_payroll <- 0
frictionless_par$collection_cost_capital <- 0
frictionless_grid <- evaluate_mixes(
  simplex_grid(0.10), d2050, frictionless_par
)
frictionless_optimum <- optimal_grid_row(frictionless_grid)

mix_2024 <- as.numeric(policy_optima[
  policy_optima$year == 2024,
  c("mix_consumption", "mix_payroll", "mix_capital")
])
names(mix_2024) <- c("consumption", "payroll", "capital")
path_reference <- solve_demographic_path(
  demographics, reference_mix, par
)
path_optimal <- solve_demographic_path(
  demographics, mix_2024, par
)
path_reference$policy <- "Reference mix"
path_optimal$policy <- "2024-optimal mix"
conditional_paths <- rbind(path_reference, path_optimal)

calibration_table <- data.frame(
  parameter = c(
    "Formal capital share", "Annual depreciation", "Annual discount factor",
    "CES elasticity", "Informal VAT coverage", "Choice scale",
    "Old-age outlay per older adult", "Formal-work fixed cost",
    "Real-resource share of social-security outlays",
    "Real bond rate", "Trend output growth", "Debt anchor",
    "Common collection-cost coefficient"
  ),
  symbol = c(
    "alpha", "delta", "beta", "eta", "xi_I", "mu", "g_bar_o",
    "kappa", "zeta", "r_b", "gamma", "q_bar", "chi"
  ),
  value = c(
    par$alpha, par$delta, par$beta, par$eta,
    par$informal_vat_coverage, par$choice_scale,
    par$old_age_outlay, par$formal_cost, par$resource_share,
    par$bond_rate, par$output_growth, par$debt_anchor,
    par$collection_cost_consumption
  ),
  status = c(
    rep("Standard / provisional", 6),
    "Calibrated to spending target",
    "Calibrated to informality target",
    "Provisional", "Provisional", "Provisional",
    "Law 2155 of 2021", "Provisional"
  )
)

write.csv(calibration_table,
          file.path(results_dir, "calibration.csv"), row.names = FALSE)
write.csv(policy_grids,
          file.path(results_dir, "tax_mix_grid.csv"), row.names = FALSE)
write.csv(policy_optima,
          file.path(results_dir, "optimal_tax_mix_by_year.csv"),
          row.names = FALSE)
write.csv(debt_optima,
          file.path(results_dir, "debt_anchor_sensitivity.csv"),
          row.names = FALSE)
write.csv(conditional_paths,
          file.path(results_dir, "conditional_stationary_paths.csv"),
          row.names = FALSE)
write.csv(frictionless_optimum,
          file.path(results_dir, "frictionless_optimum_2050.csv"),
          row.names = FALSE)

png(
  file.path(figures_dir, "tax_simplex_2050.png"),
  width = 1800, height = 1500, res = 220
)
grid2050 <- policy_searches[["2050"]]$grid
grid2050 <- grid2050[grid2050$converged, ]
x <- grid2050$mix_payroll + 0.5 * grid2050$mix_capital
y <- sqrt(3) / 2 * grid2050$mix_capital
loss <- pmin(grid2050$welfare_loss_pct, 8)
palette <- hcl.colors(100, "YlOrRd", rev = TRUE)
color_index <- 1 + floor(98 * (loss - min(loss)) /
                           pmax(max(loss) - min(loss), 1e-12))
plot(
  x, y, type = "n", asp = 1, axes = FALSE, xlab = "", ylab = "",
  main = "Consumption-equivalent loss across tax mixes, 2050"
)
polygon(c(0, 1, 0.5), c(0, 0, sqrt(3) / 2),
        border = "grey35", lwd = 1.5)
points(x, y, pch = 19, cex = 0.9, col = palette[color_index])
best <- policy_searches[["2050"]]$optimum
xb <- best$mix_payroll + 0.5 * best$mix_capital
yb <- sqrt(3) / 2 * best$mix_capital
points(xb, yb, pch = 8, cex = 2.0, lwd = 2)
text(0, -0.04, "Consumption", adj = c(0, 1), xpd = TRUE)
text(1, -0.04, "Payroll", adj = c(1, 1), xpd = TRUE)
text(0.5, sqrt(3) / 2 + 0.035, "Capital", xpd = TRUE)
legend(
  "topright",
  legend = c("Lower loss", "Higher loss", "Optimum"),
  pch = c(19, 19, 8),
  col = c(palette[1], palette[99], "black"),
  bty = "n", cex = 0.85
)
dev.off()

png(
  file.path(figures_dir, "optimal_mix_horizons.png"),
  width = 1900, height = 1100, res = 220
)
op <- par(mfrow = c(1, 2), mar = c(4.3, 4.5, 2.4, 0.8))
mix_matrix <- t(as.matrix(policy_optima[, c(
  "mix_consumption", "mix_payroll", "mix_capital"
)]))
barplot(
  mix_matrix,
  names.arg = policy_optima$year,
  col = c("#3B82F6", "#F59E0B", "#64748B"),
  border = NA, ylim = c(0, 1),
  ylab = "Share of required revenue",
  xlab = "DANE demographic year",
  main = "Optimal revenue composition"
)
legend(
  "topright",
  legend = c("Consumption", "Payroll", "Capital"),
  fill = c("#3B82F6", "#F59E0B", "#64748B"),
  bty = "n", cex = 0.85
)
tax_matrix <- t(as.matrix(policy_optima[, c(
  "tax_consumption", "tax_payroll", "tax_capital"
)])) * 100
matplot(
  policy_optima$year, t(tax_matrix), type = "b", pch = c(16, 17, 15),
  lty = 1, lwd = 2,
  col = c("#3B82F6", "#F59E0B", "#64748B"),
  xlab = "DANE demographic year", ylab = "Tax rate (percent)",
  main = "Tax rates at the optimum"
)
legend(
  "topleft",
  legend = c("Consumption", "Payroll", "Capital"),
  col = c("#3B82F6", "#F59E0B", "#64748B"),
  pch = c(16, 17, 15), lty = 1, bty = "n", cex = 0.85
)
par(op)
dev.off()

png(
  file.path(figures_dir, "conditional_aging_paths.png"),
  width = 1900, height = 1550, res = 220
)
op <- par(mfrow = c(2, 2), mar = c(4.0, 4.5, 2.3, 0.8))
cols <- c("#475569", "#2563EB")
plot(
  path_reference$year, 100 * path_reference$social_security_share,
  type = "l", lwd = 2.2, col = cols[1],
  xlab = "Year", ylab = "Percent of GDP",
  main = "Social-security outlays and revenue"
)
lines(path_reference$year, 100 * path_reference$revenue_share,
      lwd = 2.2, col = "#DC2626")
legend(
  "topleft", legend = c("Social-security outlays", "Required revenue"),
  col = c(cols[1], "#DC2626"), lty = 1, lwd = 2.2, bty = "n", cex = 0.82
)
matplot(
  path_optimal$year,
  100 * as.matrix(path_optimal[, c(
    "tax_consumption", "tax_payroll", "tax_capital"
  )]),
  type = "l", lty = 1, lwd = 2.2,
  col = c("#2563EB", "#F59E0B", "#64748B"),
  xlab = "Year", ylab = "Percent",
  main = "Tax rates under the fixed 2024-optimal mix"
)
legend(
  "topleft", legend = c("Consumption", "Payroll", "Capital"),
  col = c("#2563EB", "#F59E0B", "#64748B"),
  lty = 1, lwd = 2.2, bty = "n", cex = 0.82
)
plot(
  path_reference$year, 100 * path_reference$informality,
  type = "l", lwd = 2.2, col = cols[1],
  ylim = range(100 * c(path_reference$informality, path_optimal$informality)),
  xlab = "Year", ylab = "Percent of workers",
  main = "Labor informality"
)
lines(path_optimal$year, 100 * path_optimal$informality,
      lwd = 2.2, col = cols[2])
legend(
  "topleft", legend = c("Reference mix", "2024-optimal mix"),
  col = cols, lty = 1, lwd = 2.2, bty = "n", cex = 0.82
)
plot(
  path_reference$year, path_reference$capital,
  type = "l", lwd = 2.2, col = cols[1],
  ylim = range(c(path_reference$capital, path_optimal$capital)),
  xlab = "Year", ylab = "Model units",
  main = "Conditional stationary capital"
)
lines(path_optimal$year, path_optimal$capital,
      lwd = 2.2, col = cols[2])
legend(
  "bottomleft", legend = c("Reference mix", "2024-optimal mix"),
  col = cols, lty = 1, lwd = 2.2, bty = "n", cex = 0.82
)
par(op)
dev.off()

png(
  file.path(figures_dir, "debt_anchor_2050.png"),
  width = 1900, height = 1050, res = 220
)
op <- par(mfrow = c(1, 2), mar = c(4.2, 4.5, 2.3, 0.8))
plot(
  100 * debt_optima$debt_anchor,
  100 * debt_optima$revenue_share,
  type = "b", pch = 16, lwd = 2.2, col = "#DC2626",
  xlab = "Debt anchor (percent of GDP)",
  ylab = "Required revenue (percent of GDP)",
  main = "Long-run fiscal cost"
)
abline(v = 55, lty = 2, col = "grey50")
plot(
  100 * debt_optima$debt_anchor,
  debt_optima$composite_consumption,
  type = "b", pch = 16, lwd = 2.2, col = "#2563EB",
  xlab = "Debt anchor (percent of GDP)",
  ylab = "Composite consumption",
  main = "Long-run private consumption"
)
abline(v = 55, lty = 2, col = "grey50")
par(op)
dev.off()

fmt_pct <- function(x, digits = 1) sprintf(paste0("%.", digits, "f"), 100 * x)
fmt_num <- function(x, digits = 3) sprintf(paste0("%.", digits, "f"), x)

table_lines <- c(
  "\\begin{tabular}{lrrrrrr}",
  "\\toprule",
  "Year & $\\omega_c$ & $\\omega_p$ & $\\omega_k$ & $\\tau_c$ & $\\tau_p$ & Informality \\\\",
  "\\midrule"
)
for (i in seq_len(nrow(policy_optima))) {
  x <- policy_optima[i, ]
  table_lines <- c(
    table_lines,
    sprintf(
      "%d & %s & %s & %s & %s & %s & %s \\\\",
      x$year,
      fmt_pct(x$mix_consumption),
      fmt_pct(x$mix_payroll),
      fmt_pct(x$mix_capital),
      fmt_pct(x$tax_consumption),
      fmt_pct(x$tax_payroll),
      fmt_pct(x$informality)
    )
  )
}
table_lines <- c(table_lines, "\\bottomrule", "\\end{tabular}")
writeLines(table_lines, file.path(results_dir, "optimal_mix_table.tex"))

opt2024 <- policy_optima[policy_optima$year == 2024, ]
opt2050 <- policy_optima[policy_optima$year == 2050, ]
opt2070 <- policy_optima[policy_optima$year == 2070, ]
macros <- c(
  sprintf("\\newcommand{\\DependencyTwentyTwentyFour}{%s\\%%}", fmt_pct(demographics$old_age_dependency[demographics$year == 2024], 2)),
  sprintf("\\newcommand{\\DependencyTwentyFifty}{%s\\%%}", fmt_pct(demographics$old_age_dependency[demographics$year == 2050], 2)),
  sprintf("\\newcommand{\\DependencyTwentySeventy}{%s\\%%}", fmt_pct(demographics$old_age_dependency[demographics$year == 2070], 2)),
  sprintf("\\newcommand{\\BenchmarkInformality}{%s\\%%}", fmt_pct(calibration$equilibrium$informality, 1)),
  sprintf("\\newcommand{\\BenchmarkOutlayShare}{%s\\%%}", fmt_pct(calibration$equilibrium$social_security_share, 1)),
  sprintf("\\newcommand{\\OptimalConsumptionShareTwentyTwentyFour}{%s\\%%}", fmt_pct(opt2024$mix_consumption, 1)),
  sprintf("\\newcommand{\\OptimalPayrollShareTwentyTwentyFour}{%s\\%%}", fmt_pct(opt2024$mix_payroll, 1)),
  sprintf("\\newcommand{\\OptimalConsumptionShareTwentyFifty}{%s\\%%}", fmt_pct(opt2050$mix_consumption, 1)),
  sprintf("\\newcommand{\\OptimalPayrollShareTwentyFifty}{%s\\%%}", fmt_pct(opt2050$mix_payroll, 1)),
  sprintf("\\newcommand{\\OptimalConsumptionShareTwentySeventy}{%s\\%%}", fmt_pct(opt2070$mix_consumption, 1)),
  sprintf("\\newcommand{\\OptimalPayrollShareTwentySeventy}{%s\\%%}", fmt_pct(opt2070$mix_payroll, 1)),
  sprintf("\\newcommand{\\OptimalConsumptionTaxTwentyFifty}{%s\\%%}", fmt_pct(opt2050$tax_consumption, 1)),
  sprintf("\\newcommand{\\OptimalPayrollTaxTwentyFifty}{%s\\%%}", fmt_pct(opt2050$tax_payroll, 1)),
  sprintf("\\newcommand{\\OptimalInformalityTwentyFifty}{%s\\%%}", fmt_pct(opt2050$informality, 1)),
  sprintf("\\newcommand{\\OptimalCapitalTwentyFifty}{%s}", fmt_num(opt2050$capital, 3)),
  sprintf("\\newcommand{\\FrictionlessConsumptionShare}{%s\\%%}", fmt_pct(frictionless_optimum$mix_consumption, 1))
)
writeLines(macros, file.path(results_dir, "summary_values.tex"))

message("Paper results and figures generated in ", paper_dir)
