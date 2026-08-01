# Population Aging, Social Security Spending, and Informality

Parallel paper and R implementation by Oliver Pardo. This directory is
separate from the pension-reform paper in the repository root.

The model links demography to the government budget through

```text
old-age dependency ratio = population 65+ / population 15-64
social-security outlays per working-age adult = fixed outlay per older adult
                                                x old-age dependency ratio
```

The Colombian experiment compares the DANE old-age dependency ratios for 2024
and 2050. It holds the real outlay per older adult fixed and evaluates three
marginal financing rules:

1. broad-tax benchmark: 0% of the additional outlay on formal payroll;
2. mixed financing: 50% on formal payroll;
3. payroll-heavy financing: 75% on formal payroll.

Full marginal payroll financing is checked separately and has no interior
stationary equilibrium in the current calibration.

## Main files

- `main.tex`: paper, equations, results, limitations, and algorithm appendix.
- `data/PPED-AreaSexoEdadNac-2018-2070.xlsx`: DANE national population
  projections, updated 18 July 2025.
- `data/demographic_targets.csv`: reproducible 2024 and 2050 aggregates used
  in the calibration.
- `code/R/model.R`: equilibrium equations and stationary solver.
- `code/R/stable_path_solver.R`: stable-manifold primitives.
- `code/R/fast_solver.R`: local roots used along the equilibrium branch.
- `code/R/stable_path_solver_v2.R`: long-horizon boundary solution.
- `code/R/solver_entry.R`: common transition interface.
- `code/scripts/run_public_spending_paper.R`: tables, figures, and CSV results.
- `code/tests/test_public_spending_paper.R`: demographic identity,
  calibration, budget, Euler, and terminal-convergence tests.

## Reproduce

From this directory:

```powershell
Rscript code/tests/test_public_spending_paper.R
Rscript code/scripts/run_public_spending_paper.R
..\tmp\tectonic\bin\tectonic.exe --keep-logs --outdir build main.tex
```

## Interpretation

The results are mechanism exercises, not Colombian fiscal forecasts. The model
treats social-security outlays as a real resource claim. This is a reasonable
reduced form for health and long-term care, but cash pensions are transfers.
A literal pension interpretation requires working and retired households in an
overlapping-generations model. The broad tax is represented as lump sum, the
government balances its budget each year, informal production uses labor only,
and formal and informal goods are perfect substitutes.