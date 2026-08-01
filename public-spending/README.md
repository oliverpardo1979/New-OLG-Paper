# Public Spending, Capital Accumulation, and Informality

Parallel paper and R implementation by Oliver Pardo.

This directory is separate from the pension-reform paper in the repository
root. It studies a permanent increase in government purchases in an annual
general-equilibrium model with endogenous formal and informal employment.

The current version compares three marginal financing instruments for an
increase in government purchases equal to two percent of benchmark GDP:

1. lump-sum financing;
2. half lump-sum and half formal-payroll financing;
3. full formal-payroll financing.

All scenarios start from the same equilibrium. The benchmark matches 55.4
percent national labor informality and general-government consumption equal
to 14.72 percent of GDP. The remaining calibration is provisional.

## Main files

- `main.tex`: paper, equations, results, limitations, and algorithm appendix.
- `code/R/model.R`: equilibrium equations and stationary solver.
- `code/R/stable_path_solver.R`: stable-manifold primitives.
- `code/R/fast_solver.R`: local roots used along the equilibrium branch.
- `code/R/stable_path_solver_v2.R`: long-horizon boundary solution.
- `code/R/solver_entry.R`: common transition interface.
- `code/scripts/run_public_spending_paper.R`: tables, figures, and CSV results.
- `code/tests/test_public_spending_paper.R`: calibration, budget, Euler, and
  terminal-convergence tests.
- `results/`: reproducible numerical outputs.
- `figures/`: paper figures.

## Reproduce

From this directory:

```powershell
Rscript code/tests/test_public_spending_paper.R
Rscript code/scripts/run_public_spending_paper.R
..\tmp\tectonic\bin\tectonic.exe --keep-logs --outdir build main.tex
```

The numerical horizon is 500 years because the stable root converges slowly
under full payroll financing. The stationary comparison is the main
long-run result. The trajectories are equilibrium paths, not smoothed series.

## Interpretation

The quantitative results are mechanism exercises. They are not fiscal
forecasts for Colombia. In particular, the lump-sum instrument is a
theoretical benchmark, public spending is unproductive, the government budget
is balanced each year, and formal and informal output are perfect substitutes
in the aggregate resource constraint.
